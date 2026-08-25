#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

// #include "TPZGenGrid2D.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for
#include "TPZGenGrid2D.h"
#include <TPZGmshReader.h>
#include <TPZVTKGeoMesh.h>
#include "Poisson/TPZMatPoisson.h"      //for TPZMatLaplacian
#include "Projection/TPZL2Projection.h" //for BC in a single point
#include "pzmultiphysicscompel.h"
#include <TPZNullMaterial.h>
#include <TPZNullMaterialCS.h>
#include <Elasticity/TPZHybridElasticity2D.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h" // for Hdiv problem
#include <TPZBndCond.h>                  //for TPZBndCond
#include "TPZLinearAnalysis.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h>    //symmetric skyline matrix storage
#include <pzstepsolver.h>       //for TPZStepSolver
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzstrmatrixor.h"
#include "pzlog.h"
#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"

#include "divfree_config.h"
#include "TPZMatDivFreeBubbles.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZCompElKernelHDiv.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZKernelHdivUtils.h"
#include "TPZHDivApproxCreator.h"
#include "TPZAnalyticSolution.h"

#include "TPZAlgebraicInterface.h"
#include "TPZH1ApproxCreator.h"
#include "TPZH1HybridApproxCreator.h"
#include <pzfstrmatrix.h>
#include "TPZVTKGenerator.h"

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

/// @brief Generate a mesh with a unique element
/// @return Geometric mesh
TPZGeoMesh *GenerateGMesh();

/// @param mfmesh the input multiphysics mesh
void InsertMultiphysicsMaterials(TPZH1ApproxCreator &create);
//#define CHECK_CONSTRAINTS 1
#ifdef CHECK_CONSTRAINTS
void VerifyConstraintConsistency(TPZMultiphysicsCompMesh *mfmesh, TPZH1HybridApproxCreator &create);
#endif
int volmatid = 1;
int bcmatidl = -1;
int bcmatidr = -2;
int bcmatidbt = -3;
int skelmatid = 2;
int interfacematid = 3;

int fluxorder = 3;

TElasticity2DAnalytic gElast2d;

int main(int argc, char *argv[])
{

#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif
  TPZGeoMesh *gmesh = GenerateGMesh();
  TPZH1HybridApproxCreator H1create(gmesh);
  H1create.SetHybridType(HybridizationType::EStandardSquared);
  // H1create.HybridType() = HybridizationType::ENone;
  H1create.SetProbType(ProblemType::EElastic);
  H1create.SetDefaultOrder(fluxorder);
  H1create.SetHybridizeBoundary();
  H1create.SetExtraInternalOrder(2);
  H1create.SetShouldCondense(false);
  int fluxmatid = H1create.HybridData().fLagrangeMatId;
  InsertMultiphysicsMaterials(H1create);
  #ifdef CHECK_CONSTRAINTS
  H1create.IsRigidBodySpaces() = true;
  #else
  H1create.IsRigidBodySpaces() = false;
  #endif
  TPZMultiphysicsCompMesh *mfmesh = H1create.CreateApproximationSpace();
  if(1) {
    TPZH1HybridApproxCreator::CtoMFCel geltogel;
    H1create.ComputeOrthogonalizingRestraints(*mfmesh, geltogel, H1create.HybridData());
    #ifndef CHECK_CONSTRAINTS
    H1create.HybridizeLowOrderFluxes(*mfmesh, geltogel);
    #endif
    H1create.GroupAndCondenseElements(mfmesh);
    {
      std::ofstream out("cmesh.txt");
      mfmesh->Print(out);
      std::cout << "number of low order connects " << geltogel.size() << std::endl;
      for(auto &iter : geltogel) {
        out << "gel left " << iter.first << " gel right " << iter.second << std::endl;
      }
      H1create.HybridData().Print(out);
    }
    // H1create.GroupAndCondenseElements(mfmesh);
  }
  #ifdef CHECK_CONSTRAINTS
  VerifyConstraintConsistency(mfmesh, H1create);
  return 0;
  #endif
  TPZLinearAnalysis an(mfmesh, RenumType::ENone);
  TPZFStructMatrix<> strmat(mfmesh);
  an.SetStructuralMatrix(strmat);
  an.Assemble();
  an.Solve();
  {
    std::ofstream out("cmesh_sol_hybrid_elas.txt");
    mfmesh->Print(out);
  }
  TPZVTKGenerator vtk(mfmesh, {"displacement", "sig_x", "sig_y", "tau_xy"}, "hybrid.vtk", 1, 2);
  vtk.Do();
  return 0;
}

/// @brief Generate a mesh with a unique element
/// @return Geometric mesh
TPZGeoMesh *GenerateGMesh()
{
  TPZGenGrid2D gengrid({3, 3}, {-1., -1.}, {1., 1.});
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  gmesh->SetDimension(2);
  gengrid.Read(gmesh, volmatid);
  gengrid.SetBC(gmesh, 4, bcmatidbt);
  gengrid.SetBC(gmesh, 6, bcmatidbt);
  gengrid.SetBC(gmesh, 7, bcmatidl);
  gengrid.SetBC(gmesh, 5, bcmatidr);
  REAL angle = atan2(0.6, 0.8);
  // std::cout << "cos angle " << std::cos(angle) << " sin angle " << std::sin(angle) << std::endl;
  // gengrid.RotateGeomesh(gmesh,angle,2);
  TPZGeoEl *gel = gmesh->Element(0);
  if (gel->Dimension() != 2)
    DebugStop();
  return gmesh;
}

/// @brief Insert the hybrid elastic material and boundary material
/// @param mfmesh the input multiphysics mesh
void InsertMultiphysicsMaterials(TPZH1ApproxCreator &create)
{
  gElast2d.fProblemType = TElasticity2DAnalytic::EHomogeneous;

  REAL E = 1.;
  REAL nu = 0.;
  gElast2d.gE = E;
  gElast2d.gPoisson = nu;
  auto *mat = new TPZHybridElasticity2D(volmatid, E, nu, 1., 1.);
  mat->SetForcingFunction(gElast2d.ForceFunc(), 2);
  mat->SetExactSol(gElast2d.ExactSolution(), 2);

  create.InsertMaterialObject(mat);
  TPZFNMatrix<4, REAL> val1(2, 2, 0.);
  TPZManVector<REAL, 2> val2(2, 1.);
  auto *bndbt = mat->CreateBC(mat, bcmatidbt, 0, val1, val2);
  bndbt->SetForcingFunctionBC(gElast2d.ExactSolution(), 2);

  create.InsertMaterialObject(bndbt);
  auto *bndl = mat->CreateBC(mat, bcmatidl, 0, val1, val2);
  bndl->SetForcingFunctionBC(gElast2d.ExactSolution(), 2);
  create.InsertMaterialObject(bndl);
  auto *bndr = mat->CreateBC(mat, bcmatidr, 0, val1, val2);
  bndr->SetForcingFunctionBC(gElast2d.ExactSolution(), 2);
  create.InsertMaterialObject(bndr);
}

#ifdef CHECK_CONSTRAINTS
#include "TPZElementMatrixT.h"
void VerifyConstraintConsistency(TPZMultiphysicsCompMesh *mfmesh, TPZH1HybridApproxCreator &create) {
  // compute the number of connects before semi hybridizing the mesh
  auto &meshvec = mfmesh->MeshVector();
  int64_t original_connect_count = 0;
  for(auto mesh : meshvec) {
    original_connect_count += mesh->NConnects();
  }
  int64_t nel = mfmesh->NElements();
  for(int64_t el = 0; el < nel; el++) {
    TPZCompEl *cel = mfmesh->Element(el);
    if(!cel) continue;
    TPZCondensedCompElT<STATE> *cond = dynamic_cast<TPZCondensedCompElT<STATE> *>(cel);
    if(!cond) continue;
    TPZElementMatrixT<STATE> ek,ef;
    cond->CalcStiff(ek,ef);
    // find the connect corresponding to the rigid body space
    int firstrgb = 0;
    int nrgb = 0;
    int rgbcondindex = -1;
    int64_t ncon = cond->NConnects();
    for(int64_t icon = 0; icon < ncon; icon++) {
      TPZConnect &c = cond->Connect(icon);
      if(c.LagrangeMultiplier() == 4) {
        rgbcondindex = icon;
        nrgb = c.NShape()*c.NState();
        break;
      } else {
        firstrgb += c.NShape()*c.NState();
      }
    }
    std::cout << "element el " << el << " rgb index " << rgbcondindex << " rgb eq " << firstrgb << std::endl;
    // look for the submatrices for connects with zero lagrange multiplier
    int eq_count = 0;
    for(int64_t icon = 0; icon < ncon; icon++) {
      int64_t conindex = cel->ConnectIndex(icon);
      TPZConnect &c = cel->Connect(icon);
      if(c.LagrangeMultiplier() == 1 || conindex < original_connect_count) {
        eq_count += c.NShape()*c.NState();
        continue;
      }
      int sizeblock = c.NShape()*c.NState();
      TPZFMatrix<STATE> submat;
      ek.fMat.GetSub(firstrgb,eq_count,nrgb,sizeblock,submat);
      REAL submatnorm = Norm(submat);
      std::cout << "element el " << el << " connect " << icon << " eq count " << eq_count << " submat norm " << submatnorm << std::endl;
      eq_count += c.NShape()*c.NState();
    }
  }
}
#endif
