#ifdef HAVE_CONFIG_H
  #include <pz_config.h>
#endif

#include "TPZGenGrid2D.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZGmshReader.h>
#include <TPZVTKGeoMesh.h>
#include "Poisson/TPZMatPoisson.h" //for TPZMatLaplacian
#include <TPZBndCond.h> //for TPZBndCond
#include <pzanalysis.h> //for TPZAnalysis
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZNullMaterial.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZMultiphysicsCompMesh.h"
//#include "mixedpoisson.h"
#include "pzbuildmultiphysicsmesh.h"
//#include "TPZMixedPoissonParabolic.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"

TPZCompMesh *FluxCMesh(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMesh(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *MultiphysicCMesh(int dim, int pOrder, int *matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
TPZCompMesh *CMeshH1(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh);
void SolveProblem(TPZAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsH1(int dim, TPZAnalysis &an);
void ComputeError(TPZAnalysis &an, std::ofstream &anPostProcessFile);

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

int main(int argc, char* argv[]){
  //dimension of the problem
  constexpr int dim{2};
  constexpr int pOrder{4};
  //Materials - See the .geo file
  int matIdVec[]={1,2,3,4,5,6,7};
  //1 = Injection Well
  //2 = Production Well
  //3 = Bottom Line
  //4 = Top Line
  //5 = Left Line
  //6 = Right Line
  //7 = Domain

  //read mesh from gmsh
  TPZGeoMesh *gmesh;
  gmesh = new TPZGeoMesh();
  TPZGmshReader *reader;
  reader = new TPZGmshReader();
  reader -> GeometricGmshMesh4("../mesh/five-spot.msh",gmesh);

  //Flux mesh
  TPZCompMesh * cmeshflux= FluxCMesh(dim,pOrder,matIdVec,gmesh);

  //Pressure mesh
  TPZCompMesh * cmeshpressure= PressureCMesh(dim,pOrder,matIdVec,gmesh);

  //Multiphysics mesh
  TPZManVector< TPZCompMesh *, 2> meshvector(2);
  meshvector[0] = cmeshflux;
  meshvector[1] = cmeshpressure;
  TPZCompMesh * cmesh = MultiphysicCMesh(dim,pOrder,matIdVec,meshvector,gmesh);

  //Solve Multiphysics
  TPZAnalysis an(cmesh,true);
  SolveProblem(an,cmesh);

  //Creates H1 problem
  TPZCompMesh * cmeshH1 = CMeshH1(dim,pOrder,matIdVec,gmesh);

  //Solve H1
  TPZAnalysis anH1(cmeshH1,true);
  SolveProblem(anH1,cmeshH1);

  //Print results
  PrintResultsMultiphysic(dim,meshvector,an,cmesh);
  PrintResultsH1(dim,anH1);
  

  //Compute approximation error - same function for both H1 or multiphysics
  std::ofstream anPostProcessFile("postprocessH1.txt");
  ComputeError(anH1,anPostProcessFile);

  std::ofstream anPostProcessFile2("postprocessHdiv.txt");
  ComputeError(an,anPostProcessFile2);

  return 0;
}

//RHSfunction
constexpr int rhsPOrder{2};
const auto rhs = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u){
      const REAL &x = loc[0];
      const REAL &y = loc[1];
      u[0] = 0.;
};

//Analytical solution
constexpr int solOrder{2};
auto exactSol = [](const TPZVec<REAL> &loc,
  TPZVec<STATE>&u,
  TPZFMatrix<STATE>&gradU){
  const auto &x=loc[0];
  const auto &y=loc[1];
  const auto &d = 1.; // distance betweel injection and production wells
  u[0]= log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
  gradU(0,0) = x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2));
  gradU(1,0) = y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2));
  // gradU(2,0) = 0;//optional
};

//Flux computational mesh
TPZCompMesh *FluxCMesh(int dim, int pOrder,int *matIdVec, TPZGeoMesh *gmesh) 
{

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  TPZNullMaterial<> *mat = new TPZNullMaterial<>(matIdVec[6]);
//  mat->NStateVariables();
  cmesh->InsertMaterialObject(mat);

  cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
  mat->SetDimension(dim);
  //Insert boundary conditions
  // for(auto i = 0; i < 6; i++) {
    //TPZFMatrix<T> implements a full matrix of type T
    /* val1 and val2 are used for calculating the boundary
     * conditions. val1 goes in the matrix and val2 in the rhs.
     * for dirichlet boundary conditions, only the value of
     * val2 is used.*/
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    //dirichlet=0,neumann=1,robin=2
    constexpr int boundType{1};
    TPZBndCond * BCond = mat->CreateBC(mat, matIdVec[0], boundType, val1, val2);
    cmesh->InsertMaterialObject(BCond);

    TPZFMatrix<STATE> val3(1,1,1.);
    TPZManVector<STATE> val4(1,0.);
    TPZBndCond * BCond1 = mat->CreateBC(mat, matIdVec[1], boundType, val3, val4);
    cmesh->InsertMaterialObject(BCond1);

    TPZFMatrix<STATE> val5(1,1,1.);
    TPZManVector<STATE> val6(1,0.);
    TPZBndCond * BCond2 = mat->CreateBC(mat, matIdVec[2], boundType, val3, val4);
    TPZBndCond * BCond3 = mat->CreateBC(mat, matIdVec[3], boundType, val3, val4);
    TPZBndCond * BCond4 = mat->CreateBC(mat, matIdVec[4], boundType, val3, val4);
    TPZBndCond * BCond5 = mat->CreateBC(mat, matIdVec[5], boundType, val3, val4);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);

  // }

  cmesh->SetDefaultOrder(pOrder);
  cmesh->AutoBuild();

  // Print flux mesh
  // std::ofstream myfile("FluxMesh.txt");
  // cmesh->Print(myfile);

  return cmesh;
}

// Pressure computational mesh
TPZCompMesh *PressureCMesh(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh)
{

  //Change if not triangular elements
  bool fTriang = true;
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  TPZNullMaterial<> *mat = new TPZNullMaterial<>(matIdVec[6]);
  mat->SetDimension(dim);
//  mat->NStateVariables();
  cmesh->InsertMaterialObject(mat);

  cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetDefaultOrder(pOrder);
  cmesh->SetDimModel(dim);
  cmesh->AutoBuild();

  int ncon = cmesh->NConnects();
  for(int i=0; i<ncon; i++)
  {
      TPZConnect &newnod = cmesh->ConnectVec()[i];
      newnod.SetLagrangeMultiplier(1);
  }

  int nel = cmesh->NElements();
  for(int i=0; i<nel; i++){
      TPZCompEl *cel = cmesh->ElementVec()[i];
      TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
      celdisc->SetConstC(1.);
      celdisc->SetTrueUseQsiEta();
      if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
      {
          if(fTriang==true) celdisc->SetTotalOrderShape();
          else celdisc->SetTensorialShape();
      }
  }
  //Print pressure mesh
  // std::ofstream myfile("PressureMesh.txt");
  // cmesh->Print(myfile);

  return cmesh;
}

// Multiphysics computational mesh
TPZCompMesh *MultiphysicCMesh(int dim, int pOrder, int *matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
  gmesh->ResetReference();
  auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
  cmesh->SetDefaultOrder(pOrder);
  cmesh->SetDimModel(dim);
  // auto mat = new TPZMixedPoisson(matIdVec[6], dim);
  auto mat = new TPZMixedDarcyFlow(matIdVec[6], dim);

  mat->SetPermeabilityFunction(1.);
  // mat->SetViscosity(1.);
  // mat->SetForcingFunction(rhs,rhsPOrder); 
  cmesh->InsertMaterialObject(mat);
    
  //Boundary Conditions
  TPZFMatrix<STATE> val1(1,1,1.);
  TPZManVector<STATE> val2(1,-4.3820281971726605);
  constexpr int boundType{1};
  constexpr int boundType0{0};
  auto * BCond0 = mat->CreateBC(mat, matIdVec[0], boundType0, val1, val2);
  TPZFMatrix<STATE> val3(1,1,1.); 
  TPZManVector<STATE> val4(1,0.9288056258294208);
  auto * BCond1 = mat->CreateBC(mat, matIdVec[1], boundType0, val3, val4);
  cmesh->InsertMaterialObject(BCond1);
  cmesh->InsertMaterialObject(BCond0);

  TPZFMatrix<STATE> val5(1,1,1.);
  TPZManVector<STATE> val6(1,0.);
  TPZBndCond * BCond2 = mat->CreateBC(mat, matIdVec[2], boundType, val5, val6);
  TPZBndCond * BCond3 = mat->CreateBC(mat, matIdVec[3], boundType, val5, val6);
  TPZBndCond * BCond4 = mat->CreateBC(mat, matIdVec[4], boundType, val5, val6);
  TPZBndCond * BCond5 = mat->CreateBC(mat, matIdVec[5], boundType, val5, val6);
  cmesh->InsertMaterialObject(BCond2);
  cmesh->InsertMaterialObject(BCond3);
  cmesh->InsertMaterialObject(BCond4);
  cmesh->InsertMaterialObject(BCond5);

  TPZManVector<int> active(2,1);
  cmesh->BuildMultiphysicsSpace(active, meshvector);
  cmesh->SetAllCreateFunctionsMultiphysicElem();
  cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();

  TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
	TPZBuildMultiphysicsMesh::AddConnects(meshvector,cmesh);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);

  // Prints Multiphysics mesh
  // std::ofstream myfile("MultiPhysicsMesh.txt");
  // cmesh->Print(myfile);

  return cmesh;
}

void SolveProblem(TPZAnalysis &an, TPZCompMesh *cmesh)
{
  //sets number of threads to be used by the solver
  constexpr int nThreads{4};
  //defines storage scheme to be used for the FEM matrices
  //in this case, a symmetric skyline matrix is used
  TPZSkylineStructMatrix<STATE> matskl(cmesh);
  matskl.SetNumThreads(nThreads);
  an.SetStructuralMatrix(matskl);

  ///Setting a direct solver
  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt

  ///Setting an iterative solver
  // TPZMatrixSolver<STATE> * precond = an.BuildPreconditioner(TPZAnalysis::EBlockJacobi , true);
  // TPZCopySolve<STATE> * precond = new TPZCopySolve<STATE>( matskl.Create() );  step.ShareMatrix( *precond );
  // TPZStepSolver<STATE> * precond = new TPZStepSolver<STATE>( matskl.Create() ); step.ShareMatrix( *precond ); precond->SetJacobi(1, 0.0, 0);
  // TPZStepSolver<STATE> jac;
  // REAL tol = 1.e-10;
  // jac.SetSSOR(1,1.1,0.,0);
  // jac.ShareMatrix(step);
  // an.SetSolver(step);
  // step.SetGMRES(2000,2000, *precond, tol, 0);
  // step.SetCG(2000, *precond, tol, 0);
  an.SetSolver(step);

  //assembles the system
  an.Assemble();
	
  ///solves the system
  an.Solve();

  return;
}

void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZAnalysis &an, TPZCompMesh *cmesh)
{

  an.SetExact(exactSol,solOrder);

  TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
	TPZManVector<std::string,10> scalnames(2), vecnames(2);
    
	scalnames[0] = "Pressure";
	scalnames[1] = "ExactPressure";
	vecnames[0]= "Flux";
	vecnames[1]= "ExactFlux";

	int div = 0;
  std::string plotfile = "solution.vtk";
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
  // Print mesh properties
	// std::ofstream out("mesh.txt");
	// an.Print("nothing",out);

  return;
}

TPZCompMesh *CMeshH1(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh)
{
  //Creates cmesh object
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetAllCreateFunctionsContinuous();

  //Sets materials
  TPZMatPoisson<> *mat = new TPZMatPoisson<>(matIdVec[6],dim);
//  mat->SetPermeability(1);
  mat->SetForcingFunction(rhs,rhsPOrder);
  cmesh->InsertMaterialObject(mat);

  //Insert boundary conditions
  // for(auto i = 0; i < 6; i++)
  //   {
      //TPZFMatrix<T> implements a full matrix of type T
      /* val1 and val2 are used for calculating the boundary
       * conditions. val1 goes in the matrix and val2 in the rhs.
       * for dirichlet boundary conditions, only the value of 
       * val2 is used.*/
      TPZFMatrix<STATE> val1(1,1,1.);
      TPZManVector<STATE> val2(2,-4.3820281971726605);
      //dirichlet=0,neumann=1,robin=2
      constexpr int boundType{0};
      //TPZBndCond is a material type for boundary conditions
      TPZBndCond * BCond = mat->CreateBC(mat, matIdVec[0], boundType, val1, val2);
      TPZFMatrix<STATE> val3(1,1,1.);
      TPZManVector<STATE> val4(2,0.9288056258294208);
      TPZBndCond * BCond1 = mat->CreateBC(mat, matIdVec[1], boundType, val3, val4);
      cmesh->InsertMaterialObject(BCond);
      cmesh->InsertMaterialObject(BCond1);
    // }

  cmesh->SetDefaultOrder(pOrder);
  cmesh->AutoBuild();

  return cmesh;
}

void PrintResultsH1(int dim, TPZAnalysis &an)
{
  an.SetExact(exactSol,solOrder);
  TPZVec<std::string> scalarVars(1), vectorVars(1);
  scalarVars[0] = "Solution";
  vectorVars[0] = "Derivative";
  an.DefineGraphMesh(dim,scalarVars,vectorVars,"SolutionH1.vtk");
  constexpr int resolution{0};
  an.PostProcess(resolution,dim);	

  return;
}

void ComputeError(TPZAnalysis &an, std::ofstream &anPostProcessFile)
{
  an.SetExact(exactSol,solOrder);
  ///Calculating approximation error  
  TPZManVector<REAL,3> error;
  an.PostProcess(error,anPostProcessFile);
	
  std::cout << "\nApproximation error:\n";
  std::cout << "H1 Norm = " << error[0]<<'\n';
  std::cout << "L1 Norm = " << error[1]<<'\n'; 
  std::cout << "H1 Seminorm = " << error[2] << "\n\n";

  return;
}
