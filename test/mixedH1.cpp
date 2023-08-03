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
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZGmshReader.h>
#include <TPZVTKGeoMesh.h>
#include "Poisson/TPZMatPoisson.h" //for TPZMatLaplacian
#include "Projection/TPZL2Projection.h" //for BC in a single point
#include "pzmultiphysicscompel.h"
#include <TPZNullMaterial.h>
#include <TPZNullMaterialCS.h>
#include <TPZBndCond.h> //for TPZBndCond
#include "TPZLinearAnalysis.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
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
#include <TPZGeoMeshTools.h>

#include "divfree_config.h"
#include "TPZMatDivFreeBubbles.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZCompElKernelHDiv.h"
#include "TPZMixedDarcyH1.h"
#include "TPZKernelHdivUtils.h"
#include "TPZAnalyticSolution.h"
#include "TPZBndCondT.h"

enum EMatid  {ENone, EDomain, EBoundary, EPont, EWrap, EIntface, EPressureHyb};

template<class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcId);
TPZCompMesh *CreateFluxCMesh(TPZGeoMesh *fGeoMesh, int fDefaultPOrder);
TPZCompMesh * CreatePressureCMesh(TPZGeoMesh *gmesh, int fDefaultPOrder);
void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder) ;
TPZCompMesh * CreateConstantSpace(TPZGeoMesh *fGeoMesh);
TPZMultiphysicsCompMesh *CreateMultiphysicsCMesh(TPZGeoMesh *fGeoMesh, int fDefaultPOrder, TPZVec<TPZCompMesh *> &meshvector);
//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

//Analytical solution
constexpr int solOrder{2};
auto exactSol = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];

    // //Nabla u = 1
    // u[0] = 0.25*(x*x+y*y);
    // gradU(0,0) = -0.5*x;
    // gradU(1,0) = -0.5*y;

    // //Nabla u = 0
    u[0] = x*x*x*y - y*y*y*x;
    gradU(0,0) = (3.*x*x*y - y*y*y);
    gradU(1,0) = (x*x*x - 3.*y*y*x);
    // u[0] = x;
    // gradU(0,0) = 1.;
    // gradU(1,0) = 0.;

    // REAL a1 = 1./4;
    // REAL alpha = M_PI/2;
    // u[0] = x*a1*cos(x*alpha)*cosh(y*alpha) + y*a1*sin(x*alpha)*sinh(y*alpha);
    // gradU(0,0) = -a1*(cosh(alpha*y)*(cos(alpha*x) - alpha*x*sin(alpha*x)) + alpha*y*cos(alpha*x)*sinh(alpha*y));
    // gradU(1,0) = -a1*(alpha*y*cosh(alpha*y)*sin(alpha*x) + (alpha*x*cos(alpha*x) + sin(alpha*x))*sinh(alpha*y));

    // u[0] = x*x - y*y;
    // gradU(0,0) = 2.*x;
    // gradU(1,0) = -2.*y;

};


int main(int argc, char* argv[])
{
    //dimension of the problem
    constexpr int dim{2};
    constexpr int pOrder{2};
      

#ifdef PZ_LOG
TPZLogger::InitializePZLOG();
#endif
    
    TPZVec<int> nDivs = {10,10};
    auto gmesh = CreateGeoMesh<pzshape::TPZShapeTriang>(nDivs, EDomain, EBoundary);
    // auto gmesh = CreateGeoMesh<pzshape::TPZShapeQuad>(nDivs, EDomain, EBoundary);

    // //read mesh from gmsh
    // TPZGeoMesh *gmesh;
    // gmesh = new TPZGeoMesh();
    // {
    //     TPZGmshReader reader;
    //     // essa interface permite voce mapear os nomes dos physical groups para
    //     // o matid que voce mesmo escolher
    //     TPZManVector<std::map<std::string,int>,4> stringtoint(4);
    //     stringtoint[2]["Surface"] = 1;
    //     stringtoint[1]["Bottom"] = 2;
    //     stringtoint[1]["Right"] = 2;
    //     stringtoint[1]["Top"] = 4;
    //     stringtoint[1]["Left"] = 2;
    //     stringtoint[0]["Point"] = 3;
    //     stringtoint[1]["Top2"] = 4;
    //     reader.SetDimNamePhysical(stringtoint);
    //     reader.GeometricGmshMesh("../mesh/1element.msh",gmesh);
    //     std::ofstream out("gmesh.vtk");
    //     TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    // }


    //............................Div Free Bubbles............................   
    TPZKernelHdivUtils<STATE> util;

    util.PrintGeoMesh(gmesh);

    //Flux mesh
    TPZCompMesh * cmeshflux = CreateFluxCMesh(gmesh,pOrder);
    // ChangeInternalOrder(cmeshflux,pOrder+1);
    std::string fluxFile = "fluxcmesh";
    util.PrintCompMesh(cmeshflux,fluxFile);

    //Pressure mesh
    TPZCompMesh * cmeshpressure = CreatePressureCMesh(gmesh,pOrder-1);
    // std::string pressureFile = "PressureCMesh";
    // util.PrintCompMesh(cmeshpressure,pressureFile);

    // TPZCompMesh * cmeshpressureC = CreateConstantSpace(gmesh);

    TPZManVector< TPZCompMesh *, 2> meshvectorNew(2);
    meshvectorNew[0] = cmeshflux;
    meshvectorNew[1] = cmeshpressure;       
    // meshvectorNew[2] = cmeshpressureC;       
    auto * cmeshNew = CreateMultiphysicsCMesh(gmesh,pOrder,meshvectorNew);
    // // Group and condense the elements
    // // createSpace.Condense(cmeshNew);
    // std::string multiphysicsFile = "MultiPhysicsMeshNew";
    // util.PrintCompMesh(cmeshNew,multiphysicsFile);
    // std::cout << "Number of equations = " << cmeshNew->NEquations() << std::endl;
    // // Solve the problem
    TPZLinearAnalysis anNew(cmeshNew,false);

    // createSpace.Solve(anNew, cmeshNew, true, false);
    util.SolveProblemDirect(anNew,cmeshNew,false,false);
    // std::cout << "Number of equations = " << cmeshNew->NEquations() << std::endl;

    anNew.SetExact(exactSol,solOrder);
    // anNew.SetExact(LaplaceExact.ExactSolution());
    //Print results
    util.PrintResultsMultiphysics(meshvectorNew,anNew,cmeshNew);

    std::ofstream out4("mesh_MDFB.txt");
    anNew.Print("nothing",out4);
    std::ofstream anPostProcessFileMDFB("postprocessMDFB.txt");
    
    // util.ComputeError(anNew,anPostProcessFileMDFB);
  
    return 0;
}




TPZCompMesh *CreateFluxCMesh(TPZGeoMesh *fGeoMesh, int fDefaultPOrder)
{   
    int fDimension = fGeoMesh->Dimension();
    fGeoMesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(fDimension);

    //Inserts Null Materials
    std::set<int> allMat={EDomain};
    std::set<int> fBCMatId = {EBoundary,EWrap};

    for (std::set<int>::iterator it=allMat.begin(); it!=allMat.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        cmesh->InsertMaterialObject(mat);
        mat->SetDimension(fDimension);
        mat->SetNStateVariables(fDimension);
    } 

    for (std::set<int>::iterator it=fBCMatId.begin(); it!=fBCMatId.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it,fDimension-1);
        cmesh->InsertMaterialObject(mat);
        mat->SetNStateVariables(fDimension);
    } 

    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->AutoBuild();


    return cmesh;
}


//Create 
template <class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcId)
{
    
    MMeshType meshType;
    int dim = tshape::Dimension;

    switch (tshape::Type())
    {
    case ETriangle:
        meshType = MMeshType::ETriangular;
        break;
    case EQuadrilateral:
        meshType = MMeshType::EQuadrilateral;
        break;
    case ETetraedro:
        meshType = MMeshType::ETetrahedral;
        break;
    case ECube:
        meshType = MMeshType::EHexahedral;
        break;
        case EPrisma:
        meshType = MMeshType::EPrismatic;
        break;
    default:
        DebugStop();
    }

    TPZManVector<REAL,3> minX = {0,0,0};
    TPZManVector<REAL,3> maxX = {1,1,1};
    int nMats = 2*dim+1;

    //all bcs share the same id
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,bcId);
    matIds[0] = volId;
    matIds[1] = EBoundary;
    matIds[2] = EBoundary;
    matIds[3] = EWrap;
    matIds[4] = EBoundary;
    
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls);
    // TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(meshType,
    //                     volId,createBoundEls, bcId);
    
    bool point = false;
    for (TPZGeoEl *gel : gmesh->ElementVec()){
        if (gel->Dimension() != dim) continue;
        if (point) continue;
        TPZGeoElSide gelside(gel,0);
        TPZGeoElSide gelside2(gel,1);
        // TPZGeoElBC gelbcWrap(gelside, EPont);
        // TPZGeoElBC gelbcWrap2(gelside2, EPont);
        point = true;
    }
    

    return gmesh;
    
}



TPZCompMesh * CreatePressureCMesh(TPZGeoMesh *fGeoMesh, int fDefaultPOrder)
{
    int fDimension = fGeoMesh->Dimension();
    
    fGeoMesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);

    TPZNullMaterial<> *mat = new TPZNullMaterial<>(EDomain);
    mat->SetDimension(fDimension);
    cmesh->InsertMaterialObject(mat);

    // TPZNullMaterial<> *matBC = new TPZNullMaterial<>(EPont,fDimension-1);
    // cmesh->InsertMaterialObject(matBC);

    // cmesh->SetAllCreateFunctionsDiscontinuous();
    // cmesh->SetDefaultOrder(0);
        
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetAllCreateFunctionsContinuous();
    // cmesh->ApproxSpace().CreateDisconnectedElements(true);
        
    cmesh->SetDimModel(fDimension);
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
        if(!celdisc) continue;
        celdisc->SetConstC(1.);
        celdisc->SetTrueUseQsiEta();
        // espera-se elemento de pressao apenas para o contorno
        auto aaa = celdisc->Reference()->Dimension();
        // if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        // {
        //     DebugStop();
        // }
    }
    // } else {

    //     cmesh->InitializeBlock();

    //     if(fSpaceType == ESemiHybrid){
    //         cmesh->SetDefaultOrder(0);
    //         cmesh->SetDimModel(fDimension-1);
    //         cmesh->SetAllCreateFunctionsDiscontinuous();
    //         // cmesh->ApproxSpace().CreateDisconnectedElements(true);
    //     } else {
    //         cmesh->SetAllCreateFunctionsContinuous();
    //         cmesh->ApproxSpace().CreateDisconnectedElements(true);
    //         cmesh->SetDefaultOrder(fDefaultPOrder);
    //         cmesh->SetDimModel(fDimension-1);
    //     }
        
    //     cmesh->AutoBuild(matIdVec);
        
    //     if (fSpaceType == ESemiHybrid){
    //         hybridizer.SemiHybridizePressure(cmesh,fDefaultPOrder,fConfig.fBCHybridMatId);
    //     }

        
    //     for(auto &newnod : cmesh->ConnectVec())
    //     {
    //         newnod.SetLagrangeMultiplier(1);
    //     }
    // }
    
    cmesh->InitializeBlock();
    return cmesh;
}

auto forcefunction = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    // //Nabla u = 1
    u[0] = 0.;
};

TPZMultiphysicsCompMesh *CreateMultiphysicsCMesh(TPZGeoMesh *fGeoMesh, int fDefaultPOrder, TPZVec<TPZCompMesh *> &meshvector)
{
    int fDimension = fGeoMesh->Dimension();

    // gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(fDimension);
    
    // eh preciso criar materiais para todos os valores referenciados no enum
    auto mat = new TPZMixedDarcyH1(EDomain,fDimension);
    mat->SetConstantPermeability(1.);
    mat->SetForcingFunction(forcefunction,4);
    cmesh->InsertMaterialObject(mat);
    int nstate = mat->NStateVariables();
    mat->SetBigNumber(1.e5);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(2,2,1.);
    TPZManVector<STATE> val2(2,0.);

    //Multiphysics mesh

    //Dirichlet Boundary Conditions
    TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, EBoundary, 0, val1, val2);
    // BCond->SetForcingFunctionBC(exactSol,4);
    cmesh->InsertMaterialObject(BCond);

    //Dirichlet Boundary Conditions
    val2[0] = 1.;
    TPZBndCondT<STATE> * BCond2 = mat->CreateBC(mat, EWrap, 0, val1, val2);
    // BCond2->SetForcingFunctionBC(exactSol,4);
    cmesh->InsertMaterialObject(BCond2);
    val2[0] = 0.;
    TPZBndCondT<STATE> * BCond3 = mat->CreateBC(mat, EPont, 0, val1, val2);
    // BCond3->SetForcingFunctionBC(exactSol,4);
    cmesh->InsertMaterialObject(BCond3);
    
    TPZManVector<int> active(meshvector.size(),1);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    
    return cmesh;
}

void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder) {

    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;

        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != cmesh->Dimension()) {//Only elements with the same dimension of the mesh
            continue;
        }
        int nc = cel->NConnects();
        //Gets the volumetric connect
        int64_t conIndex = cel->ConnectIndex(nc-1);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel) DebugStop();

        TPZConnect &c = cmesh->ConnectVec()[conIndex];
        int64_t index;
        //Sets the new connect order
        c.SetOrder(pOrder,index);

        //Gets connect information to update block size (stiffness matrix data structure)
        int64_t seqnum = c.SequenceNumber();
        int nvar = 1;
        TPZMaterial * mat = cel->Material();
        if (mat) nvar = mat->NStateVariables();
        int nshape = intel->NConnectShapeF(nc-1,pOrder);
        c.SetNShape(nshape);
        // c.SetNState(nvar);
        cmesh->Block().Set(seqnum, nvar * nshape);

        cel->SetIntegrationRule(2*pOrder);
    }
    cmesh->InitializeBlock();
}

TPZCompMesh * CreateConstantSpace(TPZGeoMesh *fGeoMesh) {
    int lagLevel = 3;
    fGeoMesh->ResetReference();
    int dim = fGeoMesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDimModel(dim);

    TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(EDomain,dim,1);
    cmesh->InsertMaterialObject(nullmat);
    
    cmesh->SetDefaultOrder(0);
    cmesh->SetAllCreateFunctionsDiscontinuous();

    cmesh->AutoBuild();

    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(lagLevel);
    }

    return cmesh;
}
