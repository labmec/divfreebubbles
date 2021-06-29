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
#include "../headers/TPZMatDivFreeBubbles.h" //THE NEW MATERIAL!
#include "Poisson/TPZMatPoisson.h" //for TPZMatLaplacian
#include "Projection/TPZL2Projection.h" //for BC in a single point
#include <TPZNullMaterial.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"// for Hdiv problem
#include <TPZBndCond.h> //for TPZBndCond
#include "TPZLinearAnalysis.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzstrmatrixor.h"
#include "pzlog.h"
#include "../headers/TPZCompElKernelHdiv.h" //THE NEW MATERIAL!
#include "../headers/TPZCompElKernelHdivBC.h" //THE NEW MATERIAL!
#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"

TPZCompMesh *FluxCMesh(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *FluxCMeshNew(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZVec<int64_t> matIdNeumann, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMesh(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMeshNew(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZGeoMesh *gmesh);
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
TPZMultiphysicsCompMesh *MultiphysicCMeshNew(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZVec<int64_t> matIdNeumann, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
TPZCompMesh *CMeshDivFreeBubbles(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZGeoMesh *gmesh);
void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysicNew(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsDivFreeBubbles(int dim, TPZLinearAnalysis &an);
void ComputeError(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);
void ComputeErrorHdiv(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);

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

    u[0] = x*x*x*y - y*y*y*x;
    gradU(0,0) = (3.*x*x*y - y*y*y);
    gradU(1,0) = (x*x*x - 3.*y*y*x);

    // u[0]= x*x-y*y;
    // gradU(0,0) = 2.*x;
    // gradU(1,0) = -2.*y;

    // const auto &d = 1.5; // distance betweel injection and production wells
    // u[0]=log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
    // gradU(0,0) = x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2));
    // gradU(1,0) = y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2));
};

auto exactSol2 = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];

    u[0] = x*x*x*y - y*y*y*x;
    gradU(0,0) = -(3.*x*x*y - y*y*y);
    gradU(1,0) = -(x*x*x - 3.*y*y*x);

    // u[0]= x*x-y*y;
    // gradU(0,0) = 2.*x;
    // gradU(1,0) = -2.*y;

    // const auto &d = 1.5; // distance betweel injection and production wells
    // u[0]=log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
    // gradU(0,0) = x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2));
    // gradU(1,0) = y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2));
};

int main(int argc, char* argv[])
{
    //dimension of the problem
    constexpr int dim{2};
    constexpr int pOrder{1};
    //Materials - See the .geo file
    TPZManVector<int64_t, 7> matIdVec={1,2,3,4,5,6,7};
    //1 = Domain
    //2 = Bottom
    //3 = Right
    //4 = Top 
    //5 = Left
    //6 = Point
    // std::set<int> 
    TPZManVector<int64_t, 2> matIdNeumann = {7};
    

#ifdef PZ_LOG
TPZLogger::InitializePZLOG();
#endif
    
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    TPZGmshReader *reader;
    reader = new TPZGmshReader();
    reader -> GeometricGmshMesh4("../mesh/1element.msh",gmesh);
    
    //.................................Hdiv.................................
    // //Flux mesh
    // TPZCompMesh * cmeshflux= FluxCMesh(dim,pOrder,matIdVec,gmesh);

    // //Pressure mesh
    // TPZCompMesh * cmeshpressure= PressureCMesh(dim,pOrder,matIdVec,gmesh);

    // //Multiphysics mesh
    // TPZManVector< TPZCompMesh *, 2> meshvector(2);
    // meshvector[0] = cmeshflux;
    // meshvector[1] = cmeshpressure;
    // auto * cmesh = MultiphysicCMesh(dim,pOrder,matIdVec,meshvector,gmesh);

    // //Solve Multiphysics
    // TPZLinearAnalysis an(cmesh,true);
    // SolveProblemDirect(an,cmesh);

    // //Print results
    // PrintResultsMultiphysic(dim,meshvector,an,cmesh);
    // std::ofstream out3("mesh_Hdiv.txt");
    // an.Print("nothing",out3);

    //.........................Div Free Bubbles NEW.........................
    //Flux mesh
    TPZCompMesh * cmeshfluxNew= FluxCMeshNew(dim,pOrder,matIdVec,matIdNeumann,gmesh);

    //Pressure mesh
    TPZCompMesh * cmeshpressureNew= PressureCMeshNew(dim,pOrder,matIdNeumann,gmesh);

    //Multiphysics mesh
    TPZManVector< TPZCompMesh *, 2> meshvectorNew(2);
    meshvectorNew[0] = cmeshfluxNew;
    meshvectorNew[1] = cmeshpressureNew;
    // meshvectorNew[2] = cmeshflux;
    // meshvectorNew[3] = cmeshpressure;
    auto * cmeshNew = MultiphysicCMeshNew(dim,pOrder,matIdVec,matIdNeumann,meshvectorNew,gmesh);

    //Solve Multiphysics
    TPZLinearAnalysis anNew(cmeshNew,true);
    SolveProblemDirect(anNew,cmeshNew);

    //Print results
    PrintResultsMultiphysicNew(dim,meshvectorNew,anNew,cmeshNew);
    std::ofstream out4("mesh_MDFB.txt");
    anNew.Print("nothing",out4);

    //...........................Div Free Bubbles...........................
    // //Creates DFB problem
    // TPZCompMesh * cmeshDFB = CMeshDivFreeBubbles(dim,pOrder,matIdVec,gmesh);

    // //Solve DFB
    // TPZLinearAnalysis anDFB(cmeshDFB,true);
    // SolveProblemDirect(anDFB,cmeshDFB);

    // // //Print results
    // PrintResultsDivFreeBubbles(dim,anDFB);
    // std::ofstream out("mesh.txt");
    // anDFB.Print("nothing",out);
    
    //...........................ERROR EVALUATION...........................
    std::ofstream anPostProcessFileMDFB("postprocessMDFB.txt");
    ComputeErrorHdiv(anNew,anPostProcessFileMDFB);

    // std::ofstream anPostProcessFileHdiv("postprocessHdiv.txt");
    // ComputeErrorHdiv(an,anPostProcessFileHdiv);

    // std::ofstream anPostProcessFileDFB("postprocessDFB.txt");
    // ComputeError(anDFB,anPostProcessFileDFB);

    return 0;
}



//Flux computational mesh
TPZCompMesh *FluxCMesh(int dim, int pOrder,TPZVec<int64_t> matIdVec, TPZGeoMesh *gmesh) 
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

    for (int i = 0; i < 5; i++)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(matIdVec[i]);
        cmesh->InsertMaterialObject(mat);
        mat->SetDimension(dim);
        mat->SetBigNumber(1.e10);
    }
    
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->AutoBuild();

    // Print flux mesh
    // std::ofstream myfile("FluxMesh.txt");
    // cmesh->Print(myfile);

    return cmesh;
}

// Pressure computational mesh
TPZCompMesh *PressureCMesh(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZGeoMesh *gmesh)
{
    gmesh->ResetReference();
    //Change if not triangular elements
    bool fTriang = false;
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(matIdVec[0]);
    mat->SetDimension(dim);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
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
    // Print pressure mesh
    std::ofstream myfile("PressureMesh.txt");
    cmesh->Print(myfile);

    return cmesh;
}

// Multiphysics computational mesh
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    auto mat = new TPZMixedDarcyFlow(matIdVec[0], dim);

    mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
        
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,0.);
    constexpr int boundType{1};
    constexpr int boundType0{0};
    auto * BCond0 = mat->CreateBC(mat, matIdVec[1], 0, val1, val4);//Bottom 
    auto * BCond1 = mat->CreateBC(mat, matIdVec[2], 0, val1, val2);//Right
    BCond0->SetForcingFunctionBC(exactSol);
    BCond1->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond0);
    auto * BCond2 = mat->CreateBC(mat, matIdVec[3], 1, val1, val2);//Top
    auto * BCond3 = mat->CreateBC(mat, matIdVec[4], 0, val1, val4);//Left
    BCond2->SetForcingFunctionBC(exactSol);
    BCond3->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    // TPZBndCond * BCond4 = mat->CreateBC(mat, matIdVec[5], 0, val1, val4);//Left
    // cmesh->InsertMaterialObject(BCond4);
    // auto *mat2 = new TPZMatPoisson<>(matIdVec[5],dim);
    // cmesh->InsertMaterialObject(mat2);

    TPZManVector<int> active(2,1);
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);

    // Prints Multiphysics mesh
    std::ofstream myfile("MultiPhysicsMesh.txt");
    cmesh->Print(myfile);

    return cmesh;
}

// Multiphysics computational mesh
TPZMultiphysicsCompMesh *MultiphysicCMeshNew(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZVec<int64_t> matIdNeumann, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    auto mat = new TPZMixedDarcyFlow(matIdVec[0], dim);
    mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
        
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,0.);
    constexpr int boundType{1};
    constexpr int boundType0{0};
    auto * BCond0 = mat->CreateBC(mat, matIdVec[1], 0, val1, val4);//Bottom 
    // auto * BCond1 = mat->CreateBC(mat, matIdVec[2], 0, val1, val2);//Right
    BCond0->SetForcingFunctionBC(exactSol);
    // BCond1->SetForcingFunctionBC(exactSol);
    // cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond0);
    auto * BCond2 = mat->CreateBC(mat, matIdVec[3], 0, val1, val2);//Top
    auto * BCond3 = mat->CreateBC(mat, matIdVec[4], 0, val1, val4);//Left
    // auto * BCond4 = mat->CreateBC(mat, matIdVec[5], 0, val1, val4);//Point
    // BCond2->SetForcingFunctionBC(exactSol);
    BCond3->SetForcingFunctionBC(exactSol);
    // BCond4->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    // cmesh->InsertMaterialObject(BCond4);

    // auto mat2 = new TPZMixedDarcyFlow(matIdVec[2], dim);
    // mat2->SetPermeabilityFunction(1.);
    // mat2 -> SetBigNumber(1.e10);
    // cmesh->InsertMaterialObject(mat2);
    
    TPZManVector<int> active(2,1);
    active[0]=1;
    active[1]=1;
    // active[2]=0;
    // active[3]=0;
    cmesh->BuildMultiphysicsSpace(active, meshvector);

    auto *mat3 = new TPZLagrangeMultiplierCS<STATE>(matIdNeumann[0], dim-1, 1);
    // TPZLagrangeMultiplier<STATE> *mat3 = new TPZLagrangeMultiplier<STATE>(matIdNeumann[0], dim-1, 1);
    cmesh->InsertMaterialObject(mat3);

    gmesh->ResetReference();
    cmesh->LoadReferences();

    //Add multiphysics interface
    for (int i = 0; i < matIdNeumann.size(); i++)
    {
        for (int jel = 0; jel < gmesh->NElements();jel++)
        {
            auto *gel = gmesh -> Element(jel);
            if (gel -> MaterialId() == matIdNeumann[i]){
                int nsides= gel->NSides();
                TPZGeoElSide gelside(gel,nsides-1);
                TPZGeoElSide neighbour = gelside.Neighbour();

                for (int kel = 0; kel < gmesh->NElements(); kel++)
                { 
                    auto *gelk = gmesh -> Element(kel);
                    if (gelk -> MaterialId() == 3){
                        int nsidesk= gelk->NSides();
                        TPZGeoElSide gelsidek(gelk,nsidesk-1);
                        TPZGeoElSide neighbourk = gelsidek.Neighbour();
                
                        TPZStack<TPZCompElSide> celstack,celstack2;
                        gelside.EqualLevelCompElementList(celstack, 0, 0);
                        gelsidek.EqualLevelCompElementList(celstack2, 0, 0);
                        std::cout << "GEL SiDE = " << gelside << std::endl;
                        gel->SetMaterialId(matIdNeumann[i]);
                        int64_t index;
                        // new TPZMultiphysicsInterfaceElement(*cmesh,gel,index);
                        new TPZMultiphysicsInterfaceElement(*cmesh,gel,index,celstack2[0],celstack[0]);
                    }
                }
            }
        }
    }

    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);








    // Prints Multiphysics mesh
    std::ofstream myfile("MultiPhysicsMeshNew.txt");
    cmesh->Print(myfile);

    //Prints computational mesh properties
    std::stringstream vtk_name;
    vtk_name    << "MultiPhysicsNew" << ".vtk";
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

    return cmesh;
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    //sets number of threads to be used by the solver
    constexpr int nThreads{0};
    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);

    //assembles the system
    an.Assemble();

    ///solves the system
    an.Solve();

    return;
}



void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    //sets number of threads to be used by the solver
    constexpr int nThreads{1};
    //defines storage scheme to be used for the FEM matrices
    //in this case, a symmetric skyline matrix is used
    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    ///Setting an iterative solver
    // TPZMatrixSolver<STATE> * precond = an.BuildPreconditioner(TPZAnalysis::EBlockJacobi , true);
    // TPZCopySolve<STATE> * precond = new TPZCopySolve<STATE>( matskl.Create() );  step.ShareMatrix( *precond );
    TPZStepSolver<STATE> * precond = new TPZStepSolver<STATE>( matskl.Create() ); step.ShareMatrix( *precond ); precond->SetJacobi(1, 0.0, 0);
    TPZStepSolver<STATE> jac;
    REAL tol = 1.e-30;
    jac.SetSSOR(1,1.1,0.,0);
    jac.ShareMatrix(step);
    step.SetGMRES(2000,2000, *precond, tol, 0);
    // step.SetCG(2000, *precond, tol, 0);
    an.SetSolver(step);

    //assembles the system
    an.Assemble();

    ///solves the system
    an.Solve();

    return;
}

void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{

    an.SetExact(exactSol,solOrder);

    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
    TPZManVector<std::string,10> scalnames(2), vecnames(2);

    scalnames[0] = "Pressure";
    scalnames[1] = "ExactPressure";
    vecnames[0]= "Flux";
    vecnames[1]= "ExactFlux";

    int div = 0;
    std::string plotfile = "solutionHdiv.vtk";
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    // Print mesh properties
    // std::ofstream out("mesh.txt");
    // an.Print("nothing",out);

    return;
}

void PrintResultsMultiphysicNew(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{

    an.SetExact(exactSol,solOrder);

    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
    TPZManVector<std::string,10> scalnames(0), vecnames(2);
    
    
    for (int64_t i = 0; i < cmesh->NElements(); i++)
    {
        auto *cel = cmesh -> Element(i);
        auto type = cel -> Type();
        int64_t index;
        
        if (type == EQuadrilateral){
        auto mat = cel -> Material();
        // auto data = cel -> Data();
        }
    } 

    // scalnames[0] = "Pressure";
    // scalnames[1] = "ExactPressure";
    vecnames[0]= "Flux";
    vecnames[1]= "ExactFlux";

    int div = 0;
    std::string plotfile = "solutionMDFB.vtk";
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    // Print mesh properties
    // std::ofstream out("mesh.txt");
    // an.Print("nothing",out);

    return;
}

void PrintResultsDivFreeBubbles(int dim, TPZLinearAnalysis &an)
{
    an.SetExact(exactSol,solOrder);
    TPZVec<std::string> vectorVars(1), scalarVars(1);
    scalarVars[0] = "Solution";
    vectorVars[0] = "Derivative";
    an.DefineGraphMesh(dim,scalarVars,vectorVars,"SolutionDFB.vtk");
    constexpr int resolution{0};
    an.PostProcess(resolution,dim);	

    return;
}

void ComputeError(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile)
{
    an.SetExact(exactSol2,solOrder);
    ///Calculating approximation error  
    TPZManVector<REAL,3> error;
    an.PostProcess(error,anPostProcessFile);
        
    std::cout << "\nApproximation error:\n";
    std::cout << "H1 Norm = " << error[0]<<'\n';
    std::cout << "L1 Norm = " << error[1]<<'\n'; 
    std::cout << "H1 Seminorm = " << error[2] << "\n\n";

    return;
}

void ComputeErrorHdiv(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile)
{
    an.SetExact(exactSol2,solOrder);
    ///Calculating approximation error  
    TPZManVector<REAL,5> error;
    an.PostProcess(error,anPostProcessFile);
        
    std::cout << "\nApproximation error:\n";
    std::cout << "H1 Norm = " << error[0]<<'\n';
    std::cout << "L1 Norm = " << error[1]<<'\n'; 
    std::cout << "H1 Seminorm = " << error[2] << "\n\n";

    return;
}


TPZCompMesh *CMeshDivFreeBubbles(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZGeoMesh *gmesh)
{
    gmesh->ResetReference();
        
    //Creates cmesh object
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(pOrder);

    //Sets materials
    auto *mat = new TPZMatDivFreeBubbles<STATE>(matIdVec[0],dim);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);

    //Insert boundary conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,0.);
    constexpr int boundType{0};
    auto * BCond = mat->CreateBC(mat, matIdVec[1], 0, val1, val4);//Bottom
    auto * BCond1 = mat->CreateBC(mat, matIdVec[2], 0, val1, val2);//Right
    BCond->SetForcingFunctionBC(exactSol);
    BCond1->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond);
    cmesh->InsertMaterialObject(BCond1);
    auto * BCond2 = mat->CreateBC(mat, matIdVec[3], 1, val1, val2);//Top
    auto * BCond3 = mat->CreateBC(mat, matIdVec[4], 0, val1, val4);//Left
    BCond2->SetForcingFunctionBC(exactSol);
    BCond3->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);

    auto *mat2 = new TPZL2Projection<>(matIdVec[5],0,1);
    cmesh->InsertMaterialObject(mat2);

    //Prints computational mesh properties
    // std::stringstream text_name;
    // std::stringstream vtk_name;
    // text_name   << "geometry" << ".txt";
    // vtk_name    << "geometry" << ".vtk";
    // std::ofstream textfile(text_name.str().c_str());
    // gmesh->Print(textfile);
    // std::ofstream vtkfile(vtk_name.str().c_str());
    // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

    
    cmesh->AutoBuild();

    // Prints DFB mesh
    std::ofstream myfile("DFBMesh.txt");
    cmesh->Print(myfile);

    return cmesh;
}

//Flux computational mesh
TPZCompMesh *FluxCMeshNew(int dim,int pOrder,TPZVec<int64_t> matIdVec,TPZVec<int64_t> matIdNeumann,TPZGeoMesh *gmesh) 
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);

    for (int i = 0; i < matIdVec.size()-1; i++)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(matIdVec[i]);
        cmesh->InsertMaterialObject(mat);
        mat->SetDimension(dim);
        mat->SetBigNumber(1.e10);
    }

    for (int64_t i = 0; i < gmesh->NElements(); i++)
    {
        auto *gel = gmesh -> Element(i);
        auto type = gel -> Type();
        int64_t index;
        auto matid = gel->MaterialId();
       
        using namespace pzgeom;
        using namespace pzshape;
        if (matid != 7){
            if (type == EPoint){
                TPZCompEl *cel = new TPZCompElKernelHDivBC<TPZShapePoint>(*cmesh,gel,index);
            } else {
                if (type == EOned){
                    TPZCompEl *cel = new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,gel,index);
                } else {
                    if (type == EQuadrilateral){
                        TPZCompEl *cel = new TPZCompElKernelHDiv<TPZShapeQuad>(*cmesh,gel,index);
                    } 
                }  
            }
        }
    } 
    

    cmesh->AutoBuild();

    // Print flux mesh
    std::ofstream myfile("FluxMesh.txt");
    cmesh->Print(myfile);

    //Prints computational mesh properties
    std::stringstream vtk_name;
    vtk_name    << "FluxNew" << ".vtk";
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

    return cmesh;
}


// Pressure computational mesh
TPZCompMesh *PressureCMeshNew(int dim, int pOrder, TPZVec<int64_t> matIdVec, TPZGeoMesh *gmesh)
{
    //TPZVec(inteiros)
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    std::set<int> matIdNeumann;

    for (int i = 0; i < matIdVec.size(); i++)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(matIdVec[i]);
        mat->SetDimension(dim);
        // mat->SetMaterialId(3);
        cmesh->InsertMaterialObject(mat);
        mat->SetBigNumber(1.e10);
        matIdNeumann.insert(matIdVec[i]); 
    }
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->AutoBuild(matIdNeumann);

    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i]; 
        newnod.SetLagrangeMultiplier(1);
    }

    // Print pressure mesh
    std::ofstream myfile("PressureMeshNew.txt");
    cmesh->Print(myfile);

    //Prints computational mesh properties
    std::stringstream vtk_name;
    vtk_name    << "PressNew" << ".vtk";
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

    return cmesh;
}