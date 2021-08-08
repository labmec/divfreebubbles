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
#include "../headers/TPZL2ProjectionCS.h" //THE NEW MATERIAL!
#include "Poisson/TPZMatPoisson.h" //for TPZMatLaplacian
#include "Projection/TPZL2Projection.h" //for BC in a single point
#include "pzmultiphysicscompel.h"
#include <TPZNullMaterial.h>
#include <TPZNullMaterialCS.h>
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
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
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
#include "TPZHybridizeHDiv.h"


TPZCompMesh *FluxCMesh(int dim, int pOrder, std::set<int> matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *FluxCMeshNew(int dim, int pOrder, std::set<int> matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMeshNew(int dim, int pOrder, std::set<int> matIdVec, TPZGeoMesh *gmesh);
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, std::set<int> matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
TPZMultiphysicsCompMesh *MultiphysicCMeshNew(int dim, int pOrder, std::set<int> &matIdVec, std::set<int> &matIdNeumann, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
TPZCompMesh *CMeshDivFreeBubbles(int dim, int pOrder, std::set<int> matIdVec, TPZGeoMesh *gmesh);
void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh);
void PrintResultsMultiphysicNew(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh);
void PrintResultsDivFreeBubbles(int dim, TPZLinearAnalysis &an);
void ComputeError(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);
void ComputeErrorHdiv(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);
void CreateMultiphysicsInterfaceElements(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh, TPZVec<TPZCompMesh *> meshvector, std::set<int> &matIdNeumann);
void HybridizerGeoMesh(TPZGeoMesh *gmesh, std::set<int> &matIdBC);

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
    // gradU(0,0) = -2.*x;
    // gradU(1,0) = 2.*y;

    // const auto &d = 1.5; // distance betweel injection and production wells
    // u[0]=log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
    // gradU(0,0) = x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2));
    // gradU(1,0) = y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2));
};

enum EMatid  {ENone, EDomain, EBottom, ERight, ETop, ELeft, EPont, EWrap, EIntfaceLeft, EIntfaceRight, EPressureHyb};

int main(int argc, char* argv[])
{
    //dimension of the problem
    constexpr int dim{2};
    constexpr int pOrder{3};
      

#ifdef PZ_LOG
TPZLogger::InitializePZLOG();
#endif
    
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        stringtoint[2]["Surface"] = 1;
        stringtoint[1]["Bottom"] = 2;
        stringtoint[1]["Right"] = 3;
        stringtoint[1]["Top"] = 4;
        stringtoint[1]["Left"] = 5;
        stringtoint[0]["Point"] = 6;
        stringtoint[1]["Top2"] = 7;
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh4("../mesh/1element.msh",gmesh);
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    
    //.................................Hdiv.................................
    TPZCompMesh * cmeshflux = 0;
    TPZCompMesh * cmeshpressure = 0;
    // {
    //     // para a primeira simulacao a malha de fluxo contem todos os materiais
    //     // nao ha material EWrapBC nesta configuracao
    //     TPZManVector<int64_t, 6> matIdVecHdiv={EDomain,EBottom,ETop,ELeft,EPont,ERight};
    //     TPZManVector<int64_t, 2> matIdNeumannHdiv;

    //     //Flux mesh
    //     cmeshflux = FluxCMesh(dim,pOrder,matIdVecHdiv,gmesh);

    //     //Pressure mesh
    //     // a malha pressao sera vazio neste caso
    //     // a dimensao dos elementos de pressao eh um a menos que a dimensao do problema
    //     cmeshpressure = PressureCMesh(dim-1,pOrder,matIdVecHdiv,gmesh);

    //     //Multiphysics mesh
    //     TPZManVector< TPZCompMesh *, 2> meshvector(2,0);
    //     meshvector[0] = cmeshflux;
    //     meshvector[1] = cmeshpressure;
    //     auto * cmesh = MultiphysicCMesh(dim,pOrder,matIdVecHdiv,meshvector,gmesh);

    //     //Solve Multiphysics
    //     TPZLinearAnalysis an(cmesh,true);
    //     SolveProblemDirect(an,cmesh);

    //     // //Print results
    //     // PrintResultsMultiphysic(dim,meshvector,an,cmesh);
    //     // std::ofstream out3("mesh_Hdiv.txt");
    //     // an.Print("nothing",out3);

    //     std::ofstream anPostProcessFileHdiv("postprocessHdiv.txt");
    //     ComputeErrorHdiv(an,anPostProcessFileHdiv);

    // }

    //.........................Div Free Bubbles NEW.........................
    // {
        std::set<int> matBC={ERight};
        HybridizerGeoMesh(gmesh,matBC);

        // nesta configuracao o material ETop eh substituido por EWrap
        std::set<int> matIdVecNew={EDomain,ETop,EBottom,ELeft,EWrap,EPont};
        // criamos elementos tipo ETop para a pressao
        std::set<int> matIdNeumannNew = matBC;
        matIdNeumannNew.insert(EPressureHyb);
        
        //Flux mesh
        TPZCompMesh * cmeshfluxNew = FluxCMeshNew(dim,pOrder+1,matIdVecNew,gmesh);        

        //Pressure mesh
        TPZCompMesh * cmeshpressureNew = PressureCMeshNew(dim,pOrder,matIdNeumannNew,gmesh);

        //Multiphysics mesh
        TPZManVector< TPZCompMesh *, 2> meshvectorNew(2);
        meshvectorNew[0] = cmeshfluxNew;
        meshvectorNew[1] = cmeshpressureNew;

        auto * cmeshNew = MultiphysicCMeshNew(dim,pOrder+1,matIdVecNew,matIdNeumannNew,meshvectorNew,gmesh);

        //Solve Multiphysics
        TPZLinearAnalysis anNew(cmeshNew,false);

        SolveProblemDirect(anNew,cmeshNew);

        // anNew.Solution().Print("Sol");
        
        //Print results
        PrintResultsMultiphysicNew(dim,meshvectorNew,anNew,cmeshNew);
        std::ofstream out4("mesh_MDFB.txt");
        anNew.Print("nothing",out4);

        std::ofstream anPostProcessFileMDFB("postprocessMDFB.txt");
        ComputeErrorHdiv(anNew,anPostProcessFileMDFB);
    // }

    //...........................Div Free Bubbles...........................
    //Creates DFB problem
    // TPZCompMesh * cmeshDFB = CMeshDivFreeBubbles(dim,pOrder,matIdVec,gmesh);

    // //Solve DFB
    // TPZLinearAnalysis anDFB(cmeshDFB,true);
    // SolveProblemDirect(anDFB,cmeshDFB);

    // // //Print results
    // PrintResultsDivFreeBubbles(dim,anDFB);
    // std::ofstream out("mesh.txt");
    // anDFB.Print("nothing",out);

    return 0;
}



//Flux computational mesh
TPZCompMesh *FluxCMesh(int dim, int pOrder,std::set<int> matIdVec, TPZGeoMesh *gmesh) 
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        mat->SetDimension(dim);
        mat->SetBigNumber(1.e10);
        cmesh->InsertMaterialObject(mat);
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
TPZCompMesh *PressureCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh)
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    if(matIdVec.size() == 0) return cmesh;

    TPZNullMaterial<> *mat = new TPZNullMaterial<>(EDomain);
    mat->SetDimension(dim);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
    // distincao de ordem zero
    if(pOrder == 0)
    {
        // os elementos H1 nao tem opcao de funcao constante
        cmesh->SetAllCreateFunctionsDiscontinuous();
    }
    else
    {
        // Disconnected = true faz com que os espaços sao descontinuos
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
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
        if(!celdisc) continue;
        celdisc->SetConstC(1.);
        celdisc->SetTrueUseQsiEta();
        // espera-se elemento de pressao apenas para o contorno
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            DebugStop();
        }
    }
    // Print pressure mesh
    std::ofstream myfile("PressureMesh.txt");
    cmesh->Print(myfile);

    return cmesh;
}

// Multiphysics computational mesh
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, std::set<int> matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    auto mat = new TPZMixedDarcyFlow(EDomain, dim);

    mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
        
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,0.);
    
    auto * BCond0 = mat->CreateBC(mat, EBottom, 0, val1, val4);//Bottom
    auto * BCond1 = mat->CreateBC(mat, ERight, 0, val1, val2);//Right
    BCond0->SetForcingFunctionBC(exactSol);
    BCond1->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond0);
    auto * BCond2 = mat->CreateBC(mat, ETop, 0, val1, val2);//Top
    auto * BCond3 = mat->CreateBC(mat, ELeft, 0, val1, val4);//Left
    BCond2->SetForcingFunctionBC(exactSol);
    BCond3->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    // TPZBndCond * BCond4 = mat->CreateBC(mat, EPont, 0, val1, val4);//Left
    // cmesh->InsertMaterialObject(BCond4);
    // auto *mat2 = new TPZMatPoisson<>(EPont,dim);
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

    //Prints computational mesh properties
    std::stringstream vtk_name;
    vtk_name    << "MultiPhysics" << ".vtk";
    std::ofstream vtkfile(vtk_name.str().c_str());
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

    return cmesh;
}

void CreateMultiphysicsInterfaceElements(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvector, std::set<int> &matIdNeumann){

    cmesh->LoadReferences();
    
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->MaterialId() != EWrap) continue;
        auto nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        // here I generalized - an interface is created whenever a wrap element exists
        auto gelsidepr = gelside.HasNeighbour(matIdNeumann);
        if (!gelsidepr)
        {
            DebugStop();
        }

        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide celneigh = gelsidepr.Reference();
        if (!celside || !celneigh) {
            DebugStop();
        }

        int64_t index;
        TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(*cmesh,gel->Neighbour(2).Element(),index,celneigh,celside); // E AQUI O COMPUTACIONAL
    }

}

// Multiphysics computational mesh
TPZMultiphysicsCompMesh *MultiphysicCMeshNew(int dim, int pOrder, std::set<int> &matIdVec, std::set<int> &matIdNeumann, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
    // gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    // eh preciso criar materiais para todos os valores referenciados no enum
    auto mat = new TPZMixedDarcyFlow(EDomain, dim);
    mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
        
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,1.);
    constexpr int boundType{1};
    constexpr int boundType0{0};
    auto * BCond0 = mat->CreateBC(mat, EBottom, 0, val1, val2);
    auto * BCond1 = mat->CreateBC(mat, ERight, 1, val1, val2);
    BCond0->SetForcingFunctionBC(exactSol);
    BCond1->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond0);
    auto * BCond2 = mat->CreateBC(mat, ETop, 0, val1, val2);
    auto * BCond3 = mat->CreateBC(mat, ELeft, 0, val1, val2);
    BCond2->SetForcingFunctionBC(exactSol);
    BCond3->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);

    auto *matL2 = new TPZL2ProjectionCS<>(EPont,0,1);
    cmesh->InsertMaterialObject(matL2);

    auto * nullmat2 = new TPZNullMaterialCS<>(EWrap,1,1);
    cmesh->InsertMaterialObject(nullmat2);

    auto * nullmat3 = new TPZNullMaterialCS<>(EPressureHyb,1,1);
    cmesh->InsertMaterialObject(nullmat3);

    TPZManVector<int> active(2,1);
    active[0]=1;
    active[1]=1;
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes(); 

    auto mat3 = new TPZLagrangeMultiplierCS<STATE>(EIntfaceLeft, dim-1);
    cmesh->InsertMaterialObject(mat3);

    auto mat4 = new TPZLagrangeMultiplierCS<STATE>(EIntfaceRight, dim-1);
    cmesh->InsertMaterialObject(mat4);
    
    CreateMultiphysicsInterfaceElements(cmesh,gmesh,meshvector,matIdNeumann);
    

    // std::cout << "Nequations before = " << cmesh->NEquations() << std::endl;
    // for (auto cel:cmesh->ElementVec()) {
    //     if (!cel && cel->Dimension() != cmesh->Dimension()) continue;
    //     TPZCondensedCompEl *condense = new TPZCondensedCompEl(cel, false);
    // }
    // cmesh->CleanUpUnconnectedNodes();
    // std::cout << "Nequations after = " << cmesh->NEquations() << std::endl;

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
    constexpr int nThreads{0};
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

void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh)
{

    an.SetExact(exactSol,solOrder);

    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
    TPZManVector<std::string,10> scalnames(0), vecnames(2);

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

void PrintResultsMultiphysicNew(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh)
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


TPZCompMesh *CMeshDivFreeBubbles(int dim, int pOrder, std::set<int> matIdVec, TPZGeoMesh *gmesh)
{
    gmesh->ResetReference();
        
    //Creates cmesh object
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDefaultOrder(pOrder);

    //Sets materials
    auto *mat = new TPZMatDivFreeBubbles<STATE>(EDomain,dim);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);

    //Insert boundary conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,0.);
    constexpr int boundType{0};
    auto * BCond = mat->CreateBC(mat,EBottom, 0, val1, val4);//Bottom
    // auto * BCond1 = mat->CreateBC(mat, ERight, 0, val1, val2);//Right
    BCond->SetForcingFunctionBC(exactSol);
    // BCond1->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond);
    // cmesh->InsertMaterialObject(BCond1);
    auto * BCond2 = mat->CreateBC(mat, ETop, 1, val1, val2);//Top
    auto * BCond3 = mat->CreateBC(mat, ELeft, 0, val1, val4);//Left
    BCond2->SetForcingFunctionBC(exactSol);
    BCond3->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);

    auto *mat2 = new TPZL2Projection<>(EPont,0,1);
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
TPZCompMesh *FluxCMeshNew(int dim,int pOrder,std::set<int> matIdVec,TPZGeoMesh *gmesh) 
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    std::set<int> allmat;  

    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        allmat.insert(*it);
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        cmesh->InsertMaterialObject(mat);
        mat->SetDimension(dim);
        mat->SetBigNumber(1.e10);
    }

    // allmat.insert(EPont);
    // TPZL2Projection<> *mat = new TPZL2Projection<>(EPont,0,1);
    // cmesh->InsertMaterialObject(mat);
    // mat->SetDimension(dim);
    // mat->SetBigNumber(1.e10);

    //Disconnect elements
    // bool flag = false;
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) DebugStop();
        auto type = gel -> Type();
        int64_t index;
        auto matid = gel->MaterialId();
        if(allmat.find(matid) == allmat.end()) continue;
        // if (gel->Dimension() == 2) flag=true;
        // if (flag==false) continue;

        using namespace pzgeom;
        using namespace pzshape;

        if (type == EPoint){
            // new TPZIntelGen<TPZShapePoint>(*cmesh,gel,index);
            // TPZMaterial *mat = cmesh->FindMaterial(matid);
            // TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            // nullmat->SetDimension(0);
        } else if (type == EOned){
            // if (matid == ERight){
            //     new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,gel,index);
            //     TPZMaterial *mat = cmesh->FindMaterial(matid);
            //     TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            //     nullmat->SetDimension(1);
            // }
        } else if (type == EQuadrilateral){
            gel->ResetReference();   
            new TPZCompElKernelHDiv<TPZShapeQuad>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(2);

            for (int side=0; side<gel->NSides(); side++){
                
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();

                if (neighbour.Element()->MaterialId() == EPont){
                    new TPZIntelGen<TPZShapePoint>(*cmesh,neighbour.Element(),index);
                    TPZMaterial *mat = cmesh->FindMaterial(EPont);
                    TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                }
                
                if (neighbour.Element()->Dimension() != dim-1) continue;

                if (neighbour.Element()->MaterialId() == EWrap){
                    new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,neighbour.Element(),index);
                    TPZMaterial *mat = cmesh->FindMaterial(EWrap);
                    TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                    nullmat->SetDimension(1);
                    neighbour.Element()->ResetReference();
                } else {
                    if (allmat.find(neighbour.Element()->MaterialId()) == allmat.end()){

                    } else {
                        new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,neighbour.Element(),index);
                        TPZMaterial *mat = cmesh->FindMaterial(neighbour.Element()->MaterialId());
                        TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                        nullmat->SetDimension(1);
                        neighbour.Element()->ResetReference();
                    }
                }
                
            }
            
            for (int side = 0; side < 8; side++)
            {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                neighbour.Element()->ResetReference();
            }
            gel->ResetReference();    
        } else if(type == ETriangle) {
            new TPZCompElKernelHDiv<TPZShapeTriang>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(2);
        }
    }    

    //  cmesh->Reference()->ResetReference();
    // cmesh->LoadReferences();
    // cmesh->CleanUpUnconnectedNodes();

    // cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->ExpandSolution();

    // cmesh->AutoBuild();
    cmesh->InitializeBlock();
    cmesh->ComputeNodElCon();

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
TPZCompMesh *PressureCMeshNew(int dim, int pOrder, std::set<int> matIdVec, TPZGeoMesh *gmesh)
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    std::set<int> matIdNeumann;

    // Sets matid to BC geometric elements
    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        mat->SetDimension(1);
        cmesh->InsertMaterialObject(mat);
        mat->SetBigNumber(1.e10);
        matIdNeumann.insert(*it); 
    }

    
    // create internal geometric elements to the hybridization
    // int64_t nel = gmesh->NElements();
    // for (int64_t el = 0; el < nel; el++) {
    //     TPZGeoEl *gel = gmesh->Element(el);
    //     if(!gel) DebugStop();
    //     int matid = gel->MaterialId();
    //     // loop only over volumetric elements
    //     if(matid != EDomain) continue;
    //     if (gel->Dimension() != dim) {
    //         DebugStop();
    //     }
    //     int nsides = gel->NSides();
    //     int ncorner = gel->NCornerNodes();
    //     for (int side = ncorner; side < nsides; side++) {
    //         if(gel->SideDimension(side) != dim-1) continue;
    //         TPZGeoElSide gelside(gel,side);
    //         TPZGeoElSide neighbour = gelside.Neighbour();
            
    //         if (neighbour.Element()->MaterialId() != EWrap) continue;
    //         TPZCompElSide celside = gelside.Reference();
    //         int dir = gel->NormalOrientation(celside);
    //         if (side == 6) {
    //             TPZGeoElBC gbc(gelside, EPressureHyb);
    //         }
    //         // gel->Reference()->LoadElementReference();

    //     }
    // }

    // TPZNullMaterial<> *mat1 = new TPZNullMaterial<>(EPressureHyb);
    // mat1->SetDimension(1);
    // mat1->SetBigNumber(1.e10);
    // cmesh->InsertMaterialObject(mat1);
    // matIdNeumann.insert(EPressureHyb);

    // TPZNullMaterial<> *mat2 = new TPZNullMaterial<>(EIntfaceRight);
    // mat2->SetDimension(1);
    // mat2->SetBigNumber(1.e10);
    // cmesh->InsertMaterialObject(mat2);
    // matIdNeumann.insert(EIntfaceRight);
    
    // TPZNullMaterial<> *mat3 = new TPZNullMaterial<>(EWrap);
    // mat3->SetDimension(1);
    // mat3->SetBigNumber(1.e10);
    // cmesh->InsertMaterialObject(mat3);
    // matIdNeumann.insert(EWrap);

    cmesh->InitializeBlock();

    if(pOrder > 0)
    {
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
        cmesh->SetDefaultOrder(pOrder);
    } else
    {
        cmesh->SetDefaultOrder(pOrder);
        cmesh->SetDimModel(1);
        cmesh->SetAllCreateFunctionsDiscontinuous();
    }
    cmesh->AutoBuild(matIdNeumann);

    for(auto &newnod : cmesh->ConnectVec())
    {
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


void HybridizerGeoMesh(TPZGeoMesh * gmesh, std::set<int> &matIdBC){

    //Disconnect elements
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) DebugStop();
        auto type = gel -> Type();
        int64_t index;
        auto matid = gel->MaterialId();

        using namespace pzgeom;
        using namespace pzshape;

        if (gel->Dimension() != 2) continue;
        
        int nsides = gel->NSides();
        // gel->ResetReference();
        for (int side = 0; side < nsides; side++) {

            if(gel->SideDimension(side) != 1) continue;
            TPZGeoElSide geoside(gel,side);
            TPZGeoElSide neighbour = geoside.Neighbour();

            if (neighbour.Element()->MaterialId() == EDomain){
                TPZGeoElBC gelbcWrap(geoside, EWrap);
           
                // int dir = gel->NormalOrientation(side);
                
                TPZGeoElSide gelWrapSide(gelbcWrap.CreatedElement(),2);
                // CRIA OS ELEMENTOS GEOMÉTRICOS WRAP E INTERFACE
                if (dir == 1){
                    TPZGeoElBC gelbc(gelWrapSide, EIntfaceLeft); // AQUI CRIA O ELEMENTO DE INTERFACE GEOMÉTRICO
                    
                    bool flag = false;
                    TPZGeoElSide sidePresHyb;
                    for (int k = 4; k < 8; k++)
                    {
                        TPZGeoElSide geosideNeig(neighbour.Element(),k);
                        if (geosideNeig.Element()->Neighbour(k).Element()->MaterialId()==EWrap) {
                            flag = true;
                            sidePresHyb = geosideNeig.Element()->Neighbour(2).Element()->Neighbour(2).Element()->Neighbour(2);
                        }
                    }
                    if (flag==false){
                        TPZGeoElSide gelPresHSide(gelbc.CreatedElement(),2);
                        TPZGeoElBC gelPHyb(gelPresHSide, EPressureHyb);
                    }
                    
                }else if (dir == -1) {
                    TPZGeoElBC gelbc(gelWrapSide, EIntfaceRight); // AQUI CRIA O ELEMENTO DE INTERFACE GEOMÉTRICO
                    
                    bool flag = false;
                    TPZGeoElSide sidePresHyb;
                    for (int k = 4; k < 8; k++)
                    {
                        TPZGeoElSide geosideNeig(neighbour.Element(),k);
                        if (geosideNeig.Element()->Neighbour(k).Element()->MaterialId()==EWrap) {
                            flag = true;
                            sidePresHyb = geosideNeig.Element()->Neighbour(2).Element()->Neighbour(2).Element()->Neighbour(2);
                        }
                    }
                    if (flag==false){
                        TPZGeoElSide gelPresHSide(gelbc.CreatedElement(),2);
                        TPZGeoElBC gelPHyb(gelPresHSide, EPressureHyb);
                    }
                }
            } 
            if (matIdBC.find(neighbour.Element()->MaterialId()) != matIdBC.end()){
                TPZGeoElBC gelbcWrap(geoside, EWrap);
                TPZGeoElSide gelWrapSide(gelbcWrap.CreatedElement(),2);
                TPZGeoElBC gelbc(gelWrapSide, EIntfaceLeft);
            }
        }

        //Creates point
        TPZGeoElSide geosidePoint(gel,0);
        TPZGeoElBC gelbcPoint(geosidePoint, EPont);

        gel->ResetReference();
        for (int side = 0; side < nsides; side++) {
            TPZGeoElSide gelside(gel,side);
            gelside.Element()->ResetReference();
            TPZGeoElSide neighbour = gelside.Neighbour();
            neighbour.Element()->ResetReference();
        }   
    
    }

}