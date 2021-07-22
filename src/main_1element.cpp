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
TPZMultiphysicsCompMesh *MultiphysicCMeshNew(int dim, int pOrder, std::set<int> &matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
TPZCompMesh *CMeshDivFreeBubbles(int dim, int pOrder, std::set<int> matIdVec, TPZGeoMesh *gmesh);
void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh);
void PrintResultsMultiphysicNew(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh);
void PrintResultsDivFreeBubbles(int dim, TPZLinearAnalysis &an);
void ComputeError(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);
void ComputeErrorHdiv(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);
void HybridizeBC(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh);
void HybridizeDomain(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh);
TPZCompEl * FindFluxElement(TPZCompEl *wrapelement);
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

enum EMatid  {ENone, EDomain, EBottom, ERight, ETop, ELeft, EPont, EWrapBC, EIntfaceBC, EWrap, EIntfaceLeft, EIntfaceRight};

int main(int argc, char* argv[])
{
    //dimension of the problem
    constexpr int dim{2};
    constexpr int pOrder{1};
      

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
        // nesta configuracao o material ETop eh substituido por EWrapBC
        std::set<int> matIdVecNew={EDomain,EPont,EWrapBC};
        // criamos elementos tipo ETop para a pressao
        std::set<int> matIdNeumannNew = {ETop,ERight,EBottom,ELeft};
        
        //Flux mesh
        TPZCompMesh * cmeshfluxNew = FluxCMeshNew(dim,pOrder+1,matIdVecNew,gmesh);

        //Pressure mesh
        TPZCompMesh * cmeshpressureNew = PressureCMeshNew(dim,pOrder,matIdNeumannNew,gmesh);

        //Multiphysics mesh
        TPZManVector< TPZCompMesh *, 2> meshvectorNew(2);
        meshvectorNew[0] = cmeshfluxNew;
        meshvectorNew[1] = cmeshpressureNew;

        // TPZHybridizeHDiv hybridizer(meshvectorNew);


        auto * cmeshNew = MultiphysicCMeshNew(dim,pOrder+1,matIdVecNew,meshvectorNew,gmesh);

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

void HybridizeBC(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh * gmesh){

    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->MaterialId() != EWrapBC)
        {
            continue;
        }
        auto nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        // here I generalized - an interface is created whenever a wrap element exists
        std::set<int> neighmat = {ELeft,EBottom,ERight,ETop};
        auto gelsidepr = gelside.HasNeighbour(neighmat);
        if (!gelsidepr)
        {
            DebugStop();
        }

        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide celneigh = gelsidepr.Reference();
        if (!celside || !celneigh) {
            DebugStop();
        }
        TPZGeoElBC gelbc(gelside, EIntfaceBC); // AQUI CRIA O ELEMENTO DE INTERFACE GEOMÉTRICO
        int64_t index;
        TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(*cmesh,gelbc.CreatedElement(),index,celneigh,celside); // E AQUI O COMPUTACIONAL
    }
}

void HybridizeDomain(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh * gmesh){

    gmesh->ResetReference();
    cmesh->LoadReferences();
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();

        if(matid == EWrap)
        {
            TPZCompEl *fluxel = FindFluxElement(cel);
            TPZGeoEl *fluxgel = fluxel->Reference();
            TPZGeoElSide gelside(gel);
            TPZGeoElSide neighbour = gelside.Neighbour();
            int neighmat = neighbour.Element()->MaterialId();
            if(neighmat != EIntfaceLeft && neighmat != EIntfaceRight)
            {
                DebugStop();
            }
            // determine if the interface should be positive or negative...
            int interfacematid = neighmat;
            int64_t index;
            TPZCompElSide celwrap(cel,gel->NSides()-1);
            TPZGeoElSide fluxgelside(fluxgel);
            TPZCompElSide fluxside = fluxgelside.Reference();
//            std::cout << "Creating interface from wrap element " << gel->Index() << " using neighbour " << neighbour.Element()->Index() <<
//             " and flux element " << fluxgel->Index() << std::endl;
            if(neighbour.Element()->Reference()) DebugStop();
            new TPZMultiphysicsInterfaceElement(*cmesh,neighbour.Element(),index,celwrap,fluxside);

        }
    }
}

// Multiphysics computational mesh
TPZMultiphysicsCompMesh *MultiphysicCMeshNew(int dim, int pOrder, std::set<int> &matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
    
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();


    // eh preciso criar materiais para todos os valores referenciados no enum
    auto mat = new TPZMixedDarcyFlow(EDomain, dim);
    // auto * mat = new TPZNullMaterialCS<>(EDomain,1,1);
    mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
        
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,0.);
    constexpr int boundType{1};
    constexpr int boundType0{0};
    auto * BCond0 = mat->CreateBC(mat, EBottom, 1, val1, val4);
    auto * BCond1 = mat->CreateBC(mat, ERight, 1, val1, val2);
    BCond0->SetForcingFunctionBC(exactSol);
    BCond1->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond0);
    auto * BCond2 = mat->CreateBC(mat, ETop, 1, val1, val2);
    auto * BCond3 = mat->CreateBC(mat, ELeft, 1, val1, val4);
    BCond2->SetForcingFunctionBC(exactSol);
    BCond3->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);


    // TPZMaterial *mat8 = cmesh->FindMaterial(EDomain);
    // // if(mat8) return;
    // auto *nullmat3 = new TPZNullMaterialCS<>(EDomain);
    // nullmat3->SetDimension(dim-1);
    // nullmat3->SetNStateVariables(1);
    // cmesh->InsertMaterialObject(nullmat3);



    // the wrap material is a null material (does nothing)
    auto * nullmat = new TPZNullMaterialCS<>(EWrapBC,1,1);
    cmesh->InsertMaterialObject(nullmat);

    auto * nullmat2 = new TPZNullMaterialCS<>(EWrap,1,1);
    cmesh->InsertMaterialObject(nullmat2);

    auto mat2 = new TPZLagrangeMultiplierCS<STATE>(EIntfaceBC, dim);
    cmesh->InsertMaterialObject(mat2);

    auto mat3 = new TPZLagrangeMultiplierCS<STATE>(EIntfaceLeft, dim);
    cmesh->InsertMaterialObject(mat3);

    auto mat4 = new TPZLagrangeMultiplierCS<STATE>(EIntfaceRight, dim);
    cmesh->InsertMaterialObject(mat4);
    
    TPZManVector<int> active(2,1);
    active[0]=1;
    active[1]=1;
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes(); 

    HybridizeBC(cmesh,gmesh);
    // HybridizeDomain(cmesh,gmesh);

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

    for (int64_t i = 0; i < gmesh->NElements(); i++)
    {
        auto *gel = gmesh -> Element(i);
        auto type = gel -> Type();
        int64_t index;
        auto matid = gel->MaterialId();
        if(allmat.find(matid) == allmat.end()) continue;
       
        using namespace pzgeom;
        using namespace pzshape;
        if (type == EPoint){
            new TPZIntelGen<TPZShapePoint>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(0);
        } else if (type == EOned){
            new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(1);
        } else if (type == EQuadrilateral){
            new TPZCompElKernelHDiv<TPZShapeQuad>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(2);
        } else if(type == ETriangle) {
            new TPZCompElKernelHDiv<TPZShapeTriang>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(2);
        }
    }
    cmesh->ExpandSolution();

    // cmesh->AutoBuild();

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
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) DebugStop();
        int matid = gel->MaterialId();
        // loop only over volumetric elements
        if(matid != EDomain) continue;
        if (gel->Dimension() != dim) {
            DebugStop();
        }
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        for (int side = ncorner; side < nsides; side++) {
            if(gel->SideDimension(side) != dim-1) continue;
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            
            if (neighbour.Element()->MaterialId() != EDomain) continue;
            int dir = gel->NormalOrientation(side);
            
            // CRIA OS ELEMENTOS GEOMÉTRICOS WRAP E INTERFACE
            if (dir == 1){
                TPZGeoElBC gelbcInterface(gelside, EIntfaceLeft);
                TPZGeoElBC gelbcWrap(gelside, EWrap);
            }else if (dir == -1) {
                TPZGeoElBC gelbcInterface(gelside, EIntfaceRight); 
                TPZGeoElBC gelbcWrap(gelside, EWrap); 
            }
        }
    }

    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) DebugStop();
        int matid = gel->MaterialId();
        // loop only over volumetric elements
        if(matid != EDomain) continue;
        if (gel->Dimension() != dim) {
            DebugStop();
        }
        int nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        for (int side = ncorner; side < nsides; side++) {
            if(gel->SideDimension(side) != dim-1) continue;
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            
            std::cout << "Element = " << el << ", side = " << side 
                      << ", neighbour matid = " << neighbour.Element()->MaterialId()
                      << ", neighbour's neighbour matid = " << neighbour.Neighbour().Element() -> MaterialId() << std::endl;
        }
    }

    TPZNullMaterial<> *mat1 = new TPZNullMaterial<>(EIntfaceLeft);
    mat1->SetDimension(1);
    mat1->SetBigNumber(1.e10);
    cmesh->InsertMaterialObject(mat1);
    matIdNeumann.insert(EIntfaceLeft);

    TPZNullMaterial<> *mat2 = new TPZNullMaterial<>(EIntfaceRight);
    mat2->SetDimension(1);
    mat2->SetBigNumber(1.e10);
    cmesh->InsertMaterialObject(mat2);
    matIdNeumann.insert(EIntfaceRight);
    
    TPZNullMaterial<> *mat3 = new TPZNullMaterial<>(EWrap);
    mat3->SetDimension(1);
    mat3->SetBigNumber(1.e10);
    cmesh->InsertMaterialObject(mat3);
    matIdNeumann.insert(EWrap);

    







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

TPZCompEl * FindFluxElement(TPZCompEl *wrapelement)
{
    TPZGeoEl *gel = wrapelement->Reference();
    int nsides = gel->NSides();
    TPZGeoElSide gelside(gel,nsides-1);
    TPZStack<TPZCompElSide> celstack;
    gelside.EqualLevelCompElementList(celstack, 0, 0);
    int nelstack = celstack.size();
    for (int st = 0; st<nelstack; st++) {
        TPZCompElSide celside = celstack[st];
        TPZCompEl *cel = celside.Element();
        TPZGeoEl *gelneigh = cel->Reference();
        int matid = gelneigh->MaterialId();
        if (matid == EDomain) {
            return cel;
        }
    }
    TPZCompElSide cellarge = gelside.LowerLevelCompElementList2(0);
    if(!cellarge){
        gelside.Print(std::cout);
        cellarge = gelside.LowerLevelCompElementList2(0);
        DebugStop();
    }
    {
        TPZStack<TPZCompElSide> celstack;
        TPZGeoElSide gellarge = cellarge.Reference();
        if(gellarge.Element()->MaterialId() == EDomain)
        {
            return cellarge.Element();
        }
        gellarge.EqualLevelCompElementList(celstack, 0, 0);
        int nelst = celstack.size();
        for (int ist = 0; ist < nelst; ist++) {
            TPZCompElSide celside = celstack[ist];
            TPZGeoElSide gelside = celside.Reference();
            TPZGeoEl *candidate = gelside.Element();
            if(candidate->MaterialId() == EDomain)
            {
                return celside.Element();
            }
        }
    }
    return NULL;
}
