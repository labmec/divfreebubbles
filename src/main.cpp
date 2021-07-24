#ifdef HAVE_CONFIG_H
  #include <pz_config.h>
#endif

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
#include <TPZLagrangeMultiplierCS.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"// for Hdiv problem
#include <TPZBndCond.h> //for TPZBndCond
#include "TPZLinearAnalysis.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "TPZMultiphysicsCompMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzstrmatrixor.h"
#include "../headers/TPZCompElKernelHdiv.h" //THE NEW MATERIAL!
#include "../headers/TPZCompElKernelHdivBC.h" //THE NEW MATERIAL!
#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"
// #include "pzlog.h"

TPZCompMesh *FluxCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *FluxCMeshDFB(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMeshDFB(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
TPZMultiphysicsCompMesh *MultiphysicCMeshDFB(int dim, int pOrder, std::set<int> &matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
TPZCompMesh *CMeshH1(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *CMeshDivFreeBubbles(int dim, int pOrder, std::set<int> matIdVec, TPZGeoMesh *gmesh);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysicDFB(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsH1(int dim, TPZLinearAnalysis &an);
void PrintResultsDivFreeBubbles(int dim, TPZLinearAnalysis &an);
void ComputeError(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);
void ComputeErrorHdiv(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);
void HybridizeBC(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh * gmesh);
void HybridizeDomain(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh * gmesh);

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

enum EMatid {ENone, EDomain, EInjection, EProduction, EBottom,  ETop, ELeft, ERight, EPont, EWrapBC, EIntfaceBC, EWrap, EIntfaceLeft, EIntfaceRight};

int main(int argc, char* argv[]){
    //dimension of the problem
    constexpr int dim{2};
    constexpr int pOrder{2};

    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        stringtoint[2]["Surface"] = 1;
        stringtoint[1]["InjectionWell"] = 2;
        stringtoint[1]["ProductionWell"] = 3;
        stringtoint[1]["BottomLine"] = 4;
        stringtoint[1]["TopLine"] = 5;
        stringtoint[1]["LeftLine"] = 6;
        stringtoint[1]["RightLine"] = 7;
        stringtoint[0]["Point"] = 8;
        stringtoint[1]["BottomLine2"] = 9;
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh4("../mesh/newMesh.msh",gmesh);
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
  
    //.................................Hdiv.................................
    TPZCompMesh * cmeshflux = 0;
    TPZCompMesh * cmeshpressure = 0;
    // {
    //     std::set<int> matIdVecHdiv={EDomain,EInjection,EProduction,EBottom,ETop,ELeft,ERight};
    //     std::set<int> matIdNeumannHdiv;
        
    //     //Flux mesh
    //     cmeshflux = FluxCMesh(dim,pOrder,matIdVecHdiv,gmesh);

    //     //Pressure mesh
    //     cmeshpressure = PressureCMesh(dim,pOrder,matIdVecHdiv,gmesh);

    //     //Multiphysics mesh
    //     TPZManVector< TPZCompMesh *, 2> meshvector(2);
    //     meshvector[0] = cmeshflux;
    //     meshvector[1] = cmeshpressure;
    //     TPZCompMesh * cmesh = MultiphysicCMesh(dim,pOrder,matIdVecHdiv,meshvector,gmesh);

    //     //Solve Multiphysics
    //     TPZLinearAnalysis an(cmesh,true);
    //     SolveProblemDirect(an,cmesh);

    //     //Print results
    //     PrintResultsMultiphysic(dim,meshvector,an,cmesh);

    //     std::ofstream anPostProcessFileHdiv("postprocessHdiv.txt");
    //     ComputeErrorHdiv(an,anPostProcessFileHdiv);
    // }

    //..................................H1..................................
    TPZCompMesh * cmeshH1 = 0;
    // {
    //     std::set<int> matIdVecH1={EDomain,EInjection,EProduction,EBottom,ETop,ELeft,ERight};

    //     //Creates H1 problem
    //     cmeshH1 = CMeshH1(dim,pOrder,matIdVecH1,gmesh);

    //     //Solve H1 
    //     TPZLinearAnalysis anH1(cmeshH1,false);
    //     SolveProblemDirect(anH1,cmeshH1);

    //     //Print results
    //     PrintResultsH1(dim,anH1);

    //     std::ofstream anPostProcessFileH1("postprocessH1.txt");
    //     ComputeError(anH1,anPostProcessFileH1);

    // }

    //...........................Div Free Bubbles...........................
    TPZCompMesh * cmeshDFB = 0;
    // {
    //     std::set<int> matIdVecDFB={EDomain,EInjection,EProduction,EBottom,ETop,ELeft,ERight};

    //     //Creates DFB problem
    //     cmeshDFB = CMeshDivFreeBubbles(dim,pOrder+1,matIdVecDFB,gmesh);

    //     //Solve DFB
    //     TPZLinearAnalysis anDFB(cmeshDFB,false);
    //     SolveProblemDirect(anDFB,cmeshDFB);

    //     //Print results
    //     PrintResultsDivFreeBubbles(dim,anDFB);

    //     std::ofstream anPostProcessFileDFB("postprocessDFB.txt");
    //     ComputeError(anDFB,anPostProcessFileDFB);
    // }
  
    //.....................Div Free Bubbles Multiphysic.....................
    {
        // nesta configuracao o material ETop eh substituido por EWrap
        std::set<int> matIdVecNew={EDomain,EPont,EWrapBC};
        // criamos elementos tipo ETop para a pressao
        std::set<int> matIdNeumannNew = {EInjection,EProduction,EBottom,ETop,ELeft,ERight};

        //Flux mesh
        TPZCompMesh * cmeshfluxDFB = FluxCMeshDFB(dim,pOrder+1,matIdVecNew,gmesh);

        //Pressure mesh
        TPZCompMesh * cmeshpressureDFB = PressureCMeshDFB(dim,pOrder,matIdNeumannNew,gmesh);

        //Multiphysics mesh
        TPZManVector< TPZCompMesh *, 2> meshvectorDFB(2);
        meshvectorDFB[0] = cmeshfluxDFB;
        meshvectorDFB[1] = cmeshpressureDFB;
        auto * cmeshNew = MultiphysicCMeshDFB(dim,pOrder+1,matIdVecNew,meshvectorDFB,gmesh);
        
        //Solve Multiphysics
        TPZLinearAnalysis anNew(cmeshNew,false);
        SolveProblemDirect(anNew,cmeshNew);

        //Print results
        PrintResultsMultiphysicDFB(dim,meshvectorDFB,anNew,cmeshNew);
        std::ofstream out4("mesh_MDFB.txt");
        anNew.Print("nothing",out4);
        
        std::ofstream anPostProcessFileMDFB("postprocessMDFB.txt");
        ComputeErrorHdiv(anNew,anPostProcessFileMDFB);
    }

  return 0;
}

//Analytical solution
constexpr int solOrder{2};
auto exactSol = [](const TPZVec<REAL> &loc,
  TPZVec<STATE>&u,
  TPZFMatrix<STATE>&gradU){
  const auto &x=loc[0];
  const auto &y=loc[1];
  const auto &d = 1.; // distance between injection and production wells
  u[0]= log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
  gradU(0,0) = (x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2)));
  gradU(1,0) = (y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2)));
};
auto exactSolH1 = [](const TPZVec<REAL> &loc,
  TPZVec<STATE>&u,
  TPZFMatrix<STATE>&gradU){
  const auto &x=loc[0];
  const auto &y=loc[1];
  const auto &d = 1.; // distance between injection and production wells
  u[0]= log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
//   gradU(0,0) = (x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2)));
//   gradU(1,0) = (y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2)));
};
auto exactSolError = [](const TPZVec<REAL> &loc,
  TPZVec<STATE>&u,
  TPZFMatrix<STATE>&gradU){
  const auto &x=loc[0];
  const auto &y=loc[1];
  const auto &d = 1.; // distance between injection and production wells
  u[0]= log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
  gradU(0,0) = -(x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2)));
  gradU(1,0) = -(y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2)));
};

//Flux computational mesh
TPZCompMesh *FluxCMesh(int dim, int pOrder,std::set<int> &matIdVec, TPZGeoMesh *gmesh) 
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

//Flux computational mesh
TPZCompMesh *FluxCMeshDFB(int dim, int pOrder,std::set<int> &matIdVec, TPZGeoMesh *gmesh) 
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
        mat->SetDimension(dim);
        mat->SetBigNumber(1.e10);
        cmesh->InsertMaterialObject(mat);
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
    // // Print pressure mesh
    // std::ofstream myfile("PressureMesh.txt");
    // cmesh->Print(myfile);

    return cmesh;
}


// Pressure computational mesh
TPZCompMesh *PressureCMeshDFB(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh)
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    std::set<int> matIdNeumann;

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

    // // Print pressure mesh
    // std::ofstream myfile("PressureMeshNew.txt");
    // cmesh->Print(myfile);

    // //Prints computational mesh properties
    // std::stringstream vtk_name;
    // vtk_name    << "PressNew" << ".vtk";
    // std::ofstream vtkfile(vtk_name.str().c_str());
    // TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

    return cmesh;
}

// Multiphysics computational mesh
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
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
    auto * BCond0 = mat->CreateBC(mat, EInjection, 0, val1, val2);
    TPZFMatrix<STATE> val3(1,1,1.); 
    TPZManVector<STATE> val4(1,0.);
    auto * BCond1 = mat->CreateBC(mat, EProduction, 0, val3, val4);
    BCond0->SetForcingFunctionBC(exactSol);
    BCond1->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond0);

    TPZFMatrix<STATE> val5(1,1,1.);
    TPZManVector<STATE> val6(1,0.);
    auto * BCond2 = mat->CreateBC(mat, EBottom, 0, val5, val6);
    auto * BCond3 = mat->CreateBC(mat, ETop, 0, val5, val6);
    auto * BCond4 = mat->CreateBC(mat, ELeft, 0, val5, val6);
    auto * BCond5 = mat->CreateBC(mat, ERight, 0, val5, val6);
    BCond2->SetForcingFunctionBC(exactSol);
    BCond3->SetForcingFunctionBC(exactSol);
    BCond4->SetForcingFunctionBC(exactSol);
    BCond5->SetForcingFunctionBC(exactSol);
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
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);

    // // Prints Multiphysics mesh
    // std::ofstream myfile("MultiPhysicsMesh.txt");
    // cmesh->Print(myfile);

    //Prints computational mesh properties
    // std::stringstream vtk_name;
    // vtk_name    << "MultiPhysics" << ".vtk";
    // std::ofstream vtkfile(vtk_name.str().c_str());
    // TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

    return cmesh;
}

// Multiphysics computational mesh
TPZMultiphysicsCompMesh *MultiphysicCMeshDFB(int dim, int pOrder, std::set<int> &matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
    gmesh->ResetReference();
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
    auto * BCond0 = mat->CreateBC(mat, EInjection, 1, val1, val2);
    TPZFMatrix<STATE> val3(1,1,1.); 
    TPZManVector<STATE> val4(1,0.);
    auto * BCond1 = mat->CreateBC(mat, EProduction, 1, val3, val4);
    BCond0->SetForcingFunctionBC(exactSol);
    BCond1->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond0);

    TPZFMatrix<STATE> val5(1,1,1.);
    TPZManVector<STATE> val6(1,0.);
    auto * BCond2 = mat->CreateBC(mat, EBottom, 1, val5, val6);
    auto * BCond3 = mat->CreateBC(mat, ETop, 1, val5, val6);
    auto * BCond4 = mat->CreateBC(mat, ELeft, 1, val5, val6);
    auto * BCond5 = mat->CreateBC(mat, ERight, 1, val5, val6);
    // auto * BCond6 = mat->CreateBC(mat, EPont, 0, val5, val6);
    BCond2->SetForcingFunctionBC(exactSol);
    BCond3->SetForcingFunctionBC(exactSol);
    BCond4->SetForcingFunctionBC(exactSol);
    BCond5->SetForcingFunctionBC(exactSol);
    // BCond6->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    // cmesh->InsertMaterialObject(BCond6);

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
    HybridizeDomain(cmesh,gmesh);

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
  // TPZStepSolver<STATE> * precond = new TPZStepSolver<STATE>( matskl.Create() ); step.ShareMatrix( *precond ); precond->SetJacobi(1, 0.0, 0);
  // TPZStepSolver<STATE> jac;
  // REAL tol = 1.e-10;
  // jac.SetSSOR(1,1.1,0.,0);
  // jac.ShareMatrix(step);
  // an.SetSolver(step);
  // step.SetGMRES(20,20, *precond, tol, 0);
  // step.SetCG(2000, *precond, tol, 0);
  step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
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

TPZCompMesh *CMeshH1(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh)
{
  //Creates cmesh object
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetAllCreateFunctionsContinuous();

  //Sets materials
  TPZMatPoisson<> *mat = new TPZMatPoisson<>(EDomain,dim);
//  mat->SetPermeability(1);
  cmesh->InsertMaterialObject(mat);

  //Insert boundary conditions
  TPZFMatrix<STATE> val1(1,1,1.);
  TPZManVector<STATE> val2(2,0.);
  auto * BCond = mat->CreateBC(mat, EInjection, 0, val1, val2);
  TPZFMatrix<STATE> val3(1,1,1.);
  TPZManVector<STATE> val4(2,0.);
  auto * BCond1 = mat->CreateBC(mat, EProduction, 0, val3, val4);
  BCond->SetForcingFunctionBC(exactSolH1);
  BCond1->SetForcingFunctionBC(exactSolH1);
  cmesh->InsertMaterialObject(BCond);
  cmesh->InsertMaterialObject(BCond1);
  
  cmesh->SetDefaultOrder(pOrder);
  cmesh->AutoBuild();

  return cmesh;
}

void PrintResultsH1(int dim, TPZLinearAnalysis &an)
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
  an.SetExact(exactSolError,solOrder);
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
  an.SetExact(exactSolError,solOrder);
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
  //Creates cmesh object
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetAllCreateFunctionsContinuous();

  //Sets materials
  auto *mat = new TPZMatDivFreeBubbles<STATE>(EDomain,dim);
  cmesh->InsertMaterialObject(mat);
   
  //Insert boundary conditions
  TPZFMatrix<STATE> val1(1,1,1.);
  TPZManVector<STATE> val2(1,0.);
  TPZManVector<STATE> val5(1,0.);
  constexpr int boundType{0};
  auto * BCond = mat->CreateBC(mat, EInjection, 0, val1, val5);//Injection
  TPZFMatrix<STATE> val3(1,1,1.);
  TPZManVector<STATE> val4(1,0.);
  auto * BCond1 = mat->CreateBC(mat, EProduction, 0, val3, val5);//Production
  BCond->SetForcingFunctionBC(exactSol);
  BCond1->SetForcingFunctionBC(exactSol);
  cmesh->InsertMaterialObject(BCond);
  cmesh->InsertMaterialObject(BCond1);
  mat -> SetBigNumber(1.e10);
  
  auto * BCond2 = mat->CreateBC(mat, EBottom, 0, val3, val5);//Bottom
  auto * BCond3 = mat->CreateBC(mat, ETop, 0, val3, val5);//Top
  auto * BCond4 = mat->CreateBC(mat, ELeft, 0, val3, val5);//Left
  auto * BCond5 = mat->CreateBC(mat, ERight, 0, val3, val5);//Right
  BCond2->SetForcingFunctionBC(exactSol);
  BCond3->SetForcingFunctionBC(exactSol);
  BCond4->SetForcingFunctionBC(exactSol);
  BCond5->SetForcingFunctionBC(exactSol);
  cmesh->InsertMaterialObject(BCond2);
  cmesh->InsertMaterialObject(BCond3);
  cmesh->InsertMaterialObject(BCond4);
  cmesh->InsertMaterialObject(BCond5);

  // auto *mat2 = new TPZL2Projection<>(EPont,0,1);
  // cmesh->InsertMaterialObject(mat2);

  //Prints computational mesh properties
  // std::stringstream text_name;
  // std::stringstream vtk_name;
  // text_name   << "geometry" << ".txt";
  // vtk_name    << "geometry" << ".vtk";
  // std::ofstream textfile(text_name.str().c_str());
  // gmesh->Print(textfile);
  // std::ofstream vtkfile(vtk_name.str().c_str());
  // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

  cmesh->SetDefaultOrder(pOrder);
  cmesh->AutoBuild();

  // Prints DFB mesh
  // std::ofstream myfile("DFBMesh.txt");
  // CMeshDFB->Print(myfile);

  return cmesh;
}

void PrintResultsMultiphysicDFB(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh)
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
        std::set<int> neighmat = {EInjection, EProduction, EBottom, ETop, ELeft, ERight};
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

      // Creating interface elements
    TPZCompMesh* pressuremesh = cmesh->MeshVector()[1];
    cmesh->Reference()->ResetReference();
    pressuremesh->LoadReferences();
    const int lagrangematid = EIntfaceLeft;
    for (auto cel : pressuremesh->ElementVec()) {
        const int celmatid = cel->Material()->Id();
        if (!cel || celmatid != lagrangematid ) {
            continue;
        }
        TPZGeoEl* gel = cel->Reference();
        
        int dim = gel->Dimension()+1;
        
        cmesh->Reference()->ResetReference();
        cmesh->LoadReferences();
        
        TPZStack<TPZCompElSide> celstack;
    //    TPZGeoEl *gel = cel->Reference();
        TPZCompEl *mphysics = gel->Reference();
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZCompElSide celside = gelside.Reference();
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        int count = 0;
        for (auto &celstackside : celstack) {
            if (celstackside.Reference().Element()->Dimension() == dim - 1) {
                int matid = EIntfaceRight;
                if(count == 1) matid = EIntfaceRight;
                TPZGeoElBC gbc(gelside, matid);
                // check if the right side has a dependency
                TPZCompEl *celneigh = celstackside.Element();
                // std::cout << "COnn = " << celneigh->NConnects() << " " << std::endl;
                // if (celneigh->NConnects() != 1) {
                //     DebugStop();
                // }
                int64_t index;
                TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, celstackside);
                count++;
            }
        }
        if (count == 1)
        {
            TPZCompElSide clarge = gelside.LowerLevelCompElementList2(false);
            if(!clarge) DebugStop();
            TPZGeoElSide glarge = clarge.Reference();
            if (glarge.Element()->Dimension() == dim) {
                TPZGeoElSide neighbour = glarge.Neighbour();
                while(neighbour != glarge)
                {
                    if (neighbour.Element()->Dimension() < dim) {
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                if(neighbour == glarge) DebugStop();
                glarge = neighbour;
            }
            clarge = glarge.Reference();
            if(!clarge) DebugStop();
            TPZGeoElBC gbc(gelside, EIntfaceRight);

            int64_t index;
            TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, clarge);
            count++;
        }
        
        pressuremesh->InitializeBlock();
        
    }
}