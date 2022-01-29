#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZMultiphysicsCompMesh.h>
#include <TPZLinearAnalysis.h>
#include <TPZGmshReader.h>
#include <TPZVTKGeoMesh.h>
#include <TPZCompElDisc.h>
#include <TPZNullMaterial.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>// for Hdiv problem
#include <Poisson/TPZMatPoisson.h>
#include <pzbuildmultiphysicsmesh.h>
#include <TPZNullMaterialCS.h>
#include <pzlog.h>
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZLagrangeMultiplierCS.h>
#include <pzshapelinear.h>
#include <pzshapepoint.h>
#include <pzshapequad.h>
#include <pzshapetriang.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage

#include "divfree_config.h"
#include "TPZMatDivFreeBubbles.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZCompElKernelHDiv.h"
#include "TPZCompElKernelHDivBC.h"
#include "TPZKernelHdivHybridizer.h"
#include "TPZKernelHdivUtils.h"
#include "TPZCompElH1.h"

TPZCompMesh *FluxCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *FluxCMeshDFB(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMeshDFB(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
TPZMultiphysicsCompMesh *MultiphysicCMeshDFB(int dim, int pOrder, std::set<int> &matIdVec, std::set<int> &matIdNeumann, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
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

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

enum EMatid {ENone, EDomain, EBottom,  ETop, ELeft, ERight, EPont, EWrap, EIntface, EPressureHyb};

int main(int argc, char* argv[]){
    //dimension of the problem
    constexpr int dim{2};
    constexpr int pOrder{2};

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
        stringtoint[1]["InjectionWell"] = 2;
        stringtoint[1]["ProductionWell"] = 3;
        stringtoint[1]["BottomLine"] = 4;
        stringtoint[1]["TopLine"] = 5;
        stringtoint[1]["LeftLine"] = 6;
        stringtoint[1]["RightLine"] = 7;
        stringtoint[0]["Point"] = 8;
        stringtoint[1]["BottomLine2"] = 9;
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(string(MESHDIR)+"1element.msh",gmesh);
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
  
    //.................................Hdiv.................................
    TPZCompMesh * cmeshflux = 0;
    TPZCompMesh * cmeshpressure = 0;
    {   
        TPZKernelHdivUtils<STATE> util;
        std::set<int> matIdVecHdiv={EDomain,EBottom,ETop,ELeft,ERight};
        std::set<int> matIdNeumannHdiv;
        
        //Flux mesh
        cmeshflux = FluxCMesh(dim,pOrder,matIdVecHdiv,gmesh);
        // std::cout << "FLUX \n";
        // util.PrintCMeshConnects(cmeshflux);

        //Pressure mesh
        cmeshpressure = PressureCMesh(dim,pOrder-1,matIdVecHdiv,gmesh);
        // std::cout << "PRESSURE \n";
        // util.PrintCMeshConnects(cmeshpressure);

        //Multiphysics mesh
        TPZManVector< TPZCompMesh *, 2> meshvector(2);
        meshvector[0] = cmeshflux;
        meshvector[1] = cmeshpressure;
        TPZCompMesh * cmesh = MultiphysicCMesh(dim,pOrder,matIdVecHdiv,meshvector,gmesh);
        std::cout << "Number of equations = " << cmesh->NEquations() << std::endl;
        
        // std::cout << "MULTIPHYSICS \n";
        // util.PrintCMeshConnects(cmesh);

        //Solve Multiphysics
        TPZLinearAnalysis an(cmesh,true);
        SolveProblemDirect(an,cmesh);
        
        //Print results
        PrintResultsMultiphysic(dim,meshvector,an,cmesh);

        std::ofstream anPostProcessFileHdiv("postprocessHdiv.txt");
        ComputeErrorHdiv(an,anPostProcessFileHdiv);
    }

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
    // {
    //     TPZKernelHdivUtils util;
    //     TPZKernelHdivHybridizer hybridizer;

    //     std::set<int> matIdNeumannNew={EInjection,EProduction,EBottom,ETop,ELeft,ERight};
    //     hybridizer.SetPeriferalMaterialIds(EWrap,EPressureHyb,EIntface,EPont,EDomain);
    //     hybridizer.CreateWrapElements(gmesh,matIdNeumannNew,true);
    //     // hybridizer.PrintGeoMesh(gmesh);

    //     // nesta configuracao o material ETop eh substituido por EWrap
    //     std::set<int> matIdVecNew={EDomain,EPont,EWrap};
    //     // criamos elementos tipo ETop para a pressao
    //     matIdNeumannNew.insert(EPressureHyb);

    //     //Flux mesh
    //     TPZCompMesh * cmeshfluxDFB = FluxCMeshDFB(dim,pOrder+1,matIdVecNew,gmesh);
    //     // hybridizer.SemiHybridizeFlux(cmeshfluxDFB,matIdNeumannNew);
    //     // util.PrintCMeshConnects(cmeshfluxNew);

    //     //Pressure mesh
    //     TPZCompMesh * cmeshpressureDFB = PressureCMeshDFB(dim,pOrder,matIdNeumannNew,gmesh);

    //     //Multiphysics mesh
    //     TPZManVector< TPZCompMesh *, 2> meshvectorDFB(2);
    //     meshvectorDFB[0] = cmeshfluxDFB;
    //     meshvectorDFB[1] = cmeshpressureDFB;
    //     auto * cmeshNew = MultiphysicCMeshDFB(dim,pOrder+1,matIdVecNew,matIdNeumannNew,meshvectorDFB,gmesh);
    //     hybridizer.CreateMultiphysicsInterfaceElements(cmeshNew,gmesh,meshvectorDFB,matIdNeumannNew);
    //     // util.PrintCMeshConnects(cmeshNew);

    //     //Solve Multiphysics
    //     TPZLinearAnalysis anNew(cmeshNew,false);
    //     SolveProblemDirect(anNew,cmeshNew);

    //     //Print results
    //     PrintResultsMultiphysicDFB(dim,meshvectorDFB,anNew,cmeshNew);
    //     std::ofstream out4("mesh_MDFB.txt");
    //     anNew.Print("nothing",out4);
        
    //     std::ofstream anPostProcessFileMDFB("postprocessMDFB.txt");
    //     ComputeErrorHdiv(anNew,anPostProcessFileMDFB);
    // }

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
  u[0] = x;
gradU(0,0) = 1.;
gradU(1,0) = 0.;
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
//   u[0]= log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
//   gradU(0,0) = -(x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2)));
//   gradU(1,0) = -(y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2)));
//Harmonic2
REAL a1 = 0.25;
REAL alpha = M_PI/2;
u[0]= 
  x*a1*cos(x*alpha)*cosh(y*alpha) + y*a1*sin(x*alpha)*sinh(y*alpha);
    gradU(0,0) = a1*cos(alpha*x)*cosh(alpha*y) - 
  a1*alpha* x*cosh(alpha*y)*sin(alpha*x) + 
  a1*alpha*y*cos(alpha*x)*sinh(alpha*y);
    gradU(1,0) = a1*alpha* y*cosh(alpha*y) *sin(alpha* x) + 
  a1* alpha *x* cos(alpha* x)* sinh(alpha* y) + 
  a1* sin(alpha* x) * sinh(alpha* y);
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
    // cmesh->ApproxSpace().SetAllCreateFunctionsHDivConstant(dim);
    // cmesh->ApproxSpace().SetAllCreateFunctionsHDivKernel(dim);
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
        cmesh->InsertMaterialObject(mat);
        mat->SetDimension(dim);
        mat->SetBigNumber(1.e10);
    }

    //Creates computational elements
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) DebugStop();
        auto type = gel -> Type();
        int64_t index;
        auto matid = gel->MaterialId();
        if(allmat.find(matid) == allmat.end()) continue;

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
            new TPZCompElKernelHDiv<TPZShapeQuad>(*cmesh,gel);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(2);

            for (int side=0; side<gel->NSides(); side++){
                
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();

                if (neighbour.Element()->MaterialId() == EPont){
                    new TPZCompElH1<TPZShapePoint>(*cmesh,neighbour.Element());
                    TPZMaterial *mat = cmesh->FindMaterial(EPont);
                    TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                }
                
                if (neighbour.Element()->Dimension() != dim-1) continue;

                if (neighbour.Element()->MaterialId() == EWrap){
                    new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,neighbour.Element());
                    TPZMaterial *mat = cmesh->FindMaterial(EWrap);
                    TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                    nullmat->SetDimension(1);
                    neighbour.Element()->ResetReference();
                } else {
                    if (allmat.find(neighbour.Element()->MaterialId()) != allmat.end()){
                        new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,neighbour.Element());
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
            new TPZCompElKernelHDiv<TPZShapeTriang>(*cmesh,gel);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(2);
        }
    }    

    cmesh->ExpandSolution();
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
            // DebugStop();
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

    // Sets matid to BC geometric elements
    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        mat->SetDimension(1);
        cmesh->InsertMaterialObject(mat);
        mat->SetBigNumber(1.e10);
        matIdNeumann.insert(*it); 
    }

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

// Multiphysics computational mesh
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    auto mat = new TPZMixedDarcyFlow(EDomain, dim);
    mat->SetConstantPermeability(1.);
    
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
        
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    // auto * BCond0 = mat->CreateBC(mat, EInjection, 0, val1, val2);
    TPZFMatrix<STATE> val3(1,1,1.); 
    TPZManVector<STATE> val4(1,0.);
    // auto * BCond1 = mat->CreateBC(mat, EProduction, 0, val3, val4);
    // BCond0->SetForcingFunctionBC(exactSol);
    // BCond1->SetForcingFunctionBC(exactSol);
    // cmesh->InsertMaterialObject(BCond1);
    // cmesh->InsertMaterialObject(BCond0);

    TPZFMatrix<STATE> val5(1,1,1.);
    TPZManVector<STATE> val6(1,0.);
    auto * BCond2 = mat->CreateBC(mat, EBottom, 0, val5, val6);
    auto * BCond3 = mat->CreateBC(mat, ETop, 0, val5, val6);
    auto * BCond4 = mat->CreateBC(mat, ELeft, 0, val5, val6);
    auto * BCond5 = mat->CreateBC(mat, ERight, 0, val5, val6);
    auto * BCond6 = mat->CreateBC(mat, EPont, 0, val5, val6);
    BCond2->SetForcingFunctionBC(exactSolError);
    BCond3->SetForcingFunctionBC(exactSolError);
    BCond4->SetForcingFunctionBC(exactSolError);
    BCond5->SetForcingFunctionBC(exactSolError);
    BCond6->SetForcingFunctionBC(exactSolError);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    cmesh->InsertMaterialObject(BCond6);

    TPZManVector<int> active(2,1);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->BuildMultiphysicsSpace(active, meshvector);

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
TPZMultiphysicsCompMesh *MultiphysicCMeshDFB(int dim, int pOrder, std::set<int> &matIdVec, std::set<int> &matIdNeumann, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();


    // eh preciso criar materiais para todos os valores referenciados no enum
    auto mat = new TPZMixedDarcyFlow(EDomain, dim);
    mat->SetConstantPermeability(1.);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
    
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    // auto * BCond0 = mat->CreateBC(mat, EInjection, 1, val1, val2);
    TPZFMatrix<STATE> val3(1,1,1.); 
    TPZManVector<STATE> val4(1,0.);
    // auto * BCond1 = mat->CreateBC(mat, EProduction, 1, val3, val4);
    // BCond0->SetForcingFunctionBC(exactSolError);
    // BCond1->SetForcingFunctionBC(exactSolError);
    // cmesh->InsertMaterialObject(BCond1);
    // cmesh->InsertMaterialObject(BCond0);

    TPZFMatrix<STATE> val5(1,1,1.);
    TPZManVector<STATE> val6(1,0.);
    auto * BCond2 = mat->CreateBC(mat, EBottom, 1, val5, val6);
    auto * BCond3 = mat->CreateBC(mat, ETop, 1, val5, val6);
    auto * BCond4 = mat->CreateBC(mat, ELeft, 1, val5, val6);
    auto * BCond5 = mat->CreateBC(mat, ERight, 1, val5, val6);
    auto * BCond6 = mat->CreateBC(mat, EPont, 0, val5, val6);
    BCond2->SetForcingFunctionBC(exactSolError);
    BCond3->SetForcingFunctionBC(exactSolError);
    BCond4->SetForcingFunctionBC(exactSolError);
    BCond5->SetForcingFunctionBC(exactSolError);
    BCond6->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    cmesh->InsertMaterialObject(BCond4);
    cmesh->InsertMaterialObject(BCond5);
    cmesh->InsertMaterialObject(BCond6);

    // the wrap material is a null material (does nothing)
    auto * nullmat = new TPZNullMaterialCS<>(EWrap,1,1);
    cmesh->InsertMaterialObject(nullmat);

    auto *matL2 = new TPZL2ProjectionCS<>(EPont,0,1);
    cmesh->InsertMaterialObject(matL2);

    auto * nullmat3 = new TPZNullMaterialCS<>(EPressureHyb,1,1);
    cmesh->InsertMaterialObject(nullmat3);

    TPZManVector<int> active(2,1);
    active[0]=1;
    active[1]=1;
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes(); 

    auto mat3 = new TPZLagrangeMultiplierCS<STATE>(EIntface, dim-1);
    cmesh->InsertMaterialObject(mat3);


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
  constexpr int nThreads{10};
  //defines storage scheme to be used for the FEM matrices
  //in this case, a symmetric skyline matrix is used
//   TPZSkylineStructMatrix<STATE> matskl(cmesh);
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matskl(cmesh);
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
//   auto * BCond = mat->CreateBC(mat, EInjection, 0, val1, val2);
  TPZFMatrix<STATE> val3(1,1,1.);
  TPZManVector<STATE> val4(2,0.);
//   auto * BCond1 = mat->CreateBC(mat, EProduction, 0, val3, val4);
//   BCond->SetForcingFunctionBC(exactSolH1);
//   BCond1->SetForcingFunctionBC(exactSolH1);
//   cmesh->InsertMaterialObject(BCond);
//   cmesh->InsertMaterialObject(BCond1);
  
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
  std::cout << std::setprecision(15) << "H1 Norm = " << error[0]<<'\n';
  std::cout << std::setprecision(15) << "L1 Norm = " << error[1]<<'\n'; 
  std::cout << std::setprecision(15) << "H1 Seminorm = " << error[2] << "\n\n";

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
//   auto * BCond = mat->CreateBC(mat, EInjection, 0, val1, val5);//Injection
  TPZFMatrix<STATE> val3(1,1,1.);
  TPZManVector<STATE> val4(1,0.);
//   auto * BCond1 = mat->CreateBC(mat, EProduction, 0, val3, val5);//Production
//   BCond->SetForcingFunctionBC(exactSol);
//   BCond1->SetForcingFunctionBC(exactSol);
//   cmesh->InsertMaterialObject(BCond);
//   cmesh->InsertMaterialObject(BCond1);
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
