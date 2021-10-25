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
#include <pzshapetriang.h>
#include <TPZHybridizeHDiv.h>

#include "divfree_config.h"
#include "TPZMatDivFreeBubbles.h"
#include "TPZL2ProjectionCS.h"
#include "TPZCompElKernelHdiv.h"
#include "TPZCompElKernelHdivBC.h"

using namespace std;
// ----- Functions -----
void readGeoMesh(string& filename, TPZGeoMesh* gmesh);
TPZCompMesh *FluxCMesh(int dim, int pOrder, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMesh(int dim, int pOrder, TPZGeoMesh *gmesh);
TPZCompMesh *MultiphysicCMesh(int dim, int pOrder, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh, TPZHybridizeHDiv& hybridizer);
void HybridizeIntersection(TPZHybridizeHDiv& hybridizer, TPZVec<TPZCompMesh*>& cmesh);

TPZCompMesh *CMeshH1(int dim, int pOrder, TPZGeoMesh *gmesh);
void PrintResultsH1(int dim, TPZLinearAnalysis &an);

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh);

enum EMatid  {ENone, EDomain, EInlet, EOutlet, ENoflux, EIntersection};
// ----- End of Functions -----

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------

int main(int argc, char* argv[]){

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG("log4cxx.cfg");
#endif
    
    int dim = 2;
    int pOrder = 1;
    
    //................. Read mesh from gmsh ................................
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    string filename = MESHDIR;
    int caseSim = 2;
    switch (caseSim) {
        case 0:
            //        filename = "../mesh/Case2FracSimple.msh";
            //        filename = "../mesh/flemisch2lf_3frac_bug.msh";
            filename + "flemisch3_2frac.msh";
            break;
        case 1:
            filename + "Case1FracSimple.msh";
            break;
        case 2:
            filename + "1element.msh";
            break;
            
        default:
            break;

    }

    readGeoMesh(filename,gmesh);
    
    //.................................Hdiv.................................
    //Flux mesh
    const bool isRunHdiv = true;
    if (isRunHdiv) {
        TPZCompMesh * cmeshflux = FluxCMesh(dim,pOrder,gmesh);
        
        //Pressure mesh
        TPZCompMesh * cmeshpressure= PressureCMesh(dim,pOrder,gmesh);
            
        //Multiphysics mesh
        TPZManVector< TPZCompMesh *, 2> meshvector(2);
        meshvector[0] = cmeshflux;
        meshvector[1] = cmeshpressure;
        
        TPZHybridizeHDiv hybridizer(meshvector);
        HybridizeIntersection(hybridizer, meshvector);
        
        // Multiphysics mesh AND create interface elements at the requested intersection
        TPZCompMesh * cmesh = MultiphysicCMesh(dim,pOrder,meshvector,gmesh,hybridizer);
        
    //    ofstream outmult("multimesh.txt");
    //    cmesh->Print(outmult);
        
        //Solve Multiphysics
        TPZLinearAnalysis an(cmesh,true);
        SolveProblemDirect(an,cmesh);
        
    //    cmesh->UpdatePreviousState(-1.); NS: When do I need this??????
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
        
        //Print results
        PrintResultsMultiphysic(dim,meshvector,an,cmesh);
    }

    //..................................H1..................................
    // //Creates H1 problem
    TPZCompMesh * cmeshH1 = CMeshH1(dim,pOrder,gmesh);
    
    // // //Solve H1
    TPZLinearAnalysis anH1(cmeshH1,false);
    SolveProblemDirect(anH1,cmeshH1);
    
    // // //Print results
    PrintResultsH1(dim,anH1);
    std::ofstream out4("meshH1.txt");
    
    return 0;
}

// ----------------- End of Main -----------------

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


void HybridizeIntersection(TPZHybridizeHDiv& hybridizer, TPZVec<TPZCompMesh*>& meshvec_Hybrid) {

       
    TPZCompMesh* fluxmesh = meshvec_Hybrid[0];
    TPZGeoMesh* gmesh = fluxmesh->Reference();
    fluxmesh->LoadReferences();
    hybridizer.InsertPeriferalMaterialObjects(meshvec_Hybrid);
        
    int dim = fluxmesh->Dimension();
    for (auto gel : gmesh->ElementVec()) {
        if (gel->MaterialId() != EIntersection) {
            continue;
        }
        if (gel->Dimension() != dim - 1) {
            DebugStop();
        }
        
        // Search for first neighbor that that is domain
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neigh = gelside.Neighbour();
        
        while(neigh != gelside){
            TPZGeoEl* gelneigh = neigh.Element();
            int neighmatid = gelneigh->MaterialId();
            int neighdim = gelneigh->Dimension();
            
            if (neighmatid == EDomain && neighdim == dim) {
                cout << "\nElement with ID " << gel->Id() << " and index " << gel->Index() << " is an intersection element" << endl;
                cout << "===> Trying to split the connects of the flux mesh and create pressure element..." << endl;
                TPZCompEl* celneigh = gelneigh->Reference();
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (celneigh);
                if (!intel)
                    DebugStop();
                
                const int side = neigh.Side();
                TPZCompElSide celsideleft(intel, side);
                bool isNewInterface = hybridizer.HybridizeInterface(celsideleft,intel,side,meshvec_Hybrid);
                if (isNewInterface) {
                    break;
                }
                else{
                    DebugStop();
                }

            }
            neigh = neigh.Neighbour();
        } // while
        
    }
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh *FluxCMesh(int dim, int pOrder, TPZGeoMesh *gmesh)
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(1);
    cmesh->InsertMaterialObject(mat);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,4.);
    auto * BCond0 = mat->CreateBC(mat, ENoflux, 1, val1, val2);
    auto * BCond1 = mat->CreateBC(mat, EOutlet, 0, val1, val2);
    auto * BCond2 = mat->CreateBC(mat, EInlet, 0, val1, val4);
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);

    mat->SetDimension(dim);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
            
    
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh *PressureCMesh(int dim, int pOrder, TPZGeoMesh *gmesh)
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(1);
    mat->SetDimension(dim);
    cmesh->InsertMaterialObject(mat);
//    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);

    // NS: Do I need to create boundary cond in the pressure mesh?????
    // NO WE DONT!!!!!

    cmesh->AutoBuild();

    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++){
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    std::ofstream myfile("PressureMesh.txt");
    cmesh->Print(myfile);
    
    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh *MultiphysicCMesh(int dim, int pOrder, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh, TPZHybridizeHDiv& hybridizer)
{
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    auto mat = new TPZMixedDarcyFlow(1, dim);
    
    mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);
    
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,1.);
    auto * BCond0 = mat->CreateBC(mat, ENoflux, 1, val1, val2);
    auto * BCond1 = mat->CreateBC(mat, EOutlet, 0, val1, val2);
    auto * BCond2 = mat->CreateBC(mat, EInlet, 0, val1, val4);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond2);
    
    hybridizer.InsertPeriferalMaterialObjects(cmesh);
    
    TPZManVector<int> active(2,1);
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    //    cmesh->AdjustBoundaryElements();
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes();
    
    // Creating interface elements
    TPZCompMesh* cmeshpressure = cmesh->MeshVector()[1];
    cmesh->Reference()->ResetReference();
    cmeshpressure->LoadReferences();
    const int lagrangematid = hybridizer.lagrangeInterfaceMatId();
    for (auto cel : cmeshpressure->ElementVec()) {
        const int celmatid = cel->Material()->Id();
        if (!cel || celmatid != lagrangematid ) {
            continue;
        }
        TPZGeoEl* gel = cel->Reference();
        
        hybridizer.CreateInterfaceElementsForGeoEl(cmesh, meshvector, gel);
        
    }
        
    // Prints Multiphysics mesh
    std::ofstream myfile("MultiPhysicsMesh.txt");
    cmesh->Print(myfile);
    
    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

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

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    
    TPZManVector<std::string,10> scalnames(2), vecnames(2);
    
    scalnames[0] = "Pressure";
    scalnames[1] = "ExactPressure";
    vecnames[0]= "Flux";
    vecnames[1]= "ExactFlux";
    
    int div = 0;
    std::string plotfile = "solutionHdiv.vtk";
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

auto exactSol = [](const TPZVec<REAL> &loc,
  TPZVec<STATE>&u,
  TPZFMatrix<STATE>&gradU){
  const auto &x=loc[0];
  const auto &y=loc[1];
  u[0]= 1-x;
  gradU(0,0) = 0.;
//  gradU(1,0) = (y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2)));
};

TPZCompMesh *CMeshH1(int dim, int pOrder, TPZGeoMesh *gmesh)
{
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDimModel(dim);
    
    TPZMatPoisson<> *mat = new TPZMatPoisson<>(1,dim);
    cmesh->InsertMaterialObject(mat);
    
    //Insert boundary conditions
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,4.);
    auto * BCond = mat->CreateBC(mat, ENoflux, 1, val1, val2);
    auto * BCond1 = mat->CreateBC(mat, EOutlet, 0, val1, val2);
    auto * BCond2 = mat->CreateBC(mat, EInlet, 0, val1, val4);
    BCond->SetForcingFunctionBC(exactSol);
    BCond1->SetForcingFunctionBC(exactSol);
    BCond2->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    
    cmesh->AutoBuild();
    
    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResultsH1(int dim, TPZLinearAnalysis &an)
{
    
    TPZVec<std::string> scalarVars(1), vectorVars(1);
    scalarVars[0] = "Solution";
    vectorVars[0] = "Derivative";
    an.DefineGraphMesh(dim,scalarVars,vectorVars,"SolutionH1.vtk");
    constexpr int resolution{0};
    an.PostProcess(resolution,dim);
    
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void readGeoMesh(string& filename, TPZGeoMesh* gmesh) {
    
    TPZGmshReader reader;
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tagFine(4);
    dim_name_and_physical_tagFine[2]["frac"] = EDomain;
    dim_name_and_physical_tagFine[1]["inlet"] = EInlet;
    dim_name_and_physical_tagFine[1]["outlet"] = EOutlet;
    dim_name_and_physical_tagFine[1]["noflux"] = ENoflux;
    dim_name_and_physical_tagFine[1]["intersection"] = EIntersection;
    
    // For flemisch case 2
//    dim_name_and_physical_tagFine[2]["Fracture10"] = EDomain;
//    dim_name_and_physical_tagFine[1]["BCfrac0"] = EInlet;
//    dim_name_and_physical_tagFine[1]["BCfrac1"] = EInlet;
//    dim_name_and_physical_tagFine[1]["BCfrac2"] = EInlet;
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
//    dim_name_and_physical_tagFine[1]["fracIntersection_0_2"] = EIntersection;
    
    // For flemisch case 3
    dim_name_and_physical_tagFine[2]["Fracture14"] = EDomain;
    dim_name_and_physical_tagFine[2]["Fracture17"] = EDomain;
    dim_name_and_physical_tagFine[1]["BCfrac0"] = EInlet;
    dim_name_and_physical_tagFine[1]["BCfrac1"] = EInlet;
    dim_name_and_physical_tagFine[1]["fracIntersection_0_1"] = EIntersection;
    
    // For 1 element case
    dim_name_and_physical_tagFine[2]["Surface"] = EDomain;
    dim_name_and_physical_tagFine[1]["Bottom"] = ENoflux;
    dim_name_and_physical_tagFine[1]["Right"] = EOutlet;
    dim_name_and_physical_tagFine[1]["Top"] = ENoflux;
    dim_name_and_physical_tagFine[1]["Left"] = EInlet;
//    stringtoint[0]["Point"] = 6;
//    stringtoint[1]["Top2"] = 7;

    reader.SetDimNamePhysical(dim_name_and_physical_tagFine);
    reader.GeometricGmshMesh4(filename,gmesh,false);
    gmesh->SetDimension(2);
    ofstream outvtk("geoMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outvtk);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

