#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

// C++ includes
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>

// PZ includes
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZGmshReader.h>
#include "Poisson/TPZMatPoisson.h" //for TPZMatLaplacian
#include <TPZNullMaterial.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZLinearAnalysis.h"
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "TPZMultiphysicsCompMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"
#include "TPZHybridizeHDiv.h"
#include "pzlog.h"

// ----- Functions -----
using namespace std;

void readGeoMesh(string& filename, TPZGeoMesh* gmesh);
TPZCompMesh *FluxCMesh(int dim, int pOrder, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMesh(int dim, int pOrder, TPZGeoMesh *gmesh);
TPZCompMesh *MultiphysicCMesh(int dim, int pOrder, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
void HybridizeMiddle(TPZCompMesh* cmesh);

TPZCompMesh *CMeshH1(int dim, int pOrder, TPZGeoMesh *gmesh);
void PrintResultsH1(int dim, TPZLinearAnalysis &an);

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh);

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
TPZLogger::InitializePZLOG();
#endif
    
    int dim = 2;
    int pOrder = 1;
    
    //................. Read mesh from gmsh ................................
    TPZGeoMesh* gmesh = new TPZGeoMesh();
    string filename = "../mesh/Case1FracSimple.msh";
    readGeoMesh(filename,gmesh);
    
    //.................................Hdiv.................................
    //Flux mesh
    TPZCompMesh * cmeshflux = FluxCMesh(dim,pOrder,gmesh);
//    ofstream outvtkcmeshflux("cmeshflux.vtk");
//    TPZVTKGeoMesh::PrintCMeshVTK(cmeshflux, outvtkcmeshflux);
    
    //Pressure mesh
    TPZCompMesh * cmeshpressure= PressureCMesh(dim,pOrder,gmesh);
//    ofstream outvtkcmeshp("cmeshpressure.vtk");
//    TPZVTKGeoMesh::PrintCMeshVTK(cmeshpressure, outvtkcmeshp);

    
    //Multiphysics mesh
    TPZManVector< TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmeshflux;
    meshvector[1] = cmeshpressure;
    TPZCompMesh * cmesh = MultiphysicCMesh(dim,pOrder,meshvector,gmesh);
//    ofstream outvtkcmeshmult("cmeshmult.vtk");
//    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, outvtkcmeshmult);
    
    HybridizeMiddle(cmesh);
    
    //Solve Multiphysics
    TPZLinearAnalysis an(cmesh,true);
    SolveProblemDirect(an,cmesh);
    cmesh->UpdatePreviousState(-1.);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
//    ofstream outmult("outmult.txt");
//    cmesh->Print(outmult);
    
    //Print results
    PrintResultsMultiphysic(dim,meshvector,an,cmesh);
    
    //..................................H1..................................
    // //Creates H1 problem
    TPZCompMesh * cmeshH1 = CMeshH1(dim,pOrder,gmesh);
    
    // // //Solve H1
    TPZLinearAnalysis anH1(cmeshH1,false);
    SolveProblemDirect(anH1,cmeshH1);
    
    // // //Print results
    PrintResultsH1(dim,anH1);
    std::ofstream out4("meshH1.txt");
    
//    delete gmesh;
//    delete cmesh;
//    delete cmeshflux;
//    delete cmeshpressure;
    return 0;
}

// ----------------- End of Main -----------------

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------


void HybridizeMiddle(TPZCompMesh* cmesh) {

    TPZMultiphysicsCompMesh* mmesh = dynamic_cast<TPZMultiphysicsCompMesh*>(cmesh);
    if (!mmesh)
        DebugStop();
    
    
    TPZCompMesh* fluxmesh = mmesh->MeshVector()[0];
    TPZVec<TPZCompMesh *> &meshvec_Hybrid = mmesh->MeshVector();
    TPZHybridizeHDiv hybridizer;
    
    // now we find the intersection and hybridize it
    int dim = fluxmesh->Dimension();
    int64_t nel = fluxmesh->NElements();
    for (int64_t iel = 0; iel < nel; iel++) {
        TPZCompEl* cel = fluxmesh->Element(iel);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel || intel->Reference()->Dimension() != dim) {
            continue;
        }
        
        // loop over the side of dimension dim-1
        TPZGeoEl *gel = intel->Reference();
        for (int side = gel->NCornerNodes(); side < gel->NSides() - 1; side++) {
            if (gel->SideDimension(side) != dim - 1) {
                continue;
            }
            
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neigh = gelside.Neighbour();
            
            while(neigh != gelside){
                int neighmatid = neigh.Element()->MaterialId();
                int neighdim = neigh.Dimension();
                
                if (neighmatid == 15 && neighdim == dim-1) {
                    cout << "\nElement with ID " << gel->Id() << " and index " << gel->Index() << " has side number " << side << " with dim = " << neigh.Dimension() << " touching the requested matID" << endl;
                    cout << "===> Hybridizing the interface now..." << endl;
                    TPZCompElSide celsideleft(intel, side);
                    // hybridizer.HybridizeInterface(celsideleft,intel,side,meshvec_Hybrid);
//                    TPZCompElSide neighcomp = RightElement(intel, side);
//                    if (neighcomp) {
//                        // SplitConnects returns the geometric element index and interpolation order
//                        pressures.push_back(SplitConnects(celside, neighcomp, meshvec_Hybrid));
//                    }
                }
                neigh = neigh.Neighbour();
            } // while
            
            
        }
        
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

    //    mat -> fBigNumber = 1.e10;
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,4.);
    auto * BCond0 = mat->CreateBC(mat, -3, 1, val1, val2);
    auto * BCond1 = mat->CreateBC(mat, -2, 0, val1, val2);
    auto * BCond2 = mat->CreateBC(mat, -1, 0, val1, val4);
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);

    mat->SetDimension(dim);
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
            
    
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
    // Print flux mesh
    // std::ofstream myfile("FluxMesh.txt");
    // cmesh->Print(myfile);
    
    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh *PressureCMesh(int dim, int pOrder, TPZGeoMesh *gmesh)
{
    gmesh->ResetReference();
    //Change if not triangular elements
    bool fTriang = false;
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(1);
    mat->SetDimension(dim);
    cmesh->InsertMaterialObject(mat);
    // mat -> fBigNumber = 1.e10;
    cmesh->SetAllCreateFunctionsDiscontinuous();
//    cmesh->SetAllCreateFunctionsContinuous();
//    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);

    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,4.);
    constexpr int boundType{1};
    constexpr int boundType0{0};
    auto * BCond0 = mat->CreateBC(mat, -3, 1, val1, val2);
    auto * BCond1 = mat->CreateBC(mat, -2, 0, val1, val2);
    auto * BCond2 = mat->CreateBC(mat, -1, 0, val1, val4);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond2);

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

TPZCompMesh *MultiphysicCMesh(int dim, int pOrder, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    auto mat = new TPZMixedDarcyFlow(1, dim);
    
    mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);
//    mat -> fBigNumber = 1.e10;
    
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,4.);
    auto * BCond0 = mat->CreateBC(mat, -3, 1, val1, val2);
    auto * BCond1 = mat->CreateBC(mat, -2, 0, val1, val2);
    auto * BCond2 = mat->CreateBC(mat, -1, 0, val1, val4);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond0);
    cmesh->InsertMaterialObject(BCond2);

    
    TPZManVector<int> active(2,1);
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();

    cmesh->BuildMultiphysicsSpace(active, meshvector);
//    cmesh->AdjustBoundaryElements();
//    cmesh->CleanUpUnconnectedNodes();
    
//    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
//    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);
    
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
    // Print mesh properties
    // std::ofstream out("mesh.txt");
    // an.Print("nothing",out);
    
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh *CMeshH1(int dim, int pOrder, TPZGeoMesh *gmesh)
{
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->SetDimModel(dim);
    
    
    TPZMatPoisson<> *mat = new TPZMatPoisson<>(1,dim);
    // TPZMixedDarcyFlow *mat = new TPZMixedDarcyFlow(matIdVec[0],dim);
    // mat->SetPermeabilityFunction(1.);
//    mat->fBigNumber = 1.e12;
    cmesh->InsertMaterialObject(mat);
    
    //Insert boundary conditions
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    TPZManVector<STATE> val4(1,4.);
    auto * BCond = mat->CreateBC(mat, -3, 1, val1, val2);
    auto * BCond1 = mat->CreateBC(mat, -2, 0, val1, val2);
    auto * BCond2 = mat->CreateBC(mat, -1, 0, val1, val4);
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
    dim_name_and_physical_tagFine[2]["frac"] = 1;
    dim_name_and_physical_tagFine[1]["inlet"] = -1;
    dim_name_and_physical_tagFine[1]["outlet"] = -2;
    dim_name_and_physical_tagFine[1]["noflux"] = -3;
    dim_name_and_physical_tagFine[1]["intersection"] = 15;
    reader.SetDimNamePhysical(dim_name_and_physical_tagFine);
    reader.GeometricGmshMesh4(filename,gmesh);
    gmesh->SetDimension(2);
    ofstream outvtk("geoMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outvtk);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
