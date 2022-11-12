
// --------------------- std includes ---------------------
#include <iostream>
#include <chrono>

// --------------------- PZ includes ---------------------
#include <pzgmesh.h>
#include <TPZGenGrid2D.h>
#include <TPZVTKGeoMesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include <pzbuildmultiphysicsmesh.h>
#include <TPZNullMaterial.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZMixedDarcyFlowOrtotropic.h>
#include <TPZLinearAnalysis.h>
#include <pzskylstrmatrix.h>
#include <TPZSSpStructMatrix.h>
#include <pzstepsolver.h>
#include <pzlog.h>
#include "TPZHDivApproxCreator.h"
#include "TPZMatRedSolver.h"
#include "TPZGenGrid3D.h"
#include "TPZVTKGenerator.h"

// --------------------- Global variables ---------------------
#define PROBLEM_3D

#ifndef PROBLEM_3D
constexpr int nx = 220;
constexpr int ny = 60;
constexpr int nz = 1;//85
constexpr int dim{2};
constexpr int n_cells = nx * ny * nz;
TPZManVector<REAL, n_cells> perm_vec(n_cells, 1);
#else
constexpr int nx = 60;
constexpr int ny = 220;
constexpr int nz = 85;
constexpr int dim{3};
constexpr int n_cells = nx * ny * nz;
TPZManVector<REAL, n_cells> perm_vec(n_cells * 3, 1);
#endif

enum EMatid {ENone,EDomain,EPressureLeft,EPressureRight,ENoFlux};

// --------------------- Namespaces ---------------------
using namespace std;

// --------------------- Functions ---------------------
void ReadSPE10CellPermeabilities(TPZVec<REAL>*perm_vec, int layer);
void ReadSPE10CellPermeabilities3D(TPZVec<REAL>*perm_vec);
TPZGeoMesh *CreateSPE10CoarseGeoMesh();
// STATE PermeabilityFunction(const TPZVec<REAL> &x);
TPZManVector<REAL,3> PermeabilityFunction(const TPZVec<REAL> &x);
// STATE PermeabilityFunction3D(const TPZVec<REAL> &x);
TPZManVector<REAL,3> PermeabilityFunction3D(const TPZVec<REAL> &x);
void PrintResultsVTK(const int dim, TPZLinearAnalysis &an, const std::string &plotfile);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);


//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
int main(){
    cout << "\n--------------------- Starting SPE10 simulations ---------------------\n" << endl;
    auto start_time = std::chrono::steady_clock::now();
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    cout << "\n--------------------- Reading Permeability Data ---------------------\n" << endl;
    constexpr int layer = 36; // Same as in repo ErrorEstimation/Projects/SPE10. I suppose is a layer with interesting permeability
    if (dim == 2){
        ReadSPE10CellPermeabilities(&perm_vec, layer);
    } else if (dim == 3){
        ReadSPE10CellPermeabilities3D(&perm_vec);
    } else {
        DebugStop();
    }
    REAL min = std::numeric_limits<int>::max(), max = -std::numeric_limits<int>::max();
    for (REAL perm : perm_vec) {
        if (perm > max) {
            max = perm;
        }
        if (perm < min) {
            min = perm;
        }
    }
    cout << "\nMinimum permeability = " << min << endl;
    cout << "Maximum permeability = " << max << endl;
    
    cout << "\n--------------------- Creating GeoMesh ---------------------\n" << endl;
    TPZGeoMesh *gmesh = CreateSPE10CoarseGeoMesh();
    ofstream outvtk("spe10gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outvtk);
    
    cout << "\n--------------------- Creating CompMeshes ---------------------\n" << endl;
    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.HdivFamily() = HDivFamily::EHDivConstant;
    hdivCreator.ProbType() = ProblemType::EDarcy;
    hdivCreator.IsRigidBodySpaces() = false;
    hdivCreator.SetDefaultOrder(1);
    hdivCreator.SetExtraInternalOrder(0);
    hdivCreator.SetShouldCondense(true);
    hdivCreator.HybridType() = HybridizationType::ESemi;

    //Insert Materials
    // TPZMixedDarcyFlow* matdarcy = new TPZMixedDarcyFlow(EDomain,dim);
    TPZMixedDarcyFlowOrtotropic* matdarcy = new TPZMixedDarcyFlowOrtotropic(EDomain,dim);
    std::function<TPZManVector<REAL,3>(const TPZVec<REAL> &coord)> func;
    // std::function<STATE(const TPZVec<REAL> &coord)> func;
    if (dim == 2){
        func = PermeabilityFunction;
    } else if (dim == 3){
        func = PermeabilityFunction3D;
    } else {
        DebugStop();
    }
    matdarcy->SetPermeabilityFunction(func);
    hdivCreator.InsertMaterialObject(matdarcy);

    constexpr int neumann{1};
    constexpr int dirichlet{0};
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2_0(1,0.), val2_1(1,1.);
    // domain bcs
    auto * BCondNoFlux = matdarcy->CreateBC(matdarcy, ENoFlux, neumann, val1, val2_0);
    hdivCreator.InsertMaterialObject(BCondNoFlux);
    
    auto * BCondLeft = matdarcy->CreateBC(matdarcy,EPressureLeft, dirichlet, val1, val2_0);
    hdivCreator.InsertMaterialObject(BCondLeft);
    
    auto * BCondRight = matdarcy->CreateBC(matdarcy,EPressureRight, dirichlet, val1, val2_1);
    hdivCreator.InsertMaterialObject(BCondRight);

    //Multiphysics mesh
    TPZMultiphysicsCompMesh *mpmesh = hdivCreator.CreateApproximationSpace();
    std::string txt = "cmesh.txt";
    std::ofstream myfile(txt);
    // mpmesh->Print(myfile);

    cout << "\n--------------------- Creating Analysis ---------------------\n" << endl;
    auto start_time_anal = std::chrono::steady_clock::now();
    std::cout << "Number of equations = " << mpmesh->NEquations() << std::endl;
    TPZLinearAnalysis an(mpmesh, false);
    auto total_time_anal = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_anal).count()/1000.;
    cout << "\nTotal time opt band = " << total_time_anal << " seconds" << endl;
        
    cout << "\n--------------------- Solving system ---------------------\n" << endl;
    if (hdivCreator.HybridType() == HybridizationType::ESemi){
        std::set<int> matBCAll={ENoFlux,EPressureLeft,EPressureRight};
        TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::ESparse);
        solver.Solve(std::cout);
    } else {
        SolveProblemDirect(an, mpmesh);
    }
    
    
    cout << "\n--------------------- Post processing ---------------------\n" << endl;
    auto start_time_pp = std::chrono::steady_clock::now();
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mpmesh->MeshVector(), mpmesh);
    string outres = "HDivResults.vtk";
    PrintResultsVTK(dim, an, outres);
    auto total_time_pp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_pp).count()/1000.;
    cout << "Total time post process = " << total_time_pp << " seconds" << endl;
        
    cout << "\n--------------------- End of execution ---------------------\n" << endl;
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count()/1000.;
    cout << "Total time = " << total_time << " seconds" << endl << endl;
    
    return 0;
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
void ReadSPE10CellPermeabilities(TPZVec<REAL> *perm_vec, const int layer) {
    // Fuction copied from ErrorEstimation/Projects/SPE10
    std::cout << "Reading permeability data...\n";

    std::ifstream perm_file("test/InputData/spe_perm.dat", std::ios::in);
    if (!perm_file) {
        std::cerr << "Unable to open input file\n";
        DebugStop();
    }

    int cell_id = 0;
    const auto n_cells = perm_vec->size();
    const auto start_line = 1 + n_cells * (layer - 1) / 6;

    int line_num = 0;
    int line_num2 = 0;
    while (perm_file) {
        line_num++;
        line_num2++;
        std::string line;
        std::getline(perm_file, line, '\n');

        if (line_num < start_line) continue;

        std::stringstream stream(line);
        for (int i = 0; i < 6; i++) {
            stream >> perm_vec->operator[](cell_id);
            cell_id++;
        }
        if (cell_id == n_cells) break;
    }
    std::ofstream outmat("perm.txt");
    for (int i = 0; i < perm_vec->size(); i++)
    {
        outmat << perm_vec->operator[](i) << std::endl;
    }
    
    std::cout << "Finished reading permeability data from input file!\n";
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
void ReadSPE10CellPermeabilities3D(TPZVec<REAL> *perm_vec) {
    // Fuction copied from ErrorEstimation/Projects/SPE10
    std::cout << "Reading permeability data...\n";

    std::ifstream perm_file("test/InputData/spe_perm.dat", std::ios::in);
    if (!perm_file) {
        std::cerr << "Unable to open input file\n";
        DebugStop();
    }

    int layer = 1;
    int cell_id = 0;
    const auto n_cells = perm_vec->size();
    const auto start_line = 1 + n_cells * (layer - 1) / 6;

    int line_num = 0;
    int line_num2 = 0;
    while (perm_file) {
        line_num++;
        line_num2++;
        std::string line;
        std::getline(perm_file, line, '\n');

        if (line_num < start_line) continue;

        std::stringstream stream(line);
        for (int i = 0; i < 6; i++) {
            stream >> perm_vec->operator[](cell_id);
            cell_id++;
        }
        if (cell_id == n_cells) break;
    }
    std::ofstream outmat("perm.txt");
    for (int i = 0; i < perm_vec->size(); i++)
    {
        outmat << perm_vec->operator[](i) << std::endl;
    }
    std::cout << "Finished reading permeability data from input file!\n";
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

TPZGeoMesh *CreateSPE10CoarseGeoMesh() {
    std::cout << "Creating SPE10 initial grid...\n";

    const TPZManVector<REAL, 3> x0 = {0, 0, 0};
    const TPZManVector<REAL, 3> x1 = {(REAL)nx, (REAL)ny, (REAL)nz}; // size of domain
    
//    const TPZManVector<int, 3> ndiv = {13, 3, 0};
    const TPZManVector<int, 3> ndiv = {nx, ny, nz}; // 1 unit cell per permeability
//    const TPZManVector<int, 3> ndiv = {1, 1, 0};

    auto gmesh = new TPZGeoMesh;
    
    if (dim == 2){
        TPZGenGrid2D gen(ndiv, x0, x1);
        gen.Read(gmesh);

        gen.SetBC(gmesh, 4, ENoFlux); // bot
        gen.SetBC(gmesh, 5, EPressureRight); // right
        gen.SetBC(gmesh, 6, ENoFlux); // top
        gen.SetBC(gmesh, 7, EPressureLeft); // left
    } else if (dim == 3){
        TPZGenGrid3D gen(x0, x1, ndiv, MMeshType::EHexahedral);
        gmesh = gen.BuildVolumetricElements(1);
        gmesh = gen.BuildBoundaryElements(ENoFlux,ENoFlux,EPressureLeft,ENoFlux,EPressureRight,ENoFlux);
    }
    
    std::cout << "SPE10 initial grid created. NElem: " << gmesh->NElements() << "\n";

    return gmesh;
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

// STATE PermeabilityFunction(const TPZVec<REAL> &x) {
TPZManVector<REAL,3> PermeabilityFunction(const TPZVec<REAL> &x) {
    auto rounded_x = static_cast<int>(x[0]);
    auto rounded_y = static_cast<int>(x[1]);
    if (rounded_x == 220) rounded_x = 219;
    if (rounded_y == 60) rounded_y = 59;
    TPZManVector<REAL,3> perm = {perm_vec[rounded_x * 60 + rounded_y],0,0};

    return perm;
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

// STATE PermeabilityFunction3D(const TPZVec<REAL> &x) {
TPZManVector<REAL,3> PermeabilityFunction3D(const TPZVec<REAL> &x) {
    auto rounded_x = static_cast<int>(x[0]);
    auto rounded_y = static_cast<int>(x[1]);
    auto rounded_z = static_cast<int>(x[2]);
    if (rounded_x == 60) rounded_x = 59;
    if (rounded_y == 220) rounded_y = 219;
    if (rounded_z == 85) rounded_z = 84;
    auto permx = perm_vec[rounded_x + rounded_y * 60 +rounded_z*60*220];
    auto permy = perm_vec[n_cells + rounded_x + rounded_y * 60 +rounded_z*60*220];
    auto permz = perm_vec[n_cells * 2 + rounded_x + rounded_y * 60 +rounded_z*60*220];
    TPZManVector<REAL,3> perm = {permx,permy,permz};

    return perm;
}


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

void PrintResultsVTK(const int dim, TPZLinearAnalysis &an, const std::string &plotfile){
    TPZManVector<std::string,2> scalnames(2), vecnames(1);
    
    // scalnames[0] = "Permeability";
    // scalnames[1] = "Pressure";
    // vecnames[0]= "Flux";
    
    // int div = 0;
    // an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    // an.PostProcess(div,dim);
    constexpr int vtkRes{0};
    TPZVec<std::string> fields = {
        "Permeability", "Flux", "Pressure"
    };
    auto vtk = TPZVTKGenerator(an.Mesh(), fields, plotfile, vtkRes);

    vtk.Do();
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    constexpr int nThreads{12};
    //    TPZFStructMatrix<STATE> matskl(cmesh); // slowest - good for debugging
    //    TPZSkylineStructMatrix<STATE> matskl(cmesh); // medium speed - pz only
    TPZSSpStructMatrix<STATE> matskl(cmesh); // fast - works great with mkl

    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    auto start_time_ass = std::chrono::steady_clock::now();
    cout << "Doing assemble..." << endl;
    an.Assemble();
    auto total_time_ass = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_ass).count()/1000.;
    cout << "Total time assemble = " << total_time_ass << " seconds" << endl;

//  {
//    std::ofstream outmat("mat.nb");
//    TPZMatrixSolver<STATE>* matsol = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
//    matsol->Matrix()->Print("singmat=",outmat,EMathematicaInput);
//    std::ofstream outrhs("rhs.nb");
//      TPZFMatrix<STATE> rhs = an.Rhs();
//    rhs.Print("rhs=",outrhs,EMathematicaInput);
//  }
    
    ///solves the system
    auto start_time_solve = std::chrono::steady_clock::now();
    cout << "\nDoing solve..." << endl;
    an.Solve();
    auto total_time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_solve).count()/1000.;
    cout << "Total time solve = " << total_time_solve << " seconds" << endl;
    
//    {
//        std::ofstream outsol("sol.nb");
//        TPZFMatrix<STATE> sol = an.Solution();
//        sol.Print("sol = ",outsol,EMathematicaInput);
//    }
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
