#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <TPZGeoMeshTools.h>
#include "TPZKernelHdivUtils.h"
#include "TPZAnalyticSolution.h"
#include <TPZGmshReader.h>
#include "TPZCompMeshTools.h"
#include "pzlog.h"

#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
#include "TPZTimer.h"
#include "TPZMatRedSolver.h"
#include "fstream"
#include "TPZSimpleTimer.h"
#include "TPZVTKGenerator.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZHDivApproxCreator.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
//#include <valgrind/callgrind.h>
#include "TPZTimer.h"
#include <chrono>
#include <sys/resource.h>


double getPeakMemoryMB() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    // ru_maxrss is in KB on Linux
    return usage.ru_maxrss / 1024.0;
}


std::ofstream rprint("results_Harmonic2D.txt",std::ofstream::out);
std::ofstream printerrors("results_errors.txt",std::ofstream::out);
std::ofstream printmemoryTime("results_memory_time.txt",std::ofstream::app);

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

enum EMatid  {ENone, EDomain, EBoundary, EPont, EWrap, EIntface, EPressureHyb};

/**
   @brief Creates a geometric mesh with elements of a given type on a unit square or cube (depending on the mesh dimension).
   @param[in] meshType element type to be created.
   @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
template<class tshape>
TPZGeoMesh*
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcId);

//Analytical solution
constexpr int solOrder{4};
auto exactSol = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    // const auto &d = 1.; // distanc between injection and production wells
    // u[0]= x*x-y*y ;
    // gradU(0,0) = -2*x;
    // gradU(1,0) = 2.*y;
    // gradU(2,0) = 0.;

    // u[0] =  5. + 3. * x + 2. * y + 4. * x * y;
    // gradU(0,0) = 3. + 4. * y;
    // gradU(1,0) = 2. + 4. * x;
    // gradU(2,0) = 0.;
    
    // u[0] = x*x*x*y - y*y*y*x;
    // gradU(0,0) = (3.*x*x*y - y*y*y);
    // gradU(1,0) = (x*x*x - 3.*y*y*x);

    // u[0]= x*x - y*y;
    // gradU(0,0) = 2.*x;
    // gradU(1,0) = -2.*y;
    // // gradU(2,0) = -1;

    REAL aux = 1./sinh(sqrt(2)*M_PI);
    u[0] = sin(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2)*M_PI*z)*aux;
    gradU(0,0) = M_PI*cos(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2)*M_PI*z)*aux;
    gradU(1,0) = M_PI*cos(M_PI*y)*sin(M_PI*x)*sinh(sqrt(2)*M_PI*z)*aux;
    gradU(2,0) = sqrt(2)*M_PI*cosh(sqrt(2)*M_PI*z)*sin(M_PI*x)*sin(M_PI*y)*aux;

    // u[0]= std::sin(M_PI*x)*std::sin(M_PI*y);
    // gradU(0,0) = M_PI*cos(M_PI*x)*sin(M_PI*y);
    // gradU(1,0) = M_PI*cos(M_PI*y)*sin(M_PI*x);

    // u[0]=pow(2,2 - pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x);
    // gradU(0,0) = 0.;
    // gradU(1,0) = 0.;

    // u[0] = (x-1)*x*(y-1)*y*(z-1)*z;
    // gradU(0,0) = (x-1)*(y-1)*y*(z-1)*z + x*(y-1)*y*(z-1)*z;
    // gradU(1,0) = (x-1)*x*(y-1)*(z-1)*z + (x-1)*x*y*(z-1)*z;
    // gradU(1,0) = (x-1)*x*(y-1)*y*(z-1) + (x-1)*x*(y-1)*y*z;

    // REAL a1 = 1./4;
    // REAL alpha = M_PI/2;
    // u[0] = x*a1*cos(x*alpha)*cosh(y*alpha) + y*a1*sin(x*alpha)*sinh(y*alpha) + x*x - y*y;
    // gradU(0,0) = -a1*(cosh(alpha*y)*(cos(alpha*x) - alpha*x*sin(alpha*x)) + alpha*y*cos(alpha*x)*sinh(alpha*y));
    // gradU(1,0) = -a1*(alpha*y*cosh(alpha*y)*sin(alpha*x) + (alpha*x*cos(alpha*x) + sin(alpha*x))*sinh(alpha*y));

    // u[0] = exp(M_PI*x)*sin(M_PI*y);
    // gradU(0,0) = M_PI*exp(M_PI*x)*sin(M_PI*y);
    // gradU(1,0) = M_PI*exp(M_PI*x)*cos(M_PI*y);

    // u[0] = 0.5*(x)*x+0.5*(y)*y-(z)*z;
    // gradU(0,0) = -x;//(x-1)*(y-1)*y*(z-1)*z + x*(y-1)*y*(z-1)*z;
    // gradU(1,0) = -y;//(x-1)*x*(y-1)*(z-1)*z + (x-1)*x*y*(z-1)*z;
    // gradU(2,0) = 2.*z;//(x-1)*x*(y-1)*y*(z-1) + (x-1)*x*(y-1)*y*z;

    // u[0] = exp(M_PI*x)*sin(M_PI*y)*x;
    // gradU(0,0) = x*M_PI*exp(M_PI*x)*sin(M_PI*y)+exp(M_PI*x)*sin(M_PI*y);
    // gradU(1,0) = x*M_PI*exp(M_PI*x)*cos(M_PI*y);
};
auto forcefunction = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    // u[0] = -exp(M_PI*x)*sin(M_PI*y)*M_PI*2.;
    REAL aux = 1./sinh(sqrt(2)*M_PI);
    u[0] = -2.*sqrt(2)*M_PI*M_PI*cos(M_PI*x)*sin(M_PI*y)*cosh(sqrt(2)*M_PI*x)*aux;
};

int main(int argc, char* argv[])
{
    
    TPZTimer clockall;
    clockall.start();
    auto start = std::chrono::high_resolution_clock::now();

    const int xdiv = 20;
    const int pOrder = 1;
    const HDivFamily hdivfamily = HDivFamily::EHDivConstant;
    // const HDivFamily hdivfamily = HDivFamily::EHDivStandard;
    int DIM = 3;

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    // std::cout << "\nTest Case: \nTopology = " << MElementType_Name(tshape::Type()) << 
    //              ", xdiv = " << xdiv << ", pOrder = " << pOrder << 
    //              ", Approximation space = " << MHDivFamily_Name(hdivfamily) << "\n\n "; 
    
    std::vector<int> idivs = {2};
    // std::vector<int> idivs = {10,50,100,200};
    
    for (int iorder = 3; iorder < 4; iorder++) {
    for (auto idiv : idivs) {
        std::cout << "Running with pOrder = " << iorder << "\n";
        std::cout << "Running with idiv = " << idiv << "\n";
        rprint << "pOrder = " << iorder << " " ;

    TPZVec<int> nDivs;

    if (DIM == 2) nDivs = {idiv,idiv};
    if (DIM == 3) nDivs = {idiv,idiv,idiv};
    
    

    // Creates/import a geometric mesh 
    TPZGeoMesh *gmesh = nullptr; 
    if (DIM == 2) {
        // gmesh = CreateGeoMesh<pzshape::TPZShapeTriang>(nDivs, EDomain, EBoundary);
        gmesh = CreateGeoMesh<pzshape::TPZShapeQuad>(nDivs, EDomain, EBoundary);
    }
    if (DIM == 3) {
        gmesh = CreateGeoMesh<pzshape::TPZShapeCube>(nDivs, EDomain, EBoundary);
    }
    

    // Util for HDivKernel printing and solving
    TPZKernelHdivUtils<STATE> util;

    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.HdivFamily() = hdivfamily;
    hdivCreator.SetProbType(ProblemType::EDarcy);
    hdivCreator.IsRigidBodySpaces() = false;
    hdivCreator.SetDefaultOrder(iorder);
    hdivCreator.SetExtraInternalOrder(0);
    hdivCreator.SetShouldCondense(true);
    // hdivCreator.SetShouldCondense(false);
    hdivCreator.SetHybridType(HybridizationType::EStandard);

    //Prints gmesh mesh properties
    // std::string vtk_name = "geoMesh.vtk";
    // std::ofstream vtkfile(vtk_name.c_str());

    // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

    //Insert Materials
    TPZMixedDarcyFlow* matdarcy = new TPZMixedDarcyFlow(EDomain,DIM);
    matdarcy->SetConstantPermeability(1.);
    matdarcy->SetExactSol(exactSol,4);
    matdarcy->SetForcingFunction(forcefunction,4);

    hdivCreator.InsertMaterialObject(matdarcy);

    TPZFMatrix<STATE> val1(3,3,0.);
    TPZManVector<STATE> val2(3,0.);
    TPZBndCondT<STATE> *BCond1 = matdarcy->CreateBC(matdarcy, EBoundary, 0, val1, val2);
    BCond1->SetForcingFunctionBC(exactSol,4);
    hdivCreator.InsertMaterialObject(BCond1);

    //Multiphysics mesh
    TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
    // std::string txt = "cmesh.txt";
    // std::ofstream myfile(txt);
    // cmesh->Print(myfile);

  
    // Number of equations without condense elements
    const int nEquationsFull = cmesh->NEquations();
    std::cout << "Number of equations = " << nEquationsFull << std::endl;

    rprint << nEquationsFull << " " ;

    TPZTimer clock,clock2;
    clock.start();


    // std::string multFile = "MultiCMesh";
    // util.PrintCompMesh(cmesh,multFile);
    
    //Number of condensed problem.
    int nEquationsCondensed = cmesh->NEquations();
    std::cout << "Number of equations condensed = " << nEquationsCondensed << std::endl;
    //Create analysis environment
    TPZLinearAnalysis an(cmesh,RenumType::EMetis);
    an.SetExact(exactSol,solOrder);

    std::set<int> matBCAll = {EBoundary};
    // Solve problem
    // bool sparse = true;
    bool sparse = false;
    
    if (sparse){
    // if (approxSpace == TPZHDivApproxSpaceCreator<STATE>::EDuplicatedConnects){
        // TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::EDefault);
        // CALLGRIND_START_INSTRUMENTATION;
        // CALLGRIND_TOGGLE_COLLECT;
        TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::ESparse);
        
        solver.Solve(printmemoryTime);
        

    } else {
        TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matskl(cmesh);
        matskl.SetNumThreads(100);
        an.SetStructuralMatrix(matskl);

        // TPZBlockDiagonalStructMatrix<STATE> BDFmatrix(cmesh);
        // TPZBlockDiagonal<REAL> KBD;
        // std::cout << "Start assembling BlockDiag ...\n";
        // BDFmatrix.AssembleBlockDiagonal(KBD);
        // std::cout << "Finish assembling BlockDiag ...\n";

        // //Creates the preconditioner 
        // TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>( &KBD );
        // precond->SetDirect(ELU);
        // int64_t nMaxIter = 5000;
        // TPZVec<REAL> errors(nMaxIter);
        // errors.Fill(0.);

        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        // // step.SetCG(nMaxIter,*precond,1.e-10,0);
        // step.SetGMRES(nMaxIter,100,*precond,1.e-10,0);
        an.SetSolver(step);
        auto startassembly = std::chrono::high_resolution_clock::now();
        an.Assemble();
        auto endassembly = std::chrono::high_resolution_clock::now();
        auto durationassembly = std::chrono::duration_cast<std::chrono::milliseconds>(endassembly - startassembly);
        std::cout << "Time assembling: " << durationassembly.count() << " ms\n";
        auto mat = an.MatrixSolver<REAL>().Matrix();
        mat->operator*=(-1.0);
        mat->SetDefPositive(true);
        TPZFMatrix<REAL>& rhs = an.Rhs();
        rhs*= -1.0;
        
        auto startsolver = std::chrono::high_resolution_clock::now();
        an.Solve();
        auto endsolver = std::chrono::high_resolution_clock::now();
        auto durationsolver = std::chrono::duration_cast<std::chrono::milliseconds>(endsolver - startsolver);
        std::cout << "Time solving: " << durationsolver.count() << " ms\n";
        printmemoryTime << iorder << " " << idiv << " " << durationassembly.count() << " " << durationsolver.count() << " " << getPeakMemoryMB() << "\n";
    }

    

    clock.stop();
    clockall.stop();
    // std::cout << "Time running = " << clock.seconds() << std::endl;
    // std::cout << "Time runningall = " << clockall.seconds() << std::endl;

    auto end = std::chrono::high_resolution_clock::now();

    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Elapsed time: "
              << duration.count()
              << " ms\n";
    std::cout << "Memory usage: "
              << getPeakMemoryMB()
              << " MB\n";


    printmemoryTime << iorder << " " << idiv << " " << duration.count() << " " << getPeakMemoryMB() << "\n";
    // //Print results
    // {
    //     TPZSimpleTimer postProc("Post processing1");
        // util.PrintResultsMultiphysics(cmesh->MeshVector(),an,cmesh);
    // }
    }
    }
    {
        // TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), cmesh);
        // TPZSimpleTimer postProc("Post processing2");
        // const std::string plotfile = "myfile";//sem o .vtk no final
        // constexpr int vtkRes{0};
    

        // TPZVec<std::string> fields = {
        // "Pressure",
        // "ExactPressure",
        // "Flux",
        // "ExactFlux"};
        // auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

        // vtk.Do();
    }
    // std::string txt2 = "cmeshSol.txt";
    // std::ofstream myfile2(txt2);
    // cmesh->Print(myfile2);

    // //vamos supor que vc atualiza a solucao, roda de novo, sei la
    // vtk.Do();

    // //Compute error
    // std::ofstream anPostProcessFile("postprocess.txt");
    // TPZManVector<REAL,5> error;
    // int64_t nelem = cmesh->NElements();
    // cmesh->LoadSolution(cmesh->Solution());
    // cmesh->ExpandSolution();
    // cmesh->ElementSolution().Redim(nelem, 5);
    // an.PostProcessError(error,false,anPostProcessFile);
    
    // printerrors << xdiv << std::scientific << std::setprecision(8) << " " << error[0] << " " 
    //  << error[1] << " " << error[2] << " "  << error[3] << " "  << error[4] << std::endl;

    // //Check error
    // // REAL tolerance = 1.e-6;
    // std::cout << "ERROR[0] = " << std::scientific << std::setprecision(15) << error[0] << std::endl;
    // std::cout << "ERROR[1] = " << error[1] << std::endl;
    // std::cout << "ERROR[2] = " << error[2] << std::endl;
    // std::cout << "ERROR[3] = " << error[3] << std::endl;
    // std::cout << "ERROR[4] = " << error[4] << std::endl;
    // // // REQUIRE(error[1] < tolerance);



    return 0;
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
    // matIds[1] = bcId;
    // matIds[2] = EBoundary1;
    // matIds[3] = EBoundary1;
    // matIds[4] = EBoundary1;
    
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls);
    // TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(meshType,
    //                     volId,createBoundEls, bcId);
    
    return gmesh;
    
}
