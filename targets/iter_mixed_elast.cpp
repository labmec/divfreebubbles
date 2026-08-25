#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "TPZAnalyticSolution.h"
#include "TPZCompMeshTools.h"
#include "TPZKernelHdivUtils.h"
#include "pzlog.h"
#include <TPZGeoMeshTools.h>
#include <TPZGmshReader.h>

#include "TPZHDivSHybridApproxCreator.h"
#include "TPZMatRedSolver.h"
#include "TPZRefPattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZSimpleTimer.h"
#include "TPZTimer.h"
#include "TPZVTKGenerator.h"
#include "fstream"
#include "pzbuildmultiphysicsmesh.h"
#include "pzshapecube.h"
#include "pzshapequad.h"
#include "pzshapetetra.h"
#include "pzshapetriang.h"
#include "tpzgeoelrefpattern.h"
// #include "TPZH1HybridApproxCreator.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include <Elasticity/TPZHybridElasticity2D.h>
#include <Elasticity/TPZMixedElasticityND.h>
#include <TPZNullMaterial.h>
#include <TPZNullMaterialCS.h>

#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#ifdef PZ_USING_MUMPS
#include "TPZSSpStructMatrixMumps.h"
#endif
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include <pzstepsolver.h> //for TPZStepSolver
// #include <valgrind/callgrind.h>
#include "TPZTimer.h"
#include <chrono>
#include <sys/resource.h>

double getPeakMemoryMB() {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);

  // ru_maxrss is in KB on Linux
  return usage.ru_maxrss / 1024.0;
}

std::ofstream rprint("results_Mixed_Elastic2D.txt", std::ofstream::out);
std::ofstream printerrors("results_Mixed_Elastic2D_errors.txt", std::ofstream::out);
std::ofstream printmemoryTime("results_Mixed_Elastic2D_memory_time.txt", std::ofstream::app);

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

TElasticity2DAnalytic gElast2d;
int solOrder = 3;

/// @param mfmesh the input multiphysics mesh
void InsertMultiphysicsMaterials(TPZHDivSHybridApproxCreator &create);

enum EMatid { ENone,
              EDomain,
              EBoundary,
              EPont,
              EWrap,
              EIntface,
              EPressureHyb };

/**
   @brief Creates a geometric mesh with elements of a given type on a unit square or cube (depending on the mesh dimension).
   @param[in] meshType element type to be created.
   @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
template <class tshape>
TPZGeoMesh *
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcId);

int main(int argc, char *argv[]) {

  TPZTimer clockall;
  clockall.start();
  auto start = std::chrono::high_resolution_clock::now();

  const int xdiv = 20;
  const int pOrder = (argc > 2) ? std::atoi(argv[2]) : 2;
  const HDivFamily hdivfamily = HDivFamily::EHDivConstant;
  // const HDivFamily hdivfamily = HDivFamily::EHDivStandard;

  const int DIM = 2;

#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  // std::cout << "\nTest Case: \nTopology = " << MElementType_Name(tshape::Type()) <<
  //              ", xdiv = " << xdiv << ", pOrder = " << pOrder <<
  //              ", Approximation space = " << MHDivFamily_Name(hdivfamily) << "\n\n ";

  // std::vector<int> idivs = DIM == 3 ? std::vector<int>{2,8,12,16,32} : std::vector<int>{2,50,100,200,300,400};
  std::vector<int> idivs = DIM == 3 ? std::vector<int>{2, 8, 12, 16, 32} : std::vector<int>{3, 50, 100, 200, 300, 400};
  std::vector<int> pOrders = {1, 2, 3, 4};
  std::cout << "****************** DIM ************** " << DIM << std::endl;
  for (int iorder : pOrders) {
    for (auto idiv : idivs) {
      std::cout << "Running with pOrder = " << iorder << "\n";
      std::cout << "Running with idiv = " << idiv << "\n";
      rprint << "pOrder = " << iorder << " ";

      TPZVec<int> nDivs;

      if (DIM == 2) nDivs = {idiv, idiv};
      if (DIM == 3) nDivs = {idiv, idiv, idiv};

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

      TPZHDivSHybridApproxCreator hdivCreator(gmesh);
      // hdivCreator.HdivFamily() = hdivfamily;
      hdivCreator.SetProbType(ProblemType::EElastic);
      hdivCreator.IsRigidBodySpaces() = false;
      hdivCreator.SetDefaultOrder(iorder);
      hdivCreator.SetShouldCondense(true);
      hdivCreator.SetHybridType(HybridizationType::EStandard);
      InsertMultiphysicsMaterials(hdivCreator);
      // Multiphysics mesh
      TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();

      // Prints gmesh mesh properties
      //  std::string vtk_name = "geoMesh.vtk";
      //  std::ofstream vtkfile(vtk_name.c_str());

      // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

      // std::string txt = "cmesh.txt";
      // std::ofstream myfile(txt);
      // cmesh->Print(myfile);

      // Number of equations without condense elements
      const int nEquationsFull = cmesh->Solution().Rows();
      std::cout << "Number of equations = " << nEquationsFull << std::endl;

      rprint << nEquationsFull << " ";

      TPZTimer clock, clock2;
      clock.start();

      // std::string multFile = "MultiCMesh";
      // util.PrintCompMesh(cmesh,multFile);

      // Number of condensed problem.
      int nEquationsCondensed = cmesh->NEquations();
      std::cout << "Number of equations condensed = " << nEquationsCondensed << std::endl;
      rprint << " Number of equations condensed " << nEquationsCondensed;
      // Create analysis environment
      TPZLinearAnalysis an(cmesh, RenumType::ENone);
      an.SetExact(gElast2d.ExactSolution(), solOrder);

      std::set<int> matBCAll = {EBoundary};
      // Solve problem
      // bool sparse = true;
      bool sparse = (argc > 3) ? std::atoi(argv[3]) : true;

      if (sparse) {
        // if (approxSpace == TPZHDivApproxSpaceCreator<STATE>::EDuplicatedConnects){
        // TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::EDefault);
        // CALLGRIND_START_INSTRUMENTATION;
        // CALLGRIND_TOGGLE_COLLECT;
        TPZMatRedSolver<STATE> solver(an, TPZMatRedSolver<STATE>::EElasticityHDiv);
        printmemoryTime << "iterative porder " << iorder << " div " << idiv << " neqFull " << nEquationsFull << " neqCondensed " << nEquationsCondensed << " ";
        solver.Solve(printmemoryTime);

      } else {
#ifdef PZ_USING_MKL
        TPZSSpStructMatrix<STATE, TPZStructMatrixOR<STATE>> matskl(cmesh);
#endif
#ifdef PZ_USING_MUMPS
        TPZSSpStructMatrixMumps<STATE, TPZStructMatrixOR<STATE>> matskl(cmesh);
#endif
        matskl.SetNumThreads(32);
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
        // auto mat = an.MatrixSolver<REAL>().Matrix();
        // mat->operator*=(-1.0);
        // mat->SetDefPositive(true);
        // TPZFMatrix<REAL>& rhs = an.Rhs();
        // rhs*= -1.0;

        auto startsolver = std::chrono::high_resolution_clock::now();
        an.Solve();
        auto endsolver = std::chrono::high_resolution_clock::now();
        auto durationsolver = std::chrono::duration_cast<std::chrono::milliseconds>(endsolver - startsolver);
        std::cout << "Time solving: " << durationsolver.count() << " ms\n";
        printmemoryTime << "direct " << iorder << " " << idiv << " " << durationassembly.count() << " " << durationsolver.count() << " ";
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

      printmemoryTime << " memory " << getPeakMemoryMB() << "\n";
      printmemoryTime.flush();
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
// Create
template <class tshape>
TPZGeoMesh *
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, EMatid bcId) {

  MMeshType meshType;
  int dim = tshape::Dimension;

  switch (tshape::Type()) {
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

  TPZManVector<REAL, 3> minX = {0, 0, 0};
  TPZManVector<REAL, 3> maxX = {1, 1, 1};
  int nMats = 2 * dim + 1;

  // all bcs share the same id
  constexpr bool createBoundEls{true};
  TPZVec<int> matIds(nMats, bcId);
  matIds[0] = volId;
  // matIds[1] = bcId;
  // matIds[2] = EBoundary1;
  // matIds[3] = EBoundary1;
  // matIds[4] = EBoundary1;

  TPZGeoMesh *gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX, matIds, nDivs, meshType, createBoundEls);
  // TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(meshType,
  //                     volId,createBoundEls, bcId);

  return gmesh;
}

/// @brief Insert the hybrid elastic material and boundary material
/// @param mfmesh the input multiphysics mesh
void InsertMultiphysicsMaterials(TPZHDivSHybridApproxCreator &create) {
  gElast2d.fProblemType = TElasticity2DAnalytic::EHomogeneous;

  REAL E = 1.;
  REAL nu = 0.;
  gElast2d.gE = E;
  gElast2d.gPoisson = nu;
  int volmatid = EDomain;
  int dim = create.GeoMesh()->Dimension();
  auto *mat = new TPZMixedElasticityND(volmatid, E, nu, 1., 1., gElast2d.fPlaneStress, dim);
  mat->SetForcingFunction(gElast2d.ForceFunc(), 2);
  mat->SetExactSol(gElast2d.ExactSolution(), 2);

  create.InsertMaterialObject(mat);
  TPZFNMatrix<4, REAL> val1(2, 2, 0.);
  TPZManVector<REAL, 2> val2(2, 1.);
  auto *bndbt = mat->CreateBC(mat, EBoundary, 0, val1, val2);
  bndbt->SetForcingFunctionBC(gElast2d.ExactSolution(), 2);

  create.InsertMaterialObject(bndbt);
}
