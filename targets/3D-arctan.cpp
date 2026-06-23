/*
  This unit test verifies if the hybridization and semi hybridization techniques are working
  for any specified polynomial order and topology.

*/
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
#include <pzstepsolver.h>       //for TPZStepSolver
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
// #include <valgrind/callgrind.h>

// ----- Run tests with or without main -----
#define RUNWITHMAIN

#ifndef RUNWITHMAIN
#include <catch2/catch.hpp>
#endif

std::ofstream rprint("results.txt", std::ios_base::app);
std::ofstream printerrors("results_errors.txt", std::ios_base::app);

/** @brief Returns the name of the HDiv Family approximation space. */
inline std::string MHDivFamily_Name(HDivFamily hdivfam)
{
    switch (hdivfam)
    {
    case HDivFamily::EHDivStandard:
    {
        return "EHDivStandard";
    }
    case HDivFamily::EHDivConstant:
    {
        return "EHDivConstant";
    }
    case HDivFamily::EHDivKernel:
    {
        return "EHDivKernel";
    }
    case HDivFamily::EHDivOptimized:
    {
        return "EHDivOptimized";
    }
    default:
    {
        return "HDivFamily not found!";
    }
    }
    DebugStop();
    return "";
}

enum EMatid
{
    ENone,
    EDomain,
    EBoundary,
    EPont,
    EWrap,
    EIntface,
    EPressureHyb
};

constexpr bool printVTK{true};
constexpr bool computeError{false};

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

/**
   @brief Reads the test mesh from gmsh
   @param[in] file_name the .msh mesh file.
*/
template <class tshape>
TPZGeoMesh *
ReadMeshFromGmsh(std::string file_name);

// The test function
template <class tshape>
void TestHybridization(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily);

int main(int argc, char *argv[])
{
    int xdiv = 2;
    if (argc > 1)
    {
        xdiv = std::stoi(argv[1]);
    }
    int pOrder = 1;
    if (argc > 2)
    {
        pOrder = std::stoi(argv[2]);
    }
    HDivFamily hdivfam = HDivFamily::EHDivOptimized;
    TestHybridization<pzshape::TPZShapeCube>(xdiv, pOrder, hdivfam);
    return 0;
}

// Analytical solution
constexpr int solOrder{4};
auto exactSol = [](const TPZVec<REAL> &loc,
                   TPZVec<STATE> &u,
                   TPZFMatrix<STATE> &gradU)
{
    const auto &x = loc[0];
    const auto &y = loc[1];
    const auto &z = loc[2];

    // 3D arctan problem from Sonia's paper. Available at https://onlinelibrary.wiley.com/doi/10.1002/nme.6337
    REAL a = sqrt(pow(-1.25 + x, 2) + pow(0.25 + y, 2) + pow(0.25 + z, 2));
    u[0] = M_PI / 2.0 - atan(5.0 * (-M_PI / 3.0 + a));
    gradU(0, 0) = (-5.0 * (-1.25 + x)) / (a * (1.0 + 25.0 * pow(-M_PI / 3.0 + a, 2)));
    gradU(1, 0) = (-5.0 * (0.25 + y)) / (a * (1.0 + 25.0 * pow(-M_PI / 3.0 + a, 2)));
    gradU(2, 0) = (-5.0 * (0.25 + z)) / (a * (1.0 + 25.0 * pow(-M_PI / 3.0 + a, 2)));
};

auto forcingFunc = [](const TPZVec<REAL> &loc,
                      TPZVec<STATE> &force)
{
    const auto &x = loc[0];
    const auto &y = loc[1];
    const auto &z = loc[2];

    // Steep wave solution
    REAL b = pow(-1.25 + x, 2) + pow(0.25 + y, 2) + pow(0.25 + z, 2);
    REAL a = sqrt(b);
    force[0] = (250 * pow(-1.25 + x, 2) * (-0.3333333333333333 * M_PI + a)) /
                   (b * pow(1 + 25 * pow(-0.3333333333333333 * M_PI + a, 2), 2)) +
               (250 * pow(0.25 + y, 2) * (-0.3333333333333333 * M_PI + a)) /
                   (b * pow(1 + 25 * pow(-0.3333333333333333 * M_PI + a, 2), 2)) +
               (250 * pow(0.25 + z, 2) * (-0.3333333333333333 * M_PI + a)) /
                   (b * pow(1 + 25 * pow(-0.3333333333333333 * M_PI + a, 2), 2)) +
               (5 * pow(-1.25 + x, 2)) / (pow(b, 1.5) *
                                          (1 + 25 * pow(-0.3333333333333333 * M_PI + a, 2))) +
               (5 * pow(0.25 + y, 2)) / (pow(b, 1.5) *
                                         (1 + 25 * pow(-0.3333333333333333 * M_PI + a, 2))) +
               (5 * pow(0.25 + z, 2)) / (pow(b, 1.5) *
                                         (1 + 25 * pow(-0.3333333333333333 * M_PI + a, 2))) -
               15 / (a * (1 + 25 * pow(-0.3333333333333333 * M_PI + a, 2)));
    force[0] *= -1; // because flux is defined as -K*gradU, and K=1 in this case
};

template <class tshape>
void TestHybridization(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily)
{

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    std::cout << "\nTest Case: \nTopology = " << MElementType_Name(tshape::Type()) << ", xdiv = " << xdiv << ", pOrder = " << pOrder << ", Approximation space = " << MHDivFamily_Name(hdivfamily) << "\n\n ";

    int DIM = tshape::Dimension;
    TPZVec<int> nDivs;

    if (DIM == 2)
        nDivs = {xdiv, xdiv};
    if (DIM == 3)
        nDivs = {xdiv, xdiv, xdiv};

    // Creates/import a geometric mesh
    auto gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary);

    std::set<int> fBCMatId = {EBoundary};
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->Dimension() < gmesh->Dimension())
            continue;

        int nSides = gel->NSides();
        // For tetrahedra only, loop over the surface sides
        for (int side = 0; side < nSides; side++)
        {
            if (gel->SideDimension(side) != gel->Dimension() - 1)
                continue;

            TPZGeoElSide gelside(gel, side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            // Neighbour material id
            auto Nmatid = neighbour.Element()->MaterialId();

            /*  If the boundary has BC, delete the neighbour GeoElement and
                create another one from TPZGeoElBC with the same material id
            */
            if (fBCMatId.find(Nmatid) == fBCMatId.end())
                continue;
            gmesh->DeleteElement(neighbour.Element(), neighbour.Element()->Index());
            TPZGeoElBC gelbcWrap(gelside, Nmatid);
        }
    }

    // int dim = gmesh->Dimension();
    // TPZManVector<TPZGeoEl*,10> children;
    // int64_t nel = gmesh->NElements();
    // // for (int i = 0; i < nel; i++)
    // // {
    //     // if (gmesh->ElementVec()[0]->Dimension()==dim)
    //     gmesh->ElementVec()[4]->Divide(children);

    //     for(int64_t el = 0; el<nel; el++) {
    //         TPZGeoEl *gel = gmesh->Element(el);
    //         if(gel->Dimension() != dim-1) continue;
    //         if(gel->HasSubElement()) continue;
    //         TPZGeoElSide gelside(gel);
    //         TPZGeoElSide neighbour = gelside.Neighbour();
    //         if(neighbour.HasSubElement()) {
    //             TPZManVector<TPZGeoEl*,10> children2;
    //             gel->Divide(children2);
    //         }
    //     }
    // // }

    // Util for HDivKernel printing and solving
    TPZKernelHdivUtils<STATE> util;

    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.HdivFamily() = hdivfamily;
    hdivCreator.ProbType() = ProblemType::EDarcy;
    hdivCreator.IsRigidBodySpaces() = false;
    hdivCreator.SetDefaultOrder(pOrder);
    hdivCreator.SetExtraInternalOrder(0);
    hdivCreator.SetShouldCondense(true);
    // hdivCreator.SetShouldCondense(false);
    hdivCreator.HybridType() = HybridizationType::ESemi;
    // hdivCreator.HybridType() = HybridizationType::EStandard;

    // Prints gmesh mesh properties
    // std::string vtk_name = "geoMesh.vtk";
    // std::ofstream vtkfile(vtk_name.c_str());
    // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

    // Insert Materials
    TPZMixedDarcyFlow *matdarcy = new TPZMixedDarcyFlow(EDomain, DIM);
    matdarcy->SetConstantPermeability(1.);
    matdarcy->SetExactSol(exactSol, 4);
    matdarcy->SetForcingFunction(forcingFunc, 4);

    hdivCreator.InsertMaterialObject(matdarcy);

    TPZFMatrix<STATE> val1(1, 1, 0.);
    TPZManVector<STATE> val2(1, 0.);
    TPZBndCondT<STATE> *BCond1 = matdarcy->CreateBC(matdarcy, EBoundary, 0, val1, val2);
    BCond1->SetForcingFunctionBC(exactSol, 4);
    hdivCreator.InsertMaterialObject(BCond1);

    // Multiphysics mesh
    TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
    // std::string txt = "cmesh.txt";
    // std::ofstream myfile(txt);
    // cmesh->Print(myfile);

    if (printVTK)
    {
        TPZSimpleTimer postProc("Post processing2");
        const std::string plotfile = "myfile"; // sem o .vtk no final
        constexpr int vtkRes{5};

        TPZVec<std::string> fields = {
            "ExactPressure",
            "ExactFlux",
            "ExactDivSigma"};
        auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes, cmesh->Dimension());

        vtk.Do();
        return;
    }

    // Number of equations without condense elements
    int64_t nEquationsFull = 0;
    int64_t ncon = cmesh->NConnects();
    for (int64_t i = 0; i < ncon; i++)
    {
        TPZConnect &df = cmesh->ConnectVec()[i];
        if (df.HasDependency() || !df.NElConnected() || df.SequenceNumber() == -1)
        {
            continue;
        }

        int dofsize = df.NShape() * df.NState();
        nEquationsFull += dofsize;
    }

    std::cout << "Number of equations = " << nEquationsFull << std::endl;

    rprint << xdiv << " " << nEquationsFull << " ";

    TPZTimer clock, clock2;
    clock.start();

    // std::string multFile = "MultiCMesh";
    // util.PrintCompMesh(cmesh,multFile);

    // Number of condensed problem.
    int nEquationsCondensed = cmesh->NEquations();
    std::cout << "Number of equations condensed = " << nEquationsCondensed << std::endl;
    // Create analysis environment
    TPZLinearAnalysis an(cmesh, RenumType::EMetis);
    an.SetExact(exactSol, solOrder);

    std::set<int> matBCAll = {EBoundary};
    // Solve problem
    bool sparse = true;
    //    bool sparse = false;

    if (sparse)
    {
        // if (approxSpace == TPZHDivApproxSpaceCreator<STATE>::EDuplicatedConnects){
        // TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::EDefault);
        // CALLGRIND_START_INSTRUMENTATION;
        // CALLGRIND_TOGGLE_COLLECT;
        TPZMatRedSolver<STATE> solver(an, matBCAll, TPZMatRedSolver<STATE>::ESparse);
        clock2.start();
        solver.Solve(rprint);
        clock2.stop();
        // CALLGRIND_TOGGLE_COLLECT;
        // CALLGRIND_STOP_INSTRUMENTATION;
        // std::cout << "Time SOLVER = " << clock2 << std::endl;

        // bool filter = false;
        // if (DIM == 3 && hdivfamily == HDivFamily::EHDivKernel) filter = true;
        // createSpace.Solve(an, cmesh, true, filter);
    }
    else
    {
        bool domHyb = false;
        util.SolveProblemDirect(an, cmesh, false, domHyb);
    }
    clock.stop();

    if (printVTK || computeError)
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), cmesh);

    if (printVTK)
    {
        TPZSimpleTimer postProc("Post processing2");
        const std::string plotfile = "myfile"; // sem o .vtk no final
        constexpr int vtkRes{1};

        TPZVec<std::string> fields = {
            "Pressure",
            "ExactPressure",
            "Flux",
            "ExactFlux",
            "ExactDivSigma"};
        auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes, cmesh->Dimension());

        vtk.Do();
    }

    // Compute error
    if (computeError)
    {
        TPZManVector<REAL, 5> error;
        error.Resize(5);
        int64_t nelem = cmesh->NElements();
        cmesh->LoadSolution(cmesh->Solution());
        cmesh->ExpandSolution();
        cmesh->ElementSolution().Redim(nelem, 5);
        // an.SetThreadsForError(12);
        an.PostProcessError(error, false);

        printerrors << xdiv << std::scientific << std::setprecision(8) << " " << error[0] << " "
                    << error[1] << " " << error[2] << " " << error[3] << " " << error[4] << std::endl;

        // Check error
        //  REAL tolerance = 1.e-6;
        std::cout << "ERROR[0] = " << std::scientific << std::setprecision(15) << error[0] << std::endl;
        std::cout << "ERROR[1] = " << error[1] << std::endl;
        std::cout << "ERROR[2] = " << error[2] << std::endl;
        std::cout << "ERROR[3] = " << error[3] << std::endl;
        std::cout << "ERROR[4] = " << error[4] << std::endl;
    }
}

// Create
template <class tshape>
TPZGeoMesh *
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

    TPZGeoMesh *gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                                                             matIds, nDivs, meshType, createBoundEls);
    // TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(meshType,
    //                     volId,createBoundEls, bcId);

    return gmesh;
}

template <class tshape>
TPZGeoMesh *
ReadMeshFromGmsh(std::string file_name)
{
    // read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string, int>, 4> stringtoint(4);
        stringtoint[3]["Domain"] = 1;
        stringtoint[2]["Surfaces"] = 2;

        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name, gmesh);
    }

    return gmesh;
}
