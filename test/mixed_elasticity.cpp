/*
  This unit test verifies if the hybridization and semi hybridization techniques are working
  for any specified polynomial order and topology.
  
*/
#include <catch2/catch.hpp>
#include <TPZGeoMeshTools.h>
#include "TPZHDivApproxSpaceCreator.h"
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
#include "pzbuildmultiphysicsmesh.h"
#include "TPZAnalyticSolution.h"
#include "TPZHDivApproxCreator.h"
#include "TPZMixedElasticityND.h"

std::ofstream rprint;


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
		default:
        {
            return "HDivFamily not found!";
        }
    }
    DebugStop();
	return "";
}

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

/**
   @brief Reads the test mesh from gmsh
   @param[in] file_name the .msh mesh file.
*/
template<class tshape>
TPZGeoMesh*
ReadMeshFromGmsh(std::string file_name);


// The test function
template<class tshape>
void TestHybridization(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily);

TEST_CASE("Hybridization test")
{
    rprint.open("results_MElasticity2D.txt",std::ios_base::app);
    // const int pOrder = 1;
    const int pOrder = GENERATE(1);

    const int xdiv = 10;//GENERATE(50);
    // const int xdiv = GENERATE(2,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200);
    // const int xdiv = GENERATE(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
    // const int xdiv = GENERATE(2,3,4,5,6,7,8);
    // HDivFamily hdivfam = GENERATE(HDivFamily::EHDivConstant,HDivFamily::EHDivKernel);
    // HDivFamily hdivfam = GENERATE(HDivFamily::EHDivKernel);
    HDivFamily hdivfam = GENERATE(HDivFamily::EHDivConstant);
    // HDivFamily hdivfam = GENERATE(HDivFamily::EHDivStandard);
    // HDivFamily hdivfam = GENERATE(HDivFamily::EHDivStandard,HDivFamily::EHDivConstant);
    
    // TestHybridization<pzshape::TPZShapeTriang>(xdiv,pOrder,hdivfam);
    TestHybridization<pzshape::TPZShapeQuad>(xdiv,pOrder,hdivfam); 
    // TestHybridization<pzshape::TPZShapeTetra>(xdiv,pOrder,hdivfam); 
    // TestHybridization<pzshape::TPZShapeCube>(xdiv,pOrder,hdivfam);
}


template<class tshape>
void TestHybridization(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily)
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    std::cout << "\nTest Case: \nTopology = " << MElementType_Name(tshape::Type()) << 
                 ", xdiv = " << xdiv << ", pOrder = " << pOrder << 
                 ", Approximation space = " << MHDivFamily_Name(hdivfamily) << "\n\n "; 
    
    int DIM = tshape::Dimension;
    TPZVec<int> nDivs;

    if (DIM == 2) nDivs = {xdiv,xdiv};
    if (DIM == 3) nDivs = {xdiv,xdiv,xdiv};
    
    // Creates/import a geometric mesh
    auto gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary);
    // auto gmesh = ReadMeshFromGmsh<tshape>("../mesh/1tetra.msh");
    
    // Util for HDivKernel printing and solving
    TPZKernelHdivUtils<STATE> util;
    
    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.HdivFamily() = hdivfamily;
    hdivCreator.ProbType() = ProblemType::EElastic;
    hdivCreator.IsRigidBodySpaces() = true;
    hdivCreator.SetDefaultOrder(pOrder);
    hdivCreator.SetExtraInternalOrder(0);
    hdivCreator.SetShouldCondense(true);
    hdivCreator.HybridType() = HybridizationType::ENone;

    TPZAnalyticSolution *gAnalytic = 0;
    TPZMixedElasticityND* matelastic = 0;
    if(DIM == 2)
    {
        TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
        elas->gE = 1.e3;
        elas->gPoisson = 0.3;
        elas->fProblemType = TElasticity2DAnalytic::EThiago;
        elas->fPlaneStress = 0;
        gAnalytic = elas;
        matelastic = new TPZMixedElasticityND(EDomain, elas->gE, elas->gPoisson, 0, 0, elas->fPlaneStress, DIM);
        
    }
    else if(DIM == 3)
    {
        TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
        elas->fE = 1.;//206.8150271873455;
        elas->fPoisson = 0.0;//0.3040039545229857;
        elas->fProblemType = TElasticity3DAnalytic::EDispx;
        gAnalytic = elas;
    }


    //Insert Materials
    matelastic->SetExactSol(gAnalytic->ExactSolution(),4);
    

    hdivCreator.InsertMaterialObject(matelastic);

    TPZFMatrix<STATE> val1(DIM,DIM,0.);
    TPZManVector<STATE> val2(DIM,0.);
    TPZBndCondT<STATE> *BCond1 = matelastic->CreateBC(matelastic, EBoundary, 0, val1, val2);
    BCond1->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
    hdivCreator.InsertMaterialObject(BCond1);

    //Multiphysics mesh
    TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
    std::string txt = "cmesh.txt";
    std::ofstream myfile(txt);
    cmesh->Print(myfile);

  
    // Number of equations without condense elements
    const int nEquationsFull = cmesh->NEquations();
    std::cout << "Number of equations = " << nEquationsFull << std::endl;

    //Create analysis environment
    TPZLinearAnalysis an(cmesh,true);
    an.SetExact(gAnalytic->ExactSolution());

    std::set<int> matBCAll = {EBoundary};
    //Solve problem
    // if (approxSpace == TPZHDivApproxSpaceCreator<STATE>::EDuplicatedConnects){
        // TPZMatRedSolver<STATE> solver(&an,matBCAll,TPZMatRedSolver<STATE>::EDefault);
        TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::ESparse);
        solver.Solve(rprint);
        // std::cout << "Time SOLVER = " << clock2 << std::endl;

        // bool filter = false;
        // if (DIM == 3 && hdivfamily == HDivFamily::EHDivKernel) filter = true;
        // createSpace.Solve(an, cmesh, true, filter);

    // } else {
    //     //Equation filter (spanning trees), true if 3D and HDivKernel 
    //     bool filter = false;bool domainhybr=false;
    //     if (DIM == 3 && hdivfamily == HDivFamily::EHDivKernel) filter = true;
    // //     // createSpace.Solve(an, cmesh, true, filter);
    //     util.SolveProblemDirect(an,cmesh,filter,domainhybr);
    // }

    // std::cout << "Time running = " << clock << std::endl;

    // //Print results
    // {
    //     TPZSimpleTimer postProc("Post processing1");
    //     util.PrintResultsMultiphysics(meshvector,an,cmesh);
    // }

    // {
        
    //     TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector_HDiv, cmesh_m_HDiv);
    //     TPZSimpleTimer postProc("Post processing2");
    //     const std::string plotfile = "myfile";//sem o .vtk no final
    //     constexpr int vtkRes{0};
    

    //     TPZVec<std::string> fields = {
    //     // "ExactDisplacement",
    //     // "ExactStress",
    //     "Displacement",
    //     "SigmaX",
    //     "SigmaY",
    //     "TauXY"
    //     };
    //     auto vtk = TPZVTKGenerator(cmesh_m_HDiv, fields, plotfile, vtkRes);

    //     vtk.Do();
    //     // cmesh_m_HDiv->Solution().Print("Solution=",std::cout);
        
    // }
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
    
    // //Check error
    // REAL tolerance = 1.e-6;
    // std::cout << "ERROR[0] = " << std::scientific << std::setprecision(15) << error[0] << std::endl;
    // std::cout << "ERROR[1] = " << error[1] << std::endl;
    // std::cout << "ERROR[2] = " << error[2] << std::endl;
    // // // std::cout << "ERROR[3] = " << error[3] << std::endl;
    // // // std::cout << "ERROR[4] = " << error[4] << std::endl;
    // // REQUIRE(error[1] < tolerance);

    //Trying to design other criteria besides the approximation error
    //Contar numero de equações condensadas = numero de arestas internas * porder
    int nInternalEdges = 0;
    switch (tshape::Type()){
    case EQuadrilateral:
        nInternalEdges = xdiv*(xdiv-1)*2;
        break;
    case ETriangle:
        nInternalEdges = xdiv*(xdiv-1)*2 + xdiv*xdiv;
        break;
    
    default:
        break;
    }

    
    // REQUIRE(nInternalEdges * pOrder == nEquationsCondensed);


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

    TPZManVector<REAL,3> minX = {-1,-1,-1};
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


template <class tshape>
TPZGeoMesh*
ReadMeshFromGmsh(std::string file_name)
{
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        stringtoint[3]["Domain"] = 1;
        stringtoint[2]["Surfaces"] = 2;

        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(file_name,gmesh);
    }

    return gmesh;
}

























