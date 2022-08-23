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

std::ofstream rprint("results_Harmonic2D.txt",std::ofstream::out);

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

/** @brief Returns the name of the Hybridization type. */
inline std::string ApproxSpaceKernelHDiv_Name(TPZHDivApproxSpaceCreator<STATE>::MSpaceType approxSpace)
{
	switch (approxSpace)
	{
		case TPZHDivApproxSpaceCreator<STATE>::ENone:
		{
			return "ENone";
		}
		case TPZHDivApproxSpaceCreator<STATE>::EFullHybrid:
		{
			return "EFullHybrid";
		}
		case TPZHDivApproxSpaceCreator<STATE>::ESemiHybrid:
		{
			return "ESemiHybrid";
		}
		case TPZHDivApproxSpaceCreator<STATE>::EDuplicatedConnects:
		{
			return "EDuplicatedConnects";
		}
		default:
        {
            return "TPZHDivApproxSpaceCreator not found!";
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
void TestHybridization(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily, TPZHDivApproxSpaceCreator<STATE>::MSpaceType &approxSpace);

TEST_CASE("Hybridization test")
{
    const int pOrder = GENERATE(1);
    // const int pOrder = GENERATE(2,3,4,5);

    const int xdiv = 2;//GENERATE(50);
    // const int xdiv = GENERATE(2,5,10,15,20,25,30,35,40,45,50,60,70,80,90,100,120,140,160,180,200);
    // const int xdiv = GENERATE(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16);
    // const int xdiv = GENERATE(2,3,4,5,6,7,8);
    // HDivFamily hdivfam = GENERATE(HDivFamily::EHDivConstant,HDivFamily::EHDivKernel);
    // HDivFamily hdivfam = GENERATE(HDivFamily::EHDivKernel);
    HDivFamily hdivfam = GENERATE(HDivFamily::EHDivConstant);
    // HDivFamily hdivfam = GENERATE(HDivFamily::EHDivStandard);
    // HDivFamily hdivfam = GENERATE(HDivFamily::EHDivStandard,HDivFamily::EHDivConstant);
    // TPZHDivApproxSpaceCreator<STATE>::MSpaceType approxSpace = GENERATE(TPZHDivApproxSpaceCreator<STATE>::EFullHybrid);
    TPZHDivApproxSpaceCreator<STATE>::MSpaceType approxSpace = GENERATE(TPZHDivApproxSpaceCreator<STATE>::EDuplicatedConnects);
    // TPZHDivApproxSpaceCreator<STATE>::MSpaceType approxSpace = GENERATE(TPZHDivApproxSpaceCreator<STATE>::EDuplicatedConnects);
    
    // TestHybridization<pzshape::TPZShapeTriang>(xdiv,pOrder,hdivfam,approxSpace);
    TestHybridization<pzshape::TPZShapeQuad>(xdiv,pOrder,hdivfam,approxSpace); 
    // TestHybridization<pzshape::TPZShapeTetra>(xdiv,pOrder,hdivfam,approxSpace); 
    // TestHybridization<pzshape::TPZShapeCube>(xdiv,pOrder,hdivfam,approxSpace);
}


template<class tshape>
void TestHybridization(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily, TPZHDivApproxSpaceCreator<STATE>::MSpaceType &approxSpace)
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    std::cout << "\nTest Case: \nTopology = " << MElementType_Name(tshape::Type()) << 
                 ", xdiv = " << xdiv << ", pOrder = " << pOrder << 
                 ", Approximation space = " << MHDivFamily_Name(hdivfamily) <<
                 ", Hybridization = " << ApproxSpaceKernelHDiv_Name(approxSpace) << "\n\n "; 
    
    int DIM = tshape::Dimension;
    TPZVec<int> nDivs;

    if (DIM == 2) nDivs = {xdiv,xdiv};
    if (DIM == 3) nDivs = {xdiv,xdiv,xdiv};
    
    // Creates/import a geometric mesh
    auto gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary);
    // auto gmesh = ReadMeshFromGmsh<tshape>("../mesh/1tetra.msh");

    // Util for HDivKernel printing and solving
    TPZKernelHdivUtils<STATE> util;

    // Creates the approximation space generator
    TPZHDivApproxSpaceCreator<STATE> createSpace(gmesh, approxSpace, hdivfamily);

    //Insert here the BC material id's to be hybridized
    std::set<int> matBCHybrid={};
    std::set<int> matBCNeumann={};
    std::set<int> matBCDirichlet={EBoundary};
    std::set<int> matBCAll;
    std::set_union(matBCNeumann.begin(),matBCNeumann.end(),matBCDirichlet.begin(),matBCDirichlet.end(),std::inserter(matBCAll, matBCAll.begin()));

    //Setting material ids      
    createSpace.fConfig.fDomain = EDomain;
    createSpace.SetMaterialIds(EWrap,EPressureHyb,EIntface,EPont,matBCHybrid,matBCAll);
    createSpace.SetPOrder(pOrder);
    createSpace.SetMixedElasticity();
    createSpace.Initialize();
    
    //Setting material ids      
    // util.PrintGeoMesh(gmesh);

    TPZAnalyticSolution *gAnalytic = 0;
    if(DIM == 2)
    {
        TElasticity2DAnalytic *elas = new TElasticity2DAnalytic;
        elas->gE = 1;
        elas->gPoisson = 0.0;
        elas->fProblemType = TElasticity2DAnalytic::EDispx;
        elas->fPlaneStress = 0;
        gAnalytic = elas;
    }
    else if(DIM == 3)
    {
        TElasticity3DAnalytic *elas = new TElasticity3DAnalytic;
        elas->fE = 1.;//206.8150271873455;
        elas->fPoisson = 0.0;//0.3040039545229857;
        elas->fProblemType = TElasticity3DAnalytic::EStretchx;
        gAnalytic = elas;
    }


    unsigned int stressPOrder = pOrder + 1; //Polynomial order of the approximation
    int stressInternalPOrder = stressPOrder + 1;
    int displacementPOrder = tshape::Type() == ETriangle ? stressInternalPOrder - 1 : stressInternalPOrder;
    int rotationPOrder = displacementPOrder;
    if(hdivfamily == HDivFamily::EHDivConstant) rotationPOrder = 0;
    //Creating computational mesh:
    //Creates the computational mesh for the stress field (HDiv)
    TPZCompMesh *cmesh_S_HDiv = createSpace.CreateFluxCMesh();
    createSpace.ChangeInternalOrder(cmesh_S_HDiv, stressInternalPOrder);
    //Creates the computational mesh for the displacement field (H1 disconnected)
    TPZCompMesh *cmesh_U_HDiv = createSpace.CreatePressureCMesh();
    //Creates the computational mesh for the rotation field (Discontinuous)
    TPZCompMesh *cmesh_P_HDiv = createSpace.CreateRotationCmesh(gmesh, rotationPOrder, 1. / xdiv);
    // creates the mesh for distributed forces in each element
    TPZCompMesh *cmesh_distributedforce = createSpace.CreateConstantCmesh(gmesh, 2);
    // creates the computational mesh representing the average displacement and rotation
    TPZCompMesh *cmesh_averagedisp = createSpace.CreateConstantCmesh(gmesh, 4);

    // TPZManVector<TPZCompMesh*, 5> meshvector_HDiv(5);
    TPZManVector<TPZCompMesh*, 6> meshvector_HDiv(6);
    meshvector_HDiv[0] = cmesh_S_HDiv;
    meshvector_HDiv[1] = cmesh_U_HDiv;
    meshvector_HDiv[2] = cmesh_P_HDiv;
    meshvector_HDiv[3] = cmesh_distributedforce;
    meshvector_HDiv[4] = cmesh_averagedisp;

    meshvector_HDiv[5] = createSpace.CreatePressureCMeshHybridizedHDivConstant();

    //Creates the multi-physics computational mesh
    auto *cmesh_m_HDiv = createSpace.CreateMultiphysicsCMeshElasticity(meshvector_HDiv,gAnalytic,matBCNeumann,matBCDirichlet);
            
#ifdef PZDEBUG
    {
        //Prints the stress computational mesh in txt format
        std::ofstream filecS("MalhaC_S.txt");
        //Prints the displacement computational mesh in txt format
        std::ofstream filecU("MalhaC_U.txt");
        //Prints the rotation computational mesh in txt format
        std::ofstream filecP("MalhaC_P.txt");
        //Prints the distributed force mesh
        std::ofstream filedf("MalhaC_DistForce.txt");
        //Prints the average displacement mesh
        std::ofstream fileavdisp("MalhaC_AvDisp.txt");

        std::ofstream filehdivhybr("MalhaC_Hybrid.txt");
        cmesh_S_HDiv->Print(filecS);
        cmesh_U_HDiv->Print(filecU);
        cmesh_P_HDiv->Print(filecP);
        cmesh_distributedforce->Print(filedf);
        cmesh_averagedisp->Print(fileavdisp);
        meshvector_HDiv[5]->Print(filehdivhybr);
        
    }
#endif

    util.PrintCompMesh(cmesh_m_HDiv,"MultiCMeshBefore");
    // Group and condense the elements
    if (approxSpace == TPZHDivApproxSpaceCreator<STATE>::EDuplicatedConnects){
        createSpace.CondenseDuplicatedConnects(cmesh_m_HDiv);
    }else {
        TPZCompMeshTools::CondenseElements(cmesh_m_HDiv, 1, false);
    }
    
    // std::cout << "Multi mesh \n";
    // util.PrintCMeshConnects(cmesh);
    
    //Number of condensed problem.
    int nEquationsCondensed = cmesh_m_HDiv->NEquations();

    //Create analysis environment
    TPZLinearAnalysis an(cmesh_m_HDiv,true);
    an.SetExact(gAnalytic->ExactSolution());

    std::string multFile = "MultiCMesh";
    util.PrintCompMesh(cmesh_m_HDiv,multFile);

    //Solve problem
    if (approxSpace == TPZHDivApproxSpaceCreator<STATE>::EDuplicatedConnects){
        TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::EDefault);
        // TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::ESparse);
        solver.Solve(rprint);
        // std::cout << "Time SOLVER = " << clock2 << std::endl;

        // bool filter = false;
        // if (DIM == 3 && hdivfamily == HDivFamily::EHDivKernel) filter = true;
        // createSpace.Solve(an, cmesh, true, filter);

    } else {
        //Equation filter (spanning trees), true if 3D and HDivKernel 
        bool filter = false;bool domainhybr=false;
        if (DIM == 3 && hdivfamily == HDivFamily::EHDivKernel) filter = true;
        // createSpace.Solve(an, cmesh, true, filter);
        util.SolveProblemDirect(an,cmesh_m_HDiv,filter,domainhybr);
    }

    // std::cout << "Time running = " << clock << std::endl;

    // //Print results
    // {
    //     TPZSimpleTimer postProc("Post processing1");
    //     util.PrintResultsMultiphysics(meshvector,an,cmesh);
    // }

    {
        TPZSimpleTimer postProc("Post processing2");
        const std::string plotfile = "myfile";//sem o .vtk no final
        constexpr int vtkRes{1};
    

        TPZVec<std::string> fields = {
        // "ExactDisplacement",
        // "ExactStress",
        "Displacement",
        "SigmaX",
        "SigmaY",
        "TauXY"
        };
        auto vtk = TPZVTKGenerator(cmesh_m_HDiv, fields, plotfile, vtkRes);

        vtk.Do();
    }
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

























