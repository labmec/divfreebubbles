/*
  This unit test verifies if the hybridization and semi hybridization techniques are working
  for any specified polynomial order and topology.
  
*/
#include <catch2/catch.hpp>
#include <TPZGeoMeshTools.h>
#include "TPZApproxSpaceKernelHdiv.h"
#include "TPZKernelHdivUtils.h"
#include "TPZAnalyticSolution.h"
#include <TPZGmshReader.h>

#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"

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
		case HDivFamily::EHCurlNoGrads:
		{
			return "EHCurlNoGrads";
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
void TestHybridization(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily, TPZApproxSpaceKernelHdiv<STATE>::MSpaceType &approxSpace);

TEST_CASE("Hybridization test")
{

    const int xdiv = GENERATE(2,3);
    const int pOrder = GENERATE(2,3,4);
    HDivFamily hdivfam = GENERATE(HDivFamily::EHDivConstant,HDivFamily::EHCurlNoGrads);//HDivKernel Should not work for tetrahedra
    TPZApproxSpaceKernelHdiv<STATE>::MSpaceType approxSpace = GENERATE(TPZApproxSpaceKernelHdiv<STATE>::ENone);
                                                                        // TPZApproxSpaceKernelHdiv<STATE>::EFullHybrid);
                                                                    //    TPZApproxSpaceKernelHdiv<STATE>::ESemiHybrid);
    
    // TestHybridization<pzshape::TPZShapeTriang>(xdiv,pOrder,hdivfam,approxSpace);
    // TestHybridization<pzshape::TPZShapeQuad>(xdiv,pOrder,hdivfam,approxSpace);
    TestHybridization<pzshape::TPZShapeTetra>(xdiv,pOrder,hdivfam,approxSpace); // HDivConstant not working because of the number of equations. We need to change it to HCurlNoGrad
    TestHybridization<pzshape::TPZShapeCube>(xdiv,pOrder,hdivfam,approxSpace);
}

//Analytical solution
constexpr int solOrder{2};
auto exactSol = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    const auto &d = 1.; // distance between injection and production wells
    u[0]= x*x-y*y ;
    gradU(0,0) = -2*x;
    gradU(1,0) = 2.*y;
    gradU(2,0) = 0.;
    
    // u[0] = x*x*x*y - y*y*y*x;
    // gradU(0,0) = (3.*x*x*y - y*y*y);
    // gradU(1,0) = (x*x*x - 3.*y*y*x);

    // u[0]= x;
    // gradU(0,0) = 1.;
    // gradU(1,0) = 0.;
    // gradU(2,0) = 0.;

    // REAL aux = 1./sinh(sqrt(2)*M_PI);
    // u[0] = sin(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2)*M_PI*z)*aux;
    // gradU(0,0) = M_PI*cos(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2)*M_PI*z)*aux;
    // gradU(1,0) = M_PI*cos(M_PI*y)*sin(M_PI*x)*sinh(sqrt(2)*M_PI*z)*aux;
    // gradU(2,0) = sqrt(2)*M_PI*cosh(sqrt(2)*M_PI*z)*sin(M_PI*x)*sin(M_PI*y)*aux;
};



template<class tshape>
void TestHybridization(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily, TPZApproxSpaceKernelHdiv<STATE>::MSpaceType &approxSpace)
{

    std::cout << "\nTest Case: \nTopology = " << MElementType_Name(tshape::Type()) << 
                 ", xdiv = " << xdiv << ", pOrder = " << pOrder << ", Approximation space = " << MHDivFamily_Name(hdivfamily) << "\n\n "; 
    
    int DIM = tshape::Dimension;
    TPZVec<int> nDivs;

    if (DIM == 2) nDivs = {xdiv,xdiv};
    if (DIM == 3) nDivs = {xdiv,xdiv,xdiv};
    
    // Creates/import a geometric mesh
    auto gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary);
    // auto gmesh = ReadMeshFromGmsh<tshape>("../mesh/1element.msh");

    // Util for HDivKernel printing and solving
    TPZKernelHdivUtils<STATE> util;

    // Creates the approximation space generator
    TPZApproxSpaceKernelHdiv<STATE> createSpace(gmesh, approxSpace, hdivfamily);

    //Insert here the BC material id's to be hybridized
    std::set<int> matBCHybrid={};
    std::set<int> matBCNeumann={};
    std::set<int> matBCDirichlet={EBoundary};
    std::set<int> matBCAll;
    std::set_union(matBCNeumann.begin(),matBCNeumann.end(),matBCDirichlet.begin(),matBCDirichlet.end(),std::inserter(matBCAll, matBCAll.begin()));

    //Setting material ids      
    createSpace.fConfig.fDomain = EDomain;
    createSpace.SetPeriferalMaterialIds(EWrap,EPressureHyb,EIntface,EPont,matBCHybrid,matBCAll);
    createSpace.SetPOrder(pOrder);
    createSpace.Initialize();
    // util.PrintGeoMesh(gmesh);

    //In the case of hybridized HDivConstant, we need 2 pressure meshes, so a total of 3. Otherwise, only 2 CompMeshes are needed 
    int nMeshes = 2;
    if (approxSpace != TPZApproxSpaceKernelHdiv<STATE>::ENone && hdivfamily == HDivFamily::EHDivConstant) {
        nMeshes = 3;
    }
    TPZVec<TPZCompMesh *> meshvector;
    meshvector.Resize(nMeshes);

    //Flux mesh
    meshvector[0] = createSpace.CreateFluxCMesh();
    std::string fluxFile = "FluxCMesh";
    util.PrintCompMesh(meshvector[0],fluxFile);
    // std::cout << "Flux mesh \n";
    // util.PrintCMeshConnects(meshvector[0]);
    
    //Pressure mesh
    meshvector[1]  = createSpace.CreatePressureCMesh();
    std::string presFile = "PressureCMesh";
    util.PrintCompMesh(meshvector[1],presFile);
    // std::cout << "Pressure mesh \n";
    // util.PrintCMeshConnects(meshvector[1]);

    if (approxSpace != TPZApproxSpaceKernelHdiv<STATE>::ENone && hdivfamily == HDivFamily::EHDivConstant) {
        meshvector[2] = createSpace.CreatePressureCMeshHybridizedHDivConstant();
        util.PrintCompMesh(meshvector[1],presFile);
        // std::cout << "Pressure mesh2 \n";
        // util.PrintCMeshConnects(meshvector[2]);
    }
    
    //Multiphysics mesh
    auto * cmesh = createSpace.CreateMultiphysicsCMesh(meshvector,exactSol,matBCNeumann,matBCDirichlet);
    std::string multFile = "MultiCMesh";
    util.PrintCompMesh(cmesh,multFile);
    // std::cout << "Multi mesh \n";
    // util.PrintCMeshConnects(cmesh);

    // Number of equations without condense elements
    int nEquationsFull = cmesh->NEquations();
    std::cout << "Number of equations = " << nEquationsFull << std::endl;

    // Group and condense the elements
    if (DIM == 2){
        // createSpace.Condense(cmesh);
    }

    //Number of condensed problem.
    int nEquationsCondensed = cmesh->NEquations();

    //Create analysis environment
    TPZLinearAnalysis an(cmesh,false);

    //Solve problem
    bool filter = false;
    if (DIM == 3 && hdivfamily != HDivFamily::EHDivConstant) filter = true;
    createSpace.Solve(an, cmesh, true, filter);
    an.SetExact(exactSol,pOrder);

    //Print results
    util.PrintResultsMultiphysics(meshvector,an,cmesh);

    //Compute error
    std::ofstream anPostProcessFile("postprocess.txt");
    TPZManVector<REAL,5> error;
    int64_t nelem = cmesh->NElements();
    cmesh->LoadSolution(cmesh->Solution());
    cmesh->ExpandSolution();
    cmesh->ElementSolution().Redim(nelem, 5);
    an.PostProcessError(error,false,anPostProcessFile);
    
    //Check error
    REAL tolerance = 1.e-10;
    std::cout << "ERROR = " << error[1] << std::endl;
    REQUIRE(error[1] < tolerance);

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
    
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls);
    
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

























