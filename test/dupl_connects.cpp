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
#include <pzstepsolver.h> //for TPZStepSolver
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
//#include <valgrind/callgrind.h>

// ----- Run tests with or without main -----
#define RUNWITHMAIN

#ifndef RUNWITHMAIN
#include <catch2/catch.hpp>
#endif

std::ofstream rprint("results_Harmonic2D.txt",std::ofstream::out);
std::ofstream printerrors("results_errors.txt",std::ofstream::out);

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

#ifndef RUNWITHMAIN
TEST_CASE("Hybridization test")
{
    #define TEST
    // const int pOrder = GENERATE(9,10,11,12,13,14,15);
    const int pOrder = GENERATE(2);

    const int xdiv = 2;//GENERATE(5,10,15);
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

#else
int main() {
  const int xdiv = 2;
  const int pOrder = 2;
  HDivFamily hdivfam = HDivFamily::EHDivConstant;
  TestHybridization<pzshape::TPZShapeQuad>(xdiv,pOrder,hdivfam);
  return 0;
}
#endif

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

    // REAL aux = 1./sinh(sqrt(2)*M_PI);
    // u[0] = sin(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2)*M_PI*z)*aux;
    // gradU(0,0) = M_PI*cos(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2)*M_PI*z)*aux;
    // gradU(1,0) = M_PI*cos(M_PI*y)*sin(M_PI*x)*sinh(sqrt(2)*M_PI*z)*aux;
    // gradU(2,0) = sqrt(2)*M_PI*cosh(sqrt(2)*M_PI*z)*sin(M_PI*x)*sin(M_PI*y)*aux;

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
    // u[0] = x + y;
    // gradU(0,0) = 1.;
    // gradU(1,0) = 1.;
    u[0] = x*x - y*y;
    gradU(0,0) = 2.*x;
    gradU(1,0) = -2.*y;

    // u[0] = 0.5*(x)*x+0.5*(y)*y-(z)*z;
    // gradU(0,0) = -x;//(x-1)*(y-1)*y*(z-1)*z + x*(y-1)*y*(z-1)*z;
    // gradU(1,0) = -y;//(x-1)*x*(y-1)*(z-1)*z + (x-1)*x*y*(z-1)*z;
    // gradU(2,0) = 2.*z;//(x-1)*x*(y-1)*y*(z-1) + (x-1)*x*(y-1)*y*z;

};

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

    if (DIM == 2) nDivs = {20,20};
    if (DIM == 3) nDivs = {xdiv,xdiv,xdiv};
    
    // Creates/import a geometric mesh  
    auto gmesh = CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary);

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
//    hdivCreator.SetShouldCondense(true);
    hdivCreator.SetShouldCondense(false);
    // hdivCreator.HybridType() = HybridizationType::ESemi;
    hdivCreator.HybridType() = HybridizationType::EStandard;

    // Prints gmesh mesh properties
    std::string vtk_name = "geoMesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

    //Insert Materials
    TPZMixedDarcyFlow* matdarcy = new TPZMixedDarcyFlow(EDomain,DIM);
    matdarcy->SetConstantPermeability(1.);
    matdarcy->SetExactSol(exactSol,4);

    hdivCreator.InsertMaterialObject(matdarcy);

    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,0.);
    TPZBndCondT<STATE> *BCond1 = matdarcy->CreateBC(matdarcy, EBoundary, 0, val1, val2);
    BCond1->SetForcingFunctionBC(exactSol,4);
    hdivCreator.InsertMaterialObject(BCond1);

    //Multiphysics mesh
    TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
    std::string txt = "cmesh.txt";
    std::ofstream myfile(txt);
    cmesh->Print(myfile);

  
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
    TPZLinearAnalysis an(cmesh,true);
    an.SetExact(exactSol,solOrder);

    std::set<int> matBCAll = {EBoundary};
    // Solve problem
    bool sparse = true;
//    bool sparse = false;
    
    if (sparse){
    // if (approxSpace == TPZHDivApproxSpaceCreator<STATE>::EDuplicatedConnects){
        // TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::EDefault);
        // CALLGRIND_START_INSTRUMENTATION;
        // CALLGRIND_TOGGLE_COLLECT;
        TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::ESparse);
        clock2.start();
        solver.Solve(rprint);
        clock2.stop();
        // CALLGRIND_TOGGLE_COLLECT;
        // CALLGRIND_STOP_INSTRUMENTATION;
        // std::cout << "Time SOLVER = " << clock2 << std::endl;

        // bool filter = false;
        // if (DIM == 3 && hdivfamily == HDivFamily::EHDivKernel) filter = true;
        // createSpace.Solve(an, cmesh, true, filter);

    } else {
        bool domHyb = false;
        util.SolveProblemDirect(an,cmesh,false,domHyb);
    }

    

    clock.stop();
    // std::cout << "Time running = " << clock << std::endl;

    // //Print results
    // {
    //     TPZSimpleTimer postProc("Post processing1");
        // util.PrintResultsMultiphysics(cmesh->MeshVector(),an,cmesh);
    // }

    if(0) {
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), cmesh);
        TPZSimpleTimer postProc("Post processing2");
        const std::string plotfile = "myfile";//sem o .vtk no final
        constexpr int vtkRes{0};
    

        TPZVec<std::string> fields = {
        "Pressure",
        "ExactPressure",
        "Flux",
        "ExactFlux"};
        auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

        vtk.Do();
    }
    else {
      TPZManVector<std::string,10> scalnames(2), vecnames(2);


      scalnames[0] = "Pressure";
      scalnames[1] = "ExactPressure";
      vecnames[0]= "Flux";
      vecnames[1]= "ExactFlux";

      constexpr int resolution{0};
      std::string plotfile = "myfile.vtk";
      an.DefineGraphMesh(cmesh->Dimension(),scalnames,vecnames,plotfile);
      an.PostProcess(resolution,cmesh->Dimension());

    }
    std::string txt2 = "cmeshSol.txt";
    std::ofstream myfile2(txt2);
    cmesh->Print(myfile2);

    // //vamos supor que vc atualiza a solucao, roda de novo, sei la
    // vtk.Do();

    //Compute error
    std::ofstream anPostProcessFile("postprocess.txt");
    TPZManVector<REAL,5> error;
    int64_t nelem = cmesh->NElements();
    cmesh->LoadSolution(cmesh->Solution());
    cmesh->ExpandSolution();
    cmesh->ElementSolution().Redim(nelem, 5);
    an.PostProcessError(error,false,anPostProcessFile);
    
    printerrors << xdiv << std::scientific << std::setprecision(8) << " " << error[0] << " " 
     << error[1] << " " << error[2] << " "  << error[3] << " "  << error[4] << std::endl;

    //Check error
    // REAL tolerance = 1.e-6;
    std::cout << "ERROR[0] = " << std::scientific << std::setprecision(15) << error[0] << std::endl;
    std::cout << "ERROR[1] = " << error[1] << std::endl;
    std::cout << "ERROR[2] = " << error[2] << std::endl;
    std::cout << "ERROR[3] = " << error[3] << std::endl;
    std::cout << "ERROR[4] = " << error[4] << std::endl;
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


























