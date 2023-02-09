/*
  This unit test verifies if the hybridization and semi hybridization techniques are working
  for any specified polynomial order and topology.
  
*/
#include <catch2/catch.hpp>
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
#include <valgrind/callgrind.h>
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

std::ofstream rprint("results_Harmonic2D.txt",std::ofstream::out);
std::ofstream printerrors("results_errors.txt",std::ofstream::out);

void AssociateElementsDuplConnects(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup, int keepLagrangian);
void GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh, int keepLagrangian);
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
    #define TEST
    // const int pOrder = GENERATE(9,10,11,12,13,14,15);
    const int pOrder = GENERATE(1);

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

    if (DIM == 2) nDivs = {2,1};
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
    hdivCreator.IsRigidBodySpaces() = true;
    hdivCreator.SetDefaultOrder(pOrder);
    hdivCreator.SetExtraInternalOrder(0);
    hdivCreator.SetShouldCondense(false);
    hdivCreator.HybridType() = HybridizationType::ESemi;
    // hdivCreator.HybridType() = HybridizationType::EStandard;

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

    util.PrintCMeshConnects(cmesh);
    for (auto &con : cmesh->ConnectVec())
    {
        int64_t seqNum = con.SequenceNumber();
        if (con.IsCondensed()) continue;
        if (seqNum < 0) continue;

        if (con.LagrangeMultiplier() == 3){
            con.IncrementElConnected();
        }
    }   
    
    // TPZCompMeshTools::CreatedCondensedElements(cmesh,3,false);
    GroupAndCondenseElements(cmesh,3);

    hdivCreator.PrintAllMeshes(cmesh);


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
    // bool sparse = true;
    bool sparse = true;
    
    if (sparse){
    // if (approxSpace == TPZHDivApproxSpaceCreator<STATE>::EDuplicatedConnects){
        // TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::EDefault);
        // CALLGRIND_START_INSTRUMENTATION;
        // CALLGRIND_TOGGLE_COLLECT;
        TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::EMHMSparse);
        // TPZMatRedSolver<STATE> solver(an,matBCAll,TPZMatRedSolver<STATE>::ESparse);
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

    {
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












void GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh, int keepLagrangian) {
        
    const int dim = mcmesh->Dimension();
    
    TPZVec<int64_t> elementgroup; // it is resized inside AssociateElements()
    AssociateElementsDuplConnects(mcmesh,elementgroup,keepLagrangian);
    
    for (int i = 0; i < elementgroup.size(); i++)
    {
        std::cout << "ElGroup [" << i << "]=" << elementgroup[i] << std::endl;
    }
    

    int64_t nel = elementgroup.size();

    std::map<int64_t, TPZElementGroup *> groupmap;
    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (int64_t el = 0; el<nel; el++) {
        int64_t groupnum = elementgroup[el];
        if(groupnum == -1) continue;
        auto iter = groupmap.find(groupnum);
        if (groupmap.find(groupnum) == groupmap.end()) {
            int64_t index;
            TPZElementGroup *elgr = new TPZElementGroup(*mcmesh);
            groupmap[groupnum] = elgr;
            elgr->AddElement(mcmesh->Element(el));
        }
        else
        {
            iter->second->AddElement(mcmesh->Element(el));
        }
    }
    mcmesh->ComputeNodElCon();

    nel = mcmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = mcmesh->Element(el);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (!cel) {
            continue;
        }
        if (elgr) {
            TPZCondensedCompEl *cond = new TPZCondensedCompEl(elgr);
            cond->SetKeepMatrix(false);
        }
    }
    mcmesh->CleanUpUnconnectedNodes();
    
    
}



void AssociateElementsDuplConnects(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup, int keepLagrangian)
{
    
    int64_t nel = cmesh->NElements();
    elementgroup.Resize(nel, -1);
    elementgroup.Fill(-1);

    int64_t nconnects = cmesh->NConnects();
    TPZVec<int64_t> groupindex(nconnects, -1);
    int dim = cmesh->Dimension();
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        cel->LoadElementReference();
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != dim) {
            continue;
        }
        elementgroup[cel->Index()] = cel->Index();
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        int count = 0;
        for (auto cindex : connectlist) {
#ifdef PZDEBUG
            if (groupindex[cindex] != -1) {
                continue;
                // DebugStop();
            }
#endif
            TPZConnect &c = cel->Connect(count);
            count++;
            if (c.LagrangeMultiplier() == keepLagrangian) {
                c.IncrementElConnected();
                continue;
            }
            groupindex[cindex] = cel->Index();
        }
    }

    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        int64_t celindex = cel->Index();
        
        TPZVec<int> connectgroup(connectlist.size());
        for(int i=0; i<connectlist.size(); i++) connectgroup[i] = groupindex[connectlist[i]];
        int64_t groupfound = -1;
        //Already sets the element group to the volumetric elemet
        if (cel->Dimension() == dim) elementgroup[celindex] = celindex; 
        int count = 0;
        for (auto cindex : connectlist) {
            TPZConnect &c = cel->Connect(count);
            count++;
            if (c.LagrangeMultiplier() == keepLagrangian) {
                c.IncrementElConnected();
                continue;
            }
            //Gets only the first connect to avoid trouble with the duplicated connect
            if (groupindex[cindex] != -1 && elementgroup[celindex] == -1) {
                elementgroup[celindex] = groupindex[cindex];
                //two connects in the same element can't belong to different computational element groups,
                //but in interface elements, some connects might belong to a certain element group, while others might not be initialized.
                if(groupfound != -1 && groupfound != groupindex[cindex])
                {
                    // DebugStop();
                }
                groupfound = groupindex[cindex];
            }
        }
    }
}










