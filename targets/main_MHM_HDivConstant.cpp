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
#include "TPZMHMGeoMeshCreator.h"
#include "TPZMHMHDivApproxCreator.h"
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZLinearAnalysis.h"
#include "Common.h"
#include "TPZAnalyticSolution.h"
#include "TPZMatRedSolver.h"


enum EMatid  {ENone, EDomain, EBottom, ETop, ELeft, ERight};

// The test function
template<class tshape>
void RunMHM(const int &xdiv, const int &pOrder);

/**
   @brief Creates a geometric mesh with elements of a given type on a unit square or cube (depending on the mesh dimension).
   @param[in] meshType element type to be created.
   @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
template<class tshape>
TPZAutoPointer<TPZGeoMesh>
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenterindices);

//Analytical solution
constexpr int solOrder{4};
auto exactSol = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    u[0] = exp(M_PI*x)*sin(M_PI*y);
    gradU(0,0) = M_PI*exp(M_PI*x)*sin(M_PI*y);
    gradU(1,0) = M_PI*exp(M_PI*x)*cos(M_PI*y);
};

template<class tshape>
void RunMHM(const int &xdiv, const int &pOrder)
{
    // Util for HDivKernel printing and solving
    TPZKernelHdivUtils<STATE> util;
    TPZMHMGeoMeshCreator mhm_gcreator;
    mhm_gcreator.fSkeletonMatId = 100;
    mhm_gcreator.fDomainMatId = Emat2;
    mhm_gcreator.fBC1 = Ebc1;
    mhm_gcreator.fBC2 = Ebc2;
    mhm_gcreator.fBC3 = Ebc3;
    mhm_gcreator.fBC4 = Ebc4;
    mhm_gcreator.SetBCMatId(Ebc1);
    mhm_gcreator.SetBCMatId(Ebc2);
    mhm_gcreator.SetBCMatId(Ebc3);
    mhm_gcreator.SetBCMatId(Ebc4);
    

    int DIM = tshape::Dimension;
    TPZVec<int> nDivs;

    if (DIM == 2) nDivs = {xdiv,xdiv};
    if (DIM == 3) nDivs = {xdiv,xdiv,xdiv};
    
    // Creates/import a geometric mesh
    std::string filename;
    filename = "polygon00.txt";
    TPZVec<int64_t> elpartition;
    TPZVec<int64_t> scalingcenterindices;
    TPZAutoPointer<TPZGeoMesh> gmesh = ReadUNSWQuadtreeMesh(filename, mhm_gcreator.fElementPartition, scalingcenterindices);
    

    // // Creates/import a geometric mesh  
    // TPZAutoPointer<TPZGeoMesh> gmesh = CreateGeoMesh<pzshape::TPZShapeQuad>(nDivs, EDomain, mhm_gcreator.fElementPartition, scalingcenterindices);
    
    // std::string vtk_name = "geoMeshMHM.vtk";
    // std::ofstream vtkfile(vtk_name.c_str());
    // TPZVTKGeoMesh::PrintGMeshVTK(gmesh.operator->(), vtkfile, mhm_gcreator.fElementPartition);

    // mhm_gcreator.CreateSkeleton(gmesh);
    mhm_gcreator.AddBoundaryElements(gmesh);
    mhm_gcreator.fElementPartition.Resize(gmesh->NElements(), -1);
    scalingcenterindices.Resize(gmesh->NElements(), -1);

    // LaplaceExact.fExact = TLaplaceExample1::E2SinSin;
    LaplaceExact.fExact = TLaplaceExample1::EX;

    std::cout << "Building computational mesh\n";
    std::map<int,int> matmap;
    matmap[ESkeleton] = Emat1;
    int EPoly = 100;
    matmap[EPoly] = Emat2;

    // std::cout << "scalingcenterindices " << scalingcenterindices << std::endl;
    // std::cout << "element partition " << mhm_gcreator.fElementPartition << std::endl;

    mhm_gcreator.CreateTriangleElements(gmesh, matmap, mhm_gcreator.fElementPartition, scalingcenterindices);

    // mhm_gcreator.CreateSkeleton(gmesh);
    // mhm_gcreator.CreateSubGrids(gmesh);
    // mhm_gcreator.RefineSubGrids(gmesh);
    // mhm_gcreator.RefineSubGrids(gmesh);
    // mhm_gcreator.RefineSkeleton(gmesh);

    // std::string vtk_name = "geoMeshMHM.vtk";
    // std::ofstream vtkfile(vtk_name.c_str());
    // TPZVTKGeoMesh::PrintGMeshVTK(gmesh.operator->(), vtkfile, mhm_gcreator.fElementPartition);

    TPZMHMHDivApproxCreator mhm_ccreator(mhm_gcreator,gmesh);

    mhm_ccreator.HdivFamily() = HDivFamily::EHDivConstant;
    mhm_ccreator.ProbType() = ProblemType::EDarcy;
    mhm_ccreator.IsRigidBodySpaces() = true;
    mhm_ccreator.SetDefaultOrder(pOrder);
    mhm_ccreator.SetExtraInternalOrder(0);
    mhm_ccreator.SetShouldCondense(true);
    mhm_ccreator.HybridType() = HybridizationType::ENone;
    mhm_ccreator.SetPOrderSkeleton(pOrder);

    mhm_ccreator.InsertMaterialObjects(LaplaceExact);
    auto multiCmesh = mhm_ccreator.BuildMultiphysicsCMesh();

    if(1)
    {
        std::cout << "NEQUATIONS = " << multiCmesh ->NEquations() << std::endl;
        std::cout << "NGeoEls = " << gmesh->NElements() << std::endl;
        std::cout << "NCompEls = " << multiCmesh->NElements() << std::endl;

        
        // std::ofstream out("mphysics.txt");
        std::ofstream out("mphysics2.txt");
        multiCmesh->Print(out);
        std::ofstream out2("cmesh_multi.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(multiCmesh, out2);
        util.PrintCMeshConnects(multiCmesh);
    }
    mhm_ccreator.PutinSubstructures(*multiCmesh);
    mhm_ccreator.CondenseElements(*multiCmesh);

    bool mustOptimizeBandwidth = true;
    TPZCompMesh *SBFem = multiCmesh;
    TPZLinearAnalysis * Analysis = new TPZLinearAnalysis(SBFem,mustOptimizeBandwidth);

    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> strmat(multiCmesh);
    strmat.SetNumThreads(10);
    Analysis->SetStructuralMatrix(strmat);

    int64_t neq = multiCmesh->NEquations();
    
    if(neq > 20000)
    {
        std::cout << "Entering Assemble Equations\n";
        std::cout.flush();
    }
    
    if (mhm_ccreator.HybridType() == HybridizationType::ESemi){
        std::set<int> matBCAll={Ebc1,Ebc2,Ebc3,Ebc4};
        TPZMatRedSolver<STATE> solver(*Analysis,matBCAll,TPZMatRedSolver<STATE>::EMHMSparse);
        solver.Solve(std::cout);
    } else {
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        Analysis->SetSolver(step);
        
        Analysis->Assemble();
        
        if(neq > 20000)
        {
            std::cout << "Entering Solve\n";
            std::cout.flush();
        }
        
        Analysis->Solve();

    }

    
//    try {
//        an->Assemble();
//    } catch (...) {
//        exit(-1);
//    }
    
   
    // // Creates the approximation space generator
    // TPZHDivApproxSpaceCreator<STATE> createSpace(gmesh, approxSpace, hdivfamily);

    // //Insert here the BC material id's to be hybridized
    // std::set<int> matBCHybrid={};
    // std::set<int> matBCNeumann={};
    // std::set<int> matBCDirichlet={EBoundary};
    // std::set<int> matBCAll;
    // std::set_union(matBCNeumann.begin(),matBCNeumann.end(),matBCDirichlet.begin(),matBCDirichlet.end(),std::inserter(matBCAll, matBCAll.begin()));

    // //Setting material ids      
    // createSpace.fConfig.fDomain = EDomain;
    // createSpace.SetMaterialIds(EWrap,EPressureHyb,EIntface,EPont,matBCHybrid,matBCAll);
    // createSpace.SetPOrder(pOrder);
    // createSpace.Initialize();
    // // util.PrintGeoMesh(gmesh);

    // //In the case of hybridized HDivConstant, we need 2 pressure meshes, so a total of 3. Otherwise, only 2 CompMeshes are needed 
    // int nMeshes = 2;
    // TPZVec<TPZCompMesh *> meshvector;
    // meshvector.Resize(nMeshes);

    // //Flux mesh
    // meshvector[0] = createSpace.CreateFluxCMesh();
    // // std::string fluxFile = "FluxCMesh";
    // // util.PrintCompMesh(meshvector[0],fluxFile);
    // // std::cout << "Flux mesh \n";
    // // util.PrintCMeshConnects(meshvector[0]);
    
    // //Pressure mesh
    // meshvector[1]  = createSpace.CreatePressureCMesh();
    // // std::string presFile = "PressureCMesh";
    // // util.PrintCompMesh(meshvector[1],presFile);
    // // std::cout << "Pressure mesh \n";
    // // util.PrintCMeshConnects(meshvector[1]);

    // //Multiphysics mesh
    // auto * cmesh = createSpace.CreateMultiphysicsCMesh(meshvector,exactSol,matBCNeumann,matBCDirichlet);
    // // std::string multFile = "MultiCMesh";
    // // std::cout << "Multi mesh \n";
    // // util.PrintCMeshConnects(cmesh);

    // // Number of equations without condense elements
    // int nEquationsFull = cmesh->NEquations();
    // std::cout << "Number of equations = " << nEquationsFull << std::endl;

    // TPZTimer clock,clock2;
    // clock.start();

    // std::string multFile = "MultiCMesh";
    // util.PrintCompMesh(cmesh,multFile);
    
    // //Number of condensed problem.
    // int nEquationsCondensed = cmesh->NEquations();

    // //Create analysis environment
    // TPZLinearAnalysis an(cmesh,true);
    // an.SetExact(exactSol,solOrder);


    // //Solve problem
    // //Equation filter (spanning trees), true if 3D and HDivKernel 
    // bool filter = false;
    // if (DIM == 3 && hdivfamily == HDivFamily::EHDivKernel) filter = true;
    // createSpace.Solve(an, cmesh, true, filter);
    
    {
        // TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
        TPZSimpleTimer postProc("Post processing2");
        const std::string plotfile = "myfile";//sem o .vtk no final
        constexpr int vtkRes{0};
    

        TPZVec<std::string> fields = {
        "Pressure",
        "ExactPressure",
        "Flux",
        "ExactFlux"};
        auto vtk = TPZVTKGenerator(SBFem, fields, plotfile, vtkRes);

        vtk.Do();
    }
    
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

    

}



int main(int argc, char *argv[]){

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    #define TEST
    const int pOrder = 1;
    const int xdiv = 2;
        
    RunMHM<pzshape::TPZShapeTriang>(xdiv,pOrder); 

    return 0;
}


//Create 
template <class tshape>
TPZAutoPointer<TPZGeoMesh>
CreateGeoMesh(TPZVec<int> &nDivs, EMatid volId, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenterindices)
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

    TPZManVector<REAL,3> minX = {-1,-1,0};
    TPZManVector<REAL,3> maxX = {1,1,0};
    int nMats = 2*dim+1;

    //all bcs share the same id
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,volId);
    matIds[0] = volId;
    matIds[1] = Ebc1;
    matIds[2] = Ebc2;
    matIds[3] = Ebc3;
    matIds[4] = Ebc4;
    
    TPZAutoPointer<TPZGeoMesh> gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls);
    // TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(meshType,
    //                     volId,createBoundEls, bcId);

    int nVolumes = nDivs[0]*nDivs[1];
    if (dim == 3) nVolumes *= nDivs[2];
    scalingcenterindices.Resize(nVolumes*6,-1);

    auto nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++)
    {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if (gel->Dimension() == dim){
            int nFacets = gel->NSides(dim-1);
            int nSides = gel->NSides();
            for (int iside = 0; iside < nFacets; iside++){
                TPZGeoElSide gelside(gel,nSides-1-nFacets+iside);
                if (gelside.Neighbour().Element()->Dimension() == dim){
                    TPZGeoElBC gelbc(gel,nSides-1-nFacets+iside,ESkeleton);
                }  
            }

            TPZManVector<REAL,3> midxco(dim,0.);
            gel->CenterPoint(nSides-1,midxco);
            midxco.resize(3);
            int64_t midindex = gmesh->NodeVec().AllocateNewElement();
            gmesh->NodeVec()[midindex].Initialize(midxco, *gmesh);

            scalingcenterindices[midindex] = gel->Index();

            gmesh->DeleteElement(gel);
        }
    }
    gmesh->ResetReference();
    
    
    elpartition.Resize(gmesh->NElements(),-1);
    for (int i = 0; i < gmesh->NElements(); i++)
    {
        elpartition[i]=i;
        TPZGeoEl *gel = gmesh->ElementVec()[i];
        if (!gel) continue;
        gel->SetMaterialId(100);
        // if (gel->Dimension() == dim) DebugStop();
    }
        
    return gmesh;
    
}




















