/*
  This unit test verifies if the hybridization and semi hybridization techniques are working
  for any specified polynomial order and topology.
  
*/
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
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "pzbuildmultiphysicsmesh.h"


enum EMatid  {ENone, EDomain, EBoundary, EPont, EWrap, EIntface, EPressureHyb};

// The test function
template<class tshape>
void RunMHM(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily, TPZHDivApproxSpaceCreator<STATE>::MSpaceType &approxSpace);


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
void RunMHM(const int &xdiv, const int &pOrder, HDivFamily &hdivfamily, TPZHDivApproxSpaceCreator<STATE>::MSpaceType &approxSpace)
{
    // Util for HDivKernel printing and solving
    TPZKernelHdivUtils<STATE> util;

    int DIM = tshape::Dimension;
    TPZVec<int> nDivs;

    if (DIM == 2) nDivs = {xdiv,xdiv};
    if (DIM == 3) nDivs = {xdiv,xdiv,xdiv};
    
    // Creates/import a geometric mesh
    auto gmesh = util.CreateGeoMesh<tshape>(nDivs, EDomain, EBoundary);
    // auto gmesh = ReadMeshFromGmsh<tshape>("../mesh/1tetra.msh");   

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
    createSpace.Initialize();
    // util.PrintGeoMesh(gmesh);

    //In the case of hybridized HDivConstant, we need 2 pressure meshes, so a total of 3. Otherwise, only 2 CompMeshes are needed 
    int nMeshes = 2;
    TPZVec<TPZCompMesh *> meshvector;
    meshvector.Resize(nMeshes);

    //Flux mesh
    meshvector[0] = createSpace.CreateFluxCMesh();
    // std::string fluxFile = "FluxCMesh";
    // util.PrintCompMesh(meshvector[0],fluxFile);
    // std::cout << "Flux mesh \n";
    // util.PrintCMeshConnects(meshvector[0]);
    
    //Pressure mesh
    meshvector[1]  = createSpace.CreatePressureCMesh();
    // std::string presFile = "PressureCMesh";
    // util.PrintCompMesh(meshvector[1],presFile);
    // std::cout << "Pressure mesh \n";
    // util.PrintCMeshConnects(meshvector[1]);

    //Multiphysics mesh
    auto * cmesh = createSpace.CreateMultiphysicsCMesh(meshvector,exactSol,matBCNeumann,matBCDirichlet);
    // std::string multFile = "MultiCMesh";
    // std::cout << "Multi mesh \n";
    // util.PrintCMeshConnects(cmesh);

    // Number of equations without condense elements
    int nEquationsFull = cmesh->NEquations();
    std::cout << "Number of equations = " << nEquationsFull << std::endl;

    TPZTimer clock,clock2;
    clock.start();

    std::string multFile = "MultiCMesh";
    util.PrintCompMesh(cmesh,multFile);
    
    //Number of condensed problem.
    int nEquationsCondensed = cmesh->NEquations();

    //Create analysis environment
    TPZLinearAnalysis an(cmesh,true);
    an.SetExact(exactSol,solOrder);


    //Solve problem
    //Equation filter (spanning trees), true if 3D and HDivKernel 
    bool filter = false;
    if (DIM == 3 && hdivfamily == HDivFamily::EHDivKernel) filter = true;
    createSpace.Solve(an, cmesh, true, filter);
    
    {
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
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
    const int xdiv = 5;
    // HDivFamily hdivfam = HDivFamily::EHDivStandard;
    HDivFamily hdivfam = HDivFamily::EHDivConstant;
    TPZHDivApproxSpaceCreator<STATE>::MSpaceType approxSpace = TPZHDivApproxSpaceCreator<STATE>::ENone;
    
    RunMHM<pzshape::TPZShapeQuad>(xdiv,pOrder,hdivfam,approxSpace); 

    return 0;
}























