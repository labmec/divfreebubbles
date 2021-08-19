#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

// #include "TPZGenGrid2D.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZGmshReader.h>
#include <TPZVTKGeoMesh.h>
#include "Poisson/TPZMatPoisson.h" //for TPZMatLaplacian
#include "Projection/TPZL2Projection.h" //for BC in a single point
#include "pzmultiphysicscompel.h"
#include <TPZNullMaterial.h>
#include <TPZNullMaterialCS.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"// for Hdiv problem
#include <TPZBndCond.h> //for TPZBndCond
#include "TPZLinearAnalysis.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzstrmatrixor.h"
#include "pzlog.h"
#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"

#include "divfree_config.h"
#include "TPZMatDivFreeBubbles.h"
#include "TPZL2ProjectionCS.h"
#include "TPZCompElKernelHdiv.h"
#include "TPZCompElKernelHdivBC.h"
#include "TPZMixedDarcyFlowHybrid.h"
#include "TPZKernelHdivHybridizer.h"
#include "TPZKernelHdivUtils.h"
#include "TPZApproxSpaceKernelHdiv.h"


TPZMultiphysicsCompMesh *MultiphysicCMeshNew(int dim, int pOrder, std::set<int> &matIdVec, std::set<int> &matIdNeumann, TPZVec<TPZCompMesh *> &meshvector,TPZGeoMesh * gmesh);
TPZCompMesh *CMeshDivFreeBubbles(int dim, int pOrder, std::set<int> matIdVec, TPZGeoMesh *gmesh);
void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> &meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh);
void PrintResultsMultiphysicNew(int dim, TPZVec<TPZCompMesh *> &meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh);
void PrintResultsDivFreeBubbles(int dim, TPZLinearAnalysis &an);
void ComputeError(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);
void ComputeErrorHdiv(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

//Analytical solution
constexpr int solOrder{2};
auto exactSol = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];

    u[0] = x*x*x*y - y*y*y*x;
    gradU(0,0) = (3.*x*x*y - y*y*y);
    gradU(1,0) = (x*x*x - 3.*y*y*x);

    // u[0]= x*x-y*y;
    // gradU(0,0) = 2.*x;
    // gradU(1,0) = -2.*y;

    // const auto &d = 1.5; // distance betweel injection and production wells
    // u[0]=log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
    // gradU(0,0) = x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2));
    // gradU(1,0) = y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2));
};

auto exactSol2 = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];

    u[0] = x*x*x*y - y*y*y*x;
    gradU(0,0) = -(3.*x*x*y - y*y*y);
    gradU(1,0) = -(x*x*x - 3.*y*y*x);

    // u[0]= x*x-y*y;
    // gradU(0,0) = -2.*x;
    // gradU(1,0) = 2.*y;

    // const auto &d = 1.5; // distance betweel injection and production wells
    // u[0]=log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
    // gradU(0,0) = x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2));
    // gradU(1,0) = y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2));
};

enum EMatid  {ENone, EDomain, EBottom, ERight, ETop, ELeft, EPont, EWrap, EIntface, EPressureHyb};

int main(int argc, char* argv[])
{
    //dimension of the problem
    constexpr int dim{2};
    constexpr int pOrder{1};
      

#ifdef PZ_LOG
TPZLogger::InitializePZLOG();
#endif
    
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        stringtoint[2]["Surface"] = 1;
        stringtoint[1]["Bottom"] = 2;
        stringtoint[1]["Right"] = 3;
        stringtoint[1]["Top"] = 4;
        stringtoint[1]["Left"] = 5;
        stringtoint[0]["Point"] = 6;
        stringtoint[1]["Top2"] = 7;
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh4("../mesh/1element.msh",gmesh);
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    
    //.................................Hdiv.................................
    TPZCompMesh * cmeshflux = 0;
    TPZCompMesh * cmeshpressure = 0;
   
    //.........................Div Free Bubbles NEW.........................
    // {
        TPZKernelHdivUtils<STATE> util;
        TPZKernelHdivHybridizer hybridizer;
        //Insert here the BC material id's to be hybridized
        std::set<int> matBCHybrid={ERight};
        //Insert here the type of all boundary conditions
        std::set<int> matIDNeumann{ERight};
        std::set<int> matIDDirichlet{ETop,EBottom,ELeft};
        /// All bc's mat ID's
        std::set<int> matBC;
        std::set_union(matIDNeumann.begin(),matIDNeumann.end(),matIDDirichlet.begin(),matIDDirichlet.end(),std::inserter(matBC, matBC.begin()));

        /// Creates the approximation space - Set the type of domain hybridization
        TPZApproxSpaceKernelHdiv<STATE> createSpace(gmesh,TPZApproxSpaceKernelHdiv<STATE>::EDomainHybrid);

        //Setting material ids
        createSpace.fConfig.fDomain = EDomain;
        createSpace.SetPeriferalMaterialIds(EWrap,EPressureHyb,EIntface,EPont,matBCHybrid,matBC);
        createSpace.SetPOrder(pOrder+1);
        createSpace.Initialize();

        //Flux mesh
        TPZCompMesh * cmeshfluxNew = createSpace.CreateFluxCMesh();
        // util.PrintCMeshConnects(cmeshfluxNew);
        // std::string fluxFile = "FluxCMesh";
        // util.PrintCompMesh(cmeshfluxNew,fluxFile);

        //Pressure mesh
        TPZCompMesh * cmeshpressureNew = createSpace.CreatePressureCMesh();
        // util.PrintCMeshConnects(cmeshpressureNew);
        // std::string pressureFile = "PressureCMesh";
        // util.PrintCompMesh(cmeshpressureNew,pressureFile);

        //Multiphysics mesh
        TPZManVector< TPZCompMesh *, 2> meshvectorNew(2);
        meshvectorNew[0] = cmeshfluxNew;
        meshvectorNew[1] = cmeshpressureNew;      
        auto * cmeshNew = createSpace.CreateMultiphysicsCMesh(meshvectorNew,exactSol,matIDNeumann,matIDDirichlet);
        util.PrintCMeshConnects(cmeshNew);
        // Group and condense the elements
        // createSpace.Condense(cmeshNew);
        // std::string multiphysicsFile = "MultiPhysicsMeshNew";
        // util.PrintCompMesh(cmeshNew,multiphysicsFile);

        // Solve the problem
        TPZLinearAnalysis anNew(cmeshNew,false);
        createSpace.Solve(anNew, cmeshNew, true);

        //Print results
        PrintResultsMultiphysicNew(dim,meshvectorNew,anNew,cmeshNew);
        std::ofstream out4("mesh_MDFB.txt");
        anNew.Print("nothing",out4);

        std::ofstream anPostProcessFileMDFB("postprocessMDFB.txt");
        ComputeErrorHdiv(anNew,anPostProcessFileMDFB);
    // }

    //...........................Div Free Bubbles...........................
    //Creates DFB problem
    // TPZCompMesh * cmeshDFB = CMeshDivFreeBubbles(dim,pOrder,matIdVec,gmesh);

    // //Solve DFB
    // TPZLinearAnalysis anDFB(cmeshDFB,true);
    // SolveProblemDirect(anDFB,cmeshDFB);

    // // //Print results
    // PrintResultsDivFreeBubbles(dim,anDFB);
    // std::ofstream out("mesh.txt");
    // anDFB.Print("nothing",out);

    return 0;
}

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    //sets number of threads to be used by the solver
    constexpr int nThreads{0};
    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);

    //assembles the system
    an.Assemble();

    ///solves the system
    an.Solve();

    return;
}

void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    //sets number of threads to be used by the solver
    constexpr int nThreads{0};
    //defines storage scheme to be used for the FEM matrices
    //in this case, a symmetric skyline matrix is used
    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    ///Setting an iterative solver
    // TPZMatrixSolver<STATE> * precond = an.BuildPreconditioner(TPZAnalysis::EBlockJacobi , true);
    // TPZCopySolve<STATE> * precond = new TPZCopySolve<STATE>( matskl.Create() );  step.ShareMatrix( *precond );
    TPZStepSolver<STATE> * precond = new TPZStepSolver<STATE>( matskl.Create() ); step.ShareMatrix( *precond ); precond->SetJacobi(1, 0.0, 0);
    TPZStepSolver<STATE> jac;
    REAL tol = 1.e-30;
    jac.SetSSOR(1,1.1,0.,0);
    jac.ShareMatrix(step);
    step.SetGMRES(2000,2000, *precond, tol, 0);
    // step.SetCG(2000, *precond, tol, 0);
    an.SetSolver(step);

    //assembles the system
    an.Assemble();

    ///solves the system
    an.Solve();

    return;
}

void PrintResultsMultiphysicNew(int dim, TPZVec<TPZCompMesh *> &meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh)
{

    an.SetExact(exactSol,solOrder);

    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
    TPZManVector<std::string,10> scalnames(0), vecnames(2);


    // scalnames[0] = "Pressure";
    // scalnames[1] = "ExactPressure";
    vecnames[0]= "Flux";
    vecnames[1]= "ExactFlux";

    int div = 0;
    std::string plotfile = "solutionMDFB.vtk";
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    // Print mesh properties
    // std::ofstream out("mesh.txt");
    // an.Print("nothing",out);

    return;
}

void PrintResultsDivFreeBubbles(int dim, TPZLinearAnalysis &an)
{
    an.SetExact(exactSol,solOrder);
    TPZVec<std::string> vectorVars(1), scalarVars(1);
    scalarVars[0] = "Solution";
    vectorVars[0] = "Derivative";
    an.DefineGraphMesh(dim,scalarVars,vectorVars,"SolutionDFB.vtk");
    constexpr int resolution{0};
    an.PostProcess(resolution,dim);	

    return;
}

void ComputeError(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile)
{
    an.SetExact(exactSol2,solOrder);
    ///Calculating approximation error  
    TPZManVector<REAL,3> error;
    an.PostProcess(error,anPostProcessFile);
        
    std::cout << "\nApproximation error:\n";
    std::cout << "H1 Norm = " << error[0]<<'\n';
    std::cout << "L1 Norm = " << error[1]<<'\n'; 
    std::cout << "H1 Seminorm = " << error[2] << "\n\n";

    return;
}

void ComputeErrorHdiv(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile)
{
    an.SetExact(exactSol2,solOrder);
    ///Calculating approximation error  
    TPZManVector<REAL,5> error;

    auto cmeshNew = an.Mesh();
    int64_t nelem = cmeshNew->NElements();
    cmeshNew->LoadSolution(cmeshNew->Solution());
    cmeshNew->ExpandSolution();
    cmeshNew->ElementSolution().Redim(nelem, 5);

    an.PostProcessError(error);
        
    std::cout << "\nApproximation error:\n";
    std::cout << "H1 Norm = " << error[0]<<'\n';
    std::cout << "L1 Norm = " << error[1]<<'\n'; 
    std::cout << "H1 Seminorm = " << error[2]<<'\n'; 
    // std::cout << "error 4 = " << error[3]<<'\n'; 
    // std::cout << "error 5 = " << error[4] << "\n\n";

    return;
}
