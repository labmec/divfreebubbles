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
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZCompElKernelHDiv.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZKernelHdivUtils.h"
#include "TPZApproxSpaceKernelHdiv.h"
#include "TPZAnalyticSolution.h"

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

    // //Nabla u = 1
    // u[0] = 0.25*(x*x+y*y);
    // gradU(0,0) = -0.5*x;
    // gradU(1,0) = -0.5*y;

    //Nabla u = 0
    u[0] = x*x*x*y - y*y*y*x;
    gradU(0,0) = (3.*x*x*y - y*y*y);
    gradU(1,0) = (x*x*x - 3.*y*y*x);
    u[0] = x;
    gradU(0,0) = 1.;
    gradU(1,0) = 0.;

    REAL a1 = 1./4;
    REAL alpha = M_PI/2;
    u[0] = x*a1*cos(x*alpha)*cosh(y*alpha) + y*a1*sin(x*alpha)*sinh(y*alpha);
    gradU(0,0) = -a1*(cosh(alpha*y)*(cos(alpha*x) - alpha*x*sin(alpha*x)) + alpha*y*cos(alpha*x)*sinh(alpha*y));
    gradU(1,0) = -a1*(alpha*y*cosh(alpha*y)*sin(alpha*x) + (alpha*x*cos(alpha*x) + sin(alpha*x))*sinh(alpha*y));

};


enum EMatid  {ENone, EDomain, EBottom, ERight, ETop, ELeft, EPont, EWrap, EIntface, EPressureHyb};

int main(int argc, char* argv[])
{
    //dimension of the problem
    constexpr int dim{2};
    constexpr int pOrder{5};
      

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
        reader.GeometricGmshMesh("../mesh/1element.msh",gmesh);
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    

    //............................Div Free Bubbles............................
    TPZCompMesh * cmeshflux = 0;
    TPZCompMesh * cmeshpressure = 0;
   
    TPZKernelHdivUtils<STATE> util;

    //Insert here the BC material id's to be hybridized
    std::set<int> matBCHybrid={};
    //Insert here the type of all boundary conditions
    std::set<int> matIDNeumann{};
    std::set<int> matIDDirichlet{ERight,ETop,EBottom,ELeft};
    /// All bc's mat ID's
    std::set<int> matBC;
    std::set_union(matIDNeumann.begin(),matIDNeumann.end(),matIDDirichlet.begin(),matIDDirichlet.end(),std::inserter(matBC, matBC.begin()));

    /// Creates the approximation space - Set the type of domain hybridization
    TPZApproxSpaceKernelHdiv<STATE> createSpace(gmesh,
                                                TPZApproxSpaceKernelHdiv<STATE>::ENone,
                                                HDivFamily::EHDivConstant);

    //Setting material ids      
    createSpace.fConfig.fDomain = EDomain;
    createSpace.SetPeriferalMaterialIds(EWrap,EPressureHyb,EIntface,EPont,matBCHybrid,matBC);
    createSpace.SetPOrder(pOrder);
    createSpace.Initialize();
    // util.PrintGeoMesh(gmesh);

    //Flux mesh
    TPZCompMesh * cmeshfluxNew = createSpace.CreateFluxCMesh();
    // std::cout << "FLUX \n";
    // util.PrintCMeshConnects(cmeshfluxNew);
    // std::string fluxFile = "FluxCMesh";
    // util.PrintCompMesh(cmeshfluxNew,fluxFile);

    //Pressure mesh
    TPZCompMesh * cmeshpressureNew = createSpace.CreatePressureCMesh();
    // std::cout << "PRESSURE \n";
    // // util.PrintCMeshConnects(cmeshpressureNew);
    // std::string pressureFile = "PressureCMesh";
    // util.PrintCompMesh(cmeshpressureNew,pressureFile);

    TLaplaceExample1 LaplaceExact;
    // LaplaceExact.fExact = TLaplaceExample1::EHarmonic;
    LaplaceExact.fExact = TLaplaceExample1::EHarmonic2;
    //Multiphysics mesh
    TPZManVector< TPZCompMesh *, 2> meshvectorNew(2);
    meshvectorNew[0] = cmeshfluxNew;
    meshvectorNew[1] = cmeshpressureNew;       
    // auto * cmeshNew = createSpace.CreateMultiphysicsCMesh(meshvectorNew,LaplaceExact.ExactSolution(),matIDNeumann,matIDDirichlet);
    auto * cmeshNew = createSpace.CreateMultiphysicsCMesh(meshvectorNew,exactSol,matIDNeumann,matIDDirichlet);
    // std::cout << "MULTIPHYSICS \n";
    // util.PrintCMeshConnects(cmeshNew);
    // Group and condense the elements
    // createSpace.Condense(cmeshNew);
    std::string multiphysicsFile = "MultiPhysicsMeshNew";
    util.PrintCompMesh(cmeshNew,multiphysicsFile);
    std::cout << "Number of equations = " << cmeshNew->NEquations() << std::endl;
    // Solve the problem
    TPZLinearAnalysis anNew(cmeshNew,false);
    createSpace.Solve(anNew, cmeshNew, true, false);

    std::cout << "Number of equations = " << cmeshNew->NEquations() << std::endl;

    anNew.SetExact(exactSol,solOrder);
    // anNew.SetExact(LaplaceExact.ExactSolution());
    //Print results
    util.PrintResultsMultiphysics(meshvectorNew,anNew,cmeshNew);

    std::ofstream out4("mesh_MDFB.txt");
    anNew.Print("nothing",out4);
    std::ofstream anPostProcessFileMDFB("postprocessMDFB.txt");
    
    util.ComputeError(anNew,anPostProcessFileMDFB);
  
    return 0;
}

