#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZMultiphysicsCompMesh.h>
#include <TPZLinearAnalysis.h>
#include <TPZGmshReader.h>
#include <TPZVTKGeoMesh.h>
#include <TPZCompElDisc.h>
#include <TPZNullMaterial.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>// for Hdiv problem
#include <Poisson/TPZMatPoisson.h>
#include <pzbuildmultiphysicsmesh.h>
#include <TPZNullMaterialCS.h>
#include <pzlog.h>
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZLagrangeMultiplierCS.h>
#include <pzshapelinear.h>
#include <pzshapepoint.h>
#include <pzshapequad.h>
#include <pzshapetriang.h>

#include "divfree_config.h"
#include "TPZMatDivFreeBubbles.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZCompElKernelHDiv.h"
#include "TPZCompElKernelHDivBC.h"
#include "TPZKernelHdivHybridizer.h"
#include "TPZKernelHdivUtils.h"
#include "TPZApproxSpaceKernelHdiv.h"

//Analytical solution
constexpr int solOrder{2};
auto exactSol = [](const TPZVec<REAL> &loc,
  TPZVec<STATE>&u,
  TPZFMatrix<STATE>&gradU){
  const auto &x=loc[0];
  const auto &y=loc[1];
  const auto &d = 1.; // distance between injection and production wells
  u[0]= log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
  gradU(0,0) = (x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2)));
  gradU(1,0) = (y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2)));
};

//Analytical solution
auto exactSolError = [](const TPZVec<REAL> &loc,
  TPZVec<STATE>&u,
  TPZFMatrix<STATE>&gradU){
  const auto &x=loc[0];
  const auto &y=loc[1];
  const auto &d = 1.; // distance between injection and production wells
  u[0]= log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
  gradU(0,0) = -((x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2))));
  gradU(1,0) = -((y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2))));
};

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

enum EMatid {ENone, EDomain, EInjection, EProduction, EBottom,  ETop, ELeft, ERight, EPont, EWrap, EIntface, EPressureHyb};

int main(int argc, char* argv[]){
    //dimension of the problem
    constexpr int dim{2};
    constexpr int pOrder{2};

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
        stringtoint[1]["InjectionWell"] = 2;
        stringtoint[1]["ProductionWell"] = 3;
        stringtoint[1]["BottomLine"] = 4;
        stringtoint[1]["TopLine"] = 5;
        stringtoint[1]["LeftLine"] = 6;
        stringtoint[1]["RightLine"] = 7;
        stringtoint[0]["Point"] = 8;
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh(string(MESHDIR)+"newMesh.msh",gmesh);
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
  
    //.................................Hdiv.................................
    TPZCompMesh * cmeshflux = 0;
    TPZCompMesh * cmeshpressure = 0;
 
    TPZKernelHdivUtils<STATE> util;

    //Insert here the BC material id's to be hybridized
    std::set<int> matBCHybrid={EInjection,EProduction};
    //Insert here the type of all boundary conditions
    std::set<int> matIDNeumann{EInjection,EProduction};
    std::set<int> matIDDirichlet{EBottom,ETop,ELeft,ERight};
    /// All bc's mat ID's
    std::set<int> matBC;
    std::set_union(matIDNeumann.begin(),matIDNeumann.end(),matIDDirichlet.begin(),matIDDirichlet.end(),std::inserter(matBC, matBC.begin()));

    /// Creates the approximation space - Set the type of domain hybridization
    TPZApproxSpaceKernelHdiv<STATE> createSpace(gmesh,TPZApproxSpaceKernelHdiv<STATE>::ENone);

    //Setting material ids
    createSpace.fConfig.fDomain = EDomain;
    createSpace.SetPeriferalMaterialIds(EWrap,EPressureHyb,EIntface,EPont,matBCHybrid,matBC);
    createSpace.SetPOrder(pOrder+1);
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
    // util.PrintCMeshConnects(cmeshpressureNew);
    // std::string pressureFile = "PressureCMesh";
    // util.PrintCompMesh(cmeshpressureNew,pressureFile);

    //Multiphysics mesh
    TPZManVector< TPZCompMesh *, 2> meshvectorNew(2);
    meshvectorNew[0] = cmeshfluxNew;
    meshvectorNew[1] = cmeshpressureNew;      
    auto * cmeshNew = createSpace.CreateMultiphysicsCMesh(meshvectorNew,exactSol,matIDNeumann,matIDDirichlet);
    // std::cout << "MULTIPHYSICS \n";
    // util.PrintCMeshConnects(cmeshNew);
    // Group and condense the elements
    std::cout << "Number of equations1 = " << cmeshNew->NEquations() << std::endl;

    createSpace.Condense(cmeshNew);
    // std::string multiphysicsFile = "MultiPhysicsMeshNew";
    // util.PrintCompMesh(cmeshNew,multiphysicsFile);

    // Solve the problem
    TPZLinearAnalysis anNew(cmeshNew,false);
    createSpace.Solve(anNew, cmeshNew, true, false);

    std::cout << "Number of equations2 = " << cmeshNew->NEquations() << std::endl;
    std::cout << "Number of equations3 = " << anNew.Mesh()->NEquations() << std::endl;
    anNew.SetExact(exactSol,solOrder);
    //Print results
    util.PrintResultsMultiphysics(meshvectorNew,anNew,cmeshNew);

    anNew.SetExact(exactSol,solOrder);

    std::ofstream out4("mesh_MDFB.txt");
    anNew.Print("nothing",out4);
    std::ofstream anPostProcessFileMDFB("postprocessMDFB.txt");
    
    util.ComputeError(anNew,anPostProcessFileMDFB);

    return 0;
}
