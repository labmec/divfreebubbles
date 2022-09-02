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

#include "TPZTimer.h"
#include "fstream"
#include "TPZVTKGenerator.h"
#include "Poisson/TPZMatPoisson.h"


// The test function
void PGD(const int &xdiv, const int &pOrder);

TEST_CASE("Hybridization test")
{
    #define TEST
    const int pOrder = GENERATE(1);
    // const int pOrder = GENERATE(2,3,4,5);

    const int xdiv = 2;
    
    PGD(xdiv,pOrder); 
    
}

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

auto forcefunction = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    //Nabla u = 1
    u[0] = 1.;
};

void PGD(const int &xdiv, const int &pOrder)
{
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    int dim = 1;
    REAL minX = 0.;
    REAL maxX = 1.;
    int nMats = 2*dim+1;

    //all bcs share the same id
    constexpr bool createBoundEls{true};
    TPZVec<int> matIds(nMats,2);
    matIds[0] = 1;

    TPZGeoMesh* gmeshX = TPZGeoMeshTools::CreateGeoMesh1D(0.,2.,xdiv,matIds,createBoundEls);
    TPZGeoMesh* gmeshY = TPZGeoMeshTools::CreateGeoMesh1D(0.,1.,xdiv,matIds,createBoundEls);

    TPZCompMesh * cmeshX = new TPZCompMesh(gmeshX);
    TPZCompMesh * cmeshY = new TPZCompMesh(gmeshY);

    TPZMatPoisson<> *mat = new TPZMatPoisson<>(1,1);
    mat->SetForcingFunction(forcefunction,4);
    cmeshX->InsertMaterialObject(mat);
    cmeshY->InsertMaterialObject(mat);

    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);

    //Dirichlet Boundary Conditions
    TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, 2, 0, val1, val2);
    // BCond->SetForcingFunctionBC(exactSol,4);
    cmeshX->InsertMaterialObject(BCond);
    cmeshY->InsertMaterialObject(BCond);

    // Util for HDivKernel printing and solving
    TPZKernelHdivUtils<STATE> util;

    
    {
        const std::string plotfile = "myfile";//sem o .vtk no final
        constexpr int vtkRes{1};
    

        TPZVec<std::string> fields = {
        "Pressure",
        "ExactPressure",
        "Flux",
        "ExactFlux"};
        // auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

        // vtk.Do();
    }
    

}




























