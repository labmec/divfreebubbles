#include <TPZGeoMeshTools.h>
#include "TPZHDivApproxSpaceCreator.h"
#include "TPZKernelHdivUtils.h"
#include "TPZAnalyticSolution.h"
#include <TPZGmshReader.h>
#include "TPZCompMeshTools.h"
#include <TPZGenGrid3D.h>
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
#include "TPZHDivApproxCreator.h"
#include "Elasticity/TPZMixedElasticityND.h"

#include "divfree_config.h"

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

// GLOBAL VARIABLES FOR VARIATION OF ELASTICITY PROPERTIES
// Due to the matlab code that generates the files we always have to use ndiv = ndiv+1
//constexpr int Globnx{129}, Globny{129}, Globnz{65};
//constexpr REAL Globpartsize{78.125}; // 10000/Globnx

constexpr int Globnx{65}, Globny{65}, Globnz{33};
constexpr REAL Globpartsize{156.25}; // 10000/Globnx

//constexpr int Globnx{17}, Globny{17}, Globnz{9};
//constexpr REAL Globpartsize{625.}; // 10000/Globnx

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

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

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

enum EMatid  {ENone, EDomain, EBoundary, EPont, EWrap, EIntface, EPressureHyb, ENeumannZero};
using namespace std;

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

template<class tshape> TPZGeoMesh* CreateGeoMesh(TPZVec<int> &nDivs);
void InsertMaterials(TPZHDivApproxCreator &hdivCreator, TPZVec<int> &nDivs);

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
int main() {
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    
    std::cout << "\n----------- Starting simulation -----------" << std::endl;
    TPZSimpleTimer totaltime;
    
    int DIM = 3;
    TPZManVector<int,3> nDivs = {Globnx-1,Globny-1,Globnz-1}; // one uniform refinement is performed

    // Creates a geometric mesh
    auto gmesh = CreateGeoMesh<pzshape::TPZShapeTetra>(nDivs);
    std::cout << "Number of geo elements in mesh = " << gmesh->NElements() << std::endl;
    
    // Util for HDivKernel printing and solving
    TPZKernelHdivUtils<STATE> util;
    
    // Set simulation parameters
    TPZHDivApproxCreator hdivCreator(gmesh);
    hdivCreator.HdivFamily() = HDivFamily::EHDivStandard;
    hdivCreator.ProbType() = ProblemType::EElastic;
    hdivCreator.IsRigidBodySpaces() = false;
    hdivCreator.SetDefaultOrder(1);
    hdivCreator.SetExtraInternalOrder(0);
    hdivCreator.SetShouldCondense(true);
    hdivCreator.HybridType() = HybridizationType::EStandard;

    // Prints gmesh
    std::string vtk_name = "geomesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

    // Creates the materials, elasticity function and bcs and adds to hdivApproxCreator
    InsertMaterials(hdivCreator,nDivs);

    // Multiphysics mesh
    TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();

    // hdivCreator.PrintAllMeshes(cmesh);
  
    // Number of equations without condense elements
    const int nEquationsFull = cmesh->NEquations();
    std::cout << "Number of equations condensed = " << nEquationsFull << std::endl;
    std::cout << "Number of equations = " << cmesh->Solution().Rows() << std::endl;

    //Create analysis environment and solve
    TPZLinearAnalysis an(cmesh,true);
    bool filter = false;bool domainhybr=false;
    const int nthreads = 64;
    util.SolveProblemCholesky(an,cmesh,filter,domainhybr,nthreads);

    // Print results
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), cmesh);
    const std::string plotfile = "PostProcess"; //sem o .vtk no final
    constexpr int vtkRes{0};
    TPZManVector<std::string,6> fields = {"Displacement","SigmaX","SigmaY","TauXY","Young_Modulus","Poisson"};
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

    vtk.Do();

    
    std::cout << "Total time: " << totaltime.ReturnTimeDouble()/1000 << " seconds" << std::endl;

}
// ----------------- End of Main -----------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

template <class tshape> TPZGeoMesh* CreateGeoMesh(TPZVec<int> &nDivs) {
    
    std::cout << "\n----------- Creating gmesh -----------" << std::endl;
    
    constexpr int ndivInternal = 1;
    TPZManVector<int,10> nDivsSkel(nDivs.size());
    for (int i = 0; i < 3; i++) {
        nDivsSkel[i] = nDivs[i] / (ndivInternal*2);
    }
    
    // ----- Create Geo Mesh -----
    const TPZManVector<REAL,3> minX = {0.,0.,0.};
    const TPZManVector<REAL,3> maxX = {10000.,10000.,5000.};
    
    TPZGenGrid3D gen3d(minX,maxX,nDivsSkel,MMeshType::ETetrahedral);
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gmesh = gen3d.BuildVolumetricElements(EDomain);

    gmesh = gen3d.BuildBoundaryElements(EBoundary,ENeumannZero,ENeumannZero,ENeumannZero,ENeumannZero,ENeumannZero);
    
    TPZCheckGeom check(gmesh);
    check.UniformRefine(ndivInternal);

    std::cout << "\n----------- Finished creating gmesh -----------" << std::endl;
    
    return gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void InsertMaterials(TPZHDivApproxCreator &approxCreator, TPZVec<int> &nDivs) {
    
    if (!approxCreator.GeoMesh()) {
        cout << "\nError! Please set the geomesh before inserting materials" << endl;
        DebugStop();
    }
    std::string sourcepath(MESHDIR);
    sourcepath += "../";
    std::cout << "\npath = " << sourcepath << std::endl;
    
    const int dim = approxCreator.GeoMesh()->Dimension();
    
    // 3D volume material
    REAL E = 1., nu = 0., fx = 0., fy = 0.;
    const int plain = 0.; //* @param plainstress = 1 \f$ indicates use of plainstress
    TPZMixedElasticityND* matelas = new TPZMixedElasticityND(EDomain, E, nu, fx, fy, plain, dim);
    approxCreator.InsertMaterialObject(matelas);
    
    // Boundary conditions
    constexpr int dirType = 0, neuType = 1;
    TPZFMatrix<STATE> val1(dim,dim,0.);
    TPZManVector<STATE> val2(dim,0.);
    TPZBndCondT<STATE> *BCond1 = matelas->CreateBC(matelas, ENeumannZero, neuType, val1, val2);
    TPZBndCondT<STATE> *BCond2 = matelas->CreateBC(matelas, EBoundary, dirType, val1, val2);
    approxCreator.InsertMaterialObject(BCond1);
    approxCreator.InsertMaterialObject(BCond2);
        
    matelas->SetBodyForce(0, 0, -9.81/1.e6); // Gravity forces only

    // ----------------------- Getting E,nu data --------------------------
    // --------------------------------------------------------------------
    TPZManVector<TPZFMatrix<STATE>,Globny> edata(Globny), nudata(Globny);
    for (int iy = 0; iy < Globny; iy++) {
        edata[iy].Resize(Globnz, Globnx);
        nudata[iy].Resize(Globnz, Globnx);
    }
    REAL tempE = 0., tempNu = 0.;;

    const std::string foldername = sourcepath + "data/" + to_string(Globnx-1);
    std::string basee = foldername + "/e_", basenu = foldername + "/nu_";
    for (int iy = 0; iy < Globny; iy++) {
        std::string ename = basee + to_string(iy+1) + ".txt";
        std::string nuname = basenu + to_string(iy+1) + ".txt";
        std::ifstream inE(ename), inNu(nuname);
        if(!inE || !inNu) DebugStop(); // Files must exist
        for(int iz = 0; iz < Globnz ; iz++) {
            for(int ix = 0 ; ix < Globnx ; ix++){
                inE >> tempE;
                inNu >> tempNu;
                edata[iy](iz,ix) = tempE;
                nudata[iy](iz,ix) = tempNu;
            }
        }
    }
    
    auto ConstLawFunctionLambda = [edata,nudata] (const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv) {
        int rounded_x = static_cast<int>(floor(x[0]/Globpartsize));
        int rounded_y = static_cast<int>(floor(x[1]/Globpartsize));
        int rounded_z = static_cast<int>(floor(x[2]/Globpartsize));
        STATE eval = edata[rounded_y](rounded_z,rounded_x)/3.e9;
        STATE nuval = nudata[rounded_y](rounded_z,rounded_x);
        result[0] = eval;
        result[1] = nuval;
    };
    
    std::function<void(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)> myfunc = ConstLawFunctionLambda;
    
    matelas->SetElasticityFunction(myfunc);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
