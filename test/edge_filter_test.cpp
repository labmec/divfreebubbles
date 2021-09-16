/*
  This unit test verifies that for a given mesh, upon filtering,
  we are able to generate a HCurl approximation space such that the curl
  of our basis functions are linearly independent.

  The elements used in this target are derived from standard HCurl elements,
  and some functions are already filtered in the master element.

  Some of the edge functions, however, need to be globally filtered.
  
*/
#include <catch2/catch.hpp>
#include <TPZGeoMeshTools.h>
#include <pzcmesh.h>
#include <pzshapetriang.h>
#include <pzshapetetra.h>
#include <pzfstrmatrix.h>
#include <TPZLinearAnalysis.h>
#include <pzstepsolver.h>
#include "TPZMatCurlDotCurl.h"
#include "TPZCompElHCurlNoGrads.h"
#include "TPZKernelHdivUtils.h"
#include <TPZGmshReader.h>

/**
   @brief Creates a geometric mesh with elements of a given type on a unit cube.
   @param[in] meshType element type to be created.
   @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
TPZGeoMesh*
CreateGeoMesh(const MMeshType meshType, const TPZVec<int> &nDivs,
              const int volId, const int bcId);

/**
   @brief Creates the computational mesh.
   @param[in] gmesh Geometrical mesh.
   @param[in] pOrder Polynomial order of the elements.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
TPZAutoPointer<TPZCompMesh>
CreateCompMesh(TPZAutoPointer<TPZGeoMesh> gmesh, const int pOrder,
               const int volId, const int bcId);

/**
   @brief Reads the test mesh from gmsh
   @param[in] file_name the .msh mesh file.
*/
TPZGeoMesh*
ReadMeshFromGmsh(std::string file_name);


int CalcRank(TPZFMatrix<STATE> & S, const STATE tol);

void TestEdgeFiltering(TPZGeoMesh* gmesh, const int &volId, const int &bcId);

TEST_CASE("Edge filtering", "hcurl_els")
{
    const int xdiv = GENERATE(1,2,3);
    const int ydiv = GENERATE(1,2,3);
    const int zdiv = GENERATE(1,2,3);
    const MMeshType meshType = MMeshType::ETetrahedral;

    constexpr int volId{1};
    constexpr int bcId{-1};

    //for now this should suffice
    const TPZManVector<int,3> nDivs = {xdiv,ydiv,zdiv};
    
    auto gmesh = CreateGeoMesh(meshType, nDivs, volId, bcId);

    TestEdgeFiltering(gmesh,volId,bcId);
}

TEST_CASE("1 tetrahedron Gmsh")
{
    auto gmesh = ReadMeshFromGmsh("../mesh/1tetra.msh");
    TestEdgeFiltering(gmesh,1,2);
}

TEST_CASE("Cube Gmsh")
{
    auto gmesh = ReadMeshFromGmsh("../mesh/cube.msh");
    TestEdgeFiltering(gmesh,1,2);
}

void TestEdgeFiltering(TPZGeoMesh* gmesh, const int &volId, const int &bcId)
{
    constexpr int pOrder{1};

    auto cmesh = CreateCompMesh(gmesh, pOrder, volId, bcId);


    constexpr bool reorderEqs{true};
    TPZLinearAnalysis an(cmesh, reorderEqs);
    
    TPZAutoPointer<TPZStructMatrix> strmtrx =
        new TPZFStructMatrix<STATE>(cmesh);
    
    constexpr int nThreads{0};
    strmtrx->SetNumThreads(nThreads);

    TPZVec<int64_t> activeEqs;
    
    TPZKernelHdivUtils<STATE> util;

    if(util.FilterEdgeEquations(cmesh, activeEqs)){
        return;
    }

    const int neqs = activeEqs.size();
    
    strmtrx->EquationFilter().SetActiveEquations(activeEqs);

    an.SetStructuralMatrix(strmtrx);

    TPZStepSolver<STATE> step;

    step.SetDirect(ECholesky);
    an.SetSolver(step);

    an.Assemble();

    auto mat =
        TPZAutoPointerDynamicCast<TPZFMatrix<STATE>>(an.MatrixSolver<STATE>().Matrix());

    constexpr STATE tol = 1e-8;
    const auto rank = CalcRank(*mat, tol);
    CAPTURE(neqs,rank,neqs-rank);
    CHECK(rank == neqs);
}



TPZGeoMesh*
CreateGeoMesh(const MMeshType meshType, const TPZVec<int> &nDivs,
              const int volId, const int bcId)
{
    constexpr int dim{3};
    const TPZManVector<REAL,3> minX = {0,0,0};
    const TPZManVector<REAL,3> maxX = {1,1,1};
    constexpr bool createBoundEls{true};

    //all bcs share the same id
    const TPZManVector<int,7> matIds(7,bcId);
    matIds[0] = volId;
    
    TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
                        matIds, nDivs, meshType,createBoundEls);

    return gmesh;
}

TPZAutoPointer<TPZCompMesh>
CreateCompMesh(TPZAutoPointer<TPZGeoMesh> gmesh, const int pOrder,
               const int volId, const int bcId)
{
    constexpr bool isComplex{false};
    constexpr int dim{3};
    
    TPZAutoPointer<TPZCompMesh> cmesh =
        new TPZCompMesh(gmesh, isComplex);

    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);

    //insert volumetric material
    auto volMat = new TPZMatCurlDotCurl(volId);
    cmesh->InsertMaterialObject(volMat);
    //insert boundary material
    const int bcType = 0;//dirichlet
    TPZFNMatrix<1, STATE> val1(1, 1, 1);
    TPZManVector<STATE,1> val2(1, 0.);
    auto bcMat = volMat->CreateBC(volMat, bcId, bcType, val1, val2);
    cmesh->InsertMaterialObject(bcMat);
    
    //Creates computational elements
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) DebugStop();
            gel->ResetReference();   
            const MElementType type = gel->Type();
            const auto matid = gel->MaterialId();
            int64_t index;
            switch(type){
                case ETriangle:
                    new TPZCompElHCurlNoGrads<pzshape::TPZShapeTriang>(*cmesh,gel,index);
                    break;
                case ETetraedro:
                    new TPZCompElHCurlNoGrads<pzshape::TPZShapeTetra>(*cmesh,gel,index);
                    break;
                default:
                    const auto elName =  MElementType_Name(type);
                    CAPTURE(elName);
                    CHECK(false);
                    PZError<<__PRETTY_FUNCTION__
                            <<"\n type not yet supported. Aborting..."<<std::endl;
                    DebugStop();
        }
    }
    cmesh->SetAllCreateFunctionsHCurl();
    cmesh->AutoBuild();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}



int CalcRank(TPZFMatrix<STATE> & mat, const STATE tol){

    TPZFMatrix<STATE> S;
    TPZFMatrix<STATE> Udummy, VTdummy;
    mat.SVD(Udummy, S, VTdummy, 'N', 'N');

    int rank = 0;
    const int dimMat = S.Rows();
    for (int i = 0; i < dimMat; i++) {
        rank += S.GetVal(i, 0) > tol ? 1 : 0;
    }
    return rank;
};


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
        reader.GeometricGmshMesh4(file_name,gmesh);
    }

    return gmesh;
}