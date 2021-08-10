/*
  This target verifies that for a given tetrahedral mesh, upon filtering,
  we are able to generate a HCurl approximation space such that the curl
  of our basis functions are linearly independent.

  The elements used in this target are derived from standard HCurl elements,
  and some functions are already filtered in the master element.

  Some of the edge functions, however, need to be globally filtered.
  
*/

#include <TPZGeoMeshTools.h>
#include <pzcmesh.h>

#include "TPZMatCurlDotCurl.h"

/**
   @brief Creates a geometric mesh with tetrahedral elements
   on a unit cube.
   @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
   @param[in] volId Material identifier for the volumetric region.
   @param[in] bcId Material identifier for the boundary.
*/
TPZAutoPointer<TPZGeoMesh>
CreateGeoMesh(const TPZVec<int> &nDivs, const int volId, const int bcId);

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


TPZAutoPointer<TPZGeoMesh>
CreateGeoMesh(const TPZVec<int> &nDivs, const int volId, const int bcId)
{
  constexpr int dim{3};
  const TPZManVector<REAL,3> minX = {0,0,0};
  const TPZManVector<REAL,3> maxX = {1,1,1};
  constexpr MMeshType meshType{MMeshType::ETetrahedral};
  constexpr bool createBoundEls{true};

  //all bcs share the same id
  const TPZManVector<int,7> matIds(7,bcId);
  matIds[0] = volId;
  
  TPZAutoPointer<TPZGeoMesh> gmesh = 
    TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
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
  
  //TODO: change to the correct approx space
  cmesh->SetAllCreateFunctionsHCurl();
  cmesh->AutoBuild();
  cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;
}

int main()
{
  constexpr int volId{1};
  constexpr int bcId{-1};

  //for now this should suffice
  const TPZManVector<int,3> nDivs = {1,1,1};
  
  
  auto gmesh = CreateGeoMesh(nDivs, volId, bcId);

  constexpr int pOrder{1};

  auto cmesh = CreateCompMesh(gmesh, pOrder, volId, bcId);

  /**TODO:
     1: for each vertex: filter functions of ONE edge. 
     1.1: should we eliminate the connect itself,
     setting its order to zero? 
     or just gather the appropriate indices
     and then use the equation filter?

     2: test for p = 1.

     3: create the class for filtering the functions in the master el.
     3.1: phi_f^{e,n}: remove one for each face.

     4: test for p = 2.

     5: think about the remaining functions.
   */
  
}