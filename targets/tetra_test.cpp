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
#include <pzlog.h>
#include <pzfstrmatrix.h>
#include <TPZLinearAnalysis.h>
#include <pzstepsolver.h>
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

/**
   @brief Removes some equations associated with edges to ensure that
   the gradient of the lowest order H1 functions cannot be represented.
   @param[in] cmesh Computational mesh.
   @param[out] indices of all remaining equations.
   @return 0 if no errors were detected, 1 if a vertex was left untreated,
   2 if a vertex had all the adjacent edges removed.
*/
bool
FilterEdgeEquations(TPZAutoPointer<TPZCompMesh> cmesh,
                    TPZVec<int64_t> &activeEquations);


int CalcRank(TPZFMatrix<STATE> & S, const STATE tol);

int main()
{
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif
  constexpr int volId{1};
  constexpr int bcId{-1};

  //for now this should suffice
  const TPZManVector<int,3> nDivs = {1,1,1};
  
  
  auto gmesh = CreateGeoMesh(nDivs, volId, bcId);

  constexpr int pOrder{1};

  auto cmesh = CreateCompMesh(gmesh, pOrder, volId, bcId);


  constexpr bool reorderEqs{true};
  TPZLinearAnalysis an(cmesh, reorderEqs);
  
  TPZAutoPointer<TPZStructMatrix> strmtrx =
    new TPZFStructMatrix<STATE>(cmesh);
  
  constexpr int nThreads{0};
  strmtrx->SetNumThreads(nThreads);

  TPZVec<int64_t> activeEqs;
  
  if(FilterEdgeEquations(cmesh, activeEqs)){
    return 1;
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

  STATE tol = 1e-5;
  const auto rank = CalcRank(*mat, tol);

  std::cout<<"neq: "<<neqs<<" rank: "<<rank<<" diff: "<<neqs - rank<<std::endl;
  return 0;
}



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


bool
FilterEdgeEquations(TPZAutoPointer<TPZCompMesh> cmesh,
                    TPZVec<int64_t> &activeEquations)
{

  /**TODO:
     1: for each vertex: at least one edge must be removed

     2: test for p = 1.

     3: create the class for filtering the functions in the master el.
     3.1: phi_f^{e,n}: remove one for each face.

     4: test for p = 2.

     5: think about the remaining functions.
   */
  const auto gmesh = cmesh->Reference();
  
  const auto nnodes = gmesh->NNodes();

  /**
     The i-th position contains the indices of all the 
     connects associated with the edges adjacent to the i-th node
   */
  TPZVec<std::set<int>> vertex_edge_connects(nnodes);

  /**
     The i-th is true if the i-th node has already been dealt with.
  */
  TPZVec<bool> done_vertices(nnodes, false);
  /**
     Contains all the connects marked for removal. It is expected
     that all the connects are associated with edges.
   */
  std::set<int> removed_connects;

  constexpr int edgeDim{1};
  
  for(auto gel : gmesh->ElementVec()){
    const auto nEdges = gel->NSides(edgeDim);
    const auto firstEdge = gel->FirstSide(edgeDim);
    const auto firstFace = firstEdge + nEdges;
    
    for(auto ie = firstEdge; ie < firstFace; ie++){
      TPZGeoElSide edge(gel,ie);
      const auto con = edge.Reference().ConnectIndex();
      //check if edge has been treated already
      if (removed_connects.find(con) != removed_connects.end()) {
        continue;
      }
      
      const auto v1 = edge.SideNodeIndex(0);
      const auto v2 = edge.SideNodeIndex(1);

      vertex_edge_connects[v1].insert(con);
      vertex_edge_connects[v2].insert(con);
      if(!done_vertices[v1] || !done_vertices[v2]){
        removed_connects.insert(edge.Reference().ConnectIndex());
        done_vertices[v1] = true;
        done_vertices[v2] = true;
      }
    }
  }
  bool check{true};
  for(auto v : done_vertices){
    check = check && v;
    if(!v){
      break;
    }
  }

  if(!check){
    std::cout<<__PRETTY_FUNCTION__
             <<"\nError: could not treat all vertices"<<std::endl;
    return 1;
  }

  check = true;
  for(auto iv = 0; iv < nnodes; iv++){
    const auto all_edges = vertex_edge_connects[iv];
    bool local_check{false};
    for(auto edge: all_edges){
      if(removed_connects.find(edge) == removed_connects.end()){
        local_check = true;
        break;
      }
    }
    check = local_check && check;
  }
  
  if(!check){
    std::cout<<__PRETTY_FUNCTION__
             <<"\nError: a vertex had all the edges removed"<<std::endl;
    return 2;
  }
  
  activeEquations.Resize(0);
  for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
    if (removed_connects.find(iCon) == removed_connects.end()) {
      auto &con = cmesh->ConnectVec()[iCon];
      if (con.HasDependency()){
        continue;
      }
      const auto seqnum = con.SequenceNumber();
      const auto pos = cmesh->Block().Position(seqnum);
      const auto blocksize = cmesh->Block().Size(seqnum);
      if (blocksize == 0){
        continue;
      }
      
      const auto vs = activeEquations.size();
      activeEquations.Resize(vs + blocksize);
      for (auto ieq = 0; ieq < blocksize; ieq++) {
        activeEquations[vs + ieq] = pos + ieq;
      }
    }
  }
  
  return 0;
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