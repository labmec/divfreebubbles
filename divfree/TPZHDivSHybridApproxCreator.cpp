#include "TPZHDivSHybridApproxCreator.h"
#include "pzgeoelbc.h"
#include "TPZNullMaterial.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzintel.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("TPZHDivSHybridApproxCreator");
#endif

TPZHDivSHybridApproxCreator::TPZHDivSHybridApproxCreator(TPZGeoMesh *gmesh)
    : TPZHDivApproxCreator(gmesh)
{
}

void TPZHDivSHybridApproxCreator::AddHybridizationGeoElements(){
    /// For fHybridType == HybridizationType::EStandard, wrap, interface and lagrange geometric elements are created;
    ///     if fHybridizeBCLevel == 0, the boundary is not hybridized, if fHybridizeBCLevel == 1, it is.
    /// For fHybridType == HybridizationType::EStandardSquared, wrap, interface, lagrange, second interface and second lagrange
    /// geometric elements are created;
    ///     BC can be hybridized 0, 1 or 2 times in this set up;
    /// fHybridType == HybridizationType::Semi has yet to be implemented.

    if(fHybridType == HybridizationType::ENone) {
        std::cout << "Hybridization type is none, this method should not be called\n";
        DebugStop();
    }
    if(fHybridizationData.fWrapMatId == -123456){
        std::cout << "\nERROR! Please call TPZApproxCreator::ComputePeriferalMaterialIds() before TPZApproxCreator::AddHybridizationGeoElements()" << std::endl;
        DebugStop();
    }

#ifdef PZ_LOG
    std::map<int,int> numcreated;
#endif
    // create the wrap and interface geometric elements
    // the creation of wrap and interface needs to happen in separate loops when there are hanging nodes
    int64_t nel = fGeoMesh->NElements();
    int dim = fGeoMesh->Dimension();
    std::set<int> bcMatIds = GetBCMatIds();
    // Wrap and interface creation for standard hyb.
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel || gel->HasSubElement() || gel->Dimension() != dim) continue;
        int nsides = gel->NSides();
        int side = gel->FirstSide(dim-1);
        // loop over the sides of dimension dim-1
        for(; side < nsides-1; side++)
        {
            TPZGeoElSide gelside(gel,side);
            // we want to create side elements of type
            // first fMatWrapId
            TPZGeoElSide neighbour = gelside.Neighbour();
#ifdef PZDEBUG
            {
                int neighMatId = neighbour.Element()->MaterialId();
                if(neighMatId == fHybridizationData.fWrapMatId)
                {
                    std::cout << __PRETTY_FUNCTION__ << " should be called only once!\n";
                    DebugStop();
                }
            }
#endif
            // if there is a neighbour that is a boundary condition
            // do not create the wrap layers
            bool hasBCNeighbour = gelside.HasNeighbour(bcMatIds);
            if(hasBCNeighbour && fHybridizationData.fHybridizeBCLevel == 0)
            {
                // no interface will be created between the element and a flux space
                continue;
            }
            if(gel->NormalOrientation(side) == -1)
            {
              continue; // we only create the wrap and interface elements for one of the two neighbouring elements
            }
            TPZGeoElBC(gelside, fHybridizationData.fWrapMatId);
#ifdef PZ_LOG
            numcreated[fHybridizationData.fWrapMatId]++;
#endif
#ifdef PZDEBUG
            neighbour = gelside.Neighbour();
            if(neighbour.Element()->MaterialId() != fHybridizationData.fWrapMatId)
            {
                DebugStop();
            }
#endif
        }
    }
    
  if(fHybridType == HybridizationType::EStandardSquared) {
    DebugStop(); // this is not implemented yet
  }
}

TPZCompMesh * TPZHDivSHybridApproxCreator::CreateHDivSpace(){
  auto hybrid = fHybridType;
  fHybridType = HybridizationType::ENone; // we do not want to create the hybridization geometric elements in the HDiv space
  // @TODO this will not create the wrap hdiv bound elements
  auto cmesh = TPZHDivApproxCreator::CreateHDivSpace();
  {
    int dim = cmesh->Dimension();
    TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial<STATE>(fHybridizationData.fWrapMatId,dim-1,dim);
    cmesh->InsertMaterialObject(nullmat);
    std::set<int> matids = {fHybridizationData.fWrapMatId};
    cmesh->AutoBuild(matids);
  }
  fHybridType = hybrid; // restore the original hybridization type
  return cmesh;
}

    /// Insert interface periferal material objects related to geometric objects created during hybridization
void TPZHDivSHybridApproxCreator::InsertInterfaceMaterialObjects(TPZMultiphysicsCompMesh *mphys){}

    /// Create interface elements on hybridizes spaces
    /// @param mphys multiphysics compmesh
void TPZHDivSHybridApproxCreator::AddInterfaceComputationalElements(TPZMultiphysicsCompMesh *mfmesh) {
  // return;
  int64_t nel = mfmesh->NElements();
  for(int64_t el = 0; el<nel; el++)
  {
    TPZCompEl *cel = mfmesh->Element(el);
    if(!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    if(!gel) continue;
    int matid = gel->MaterialId();
    if(matid != fHybridizationData.fWrapMatId)
    {
      continue;
    }
    TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if(!mcel)
    {
      std::cout << "TPZHDivSHybridApproxCreator::AddInterfaceComputationalElements error. The element is not a multiphysics element\n";
      DebugStop();
    }
    TPZInterpolatedElement *wrapel = dynamic_cast<TPZInterpolatedElement *>(mcel->Element(0));
    if(!wrapel)
    {
      std::cout << "TPZHDivSHybridApproxCreator::AddInterfaceComputationalElements error. The first element of the multiphysics element is not an interpolated element\n";
      DebugStop();
    }
    if(mcel->NConnects() != 1)
    {
      std::cout << "TPZHDivSHybridApproxCreator::AddInterfaceComputationalElements error. The multiphysics element has more than one connect\n";
      DebugStop();
    }
    TPZConnect &cleft = mcel->Connect(0);
    // restrain the connect of the left lagrange multiplier to a connect representing the low order flux space and the remainder of dofs of the lagrange multiplier connect
    int64_t newind1 = mfmesh->AllocateNewConnect(cleft);
    TPZConnect &cnew1 = mfmesh->ConnectVec()[newind1];
  
    int64_t newind2 = mfmesh->AllocateNewConnect(cleft);
    TPZConnect &cnew2 = mfmesh->ConnectVec()[newind2];

    // std::cout << "cleftindex " << cleftindex <<
    // " crightindex " << crightindex << " newind1 " << newind1 << " newind2 " << newind2 << std::endl;
    cnew1.SetNState(1);
    cnew1.SetNShape(3);
    int64_t seq1 = cnew1.SequenceNumber();
    mfmesh->Block().Set(seq1, 3);
    cnew2.SetNState(1);
    cnew2.SetNShape(cleft.NDof() - 3);
    int64_t seq2 = cnew2.SequenceNumber();
    mfmesh->Block().Set(seq2, cleft.NDof() - 3);

    RestraintConnect(wrapel, mcel, newind1, newind2);
    TPZConnect &c = mcel->Connect(0);
    TPZConnect::TPZDependBase *dep = c.FirstDepend();
    while(dep)
    {
      TPZConnect::TPZDepend<STATE> *depcast = dynamic_cast<TPZConnect::TPZDepend<STATE> *>(dep);
      if(!depcast)
      {
        std::cout << "TPZHDivSHybridApproxCreator::AddInterfaceComputationalElements error. The depend is not of type TPZConnect::TPZDepend<STATE>\n";
        DebugStop();
      }
      TPZFMatrix<STATE> depmat = depcast->GetDepMatrix();
      int64_t deprows = depmat.Rows();
      int64_t depcols = depmat.Cols();
      depmat.Redim(deprows, depcols);
      if(depcols == 3) {
        for(int i=0; i<3; i++) {
          depmat(i,i) = 1.;
        }
      } else if(depcols == 5) {
        for(int i=0; i<5; i++) {
          depmat(i+3,i) = 1.;
        }
      } else {
        std::cout << "TPZHDivSHybridApproxCreator unexpected number of columns in the depend matrix. Expected 3 or 5, got " << depcols << std::endl;
        DebugStop();
      }
      depcast->SetDepMatrix(depmat);
      depmat.Print("depmat");
      dep = dep->fNext;
    }
  }
  mfmesh->ExpandSolution();
  mfmesh->CleanUpUnconnectedNodes();
}
