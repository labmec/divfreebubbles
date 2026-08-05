#include "TPZH1HybridApproxCreator.h"
#include "pzcmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzmultiphysicselement.h"
#include "pzelementgroup.h"
#include "TPZAlgebraicInterface.h"
#include "TPZMaterialDataT.h"
#include "pzfmatrix.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("TPZH1HybridApproxCreator");
#endif

TPZH1HybridApproxCreator::TPZH1HybridApproxCreator(TPZGeoMesh * gmesh) : TPZH1ApproxCreator(gmesh),
fLowOrderFluxHybridized(false) {
}

TPZMultiphysicsCompMesh *TPZH1HybridApproxCreator::CreateHybridApproximationSpace() {
    if (HybridType() == HybridizationType::ENone) {
        return TPZH1ApproxCreator::CreateApproximationSpace();
    }

    return TPZH1ApproxCreator::CreateApproximationSpace();
}

/// Add a hybrid geometric configuration
void TPZH1HybridApproxCreator::AddHybrid(std::map<int64_t,TPZHybrid> &tree, int64_t interface_id, TPZHybrid &hybrid) {
//    std::cout << "Hybrid adding interface " << interface_id << " left " << hybrid.fLeft << " right " << hybrid.fRight << std::endl;
#ifdef PZDEBUG
    if(tree.find(interface_id) != tree.end()) {
        DebugStop();
    }
#endif
    tree[interface_id] = hybrid;
}

    /// Add geometric elements to represent squared hybridization
void TPZH1HybridApproxCreator::AddHybridSquareGeoElements() {
    // Implementation for adding squared hybrid geometric elements
    
  if(fHybridType != HybridizationType::EStandardSquared) DebugStop();
    if(0) {
        std::cout << "Hybridization structure before hybrid square\n";
        for(auto &it : fHybridizationData.fInterfaces) {
            TPZHybrid &hybridlarge = it.second;
            std::cout << "lag " << it.first << " left " << hybridlarge.fLeft << " right " << hybridlarge.fRight << std::endl;
        }
    }
#ifdef PZDEBUG
    std::set<int64_t> largeflux;
#endif


    std::map<int64_t,TPZHybrid> hybridsq;
    int bchybridlevel = fHybridizationData.fHybridizeBCLevel;
    for(auto &it : fHybridizationData.fInterfaces) {
        TPZGeoEl *gelinterface = fGeoMesh->Element(it.first);
        TPZGeoEl *gelwrap = fGeoMesh->Element(it.second.fLeft);
        if(gelwrap->MaterialId() != fHybridizationData.fWrapMatId) DebugStop();
        TPZGeoElSide wrapside(gelwrap);
        
        bool hasLargeElementNeighbour = wrapside.HasLowerLevelNeighbour(fHybridizationData.fWrapMatId);
        TPZGeoElSide wrapneighbour = wrapside.Neighbour().HasNeighbour(fHybridizationData.fWrapMatId);
        TPZGeoElSide lagrangeneighbour = wrapside.HasNeighbour(fHybridizationData.fLagrangeMatId);
        bool haswrapNeighbour = (wrapneighbour != wrapside);
        bool hasbcNeighbour = wrapside.HasNeighbour(GetBCMatIds());
        bool isLargeElement = false;
        if(!haswrapNeighbour && !hasLargeElementNeighbour && !hasbcNeighbour) {
            // we are in the case where there is no wrap neighbour, no large element neighbour and no bc neighbour, so we are the large element and we need to create the second interface and second lagrange elements
            isLargeElement = true;
        }
        if(isLargeElement || hasLargeElementNeighbour) {
          // We dont implement semi hybridization
          std::cout << "Large element neighbour " << hasLargeElementNeighbour << " is large element " << isLargeElement << std::endl;
          std::cout << "Semi hybridization is not implemented\n";
          std::cout << "Element " << gelinterface->Index() << " interface material id " << gelinterface->MaterialId() << std::endl;
          DebugStop();
        }
        // find the element configuration :
        // - it is an equal level element
        // - it is a small element linked to a large element
        // - it is a large element
        if(hasbcNeighbour) {
          // we do not doubly hybridize boundary conditions
          // copy the TPZHybrid object
          TPZHybrid bchybrid(it.second);
          AddHybrid(hybridsq,it.first,bchybrid);
          continue;
        }
        if(haswrapNeighbour) {
          TPZManVector<int64_t> geoelem(5,-1);
          geoelem[0] = gelwrap->Index();
          geoelem[1] = gelinterface->Index();
          int interfacematid = gelinterface->MaterialId();
          int secondinterfacematid = 0;
          // establish the second interface material id based on the first interface material id
          if(interfacematid == fHybridizationData.fLeftInterfaceMatId) {
            secondinterfacematid = fHybridizationData.fSecondLeftInterfaceMatId;
          } else if(interfacematid == fHybridizationData.fRightInterfaceMatId) {
            secondinterfacematid = fHybridizationData.fSecondRightInterfaceMatId;
          } else {
            DebugStop();
          }
          if(interfacematid == fHybridizationData.fLeftInterfaceMatId) {
            // we are on the left interface, we need to create a lagrange element and use it as the lagrange neighbour
            int nsides = gelinterface->NSides();
            auto newflux = gelinterface->CreateBCGeoEl(nsides-1,fHybridizationData.fLagrangeMatId);
            geoelem[2] = newflux->Index();
            lagrangeneighbour = TPZGeoElSide(newflux);      
          } else {
            // if we are on the right, we use the existing lagrange element as the lagrange neighbour
              int64_t lagrangeindex = it.second.fRight;
              auto gel = fGeoMesh->Element(lagrangeindex);
              if(gel->MaterialId() != fHybridizationData.fLagrangeMatId) DebugStop();
              lagrangeneighbour = TPZGeoElSide(gel);
              if(!lagrangeneighbour) DebugStop();
              geoelem[2] = lagrangeneighbour.Element()->Index();
          }
          TPZHybrid hybr1(geoelem[1],geoelem[0],geoelem[2]);
          /// we add the first hybrid element with the same interface id as the original one
          AddHybrid(hybridsq,geoelem[1],hybr1);
          continue;
        }
        if(hasLargeElementNeighbour) {
          /// this case is not implemented. Call DebugStop to avoid silent errors
          std::cout << "Large element neighbour " << hasLargeElementNeighbour << std::endl;
          std::cout << "Semi hybridization is not implemented\n";
          std::cout << "Element " << gelinterface->Index() << " interface material id " << gelinterface->MaterialId() << std::endl;
          DebugStop();
        }
    }
#ifdef PZDEBUG
    {
        // number of times an element appears in a position (only once)
        std::set<int64_t> asleft, asright, aslagrange;
        // number of times an element is referenced (max 2)
        std::map<int64_t,int> numref;
        int error = 0;
        for(auto &it : hybridsq) {
            int64_t lag = it.first;
            int64_t left = it.second.fLeft;
            int64_t right = it.second.fRight;
            bool leftlarge = largeflux.find(left) != largeflux.end();
            bool rightlarge = largeflux.find(right) != largeflux.end();
            if(!leftlarge && asleft.find(left) != asleft.end()) {
                std::cout << "left " << left << " appears more than once\n";
                error++;
            }
            if(rightlarge && asright.find(right) != asright.end()) {
                TPZGeoEl *gel = fGeoMesh->Element(right);
                int matid = gel->MaterialId();
                if(matid != fHybridizationData.fSecondLagrangeMatId) {
                    std::cout << "right " << right << " appears more than once\n";
                    error++;
                }
            }
            if(aslagrange.find(lag) != aslagrange.end()) {
                std::cout << "lagrange " << lag << " appears more than once\n";
                error++;
            }
            asleft.insert(left);
            asright.insert(right);
            aslagrange.insert(lag);
            numref[left]++;
            numref[right]++;
            numref[lag]++;
        }
        for(auto &it : numref) {
            bool islarge = largeflux.find(it.first) != largeflux.end();
            if(!islarge && it.second > 2) {
                std::cout << "element " << it.first << " appears more than twice\n";
                error++;
            }
        }
        if(error > 0) {
            for(auto &it : hybridsq) {
                TPZHybrid &hybridlarge = it.second;
                std::cout << "lag " << it.first << " left " << hybridlarge.fLeft << " right " << hybridlarge.fRight << std::endl;
            }
            DebugStop();
        }
    }
#endif
#ifdef PZ_LOG
            if(logger.isDebugEnabled())
            {
                std::stringstream out;
                for(auto &it : hybridsq) {
                    TPZHybrid &hybridlarge = it.second;
                    out << "lag " << it.first << " left " << hybridlarge.fLeft << " right " << hybridlarge.fRight << std::endl;
                }
                LOGPZ_DEBUG(logger, out.str())
            }
#endif
    fHybridizationData.fInterfaces = hybridsq;
}

    /// Put elements in element groups
void TPZH1HybridApproxCreator::GroupElements(TPZMultiphysicsCompMesh *mcmesh) {
  if(HybridType() != HybridizationType::EStandardSquared) {
      DebugStop();
  }
  // verify that the lagrange elements do not share connects.
  {
    std::map<int64_t,int> connectcount;
    int64_t nel = mcmesh->NElements();
    for(int64_t el = 0; el < nel; el++) {
      TPZCompEl *cel = mcmesh->Element(el);
      if(!cel) continue;
      TPZGeoEl *gel = cel->Reference();
      if(!gel) continue;
      if(gel->MaterialId() == fHybridizationData.fLagrangeMatId) {
        int nc = cel->NConnects();
        for(int i = 0; i < nc; i++) {
          int64_t connectindex = cel->ConnectIndex(i);
          connectcount[connectindex]++;
        }
      }
    }
    for(const auto &pair : connectcount) {
      if(pair.second > 1) {
        std::cout << "Error: connect " << pair.first << " is shared by " << pair.second << " lagrange elements\n";
        DebugStop();
      }
    }
  }
  fHybridType = HybridizationType::EStandard;
  TPZH1ApproxCreator::GroupElements(mcmesh);
  fHybridType = HybridizationType::EStandardSquared;
  // add the algebraic elements to the groups
  TPZVec<int64_t> connectTOgroup(mcmesh->NConnects(), -1);
  int64_t nel = mcmesh->NElements();
  for(int64_t el = 0; el < nel; el++) {
      TPZCompEl *cel = mcmesh->Element(el);
      if(!cel) continue;
      TPZElementGroup *group = dynamic_cast<TPZElementGroup *>(cel);
      if(group) {
        std::set<int64_t> groupconnects;
        group->BuildConnectList(groupconnects);
        for(auto i : groupconnects) {
              int64_t connectindex = i;
              if(connectTOgroup[connectindex] >= 0) {
                  connectTOgroup[connectindex] = -2; // mark as invalid.
              } else {
                  connectTOgroup[connectindex] = group->Index();
              }
          }
      }
  }
  if(fLowOrderFluxHybridized == false) {
    return;
  }
  // put the lagrange elements in the group to which they belong
  std::set<int> lagrangeids = this->GetBCMatIds();
  lagrangeids.insert(fHybridizationData.fLagrangeMatId);
  for(int64_t el = 0; el<nel; el++) {
    TPZCompEl *cel = mcmesh->Element(el);
    if(!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    if(!gel) continue;
    if(lagrangeids.find(gel->MaterialId()) != lagrangeids.end()) {
      int nc = cel->NConnects();
      int64_t groupindex = -1;
      for(int i = 0; i < nc; i++) {
        int64_t connectindex = cel->ConnectIndex(i);
        int64_t groupid = connectTOgroup[connectindex];
        if(groupid < 0) {
          std::cout << "Error: lagrange element " << cel->Index() << " has connect " << connectindex << " which does not belong to any group\n";
          DebugStop();
        }
        if(groupindex < 0) {
          groupindex = groupid;
        } else if(groupindex != groupid) {
          std::cout << "Error: lagrange element " << cel->Index() << " has connects belonging to different groups " << groupindex << " and " << groupid << std::endl;
          DebugStop();
        }
      }
      if(groupindex < 0) {
        std::cout << "Error: lagrange element " << cel->Index() << " has no connects belonging to any group\n";
        DebugStop();
      }
      TPZElementGroup *group = dynamic_cast<TPZElementGroup *>(mcmesh->Element(groupindex));
      if(!group) DebugStop();
      group->AddElement(cel);
    }
  }
  if(0)
  {
    mcmesh->ComputeNodElCon();
    std::ofstream out("cmesh_before_algebraic_group.txt");
    mcmesh->Print(out);
    out << "connect to group\n";
    for(size_t i = 0; i < connectTOgroup.size(); i++) {
      out << "connect " << i << " group " << connectTOgroup[i] << std::endl;
    }
  }
  for(int64_t el = 0; el < nel; el++) {
    TPZCompEl *cel = mcmesh->Element(el);
    if(!cel) continue;
    TPZAlgebraicInterface *alg = dynamic_cast<TPZAlgebraicInterface *>(cel);
    if(alg) {
      int64_t algGroupIndex = -1;
      int nc = alg->NConnects();
      for(int i = 0; i < nc; i++) {
        int64_t connectindex = alg->ConnectIndex(i);
        int64_t groupindex = connectTOgroup[connectindex];
        if(groupindex < 0) continue;
        if(algGroupIndex < 0) {
          algGroupIndex = groupindex;
        } else if(algGroupIndex != groupindex) {
          std::cout << "Error: algebraic element " << cel->Index() << " has connects in different groups " << algGroupIndex << " and " << groupindex << std::endl;
          DebugStop();
        }
      }
      if(algGroupIndex < 0) {
        std::cout << "Error: algebraic element " << cel->Index() << " has no connects in any group\n";
        for(int i = 0; i < nc; i++) {
          int64_t connectindex = alg->ConnectIndex(i);
          std::cout << "connect " << connectindex << " group " << connectTOgroup[connectindex] << std::endl;
        }
        DebugStop();
      }
      TPZElementGroup *group = dynamic_cast<TPZElementGroup *>(mcmesh->Element(algGroupIndex));
      if(!group) DebugStop();
      group->AddElement(cel);
    }
  }
  mcmesh->ComputeNodElCon();
  if(0)
   {
    std::ofstream out("cmesh_after_group.txt");
    mcmesh->Print(out);
  }
}
    
/// Create condensed elements around group elements
/// this method will adjust the connect count of DOF's that should not be condensed
void TPZH1HybridApproxCreator::CondenseElements(TPZMultiphysicsCompMesh *mcmesh) {
  if(HybridType() != HybridizationType::EStandardSquared) {
      DebugStop();
  }
  // if the lower order fluxes are not hybridized, the rigid body modes cannot be condensed.
  if(fLowOrderFluxHybridized == false) {
    fHybridType = HybridizationType::EStandard;
  }
  TPZH1ApproxCreator::CondenseElements(mcmesh);
  fHybridType = HybridizationType::EStandardSquared;

}
#include "pzcompel.h"
/// @brief hybridize the low order fluxes
void TPZH1HybridApproxCreator::HybridizeLowOrderFluxes(TPZMultiphysicsCompMesh &mfmesh, CtoMFCel &geltogel)
{
  for(auto iter : geltogel) {
    int locindex1(-1),locindex2(-1);
    int64_t gel1 = iter.first;
    int64_t gel2 = iter.second;
    if(gel2 < 0) DebugStop();
    TPZGeoMesh *gmesh = mfmesh.Reference();
    TPZCompEl *cel1 = gmesh->Element(gel1)->Reference();
    TPZCompEl *cel2 = gmesh->Element(gel2)->Reference();
    if(!cel1 || !cel2) DebugStop();
    if(cel1->NConnects() != 1 || cel2->NConnects() != 1) DebugStop();
    locindex1 = cel1->ConnectIndex(0);
    locindex2 = cel2->ConnectIndex(0);
    if(locindex1 < 0 || locindex2 < 0) DebugStop();
    TPZConnect &c1 = cel1->Connect(0);
    TPZConnect &c2 = cel2->Connect(0);
    if(!c1.HasDependency() || !c2.HasDependency()) DebugStop();
    TPZConnect::TPZDependBase *dep1 = c1.FirstDepend();
    TPZConnect::TPZDependBase *dep2 = c2.FirstDepend();
    if(!dep1 || !dep2) DebugStop();
    int64_t depindex1 = dep1->fDepConnectIndex;
    int64_t depindex2 = dep2->fDepConnectIndex;
    if(depindex1 != depindex2) DebugStop();
    TPZConnect &c = mfmesh.ConnectVec()[depindex1];
    int64_t newindex = mfmesh.AllocateNewConnect(c);
    dep2->fDepConnectIndex = newindex;
    depindex2 = newindex;
    int64_t disp_connectIndex = mfmesh.AllocateNewConnect(c);
    TPZConnect &cdisp = mfmesh.ConnectVec()[disp_connectIndex];
    // create two lagrange multiplier elements to link the hybridized flux to the two low order fluxes
    TPZAlgebraicInterface *lag1 = new TPZAlgebraicInterface(mfmesh, NULL);
    lag1->SetConnectIndex(0, depindex1);
    lag1->SetConnectIndex(1, disp_connectIndex);
    TPZAlgebraicInterface *lag2 = new TPZAlgebraicInterface(mfmesh, NULL);
    lag2->SetConnectIndex(0, depindex2);
    lag2->SetConnectIndex(1, disp_connectIndex);
    lag2->SetMultiplier(-1.);
  }
  mfmesh.ExpandSolution();
  mfmesh.ComputeNodElCon();
  mfmesh.CleanUpUnconnectedNodes();
  fLowOrderFluxHybridized = true;
}
