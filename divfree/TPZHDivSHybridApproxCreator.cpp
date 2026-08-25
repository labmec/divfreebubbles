#include "TPZHDivSHybridApproxCreator.h"
#include "TPZAlgebraicInterface.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterial.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"
#include "pzgeoelbc.h"
#include "pzintel.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("TPZHDivSHybridApproxCreator");
#endif

TPZHDivSHybridApproxCreator::TPZHDivSHybridApproxCreator(TPZGeoMesh *gmesh)
    : TPZHDivApproxCreator(gmesh) {
}

void TPZHDivSHybridApproxCreator::AddHybridizationGeoElements() {
  /// For fHybridType == HybridizationType::EStandard, wrap, interface and lagrange geometric elements are created;
  ///     if fHybridizeBCLevel == 0, the boundary is not hybridized, if fHybridizeBCLevel == 1, it is.
  /// For fHybridType == HybridizationType::EStandardSquared, wrap, interface, lagrange, second interface and second lagrange
  /// geometric elements are created;
  ///     BC can be hybridized 0, 1 or 2 times in this set up;
  /// fHybridType == HybridizationType::Semi has yet to be implemented.

  if (fHybridType == HybridizationType::ENone) {
    std::cout << "Hybridization type is none, this method should not be called\n";
    DebugStop();
  }
  if (fHybridizationData.fWrapMatId == -123456) {
    std::cout << "\nERROR! Please call TPZApproxCreator::ComputePeriferalMaterialIds() before TPZApproxCreator::AddHybridizationGeoElements()" << std::endl;
    DebugStop();
  }

#ifdef PZ_LOG
  std::map<int, int> numcreated;
#endif
  // create the wrap and interface geometric elements
  // the creation of wrap and interface needs to happen in separate loops when there are hanging nodes
  int64_t nel = fGeoMesh->NElements();
  int dim = fGeoMesh->Dimension();
  std::set<int> bcMatIds = GetBCMatIds();
  // Wrap and interface creation for standard hyb.
  for (int64_t el = 0; el < nel; el++) {
    TPZGeoEl *gel = fGeoMesh->Element(el);
    if (!gel || gel->HasSubElement() || gel->Dimension() != dim) continue;
    int nsides = gel->NSides();
    int side = gel->FirstSide(dim - 1);
    // loop over the sides of dimension dim-1
    for (; side < nsides - 1; side++) {
      TPZGeoElSide gelside(gel, side);
      // we want to create side elements of type
      // first fMatWrapId
      TPZGeoElSide neighbour = gelside.Neighbour();
#ifdef PZDEBUG
      {
        int neighMatId = neighbour.Element()->MaterialId();
        if (neighMatId == fHybridizationData.fWrapMatId) {
          std::cout << __PRETTY_FUNCTION__ << " should be called only once!\n";
          DebugStop();
        }
      }
#endif
      // if there is a neighbour that is a boundary condition
      // do not create the wrap layers
      bool hasBCNeighbour = gelside.HasNeighbour(bcMatIds);
      if (hasBCNeighbour && fHybridizationData.fHybridizeBCLevel == 0) {
        // no interface will be created between the element and a flux space
        continue;
      }
      if (gel->NormalOrientation(side) == -1) {
        continue; // we only create the wrap and interface elements for one of the two neighbouring elements
      }
      TPZGeoElBC(gelside, fHybridizationData.fWrapMatId);
#ifdef PZ_LOG
      numcreated[fHybridizationData.fWrapMatId]++;
#endif
#ifdef PZDEBUG
      neighbour = gelside.Neighbour();
      if (neighbour.Element()->MaterialId() != fHybridizationData.fWrapMatId) {
        DebugStop();
      }
#endif
    }
  }

  if (fHybridType == HybridizationType::EStandardSquared) {
    DebugStop(); // this is not implemented yet
  }
}

TPZCompMesh *TPZHDivSHybridApproxCreator::CreateHDivSpace() {
  auto hybrid = fHybridType;
  fHybridType = HybridizationType::EStandard; // we do not want to create the hybridization geometric elements in the HDiv space
  // @TODO this will not create the wrap hdiv bound elements
  auto cmesh = TPZHDivApproxCreator::CreateHDivSpace();
  fHybridType = hybrid; // restore the original hybridization type
  return cmesh;
}

/// Insert interface periferal material objects related to geometric objects created during hybridization
void TPZHDivSHybridApproxCreator::InsertInterfaceMaterialObjects(TPZMultiphysicsCompMesh *mphys) {}

/// Create interface elements on hybridizes spaces
/// @param mphys multiphysics compmesh
void TPZHDivSHybridApproxCreator::AddInterfaceComputationalElements(TPZMultiphysicsCompMesh *mfmesh) {
  // return;
  mfmesh->LoadReferences();
  int64_t nel = mfmesh->NElements();
  std::set<int> bcMatIds = GetBCMatIds();
  bcMatIds.insert(fHybridizationData.fWrapMatId);
  for (int64_t el = 0; el < nel; el++) {
    TPZCompEl *cel = mfmesh->Element(el);
    if (!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    if (!gel) continue;
    int matid = gel->MaterialId();
    if (bcMatIds.find(matid) == bcMatIds.end()) {
      continue;
    }
    TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if (!mcel) {
      std::cout << "TPZHDivSHybridApproxCreator::AddInterfaceComputationalElements error. The element is not a multiphysics element\n";
      DebugStop();
    }
    TPZInterpolatedElement *wrapel = dynamic_cast<TPZInterpolatedElement *>(mcel->Element(0));
    if (!wrapel) {
      std::cout << "TPZHDivSHybridApproxCreator::AddInterfaceComputationalElements error. The first element of the multiphysics element is not an interpolated element\n";
      DebugStop();
    }
    if (mcel->NConnects() != 1) {
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
    // restrain the connect of the neighbour element to the same connects
    // any neigbour element that have a non constrained connect will be constrained
    TPZConnect &c = mcel->Connect(0);
    {
      TPZGeoElSide gelside(gel);
      TPZGeoElSide neighbour = gelside.Neighbour();
      while(neighbour != gelside) {
        TPZGeoEl *neiggel = neighbour.Element();
        TPZCompEl *neigcel = neiggel->Reference();
        if (neigcel) {
          TPZMultiphysicsElement *neigmcel = dynamic_cast<TPZMultiphysicsElement *>(neigcel);
          if (!neigmcel) {
            std::cout << "TPZHDivSHybridApproxCreator::AddInterfaceComputationalElements error. The neighbour element is not a multiphysics element\n";
            DebugStop();
          }
          // get the first element of the multiphysics element.
          TPZInterpolatedElement *neigdivel = dynamic_cast<TPZInterpolatedElement *>(neigmcel->Element(0));
          int neigconlocindex = neigdivel->SideConnectLocId(0,neighbour.Side());
          TPZConnect &neigc = neigmcel->Connect(neigconlocindex);
          if (neigc.NDof() != c.NDof()) {
            std::cout << "TPZHDivSHybridApproxCreator::AddInterfaceComputationalElements error. The neighbour connect has a different number of dofs than the original connect\n";
            DebugStop();
          }
          if (!neigc.HasDependency()) {
            neigdivel->SetSideOrient(neighbour.Side(), -1);
            CopyRestraintConnect(c, neigc);
          }
        }
        neighbour = neighbour.Neighbour();
      }
    }
    TPZConnect::TPZDependBase *dep = c.FirstDepend();
    while (dep) {
      TPZConnect::TPZDepend<STATE> *depcast = dynamic_cast<TPZConnect::TPZDepend<STATE> *>(dep);
      if (!depcast) {
        std::cout << "TPZHDivSHybridApproxCreator::AddInterfaceComputationalElements error. The depend is not of type TPZConnect::TPZDepend<STATE>\n";
        DebugStop();
      }
      TPZFMatrix<STATE> depmat = depcast->GetDepMatrix();
      // depmat.Print("depmat");

      /*
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
      */
      dep = dep->fNext;
    }
  }
  mfmesh->LoadReferences();
  if (fShouldHybridizeLowOrderFluxes) {
    HybridizeLowOrderFluxes(*mfmesh);
  }
  mfmesh->ExpandSolution();
  mfmesh->CleanUpUnconnectedNodes();
  if(0)
  {
    std::ofstream out("hdiv mesh_after_interface.txt");
    mfmesh->MeshVector()[0]->Print(out);
  }
}

/// Groups the elements in data structure to be condensed
/// @param mcmesh multiphysics compmesh with elements to be condensed
void TPZHDivSHybridApproxCreator::GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh) {
  TPZGeoMesh *gmesh = mcmesh->Reference();
  gmesh->ResetReference();
  mcmesh->LoadReferences();
  int64_t nel = mcmesh->NElements();
  int dim = mcmesh->Dimension();
  int64_t ncon = mcmesh->NConnects();
  TPZManVector<int64_t,300> connectgroup(ncon,-1);
  std::set<int> bcMatIds = GetBCMatIds();
  std::list<TPZCompEl *> elgrouplist;
  for (int64_t el = 0; el < nel; el++) {
    TPZCompEl *cel = mcmesh->Element(el);
    if (!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    if (!gel) continue;
    if (gel->Dimension() != dim) continue;
    TPZElementGroup *elgroup = new TPZElementGroup(*mcmesh);
    elgroup->AddElement(cel);
    int64_t groupindex = elgroup->Index();
    std::set<int64_t> connectlist;
    cel->BuildConnectList(connectlist);
    int ncon = cel->NConnects();
    for (auto conindex : connectlist) {
      if (connectgroup[conindex] != -1) {
        connectgroup[conindex] = -2; // this connect is associated with more than one element group
      } else if(connectgroup[conindex] == -2) {
        DebugStop(); // this should not happen, a connect should not be associated with more than two element groups
      } else {
        connectgroup[conindex] = groupindex;
      }
    }
    int firstside = gel->FirstSide(dim - 1);
    int lastside = gel->NSides();
    for (int side = firstside; side < lastside - 1; side++) {
      TPZGeoElSide gelside(gel, side);

      TPZGeoElSide neighbour = gelside.Neighbour();
      int neighmatid = neighbour.Element()->MaterialId();
      if(neighmatid == fHybridizationData.fWrapMatId) {
        TPZGeoEl *gelwrap = neighbour.Element();
        TPZCompEl *celwrap = gelwrap->Reference();
        if (!celwrap) {
          std::cout << "TPZHDivSHybridApproxCreator::GroupAndCondenseElements error. The wrap element does not have a computational element\n";
          DebugStop();
        }
        int64_t wrapconindex = celwrap->ConnectIndex(0);
        if (connectgroup[wrapconindex] != groupindex) {
          std::cout << "TPZHDivSHybridApproxCreator::GroupAndCondenseElements error. The wrap connect is not associated with the same element group as the original element\n";
          DebugStop();
        }
        elgroup->AddElement(celwrap);
      }
      auto bcneighbor = gelside.HasNeighbour(bcMatIds);
      if (bcneighbor) {
        TPZGeoEl *gel = bcneighbor.Element();
        TPZCompEl *cel = gel->Reference();
        if (!cel) {
          std::cout << "TPZHDivSHybridApproxCreator::GroupAndCondenseElements error. The neighbour element does not have a computational element\n";
          DebugStop();
        }
        int64_t bcconindex = cel->ConnectIndex(0);
        if (connectgroup[bcconindex] != groupindex) {
          std::cout << "TPZHDivSHybridApproxCreator::GroupAndCondenseElements error. The neighbour connect is not associated with the same element group as the original element\n";
          DebugStop();
        }
        elgroup->AddElement(cel);
      }
    }
    elgrouplist.push_back(elgroup);
  }
  // add all elements that have a connect that is associated with an element group to that element group
  for (int64_t el = 0; el < nel; el++) {
    TPZCompEl *cel = mcmesh->Element(el);
    if (!cel) continue;
    int ncon = cel->NConnects();
    std::set<int64_t> connectlist;
    cel->BuildConnectList(connectlist);
    int64_t elgroupindex = -1;
    for (auto conindex : connectlist) {
      int64_t groupindex = connectgroup[conindex];
      if (groupindex == -1 || groupindex == -2) continue;
      if(elgroupindex == -1) {
        elgroupindex = groupindex;
      } else if(groupindex != elgroupindex) {
        std::cout << "TPZHDivSHybridApproxCreator::GroupAndCondenseElements error. The connect is associated with more than one element group\n";
        DebugStop();
      }
      TPZElementGroup *elgroup = dynamic_cast<TPZElementGroup *>(mcmesh->Element(groupindex));
      if (!elgroup) {
        std::cout << "TPZHDivSHybridApproxCreator::GroupAndCondenseElements error. The element group is not a TPZElementGroup\n";
        DebugStop();
      }
      elgroup->AddElement(cel);
    }
  }
  mcmesh->ComputeNodElCon();
  auto &meshvec = mcmesh->MeshVector();
  int nmeshes = meshvec.size();
  auto &activevec = mcmesh->GetActiveApproximationSpaces();
  TPZManVector<int64_t, 7> firstconnect(meshvec.size() + 1, 0);
  for (int m = 0; m < nmeshes; m++) {
    if (activevec[m] == 0) {
      firstconnect[m + 1] = firstconnect[m];
    } else {
      firstconnect[m + 1] = firstconnect[m] + meshvec[m]->NConnects();
    }
  }
  if (this->fIsRBSpaces) {
    int rgbspace = 4;
    for (int64_t ic = firstconnect[rgbspace]; ic < firstconnect[rgbspace + 1]; ic++) {
      TPZConnect &c = mcmesh->ConnectVec()[ic];
      c.IncrementElConnected();
    }
  }

  for (auto elgroup : elgrouplist) {
    TPZCondensedCompElT<STATE> *condensed = new TPZCondensedCompElT<STATE>(elgroup, false);
  }
  mcmesh->ComputeNodElCon();
  mcmesh->CleanUpUnconnectedNodes();
}

/// @brief hybridize the low order fluxes
void TPZHDivSHybridApproxCreator::HybridizeLowOrderFluxes(TPZMultiphysicsCompMesh &mfmesh) {
  int64_t nel = mfmesh.NElements();
  for (int64_t el = 0; el < nel; el++) {
    TPZCompEl *cel = mfmesh.Element(el);
    if (!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    if (!gel) continue;
    int matid = gel->MaterialId();
    if (matid != fHybridizationData.fWrapMatId) continue;
    TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if (!mcel) DebugStop();
    TPZGeoElSide leftside, rightside;
    TPZGeoElSide gelside(gel);
    TPZGeoElSide neighbour = gelside.Neighbour();
    if (neighbour.Element()->NormalOrientation(gelside.Side()) == 1) {
      leftside = neighbour;
      rightside = neighbour.Neighbour();
    } else {
      rightside = neighbour;
      leftside = neighbour.Neighbour();
    }
    TPZGeoEl *leftgel = leftside.Element();
    TPZGeoEl *rightgel = rightside.Element();
    if (leftgel->Dimension() != mfmesh.Dimension() || rightgel->Dimension() != mfmesh.Dimension()) DebugStop();
    if (!leftgel->Reference() || !rightgel->Reference()) DebugStop();
    int locindex1(-1), locindex2(-1);
    TPZGeoMesh *gmesh = mfmesh.Reference();
    auto *mfcel1 = dynamic_cast<TPZMultiphysicsElement *>(leftgel->Reference());
    auto *mfcel2 = dynamic_cast<TPZMultiphysicsElement *>(rightgel->Reference());
    if (!mfcel1 || !mfcel2) DebugStop();
    auto *cel1 = dynamic_cast<TPZInterpolatedElement *>(mfcel1->Element(0));
    auto *cel2 = dynamic_cast<TPZInterpolatedElement *>(mfcel2->Element(0));
    if (!cel1 || !cel2) DebugStop();
    locindex1 = cel1->MidSideConnectLocId(leftside.Side());
    locindex2 = cel2->MidSideConnectLocId(rightside.Side());
    if (locindex1 < 0 || locindex2 < 0) DebugStop();
    TPZConnect &c1 = mfcel1->Connect(locindex1);
    TPZConnect &c2 = mfcel2->Connect(locindex2);
    // both connects should be the same
    if (!c1.HasDependency() || !c2.HasDependency()) DebugStop();
    TPZConnect::TPZDependBase *dep1 = c1.FirstDepend();
    TPZConnect::TPZDependBase *dep2 = c2.FirstDepend();
    if (!dep1 || !dep2) DebugStop();
    int64_t depindex1 = dep1->fDepConnectIndex;
    int64_t depindex2 = dep2->fDepConnectIndex;
    if (depindex1 != depindex2) DebugStop();
    // the connect depindex1 is the one that will be used to link the hybridized flux to the low order displacment.
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
}
