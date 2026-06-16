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

/// @brief Compute the constraints to orthogonalize the restraints
// this function will compute restraints for all boundary flux connects
void TPZH1HybridApproxCreator::ComputeOrthogonalizingRestraints(TPZMultiphysicsCompMesh &mfmesh, CtoMFCel &geltogel, TPZApproxCreator::HybridizationData hybridData)
{
  int64_t nel = mfmesh.NElements();
  int meshdim = mfmesh.Dimension();
  int matidleft = hybridData.fLagrangeMatId;
  int matidright = hybridData.fLagrangeMatId;
  for (int64_t el = 0; el < nel; el++)
  {
    TPZCompEl *celleft = mfmesh.Element(el);
    if (!celleft)
      continue;
    TPZGeoEl *gelleft = celleft->Reference();
    if (!gelleft)
      continue;
    int matid = gelleft->MaterialId();
    if(matid != matidleft) continue;
    if (gelleft->Dimension() != meshdim - 1)
      DebugStop();
    if (celleft->NConnects() != 1)
      DebugStop();
    TPZConnect &cleft = celleft->Connect(0);
    int64_t cleftindex = celleft->ConnectIndex(0);
    if (cleft.HasDependency()) {
      continue;
    }
    TPZGeoEl *gelright = 0;
    {
      TPZGeoElSide gelside(gelleft);
      TPZGeoElSide neighbour = gelside.Neighbour().HasNeighbour(matidright);
      if (!neighbour || neighbour == gelside) {
        DebugStop();
      }
      gelright = neighbour.Element();
    }
    TPZCompEl *rightcel = gelright->Reference();
    if (!rightcel) {
      DebugStop();
    }
    if(rightcel->NConnects() != 1) {
      DebugStop();
    }
    TPZConnect &cright = rightcel->Connect(0);
    int64_t crightindex = rightcel->ConnectIndex(0);
    if (cright.HasDependency()) {
      DebugStop();
    }
    geltogel[gelleft->Index()] = gelright->Index();
  

    // restrain the connect of the left lagrange multiplier to a connect representing the low order flux space and the remainder of dofs of the lagrange multiplier connect
    int64_t newind1 = mfmesh.AllocateNewConnect(cleft);
    TPZConnect &cnew1 = mfmesh.ConnectVec()[newind1];
  
    int64_t newind2 = mfmesh.AllocateNewConnect(cleft);
    TPZConnect &cnew2 = mfmesh.ConnectVec()[newind2];

    std::cout << "cleftindex " << cleftindex <<
    " crightindex " << crightindex << " newind1 " << newind1 << " newind2 " << newind2 << std::endl;
    cnew1.SetNState(1);
    cnew1.SetNShape(3);
    int64_t seq1 = cnew1.SequenceNumber();
    mfmesh.Block().Set(seq1, 3);
    cnew2.SetNState(1);
    cnew2.SetNShape(cleft.NDof() - 3);
    int64_t seq2 = cnew2.SequenceNumber();
    mfmesh.Block().Set(seq2, cleft.NDof() - 3);
    mfmesh.Block().Resequence();
    TPZMultiphysicsElement *mfel = dynamic_cast<TPZMultiphysicsElement *>(celleft);
    if (!mfel)
      DebugStop();
    // we assume the fluxmesh is the first mesh in the multiphysics computational mesh
    TPZCompEl *fluxel = mfel->Element(0);
    TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(fluxel);
    RestraintConnect(intel, mfel, newind1, newind2);

    // restrain the connect of the right lagrange multiplier to the same low order flux connect and the remainder of dofs of the lagrange multiplier connect
    TPZMultiphysicsElement *mfelright = dynamic_cast<TPZMultiphysicsElement *>(rightcel);
    if (!mfelright)
      DebugStop();
    if (!cleft.HasDependency()) {
      DebugStop();
    }
    TPZConnect::TPZDependBase *dep = cleft.FirstDepend();
    TPZConnect::TPZDependBase *depnext = dep->fNext;
    if(!depnext || depnext->fNext) {
      std::cout << "We expect only one dependency for the left lagrange connect\n";
      DebugStop();
    }
    TPZManVector<TPZConnect::TPZDepend<STATE> *,2> depvec(2);
    depvec[0] = dynamic_cast<TPZConnect::TPZDepend<STATE> *>(depnext);
    depvec[1] = dynamic_cast<TPZConnect::TPZDepend<STATE> *>(dep);
    for(int i = 0; i < 2; i++) {
      TPZConnect::TPZDepend<STATE> *dep = depvec[i];
      if(!dep) DebugStop();
      int64_t nr = dep->fDepMatrix.Rows();
      int64_t nc = dep->fDepMatrix.Cols();
      cright.AddDependency(crightindex, dep->fDepConnectIndex, dep->fDepMatrix, 0,0,nr,nc);
    }
  }
  mfmesh.ExpandSolution();
  mfmesh.ComputeNodElCon();
  mfmesh.CleanUpUnconnectedNodes();
}

/// @brief compute the projection directions for a geometric element
// each direction is a line of the matrix
void TPZH1HybridApproxCreator::ProjectionDirections(TPZGeoEl *gel, TPZFMatrix<REAL> &projdir)
{
  int geldim = gel->Dimension();
  TPZGeoElSide gelside(gel);
  projdir.Redim(geldim + 1, geldim + 1);
  TPZFNMatrix<9, REAL> axes(geldim, 3), jac(geldim, geldim), jacinv(geldim, geldim), gradx(3, geldim);
  TPZManVector<REAL, 3> center(geldim, 0.);
  gelside.CenterPoint(center);
  gel->GradX(center, gradx);
  REAL detjac;
  gel->Jacobian(gradx, jac, axes, detjac, jacinv);
  if (geldim == 1)
  {
    // first line is in the direction of the normal
    projdir(0, 0) = axes(0, 1);
    projdir(0, 1) = -axes(0, 0);
    // second line is in the direction of the axes
    projdir(1, 0) = axes(0, 0);
    projdir(1, 1) = axes(0, 1);
  }
  else if (geldim == 2)
  {
    // the first direction is in the direction of the vector product of the axes
    DebugStop();
    // the second and third directions are the directions of the axes
    for (int i = 0; i < 3; i++)
    {
      projdir(1, i) = axes(0, i);
      projdir(2, i) = axes(1, i);
    }
  }
  // std::cout << "projection directions\n";
  // projdir.Print("projdir = ", std::cout, EMathematicaInput);
}

/// @brief Values of the functions onto which the boundary space will be _Projected
/// out : funcval (nxdim) vector values for each component (force, momentum, others)
void TPZH1HybridApproxCreator::ProjectionValues(TPZVec<REAL> &xrelative, TPZFMatrix<REAL> &projdir, int ncorner, TPZFMatrix<REAL> &funcval)
{
  if (ncorner == 2)
  {
    int dim = 2;
    // project xrelative onto the tangential direction
    REAL relativeT = xrelative[0] * projdir(1, 0) + xrelative[1] * projdir(1, 1);
    for (int i = 0; i < 2; i++)
    {
      funcval(i, 0) = projdir(0, i);
      funcval(i, 1) = projdir(1, i);
    }
    // moment value
    funcval(0, 2) = -relativeT * projdir(1, 1);
    funcval(1, 2) = relativeT * projdir(1, 0);
    // function from the center out
    funcval(0, 3) = relativeT * projdir(1, 0);
    funcval(1, 3) = relativeT * projdir(1, 1);
  }
  else
  {
    std::cout << "Please implement me\n";
    DebugStop();
  }
  // funcval.Print("functionVal = ", std::cout, EMathematicaInput);
}

#include "pzvec_extras.h"
/// @brief compute the projection matrix
void TPZH1HybridApproxCreator::ComputeProjectionMatrix(TPZInterpolationSpace *intel, TPZFMatrix<REAL> &projection)
{
  TPZGeoEl *gel = intel->Reference();
  TPZConnect &c = intel->Connect(0);
  int intorder = 2 * c.Order();
  TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gel->NSides() - 1, intorder);
  projection.Zero();
  TPZGeoElSide gelside(gel);
  TPZManVector<REAL, 3> xcenter(3, 0.), xval(3, 0.);
  gelside.CenterX(xcenter);
  // std::cout << "x center " << xcenter << std::endl;
  int geldim = gel->Dimension();
  int meshdim = geldim + 1;
  int ncorner = gel->NCornerNodes();
  // compute the projection directons (2x2) or (3x3)
  TPZFNMatrix<9, REAL> projdir(geldim + 1, geldim + 1, 0.);
  ProjectionDirections(gel, projdir);
  TPZMaterialDataT<STATE> matdata;
  intel->InitMaterialData(matdata);
  int nshape = ncorner;
  // compute the L2 projection of the shape functions on forces and moments
  TPZFNMatrix<144, REAL> L2(ncorner * meshdim, ncorner * meshdim, 0.), funcval(meshdim, ncorner * meshdim), rhs(nshape * meshdim, nshape * meshdim, 0.);
  int np = intrule->NPoints();
  REAL weight;
  TPZManVector<REAL, 3> point(geldim);
  TPZFNMatrix<9, REAL> jac(geldim, geldim), jacinv(geldim, geldim), axes(geldim, 3), gradx(3, geldim);
  REAL detjac;
  for (int ip = 0; ip < np; ip++)
  {
    intrule->Point(ip, point, weight);
    intel->ComputeRequiredData(matdata, point);
    gel->X(point, xval);
    gel->GradX(point, gradx);
    gel->Jacobian(gradx, jac, axes, detjac, jacinv);
    xval = xval - xcenter;
    ProjectionValues(xval, projdir, ncorner, funcval);
    for (int ish = 0; ish < nshape; ish++)
    {
      for (int jsh = 0; jsh < nshape; jsh++)
      {
        for (int d = 0; d < meshdim; d++)
        {
          L2(ish * meshdim + d, jsh * meshdim + d) += matdata.phi(ish, 0) * matdata.phi(jsh, 0) * detjac * weight;
        }
      }
      for (int d = 0; d < meshdim; d++)
      {
        for (int f = 0; f < nshape * meshdim; f++)
        {
          rhs(ish * meshdim + d, f) += matdata.phi(ish, 0) * funcval(d, f) * detjac * weight;
        }
      }
    }
  }
  L2.SolveDirect(rhs, ECholesky);
  projection = rhs;
  delete intrule;
}

/// @brief compute the restraints of a connect for a geometric element
void TPZH1HybridApproxCreator::RestraintConnect(TPZInterpolationSpace *intel, TPZMultiphysicsElement *mfcel, int64_t newind1, int64_t newind2)
{
  TPZConnect &c = mfcel->Connect(0);
  TPZCompMesh *mfmesh = mfcel->Mesh();
  TPZConnect &cnew1 = mfmesh->ConnectVec()[newind1];
  TPZConnect &cnew2 = mfmesh->ConnectVec()[newind2];
  TPZFNMatrix<36, STATE> projection(4, 4);
  ComputeProjectionMatrix(intel, projection);
  // projection.Print("proj = ", std::cout, EMathematicaInput);
  int ndof = c.NDof();
  TPZFNMatrix<64, STATE> rest(ndof, ndof, 0.);
  rest.Identity();

  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      rest(i, j) = projection(i, j);
    }
  }

  TPZFNMatrix<54, REAL> B(3, ndof), BRes(3, ndof), BRes2(3, ndof);
  ComputeBMatrix(intel, B);
  B.Multiply(rest, BRes);
  // B.Print("B =", std::cout, EMathematicaInput);
  // BRes.Print("BRes = ", std::cout, EMathematicaInput);
  TPZFNMatrix<36, REAL> rest2, rest3;
  ComputeOrthogonalizingR(BRes, rest2);
  BRes.Multiply(rest2, BRes2);
  // BRes2.Print("BRes2 = ", std::cout, EMathematicaInput);
  rest.Multiply(rest2, rest3);
  c.AddDependency(mfcel->ConnectIndex(0), newind2, rest3, 0, 3, ndof, ndof - 3);
  c.AddDependency(mfcel->ConnectIndex(0), newind1, rest3, 0, 0, ndof, 3);
}

/// @brief compute an equivalent B matrix
void TPZH1HybridApproxCreator::ComputeBMatrix(TPZInterpolationSpace *intel, TPZFMatrix<STATE> &B)
{
  if (intel->NConnects() != 1)
    DebugStop();
  TPZConnect &c = intel->Connect(0);
  int order = c.Order();
  int nstate = c.NState();
  int nshape = c.NShape();
  int ndof = c.NDof();
  B.Redim(3, ndof);
  TPZGeoEl *gel = intel->Reference();
  int geldim = gel->Dimension();
  int dim = geldim + 1;
  if (geldim != 1)
  {
    std::cout << "Please implement me\n";
    DebugStop();
  }
  TPZGeoElSide gelside(gel);
  TPZManVector<REAL, 3> xcenter(3);
  gelside.CenterX(xcenter);
  TPZMaterialDataT<STATE> matdata;
  intel->InitMaterialData(matdata);
  TPZIntPoints *intrule = gelside.CreateIntegrationRule(2 * order);
  int np = intrule->NPoints();
  REAL weight, detjac;
  TPZFNMatrix<4, REAL> jac(geldim, geldim), jacinv(geldim, geldim), gradx(3, geldim), axes(geldim, 3);
  TPZManVector<REAL, 3> pos(geldim);
  for (int ip = 0; ip < np; ip++)
  {
    intrule->Point(ip, pos, weight);
    gel->GradX(pos, gradx);
    gel->Jacobian(gradx, jac, axes, detjac, jacinv);
    intel->ComputeRequiredData(matdata, pos);
    TPZManVector<REAL, 3> delx(dim);
    for (int d = 0; d < dim; d++)
      delx[d] = matdata.x[d] - xcenter[d];
    for (int ish = 0; ish < nshape; ish++)
    {
      for (int d = 0; d < dim; d++)
      {
        TPZManVector<REAL, 3> shapevec(dim, 0.);
        shapevec[d] = matdata.phi(ish);
        B(d, ish * dim + d) += matdata.phi(ish, 0) * weight * detjac;
        REAL mom = -delx[1] * shapevec[0] + delx[0] * shapevec[1];
        B(2, ish * dim + d) += weight * detjac * mom;
      }
    }
  }
  delete intrule;
}

/// @brief  compute orthogonalizing restraint
void TPZH1HybridApproxCreator::ComputeOrthogonalizingR(TPZFMatrix<REAL> &B, TPZFMatrix<REAL> &Restraint)
{
  if (B.Rows() != 3)
    DebugStop();
  TPZFNMatrix<9, REAL> diag(3, 3);
  int64_t ncol = B.Cols();
  TPZFMatrix<REAL> b2(3, ncol - 3);
  B.GetSub(0, 3, 3, ncol - 3, b2);
  B.GetSub(0, 0, 3, 3, diag);
  diag.Solve_LU(&b2);
  b2 *= -1.;
  Restraint.Redim(ncol, ncol);
  Restraint.Identity();
  Restraint.PutSub(0, 3, b2);
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
  for(int64_t el = 0; el<nel; el++) {
    TPZCompEl *cel = mcmesh->Element(el);
    if(!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    if(!gel) continue;
    if(gel->MaterialId() == fHybridizationData.fLagrangeMatId) {
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
