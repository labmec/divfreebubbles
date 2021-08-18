#include "TPZCompElHCurlNoGrads.h"

#include <TPZMaterial.h>
#include <pzconnect.h>
#include <pzcmesh.h>

template<class TSHAPE>
TPZCompElHCurlNoGrads<TSHAPE>::TPZCompElHCurlNoGrads() : TPZCompElHCurlFull<TSHAPE>()
{
  this->AdjustConnects();
}

template<class TSHAPE>
TPZCompElHCurlNoGrads<TSHAPE>::TPZCompElHCurlNoGrads(
  TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
  TPZCompElHCurlFull<TSHAPE>(mesh,gel,index)
{
  this->AdjustConnects();
}

template<class TSHAPE>
void TPZCompElHCurlNoGrads<TSHAPE>::AdjustConnects()
{
  constexpr auto nNodes = TSHAPE::NCornerNodes;
  constexpr auto ncon = TSHAPE::NSides - nNodes;
  for(int icon = 0; icon < ncon; icon++){
    const int connect = this->MidSideConnectLocId(icon+nNodes);
    TPZConnect &c = this->Connect(connect);
    const int nshape =this->NConnectShapeF(connect,c.Order());
    c.SetNShape(nshape);
    const auto seqnum = c.SequenceNumber();
    const int nStateVars = [&](){
        TPZMaterial * mat =this-> Material();
        if(mat) return mat->NStateVariables();
        else {
            return 1;
        }
    }();
    this-> Mesh()->Block().Set(seqnum,nshape*nStateVars);
  }
}
  
template<class TSHAPE>
void TPZCompElHCurlNoGrads<TSHAPE>::InitMaterialData(TPZMaterialData &data){
	TPZIntelGen<TSHAPE>::InitMaterialData(data);
  data.fMasterDirections = this->fMasterDirections;

  /*
    we first compute ALL the hcurl traditional functions(scalar+vector)
    then, at each integration point we will filter them.
    
    some of this code is copied from the TPZCompElHCurlFull since we want to avoid
    calls for the virtual methods.
  */
  
  //computes the index that will associate each scalar function to a constant vector field
  constexpr auto nConnects = TSHAPE::NSides - TSHAPE::NCornerNodes;
  TPZManVector<int,nConnects> connectOrders(nConnects,-1);
  int unfiltnshape = 0;
  for(auto i = 0; i < nConnects; i++){
    const auto conorder = this->EffectiveSideOrder(i + TSHAPE::NCornerNodes);
    connectOrders[i] = conorder;
    unfiltnshape += TPZCompElHCurlFull<TSHAPE>::NConnectShapeF(i, conorder);
  }


  const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
  const auto nEdges = TSHAPE::NumSides(1);
  constexpr auto nNodes = TSHAPE::NCornerNodes;

  TPZManVector<int64_t,nNodes> nodes(nNodes, 0);

  for (auto iNode = 0; iNode < nNodes; iNode++){
    nodes[iNode] = this->Reference()->NodeIndex(iNode);
  }

  
  TPZManVector<int64_t, TSHAPE::NSides - nNodes>
    firstH1ShapeFunc(TSHAPE::NSides - nNodes,0);
  
  //calculates the first shape function associated with each side of dim > 0
  TPZManVector<int,TSHAPE::NSides-nNodes> sidesH1Ord(TSHAPE::NSides - nNodes,-1);
  this->CalcH1ShapeOrders(sidesH1Ord);
  firstH1ShapeFunc[0] = nNodes;
  
  for (int iSide = nNodes + 1; iSide < TSHAPE::NSides; iSide++) {
    const int iCon = iSide - nNodes;
    firstH1ShapeFunc[iCon] =
      firstH1ShapeFunc[iCon - 1] +
      TSHAPE::NConnectShapeF(iSide - 1, sidesH1Ord[iCon-1]);

  }

  auto &indexVecShape = data.fVecShapeIndex;


  
  
  indexVecShape.Resize(unfiltnshape);
  TPZVec<unsigned int> shapeCountVec(TSHAPE::NSides - nNodes, 0);
  this->StaticIndexShapeToVec(indexVecShape, connectOrders,
                              firstH1ShapeFunc,sidesH1Ord, shapeCountVec, nodes);

  //setting the type of shape functions as vector shape functions
  data.fShapeType = TPZMaterialData::EVecShape;
}

template<class TSHAPE>
template<class TVar>
void TPZCompElHCurlNoGrads<TSHAPE>::ComputeRequiredDataT(
  TPZMaterialDataT<TVar> &data, TPZVec<REAL> &qsi)
{
  
  {
    const bool needsSol = data.fNeedsSol;
    data.fNeedsSol = false;
    TPZIntelGen<TSHAPE>::ComputeRequiredData(data,qsi);//in this method, Shape will be called
    data.fNeedsSol = needsSol;
  }

  this->ComputeDeformedDirections(data);    
  /******************************************************************************************************************
   * at this point, we already have the basis functions on the deformed element, since we have data.phi,
   * data.fVecShapeIndex and data.fDeformedDirections. Now it is time to compute the curl, which will be stored in
   * data.curlphi.
   *******************************************************************************************************************/
  TPZFMatrix<REAL> phiHCurl;
  this->ComputeShape(data, phiHCurl);
  

  constexpr auto dim{TSHAPE::Dimension};
  data.curlphi.Redim(2*dim - 3 > 0 ? 2*dim - 3 : 1, this->NShapeF());
  this->ComputeCurl<dim>(data);
  
  data.phi = phiHCurl;
  if (data.fNeedsSol) {
    this->ReallyComputeSolution(data);
  }
}

template<class TSHAPE>
template<class TVar>
void TPZCompElHCurlNoGrads<TSHAPE>::ReallyComputeSolutionT(TPZMaterialDataT<TVar> &data)
{
    this->ComputeSolutionHCurlT(data.phi, data.curlphi,
                         data.sol, data.curlsol);
}


template<class TSHAPE>
int TPZCompElHCurlNoGrads<TSHAPE>::NConnectShapeF(int icon, int order) const
{
  const int side = icon + TSHAPE::NCornerNodes;
#ifdef PZDEBUG
  if (side < TSHAPE::NCornerNodes || side >= TSHAPE::NSides) {
    DebugStop();
  }
#endif
  if(order == 0) {
    PZError<<__PRETTY_FUNCTION__
           <<"\nERROR: polynomial order not compatible.\nAborting..."
           <<std::endl;
    DebugStop();
    return 0;
  }
  const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
  const auto nEdges = TSHAPE::NumSides(1);
  const int nShapeF = [&](){
    if (side < TSHAPE::NCornerNodes + nEdges) {//edge connect
      return 1;
    }
    else if(side < TSHAPE::NCornerNodes + nEdges + nFaces){//face connect
      return 0;
      // switch(TSHAPE::Type(side)){
      // case ETriangle://triangular face
      //   return (order - 1) * (order+1);
      // case EQuadrilateral://quadrilateral face
      //   return 2 * order * (order+1);
      // default:
      //   PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
      //   DebugStop();
      //   return 0;
      // }
    }
    else{//internal connect (3D element only)
      return 0;
      // int count = 0;
      // //first we count the face-based interior functions \phi^{K,F}
      // for(int iFace = 0; iFace < nFaces; iFace++){
      //   switch(TSHAPE::Type(TSHAPE::NCornerNodes+nEdges + iFace)){
      //   case ETriangle://triangular face
      //     count +=  (order - 1) * ( order - 2)/2;
      //     break;
      //   case EQuadrilateral://quadrilateral face
      //     //we need the functions of k+1
      //     count +=  order * order;
      //     break;
      //   default:
      //     PZError<<__PRETTY_FUNCTION__<<" error."<<std::endl;
      //     DebugStop();
      //   }
      // }

      // const int nVkf = count;
      // //now we count the interior bubble functions
            
      // //number of H1 funcs
      // const auto orderH1 = TSHAPE::Type() == ECube ?
      //   order +1 : order;
            
      // const auto nFuncsH1 = TSHAPE::NConnectShapeF(side, orderH1);
      // if constexpr (TSHAPE::Type() == ECube) {
      //   //shapeorders[k] contains the polynomial orders
      //   //in xi, eta and zeta for each of the internal funcs
      //   TPZGenMatrix<int> shapeorders(nFuncsH1,3);
      //   //not really used since we are interested in internal funcs
      //   TPZManVector<int64_t,0> idvec(0);
      //   TSHAPE::SideShapeOrder(side,idvec,orderH1,shapeorders);              
      //   /*for compatibility, we need the spaces
      //     Q^{k,k+1,k+1} \times Q^{k+1,k,k+1} \times Q^{k+1,k+1,k}
      //     I am supposing that #funcs in xi-dir is equal to 
      //     #funcs in eta-dir and in zeta-dir. So I will just
      //     see how many we will have in xi-dir after skipping
      //     internals that were included in the k-1 set of funcs*/
      //   int countInt = 0;
      //   for(int ifunc = 0; ifunc < nFuncsH1; ifunc++){
      //     if(shapeorders(ifunc,0) <= order) {countInt++;}
      //   }
      //   count += 3* countInt;

      //   const auto nVki = count - nVkf;
      //   // std::stringstream sout;
      //   // sout << __PRETTY_FUNCTION__<<'\n'
      //   //      <<"\tside "<<side<<"\tcon "<<icon<<"\torder "<<order<<'\n'
      //   //      <<"\tnfuncs "<<count<<"\tnvkf "<<nVkf<<"\tnvki "<<nVki<<'\n';
      //   // std::cout<<sout.str()<<std::flush;
      // }
      // else{
      //   count += 3 * nFuncsH1;
      // }
      // return count;
    }
  }();
  return nShapeF;
}

template<class TSHAPE>
void TPZCompElHCurlNoGrads<TSHAPE>::ComputeShape(TPZMaterialData &data,
                                                 TPZFMatrix<REAL> &phiHCurl)
{
  constexpr auto dim = TSHAPE::Dimension;
  constexpr auto nNodes = TSHAPE::NCornerNodes;
  constexpr auto nConnects = TSHAPE::NSides - nNodes;
  const auto nFaces = TSHAPE::NumSides(2);
  const auto nEdges = TSHAPE::NumSides(1);
  
  //unfiltered shape functions count for each connect
  TPZManVector<int,nConnects> firstHCurlFunc(nConnects,0);
  for(auto icon = 1; icon < nConnects; icon++){
    const auto conorder = this->ConnectOrder(icon);
    firstHCurlFunc[icon] = firstHCurlFunc[icon-1] +
      TPZCompElHCurlFull<TSHAPE>::NConnectShapeF(icon,conorder);
  }

  //number of FILTERED hcurl functions
  const auto nshape = this->NShapeF();
  
  const auto &vecShapeIndex = data.fVecShapeIndex;
  const auto &deformedDirs = data.fDeformedDirections;
  const auto &phiH1 = data.phi;

        
  phiHCurl.Resize(nshape, dim);
  int fcount = 0;
        
  /*
    edge connects: we combine the lowest order edge functions in order to
    create a function with constant tangential trace
  */

  for(auto ie = 0; ie < nEdges; ie++){
    const auto firstSideShape = firstHCurlFunc[ie];
    const auto vIndex1 = vecShapeIndex[firstSideShape].first;
    const auto sIndex1 = vecShapeIndex[firstSideShape].second;
    const auto vIndex2 = vecShapeIndex[firstSideShape+1].first;
    const auto sIndex2 = vecShapeIndex[firstSideShape+1].second;
    for(auto x = 0; x < dim; x++){
      phiHCurl(fcount,x) = phiH1.GetVal(sIndex1,0) *
        deformedDirs.GetVal(x,vIndex1) +
        phiH1.GetVal(sIndex2,0) *
        deformedDirs.GetVal(x,vIndex2);
    }
    fcount++;
  }
  if constexpr (dim < 2) return;
  /*
    face connects: 
  */
  if constexpr (dim < 3) return;
  /*
    interior connects:
  */
    
}

template<class TSHAPE>
template<int DIM>
void TPZCompElHCurlNoGrads<TSHAPE>::ComputeCurl(TPZMaterialData &data)
{

  
  TPZFMatrix<REAL> unfiltCurl;
  //compute curl of ALL HCurl functions
  TPZCompElHCurlFull<TSHAPE>::ComputeCurl(
    data.fVecShapeIndex,data.dphi,this->fMasterDirections,
    data.jacobian,data.detjac,data.axes,unfiltCurl);

  constexpr auto dim = TSHAPE::Dimension;
  constexpr auto nNodes = TSHAPE::NCornerNodes;
  constexpr auto nConnects = TSHAPE::NSides - nNodes;
  const auto nFaces = TSHAPE::NumSides(2);
  const auto nEdges = TSHAPE::NumSides(1);

  //number of FILTERED hcurl functions
  const int nshape = this->NShapeF();
  //unfiltered shape functions count for each connect
  TPZManVector<int,nConnects> firstHCurlFunc(nConnects,0);
  for(auto icon = 1; icon < nConnects; icon++){
    const auto conorder = this->ConnectOrder(icon);
    firstHCurlFunc[icon] = firstHCurlFunc[icon-1] +
      TPZCompElHCurlFull<TSHAPE>::NConnectShapeF(icon,conorder);
  }
  
  constexpr auto curlDim= 2*dim - 3 > 0 ? 2*dim -3 : 1;

  auto &curlPhi = data.curlphi;
  curlPhi.Redim(curlDim,nshape);

  int fcount = 0;
  //edges
  for(auto ie = 0; ie < nEdges; ie++){
    //fss = first side shape
    const auto fss = firstHCurlFunc[ie];
    for(auto x = 0; x < curlDim; x++){
      curlPhi(x,fcount) += unfiltCurl(x,fss) + unfiltCurl(x,fss+1);
    }
    fcount++;
  }
  if constexpr (dim < 2) return;
  //face
  if constexpr (dim < 3) return;
  //interior
}

#include <pzshapetriang.h>
#include <pzshapetetra.h>


#define IMPLEMENTHCURLNOGRADS(TSHAPE)                           \
                                                                \
  template class                                                \
  TPZRestoreClass< TPZCompElHCurlNoGrads<TSHAPE> >;             \
  template class TPZCompElHCurlNoGrads<TSHAPE>;                 \
                                                                \
  template void                                                 \
  TPZCompElHCurlNoGrads<TSHAPE>::ComputeRequiredDataT<STATE>(   \
    TPZMaterialDataT<STATE> &data,TPZVec<REAL> &qsi);           \
  template void                                                 \
  TPZCompElHCurlNoGrads<TSHAPE>::ComputeRequiredDataT<CSTATE>(  \
    TPZMaterialDataT<CSTATE> &data, TPZVec<REAL> &qsi);

IMPLEMENTHCURLNOGRADS(pzshape::TPZShapeTriang)
IMPLEMENTHCURLNOGRADS(pzshape::TPZShapeTetra)

#undef IMPLEMENTHCURLNOGRADS