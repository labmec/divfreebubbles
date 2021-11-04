/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHdivBC3D.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix");
#endif

template<class TSHAPE>
TPZCompElKernelHDivBC3D<TSHAPE>::TPZCompElKernelHDivBC3D(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElKernelHDivBC3D::ClassId), TPZCompElHCurlNoGrads<TSHAPE>(mesh,gel,index)  {

}

template<class TSHAPE>
void TPZCompElKernelHDivBC3D<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZCompElHCurlNoGrads<TSHAPE>::InitMaterialData(data);

    
    // int nshape = this->NShapeF();    
    // int64_t size = nshape*3;//(TSHAPE::Dimension);
    // data.fVecShapeIndex.Resize(size);
    // // auto size = data.fVecShapeIndex.size();
    
    // for (int i=0; i<size; i++) {
	// 	data.fVecShapeIndex[i] = std::make_pair(i,1);
    // }

}

template<class TSHAPE>
void TPZCompElKernelHDivBC3D<TSHAPE>::ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi){

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    TPZCompElHCurlNoGrads<TSHAPE>::ComputeRequiredData(data,qsi);
    data.fNeedsSol = needsol;
    
    // TPZFNMatrix<220,REAL> dphix(3,data.dphix.Cols());
    // TPZFMatrix<REAL> &dphi = data.dphix;
    // TPZAxesTools<REAL>::Axes2XYZ(dphi, dphix, data.axes);
    

    data.phi.Zero();
    if (data.phi.Rows()>1){
      for (int i = 0; i < data.phi.Rows(); i++){
		data.phi(i,0) = -data.curlphi(0,i) * fSideOrient;
	  }
    }

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        //	this->Print(sout);
        sout << "\nVecshape = " << data.fVecShapeIndex << std::endl;
        sout << "MASTER = " << data.fMasterDirections << std::endl;
        sout << "Phi = " << data.phi << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
        
    }
#endif

}//void

template<class TSHAPE>
void TPZCompElKernelHDivBC3D<TSHAPE>::SetSideOrient(int orient){
    fSideOrient = orient;
}

template<class TSHAPE>
int TPZCompElKernelHDivBC3D<TSHAPE>::GetSideOrient(){
    return fSideOrient;
}

// template<class TSHAPE>
// int TPZCompElKernelHDivBC3D<TSHAPE>::NConnectShapeF(int icon, int order) const{
  
//         const int side = icon;// + TSHAPE::NCornerNodes;
//     // #ifdef PZDEBUG
//     // if (side < TSHAPE::NCornerNodes || side >= TSHAPE::NSides) {
//     //     DebugStop();
//     // }
//     // #endif
//     if(order == 0) {
//         PZError<<__PRETTY_FUNCTION__
//             <<"\nERROR: polynomial order not compatible.\nAborting..."
//             <<std::endl;
//         DebugStop();
//         return 0;
//     }
//     const auto nFaces = TSHAPE::Dimension < 2 ? 0 : TSHAPE::NumSides(2);
//     const auto nEdges = TSHAPE::NumSides(1);
//     const int nShapeF = [&](){
//         if (side < TSHAPE::NCornerNodes + nEdges) {//edge connect
//         return 1;
//         }
//         else if(side < TSHAPE::NCornerNodes + nEdges + nFaces){//face connect
//         switch(TSHAPE::Type(side)){
//         case ETriangle://triangular face
//             /**
//              we remove one internal function for each h1 face function of order k+1
//             since there are (k-1)(k-2)/2 functions per face in a face with order k,
//             we remove k(k-1)/2.
//             so:
//             (k-1)*(k+1)-k*(k-1)/2
//             */
//             return (order - 1) * (order+2) / 2;
//         default:
//             PZError<<__PRETTY_FUNCTION__<<" error. Not yet implemented"<<std::endl;
//             DebugStop();
//             return 0;
//         }
//         }
//         else{//internal connect (3D element only)
//         if constexpr (TSHAPE::Type() == ETetraedro){
//             return (order-1)*(order-2)*(2*order+3)/6;
//         }
//         return 0;
//         }
//     }();
//     return nShapeF;


// }



#include "pzshapelinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgraphelq2dd.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelt2dmapped.h"

using namespace pztopology;

// #include "tpzpoint.h"
// #include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"

// #include "TPZCompElHCurl.h"

using namespace pzgeom;
using namespace pzshape;

template class TPZCompElKernelHDivBC3D<TPZShapeTriang>;
// template class TPZCompElKernelHDivBC3D<TPZShapeQuad>;
