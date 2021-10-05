/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHdivBC3D.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"


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

    // std::cout << "Phi = " << data.phi << std::endl;
    // std::cout << "ELEMENT MATERIAL ID = " << this->Reference()->MaterialId() << std::endl; 

}//void

template<class TSHAPE>
void TPZCompElKernelHDivBC3D<TSHAPE>::SetSideOrient(int orient){
    fSideOrient = orient;
}

template<class TSHAPE>
int TPZCompElKernelHDivBC3D<TSHAPE>::GetSideOrient(){
    return fSideOrient;
}




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
// #include "tpzquadrilateral.h"
#include "tpztriangle.h"

// #include "TPZCompElHCurl.h"

using namespace pzgeom;
using namespace pzshape;

// template class TPZCompElKernelHDivBC3D<TPZShapePoint>;
// template class TPZCompElKernelHDivBC3D<TPZShapeLinear>;
template class TPZCompElKernelHDivBC3D<TPZShapeTriang>;
// template class TPZCompElKernelHDivBC3D<TPZShapeQuad>;
