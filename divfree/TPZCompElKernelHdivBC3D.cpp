/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHdivBC3D.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"


template<class TSHAPE>
TPZCompElKernelHDivBC3D<TSHAPE>::TPZCompElKernelHDivBC3D(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElKernelHDivBC3D::ClassId), TPZIntelGen<TSHAPE>(mesh,gel,index)  {

}

template<class TSHAPE>
TPZCompElKernelHDivBC3D<TSHAPE>::~TPZCompElKernelHDivBC3D(){
    this->~TPZIntelGen<TSHAPE>();
}

template<class TSHAPE>
void TPZCompElKernelHDivBC3D<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZIntelGen<TSHAPE>::InitMaterialData(data);

    int nshape = this->NShapeF();
    data.fVecShapeIndex.Resize(nshape);
    for (int i=0; i<nshape; i++) {
		data.fVecShapeIndex[i] = std::make_pair(i,1);
    }

}

template<class TSHAPE>
void TPZCompElKernelHDivBC3D<TSHAPE>::ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi){

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    TPZIntelGen<TSHAPE>::ComputeRequiredData(data,qsi);
    data.fNeedsSol = needsol;
    
    TPZFNMatrix<220,REAL> dphix(3,data.dphix.Cols());
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZAxesTools<REAL>::Axes2XYZ(dphi, dphix, data.axes);

    if (data.phi.Rows()>1){
      for (int i = 0; i < data.phi.Rows(); i++){
		data.phi(i,0) = dphi(1,i)-dphi(0,i);
	  }
    }

}//void




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

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"

#include "pzelchdivbound2.h"

using namespace pzgeom;
using namespace pzshape;

template class TPZCompElKernelHDivBC3D<TPZShapePoint>;
template class TPZCompElKernelHDivBC3D<TPZShapeLinear>;
template class TPZCompElKernelHDivBC3D<TPZShapeTriang>;
template class TPZCompElKernelHDivBC3D<TPZShapeQuad>;
