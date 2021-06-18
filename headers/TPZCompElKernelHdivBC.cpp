/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHdivBC.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"


template<class TSHAPE>
TPZCompElKernelHDivBC<TSHAPE>::TPZCompElKernelHDivBC(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElKernelHDivBC::ClassId), TPZIntelGen<TSHAPE>(mesh,gel,index)  {

}

template<class TSHAPE>
TPZCompElKernelHDivBC<TSHAPE>::~TPZCompElKernelHDivBC(){
    this->~TPZIntelGen<TSHAPE>();
}

template<class TSHAPE>
void TPZCompElKernelHDivBC<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZIntelGen<TSHAPE>::InitMaterialData(data);

    int nshape = this->NShapeF();
    data.fVecShapeIndex.Resize(nshape);
    for (int i=0; i<nshape; i++) {
		data.fVecShapeIndex[i] = std::make_pair(i,1);
    }

}

template<class TSHAPE>
void TPZCompElKernelHDivBC<TSHAPE>::ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi){

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    TPZIntelGen<TSHAPE>::ComputeRequiredData(data,qsi);
    data.fNeedsSol = needsol;

     for (int i = 0; i < data.phi.Rows(); i++){
		data.phi(i,0) = -data.dphix(0,i);
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

template class TPZCompElKernelHDivBC<TPZShapePoint>;
template class TPZCompElKernelHDivBC<TPZShapeLinear>;
template class TPZCompElKernelHDivBC<TPZShapeTriang>;
template class TPZCompElKernelHDivBC<TPZShapeQuad>;
