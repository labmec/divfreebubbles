/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElHDivConstantBC.h"
#include "TPZMaterialData.h"
#include "TPZMaterialDataT.h"
#include "TPZShapeHDivKernel2DBound.h"
#include "TPZShapeHDivConstantBound.h"

#include "pzcmesh.h"

template<class TSHAPE>
TPZCompElHDivConstantBC<TSHAPE>::TPZCompElHDivConstantBC(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index, int shapetype) :
TPZRegisterClassId(&TPZCompElHDivConstantBC::ClassId), TPZCompElHDivBound2<TSHAPE>(mesh,gel,index), fShapeType(shapetype)  {
    this->AdjustConnects();
}

template<class TSHAPE>
TPZCompElHDivConstantBC<TSHAPE>::~TPZCompElHDivConstantBC(){
    this->~TPZCompElHDivBound2<TSHAPE>();
}

template<class TSHAPE>
void TPZCompElHDivConstantBC<TSHAPE>::AdjustConnects()
{
    constexpr auto nNodes = TSHAPE::NCornerNodes;
    constexpr auto ncon = TSHAPE::NSides - nNodes;
    for(int icon = 0; icon < ncon; icon++){
        const int connect = icon;//this->MidSideConnectLocId(icon+1);
        TPZConnect &c = this->Connect(connect);
        TPZVec<int> pOrder;
        this->GetInterpolationOrder(pOrder);
        int64_t index;
        c.SetOrder(pOrder[0],index);
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
        std::cout << "Connect = " << c << std::endl;
        std::cout << "NSHAPE = " << c.NShape() << std::endl;
    }
}


template<class TSHAPE>
void TPZCompElHDivConstantBC<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZCompElHDivBound2<TSHAPE>::InitMaterialData(data);
    // this->AdjustConnects();
    int nshape = this->NShapeF();
    data.fVecShapeIndex.Resize(nshape);
    for (int i=0; i<nshape; i++) {
		data.fVecShapeIndex[i] = std::make_pair(i,1);
    }

}

template<class TSHAPE>
void TPZCompElHDivConstantBC<TSHAPE>::ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data){

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    data.fNeedsSol = needsol;
    // TPZCompElHDivBound2<TSHAPE>::ComputeShape(qsi,data);

    TPZShapeData &shapedata = data;
    int nshape = this->NShapeF();
    data.phi.Resize(nshape,1);

    TPZFMatrix<REAL> auxPhi(nshape,1);
    auxPhi.Zero();

    TPZShapeHDivConstantBound<TSHAPE>::Shape(qsi, shapedata, auxPhi);

    for (int i = 0; i < data.phi.Rows(); i++){
        data.phi(i,0) = auxPhi(i,0) / data.detjac;
    }
 
}//void

template<class TSHAPE>
int TPZCompElHDivConstantBC<TSHAPE>::NConnectShapeF(int connect, int connectorder) const
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif
    return TPZShapeHDivConstantBound<TSHAPE>::ComputeNConnectShapeF(connect,connectorder);
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

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"

#include "pzelchdivbound2.h"

using namespace pzgeom;
using namespace pzshape;

// template class TPZCompElKernelHDivBC<TPZShapePoint>;
template class TPZCompElHDivConstantBC<TPZShapeLinear>;
// template class TPZCompElKernelHDivBC<TPZShapeTriang>;
// template class TPZCompElKernelHDivBC<TPZShapeQuad>;
