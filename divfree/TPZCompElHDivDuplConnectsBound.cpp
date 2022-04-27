#include "TPZCompElHDivDuplConnectsBound.h"
#include "TPZMaterial.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeHDivConstantBound.h"
#include "pzlog.h"
#include "pzcmesh.h"
#include <sstream>

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHDiv");
static TPZLogger loggerdiv("pz.mesh.tpzinterpolatedelement.divide");
#endif

template<class TSHAPE>
TPZCompElHDivDuplConnectsBound<TSHAPE>::TPZCompElHDivDuplConnectsBound(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElHDivDuplConnectsBound::ClassId), TPZCompElHDivBound2<TSHAPE>(mesh,gel,hdivfam), fSideOrient(TSHAPE::NFacets,1) {
    
    this->fConnectIndexes.Resize(2);//Change it to 3D

    if (TSHAPE::Dimension == 2) DebugStop();
}

template<class TSHAPE>
int TPZCompElHDivDuplConnectsBound<TSHAPE>::NConnects() const {
	return 2;
}

template<class TSHAPE>
int TPZCompElHDivDuplConnectsBound<TSHAPE>::NSideConnects(int side) const{
	if(side == TSHAPE::NSides-1)
	{
		return 2;
	}
	return 0;
}

template<class TSHAPE>
int TPZCompElHDivDuplConnectsBound<TSHAPE>::NConnectShapeF(int connect, int connectorder)const
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif

    // if (connect > 0) return 0;

	switch (this->fhdivfam)
    {
    case HDivFamily::EHDivStandard:
        {
            if(connectorder == 0) return 1;
            TPZManVector<int,22> order(TSHAPE::NSides-TSHAPE::NCornerNodes,connectorder);
            int nshape = TSHAPE::NShapeF(order);
            int conCorrect = connect/2;
            int res = connect % 2;
            if (res == 1){ 
                // nshape = 0;
                nshape -= 1;
            } else {
                nshape = 1;
            }

            return nshape;
        }    
        break;
    case HDivFamily::EHDivConstant:
        {
            int conCorrect = connect/2;
            int res = connect % 2;
            int nshape = TPZShapeHDivConstantBound<TSHAPE>::ComputeNConnectShapeF(connect,connectorder);
            // if (res == 1) nshape = 0;
            if (res == 1){ 
                // nshape = 0;
                nshape -= 1;
            } else {
                nshape = 1;
            }

            return nshape;
        }
        break;
    
    default:
        DebugStop();//You shoud choose an approximation space
        break;
    }
    
    return -1;
 }

template<class TSHAPE>
int64_t TPZCompElHDivDuplConnectsBound<TSHAPE>::ConnectIndex(int con) const{
#ifndef PZNODEBUG
	if(con<0 || con > 2) {
		std::cout << "TPZCompElHDivDuplConnectsBound::ConnectIndex wrong parameter connect " << con <<
		" NConnects " << TSHAPE::NFacets << std::endl;
		DebugStop();
		return -1;
	}

#endif

	return this->fConnectIndexes[con];
}

template<class TSHAPE>
int TPZCompElHDivDuplConnectsBound<TSHAPE>::SideConnectLocId(int node, int side) const
{
	if(side == TSHAPE::NSides-1 && node <2)
	{
		return node;
	}
	else{
	return -1;
	}
	
}

template<class TSHAPE>
void TPZCompElHDivDuplConnectsBound<TSHAPE>::SetConnectIndex(int i, int64_t connectindex)
{
	this->fConnectIndexes[i] = connectindex;
}


#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
using namespace pzshape;

template class TPZCompElHDivDuplConnectsBound<TPZShapeLinear>;



