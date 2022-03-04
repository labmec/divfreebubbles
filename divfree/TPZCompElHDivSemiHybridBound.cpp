#include "TPZCompElHDivSemiHybridBound.h"
#include "TPZMaterial.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeHDivConstantBound.h"
#include "pzlog.h"
#include "pzcmesh.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHDiv");
static TPZLogger loggerdiv("pz.mesh.tpzinterpolatedelement.divide");
#endif

template<class TSHAPE>
TPZCompElHDivSemiHybridBound<TSHAPE>::TPZCompElHDivSemiHybridBound(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElHDivSemiHybridBound::ClassId), TPZCompElHDivBound2<TSHAPE>(mesh,gel,hdivfam), fSideOrient(TSHAPE::NFacets,1) {
    
    std::cout << "ConIndex " << this->fConnectIndexes << std::endl;
    this->fConnectIndexes.Resize(2);
    //Allocate new connects for the faces
    this->fConnectIndexes[1] =  this->fConnectIndexes[0] + 5;//mesh.AllocateNewConnect(0,1,pOrder);
    // this->fConnectIndexes[1]= this->CreateMidSideConnect(TSHAPE::NSides-1);;
    // mesh.ConnectVec()[this->fConnectIndexes[1]].IncrementElConnected();

    
    std::cout << "Connect " << this->ConnectVec() << std::endl;

    constexpr auto nNodes = TSHAPE::NCornerNodes;
    auto ncon = 2;//this->NConnects();
    for(int icon = 0; icon < ncon; icon++){
        const int connect = icon;//this->MidSideConnectLocId(icon+1);
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


    for (int i = 0; i<NConnects(); i++)
    {
        auto con = this->ConnectVec()[i];
        std::cout << "C = " << con << std::endl;
    }
    
    // for (int i = 0; i < this->Mesh()->ConnectVec().Size(); i++)
    // {
    //     std::cout << "Connect " << i << " " << this->Mesh()->Connect(i) << std::endl;
    // }
    

}



template<class TSHAPE>
int TPZCompElHDivSemiHybridBound<TSHAPE>::NConnects() const {
	return 2;
}

template<class TSHAPE>
int TPZCompElHDivSemiHybridBound<TSHAPE>::NSideConnects(int side) const{
	if(side == TSHAPE::NSides-1)
	{
		return 2;
	}
	return 0;
}

template<class TSHAPE>
int TPZCompElHDivSemiHybridBound<TSHAPE>::NConnectShapeF(int connect, int connectorder)const
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif

    if (connect > 0) return 0;

	switch (this->fhdivfam)
    {
    case HDivFamily::EHDivStandard:
        if(connect == 0)
        {
            if(connectorder == 0) return 1;
            TPZManVector<int,22> order(TSHAPE::NSides-TSHAPE::NCornerNodes,connectorder);
            return TSHAPE::NShapeF(order);
        }    
        break;
    case HDivFamily::EHDivConstant:
        return TPZShapeHDivConstantBound<TSHAPE>::ComputeNConnectShapeF(connect,connectorder);
        break;
    
    default:
        DebugStop();//You shoud choose an approximation space
        break;
    }
    
    return -1;
 }


// NAO TESTADO
template<class TSHAPE>
void TPZCompElHDivSemiHybridBound<TSHAPE>::SetSideOrient(int side, int sideorient)
{
    this->fSideOrient = sideorient;
}

// NAO TESTADO
template<class TSHAPE>
int TPZCompElHDivSemiHybridBound<TSHAPE>::GetSideOrient(int side)
{
    return this->fSideOrient[0];
}




template<class TSHAPE>
int64_t TPZCompElHDivSemiHybridBound<TSHAPE>::ConnectIndex(int con) const{
#ifndef PZNODEBUG
	if(con<0 || con > 2) {
		std::cout << "TPZCompElHDivSemiHybridBound::ConnectIndex wrong parameter connect " << con <<
		" NConnects " << TSHAPE::NFacets << std::endl;
		DebugStop();
		return -1;
	}

#endif

	return this->fConnectIndexes[con];
}

template<class TSHAPE>
int TPZCompElHDivSemiHybridBound<TSHAPE>::SideConnectLocId(int node, int side) const
{
	if(side == TSHAPE::NSides-1 && node <2)
	{
		return node;
	}
	else{
	return -1;
	}
	
}


#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
using namespace pzshape;

template class TPZCompElHDivSemiHybridBound<TPZShapeLinear>;



