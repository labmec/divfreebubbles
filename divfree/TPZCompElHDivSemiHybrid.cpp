#include "TPZCompElHDivSemiHybrid.h"
#include "TPZCompElHDivSemiHybridBound.h"
#include "TPZMaterial.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeHDivConstant.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHDiv");
static TPZLogger loggerdiv("pz.mesh.tpzinterpolatedelement.divide");
#endif

template<class TSHAPE>
TPZCompElHDivSemiHybrid<TSHAPE>::TPZCompElHDivSemiHybrid(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElHDivSemiHybrid::ClassId), TPZCompElHDiv<TSHAPE>(mesh,gel,hdivfam), fSideOrient(TSHAPE::NFacets,1) {
    
    this->fConnectIndexes.Resize(TSHAPE::NFacets*2+1);
    //Allocate new connects for the faces
    for (int i = 0; i < TSHAPE::NFacets; i++)
    {
        auto pOrder = this->ConnectOrder(i);
        auto c = this->Connect(i);
        // this->fConnectIndexes[TSHAPE::NFacets+1+i]= this->Mesh()->AllocateNewConnect(0,1,pOrder);
        this->fConnectIndexes[TSHAPE::NFacets+1+i]= this->Mesh()->AllocateNewConnect(c);
        // this->Connect(ConnectVec()[TSHAPE::NFacets+1+i]) = this->Connect(ConnectVec()[i]);
        // this->Connect(ConnectVec()[TSHAPE::NFacets+1+i]).SetNShape(0);
        // mesh.ConnectVec()[this->fConnectIndexes[TSHAPE::NFacets+1+i]].IncrementElConnected();
    } 
    
    std::cout << "ConIndex " << this->fConnectIndexes << std::endl;
    std::cout << "Connect " << this->ConnectVec() << std::endl;

    constexpr auto nNodes = TSHAPE::NCornerNodes;
    auto ncon = this->NConnects();
    // if (TSHAPE::Dimension == 2) ncon =TSHAPE::NSides - nNodes;
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
        TPZConnect con = this->Connect(i);
        std::cout << "C = " << con << std::endl;
        std::cout << "NShape = " << con.NShape() << std::endl;
    }
    
    // for (int i = 0; i < this->Mesh()->ConnectVec().Size(); i++)
    // {
    //     std::cout << "Connect " << i << " " << this->Mesh()->Connect(i) << std::endl;
    // }
    

}



template<class TSHAPE>
int TPZCompElHDivSemiHybrid<TSHAPE>::NConnects() const {
	return TSHAPE::NFacets*2 + 1;
}

template<class TSHAPE>
int TPZCompElHDivSemiHybrid<TSHAPE>::NSideConnects(int side) const{
	if(TSHAPE::SideDimension(side)<= this->Dimension()-2) return 0;
	if(TSHAPE::SideDimension(side)== this->Dimension()-1) return 1;
	if(TSHAPE::SideDimension(side)== this->Dimension()) {
        int ncon = 1;
        return ncon;
    }
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << "Side: " << side <<"unhandled case ";
		LOGPZ_ERROR(logger,sout.str())
	}
#endif
	return -1;

}

template<class TSHAPE>
int TPZCompElHDivSemiHybrid<TSHAPE>::NConnectShapeF(int connect, int order)const
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets*2) {
        DebugStop();
    }
#endif
    if (connect >= TSHAPE::NFacets+1) return 0;

    switch (this->fhdivfam)
    {
    case HDivFamily::EHDivStandard:
        return TPZShapeHDiv<TSHAPE>::ComputeNConnectShapeF(connect,order);    
        break;
    case HDivFamily::EHDivConstant:
        return TPZShapeHDivConstant<TSHAPE>::ComputeNConnectShapeF(connect,order);
        break;
    
    default:
        return -1;
        break;
    }
    return -1;
 }


template<class TSHAPE>
int64_t TPZCompElHDivSemiHybrid<TSHAPE>::ConnectIndex(int con) const{
#ifndef PZNODEBUG
	if(con<0 || con > TSHAPE::NFacets*2) {
		std::cout << "TPZCompElHDivSemiHybrid::ConnectIndex wrong parameter connect " << con <<
		" NConnects " << TSHAPE::NFacets << std::endl;
		DebugStop();
		return -1;
	}

#endif

	return this->fConnectIndexes[con];
}







#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
using namespace pzshape;

// template class TPZRestoreClass< TPZCompElHDivSemiHybrid<TPZShapeLinear>>;
// template class TPZRestoreClass< TPZCompElHDivSemiHybrid<TPZShapeTriang>>;
// template class TPZRestoreClass< TPZCompElHDivSemiHybrid<TPZShapeQuad>>;
// template class TPZRestoreClass< TPZCompElHDivSemiHybrid<TPZShapeCube>>;
// template class TPZRestoreClass< TPZCompElHDivSemiHybrid<TPZShapeTetra>>;

// template class TPZCompElHDivSemiHybrid<TPZShapeLinear>;
template class TPZCompElHDivSemiHybrid<TPZShapeTriang>;
template class TPZCompElHDivSemiHybrid<TPZShapeQuad>;
template class TPZCompElHDivSemiHybrid<TPZShapeTetra>;
template class TPZCompElHDivSemiHybrid<TPZShapeCube>;


TPZCompEl * CreateHDivSemiHybridLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElHDivSemiHybridBound< TPZShapeLinear>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivSemiHybridQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivSemiHybrid< TPZShapeQuad>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivSemiHybridTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivSemiHybrid< TPZShapeTriang >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivSemiHybridCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivSemiHybrid< TPZShapeCube >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivSemiHybridTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivSemiHybrid< TPZShapeTetra >(mesh,gel,hdivfam);
}


