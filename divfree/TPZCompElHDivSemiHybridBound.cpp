#include "TPZCompElHDivSemiHybridBound.h"
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
TPZCompElHDivSemiHybridBound<TSHAPE>::TPZCompElHDivSemiHybridBound(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElHDivSemiHybridBound::ClassId), TPZCompElHDivBound2<TSHAPE>(mesh,gel,hdivfam), fSideOrient(TSHAPE::NFacets,1) {
    
    this->fConnectIndexes.Resize(2);//Change it to 3D

    //Get the index of the additional created connect
    this->fConnectIndexes[1] = this->CreateMidSideConnect(TSHAPE::NSides-1);

    // for (int i = 0; i<NConnects(); i++)
    // {
    //     TPZConnect con = this->Connect(i);
    //     std::cout << "C = " << con << std::endl;
    //     std::cout << "NShape = " << con.NShape() << std::endl;
    // }
    // std::cout << "\n\n\n";
    

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
        {
            int conCorrect = connect/2;
            int res = connect % 2;
            int nshape = TPZShapeHDivConstantBound<TSHAPE>::ComputeNConnectShapeF(connect,connectorder);
            if (res == 1) nshape = 0;
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

template<class TSHAPE>
int64_t TPZCompElHDivSemiHybridBound<TSHAPE>::CreateMidSideConnect(int side) {
    TPZCompMesh *cmesh = this->Mesh();
    TPZMaterial * mat = this->Material();
#ifdef PZDEBUG
    if (!mat) {
        std::cout << __PRETTY_FUNCTION__ << " no material associated with matid " << this->Reference()->MaterialId() << std::endl;
    }
#endif
    int nvar = 1;
    if (mat) nvar = mat->NStateVariables();
    int64_t newnodeindex = -1;
    int64_t il;
    int64_t nodloc = this->MidSideConnectLocId(side);

    TPZStack<TPZCompElSide> elvec;
    TPZCompElSide thisside(this, side);

    // Connect looks for a connecting element of equal or lower level
    TPZInterpolatedElement *cel = 0;
    int side_neig = 0;
    thisside.EqualLevelElementList(elvec, 1, 0);
    int64_t nelem = elvec.NElements();
    // find an element in the list which is interpolated
    if (nelem) {
        cel = dynamic_cast<TPZInterpolatedElement *> (elvec[0].Element());
        side_neig = elvec[0].Side();
    }
    int64_t newnodecreated = 0;
    if (cel) {
        auto cind = cel->MidSideConnectLocId(side_neig);
        newnodeindex = cel->ConnectIndex(cind+1);
        return newnodeindex;    
    } else {
        DebugStop();
    }
    return newnodeindex;
}


#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
using namespace pzshape;

template class TPZCompElHDivSemiHybridBound<TPZShapeLinear>;



