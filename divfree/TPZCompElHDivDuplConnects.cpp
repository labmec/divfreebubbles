#include "TPZCompElHDivDuplConnects.h"
#include "TPZCompElHDivDuplConnectsBound.h"
#include "TPZMaterial.h"
#include "TPZShapeHDiv.h"
#include "TPZShapeHDivConstant.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElHDiv");
static TPZLogger loggerdiv("pz.mesh.tpzinterpolatedelement.divide");
#endif

template<class TSHAPE>
TPZCompElHDivDuplConnects<TSHAPE>::TPZCompElHDivDuplConnects(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam) :
TPZRegisterClassId(&TPZCompElHDivDuplConnects::ClassId), TPZCompElHDiv<TSHAPE>(mesh,gel,hdivfam), fSideOrient(TSHAPE::NFacets,1) {
    
    std::cout << "Connects before = " << this->fConnectIndexes << std::endl;
    this->fConnectIndexes.Resize(TSHAPE::NFacets*2+1);

    //Reorder the connects
    auto prevCon = this->fConnectIndexes;
    for (int i = 0; i < TSHAPE::NFacets; i++)
    {
        this->fConnectIndexes[2*i  ] = prevCon[i];    
        this->fConnectIndexes[2*i+1] = prevCon[i+TSHAPE::NFacets+1];
    }
    this->fConnectIndexes[TSHAPE::NFacets*2] = prevCon[TSHAPE::NFacets];
    std::cout << "Connects after = " << this->fConnectIndexes << std::endl;
}


template<class TSHAPE>
int TPZCompElHDivDuplConnects<TSHAPE>::NConnects() const {
	return this->fConnectIndexes.size();
}

template<class TSHAPE>
int TPZCompElHDivDuplConnects<TSHAPE>::NSideConnects(int side) const{
    // o erro provavelmente esta aqui
	if(TSHAPE::SideDimension(side)<= this->Dimension()-2) return 0;
	if(TSHAPE::SideDimension(side)== this->Dimension()-1) return 2;
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
int TPZCompElHDivDuplConnects<TSHAPE>::NConnectShapeF(int connect, int order)const
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets*2) {
        DebugStop();
    }
#endif
    if (connect >= 2*TSHAPE::NFacets+1) return 0;

    switch (this->fhdivfam)
    {
    case HDivFamily::EHDivStandard:
        {
            int conCorrect = connect/2;
            int res = connect % 2;
            int nshape = TPZShapeHDiv<TSHAPE>::ComputeNConnectShapeF(conCorrect,order); 
            if (res == 1){ 
                nshape -= 1;
            } else {
                if (connect != 2*TSHAPE::NFacets){
                    nshape = 1;
                }
            }
            return nshape;   
        }
        break;
    case HDivFamily::EHDivConstant:
        {
            int conCorrect = connect/2;
            int res = connect % 2;
            int nshape = TPZShapeHDivConstant<TSHAPE>::ComputeNConnectShapeF(conCorrect,order);
            if (res == 1){ 
                nshape -= 1;
            } else {
                if (connect != 2*TSHAPE::NFacets){
                    nshape = 1;
                }
            }
            return nshape;
        }
        break;
    
    default:
        return -1;
        break;
    }
    return -1;
 }


template<class TSHAPE>
int64_t TPZCompElHDivDuplConnects<TSHAPE>::ConnectIndex(int con) const{
#ifndef PZNODEBUG
	if(con<0 || con > TSHAPE::NFacets*2) {
		std::cout << "TPZCompElHDivDuplConnects::ConnectIndex wrong parameter connect " << con <<
		" NConnects " << TSHAPE::NFacets << std::endl;
		DebugStop();
		return -1;
	}

#endif

	return this->fConnectIndexes[con];
}


template<class TSHAPE>
int TPZCompElHDivDuplConnects<TSHAPE>::SideConnectLocId(int node,int side) const {
    if (TSHAPE::Dimension == 2){
        return 2*(side-TSHAPE::NCornerNodes);
    } else if (TSHAPE::Dimension == 3){
        return 2*(side-(TSHAPE::NSides-TSHAPE::NumSides(TSHAPE::Dimension-1)-1));
    } else {
        DebugStop();
    }
    return -1;
}


template<class TSHAPE>
void TPZCompElHDivDuplConnects<TSHAPE>::SetConnectIndex(int i, int64_t connectindex)
{
	this->fConnectIndexes[i] = connectindex;
}

#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzshapetetra.h"
using namespace pzshape;

// template class TPZRestoreClass< TPZCompElHDivDuplConnects<TPZShapeLinear>>;
// template class TPZRestoreClass< TPZCompElHDivDuplConnects<TPZShapeTriang>>;
// template class TPZRestoreClass< TPZCompElHDivDuplConnects<TPZShapeQuad>>;
// template class TPZRestoreClass< TPZCompElHDivDuplConnects<TPZShapeCube>>;
// template class TPZRestoreClass< TPZCompElHDivDuplConnects<TPZShapeTetra>>;

// template class TPZCompElHDivDuplConnects<TPZShapeLinear>;
template class TPZCompElHDivDuplConnects<TPZShapeTriang>;
template class TPZCompElHDivDuplConnects<TPZShapeQuad>;
template class TPZCompElHDivDuplConnects<TPZShapeTetra>;
template class TPZCompElHDivDuplConnects<TPZShapeCube>;

//BC elements
TPZCompEl * CreateHDivDuplConnectsBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElHDivDuplConnectsBound< TPZShapeLinear>(mesh,gel,hdivfam);
}
TPZCompEl * CreateHDivDuplConnectsBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElHDivDuplConnectsBound< TPZShapeQuad>(mesh,gel,hdivfam);
}
TPZCompEl * CreateHDivDuplConnectsBoundTriangEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
    return new TPZCompElHDivDuplConnectsBound< TPZShapeTriang>(mesh,gel,hdivfam);
}



//Volumetric elements
TPZCompEl * CreateHDivDuplConnectsQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivDuplConnects< TPZShapeQuad>(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivDuplConnectsTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivDuplConnects< TPZShapeTriang >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivDuplConnectsCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivDuplConnects< TPZShapeCube >(mesh,gel,hdivfam);
}

TPZCompEl * CreateHDivDuplConnectsTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam) {
	return new TPZCompElHDivDuplConnects< TPZShapeTetra >(mesh,gel,hdivfam);
}


