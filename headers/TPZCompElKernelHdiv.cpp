/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHdiv.h"

#include "pzcmesh.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "TPZMaterial.h"
#include "TPZMatSingleSpace.h"
#include "pzlog.h"
#include "pzgeoquad.h"
#include "TPZShapeDisc.h"
#include "TPZCompElDisc.h"
#include "TPZMaterialDataT.h"
#include "pzshapepiram.h"
#include "tpzline.h"
#include "tpztriangle.h"


#include "pzshtmat.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.TPZCompElKernelHDiv");
static TPZLogger loggerdiv("pz.mesh.tpzinterpolatedelement.divide");
#endif

using namespace std;


template<class TSHAPE>
TPZCompElKernelHDiv<TSHAPE>::TPZCompElKernelHDiv(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElKernelHDiv::ClassId), TPZIntelGen<TSHAPE>(mesh,gel,index) {

}

template<class TSHAPE>
TPZCompElKernelHDiv<TSHAPE>::TPZCompElKernelHDiv(TPZCompMesh &mesh, const TPZCompElKernelHDiv<TSHAPE> &copy) :
TPZRegisterClassId(&TPZCompElKernelHDiv::ClassId),
TPZIntelGen<TSHAPE>(mesh,copy)
{
	this-> fPreferredOrder = copy.fPreferredOrder;
    this->fConnectIndexes = copy.fConnectIndexes;
}

template<class TSHAPE>
TPZCompElKernelHDiv<TSHAPE>::TPZCompElKernelHDiv(TPZCompMesh &mesh,
									             const TPZCompElKernelHDiv<TSHAPE> &copy,
									             std::map<int64_t,int64_t> & gl2lcConMap,
									             std::map<int64_t,int64_t> & gl2lcElMap) :
TPZRegisterClassId(&TPZCompElKernelHDiv::ClassId),
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap){

}

template<class TSHAPE>
TPZCompElKernelHDiv<TSHAPE>::TPZCompElKernelHDiv() :
TPZRegisterClassId(&TPZCompElKernelHDiv::ClassId), TPZIntelGen<TSHAPE>()
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NSides;i++) {
		this-> fConnectIndexes[i] = -1;
	}
}

template<class TSHAPE>
TPZCompElKernelHDiv<TSHAPE>::~TPZCompElKernelHDiv(){
    this->~TPZIntelGen<TSHAPE>();
}

template<class TSHAPE>
MElementType TPZCompElKernelHDiv<TSHAPE>::Type() {
	return TSHAPE::Type();
}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::SetConnectIndex(int i, int64_t connectindex){
	this-> SetConnectIndex(i,connectindex);
}

template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::NConnectShapeF(int connect, int order)const
{
    if(connect < TSHAPE::NCornerNodes) return TSHAPE::NConnectShapeF(connect,0);
	if(order < 0) return 0;
    int nshape = TSHAPE::NConnectShapeF(connect, order);
#ifdef PZDEBUG
    if(nshape < 0 )
    {
        nshape = TSHAPE::NConnectShapeF(connect, order);
        DebugStop();
    }
#endif
	return nshape;
 }

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::SetIntegrationRule(int ord) {
	TPZManVector<int,3> order(TSHAPE::Dimension,ord);
	this->fIntRule.SetOrder(order);
}

template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::NSideConnects(int side) const{
	return TSHAPE::NContainedSides(side);
}

template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::SideConnectLocId(int node,int side) const {
    return TSHAPE::ContainedSideLocId(side,node);
}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
    this->GetInterpolationOrder(ord);
}


template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::PreferredSideOrder(int side) {
	if(side < TSHAPE::NCornerNodes) return 0;
	if(side<TSHAPE::NSides) {
		int order =this->fPreferredOrder;
		return this->AdjustPreferredSideOrder(side,order);
	}
	PZError << "TPZIntelgen::PreferredSideOrder called for side = " << side << "\n";
	return 0;
}

template<class TSHAPE>
int64_t TPZCompElKernelHDiv<TSHAPE>::ConnectIndex(int con) const{
#ifndef NODEBUG
	if(con<0 || con>= this->NConnects()) {
		std::cout << "TPZIntelgen::ConnectIndex wrong parameter con " << con <<
		" NSides " << TSHAPE::NSides << " NConnects " << this->NConnects() << std::endl;
		DebugStop();
	}
#endif
	return this->fConnectIndexes[con];
}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::SetPreferredOrder(int order)
{
    TPZIntelGen<TSHAPE>:: SetPreferredOrder(order);
	//this->fPreferredOrder = order;
}

template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::ConnectOrder(int connect) const{
	return this->Connect(connect).Order();
}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
{	
    TPZManVector<int64_t,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);
	TPZManVector<int, TSHAPE::NSides-TSHAPE::NCornerNodes+1> ord(TSHAPE::NSides-TSHAPE::NCornerNodes,0);
	int i;
	TPZGeoEl *ref = this->Reference();
	for(i=0; i<TSHAPE::NCornerNodes; i++) {
		id[i] = ref->NodePtr(i)->Id();
	}
	for(i=0; i<TSHAPE::NSides-TSHAPE::NCornerNodes; i++) {
		ord[i] = this->Connect(i+TSHAPE::NCornerNodes).Order();
	}
	TSHAPE::Shape(pt,id,ord,phi,dphi);
}

template<class TSHAPE>
template<class TVar>
void TPZCompElKernelHDiv<TSHAPE>::ComputeRequiredDataT(TPZMaterialDataT<TVar> &data,
                                                TPZVec<REAL> &qsi){

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    TPZIntelGen<TSHAPE>::ComputeRequiredData(data,qsi);
    data.fNeedsSol = needsol;

#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        data.fDeformedDirections.Print("Normal Vectors " , sout,EMathematicaInput);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif


}//void

/** Initialize a material data and its attributes based on element dimension, number
 * of state variables and material definitions
 */

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	TPZIntelGen<TSHAPE>::InitMaterialData(data);
}


// Save the element data to a stream
template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::Write(TPZStream &buf, int withclassid) const
{
    this->Write(buf,withclassid);
}


#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzgraphelq2dd.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt2dmapped.h"

using namespace pztopology;

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"
#include "tpzcube.h"
#include "tpztetrahedron.h"
#include "tpzprism.h"

#include "pzelchdivbound2.h"

using namespace pzgeom;
using namespace pzshape;

template<>
void TPZCompElKernelHDiv<pzshape::TPZShapePoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == 0) std::cout << "A point element has no graphical representation\n";
}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == TSHAPE::Dimension /* && Material()->Id() > 0 */) {
		new typename TSHAPE::GraphElType(this,&grafgrid);
	}
}


// template<class TSHAPE>
// void TPZCompElKernelHDiv<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
// 	// if(dimension == TSHAPE::Dimension && this->Material()->Id() > 0) {
// 	// 	new typename TSHAPE::GraphElType(this,&grafgrid);
// 	// }
//     this->CreateGraphicalElement(grafgrid,dimension);
// }



template class TPZRestoreClass< TPZCompElKernelHDiv<TPZShapeLinear>>;
template class TPZRestoreClass< TPZCompElKernelHDiv<TPZShapeTriang>>;
template class TPZRestoreClass< TPZCompElKernelHDiv<TPZShapeQuad>>;
template class TPZRestoreClass< TPZCompElKernelHDiv<TPZShapeCube>>;
template class TPZRestoreClass< TPZCompElKernelHDiv<TPZShapeTetra>>;
template class TPZRestoreClass< TPZCompElKernelHDiv<TPZShapePrism>>;
template class TPZRestoreClass< TPZCompElKernelHDiv<TPZShapePiram>>;

template class TPZCompElKernelHDiv<TPZShapeLinear>;
template class TPZCompElKernelHDiv<TPZShapeTriang>;
template class TPZCompElKernelHDiv<TPZShapeQuad>;
template class TPZCompElKernelHDiv<TPZShapeTetra>;
template class TPZCompElKernelHDiv<TPZShapePrism>;
template class TPZCompElKernelHDiv<TPZShapePiram>;
template class TPZCompElKernelHDiv<TPZShapeCube>;

TPZCompEl * CreateKernelHDivPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElKernelHDiv< TPZShapePoint>(mesh,gel,index);
}

TPZCompEl * CreateKernelHDivLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    return new TPZCompElKernelHDiv< TPZShapeLinear>(mesh,gel,index);
}

TPZCompEl * CreateKernelHDivQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	return new TPZCompElKernelHDiv< TPZShapeQuad>(mesh,gel,index);
}

TPZCompEl * CreateKernelHDivTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	return new TPZCompElKernelHDiv< TPZShapeTriang >(mesh,gel,index);
}

TPZCompEl * CreateKernelHDivCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	return new TPZCompElKernelHDiv< TPZShapeCube >(mesh,gel,index);
}

TPZCompEl * CreateKernelHDivPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	return new TPZCompElKernelHDiv< TPZShapePrism>(mesh,gel,index);
}

TPZCompEl * CreateKernelHDivPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	return new TPZCompElKernelHDiv< TPZShapePiram >(mesh,gel,index);
}

TPZCompEl * CreateKernelHDivTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
	return new TPZCompElKernelHDiv< TPZShapeTetra >(mesh,gel,index);
}

