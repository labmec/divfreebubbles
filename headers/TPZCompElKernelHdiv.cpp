/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHdiv.h"
#include "TPZCompElKernelHdivBC.h"
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
TPZRegisterClassId(&TPZCompElKernelHDiv::ClassId), TPZIntelGen<TSHAPE>(mesh,gel,index)  {

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
TPZCompElKernelHDiv<TSHAPE>::~TPZCompElKernelHDiv(){
    this->~TPZIntelGen<TSHAPE>();
}
 

template<class TSHAPE>
template<class TVar>
void TPZCompElKernelHDiv<TSHAPE>::ComputeRequiredDataT(TPZMaterialDataT<TVar> &data,
                                                TPZVec<REAL> &qsi){

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    TPZIntelGen<TSHAPE>::ComputeRequiredData(data,qsi);
    data.fNeedsSol = needsol;


    for (int i = 0; i < data.dphix.Cols(); i++){
        if (data.dphix.Rows()>1){
            data.fDeformedDirections(0,i) =  data.dphix(1,i);
            data.fDeformedDirections(1,i) = -data.dphix(0,i);
        }
    }

    for (int i = 0; i < data.phi.Rows(); i++){
		data.phi(i,0) = 1.;
	}
	// for (int i = 0; i < data.dphix.Rows(); i++)
    //     for (int j = 0; j < data.dphix.Cols(); j++)
    // 	    	data.dphix(i,j) = 1.;

}//void

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	data.fNeedsSol = true;
	TPZIntelGen<TSHAPE>::InitMaterialData(data);

	int nshape = this->NShapeF();
    // int64_t numvec = TSHAPE::Dimension*TSHAPE::NSides;
    data.fMasterDirections.Resize(3, nshape);
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < nshape; j++)
			data.fMasterDirections(i,j) = 1;

    data.divphi.Zero();
    
    data.fVecShapeIndex.Resize(nshape);
    for (int i=0; i<nshape; i++) {
		data.fVecShapeIndex[i] = std::make_pair(i,1);
    }
    data.fDeformedDirections.Resize(3,nshape);

	// std::cout << "Connect Index = " << this->SequenceNumber() << std::endl;
   
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

