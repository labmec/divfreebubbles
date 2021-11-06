/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHDiv3D.h"
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
// #include "tpzline.h"
#include "tpztriangle.h"


#include "pzshtmat.h"


#ifdef PZ_LOG
static TPZLogger logger("pz.strmatrix");
#endif

using namespace std;

template<class TSHAPE>
TPZCompElKernelHDiv3D<TSHAPE>::TPZCompElKernelHDiv3D(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index, int shapetype) :
TPZRegisterClassId(&TPZCompElKernelHDiv3D::ClassId), TPZCompElHCurlNoGrads<TSHAPE>(mesh,gel,index), fSideOrient(TSHAPE::NFacets,1) {
    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    for(int side = firstside ; side < TSHAPE::NSides-1; side++ )
    {
        fSideOrient[side-firstside] = this->Reference()->NormalOrientation(side);
    }
}

template<class TSHAPE>
TPZCompElKernelHDiv3D<TSHAPE>::~TPZCompElKernelHDiv3D(){
    this->~TPZCompElHCurlNoGrads<TSHAPE>();
}
 

template<class TSHAPE>
template<class TVar>
void TPZCompElKernelHDiv3D<TSHAPE>::ComputeRequiredDataT(TPZMaterialDataT<TVar> &data,
                                                TPZVec<REAL> &qsi){

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = true;
    TPZCompElHCurlNoGrads<TSHAPE>::ComputeRequiredData(data,qsi);
    data.fNeedsSol = needsol;

    int nshape = this->NShapeF();
    data.fDeformedDirections.Resize(3,nshape);

    data.fDeformedDirections=data.curlphi;

#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        //	this->Print(sout);
        // sout << "\nVecshape = " << data.fVecShapeIndex << std::endl;
        sout << "Phi = " << data.curlphi << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
        
    }
#endif

    data.phi.Resize(nshape,3);
    REAL val = 1.;

    for (int i = 0; i < data.phi.Rows(); i++){
		data.phi(i,0) = val;
        data.phi(i,1) = val;
        data.phi(i,2) = val;
	}

	for (int i = 0; i < data.dphix.Rows(); i++)
        for (int j = 0; j < data.dphix.Cols(); j++)
    	    	data.dphix(i,j) = 1.;
    
    // for (int i=0; i<data.fVecShapeIndex.size(); i++) {
	// 	data.fVecShapeIndex[i] = std::make_pair(i,1);
    // }
    if (data.fNeedsSol) {
        this->ReallyComputeSolution(data);
    }

}//void

template<class TSHAPE>
void TPZCompElKernelHDiv3D<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	data.fNeedsSol = true;
	TPZCompElHCurlNoGrads<TSHAPE>::InitMaterialData(data);

    data.fShapeType = data.EVecandShape;
    
}

/**
 * @brief It returns the normal orientation of the reference element by the side.
 * Only side that has dimension larger than zero and smaller than me.
 * @param side: side of the reference elemen
 */
template<class TSHAPE>
int TPZCompElKernelHDiv3D<TSHAPE>::GetSideOrient(int side){

    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    if (side < firstside || side >= TSHAPE::NSides - 1) {
        DebugStop();
    }
    return fSideOrient[side-firstside];
}

/**
 * @brief It set the normal orientation of the element by the side.
 * Only side that has dimension equal to my dimension minus one.
 * @param side: side of the reference elemen
 */
template<class TSHAPE>
void TPZCompElKernelHDiv3D<TSHAPE>::SetSideOrient(int side, int sideorient){

    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    if (side < firstside || side >= TSHAPE::NSides - 1) {
        DebugStop();
    }
    fSideOrient[side-firstside] = sideorient;
}

template<class TSHAPE>
template<class TVar>
void TPZCompElKernelHDiv3D<TSHAPE>::ComputeSolutionKernelHdivT(TPZMaterialDataT<TVar> &data)
{
    TPZCompElHCurlNoGrads<TSHAPE>::ReallyComputeSolution(data);
    // data.fDeformedDirections=data.curlphi;

    const int dim = 3; // Hdiv vectors are always in R3
    const int nstate = this->Material()->NStateVariables();
    const int ncon = this->NConnects();

    TPZFMatrix<TVar> &MeshSol = this->Mesh()->Solution();

    int64_t numbersol = MeshSol.Cols();

    if(numbersol != 1)
    {
        DebugStop();
    }
    data.sol.Resize(numbersol);
    data.dsol.Resize(numbersol);
    data.divsol.Resize(numbersol);

    for (int64_t is=0; is<numbersol; is++)
    {
        data.sol[is].Resize(dim*nstate);
        data.sol[is].Fill(0);
        data.dsol[is].Redim(dim*nstate, dim);
        data.divsol[is].Resize(nstate);
        data.divsol[is].Fill(0.);
    }
    data.sol = data.curlsol;

}


#include "pzshapecube.h"
#include "pzshapetetra.h"

using namespace pztopology;

#include "tpzcube.h"
#include "tpztetrahedron.h"

using namespace pzgeom;
using namespace pzshape;


template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeCube>>;
template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeTetra>>;

template class TPZCompElKernelHDiv3D<TPZShapeTetra>;
template class TPZCompElKernelHDiv3D<TPZShapeCube>;
