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
static TPZLogger logger("pz.mesh.TPZCompElKernelHDiv3D");
#endif

using namespace std;


template<class TSHAPE>
TPZCompElKernelHDiv3D<TSHAPE>::TPZCompElKernelHDiv3D(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index) :
TPZRegisterClassId(&TPZCompElKernelHDiv3D::ClassId), TPZCompElHCurlNoGrads<TSHAPE>(mesh,gel,index), fSideOrient(TSHAPE::NFacets,1) {

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
    data.fDeformedDirections.Zero();

    data.fDeformedDirections=data.curlphi;
    // GetCurl(data);
    // // TPZFNMatrix<220,REAL> dphix(3,data.dphix.Cols());
    // // TPZFMatrix<REAL> &dphi = data.dphix;;
    // // TPZAxesTools<REAL>::Axes2XYZ(dphi, dphix, data.axes);

    
    // for (int i = 0; i < data.curlphi.Cols(); i++)
    // {
    //     for (int j = 0; j < data.dphix.Cols(); j++)
    //     {
    //         data.fDeformedDirections(0,i) += data.curlphi(2,i)*data.dphix(1,j)-data.curlphi(1,i)*data.dphix(2,j);
    //         data.fDeformedDirections(1,i) += data.curlphi(0,i)*data.dphix(2,j)-data.curlphi(2,i)*data.dphix(0,j);
    //         data.fDeformedDirections(2,i) += data.curlphi(1,i)*data.dphix(0,j)-data.curlphi(0,i)*data.dphix(1,j);
    //     }
    // }


    for (int i = 0; i < data.phi.Rows(); i++){
		data.phi(i,0) = 1.;
	}
	for (int i = 0; i < data.dphix.Rows(); i++)
        for (int j = 0; j < data.dphix.Cols(); j++)
    	    	data.dphix(i,j) = 1.;
    
    // for (int i=0; i<data.fVecShapeIndex.size(); i++) {
	// 	data.fVecShapeIndex[i] = std::make_pair(i,1);
    // }

}//void

template<class TSHAPE>
void TPZCompElKernelHDiv3D<TSHAPE>::InitMaterialData(TPZMaterialData &data)
{
	data.fNeedsSol = true;
	TPZCompElHCurlNoGrads<TSHAPE>::InitMaterialData(data);

	// int nshape = this->NShapeF();    
    // int64_t size = nshape*3;//(TSHAPE::Dimension);
    // data.fVecShapeIndex.Resize(size);
    // // auto size = data.fVecShapeIndex.size();
    
    // for (int i=0; i<size; i++) {
	// 	data.fVecShapeIndex[i] = std::make_pair(1,i);
    // }
    
   
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
void TPZCompElKernelHDiv3D<TSHAPE>:: Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol)
{
    
    TPZMaterialDataT<STATE> data;
    constexpr bool hasPhi{false};
    this->ComputeSolution(qsi,data,hasPhi);

    sol.Resize(3);
    
    // REAL Sol = data.sol[0];
    // data.sol.resize(3);
    // data.sol[0] = Sol;
    sol[0] = data.sol[0][0];
    data.sol[0].Resize(3);
    // sol = std::move(data.sol[0]);
}


template<class TSHAPE>
template<class TVar>
void TPZCompElKernelHDiv3D<TSHAPE>::ComputeSolutionKernelHdivT(TPZMaterialDataT<TVar> &data)
{
    
    const int dim = 3; 
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
    TPZFNMatrix<220,REAL> dphix(3,data.dphix.Cols());
    TPZFMatrix<REAL> &dphi = data.dphix;;

    TPZAxesTools<REAL>::Axes2XYZ(dphi, dphix, data.axes);

    TPZBlock &block =this->Mesh()->Block();
    int ishape=0,ivec=0,counter=0;



    // TPZFMatrix<TVar> dsolX(3,1);
    // // const auto &sol = data.sol[0];
    // const auto &dsol = data.dsol[0];
    // TPZAxesTools<TVar>::Axes2XYZ(dsol,dsolX,data.axes);

//     int nshapeV = data.fVecShapeIndex.NElements();

    for(int in=0; in<ncon; in++)
    {
        TPZConnect *df = &this->Connect(in);
        int64_t dfseq = df->SequenceNumber();
        int dfvar = block.Size(dfseq);
//         // pos : position of the block in the solution matrix
        int64_t pos = block.Position(dfseq);

//         /// ish loops of the number of shape functions associated with the block
        for(int ish=0; ish<dfvar/nstate; ish++)
        {
            ishape  = data.fVecShapeIndex[counter].first;
            // std::cout << "ISHAPE = " << ishape << std::endl;
            for(int idf=0; idf<nstate; idf++)
            {
                TVar meshsol = MeshSol(pos+ish*nstate+idf,0);
                // REAL phival = data.phi(ishape,0);
                //Computes sol and dsol
                // data.sol[0][dim*idf] += phival*meshsol;
                // data.dsol[0](dim*idf,0)+= meshsol * dphix(0,ishape);
                // data.dsol[0](dim*idf,1)+= meshsol * dphix(1,ishape);

                //Compute rotated flux
                data.sol[0][dim*idf+0] += meshsol;
                data.sol[0][dim*idf+1] += meshsol;
                data.sol[0][dim*idf+2] += meshsol;
            }
            counter++;
        }
    }
    // data.sol[1][0] = 0.;
    // data.sol[0][0] = -data.dsol[0](0,1);
    // data.sol[0][1] = data.dsol[0](0,0);
    // data.sol[0][2] = 0.;
}

template<class TSHAPE>
template<class TVar>
void TPZCompElKernelHDiv3D<TSHAPE>::GetCurl(TPZMaterialDataT<TVar> &data)
{
    data.fDeformedDirections=data.curlphi;

    // constexpr auto dim = 3;
    // const auto nShapeFuncs = data.fVecShapeIndex.size();
    
    // const REAL jacInv = 1/data.detjac;
    // TPZFNMatrix<dim,REAL> tempCurl(dim, 1, 0),gradPhiCrossDirections(dim, 1, 0);
    
    // for(auto iShapeFunc = 0; iShapeFunc < nShapeFuncs; iShapeFunc++) {
    //     const auto iVec = data.fVecShapeIndex[iShapeFunc].first;
    //     const auto iShape = data.fVecShapeIndex[iShapeFunc].second;
        
    //     for(auto ix = 0; ix < dim; ix++) {
    //         const auto i = (ix+1)%dim;
    //         const auto j = (ix+2)%dim;
    //         gradPhiCrossDirections(ix,0) =
    //             data.curlphi.GetVal(i,iShape) * data.fMasterDirections.GetVal(j,iVec)-
    //             data.curlphi.GetVal(j,iShape) * data.fMasterDirections.GetVal(i,iVec);
    //     }
        

    //     tempCurl = data.jacobian * gradPhiCrossDirections;
    //     tempCurl *= jacInv;
    //     for (auto ix = 0; ix < dim; ix++) {
    //         data.fDeformedDirections.PutVal(ix, iShapeFunc,tempCurl.GetVal(ix,0));
    //     }
    // }

    double tol = 1.e-10;
    for (int i = 0; i < data.curlphi.Cols(); i++)
    {
        if ((fabs(data.curlphi(0,i)) < tol) && (fabs(data.curlphi(1,i)) < tol) && (fabs(data.curlphi(2,i)) < tol)){
            std::cout << "PROBLEM WITH CURL = 0 !!" << std::endl; 
            std::cout << "Curl = \n " << data.curlphi(0,i) << " " << data.curlphi(1,i) << " " << data.curlphi(2,i) << std::endl;
        }
        if ((fabs(data.fDeformedDirections(0,i)) < tol) && (fabs(data.fDeformedDirections(1,i)) < tol) && (fabs(data.fDeformedDirections(2,i)) < tol)){
            std::cout << "PROBLEM WITH fDeformedDirections = 0 !!" << std::endl; 
            std::cout << "Curl = \n " << data.fDeformedDirections(0,i) << " " << data.fDeformedDirections(1,i) << " " << data.fDeformedDirections(2,i) << std::endl;
        }
    }
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
// #include "pzrefpoint.h"
// #include "pzgeopoint.h"
// #include "pzshapepoint.h"
#include "pzgraphelq2dd.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt2dmapped.h"

using namespace pztopology;

// #include "tpzpoint.h"
// #include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"
#include "tpzcube.h"
#include "tpztetrahedron.h"
#include "tpzprism.h"

#include "pzelchdivbound2.h"

using namespace pzgeom;
using namespace pzshape;


// template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeLinear>>;
// template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeTriang>>;
// template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeQuad>>;
// template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeCube>>;
template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapeTetra>>;
// template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapePrism>>;
// template class TPZRestoreClass< TPZCompElKernelHDiv3D<TPZShapePiram>>;

// template class TPZCompElKernelHDiv3D<TPZShapeLinear>;
// template class TPZCompElKernelHDiv3D<TPZShapeTriang>;
// template class TPZCompElKernelHDiv3D<TPZShapeQuad>;
template class TPZCompElKernelHDiv3D<TPZShapeTetra>;
// template class TPZCompElKernelHDiv3D<TPZShapePrism>;
// template class TPZCompElKernelHDiv3D<TPZShapePiram>;
// template class TPZCompElKernelHDiv3D<TPZShapeCube>;
