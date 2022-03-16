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
        int nshape = 0;
        int nstate = 1;
        this->fConnectIndexes[TSHAPE::NFacets+1+i]= this->Mesh()->AllocateNewConnect(nshape,nstate,pOrder);//Atualizar fConnectIndexes para ficar na ordem desejada.        
        // this->Connect(TSHAPE::NFacets+1+i).IncrementElConnected();
    } 
    //Reorder the connects
    auto prevCon = this->fConnectIndexes;
    for (int i = 0; i < TSHAPE::NFacets; i++)
    {
        this->fConnectIndexes[2*i  ] = prevCon[i];    
        this->fConnectIndexes[2*i+1] = prevCon[i+TSHAPE::NFacets+1];
    }
    this->fConnectIndexes[TSHAPE::NFacets*2] = prevCon[TSHAPE::NFacets];
    
    // std::cout << "ConIndex " << this->fConnectIndexes << std::endl;
    // std::cout << "Connect " << this->ConnectVec() << std::endl;

    // for (int i = 0; i<NConnects(); i++)
    // {
    //     TPZConnect con = this->Connect(i);
    //     std::cout << "C = " << con << std::endl;
    //     std::cout << "NShape = " << con.NShape() << std::endl;
    // }
    // std::cout << "\n\n\n";
    
    
    //Updates the number of shape functions and also the integration rule
    for (int icon = 0; icon < this->NConnects(); icon++)
    {
        TPZConnect &c = this->Connect(icon);
        int nShapeF = NConnectShapeF(icon,c.Order());
        c.SetNShape(nShapeF);
        int64_t seqnum = c.SequenceNumber();
        int nvar = 1;
        TPZMaterial * mat = this->Material();
        if (mat) nvar = mat->NStateVariables();
        this->Mesh()->Block().Set(seqnum, nvar * nShapeF);
    }
    this->Mesh()->ExpandSolution();
    // this->Mesh()->Block().Resequence();
}


template<class TSHAPE>
int TPZCompElHDivSemiHybrid<TSHAPE>::NConnects() const {
	return this->fConnectIndexes.size();
}

template<class TSHAPE>
int TPZCompElHDivSemiHybrid<TSHAPE>::NSideConnects(int side) const{
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
int TPZCompElHDivSemiHybrid<TSHAPE>::NConnectShapeF(int connect, int order)const
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
        return TPZShapeHDiv<TSHAPE>::ComputeNConnectShapeF(connect,order);    
        break;
    case HDivFamily::EHDivConstant:
        // ARRUMAR ISSO AQUI. ESTA PEGANDO OS NUMEROS DE SHAPE FUNCTIONS
        {
            int conCorrect = connect/2;
            int res = connect % 2;
            int nshape = TPZShapeHDivConstant<TSHAPE>::ComputeNConnectShapeF(conCorrect,order);
            if (res == 1) nshape = 0;
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


template<class TSHAPE>
int TPZCompElHDivSemiHybrid<TSHAPE>::SideConnectLocId(int node,int side) const {
    if (TSHAPE::Dimension == 3){
        std::cout << "this will not work for 3D\n"; 
        DebugStop();
    }

    return 2*(side-TSHAPE::NCornerNodes);
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


