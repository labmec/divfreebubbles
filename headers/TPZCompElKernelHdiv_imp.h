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
#include "pzhdivpressure.h"
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
TPZRegisterClassId(&TPZCompElKernelHDiv::ClassId),
TPZIntelGen<TSHAPE>(mesh,gel,index,1), fSideOrient(TSHAPE::NFacets,1) {
	this->TPZInterpolationSpace::fPreferredOrder = mesh.GetDefaultOrder();
	int nconflux= TPZCompElKernelHDiv::NConnects();
    this->fConnectIndexes.Resize(nconflux);
	gel->SetReference(this);

//    int nfaces = TSHAPE::NumSides(TSHAPE::Dimension-1);
    TPZStack<int> facesides;
    TSHAPE::LowerDimensionSides(TSHAPE::NSides-1,facesides,TSHAPE::Dimension-1);
    facesides.Push(TSHAPE::NSides-1);
	for(int i=0;i< facesides.size(); i++)
	{
        int sideaux= facesides[i];
		this->fConnectIndexes[i] = this->CreateMidSideConnect(sideaux);
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "After creating last flux connect " << i << std::endl;
            //	this->Print(sout);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif

		mesh.ConnectVec()[this->fConnectIndexes[i]].IncrementElConnected();
		this->IdentifySideOrder(sideaux);
    }


    int sideorder = EffectiveSideOrder(TSHAPE::NSides-1);
//    if(TSHAPE::Type()==EQuadrilateral)
//    {
//        sideorder++;
//    }

    sideorder++;

	sideorder = 2*sideorder;
	if (sideorder > this->fIntRule.GetMaxOrder()) sideorder = this->fIntRule.GetMaxOrder();
	TPZManVector<int,3> order(3,sideorder);
	this->fIntRule.SetOrder(order);
    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    for(int side = firstside ; side < TSHAPE::NSides-1; side++ )
    {
        fSideOrient[side-firstside] = this->Reference()->NormalOrientation(side);
    }
    auto *mat =
        dynamic_cast<TPZMatSingleSpace *>(this->Material());
    if (mat)
    {
        int order = mat->IntegrationRuleOrder(MaxOrder());
        TPZManVector<int,3> ord(gel->Dimension(),order);
        this->fIntRule.SetOrder(ord);
    }


}

template<class TSHAPE>
TPZCompElKernelHDiv<TSHAPE>::TPZCompElKernelHDiv(TPZCompMesh &mesh, const TPZCompElKernelHDiv<TSHAPE> &copy) :
TPZRegisterClassId(&TPZCompElKernelHDiv::ClassId),
TPZIntelGen<TSHAPE>(mesh,copy), fSideOrient(copy.fSideOrient)
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
TPZIntelGen<TSHAPE>(mesh,copy,gl2lcConMap,gl2lcElMap), fSideOrient(copy.fSideOrient)
{
	this-> fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<NConnects();i++)
	{
		int lcIdx = -1;
		int glIdx = copy.fConnectIndexes[i];
		if (gl2lcConMap.find(glIdx) != gl2lcConMap.end()) lcIdx = gl2lcConMap[glIdx];
		else
		{
			std::stringstream sout;
			sout << "ERROR in : " << __PRETTY_FUNCTION__
			<< " trying to clone the connect index: " << glIdx
			<< " wich is not in mapped connect indexes!";
			LOGPZ_ERROR(logger, sout.str().c_str());
			this-> fConnectIndexes[i] = -1;
			return;
		}
		this-> fConnectIndexes[i] = lcIdx;
	}
}

template<class TSHAPE>
TPZCompElKernelHDiv<TSHAPE>::TPZCompElKernelHDiv() :
TPZRegisterClassId(&TPZCompElKernelHDiv::ClassId),
TPZIntelGen<TSHAPE>()
{
	this->fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NSides;i++) {
		this-> fConnectIndexes[i] = -1;
	}

}

template<class TSHAPE>
TPZCompElKernelHDiv<TSHAPE>::~TPZCompElKernelHDiv(){
    TPZGeoEl *gel = this->Reference();
    if (gel && gel->Reference() != this) {
        return;
    }
    for (int side=TSHAPE::NCornerNodes; side < TSHAPE::NSides; side++) {
        if (TSHAPE::SideDimension(side) != TSHAPE::Dimension-1) {
            continue;
        }
        TPZGeoElSide gelside(this->Reference(),side);
        TPZStack<TPZCompElSide> celstack;
        TPZCompElSide largecel = gelside.LowerLevelCompElementList2(0);
        if (largecel) {
            int cindex = SideConnectLocId(0, side);
            TPZConnect &c = this->Connect(cindex);
            c.RemoveDepend();
        }
        if (gelside.Element()){
            gelside.HigherLevelCompElementList3(celstack, 0, 1);
        }
        int64_t ncel = celstack.size();
        for (int64_t el=0; el<ncel; el++) {
            TPZCompElSide celside = celstack[el];
            TPZCompEl *celsmall = celside.Element();
            TPZGeoEl *gelsmall = celsmall->Reference();
            if (gelsmall->SideDimension(celside.Side()) != gel->Dimension()-1) {
                continue;
            }
            TPZInterpolatedElement *intelsmall = dynamic_cast<TPZInterpolatedElement *>(celsmall);
            if (!intelsmall) {
                DebugStop();
            }
            int cindex = intelsmall->SideConnectLocId(0, celside.Side());
            TPZConnect &c = intelsmall->Connect(cindex);
            c.RemoveDepend();
        }
    }
    if (gel){
        gel->ResetReference();
    }
}

template<class TSHAPE>
MElementType TPZCompElKernelHDiv<TSHAPE>::Type() {
	return TSHAPE::Type();
}


template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::NConnects() const {
	return TSHAPE::NFacets + 1;
}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::SetConnectIndex(int i, int64_t connectindex){
#ifndef NODEBUG
	if(i<0 || i>= this->NConnects()) {
		std::cout << " TPZCompElKernelHDiv<TSHAPE>::SetConnectIndex index " << i <<
		" out of range\n";
		DebugStop();
		return;
	}
#endif
	this-> fConnectIndexes[i] = connectindex;
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << endl<<"Setting Connect : " << i << " to connectindex " << connectindex<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::NConnectShapeF(int connect, int order)const
{
#ifdef DEBUG
    if (connect < 0 || connect > TSHAPE::NFacets) {
        DebugStop();
    }
#endif
    MElementType thistype = TSHAPE::Type();
    if(thistype == EOned)
    {
        if(connect < 2) return 1;
        else return order;
    }
    else if(thistype == ETriangle)
    {
        if(connect < TSHAPE::NFacets) return (order+1);
        else return (order+1)*(order+1)-1;
    }
    else if(thistype == EQuadrilateral)
    {
        if(connect < TSHAPE::NFacets) return (order+1);
        else return 2*order*(order+1);
    }
    else if(thistype == ETetraedro)
    {
        if(connect < TSHAPE::NFacets) return (order+1)*(order+2)/2;
        else return order*(order+2)*(order+3)/2;
    }
    else if(thistype == EPrisma)
    {
        if(connect == 0 || connect == 4) return (order+1)*(order+2)/2;
        else if(connect < TSHAPE::NFacets) return (order+1)*(order+1);
        else return order*order*(3*order+5)/2+7*order-2;
    }
    else if(thistype == ECube)
    {
        if(connect < TSHAPE::NFacets) return (order+1)*(order+1);
        else return 3*order*(order+1)*(order+1);
    }
    DebugStop();
    // @TODO put in the analytic values
    if (connect < TSHAPE::NFacets) {
         int64_t connectindex = ConnectIndex(connect);
//         int order = 0;
//         if (connectindex >= 0) {
//             order = this->Connect(connect).Order();
//         }
//         else
//         {
//             order = this->fPreferredOrder;
//         }
         const int nfaces = TSHAPE::NFacets;
         int face = TSHAPE::NSides-nfaces+connect-1;
         TPZStack<int> lowerdimensionsides;
         TSHAPE::LowerDimensionSides(face,lowerdimensionsides);
         int nshape = TSHAPE::NConnectShapeF(face,order);
         for (int is=0; is<lowerdimensionsides.size(); is++) {
             nshape += TSHAPE::NConnectShapeF(lowerdimensionsides[is],order);
         }
         return nshape;
     }

    // connect == NFacets

     const int nfaces = TSHAPE::NFacets;
     int face = TSHAPE::NSides-nfaces-1;
    // nvecignore is the first vector index associated with the internal vectors
     int nvecignore = 0;
     for (int side = face; side < TSHAPE::NSides-1; side++) {
         nvecignore += TSHAPE::NContainedSides(side);
     }


     TPZManVector<int,TSHAPE::Dimension*TSHAPE::NSides+1> vecside(TSHAPE::Dimension*TSHAPE::NSides),bilinear(TSHAPE::Dimension*TSHAPE::NSides),directions(TSHAPE::Dimension*TSHAPE::NSides);
     TSHAPE::GetSideHDivDirections(vecside,directions,bilinear);
//     if (TSHAPE::Type()==ETriangle||TSHAPE::Type()==ETetraedro) {
////         pressureorder=this->fPreferredOrder-1;
//     }
//     else if (TSHAPE::Type()==ECube||TSHAPE::Type()==EQuadrilateral||TSHAPE::Type()==EPrisma) {
//         pressureorder=this->fPreferredOrder;
//     }
//     else
//     {
//         // Tipo nao implementado
//         DebugStop();
//     }
     TPZManVector<int,27> orders(TSHAPE::NSides-TSHAPE::NCornerNodes,0);
     FillOrder(orders);
     int nshape = TSHAPE::NShapeF(orders);

     TPZManVector<int64_t, TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes);
     for (int i=0; i<id.size(); i++) {
         id[i] = this->Reference()->NodePtr(i)->Id();
     }


     int nexternalvectors = 0;
     nexternalvectors = nvecignore;

     TPZGenMatrix<int> shapeorders(nshape,3);
     TSHAPE::ShapeOrder(id, orders, shapeorders);
    {
        static int first = 0;
        if (first==0) {
            shapeorders.Print("ShapeOrders");
            first++;
        }
    }
     // VectorSide indicates the side associated with each vector entry
     TPZManVector<int64_t,30> FirstIndex(TSHAPE::NSides+1);
     // the first index of the shape functions
     FirstShapeIndex(FirstIndex);
     //FirstIndex.Print();

#ifdef PZ_LOG
     if (logger.isDebugEnabled())
     {
         std::stringstream sout;
         sout << "FirstIndex "<<FirstIndex << std::endl;
         LOGPZ_DEBUG(logger,sout.str())
     }
#endif

     int count = 0;

    int internalorder = this->Connect(connect).Order();

     int64_t nvec = vecside.NElements();
     for (int locvec = nexternalvectors; locvec<nvec; locvec++)
     {
         int side = vecside[locvec];
         int bil = bilinear[locvec];
         int dir = directions[locvec];

         MElementType tipo = TSHAPE::Type(side);

         int firstshape = FirstIndex[side];
         int lastshape = FirstIndex[side+1];

         for (int ish = firstshape; ish<lastshape; ish++)
         {
             int sidedimension = TSHAPE::SideDimension(side);
             int maxorder[3] = {internalorder,internalorder,internalorder};
             if (bil) {
                 maxorder[dir]++;
             }
             int include=true;
             for (int d=0; d<sidedimension; d++)
             {
                 if (tipo==ETriangle||tipo==ETetraedro)//
                 {
                     if (shapeorders(ish,d) > maxorder[d]+1) {
                         include = false;
                     }
                 }
                 else if(tipo==EQuadrilateral)
                 {
                     if (shapeorders(ish,d) > maxorder[d]) {
                         include = false;
                     }
                 }
                 else if(tipo==ECube)
                 {
                     if (shapeorders(ish,d) > maxorder[d]) {
                         include = false;
                     }
                 }
                 else if (tipo==EPrisma)
                 {
                     //DebugStop();
                     if (shapeorders(ish,d) > maxorder[d]) {
                         include = false;
                     }
                 }
                 else if (tipo == EOned)
                 {
                     if (shapeorders(ish,0) > maxorder[d]) {
                         include = false;
                     }
                 }
                 else if (tipo == EPoint)
                 {
                     DebugStop();
                 }
                 else if (tipo == EPiramide)
                 {
                     if (shapeorders(ish,d) > maxorder[d]) {
                         include = false;
                     }
                 }
                 else
                 {
                     DebugStop();
                 }
             }
             if (include)
             {
                 count++;
             }
         }
     }
    return count;


 #ifdef PZ_LOG
     {
         std::stringstream sout;
         sout <<__PRETTY_FUNCTION__<< "unhandled case ";
         LOGPZ_ERROR(logger,sout.str())
     }
 #endif
 return -1;

 }




////
template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::SetIntegrationRule(int ord) {
	TPZManVector<int,3> order(TSHAPE::Dimension,ord);
	this->fIntRule.SetOrder(order);
}


template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::NSideConnects(int side) const{
	if(TSHAPE::SideDimension(side)<= Dimension()-2) return 0;
	if(TSHAPE::SideDimension(side)==Dimension()-1) return 1;
	if(TSHAPE::SideDimension(side)== Dimension()) {
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
int TPZCompElKernelHDiv<TSHAPE>::SideConnectLocId(int node,int side) const {
#ifdef PZDEBUG
	if(TSHAPE::SideDimension(side)<= TSHAPE::Dimension - 2 || node >= NSideConnects(side)) {
		PZError << "TPZCompElKernelHDiv<TSHAPE>::SideConnectLocId no connect associate " <<  endl;
		return -1;
	}
#endif

    return node+side-(TSHAPE::NSides-TSHAPE::NumSides(TSHAPE::Dimension-1)-1);
}

template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::ConnectSideLocId(int connect) const{

    int side = connect+TSHAPE::NSides-TSHAPE::NumSides(TSHAPE::Dimension-1)-1 ;
    return side;
}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
	ord.Resize(NConnects());
	int i;
	for(i=0; i<NConnects(); i++) {
		ord[i] = ConnectOrder(i);
	}
}


template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::PreferredSideOrder(int side) {
	if(TSHAPE::SideDimension(side) < Dimension()-1)
	{
		PZError << __PRETTY_FUNCTION__ << " side " << side << std::endl;
	}
	int connect= SideConnectLocId(0,side);
	if(connect<0 || connect > NConnects()) {
		PZError << "TPZCompElKernelHDiv<TSHAPE>::PreferredSideOrder no polynomial associate " <<  endl;
		return -1;
	}
	if(connect<NConnects()) {
			int order =this->fPreferredOrder;
			return order;//this->AdjustPreferredSideOrder(side,order);
	}
	PZError << "TPZCompElKernelHDiv<TSHAPE>::PreferredSideOrder called for connect = " << connect << "\n";
	return 0;

}

template<class TSHAPE>
int64_t TPZCompElKernelHDiv<TSHAPE>::ConnectIndex(int con) const{
#ifndef NODEBUG
	if(con<0 || con > TSHAPE::NFacets) {
		std::cout << "TPZCompElKernelHDiv::ConnectIndex wrong parameter connect " << con <<
		" NConnects " << TSHAPE::NFacets << std::endl;
		DebugStop();
		return -1;
	}

#endif

//    #ifndef NODEBUG
//    	if(con<0) {
//    		std::cout << "TPZCompElHDiv::ConnectIndex wrong parameter connect " << con <<
//    		" NConnects " << this-> NConnects() << std::endl;
//    		DebugStop();
//    		return -1;
//    	}
//    	
//    #endif
//    
//   
//    if(con>= this->NConnects())
//    {
//        int con2= con-TSHAPE::NCornerNodes;
//        return this->fConnectIndexes[con2];
//    }

	return this->fConnectIndexes[con];
}


template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::SetPreferredOrder(int order)
{
		TPZIntelGen<TSHAPE>:: SetPreferredOrder(order);
	//this->fPreferredOrder = order;
}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::SetSideOrder(int side, int order){
	int connectaux= SideConnectLocId(0,side);
	if(connectaux<0 || connectaux > this-> NConnects()) {
		PZError << "TPZCompElKernelHDiv::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
#ifdef PZ_LOG
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
		LOGPZ_DEBUG(logger,sout.str())
#endif
		return;
	}
	TPZConnect &c = this->Connect(connectaux);
    c.SetOrder(order,this->fConnectIndexes[connectaux]);
    int64_t seqnum = c.SequenceNumber();
    int nvar = 1;
    TPZMaterial * mat =this-> Material();
    if(mat) nvar = mat->NStateVariables();
    c.SetNState(nvar);
    int nshape =this->NConnectShapeF(connectaux,order);
    c.SetNShape(nshape);
	this-> Mesh()->Block().Set(seqnum,nshape*nvar);
}


template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::ConnectOrder(int connect) const{
	if (connect < 0 || connect >= this->NConnects()){
#ifdef PZ_LOG
		{
			std::stringstream sout;
			sout << "Connect index out of range connect " << connect <<
			" nconnects " << NConnects();
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		return -1;
	}

	if (this->fConnectIndexes[connect] == -1) {
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " connect " << connect
		<< " is not initialized" << std::endl;
#ifdef PZ_LOG
		LOGPZ_ERROR(logger,sout.str());
#else
		std::cout << sout.str() << std::endl;
#endif
		return 0;
	}

    TPZConnect &c = this-> Connect(connect);
    return c.Order();
}

template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::EffectiveSideOrder(int side) const
{
	if(!NSideConnects(side)) return -1;
	int cindex =SideConnectLocId(0, side);
#ifdef PZDEBUG
    if(cindex<0 || cindex >= NConnects()) DebugStop();
#endif
    return ConnectOrder(cindex);

}

/**return the first shape associate to each side*/
template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::FirstShapeIndex(TPZVec<int64_t> &Index) const {

    TPZManVector<int> orders(TSHAPE::NSides-TSHAPE::NCornerNodes);
    FillOrder(orders);
    Index[0] = 0;
	for(int iside=0;iside<TSHAPE::NSides;iside++)
	{
        int sideorder = 1;
        if (iside >= TSHAPE::NCornerNodes) {
            sideorder = orders[iside-TSHAPE::NCornerNodes];
        }
        int temp = Index[iside] + TSHAPE::NConnectShapeF(iside,sideorder);
        Index[iside+1] = temp;
	}

#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "First  Index " << Index;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
}


template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::NFluxShapeF() const{
    int in,result=0;
    int nn=TPZCompElKernelHDiv::NConnects();
    for(in=0;in<nn;in++){
//#ifdef PZ_LOG
//				std::stringstream sout;
//				sout << "conect " << in<< " seq number "<<seqnum<<" num func "<<TPZCompElHDiv::NConnectShapeF(in);
//				LOGPZ_DEBUG(logger,sout.str())
//#endif
        int order = this->Connect(in).Order();
        result += TPZCompElKernelHDiv::NConnectShapeF(in,order);
    }


//#ifdef PZ_LOG
//    std::stringstream sout;
//    sout << "Num funcoes associada ao fluxo " << result;
//    LOGPZ_DEBUG(logger,sout.str())
//#endif
    return result;


}


/**
 * @brief It returns the normal orientation of the reference element by the side.
 * Only side that has dimension larger than zero and smaller than me.
 * @param side: side of the reference elemen
 */
template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::GetSideOrient(int side){

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
void TPZCompElKernelHDiv<TSHAPE>::SetSideOrient(int side, int sideorient){

    int firstside = TSHAPE::NSides-TSHAPE::NFacets-1;
    if (side < firstside || side >= TSHAPE::NSides - 1) {
        DebugStop();
    }
    fSideOrient[side-firstside] = sideorient;
}

//compute the values of the shape function of the side
// this method is used by the method RestrainSide. It does not consider the side orientation
template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {

    if(side==TSHAPE::NSides || point.size() != TSHAPE::Dimension-1){
        std::cout<<"Don't have side shape associated to this side";
        DebugStop();
    }
	if(TSHAPE::SideDimension(side)!= TSHAPE::Dimension -1 ){
		return ;
	}
    int ncontained = TSHAPE::NContainedSides(side);
    int nsideshape = 0;
    int connectlocid = SideConnectLocId(0, side);
    int order = this->Connect(connectlocid).Order();
    int is;
    for (is=0; is<ncontained; is++) {
        int ic = TSHAPE::ContainedSideLocId(side,is);
        nsideshape += TSHAPE::NConnectShapeF(ic,order);
    }
#ifdef PZDEBUG
    if (nsideshape != this->NSideShapeF(side)) {
        DebugStop();
    }
#endif

    TPZGeoEl *gel = this->Reference();
    //int nc = gel->NCornerNodes();
    int nsn = TSHAPE::NSideNodes(side);
    TPZManVector<int64_t,8> id(nsn);
    for (int ic=0; ic<nsn; ic++) {
        int locid = TSHAPE::SideNodeLocId(side,ic);
        id[ic] = gel->Node(locid).Id();
    }

    //int idsize = id.size();
    TPZManVector<int,9> permutegather(ncontained);
    int transformid;


    MElementType sidetype = TSHAPE::Type(side);
    switch (sidetype) {
        case EOned:
            transformid = pztopology::TPZLine::GetTransformId(id);
            pztopology::TPZLine::GetSideHDivPermutation(transformid, permutegather);
            break;
        case EQuadrilateral:
            transformid = pztopology::TPZQuadrilateral::GetTransformId(id);
            pztopology::TPZQuadrilateral::GetSideHDivPermutation(transformid, permutegather);
            break;
        case ETriangle:
            transformid = pztopology::TPZTriangle::GetTransformId(id);
            pztopology::TPZTriangle::GetSideHDivPermutation(transformid, permutegather);
            break;
        default:
            DebugStop();
            break;
    }

    TPZManVector<int,TSHAPE::NSides> ord(TSHAPE::NSides,order);

    int sidedimension = TSHAPE::SideDimension(side);
    TPZFNMatrix<50,REAL> philoc(nsideshape,1),dphiloc(sidedimension,nsideshape);

    TSHAPE::SideShape(side,point,id,ord,philoc,dphiloc);

    int ncs = TSHAPE::NContainedSides(side);
    TPZManVector<int64_t,28> FirstIndex(ncs+1,0);
    for (int ls=0; ls<ncs; ls++) {
        int localside = TSHAPE::ContainedSideLocId(side,ls);
        FirstIndex[ls+1] = FirstIndex[ls]+TSHAPE::NConnectShapeF(localside,order);
    }

    REAL detjac = 1.;
    {
        TPZGeoElSide gelside = TPZGeoElSide(this->Reference(),side);
        int dim = gel->SideDimension(side);
        TPZFNMatrix<9,REAL> jac(dim,dim),jacinv(dim,dim),axes(dim,3);
        gelside.Jacobian(point, jac, axes, detjac, jacinv);
    }

    for (int side=0; side < ncs; side++) {
        int ifirst = FirstIndex[side];
        int kfirst = FirstIndex[permutegather[side]];
        int nshape = FirstIndex[side+1]-FirstIndex[side];
        for (int i=0; i<nshape; i++) {
            phi(ifirst+i,0) = philoc(kfirst+i,0)/detjac;
            for (int d=0; d< sidedimension; d++) {
                dphi(d,ifirst+i) = dphiloc(d,kfirst+i)/detjac;
            }
        }
    }

}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>:: Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol)
{
    //TODOCOMPLEX
    if (var == 99) {
        return TPZIntelGen<TSHAPE>::Solution(qsi,var,sol);
    }
    TPZMaterialDataT<STATE> data;
    constexpr bool hasPhi{false};
    this->ComputeSolution(qsi,data,hasPhi);
    sol = std::move(data.sol[0]);
}

template<class TSHAPE>
template<class TVar>
void TPZCompElKernelHDiv<TSHAPE>::ComputeSolutionHDivT(TPZMaterialDataT<TVar> &data)
{
    
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
    TPZFNMatrix<220,REAL> dphix(3,data.dphix.Cols());
    TPZFMatrix<REAL> &dphi = data.dphix;;

    TPZAxesTools<REAL>::Axes2XYZ(dphi, dphix, data.axes);

    TPZFMatrix<TVar> GradOfPhiHdiv(dim,dim);
    GradOfPhiHdiv.Zero();


    int normvecRows = data.fDeformedDirections.Rows();
    int normvecCols = data.fDeformedDirections.Cols();
    TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);
    TPZManVector<TPZFNMatrix<9,REAL>,18> GradNormalvec(normvecCols);
    for (int i=0; i<GradNormalvec.size(); i++) {
        GradNormalvec[i].Redim(dim,dim);
    }

    if (data.fNeedsDeformedDirectionsFad) {
        for (int e = 0; e < normvecRows; e++) {
            for (int s = 0; s < normvecCols; s++) {
                Normalvec(e,s)=data.fDeformedDirectionsFad(e,s).val();
            }
        }

        TPZFNMatrix<4,REAL> Grad0(3,3,0.);
	TPZGeoEl *ref = this->Reference();
	const int gel_dim = ref->Dimension();

	for (int s = 0; s < normvecCols; s++) {
            for (int i = 0; i < gel_dim; i++) {
                for (int j = 0; j < gel_dim; j++) {
                    Grad0(i,j)=data.fDeformedDirectionsFad(i,s).fastAccessDx(j);
                }
            }
            GradNormalvec[s] = Grad0;
        }

    }else{
        Normalvec=data.fDeformedDirections;
    }

    TPZBlock &block =this->Mesh()->Block();
    int ishape=0,ivec=0,counter=0;

    int nshapeV = data.fVecShapeIndex.NElements();

    for(int in=0; in<ncon; in++)
    {
        TPZConnect *df = &this->Connect(in);
        int64_t dfseq = df->SequenceNumber();
        int dfvar = block.Size(dfseq);
        // pos : position of the block in the solution matrix
        int64_t pos = block.Position(dfseq);

        /// ish loops of the number of shape functions associated with the block
        for(int ish=0; ish<dfvar/nstate; ish++)
        {
            ivec    = data.fVecShapeIndex[counter].first;
            ishape  = data.fVecShapeIndex[counter].second;


            // portion of the gradient coming from the gradient of the scalar function
            for (int e = 0; e < dim; e++) {
                for (int f = 0; f< dim; f++) {
                    GradOfPhiHdiv(e,f) = Normalvec(e,ivec)*dphix(f,ishape);
                }
            }

            for (int64_t is=0; is<numbersol; is++)
            {
                for(int idf=0; idf<nstate; idf++)
                {
                    TVar meshsol = MeshSol(pos+ish*nstate+idf,is);
                    REAL phival = data.phi(ishape,0);
                    TPZManVector<REAL,3> normal(3);

                    for (int i=0; i<3; i++)
                    {
                        if (data.fNeedsDeformedDirectionsFad) {
                            normal[i] = data.fDeformedDirectionsFad(i,ivec).val();
                        }else{
                            normal[i] = data.fDeformedDirections(i,ivec);
                        }
                    }

#ifdef PZ_LOG
                    if(logger.isDebugEnabled() && abs(meshsol) > 1.e-6)
                    {
                        std::stringstream sout;
                        sout << "meshsol = " << meshsol << " ivec " << ivec << " ishape " << ishape << " x " << data.x << std::endl;
                        sout << " phi = " << data.phi(ishape,0) << " dphix " << dphix(0,ishape) << " " << dphix(1,ishape) << std::endl;
                        sout << "normal = " << normal << std::endl;
                        sout << "GradOfPhiHdiv " << GradOfPhiHdiv << std::endl;
                        sout << "GradNormalVec " << GradNormalvec[ivec] << std::endl;
                        LOGPZ_DEBUG(logger,sout.str())
                    }
#endif

                    data.divsol[is][idf] += data.divphi(counter,0)*meshsol;
                    for (int ilinha=0; ilinha<dim; ilinha++) {
                        data.sol[is][ilinha+dim*idf] += normal[ilinha]*phival*meshsol;
                        for (int kdim = 0 ; kdim < dim; kdim++) {
                            data.dsol[is](ilinha+dim*idf,kdim)+= meshsol * GradOfPhiHdiv(ilinha,kdim);
                            if(data.fNeedsDeformedDirectionsFad){
                                data.dsol[is](ilinha+dim*idf,kdim)+=meshsol *GradNormalvec[ivec](ilinha,kdim)*data.phi(ishape,0);
                            }
                        }

                    }

                }
            }
            counter++;
        }
    }






#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "x " << data.x << " sol " << data.sol[0] << std::endl;
        data.dsol[0].Print("dsol",sout);
        sout << "divsol" << data.divsol[0] << std::endl;
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif

}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12)
{

	bool Is_u1PHI = (u1.Cols() == 1) ? true : false;
	bool Is_u2PHI = (u2.Cols() == 1) ? true : false;

	if(Is_u1PHI && Is_u2PHI)
	{
		int64_t nu1 = u1.Rows(),nu2 = u2.Rows();
		u12.Redim(nu1+nu2,1);
		int64_t i;
		for(i=0; i<nu1; i++) u12(i,0) = u1(i,0);
		for(i=0; i<nu2; i++) u12(i+nu1,0) = u2(i,0);


	}
	else if(!Is_u1PHI || !Is_u2PHI)
	{
		int64_t ru1 = u1.Rows(), cu1 = u1.Cols(), ru2 = u2.Rows(), cu2 = u2.Cols();
		int64_t ru12 = ru1 < ru2 ? ru2 : ru1;
		int64_t cu12 = cu1+cu2;
		u12.Redim(ru12,cu12);
		int64_t i,j;
		for(i=0; i<ru1; i++) for(j=0; j<cu1; j++) u12(i,j) = u1(i,j);
		for(i=0; i<ru2; i++) for(j=0; j<cu2; j++) u12(i,j+cu1) = u2(i,j);
	}
	else
	{
		PZError << "TPZCompElKernelHDiv::Append. Bad input parameters " << std::endl;

	}

}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::FillOrder(TPZVec<int> &order) const
{
    order.resize(TSHAPE::NSides-TSHAPE::NCornerNodes);
    int ncon = TSHAPE::NFacets+1;

    TPZConnect &c = this->Connect(ncon-1);
    int internalorder = c.Order();
    order.Fill(internalorder+1);
#ifdef PZDEBUG
    for(int ic=0; ic<ncon-1; ic++)
    {
        if(ConnectOrder(ic) > internalorder) DebugStop();
    }
#endif
    return;

    int nvecs = TSHAPE::Dimension*TSHAPE::NSides;
    TPZManVector<int,3*27> associated_side(nvecs),bilinear(nvecs),direction(nvecs);
    TSHAPE::GetSideHDivDirections(associated_side,direction,bilinear);
    TPZManVector<int,27> sideinc(TSHAPE::NSides,0);
    for (int iv=0; iv<nvecs; iv++) {
        int side = associated_side[iv];
        int bil = bilinear[iv];
        if (bil) {
            sideinc[side] = 1;
        }
    }
    int nsides = TSHAPE::NSides;
    for (int is=0; is<nsides; is++) {
        if (TSHAPE::SideDimension(is) < TSHAPE::Dimension-1) {
            continue;
        }
        else if(TSHAPE::SideDimension(is) == TSHAPE::Dimension -1)
        {
            int intorder = internalorder;
            int connectindex = SideConnectLocId(0, is);
            if (connectindex < 0) {
                connectindex = SideConnectLocId(0,is);
                DebugStop();
            }
            TPZConnect &c = this->Connect(connectindex);
            if (c.Order() > intorder) {
                DebugStop();
                intorder = c.Order();
            }
            if (sideinc[is]) {
                intorder++;
            }
            order[is-TSHAPE::NCornerNodes] = intorder;
        }
        else
        {
            int intorder = internalorder;
            if (sideinc[is]) {
                intorder++;
            }
            order[is-TSHAPE::NCornerNodes] = intorder;
        }
    }
}


template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
	TPZManVector<int64_t,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);
	TPZManVector<int, TSHAPE::NSides-TSHAPE::NCornerNodes+1> ord(TSHAPE::NSides-TSHAPE::NCornerNodes,0);
    int i;
    TPZGeoEl *ref = this->Reference();
    for(i=0; i<TSHAPE::NCornerNodes; i++) {
        id[i] = ref->NodePtr(i)->Id();
    }


    FillOrder(ord);
    int nshape= this->NShapeContinuous(ord);

//    phi.Redim(nshape, 1);
//    dphi.Redim(TSHAPE::Dimension, nshape);
    TSHAPE::Shape(pt,id,ord,phi,dphi);

}

template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::NShapeContinuous(TPZVec<int> &order ){

    return TSHAPE::NShapeF(order);
//#ifdef PZ_LOG
//		{
//				std::stringstream sout;
//				sout << "ordem max "<<maxorder<< " vec order " << order<<" num func cont "<< nshape<<std::endl;
//				LOGPZ_DEBUG(logger,sout.str())
//		}
//#endif


}


template<class TSHAPE>
TPZTransform<> TPZCompElKernelHDiv<TSHAPE>::TransformSideToElement(int side){
	return TSHAPE::TransformSideToElement(side);
}

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int64_t> &shapeindex){

    DebugStop();
	TPZManVector<int64_t> firstshapeindex;
	FirstShapeIndex(firstshapeindex);
	int nshape = TPZIntelGen<TSHAPE>::NShapeF();
	shapeindex.Resize(nshape);
	int64_t nsides = sides.NElements();
	int64_t is, count=0;
	for(is=0 ; is<nsides; is++)
	{
		int side = sides[is];
        int conind = SideConnectLocId(0, side);
		int sideorder = this->Connect(conind).Order();
		int NShapeFace = TSHAPE::NConnectShapeF(side,sideorder);
		int ishapeface;
		for(ishapeface=0; ishapeface<NShapeFace; ishapeface++)
		{
			shapeindex[count++] = is;
		}
	}
	shapeindex.Resize(count);
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "count = " << count << " nshape " << nshape;
		sout << endl<<"sides associated with the normals "<< sides <<
		"\nnormal associated with each shape function : shape function indexes " << shapeindex;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

template<class TSHAPE>
template<class TVar>
void TPZCompElKernelHDiv<TSHAPE>::ComputeRequiredDataT(TPZMaterialDataT<TVar> &data,
                                                TPZVec<REAL> &qsi){

//    TPZManVector<int,TSHAPE::NSides*TSHAPE::Dimension> normalsidesDG(TSHAPE::Dimension*TSHAPE::NSides);

    bool needsol = data.fNeedsSol;
    data.fNeedsSol = false;
    TPZIntelGen<TSHAPE>::ComputeRequiredData(data,qsi);
    data.fNeedsSol = needsol;



    if(data.fNeedsDeformedDirectionsFad){
        TPZIntelGen<TSHAPE>::Reference()->HDivDirections(qsi,data.fDeformedDirectionsFad);
        int cont = 0;
        int firstface = TSHAPE::NSides - TSHAPE::NFacets - 1;
        int lastface = TSHAPE::NSides - 1;
        for(int side = firstface; side < lastface; side++)
        {
            int nvec = TSHAPE::NContainedSides(side);
            for (int ivet = 0; ivet<nvec; ivet++)
            {
                for (int il = 0; il<3; il++)
                {
                    data.fDeformedDirectionsFad(il,ivet+cont) *= fSideOrient[side-firstface];
                }

            }
            cont += nvec;
        }
    }else{
        TPZIntelGen<TSHAPE>::Reference()->HDivDirections(qsi,data.fDeformedDirections);
        int cont = 0;
        int firstface = TSHAPE::NSides - TSHAPE::NFacets - 1;
        int lastface = TSHAPE::NSides - 1;

        for(int side = firstface; side < lastface; side++)
        {
            int nvec = TSHAPE::NContainedSides(side);
            for (int ivet = 0; ivet<nvec; ivet++)
            {
                for (int il = 0; il<3; il++)
                {
                    data.fDeformedDirections(il,ivet+cont) *= fSideOrient[side-firstface];
                }
            }
            cont += nvec;
        }
    }

    data.ComputeFunctionDivergence();
    if (data.fNeedsSol) {
        constexpr bool hasPhi{true};
        ReallyComputeSolution(data);
    }


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
//	if (TSHAPE::Type()==EQuadrilateral) {
//        int maxorder = this->MaxOrder();
//        data.p = maxorder+1;
//    }
    {
        TPZManVector<int> orders;
        FillOrder(orders);
        int nshapescalar = TSHAPE::NShapeF(orders);
        data.phi.Resize(nshapescalar, 1);
        data.dphi.Resize(TSHAPE::Dimension, nshapescalar);
        data.dphix.Resize(TSHAPE::Dimension, nshapescalar);
    }
#ifdef PZ_LOG
        if(logger.isDebugEnabled())
		{
				LOGPZ_DEBUG(logger,"Initializing MaterialData of TPZCompElKernelHDiv")
		}
#endif
    TPZManVector<int,TSHAPE::Dimension*TSHAPE::NSides> vecside(TSHAPE::Dimension*TSHAPE::NSides),bilinear(TSHAPE::Dimension*TSHAPE::NSides),directions(TSHAPE::Dimension*TSHAPE::NSides);

	TPZManVector<int,TSHAPE::NSides*TSHAPE::Dimension> normalsides(TSHAPE::Dimension*TSHAPE::NSides);

    {
        int nconnects = TSHAPE::NFacets+1;
        int internalorder = this->Connect(nconnects-1).Order();
        TPZVec<std::pair<int,int64_t> > IndexVecShape;
        TSHAPE::GetSideHDivDirections(vecside,directions,bilinear,normalsides);
        int64_t numvec = TSHAPE::Dimension*TSHAPE::NSides;

        data.fMasterDirections.Resize(3, numvec);
        this->Reference()->HDivDirectionsMaster(data.fMasterDirections);
        // Acerta o vetor data.fDeformedDirections para considerar a direcao do campo. fSideOrient diz se a orientacao e de entrada
        // no elemento (-1) ou de saida (+1), dependedo se aquele lado eh vizinho pela direita (-1) ou pela esquerda(+1)
        int firstface = TSHAPE::NSides - TSHAPE::NFacets - 1;
        int lastface = TSHAPE::NSides - 1;
        int cont = 0;
        for(int side = firstface; side < lastface; side++)
        {
            int nvec = TSHAPE::NContainedSides(side);
            for (int ivet = 0; ivet<nvec; ivet++)
            {
                for (int il = 0; il<3; il++)
                {
                    data.fMasterDirections(il,ivet+cont) *= SideOrient(side-firstface);
                }
            }
            cont += nvec;
        }

        data.fDeformedDirectionsFad.Resize(3, numvec);

        data.fDeformedDirections.Resize(3, numvec);


        // IndexShapeToVec2(normalsides, bilinear, directions,data.fVecShapeIndex,internalorder);
        data.divphi.Resize(data.fVecShapeIndex.size(), 1);
    }
    data.fShapeType = TPZMaterialData::EVecandShape;

//    cout << "vecShape " << endl;
//    cout << data.fVecShapeIndex << endl;;


#ifdef PZ_LOG
    if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		data.fDeformedDirections.Print("Normal vector ", sout,EMathematicaInput);
        for (int i=0; i<TSHAPE::NCornerNodes; i++) {
            sout << "Id[" << i << "] = " << this->Reference()->NodePtr(i)->Id() << " ";
        }
        sout << "\n\nSides associated with the normals\n";
        for (int i=0; i<normalsides.size(); i++) {
            sout << i << '|' << normalsides[i] << " ";
        }
        sout << std::endl;

        sout << std::endl;
		sout << "NormalVector/Shape indexes \n";
        for (int i=0; i<data.fVecShapeIndex.size(); i++) {
            sout << i << '|' << data.fVecShapeIndex[i] << " ";
        }
        sout << std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif

}


// Save the element data to a stream
template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::Write(TPZStream &buf, int withclassid) const
{
	TPZInterpolatedElement::Write(buf,withclassid);
	TPZManVector<int,3> order(3,0);
	this->fIntRule.GetOrder(order);
	buf.Write(order);
    buf.Write(fSideOrient);

	buf.Write(this->fConnectIndexes.begin(),TSHAPE::NSides);
	buf.Write(&this->fPreferredOrder,1);
    buf.Write(fSideOrient);
    int sz = fRestraints.size();
    buf.Write(&sz);
    for (std::list<TPZOneShapeRestraint>::const_iterator it = fRestraints.begin(); it != fRestraints.end(); it++) {
        it->Write(buf);
    }
	int classid = this->ClassId();
	buf.Write ( &classid, 1 );
}


// Read the element data from a stream
template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZInterpolatedElement::Read(buf,context);
	TPZManVector<int,3> order;
	buf.Read(order);
	this-> fIntRule.SetOrder(order);
    TPZManVector<int, TSHAPE::NFacets> SideOrient;
    buf.Read(SideOrient);
    fSideOrient = SideOrient;
	buf.Read(this->fConnectIndexes.begin(),TSHAPE::NSides);
	buf.Read(&this->fPreferredOrder,1);
    buf.Read(fSideOrient);
    int sz;
    buf.Read(&sz);
    for (int i=0; i<sz; i++) {
        TPZOneShapeRestraint one;
        one.Read(buf);
        fRestraints.push_back(one);
    }
	int classid = -1;
	buf.Read( &classid, 1 );
	if ( classid != this->ClassId())
	{
		std::stringstream sout;
		sout << "ERROR - " << __PRETTY_FUNCTION__
        << " trying to restore an object id " << this->ClassId() << " and classid read = " << classid;
		LOGPZ_ERROR ( logger, sout.str().c_str() );
	}
}
//refinamento
template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::PRefine(int order)
{
    this->SetPreferredOrder(order);
    int side;
    int icon;
    int ncon=NConnects();
    TPZCompElHDivPressure<TSHAPE> *hdivpressure = dynamic_cast<TPZCompElHDivPressure<TSHAPE> *>(this);

    if (hdivpressure) {
        ncon--;
    }
    int nnodes = this->Reference()->NNodes();
    for(icon=0; icon<ncon; icon++)
    {//somente para os conects de fluxo
//        TPZConnect &con = this->Connect(icon);
//        con.SetOrder(order);
        side= ConnectSideLocId(icon);

#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
                std::stringstream sout;
                sout << "side " << side << " order " << this->PreferredSideOrder(side)<<std::endl;
                LOGPZ_DEBUG(logger,sout.str())
        }
#endif

        this->IdentifySideOrder(side);
    }
    #ifdef PZ_LOG
    if (loggerdiv.isDebugEnabled()) {
        std::stringstream sout;
        sout << (void*) this->Mesh() << "PRefine elindex " << this->Index() << " gel index " << this->Reference()->Index() << " " << order;
        sout << "\nPRefine connect orders ";
        int nc = this->NConnects();
        for(int ic=0; ic<nc; ic++) sout << (int)this->Connect(ic).Order() << " ";
        LOGPZ_DEBUG(loggerdiv, sout.str())
    }
#endif

		// conect da pressao

    if(ncon>nnodes+1)
    {
		TPZCompElHDivPressure<TSHAPE> *hdivpressure = dynamic_cast<TPZCompElHDivPressure<TSHAPE> *>(this);
		TPZConnect &con = this->Connect(ncon-1);

		if (TSHAPE::Type()==EQuadrilateral) {
				hdivpressure->SetPressureOrder(order);
				con.SetOrder(order,this->fConnectIndexes[ncon-1]);

		}
		else {
				hdivpressure->SetPressureOrder(order-1);
				con.SetOrder(order-1,this->fConnectIndexes[ncon-1]);

		}
		int nshape = hdivpressure-> NConnectShapeF(ncon-1,con.Order());
		con.SetNShape(nshape);
		int64_t seqnum = con.SequenceNumber();
		this->Mesh()->Block().Set(seqnum,nshape);
    }


}

/** @brief Prints the relevant data of the element to the output stream */
template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::Print(std::ostream &out) const
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZIntelGen<TSHAPE>::Print(out);
    out << "Side orientation " << fSideOrient << std::endl;
    if (fRestraints.size()) {
        out << "One shape restraints associated with the element\n";
        for (std::list<TPZOneShapeRestraint>::const_iterator it = fRestraints.begin(); it != fRestraints.end(); it++)
        {
            it->Print(out);
        }
    }



}

template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::MaxOrder(){

    int maxorder = TPZInterpolationSpace::MaxOrder();
    return maxorder+1;
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


template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == TSHAPE::Dimension && this->Material()->Id() > 0) {
		new typename TSHAPE::GraphElType(this,&grafgrid);
	}
}

/// return the first one dof restraint
template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::RestrainedFace()
{
    return -1;
}

template<>
int TPZCompElKernelHDiv<TPZShapePiram>::RestrainedFace()
{
    if (fRestraints.size() == 0) {
        return -1;
        DebugStop(); //AQUIPHIL
    }
    std::list<TPZOneShapeRestraint>::iterator it = fRestraints.begin();
    int foundis = -1;
    bool found = false;
    while (found == false && it != fRestraints.end()) {
        int64_t connectindex = it->fFaces[3].first;
        int64_t cindex = -1;
        for (int is = 14; is<18; is++) {
            cindex = ConnectIndex(is-13);
            if (connectindex == cindex) {
                found = true;
                foundis = is+1;
                if (foundis == 18) {
                    foundis = 14;
                }
            }
        }
        it++;
    }
    if (found == false) {
        DebugStop();
    }
    return foundis;
}

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


// TPZCompEl * CreateHDivBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
// 	return new TPZCompElHDivBound2<TPZShapePoint>(mesh,gel,index);
// }

// TPZCompEl * CreateHDivBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
// 	return new TPZCompElHDivBound2< TPZShapeLinear>(mesh,gel,index);
// }

// TPZCompEl * CreateHDivBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
//     return new TPZCompElHDivBound2< TPZShapeQuad>(mesh,gel,index);
// }

// TPZCompEl * CreateHDivBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
//     return new TPZCompElHDivBound2< TPZShapeTriang >(mesh,gel,index);
// }

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

