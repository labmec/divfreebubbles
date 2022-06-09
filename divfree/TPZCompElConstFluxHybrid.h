
#ifndef PZELC_CONSTFLUX_HYBRID
#define PZELC_CONSTFLUX_HYBRID

#include "TPZCompElDisc.h"
#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatLoadCases.h"
#include "TPZMatErrorSingleSpace.h"

class TPZCompElConstFluxHybrid : public TPZCompElDisc{

public:

    TPZCompElConstFluxHybrid(TPZCompMesh &mesh, TPZGeoEl *reference);

    void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override{
        CalcStiffInternal(ek,ef);
    }
    void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef) override;

    virtual int NConnects() const override;

     /**
     * @brief Number of shapefunctions of the connect associated
     * @param connect connect number
     * @return number of shape functions
     */
	virtual int NConnectShapeF(int connect, int order) const override;
	
    virtual int NSideConnects(int side) const override;

    virtual int64_t ConnectIndex(int con) const override;

    /** 
     * @brief return the local index for connect
	 **/
	virtual int SideConnectLocId(int node, int side) const override;
    
    void SetConnectIndex(int i, int64_t connectindex) override;    

    void SetSideOrient(int sideorient){
        fSideOrient = sideorient;
    };

    int GetSideOrient();

protected:
	
	/** @brief It preserves index of connect associated to the element */
	TPZManVector<int64_t,2> fConnectIndexes;

    int fSideOrient = 1;

    template<class TVar>
    void CalcStiffInternal(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef);
};

#endif