/**
 * @file
 * @brief Contains declaration of TPZCompElHDivSemiHybridBound
 */

#ifndef PZELCHDIVBOUND_SEMI_HYBRID
#define PZELCHDIVBOUND_SEMI_HYBRID

#include "pzelchdivbound2.h"

/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivSemiHybridBound : public TPZCompElHDivBound2<TSHAPE> {

    /// vector which defines whether the normal is outward or not
    TPZManVector<int, TSHAPE::NFacets> fSideOrient;

public:

    TPZCompElHDivSemiHybridBound(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam = DefaultFamily::fHDivDefaultValue);

    virtual int NConnects() const override;

     /**
     * @brief Number of shapefunctions of the connect associated
     * @param connect connect number
     * @return number of shape functions
     */
	virtual int NConnectShapeF(int connect, int order) const override;
	

    virtual int NSideConnects(int side) const override;

    virtual int64_t ConnectIndex(int con) const override;

    virtual int SideConnectLocId(int node, int side) const override;

    void SetConnectIndex(int i, int64_t connectindex) override;
protected:
    
	
	
};



#endif
