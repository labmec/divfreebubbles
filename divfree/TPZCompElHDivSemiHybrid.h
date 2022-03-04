/**
 * @file
 * @brief Contains declaration of TPZCompElHDivSemiHybrid
 */

#ifndef PZELCHDIV_SEMI_HYBRID
#define PZELCHDIV_SEMI_HYBRID

#include "pzelchdiv.h"

/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivSemiHybrid : public TPZCompElHDiv<TSHAPE> {

    /// vector which defines whether the normal is outward or not
    TPZManVector<int, TSHAPE::NFacets> fSideOrient;

public:

    TPZCompElHDivSemiHybrid(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam = DefaultFamily::fHDivDefaultValue);

    virtual int NConnects() const override;

     /**
     * @brief Number of shapefunctions of the connect associated
     * @param connect connect number
     * @return number of shape functions
     */
	virtual int NConnectShapeF(int connect, int order) const override;
	

    virtual int NSideConnects(int side) const override;

    virtual int64_t ConnectIndex(int con) const override;

protected:
    
	
	
};






/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivSemiHybridLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivSemiHybridQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivSemiHybridTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational cube element for HDiv approximate space */
TPZCompEl *CreateHDivSemiHybridCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational tetrahedral element for HDiv approximate space */
TPZCompEl *CreateHDivSemiHybridTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);


#endif
