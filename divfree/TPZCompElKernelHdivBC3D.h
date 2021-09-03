/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef TPZCompElKernelHDivBC3D_H
#define TPZCompElKernelHDivBC3D_H

#include "pzelctemp.h"
#include "TPZBndCond.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"


/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElKernelHDivBC3D : public TPZIntelGen<TSHAPE>  {

public:
	    
	TPZCompElKernelHDivBC3D();
    
    TPZCompElKernelHDivBC3D(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
	virtual ~TPZCompElKernelHDivBC3D();

    virtual void InitMaterialData(TPZMaterialData &data) override;

    void ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) override;
	
};


#endif
