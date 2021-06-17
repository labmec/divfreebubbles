/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef TPZCOMPELKERNELHDIVBC_H
#define TPZCOMPELKERNELHDIVBC_H

#include "TPZBndCond.h"
#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatLoadCases.h"
#include "TPZMatErrorSingleSpace.h"
#include "TPZMaterialDataT.h"
#include "TPZMatCombinedSpaces.h"



/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TVar=STATE>
class TPZCompElKernelHDivBC : public TPZMatCombinedSpacesBC<TVar> {
    using TBase = TPZMatBase<TVar,
                             TPZMatSingleSpaceT<TVar>,
                             TPZMatErrorSingleSpace<TVar>,
                             TPZMatLoadCases<TVar>>;

public:
	    
	TPZCompElKernelHDivBC();
    
    TPZCompElKernelHDivBC(int id, int dim){};

    TPZCompElKernelHDivBC(TPZMaterial * material, int matid, int type, TPZFMatrix<TVar> &val1,TPZFMatrix<TVar> &val2);
	
	virtual ~TPZCompElKernelHDivBC();

    void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<TVar> &ek,TPZFMatrix<TVar> &ef, TPZBndCondT<TVar> &bc);

    /** @brief Returns the unique identifier for reading/writing objects to streams */
    int ClassId() const override;
	
};

template class TPZCompElKernelHDivBC<STATE>;
template class TPZCompElKernelHDivBC<CSTATE>;

#endif
