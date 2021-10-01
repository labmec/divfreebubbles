/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef TPZCOMPELKERNELHDIV_H
#define TPZCOMPELKERNELHDIV_H

#include "pzelctemp.h"
#include "TPZOneShapeRestraint.h"
#include "TPZCompElH1.h"

/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElKernelHDiv : public TPZCompElH1<TSHAPE> {

    /// vector which defines whether the normal is outward or not
    TPZManVector<int, TSHAPE::NFacets> fSideOrient;
    
    /// Data structure which defines the restraints
    std::list<TPZOneShapeRestraint> fRestraints;

public:
	    
	TPZCompElKernelHDiv(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
	
	TPZCompElKernelHDiv(TPZCompMesh &mesh, const TPZCompElKernelHDiv<TSHAPE> &copy);
	
	/**
	 * @brief Constructor used to generate patch mesh... generates a map of connect index from
	 * global mesh to clone mesh
	 */
	TPZCompElKernelHDiv(TPZCompMesh &mesh,
				        const TPZCompElKernelHDiv<TSHAPE> &copy,
				        std::map<int64_t,int64_t> & gl2lcConMap,
				        std::map<int64_t,int64_t> & gl2lcElMap);
	
	TPZCompElKernelHDiv(){};
	
	virtual ~TPZCompElKernelHDiv();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const  override {
		return new TPZCompElKernelHDiv<TSHAPE> (mesh, *this);
	}
	
    /** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh) override;
		
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	virtual void InitMaterialData(TPZMaterialData &data) override;

    //@{
	/** @brief Compute and fill data with requested attributes */
	void ComputeRequiredData(TPZMaterialDataT<STATE> &data,
                             TPZVec<REAL> &qsi) override{
        ComputeRequiredDataT(data,qsi);
    }
     void ComputeRequiredData(TPZMaterialDataT<CSTATE> &data,
                              TPZVec<REAL> &qsi) override{
        ComputeRequiredDataT(data,qsi);
    }
    //@}
		
	/** @brief Compute the solution for a given variable */
	virtual void Solution( TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override;

	/** @brief Returns the unique identifier for reading/writing objects to streams */
    int ClassId() const override;

    /**
     * @brief It returns the normal orientation of the reference element by the side.
     * Only side that has dimension larger than zero and smaller than me.
     * @param side: side of the reference elemen
     */
    virtual int GetSideOrient(int side) override;
    
    /**
     * @brief It set the normal orientation of the element by the side.
     * Only side that has dimension equal to my dimension minus one.
     * @param side: side of the reference elemen
     */
    virtual void SetSideOrient(int side, int sideorient) override;

    /// the orientation of the face
    int SideOrient(int face)
    {
#ifdef PZDEBUG
        if (face < 0 || face >= TSHAPE::NFacets) {
            DebugStop();
        }
#endif
        return fSideOrient[face];
    }
	


protected:

	 //@{
    /** @brief Compute the solution using Hdiv structure */
	void ReallyComputeSolution(TPZMaterialDataT<STATE> &data) override{
        ComputeSolutionKernelHdivT(data);
    }
    void ReallyComputeSolution(TPZMaterialDataT<CSTATE> &data) override{
        ComputeSolutionKernelHdivT(data);
    }

	template<class TVar>
    void ComputeSolutionKernelHdivT(TPZMaterialDataT<TVar> &data);
    template<class TVar>
    void ComputeRequiredDataT(TPZMaterialDataT<TVar> &data, TPZVec<REAL>&qsi);
};


template<class TSHAPE>
int TPZCompElKernelHDiv<TSHAPE>::ClassId() const{
    return Hash("TPZCompElKernelHDiv") ^ TPZIntelGen<TSHAPE>::ClassId() << 1;
}

#include "pzcmesh.h"

template<class TSHAPE>
void TPZCompElKernelHDiv<TSHAPE>::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->SetAllCreateFunctionsContinuous();
}


#endif
