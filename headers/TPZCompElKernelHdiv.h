/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef TPZCOMPELKERNELHDIV_H
#define TPZCOMPELKERNELHDIV_H

#include "pzelctemp.h"
#include "TPZOneShapeRestraint.h"


/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElKernelHDiv : public TPZIntelGen<TSHAPE> {
	
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
	
	TPZCompElKernelHDiv();
	
	virtual ~TPZCompElKernelHDiv();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const  override {
		return new TPZCompElKernelHDiv<TSHAPE> (mesh, *this);
	}
	
	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes
	 * mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap) const override
	{
		return new TPZCompElKernelHDiv<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
    /** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh) override;
	

	
	virtual MElementType Type() override;
		
	virtual void SetConnectIndex(int i, int64_t connectindex) override;
    
    /**
     * @brief Number of shapefunctions of the connect associated
     * @param connect connect number
     * @return number of shape functions
     */
	virtual int NConnectShapeF(int connect, int order) const override;
	
	virtual int Dimension() const  override {
		return TSHAPE::Dimension;
	}
	
	virtual int NCornerConnects() const override {
		return 0;
	}
	
	virtual int NSideConnects(int side) const override;
    
	/** 
     * @brief return the local index for connect
	 **/
	virtual int SideConnectLocId(int node, int side) const override;
    
	
	virtual int64_t ConnectIndex(int con) const override;
    
    /// Add a shape restraint (meant to fit the pyramid to restraint
    virtual void AddShapeRestraint(TPZOneShapeRestraint restraint) override
    {
        fRestraints.push_back(restraint);
    }
    
    /// Return a list with the shape restraints
    virtual std::list<TPZOneShapeRestraint> GetShapeRestraints() const override
    {
        return fRestraints;
    }
    
    /// Return a list with the shape restraints
    virtual void ResetShapeRestraints() override
    {
        fRestraints.clear();
    }
	
	virtual void SetIntegrationRule(int ord) override;
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord) override;
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
	virtual int PreferredSideOrder(int iside) override;
	
	/*
     * @brief Sets the preferred interpolation order along a side \n
	 * This method only updates the datastructure of the element
	 * In order to change the interpolation order of an element, use the method PRefine
	 */
	virtual void SetPreferredOrder(int order) override;
	
	/** @brief Sets the interpolation order of side to order*/
	// virtual void SetSideOrder(int side, int order) override;
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
	// virtual int EffectiveSideOrder(int side) const override;
	
    /**
     * @brief return the interpolation order of the polynomial for connect
     **/
	virtual int ConnectOrder(int connect) const override;
	
    
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

	/** @brief Compute the correspondence between the normal vectors and the shape functions */
	void ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int64_t> &shapeindex);
	
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) override;

	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) override;
		
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
    int ClassId() const override;

	/** @brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
protected:
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
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"

/** @brief Creates computational point element for HDiv approximate space */
TPZCompEl *CreateKernelHDivPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateKernelHDivLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateKernelHDivQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateKernelHDivTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational cube element for HDiv approximate space */
TPZCompEl *CreateKernelHDivCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational prismal element for HDiv approximate space */
TPZCompEl *CreateKernelHDivPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational pyramidal element for HDiv approximate space */
TPZCompEl *CreateKernelHDivPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
/** @brief Creates computational tetrahedral element for HDiv approximate space */
TPZCompEl *CreateKernelHDivTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);

TPZCompEl * CreateRefKernelHDivLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefKernelHDivQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefKernelHDivTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefKernelHDivCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefKernelHDivPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefKernelHDivPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
TPZCompEl * CreateRefKernelHDivTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);


#endif
