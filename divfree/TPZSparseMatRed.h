/**
 * @file
 * @brief Contains TPZSparseMatRed class which implements a simple substructuring of a linear system of equations, composed of 4 submatrices.
 */

#ifndef _TSPARSEMATRED_
#define _TSPARSEMATRED_


#include "pzmatred.h"

/**
 * @brief Implements a simple substructuring of a linear system of equations, composed of 4 submatrices. \ref matrix "Matrix"
 * @ingroup matrix
 */
/**
 * Implements a matrix composed of 4 submatrices:
 *			\f[ [K00] [U0] + [K01] [U1] = [F0] \f]
 *			\f[ [K10] [U0] + [K11] [U1] = [F1] \f]
 */

template<class TVar, class TSideMatrix >
class TPZSparseMatRed: public TPZMatRed<TVar,TSideMatrix>
{
public:
	
	// friend class TPZSparseMatRed<TVar, TPZFMatrix<TVar> >;
	// friend class TPZSparseMatRed<TVar ,TPZVerySparseMatrix<TVar> >;
	/** @brief Simple constructor */
	TPZSparseMatRed();
	
	/**
	 * @brief Constructor with 2 parameters
	 * @param dim assumes the value of n1+n2
	 * @param dim00 equals n1
	 */
	// TPZSparseMatRed(const int64_t dim, const int64_t dim00);
	
	template<class TSideCopy>
	TPZSparseMatRed<TVar ,TSideMatrix>(const TPZSparseMatRed<TVar, TSideCopy> &cp): TPZMatrix<TVar>(cp)
	{
        this->fK11 = cp.fK11;
        this->fK01 = cp.fK01;
        this->fK10 = cp.fK10;
        this->fF0 = cp.fF0;
        this->fF1 = cp.fF1;
        this->fMaxRigidBodyModes = cp.fMaxRigidBodyModes;
        this->fNumberRigidBodyModes = cp.fNumberRigidBodyModes;
        this->fF0IsComputed = cp.fF0IsComputed;
		this->fDim0=cp.fDim0;
		this->fDim1=cp.fDim1;
		this->fK01IsComputed = cp.fK01IsComputed;
		this->fIsReduced = cp.fIsReduced;
		this->fSolver = cp.fSolver;
		
		if(cp.fK00) this->fK00 = cp.fK00;
	}
	inline TPZSparseMatRed<TVar,TSideMatrix>*NewMatrix() const override {return new TPZSparseMatRed<TVar,TSideMatrix>{};}
	CLONEDEF(TPZSparseMatRed)

    /** @brief Saveable methods */
    int ClassId() const override;

    
  
};

template<class TVar, class TSideMatrix>
int TPZSparseMatRed<TVar,TSideMatrix>::ClassId() const{
    return Hash("TPZSparseMatRed") ^ TPZMatrix<TVar>::ClassId() << 1 ^ TSideMatrix().ClassId() << 2;
}


#endif
