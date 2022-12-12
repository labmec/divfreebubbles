//
// Created by Jeferson Fernandes on 12/05/22.
//

#ifndef TPZ_MATRED_SOLVER_H
#define TPZ_MATRED_SOLVER_H

#include "pzmatred.h"
#include "TPZSparseMatRed.h"
#include "TPZLinearAnalysis.h"

template <class TVar>
class TPZMatRedSolver {
public:
    enum SolverType{EDefault, ESparse}; 

    TPZMatRedSolver() = default;

    TPZMatRedSolver(TPZLinearAnalysis &an, std::set<int> &matIdBC, SolverType sType = EDefault, std::function<TPZManVector<STATE,3>(const TPZVec<REAL> &coord)> permFunction = nullptr){
        fAnalysis = &an;
        fBCMaterialID = &matIdBC;
        fSolverType = sType;
        fPermFunction = permFunction;
    };

    void Solve(std::ostream &out = std::cout);

    void SolveProblemDefault(std::ostream &out);

    void SolveProblemSparse(std::ostream &out); 
    
    void ComputeConditionNumber(TPZSparseMatRed<STATE> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond);
    void ComputeConditionNumber(TPZMatRed<STATE,TPZFMatrix<STATE>> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond);
    void ThresholdPermeability(REAL threshold);

protected:
    SolverType fSolverType;

    TPZLinearAnalysis *fAnalysis;

    std::set<int> *fBCMaterialID;

    TPZVec<int64_t> fActiveEquations;

    std::function<TPZManVector<STATE,3>(const TPZVec<REAL> &coord)> fPermFunction;
};

#endif