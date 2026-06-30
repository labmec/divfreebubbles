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
    enum ProblemOrigin {EDarcyHDiv, EElasticityHDiv, EDarcyH1Hybrid, EElasticityH1Hybrid};

    TPZMatRedSolver() = default;

    TPZMatRedSolver(TPZLinearAnalysis &an, ProblemOrigin pOrigin){
        fAnalysis = &an;
        fProblemOrigin = pOrigin;
    }


    void Solve(std::ostream &out = std::cout);

    void ComputeConditionNumber(TPZSparseMatRed<STATE> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond);
    void ComputeConditionNumber(TPZMatRed<STATE,TPZFMatrix<STATE>> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond);

protected:
    ProblemOrigin fProblemOrigin;

    TPZLinearAnalysis *fAnalysis;


};

#endif