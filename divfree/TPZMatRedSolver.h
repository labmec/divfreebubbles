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

    TPZMatRedSolver(TPZLinearAnalysis &an, std::set<int> &matIdBC, SolverType sType = EDefault){
        fAnalysis = &an;
        fBCMaterialID = &matIdBC;
        fSolverType = sType;
    };

    void Solve(std::ostream &out = std::cout);

    void SolveProblemDefault(int64_t &nEqLinr, int64_t &nEqHigh, TPZStructMatrix &Stiffness, TPZFMatrix<TVar> &rhsFull, TPZFMatrix<TVar> &rhsHigh, std::ostream &out);

    void SolveProblemSparse(int64_t &nEqLinr, int64_t &nEqHigh, TPZStructMatrix &Stiffness, TPZFMatrix<TVar> &rhsFull, TPZFMatrix<TVar> &rhsHigh, std::ostream &out); 
    
    void SolveProblemSparseNew(int64_t &nEqLinr, int64_t &nEqHigh, TPZStructMatrix &Stiffness, TPZFMatrix<TVar> &rhsFull, TPZFMatrix<TVar> &rhsHigh, std::ostream &out); 



protected:
    SolverType fSolverType;

    TPZLinearAnalysis *fAnalysis;

    std::set<int> *fBCMaterialID;
};

#endif