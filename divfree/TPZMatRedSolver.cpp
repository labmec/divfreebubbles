#include "TPZMatRedSolver.h"
#include "pzcmesh.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "TPZTimer.h"
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "TPZEigenSolver.h"

template<class TVar>
void TPZMatRedSolver<TVar>::Solve(std::ostream &out){
    
    switch (fSolverType)
    {
    case EDefault:
        SolveProblemDefault(out);
        break;
    case ESparse:
        SolveProblemSparse(out);
        break;
    
    default:
        DebugStop();
        break;
    }

}


template<class TVar>
void TPZMatRedSolver<TVar>::SolveProblemDefault(std::ostream &out){

    //HERE STARTS THE ITERATIVE SOLVER SET
    auto cmesh = fAnalysis->Mesh();

    //Compute the number of equations in the system
    int64_t nEqFull = cmesh->NEquations();
    int64_t nEqLinr, nEqHigh;

    //Cria a matriz esparsa
    TPZSparseMatRed<STATE> *matRed2 = new TPZSparseMatRed<STATE>(1,1);

    std::set<int> lag={1};
    matRed2->ReorderEquations(cmesh,lag,nEqFull,nEqLinr);

    nEqHigh = nEqFull-nEqLinr;
    out << nEqHigh << " " << nEqLinr << " ";

    std::cout << "NUMBER OF EQUATIONS:\n " << 
                 "Full problem = " << nEqFull << 
                 ", High Order Flux = " << nEqHigh << 
                 ", Linear Flux = " << nEqLinr << std::endl;

    //Sets number of threads to be used by the solver
    constexpr int nThreads{12};

    // Create the RHS vectors
    TPZFMatrix<STATE> rhsFull(nEqLinr+nEqHigh,1,0.);
    TPZFMatrix<STATE> rhsHigh(nEqHigh,1,0.);

    //Creates the problem matrix    
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Stiffness(fAnalysis->Mesh());
    Stiffness.SetNumThreads(nThreads);

    TPZAutoPointer<TPZGuiInterface> guiInterface;

    //Cria duas matrizes, para inverter a ordem das matrizes em bloco
    TPZMatRed<STATE, TPZFMatrix<STATE>> *matRed = new TPZMatRed<STATE, TPZFMatrix<STATE>>(nEqLinr+nEqHigh,nEqLinr);
    // TPZMatRed<STATE, TPZFMatrix<STATE>> *K00red = new TPZMatRed<STATE, TPZFMatrix<STATE>>(nEqLinr,nEqLinr-25);

    //Primeiro cria a matriz auxiliar
    TPZFMatrix<REAL> K00(nEqLinr,nEqLinr,0.);
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    fAnalysis->SetSolver(step);
    
    //Altera range da matriz stiffness.
    Stiffness.SetEquationRange(0,nEqLinr);
    
    //Transfere as submatrizes da matriz auxiliar para a matriz correta.
    step.SetMatrix(&K00);
    matRed->SetSolver(&step);
  
    Stiffness.EquationFilter().Reset();
    fAnalysis->SetStructuralMatrix(Stiffness);

    TPZTimer clock;
    clock.start();
    //Monta a matriz auxiliar
    rhsFull.Zero();
    std::cout << "Start assembling matRed ...\n";
    Stiffness.Assemble(*matRed,rhsFull,guiInterface);
    std::cout << "Finish assembling matRed ...\n";
    clock.stop();
    // std::cout << "Time Assemble " << clock << std::endl;
    
    

    // TPZFMatrix<REAL> *K11Red = new TPZFMatrix<REAL>(nEqHigh,nEqHigh);

    matRed->SetF(rhsFull);
    matRed->SetReduced();
    // matRed->Print("MATRED",out2,EMathematicaInput);

    // matRed->K11Reduced(*K11Red,rhsHigh);
    
    matRed->F1Red(rhsHigh);
    // rhsFull.Zero();


    TPZBlockDiagonalStructMatrix<STATE> BDFmatrix(fAnalysis->Mesh());
    BDFmatrix.SetEquationRange(nEqLinr,nEqLinr+nEqHigh);
    TPZBlockDiagonal<REAL> KBD;
    BDFmatrix.AssembleBlockDiagonal(KBD);

    //Creates the preconditioner 
    TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>( &KBD );
    precond->SetDirect(ELU);

    int64_t nMaxIter = 50;
    TPZVec<REAL> errors(nMaxIter);
    errors.Fill(0.);

    // for (int64_t iter = 1; iter < nMaxIter; iter++){
        
        // std::cout << "ITER = " << iter << std::endl;
        TPZFMatrix<STATE> residual(nMaxIter,1,0.);
        REAL tol = 1.e-10;
        TPZFMatrix<STATE> solution(nEqHigh,1,0.);
        
        // std::cout << "rhsHigh = " << rhsHigh << std::endl;

        clock.start();
        // K11Red->SolveCG(nMaxIter,*precond,rhsHigh,solution,&residual,tol);
        std::cout << "Start CG ...\n";
        
        matRed->SolveCG(nMaxIter,*precond,rhsHigh,solution,&residual,tol);
        std::cout << "Finish CG ...\n";

        clock.stop();
        // std::cout << "Time CG " << clock << std::endl;
        
        
        // solution.Print("Solution",out,EMathematicaInput);

        REAL norm = 0.;
        std::cout << "Number of CG iterations = " << nMaxIter << " , residual = " << tol << std::endl;
        out << nMaxIter << "\n";
        TPZFMatrix<STATE> result(nEqLinr+nEqHigh,1,0.);
        matRed->UGlobal(solution,result);
        

        //Update solution in Analysis
        // for (int i = 0; i < nEqHigh; i++){
        //     rhsFull(nEqLinr+i,0) = solution(i,0);
        // }
        // for (int i = 0; i < nEqLinr; i++){
        //     rhsFull(i,0) = press(i,0);
        // }

        fAnalysis->Solution()=result;
        fAnalysis->LoadSolution();
        // rhsFull.Print("Solution = ");
                        
        ///Calculating approximation error  
        TPZManVector<REAL,5> error;

        // auto cmeshAux = an.Mesh();
        // int64_t nelem = cmeshAux->NElements();
        // cmeshAux->LoadSolution(cmeshAux->Solution());
        // cmeshAux->ExpandSolution();
        // an.Mesh()->ElementSolution().Redim(an.Mesh()->NElements(), 5);

        // an.PostProcessError(error);
            
        // std::cout << "\nApproximation error:\n";
        // std::cout << "H1 Norm = " << std::scientific << std::setprecision(15) << error[0]<<'\n';
        // std::cout << "L1 Norm = " << std::scientific << std::setprecision(15) << error[1]<<'\n'; 
        // std::cout << "H1 Seminorm = " << std::scientific << std::setprecision(15) << error[2]<<'\n'; 

        // if (tol < 1.e-10) break;
        // errors[iter-1] = error[1];
    // }

    // std::cout << "Number of iterations = " << residual << std::endl;
}

template<class TVar>
void TPZMatRedSolver<TVar>::SolveProblemSparse(std::ostream &out){
    
    //HERE STARTS THE ITERATIVE SOLVER SET
    auto cmesh = fAnalysis->Mesh();

    //Primeiro cria a matriz auxiliar K00 - que será decomposta
    TPZSYsmpMatrix<REAL> K00;
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    fAnalysis->SetSolver(step);
    step.SetMatrix(&K00);

    //Cria a matriz esparsa
    std::set<int> lag = {10};
    TPZSparseMatRed<STATE> *matRed = new TPZSparseMatRed<STATE>(cmesh,lag);
    std::cout << "Allocating Sub Matrices ...\n";
    //Transfere as submatrizes da matriz auxiliar para a matriz correta.
    matRed->SetSolver(&step);
    K00.Resize(matRed->Dim0(),matRed->Dim0());
    matRed->AllocateSubMatrices(cmesh);

    //Compute the number of equations in the system
    int64_t nEqFull = cmesh->NEquations();
    int64_t nEqLinr = matRed->Dim0();
    int64_t nEqHigh = matRed->Dim1();

    out << nEqHigh << " " << nEqLinr << " ";

    std::cout << "NUMBER OF EQUATIONS:\n " << 
                 "Full problem = " << nEqFull << 
                 ", High Order Flux = " << nEqHigh << 
                 ", Linear Flux = " << nEqLinr << std::endl;

    //Sets number of threads to be used by the solver
    constexpr int nThreads{12};
    
    // Create the RHS vectors
    TPZFMatrix<STATE> rhsFull(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsHigh(nEqHigh,1,0.);

    //Creates the problem matrix    
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Stiffness(fAnalysis->Mesh());
    Stiffness.SetNumThreads(nThreads);

    TPZAutoPointer<TPZGuiInterface> guiInterface;
  
    Stiffness.EquationFilter().Reset();
    fAnalysis->SetStructuralMatrix(Stiffness);

    //Monta a matriz auxiliar
    rhsFull.Zero();
    TPZTimer clock;
    clock.start();
    std::cout << "Start assembling matRed ...\n";
    Stiffness.Assemble(*matRed,rhsFull,guiInterface);
    std::cout << "Finish assembling matRed ...\n";
    clock.stop();
    // std::cout << "Time Assemble " << clock << std::endl;

    TPZBlockDiagonalStructMatrix<STATE> BDFmatrix(fAnalysis->Mesh());
    BDFmatrix.SetEquationRange(nEqLinr,nEqLinr+nEqHigh);
    TPZBlockDiagonal<REAL> KBD;
    
    std::cout << "Start assembling BlockDiag ...\n";
    BDFmatrix.AssembleBlockDiagonal(KBD);
    std::cout << "Finish assembling BlockDiag ...\n";
    // std::ofstream out3("out.txt");
    
    matRed->SetF(rhsFull);
    matRed->SetReduced();
    matRed->F1Red(rhsHigh);

    // rhsHigh.Print("RhsHigh=");

    //Creates the preconditioner 
    TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>( &KBD );
    precond->SetDirect(ELU);
    int64_t nMaxIter = 500;
    TPZVec<REAL> errors(nMaxIter);
    errors.Fill(0.);
    // for (int64_t iter = 1; iter < nMaxIter; iter++){
        
        // std::cout << "ITER = " << iter << std::endl;
        TPZFMatrix<STATE> residual(nMaxIter,1,0.);
        REAL tol = 1.e-10;
        TPZFMatrix<STATE> solution(nEqHigh,1,0.);
        clock.start();
        // K11Red->SolveCG(nMaxIter,*precond,rhsHigh,solution,&residual,tol);

        // std::ofstream out2("out3.txt");
        // matRed->Print("MATRED",out2,EMathematicaInput);
        // KBD.Print("BDiag",out2,EMathematicaInput);

        std::cout << "Start CG ...\n";
        matRed->SolveCG(nMaxIter,*precond,rhsHigh,solution,&residual,tol);
        std::cout << "Finish CG ...\n";
        
        clock.stop();
        // std::cout << "Time CG " << clock << std::endl;

        REAL norm = 0.;
        std::cout << "Number of CG iterations = " << nMaxIter << " , residual = " << tol << std::endl;
        out << nMaxIter << "\n";
        TPZFMatrix<STATE> result(nEqLinr+nEqHigh,1,0.);

        matRed->UGlobal2(solution,result);

        fAnalysis->Solution()=result;
        fAnalysis->LoadSolution();
        // rhsFull.Print("Solution = ");
                        
        ///Calculating approximation error  
        TPZManVector<REAL,5> error;

        // auto cmeshAux = an.Mesh();
        // int64_t nelem = cmeshAux->NElements();
        // cmeshAux->LoadSolution(cmeshAux->Solution());
        // cmeshAux->ExpandSolution();
        // an.Mesh()->ElementSolution().Redim(an.Mesh()->NElements(), 5);

        // an.PostProcessError(error);
            
        // std::cout << "\nApproximation error:\n";
        // std::cout << "H1 Norm = " << std::scientific << std::setprecision(15) << error[0]<<'\n';
        // std::cout << "L1 Norm = " << std::scientific << std::setprecision(15) << error[1]<<'\n'; 
        // std::cout << "H1 Seminorm = " << std::scientific << std::setprecision(15) << error[2]<<'\n'; 

        // if (tol < 1.e-10) break;
        // errors[iter-1] = error[1];
    // }

    // std::cout << "Number of iterations = " << residual << std::endl;

}


template class TPZMatRedSolver<STATE>;