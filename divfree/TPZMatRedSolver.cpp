#include "TPZMatRedSolver.h"
#include "pzcmesh.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "TPZTimer.h"
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "TPZEigenSolver.h"
#include "TPZLapackEigenSolver.h"
#include "pzspblockdiagpivot.h"
#include "pzintel.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZDoubleMatRed.h"

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
    case EMHMSparse:
        SolveProblemMHMSparse(out);
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

    std::ofstream myfile("MultiCMesh1.txt");
    fAnalysis->Mesh()->Print(myfile);

    //Cria a matriz esparsa
    TPZSparseMatRed<STATE> *matRed2 = new TPZSparseMatRed<STATE>(1,1);

    std::ofstream myfile2("MultiCMesh2.txt");
    fAnalysis->Mesh()->Print(myfile2);

    std::set<int> lag={1};
    matRed2->ReorderEquations(cmesh,lag,nEqFull,nEqLinr);

    std::ofstream myfile3("MultiCMesh3.txt");
    fAnalysis->Mesh()->Print(myfile3);

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

    std::ofstream myfile4("MultiCMesh4.txt");
    fAnalysis->Mesh()->Print(myfile4);

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
    
    // matRed->Print("MatRed = ",std::cout,EMathematicaInput);

    // TPZFMatrix<REAL> *K11Red = new TPZFMatrix<REAL>(nEqHigh,nEqHigh);

    matRed->SetF(rhsFull);
    matRed->SetReduced();
    // matRed->Print("MATRED",out2,EMathematicaInput);

    // matRed->K11Reduced(*K11Red,rhsHigh);
    
    matRed->F1Red(rhsHigh);
    // rhsFull.Zero();

    std::cout << "Assembling block diagonal " << std::endl;
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

    // ComputeConditionNumber(*matRed,precond->Matrix());

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
    step.SetDirect(ECholesky);//ELU //ECholesky // ELDLt
    fAnalysis->SetSolver(step);
    step.SetMatrix(&K00);

    //Cria a matriz esparsa
    std::set<int> lag = {1};
    TPZSparseMatRed<STATE> *matRed = new TPZSparseMatRed<STATE>(cmesh,lag);
    matRed->SetK00IsNegativeDefinite();
    std::cout << "Allocating Sub Matrices ...\n";
    auto start_time_allocating = std::chrono::steady_clock::now();
    //Transfere as submatrizes da matriz auxiliar para a matriz correta.
    matRed->SetSolver(&step);
    K00.Resize(matRed->Dim0(),matRed->Dim0());
    matRed->AllocateSubMatrices(cmesh);
    auto total_time_allocating = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_allocating).count()/1000.;
    std::cout << "Time Allocating Sub Matrices = " << total_time_allocating << " seconds" << std::endl;


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
//    constexpr int nThreads{50};
    constexpr int nThreads{0};
    
    // Create the RHS vectors
    TPZFMatrix<STATE> rhsFull(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsHigh(nEqHigh,1,0.);

    // ThresholdPermeability(5.e-3);

    //Creates the problem matrix    
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Stiffness(fAnalysis->Mesh());
    Stiffness.SetNumThreads(nThreads);
    

    TPZAutoPointer<TPZGuiInterface> guiInterface;
  
    Stiffness.EquationFilter().Reset();
    // Stiffness.EquationFilter().SetActiveEquations(fActiveEquations);
    fAnalysis->SetStructuralMatrix(Stiffness);

    //Monta a matriz auxiliar
    rhsFull.Zero();

    auto start_time_assemble = std::chrono::steady_clock::now();
    Stiffness.Assemble(*matRed,rhsFull,guiInterface);
    auto total_time_assemble = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_assemble).count()/1000.;
    std::cout << "Time Assembling SparseMatRed " << total_time_assemble << " seconds" << std::endl;

    //Block Diagonal
    auto start_time_bd = std::chrono::steady_clock::now();
    TPZBlockDiagonal<REAL> KBD;
    int ord = cmesh->GetDefaultOrder();
    if (ord == 0) ord = 1;
    TPZVec<int> blocksize(nEqHigh/ord,ord);

    KBD.Initialize(blocksize);    
    KBD.UpdateFrom(matRed->K11());
    
    auto total_time_bd = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_bd).count()/1000.;
    std::cout << "Time Assembling block diagonal " << total_time_bd << " seconds" << std::endl;
    
    
    // std::ofstream out2("out3.txt");
    // matRed->Print("MATRED",out2,EMathematicaInput);

    matRed->SetF(rhsFull);
    matRed->SetReduced();
    
    //Decomposes the reduced matrix
    auto start_time_decomp = std::chrono::steady_clock::now();
    matRed->F1Red(rhsHigh);
    auto total_time_decomp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_decomp).count()/1000.;
    std::cout << "Time decomposing k00 " << total_time_decomp << " seconds" << std::endl;

    // std::ofstream out3("out.txt");
    
    // matRed->Print("MatRed = ",out3,EMathematicaInput);

    //Creates the preconditioner 
    TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>( &KBD );
    precond->SetDirect(ELU);
    int64_t nMaxIter = 500;
    TPZVec<REAL> errors(nMaxIter);
    errors.Fill(0.);

    // ComputeConditionNumber(*matRed,precond->Matrix());
    
    // for (int64_t iter = 1; iter < nMaxIter; iter++){
        
        // std::cout << "ITER = " << iter << std::endl;
        TPZFMatrix<STATE> residual(nMaxIter,1,0.);
        REAL tol = 1.e-10;
        TPZFMatrix<STATE> solution(nEqHigh,1,0.);
        auto start_time_solve = std::chrono::steady_clock::now();
        // K11Red->SolveCG(nMaxIter,*precond,rhsHigh,solution,&residual,tol);

        // std::ofstream out2("out3.txt");
        // matRed->Print("MATRED",out2,EMathematicaInput);
        // KBD.Print("BDiag",out2,EMathematicaInput);

        // std::cout << "Start CG ...\n";
        matRed->SolveCG(nMaxIter,*precond,rhsHigh,solution,&residual,tol);
        // matRed->K11().SolveCG(nMaxIter,*precond,rhsHigh,solution,&residual,tol);
        // std::cout << "Finish CG ...\n";
        
        auto total_time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_solve).count()/1000.;
        std::cout << "Time CG " << total_time_solve << std::endl;

        REAL norm = 0.;
        std::cout << "Number of CG iterations = " << nMaxIter << " , residual = " << tol << std::endl;
        out << nMaxIter << "\n";
        TPZFMatrix<STATE> result(nEqLinr+nEqHigh,1,0.);

        matRed->UGlobal(solution,result);

        fAnalysis->Solution() = result;
        fAnalysis->LoadSolution();
        // rhsFull.Print("Solution = ");
        // std::cout << "Result = " << result << std::endl;
        // std::cout << "solution = " << solution << std::endl;
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
void TPZMatRedSolver<TVar>::SolveProblemMHMSparse(std::ostream &out){
    
    //HERE STARTS THE ITERATIVE SOLVER SET
    auto cmesh = fAnalysis->Mesh();

    //Primeiro cria a matriz auxiliar K00 - que será decomposta
    TPZSYsmpMatrix<REAL> K00;
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);//ELU //ECholesky // ELDLt
    fAnalysis->SetSolver(step);
    step.SetMatrix(&K00);

    //Cria a matriz esparsa
    std::set<int> lag = {1};
    TPZDoubleMatRed<STATE> *matRed = new TPZDoubleMatRed<STATE>(cmesh,lag);
    matRed->SetK00IsNegativeDefinite();
    std::cout << "Allocating Sub Matrices ...\n";
    auto start_time_allocating = std::chrono::steady_clock::now();
    //Transfere as submatrizes da matriz auxiliar para a matriz correta.
    matRed->SetSolver(&step);
    K00.Resize(matRed->Dim0(),matRed->Dim0());
    matRed->AllocateSubMatrices(cmesh);
    auto total_time_allocating = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_allocating).count()/1000.;
    std::cout << "Time Allocating Sub Matrices = " << total_time_allocating << " seconds" << std::endl;


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
    constexpr int nThreads{50};
    
    // Create the RHS vectors
    TPZFMatrix<STATE> rhsFull(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsHigh(nEqHigh,1,0.);

    // ThresholdPermeability(5.e-3);

    //Creates the problem matrix    
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Stiffness(fAnalysis->Mesh());
    Stiffness.SetNumThreads(nThreads);
    

    TPZAutoPointer<TPZGuiInterface> guiInterface;
  
    Stiffness.EquationFilter().Reset();
    // Stiffness.EquationFilter().SetActiveEquations(fActiveEquations);
    fAnalysis->SetStructuralMatrix(Stiffness);

    //Monta a matriz auxiliar
    rhsFull.Zero();

    auto start_time_assemble = std::chrono::steady_clock::now();
    Stiffness.Assemble(*matRed,rhsFull,guiInterface);
    auto total_time_assemble = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_assemble).count()/1000.;
    std::cout << "Time Assembling SparseMatRed " << total_time_assemble << " seconds" << std::endl;

    //Block Diagonal
    auto start_time_bd = std::chrono::steady_clock::now();
    TPZBlockDiagonal<REAL> KBD;
    int ord = cmesh->GetDefaultOrder();
    if (ord == 0) ord = 1;
    TPZVec<int> blocksize(nEqHigh/ord,ord);

    KBD.Initialize(blocksize);    
    KBD.UpdateFrom(matRed->K11());
    
    auto total_time_bd = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_bd).count()/1000.;
    std::cout << "Time Assembling block diagonal " << total_time_bd << " seconds" << std::endl;
    
    
    // std::ofstream out2("out3.txt");
    // matRed->Print("MATRED",out2,EMathematicaInput);

    matRed->SetF(rhsFull);
    matRed->SetReduced();
    
    //Decomposes the reduced matrix
    auto start_time_decomp = std::chrono::steady_clock::now();
    matRed->F1Red(rhsHigh);
    auto total_time_decomp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_decomp).count()/1000.;
    std::cout << "Time decomposing k00 " << total_time_decomp << " seconds" << std::endl;

    // std::ofstream out3("out.txt");
    
    // matRed->Print("MatRed = ",out3,EMathematicaInput);

    //Creates the preconditioner 
    TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>( &KBD );
    precond->SetDirect(ELU);
    int64_t nMaxIter = 500;
    TPZVec<REAL> errors(nMaxIter);
    errors.Fill(0.);

    // ComputeConditionNumber(*matRed,precond->Matrix());
    
    // for (int64_t iter = 1; iter < nMaxIter; iter++){
        
        // std::cout << "ITER = " << iter << std::endl;
        TPZFMatrix<STATE> residual(nMaxIter,1,0.);
        REAL tol = 1.e-10;
        TPZFMatrix<STATE> solution(nEqHigh,1,0.);
        auto start_time_solve = std::chrono::steady_clock::now();
        // K11Red->SolveCG(nMaxIter,*precond,rhsHigh,solution,&residual,tol);

        // std::ofstream out2("out3.txt");
        // matRed->Print("MATRED",out2,EMathematicaInput);
        // KBD.Print("BDiag",out2,EMathematicaInput);

        // std::cout << "Start CG ...\n";
        matRed->SolveCG(nMaxIter,*precond,rhsHigh,solution,&residual,tol);
        // matRed->K11().SolveCG(nMaxIter,*precond,rhsHigh,solution,&residual,tol);
        // std::cout << "Finish CG ...\n";
        
        auto total_time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_solve).count()/1000.;
        std::cout << "Time CG " << total_time_solve << std::endl;

        REAL norm = 0.;
        std::cout << "Number of CG iterations = " << nMaxIter << " , residual = " << tol << std::endl;
        out << nMaxIter << "\n";
        TPZFMatrix<STATE> result(nEqLinr+nEqHigh,1,0.);

        matRed->UGlobal(solution,result);

        fAnalysis->Solution() = result;
        fAnalysis->LoadSolution();
        // rhsFull.Print("Solution = ");
        // std::cout << "Result = " << result << std::endl;
        // std::cout << "solution = " << solution << std::endl;
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
void TPZMatRedSolver<TVar>::ComputeConditionNumber(TPZSparseMatRed<STATE> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond){
    
    TPZFMatrix<REAL> KBDInv;
    TPZAutoPointer<TPZFMatrix<REAL>> Res = new TPZFMatrix<REAL>;
    auto dim = precond->Rows();
    Res->Redim(dim,dim);
    precond->Inverse(KBDInv,ELU);
    // KBDInv.Identity();
    // KBDInv.Print("KBDInv=",std::cout,EMathematicaInput);
    // matRed.K11().Print("K11=",std::cout,EMathematicaInput);
    // matRed.MultAdd(KBDInv,*Res,*Res,1.,0.);
    matRed.K11().MultAdd(KBDInv,*Res,*Res,1.,0.);
    // KBDInv.Multiply(*matRed,Res);

    
    
    // Res->Print("Res=",std::cout,EMathematicaInput);

    TPZLapackEigenSolver<REAL> eigSolver;
    
    TPZVec<std::complex<REAL>> eigenvalues;
    eigSolver.SetMatrixA(Res);
    auto a1 = eigSolver.SolveEigenProblem(eigenvalues);

    std::ofstream rprint3,rprint4;
    rprint3.open("REAL_EIGEN_ALL.txt",std::ios_base::app);
    rprint4.open("REAL_EIGEN.txt",std::ios_base::app);

    REAL maxEig = 0.;
    REAL minEig = 1e3;
    REAL minAbs = 1e3;
    REAL maxAbs = 0.;
    int nonzeroEigenvalues = 0;
    REAL tol = 1e-10;
    for (int i = 0; i < eigenvalues.size(); i++)
    {
        rprint3 << eigenvalues[i].real() << std::endl;
        if (eigenvalues[i].real() > maxEig) maxEig = eigenvalues[i].real();
        if (eigenvalues[i].real() < minEig) minEig = eigenvalues[i].real();
        if (fabs(eigenvalues[i].real()) < minAbs) minAbs = fabs(eigenvalues[i].real());
        if (fabs(eigenvalues[i].real()) > maxAbs) maxAbs = fabs(eigenvalues[i].real());
        if (fabs(eigenvalues[i].real()) > tol) nonzeroEigenvalues++;
    }
    rprint3 << std::endl;

    rprint4 << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim<< std::endl;
        std::cout << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim<< std::endl;


}

template<class TVar>
void TPZMatRedSolver<TVar>::ComputeConditionNumber(TPZMatRed<STATE,TPZFMatrix<STATE>> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond){
    
    TPZFMatrix<REAL> KBDInv;
    TPZAutoPointer<TPZFMatrix<REAL>> Res = new TPZFMatrix<REAL>;
    auto dim = precond->Rows();
    // Res->(dim,dim,true);
    precond->Inverse(KBDInv,ELU);
    // KBDInv.Identity();
    KBDInv.Print("KBDInv=",std::cout,EMathematicaInput);
    matRed.Print("MatRed=",std::cout,EMathematicaInput);

    matRed.Multiply(KBDInv,Res);
    // KBDInv.Multiply(*matRed,Res);

    
    
    Res->Print("Res=",std::cout,EMathematicaInput);

    TPZLapackEigenSolver<REAL> eigSolver;
    
    TPZVec<std::complex<REAL>> eigenvalues;
    eigSolver.SetMatrixA(Res);
    
    // auto a1 = eigSolver.SolveEigenProblem(eigenvalues);

    std::ofstream rprint3,rprint4;
    rprint3.open("REAL_EIGEN_ALL.txt",std::ios_base::app);
    rprint4.open("REAL_EIGEN.txt",std::ios_base::app);

    REAL maxEig = 0.;
    REAL minEig = 1e3;
    REAL minAbs = 1e3;
    REAL maxAbs = 0.;
    int nonzeroEigenvalues = 0;
    REAL tol = 1e-10;
    for (int i = 0; i < eigenvalues.size(); i++)
    {
        rprint3 << eigenvalues[i].real() << std::endl;
        if (eigenvalues[i].real() > maxEig) maxEig = eigenvalues[i].real();
        if (eigenvalues[i].real() < minEig) minEig = eigenvalues[i].real();
        if (fabs(eigenvalues[i].real()) < minAbs) minAbs = fabs(eigenvalues[i].real());
        if (fabs(eigenvalues[i].real()) > maxAbs) maxAbs = fabs(eigenvalues[i].real());
        if (fabs(eigenvalues[i].real()) > tol) nonzeroEigenvalues++;
    }
    rprint3 << std::endl;

    rprint4 << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim<< std::endl;
        std::cout << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim<< std::endl;


}



template<class TVar>
void TPZMatRedSolver<TVar>::ThresholdPermeability(REAL threshold){

    auto cmeshanalysis = fAnalysis->Mesh();
    TPZMultiphysicsCompMesh *mCmesh = dynamic_cast<TPZMultiphysicsCompMesh *> (cmeshanalysis);
    auto cmesh = mCmesh->MeshVector()[0];
    const int64_t nEquations = fAnalysis->Mesh()->NEquations();

    
    std::set<int64_t> remove_eq;

    for (TPZCompEl* cel : cmesh->ElementVec()){
        if (!cel) continue;
        int dim = cel->Dimension();
        if (dim != cmesh->Dimension()) continue;

        TPZGeoEl* gel = cel->Reference();
        if (!gel) DebugStop();

        int nFaces = gel->NSides(dim-1);
        int nSides = gel->NSides();

        for (int iFace = 0; iFace < nFaces; iFace++){
            // Pega o geoElSide de cada face, depois pega CenterX.
            
            TPZGeoElSide gelside(gel,nSides-1-nFaces+iFace);
            TPZManVector<REAL,3> xCenter(3,0.), perm(3,0.);
            gelside.CenterX(xCenter);

            perm = fPermFunction(xCenter);

            if (perm[0] < threshold){
                
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                if (!intel) DebugStop();
                int64_t conindex1 = cel->ConnectIndex(2*iFace);
                int64_t conindex2 = cel->ConnectIndex(2*iFace+1);

                TPZConnect &c1 = cmesh->ConnectVec()[conindex1];
                TPZConnect &c2 = cmesh->ConnectVec()[conindex2];

                // std::cout << "Dim = " << dim << std::endl;

                int64_t seqnum = c1.SequenceNumber();
                const auto pos = cmesh->Block().Position(seqnum);
                const auto blocksize = cmesh->Block().Size(seqnum);
                // if (blocksize == 0){
                //     continue;
                // }
                // const auto vs = fActiveEquations.size();
                if (pos < nEquations){
                    // fActiveEquations.Resize(vs + blocksize);
                    for (auto ieq = 0; ieq < blocksize; ieq++) {
                        remove_eq.insert(pos + ieq);
                    }
                }
                
                // std::cout << "Connect index = " << conindex1 << " " << blocksize << std::endl;

                // int64_t seqnum2 = c2.SequenceNumber();
                // const auto pos2 = cmesh->Block().Position(seqnum2);
                // const auto blocksize2 = cmesh->Block().Size(seqnum2);
                // if (blocksize == 0){
                //     continue;
                // }
                // if (c2.IsCondensed()){
                //     std::cout << "Is Condensed " << seqnum2 << std::endl;
                //     continue;
                // }
                // // std::cout << "Connect index = " << conindex2 << " " << blocksize2 << std::endl;
                // const auto vs2 = fActiveEquations.size();
                // if (pos2 < nEquations){
                //     fActiveEquations.Resize(vs2 + blocksize2);
                //     for (auto ieq = 0; ieq < blocksize2; ieq++) {
                //         fActiveEquations[vs2 + ieq] = pos2 + ieq;
                //     }
                // }
                
                
                // int64_t index;
                // //Sets the new connect order
                // c2.SetOrder(0,index);//Div constant function
                // c1.SetOrder(0,index);//High order div-free functions

                // //Gets connect information to update block size (stiffness matrix data structure)
                // int64_t seqnum = c2.SequenceNumber();
                // int nvar = 1;
                // TPZMaterial * mat = cel->Material();
                // if (mat) nvar = mat->NStateVariables();
                // int nshape = intel->NConnectShapeF(2*iFace,c1.Order());
                // int nshape2 = intel->NConnectShapeF(2*iFace+1,c2.Order());
                // c2.SetNShape(nshape2);
                // // c.SetNState(nvar);
                // cmesh->Block().Set(seqnum, nvar * nshape2);


                // seqnum = c1.SequenceNumber();
                // // nshape = 1;
                // cmesh->Block().Set(seqnum, nvar * nshape);

            }
        }
    }
    cmesh->InitializeBlock();
    // cmesh->ExpandSolution();

    std::cout << "remove_eq = " ;
    for (auto it:remove_eq)
    {
        std::cout << it << " " ;
    }
    
    fActiveEquations.Resize(nEquations - remove_eq.size());
    std::cout << "NEquations = " << nEquations << std::endl;
    std::cout << "remove_eq = " << remove_eq.size() << std::endl;
    std::cout << "fActiveEquations = " << fActiveEquations.size() << std::endl;

    int64_t count = 0;
    for (int i = 0; i < nEquations; i++)
    {
        if (remove_eq.find(i) == remove_eq.end()) {
            fActiveEquations[count] = i;
            std::cout << "Active equation [" << count << "] = " << fActiveEquations[count] << std::endl;
            count++;
        }
    }

}



template class TPZMatRedSolver<STATE>;
