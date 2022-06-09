#include "TPZMatRedSolver.h"
#include "pzcmesh.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "TPZTimer.h"
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"

template<class TVar>
void TPZMatRedSolver<TVar>::ReorderEquations(int64_t &nEqLinr, int64_t &nEqHigh){
    
    auto cmesh = fAnalysis->Mesh();
    cmesh->LoadReferences();
    std::set<int64_t> auxConnectsP, auxConnectsF;
    // TPZSkylineStructMatrix<STATE> Pmatrix(cmesh);
    // TPZSkylineStructMatrix<STATE> Fmatrix(cmesh);
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Pmatrix(cmesh);
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Fmatrix(cmesh);
    auto neqsTotal = cmesh->NEquations();

    auto gmesh = cmesh->Reference();
    nEqLinr = 0;
    nEqHigh = 0;

    // PrintCompMesh(cmesh,"CMESH_0");

    // loop over the connects and create two std::sets one to the "pressure" ones, representing the 
    // degrees of freedom to be condensed. The second set contains the "flux" connects, which will not be condensed
    // This can change depending on the problem.
    for (auto gel:gmesh->ElementVec()){
        
        if (!gel) continue;
        //Looking for the hybridized pressure
        if (gel->Dimension() == gmesh->Dimension()-1 && gel->MaterialId() == 6){
            auxConnectsP.insert(gel->Reference()->ConnectIndex(0));
        }

        if (gel->Dimension() != gmesh->Dimension()) continue;//Only looks for the volumetric elements
        auto nconnects = gel->Reference()->NConnects();
        int nFacets = gel->NSides(gmesh->Dimension()-1);

        for (size_t i = 0; i < nFacets; i++)
        {
            //High order edge functions will not be condensed
            auxConnectsF.insert(gel->Reference()->ConnectIndex(2*i+1));
            //Lower order edge functions will be condensed
            auxConnectsP.insert(gel->Reference()->ConnectIndex(2*i  ));
        }
        //The internal connect will always be condensed
        auxConnectsP.insert(gel->Reference()->ConnectIndex(2*nFacets));
        //Pressure degrees of freedom will not be condensed
        for (size_t i = 2*nFacets+1; i < nconnects; i++){
            auxConnectsP.insert(gel->Reference()->ConnectIndex(i));
        }
    }
    // // Now loop over the boundary elements to condense it and correct the map structure
    // for (auto gel:gmesh->ElementVec()){
        
    //     if (gel->Dimension() == gmesh->Dimension()) continue;//Only looks for the BC elements
    //     auto nconnects = gel->Reference()->NConnects();
    //     int matId = gel->MaterialId();

    //     //Loops over all the BC element connects
    //     for (size_t icon = 0; icon < nconnects; icon++)
    //     {
    //         //Sets all connects to be condensed and remove the entry from auxConnectsF.
    //         auxConnectsP.insert(gel->Reference()->ConnectIndex(icon));
    //         auxConnectsF.erase(gel->Reference()->ConnectIndex(icon));           
    //     }
        
    // }

    // std::cout << "AuxConnects P = " << auxConnectsP << std::endl;
    // std::cout << "AuxConnects F = " << auxConnectsF << std::endl;

    //If the previous structure was properly filled, there is no need to change from here.
    //First - set the sequence number for the "pressure" variables, i.e., the variables to be condensed 
    int64_t seqNumP = 0;
    cmesh->Block().Resequence();

    for (int icon = 0; icon < cmesh->NConnects(); icon++){

        TPZConnect &con = cmesh->ConnectVec()[icon];
        if (auxConnectsP.find(icon) != auxConnectsP.end()) {
            int64_t seqNum = con.SequenceNumber();
            if (con.IsCondensed()) continue;
            if (seqNum < 0) continue;

            con.SetSequenceNumber(seqNumP);
            
            //Em cada caso precisa atualizar o tamanho do bloco
            int neq = con.NDof()*con.NState();
            seqNumP++;
            nEqLinr += neq;
            seqNum=con.SequenceNumber();
            cmesh->Block().Set(seqNum,neq);
        }
    }

    int64_t seqNumF = seqNumP;

    //Second - Set the sequence number to the flux variables - which will not be condensed 
    for (int icon = 0; icon < cmesh->NConnects(); icon++){

        TPZConnect &con = cmesh->ConnectVec()[icon];
        if (auxConnectsF.find(icon) != auxConnectsF.end()) {
            int64_t seqNum = con.SequenceNumber();
            if (con.IsCondensed()) continue;
            if (seqNum < 0) continue;

            con.SetSequenceNumber(seqNumF);
            
            //Em cada caso precisa atualizar o tamanho do bloco
            int neq = con.NDof()*con.NState();
            seqNumF++;
            nEqHigh += neq;
            seqNum=con.SequenceNumber();
            cmesh->Block().Set(seqNum,neq);
        }
    }   
    cmesh->ExpandSolution();

}


template<class TVar>
void TPZMatRedSolver<TVar>::Solve(std::ostream &out){
    
    //HERE STARTS THE ITERATIVE SOLVER SET
    //Sets number of threads to be used by the solver
    constexpr int nThreads{10};
    auto cmesh = fAnalysis->Mesh();

    //Compute the number of equations in the system
    int64_t nEqFull = cmesh->NEquations();
    int64_t nEqLinr, nEqHigh;

    ReorderEquations(nEqLinr,nEqHigh);

    out << nEqHigh << " " << nEqLinr << " ";

    std::cout << "NUMBER OF EQUATIONS:\n " << 
                 "Full problem = " << nEqFull << 
                 ", High Order Flux = " << nEqHigh << 
                 ", Linear Flux = " << nEqLinr << std::endl;

    // Create the RHS vectors
    TPZFMatrix<STATE> rhsFull(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsHigh(nEqHigh,1,0.);

    //Creates the problem matrix    
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Stiffness(cmesh);
    Stiffness.SetNumThreads(nThreads);

    TPZMatrix<TVar> *K11Red, *matRed;

    switch (fSolverType)
    {
    case EDefault:
        SolveProblemDefault(nEqLinr,nEqHigh,Stiffness,rhsFull,rhsHigh,out);
        break;
    case ESparse:
        SolveProblemSparse(nEqLinr,nEqHigh,Stiffness,rhsFull,rhsHigh,out);
        break;
    
    default:
        DebugStop();
        break;
    }

}


template<class TVar>
void TPZMatRedSolver<TVar>::SolveProblemDefault(int64_t &nEqLinr, int64_t &nEqHigh, TPZStructMatrix &Stiffness, TPZFMatrix<TVar> &rhsFull, TPZFMatrix<TVar> &rhsHigh, std::ostream &out){

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
    
    // std::ofstream out2("out2.txt");

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
void TPZMatRedSolver<TVar>::SolveProblemSparse(int64_t &nEqLinr, int64_t &nEqHigh, TPZStructMatrix &Stiffness, TPZFMatrix<TVar> &rhsFull, TPZFMatrix<TVar> &rhsHigh, std::ostream &out){
    
    TPZAutoPointer<TPZGuiInterface> guiInterface;

    //Cria a matriz esparsa
    TPZSparseMatRed<STATE> *matRed = new TPZSparseMatRed<STATE>(nEqLinr+nEqHigh,nEqLinr);



    //Primeiro cria a matriz auxiliar
    TPZSYsmpMatrix<REAL> K00(nEqLinr,nEqLinr);
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    fAnalysis->SetSolver(step);
    
    std::cout << "Allocating Sub Matrices ...\n";
    AllocateSubMatrices(nEqLinr,nEqHigh,Stiffness,K00,matRed);
    
    //Transfere as submatrizes da matriz auxiliar para a matriz correta.
    step.SetMatrix(&K00);
    matRed->SetSolver(&step);
  
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
    
    // TPZFMatrix<REAL> *K11Red = new TPZFMatrix<REAL>(nEqHigh,nEqHigh,0.);

    matRed->SetF(rhsFull);
    matRed->SetReduced();
    // matRed->Print("MATRED",out3,EMathematicaInput);
    // matRed->K11Reduced(*K11Red,rhsHigh);
    matRed->F1Red(rhsHigh);

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
        std::cout << "Start CG ...\n";
        matRed->SolveCG(nMaxIter,*precond,rhsHigh,solution,&residual,tol);
        std::cout << "Finish CG ...\n";
        
        clock.stop();
        // std::cout << "Time CG " << clock << std::endl;

        REAL norm = 0.;
        std::cout << "Number of CG iterations = " << nMaxIter << " , residual = " << tol << std::endl;
        out << nMaxIter << "\n";
        TPZFMatrix<STATE> result(nEqLinr+nEqHigh,1,0.);

        matRed->UGlobal(solution,result);

        // //Update solution in Analysis
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
void TPZMatRedSolver<TVar>::AllocateSubMatrices(int64_t &nEqLinr, int64_t &nEqHigh, TPZStructMatrix &Stiffness, TPZSYsmpMatrix<TVar> &K00, TPZSparseMatRed<TVar> *matRed) {
    
    //Aloca as submatrizes no formato esparso.
    Stiffness.EquationFilter().Reset();
    TPZSYsmpMatrix<REAL> *StiffK11 = dynamic_cast<TPZSYsmpMatrix<REAL> *>(Stiffness.Create());

    //Fazer uma rotina para separar IA, JA e A de K01, K10 e K11;
    TPZVec<int64_t> IA_K00(nEqLinr+1,0), IA_K01(nEqLinr+1,0), IA_K10(nEqHigh+1,0), IA_K11(nEqHigh+1,0);
    
    std::vector<int64_t> auxK00, auxK01, auxK11;
    auxK00.reserve(StiffK11->JA().size());
    auxK01.reserve(StiffK11->JA().size());
    auxK11.reserve(StiffK11->JA().size());
    
    for (int i = 0; i < nEqLinr; i++){
        for (int j = StiffK11->IA()[i]; j < StiffK11->IA()[i+1]; j++){
            if (StiffK11->JA()[j] < nEqLinr){
                // Faz parte da matriz K00
                auxK00.push_back(StiffK11->JA()[j]);
            } else {
                // Faz parte da matriz K01
                auxK01.push_back(StiffK11->JA()[j]-nEqLinr);
            }
        }
        IA_K00[i+1] = auxK00.size();
        IA_K01[i+1] = auxK01.size();
    }
    for (int i = nEqLinr; i < nEqLinr+nEqHigh; i++){
        for (int j = StiffK11->IA()[i]; j < StiffK11->IA()[i+1]; j++){
            if (StiffK11->JA()[j] >= nEqLinr){
                // Faz parte da matriz K11
                auxK11.push_back(StiffK11->JA()[j]-nEqLinr);
            }
        }
        IA_K11[i-nEqLinr+1] = auxK11.size();
    }
    
    //Do the transpose - Matriz K10
    IA_K10[0]=0;
    for (int i = 0  ; i < nEqHigh; i++){
        int nNonZeros = std::count(auxK01.begin(),auxK01.end(),i);
        IA_K10[i+1] = IA_K10[i] + nNonZeros;
    }

    // Resize the CRS structure with the correct size
    TPZVec<int64_t> JA_K00(auxK00.size(),0), JA_K01(auxK01.size(),0), JA_K10(auxK01.size(),0), JA_K11(auxK11.size(),0);
    TPZVec<double> A_K00(auxK00.size(),0.), A_K01(auxK01.size(),0.), A_K10(auxK01.size(),0.), A_K11(auxK11.size(),0.);

    // Sets values to the nonzero values columns entries
    for (int i = 0; i < JA_K00.size(); i++) JA_K00[i] = auxK00[i];
    for (int i = 0; i < JA_K01.size(); i++) JA_K01[i] = auxK01[i];
    // K10 is skiped because the transposition is performed inside TPZSparseMatRed, so here we insert a zero vector.
    for (int i = 0; i < JA_K11.size(); i++) JA_K11[i] = auxK11[i];
    
    //Aloca estrutura das matrizes esparsas
    K00.SetData(IA_K00,JA_K00,A_K00);
    matRed->K01().SetData(IA_K01,JA_K01,A_K01);
    matRed->K01().SetData(IA_K01,JA_K01,A_K01);
    matRed->K10().SetData(IA_K10,JA_K10,A_K10);
    matRed->K11().SetData(IA_K11,JA_K11,A_K11);
}

template class TPZMatRedSolver<STATE>;