#include "TPZMatRedSolver.h"
#include "pzcmesh.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzstepsolver.h>       //for TPZStepSolver
#include "TPZTimer.h"
#include "pzblockdiag.h"
#include "pzbdstrmatrix.h"
#include "TPZEigenSolver.h"
#include "TPZLapackEigenSolver.h"
#include "pzspblockdiagpivot.h"
#include "pzintel.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZDoubleMatRed.h"
#include "TPZGuiInterface.h"
#ifdef PZ_USING_MKL
#include "TPZSYSMPPardiso.h"
#endif
#ifdef PZ_USING_MUMPS
#include "TPZSSpStructMatrixMumps.h"
#include "TPZSYSMPMatrix.h"
#include "TPZMumpsSolver.h"
#include "TPZSYSMPMumps.h"
#endif

template <class TVar>
void TPZMatRedSolver<TVar>::Solve(std::ostream &out)
{

  // HERE STARTS THE ITERATIVE SOLVER SET
  auto cmesh = fAnalysis->Mesh();
  int dimension = cmesh->Dimension();
// Primeiro cria a matriz auxiliar K00 - que será decomposta
//     TPZSYsmpMatrix<REAL> K00;
#ifdef PZ_USING_MKL
  TPZSYsmpMatrixPardiso<REAL> K00;
#else
  TPZSYsmpMatrixMumps<REAL> K00;
#endif

  TPZStepSolver<STATE> step;
  K00.SetSymmetry(SymProp::Sym);
  // step.SetDirect(ELU);//ELU //ECholesky // ELDLt
  step.SetDirect(ECholesky); // ELU //ECholesky // ELDLt
  // step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt

  fAnalysis->SetSolver(step);
  step.SetMatrix(&K00);

  // Cria a matriz esparsa
  std::set<int> lag = {1};
  TPZSparseMatRed<STATE> *matRed = new TPZSparseMatRed<STATE>(cmesh, lag);
  matRed->SetK00IsNegativeDefinite();
  std::cout << "Allocating Sub Matrices ...\n";
  auto start_time_allocating = std::chrono::steady_clock::now();
  // Transfere as submatrizes da matriz auxiliar para a matriz correta.
  matRed->SetSolver(&step);
#ifdef PZ_USING_MKL
  K00 = TPZSYsmpMatrixPardiso<REAL>(matRed->Dim0(), matRed->Dim0());
#endif
#ifdef PZ_USING_MUMPS
  K00 = TPZSYsmpMatrixMumps<REAL>(matRed->Dim0(), matRed->Dim0());
#endif
  // K00.Resize(matRed->Dim0(),matRed->Dim0());
  matRed->AllocateSubMatrices(cmesh);
  auto total_time_allocating = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_allocating).count() / 1000.;
  std::cout << "Time Allocating Sub Matrices = " << total_time_allocating << " seconds" << std::endl;

  // Compute the number of equations in the system
  int64_t nEqFull = cmesh->NEquations();
  int64_t nEqLinr = matRed->Dim0();
  int64_t nEqHigh = matRed->Dim1();

  // out << nEqHigh << " " << nEqLinr << " ";

  std::cout << "NUMBER OF EQUATIONS:\n " << "Full problem = " << nEqFull << ", High Order Flux = " << nEqHigh << ", Linear Flux = " << nEqLinr << std::endl;

  // Sets number of threads to be used by the solver
  constexpr int nThreads{32};
    // constexpr int nThreads{0};

  // Create the RHS vectors
  TPZFMatrix<STATE> rhsFull(nEqFull, 1, 0.);
  TPZFMatrix<STATE> rhsHigh(nEqHigh, 1, 0.);

// ThresholdPermeability(5.e-3);

// Creates the problem matrix
#ifdef PZ_USING_MKL
  TPZSSpStructMatrix<STATE, TPZStructMatrixOR<STATE>> Stiffness(fAnalysis->Mesh());
#endif
#ifdef PZ_USING_MUMPS
  TPZSSpStructMatrixMumps<STATE, TPZStructMatrixOR<STATE>> Stiffness(fAnalysis->Mesh());
#endif
  Stiffness.SetNumThreads(nThreads);

  fAnalysis->SetStructuralMatrix(Stiffness);

  // Monta a matriz
  rhsFull.Zero();
  auto start_time_assemble = std::chrono::steady_clock::now();
  Stiffness.Assemble(*matRed, rhsFull);
  if (this->fProblemOrigin == ProblemOrigin::EElasticityH1Hybrid || this->fProblemOrigin == ProblemOrigin::EDarcyH1Hybrid)
  {
    matRed->K11() *= (-1.);
    matRed->K10() *= (-1.);
    matRed->K01() *= (-1.);
    matRed->SetK00IsUpdated(true);
    rhsFull *= -1.;
  }
  auto total_time_assemble = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_assemble).count();
  std::cout << "Time Assembling SparseMatRed " << total_time_assemble << " ms" << std::endl;

  // Block Diagonal
  auto start_time_bd = std::chrono::steady_clock::now();
  TPZBlockDiagonal<REAL> KBD;
  int ord = cmesh->GetDefaultOrder();
  int nstate = -1;
  if (ord == 0)
  {
    DebugStop();
  }
  int bsize = 0;
  if (dimension == 2)
  {
    if (fProblemOrigin == ProblemOrigin::EDarcyHDiv || fProblemOrigin == ProblemOrigin::EDarcyH1Hybrid)
    {
      bsize = ord; // the size of the block is ord+1 and we subtract 1 because we are not considering the lagrange multiplier
      nstate = 1;
    }
    else if (fProblemOrigin == ProblemOrigin::EElasticityH1Hybrid || fProblemOrigin == ProblemOrigin::EElasticityHDiv)
    {
      bsize = 2 * (ord + 1) - 3;
      nstate = 2;
    }
    else
    {
      DebugStop();
    }
  }
  else if (dimension == 3 && fProblemOrigin == ProblemOrigin::EDarcyHDiv)
  {
    bsize = (ord + 1) * (ord + 1) - 1;
  }
  else
  {
    DebugStop();
  }
  TPZVec<int> blocksize(nEqHigh / bsize, bsize);

  KBD.Initialize(blocksize);
  KBD.UpdateFrom(matRed->K11());

  //   KBD.Print("KBD",std::cout,EMathematicaInput);

  auto total_time_bd = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_bd).count();
  std::cout << "Time Assembling block diagonal " << total_time_bd << " ms" << std::endl;

  matRed->SetF(rhsFull);
  matRed->SetReduced();

  // Decomposes the reduced matrix
  auto start_time_decomp = std::chrono::steady_clock::now();
  matRed->F1Red(rhsHigh);
  auto total_time_decomp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_decomp).count();
  std::cout << "Time decomposing k00 " << total_time_decomp << " ms" << std::endl;

  // Creates the preconditioner
  TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>(&KBD);
  precond->SetDirect(ELU);
  int64_t nMaxIter = 500;
  TPZVec<REAL> errors(nMaxIter);
  errors.Fill(0.);

  //   ComputeConditionNumber(*matRed,precond->Matrix());
  //   return;
  // for (int64_t iter = 1; iter < nMaxIter; iter++){

  // std::cout << "ITER = " << iter << std::endl;
  TPZFMatrix<STATE> residual(nMaxIter, 1, 0.);
  REAL tol = 1.e-10;
  TPZFMatrix<STATE> solution(nEqHigh, 1, 0.);
  auto start_time_solve = std::chrono::steady_clock::now();

  // std::ofstream out2("out3.txt");
  // matRed->Print("MATRED",out2,EMathematicaInput);
  //   KBD.Print("BDiag",std::cout,EMathematicaInput);

  // std::cout << "Start CG ...\n";
  //   precond->Matrix()->Print("Precond",std::cout,EMathematicaInput);
  //   rhsHigh.Print("rhsHigh",std::cout,EMathematicaInput);
  //   solution.Print("solution",std::cout,EMathematicaInput);
  //   residual.Print("residual",std::cout,EMathematicaInput);
  matRed->SolveCG(nMaxIter, *precond, rhsHigh, solution, &residual, tol);
  // std::cout << "Finish CG ...\n";

  auto total_time_solve = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time_solve).count();
  std::cout << "Time CG " << total_time_solve << std::endl;

  out << total_time_assemble << " " << total_time_decomp + total_time_bd + total_time_solve << " ";

  REAL norm = 0.;
  std::cout << "Number of CG iterations = " << nMaxIter << " , residual = " << tol << std::endl;
  // out << nMaxIter << "\n";
  TPZFMatrix<STATE> result(nEqLinr + nEqHigh, 1, 0.);

  matRed->UGlobal(solution, result);

  fAnalysis->Solution() = result;
  fAnalysis->LoadSolution();

  /// Calculating approximation error
  TPZManVector<REAL, 5> error;

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

template <class TVar>
void TPZMatRedSolver<TVar>::ComputeConditionNumber(TPZSparseMatRed<STATE> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond)
{

  TPZFMatrix<REAL> KBDInv;
  TPZAutoPointer<TPZFMatrix<REAL>> Res = new TPZFMatrix<REAL>;
  auto dim = precond->Rows();
  Res->Redim(dim, dim);
  precond->Inverse(KBDInv, ELU);
  //   KBDInv.Print("KBD=",std::cout,EMathematicaInput);
  // KBDInv.Identity();
  // KBDInv.Print("KBDInv=",std::cout,EMathematicaInput);
  // matRed.K11().Print("K11=",std::cout,EMathematicaInput);
  // matRed.MultAdd(KBDInv,*Res,*Res,1.,0.);
  TPZFMatrix<REAL> F1, Aux(dim, dim, 0);
  matRed.K11Reduced(Aux, F1);
  Aux.Multiply(KBDInv, Res);

  //   Res->Print("Res=",std::cout,EMathematicaInput);

  TPZLapackEigenSolver<REAL> eigSolver;

  TPZVec<std::complex<REAL>> eigenvalues;
  eigSolver.SetMatrixA(Res);
  auto a1 = eigSolver.SolveEigenProblem(eigenvalues);

  std::ofstream rprint3, rprint4;
  rprint3.open("REAL_EIGEN_ALL.txt", std::ios_base::app);
  rprint4.open("REAL_EIGEN.txt", std::ios_base::app);

  REAL maxEig = 0.;
  REAL minEig = 1e3;
  REAL minAbs = 1e3;
  REAL maxAbs = 0.;
  int nonzeroEigenvalues = 0;
  REAL tol = 1e-8;
  for (int i = 0; i < eigenvalues.size(); i++)
  {
    rprint3 << eigenvalues[i].real() << std::endl;
    if (eigenvalues[i].real() > maxEig)
      maxEig = eigenvalues[i].real();
    if (eigenvalues[i].real() < minEig)
      minEig = eigenvalues[i].real();
    if (fabs(eigenvalues[i].real()) < minAbs)
      minAbs = fabs(eigenvalues[i].real());
    if (fabs(eigenvalues[i].real()) > maxAbs)
      maxAbs = fabs(eigenvalues[i].real());
    if (fabs(eigenvalues[i].real()) > tol)
      nonzeroEigenvalues++;
  }
  rprint3 << std::endl;

  rprint4 << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim << std::endl;
  std::cout << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim << std::endl;
}

template <class TVar>
void TPZMatRedSolver<TVar>::ComputeConditionNumber(TPZMatRed<STATE, TPZFMatrix<STATE>> &matRed, TPZAutoPointer<TPZMatrix<REAL>> precond)
{

  TPZFMatrix<REAL> KBDInv;
  TPZAutoPointer<TPZFMatrix<REAL>> Res = new TPZFMatrix<REAL>;
  auto dim = precond->Rows();
  // Res->(dim,dim,true);
  precond->Inverse(KBDInv, ELU);
  // KBDInv.Identity();
  //   KBDInv.Print("KBDInv=",std::cout,EMathematicaInput);
  //   matRed.Print("MatRed=",std::cout,EMathematicaInput);

  matRed.Multiply(KBDInv, Res);
  // KBDInv.Multiply(*matRed,Res);

  //   Res->Print("Res=",std::cout,EMathematicaInput);

  TPZLapackEigenSolver<REAL> eigSolver;

  TPZVec<std::complex<REAL>> eigenvalues;
  eigSolver.SetMatrixA(Res);

  // auto a1 = eigSolver.SolveEigenProblem(eigenvalues);

  std::ofstream rprint3, rprint4;
  rprint3.open("REAL_EIGEN_ALL.txt", std::ios_base::app);
  rprint4.open("REAL_EIGEN.txt", std::ios_base::app);

  REAL maxEig = 0.;
  REAL minEig = 1e3;
  REAL minAbs = 1e3;
  REAL maxAbs = 0.;
  int nonzeroEigenvalues = 0;
  REAL tol = 1e-10;
  for (int i = 0; i < eigenvalues.size(); i++)
  {
    rprint3 << eigenvalues[i].real() << std::endl;
    if (eigenvalues[i].real() > maxEig)
      maxEig = eigenvalues[i].real();
    if (eigenvalues[i].real() < minEig)
      minEig = eigenvalues[i].real();
    if (fabs(eigenvalues[i].real()) < minAbs)
      minAbs = fabs(eigenvalues[i].real());
    if (fabs(eigenvalues[i].real()) > maxAbs)
      maxAbs = fabs(eigenvalues[i].real());
    if (fabs(eigenvalues[i].real()) > tol)
      nonzeroEigenvalues++;
  }
  rprint3 << std::endl;

  rprint4 << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim << std::endl;
  std::cout << maxEig << " " << minEig << " " << maxAbs << " " << minAbs << " " << nonzeroEigenvalues << "/" << dim << std::endl;
}

template class TPZMatRedSolver<STATE>;
