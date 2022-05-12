//
// Created by Jeferson Fernandes on 11/08/21.
//

#include "TPZKernelHdivUtils.h"

#include "pzcmesh.h"
#include "pzgmesh.h"
#include <TPZVTKGeoMesh.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZHCurlEquationFilter.h"
#include "pzstrmatrixflowtbb.h"
#include "pzstrmatrixot.h"
#include <chrono>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include "pzblockdiag.h"
#include "tpzsparseblockdiagonal.h"
#include "pzmatred.h"
#include "TPZSparseMatRed.h"
#include "TPZSpStructMatrix.h"
#include "pzbdstrmatrix.h"
#include "tpzverysparsematrix.h"
#include "TPZPardisoSolver.h"
#include "pzsysmp.h"
#include "pzysmp.h"
#include "TPZTimer.h"

// Util to print a summary of element information (mainly the connects) of a computational mesh
template <class TVar>
void TPZKernelHdivUtils<TVar>::PrintCMeshConnects(TPZCompMesh *cmesh){
    
    for (int i = 0; i < cmesh->NElements(); i++)
    {
        TPZCompEl *cel = cmesh->Element(i);
        cel->LoadElementReference();
        int matid = cel->Reference()->MaterialId();
        auto nconnects = cel->NConnects();
        std::cout << "Element = " << i << ", dim= " << cel->Dimension() << ",mat = " << cel->Reference()->MaterialId() << ", nconnects= " << nconnects << ": ";
        for (int j = 0; j < nconnects; j++)
        {
            std::cout << cel->ConnectIndex(j) << ", ";
        }
        std::cout << std::endl;
    
        // std::cout << cel->Connect() << std::endl;

        // // loop only over volumetric elements
        // if(matid != EDomain) continue;
        // if (cel->Reference()->Dimension() != dim) {
        //     DebugStop();
        // }

        // int nsides = cel->Reference()->NSides();
        // int ncorner = cel->Reference()->NCornerNodes();
        // for (int side = 0; side < nsides; side++) {
        //     if(cel->Reference()->SideDimension(side) != dim-1) continue;
        //     TPZGeoElSide gelside(cel->Reference(),side);
        //     TPZGeoElSide neighbour = gelside.Neighbour();
            
        //     std::cout << "Element = " << i << ", side = " << side  
        //             << ", NEl = " << neighbour.Element()->Index()
        //             << ", Nmatid = " << neighbour.Element()->MaterialId()
        //             << ", NNEl = " << neighbour.Neighbour().Element()->Index()
        //             << ", NNmatid = " << neighbour.Neighbour().Element() -> MaterialId() << std::endl;
        //     std::cout << "Neigh connect : " ;
        //     nconnects = neighbour.Element()->Reference()->NConnects();
        //     for (int j = 0; j < nconnects; j++)
        //     {
        //         std::cout << neighbour.Element()->Reference()->ConnectIndex(j) << ", ";
        //     }
        //     std::cout << std::endl;
        // }
    }
}

//Util to print the element properties of a geometric mesh
template <class TVar>
void TPZKernelHdivUtils<TVar>::PrintGeoMesh(TPZGeoMesh *gmesh){
    

    for (int i = 0; i < gmesh->NElements(); i++)
    {
        auto *gel = gmesh->Element(i);
        if (!gel) continue;

        int matid = gel->MaterialId();
        auto nsides = gel->NSides();
        // auto nconnects = gel->Reference()->NConnects();
        std::cout << "ELGeometric = " << i << ", dim= " << gel->Dimension() << ",mat = " << gel->MaterialId() << std::endl;
  
        nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        for (int side = 0; side < nsides; side++) {
            // if(gel->SideDimension(side) != 1) continue;
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            
            std::cout << "Element = " << i << ", side = " << side  
                    << ", NEL = " << neighbour.Element()->Index() 
                    << ", Nmatid = " << neighbour.Element()->MaterialId()
                    << ", NNEL = " << neighbour.Neighbour().Element()->Index() 
                    << ", NNmatid = " << neighbour.Neighbour().Element() -> MaterialId() << std::endl;
        }
    }

    //Prints gmesh mesh properties
    std::string vtk_name = "geoMesh.vtk";
    std::ofstream vtkfile(vtk_name.c_str());

    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

}

// Util to print the computational mesh
template <class TVar>
void TPZKernelHdivUtils<TVar>::PrintCompMesh(TPZCompMesh *cmesh,std::string file_name){

    // Print pressure mesh
    std::string txt = file_name + ".txt";
    std::ofstream myfile(txt);
    cmesh->Print(myfile);

    //Prints computational mesh properties
    std::string vtk_name = file_name + ".vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

}

// Util to solve the arising linear sistem by means of a direct method
template <class TVar>
void TPZKernelHdivUtils<TVar>::SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh, bool filterEquations, bool &domainHybridization)
{
    //sets number of threads to be used by the solver
    constexpr int nThreads{10};
    // TPZSkylineStructMatrix<REAL> matskl(cmesh);
    // TPZSSpStructMatrix<STATE> matskl(cmesh);
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matskl(cmesh);   
    
    // 
    // TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matskl(cmesh);
    // TPZSSpStructMatrix<STATE,TPZStructMatrixOT<STATE>> matskl(cmesh);
    // TPZSSpStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>> matskl(cmesh);
    matskl.SetNumThreads(nThreads);
    TPZHCurlEquationFilter<TVar> filter;
    //-----------------------
    if (filterEquations){    
        

        TPZVec<int64_t> activeEqs;
    
        if(filter.FilterEdgeEquations(cmesh, activeEqs, domainHybridization)){
            return;
        }
        // vertexData = filter.GetVertexDataStructure();
        // edgeData = filter.GetEdgeDataStructure();
        const int neqs = activeEqs.size();
        // std::cout << "ACtiveEqu - " << activeEqs << std::endl;
        matskl.EquationFilter().SetActiveEquations(activeEqs);
        std::cout << "Active equations = " << activeEqs.size() << std::endl;
    }
    //----------------------

    an.SetStructuralMatrix(matskl);
    

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    
    // TPZStepSolver<STATE> jac;
    // REAL tol = 1.e-30;
    // jac.SetJacobi(100,tol,0);
    // jac.ShareMatrix(step);
// #ifdef USING_MKL
    an.SetSolver(step);
// #endif
    //assembles the system

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    an.Assemble();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time Assemble = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;

    ///solves the system
    std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
    an.Solve();
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    std::cout << "Time Solve = " << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - begin2).count() << "[ms]" << std::endl;

    return;
}

// An util to solve the arising linear system employing an iterative method
template <class TVar>
void TPZKernelHdivUtils<TVar>::SolveProblemIterative(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    //sets number of threads to be used by the solver
    constexpr int nThreads{0};
    //defines storage scheme to be used for the FEM matrices
    //in this case, a symmetric skyline matrix is used
    // TPZSkylineStructMatrix<STATE> matskl(cmesh);
    // TPZSSpStructMatrix<STATE> matskl(cmesh);
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matskl(cmesh);
    // TPZSSpStructMatrix<STATE,TPZStructMatrixOT<STATE>> matskl(cmesh);
    // TPZSSpStructMatrix<STATE,TPZStructMatrixTBBFlow<STATE>> matskl(cmesh);
    

    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    ///Setting an iterative solver
    // TPZCopySolve<STATE> * precond = new TPZCopySolve<STATE>( matskl.Create() );  step.ShareMatrix( *precond );
    TPZStepSolver<STATE> * precond = new TPZStepSolver<STATE>( matskl.Create() ); step.ShareMatrix( *precond ); precond->SetJacobi(1, 0.0, 0);
    TPZStepSolver<STATE> jac;
    REAL tol = 1.e-30;
    jac.SetSSOR(1,1.1,0.,0);
    jac.ShareMatrix(step);
    // step.SetGMRES(2000,2000, *precond, tol, 0);
    step.SetCG(2000, *precond, tol, 0);
    an.SetSolver(step);

    //assembles the system
    an.Assemble();

    ///solves the system
    an.Solve();

    return;
}

//An Util to print the results of a multiphysics mesh to a .vtk file
template <class TVar>
void TPZKernelHdivUtils<TVar>::PrintResultsMultiphysics(TPZVec<TPZCompMesh *> &meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh)
{
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
    TPZManVector<std::string,10> scalnames(2), vecnames(2);


    scalnames[0] = "Pressure";
    scalnames[1] = "ExactPressure";
    vecnames[0]= "Flux";
    vecnames[1]= "ExactFlux";

    constexpr int resolution{0};
    std::string plotfile = "solutionMDFB.vtk";
    an.DefineGraphMesh(cmesh->Dimension(),scalnames,vecnames,plotfile);
    an.PostProcess(resolution,cmesh->Dimension());
    // Print mesh properties
    // std::ofstream out("mesh.txt");
    // an.Print("nothing",out);

}

// An Util to compute the error on Kernel Hdiv simulations
template <class TVar>
void TPZKernelHdivUtils<TVar>::ComputeError(TPZLinearAnalysis &an, std::ostream &anPostProcessFile)
{
    ///Calculating approximation error  
    TPZManVector<REAL,5> error;

    auto cmeshNew = an.Mesh();
    int64_t nelem = cmeshNew->NElements();
    cmeshNew->LoadSolution(cmeshNew->Solution());
    cmeshNew->ExpandSolution();
    cmeshNew->ElementSolution().Redim(nelem, 5);

    an.PostProcessError(error,false,anPostProcessFile);
        
    std::cout << "\nApproximation error:\n";
    std::cout << "H1 Norm = " << std::scientific << std::setprecision(15) << error[0]<<'\n';
    std::cout << "L1 Norm = " << std::scientific << std::setprecision(15) << error[1]<<'\n'; 
    std::cout << "H1 Seminorm = " << std::scientific << std::setprecision(15) << error[2]<<'\n'; 
    // std::cout << "H1 Seminorm = " << std::scientific << std::setprecision(15) << error[3]<<'\n'; 
    // std::cout << "H1 Seminorm = " << std::scientific << std::setprecision(15) << error[4]<<'\n'; 
    // std::cout << "error 4 = " << error[3]<<'\n'; 
    // std::cout << "error 5 = " << error[4] << "\n\n";
}


// An Util to compute the error on Kernel Hdiv simulations
template <class TVar>
void TPZKernelHdivUtils<TVar>::SolveProblemMatRed(TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh, std::set<int> &matIdBC)
{

    //HERE STARTS THE ITERATIVE SOLVER SET
    //sets number of threads to be used by the solver
    constexpr int nThreads{10};

    // Compute the number of equations in the system
    int64_t nEqFull = cmesh->NEquations();
    int64_t nEqPres, nEqFlux;

    ReorderEquations(cmesh,nEqPres,nEqFlux,matIdBC);

    PrintCompMesh(cmesh,"CMESH_reordered");
    std::cout << "NUMBER OF EQUATIONS:\n " << 
                 "\n Full problem = " << nEqFull << 
                 "\n Flux = " << nEqFlux << 
                 " \n Pressure = " << nEqPres << std::endl;

    // Create the RHS vectors
    TPZFMatrix<STATE> rhsFull(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsAux(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsFlux(nEqFlux,1,0.);
    TPZAutoPointer<TPZGuiInterface> guiInterface;

    //Creates the problem matrix    
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Stiffness(cmesh);
    Stiffness.SetNumThreads(nThreads);

    //Cria duas matrizes, para inverter a ordem das matrizes em bloco
    TPZMatRed<STATE, TPZFMatrix<STATE>> *matRed = new TPZMatRed<STATE, TPZFMatrix<STATE>>(nEqFull,nEqPres);

    std::ofstream out("out2.txt");

    //Primeiro cria a matriz auxiliar
    TPZFMatrix<REAL> K00(nEqPres,nEqPres,0.);
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //Altera range da matriz stiffness.
    Stiffness.SetEquationRange(0,nEqPres);
    
    //Transfere as submatrizes da matriz auxiliar para a matriz correta.
    step.SetMatrix(&K00);
    matRed->SetSolver(&step);
  
    Stiffness.EquationFilter().Reset();
    an.SetStructuralMatrix(Stiffness);

    TPZTimer clock;
    clock.start();
    //Monta a matriz auxiliar
    rhsAux.Zero();
    Stiffness.Assemble(*matRed,rhsAux,guiInterface);
    clock.stop();
    std::cout << "Time Assemble " << clock << std::endl;

    rhsFull=rhsAux;

    TPZBlockDiagonalStructMatrix<STATE> BDFmatrix(cmesh);
    BDFmatrix.SetEquationRange(nEqPres,nEqFull);
    TPZBlockDiagonal<REAL> KBD;
    BDFmatrix.AssembleBlockDiagonal(KBD);
    
    TPZFMatrix<REAL> *K11Red = new TPZFMatrix<REAL>(nEqFlux,nEqFlux);

    matRed->SetF(rhsFull);
    matRed->K00()->Print("K00=",out,EMathematicaInput);

    matRed->K11Reduced(*K11Red,rhsFlux);
    K11Red->Print("K11Red=",out,EMathematicaInput);
    // matRed->F1Red(rhsFlux);
    matRed->SetF(rhsFull);
    matRed->K11().Print("K11=",out,EMathematicaInput);
    KBD.Print("KBD=",out,EMathematicaInput);
    matRed->F1().Print("F1=",out,EMathematicaInput);
    matRed->K00()->Print("K00=",out,EMathematicaInput);
    matRed->K01().Print("K01=",out,EMathematicaInput);
    matRed->K10().Print("K10=",out,EMathematicaInput);
    rhsFlux.Print("RHSFlux = ",out,EMathematicaInput);

    //Creates the preconditioner 
    // TPZCopySolve<STATE> *precond = new TPZCopySolve<STATE>( &matRed->K11() );
    // TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>( &matRed->K11() );
    TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>( &KBD );
    precond->SetDirect(ELU);
    // precond->SetJacobi(1,1.e-6,0);
    // precond->SetCopySolve();
    // precond->Solve(matRed->F1(),rhsAux);
    // std::cout << "matRed->F1() " << matRed->F1() << std::endl;
    // std::cout << "rhsAux " << rhsAux << std::endl;
    int64_t nMaxIter = 50;
    TPZVec<REAL> errors(nMaxIter);
    errors.Fill(0.);

    // for (int64_t iter = 1; iter < nMaxIter; iter++){
        
        // std::cout << "ITER = " << iter << std::endl;
        TPZFMatrix<STATE> residual(nMaxIter,1,0.);
        REAL tol = 1.e-10;
        TPZFMatrix<STATE> solution(nEqFlux,1,0.);

        clock.start();
        K11Red->SolveCG(nMaxIter,*precond,rhsFlux,solution,&residual,tol);
        clock.stop();
        std::cout << "Time CG " << clock << std::endl;

        REAL norm = 0.;
        std::cout << "Number of CG iterations = " << nMaxIter << " , residual = " << tol << std::endl;
        TPZFMatrix<STATE> press(nEqPres,1,0.);
        matRed->UGlobal(solution,press);

        //Update solution in Analysis
        for (int i = 0; i < nEqFlux; i++){
            rhsFull(nEqPres+i,0) = solution(i,0);
        }
        for (int i = 0; i < nEqPres; i++){
            rhsFull(i,0) = press(i,0);
        }

        an.Solution()=rhsFull;
        an.LoadSolution();
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

template <class TVar>
void TPZKernelHdivUtils<TVar>::ReorderEquations(TPZMultiphysicsCompMesh* cmesh, int64_t &nEqPres, int64_t &nEqFlux, std::set<int> &matIdBC)
{
    cmesh->LoadReferences();
    std::set<int64_t> auxConnectsP, auxConnectsF;
    // TPZSkylineStructMatrix<STATE> Pmatrix(cmesh);
    // TPZSkylineStructMatrix<STATE> Fmatrix(cmesh);
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Pmatrix(cmesh);
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Fmatrix(cmesh);
    auto neqsTotal = cmesh->NEquations();

    auto gmesh = cmesh->Reference();
    nEqPres = 0;
    nEqFlux = 0;

    // PrintCompMesh(cmesh,"CMESH_0");


    // loop over the connects and create two std::sets one to the "pressure" ones, representing the 
    // degrees of freedom to be condensed. The second set contains the "flux" connects, which will not be condensed
    // This can change depending on the problem.
    for (auto gel:gmesh->ElementVec()){
        
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
            nEqPres += neq;
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
            nEqFlux += neq;
            seqNum=con.SequenceNumber();
            cmesh->Block().Set(seqNum,neq);
        }
    }   
    cmesh->ExpandSolution();

    
}



// An Util to compute the error on Kernel Hdiv simulations
template <class TVar>
void TPZKernelHdivUtils<TVar>::SolveProblemMatRedSparse(TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh, std::set<int> &matIdBC)
{

    //HERE STARTS THE ITERATIVE SOLVER SET
    //sets number of threads to be used by the solver
    constexpr int nThreads{0};

    // Compute the number of equations in the system
    int64_t nEqFull = cmesh->NEquations();
    int64_t nEqPres, nEqFlux;

    ReorderEquations(cmesh,nEqPres,nEqFlux,matIdBC);

    PrintCompMesh(cmesh,"CMESH_reordered");
    std::cout << "NUMBER OF EQUATIONS:\n " << 
                 "\n Full problem = " << nEqFull << 
                 "\n Flux = " << nEqFlux << 
                 " \n Pressure = " << nEqPres << std::endl;

    // Create the RHS vectors
    TPZFMatrix<STATE> rhsFull(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsAux(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsFlux(nEqFlux,1,0.);
    TPZAutoPointer<TPZGuiInterface> guiInterface;

    //Creates the problem matrix    
    TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> Stiffness(cmesh);
    Stiffness.SetNumThreads(nThreads);

    //Cria duas matrizes, para inverter a ordem das matrizes em bloco
    // TPZSparseMatRed<STATE, TPZVerySparseMatrix<STATE>> *matRed = new TPZSparseMatRed<STATE, TPZVerySparseMatrix<STATE>>(nEqFull,nEqPres);
    // TPZSparseMatRed<STATE, TPZSYsmpMatrix<STATE>> *matRed = new TPZSparseMatRed<STATE, TPZSYsmpMatrix<STATE>>(nEqFull,nEqPres);
    TPZSparseMatRed<STATE> *matRed = new TPZSparseMatRed<STATE>(nEqFull,nEqPres);

    std::ofstream out("out.txt");


    //Primeiro cria a matriz auxiliar
    TPZSYsmpMatrix<REAL> K00(nEqPres,nEqPres);
    TPZSYsmpMatrix<REAL> K11(nEqFlux,nEqFlux);
    // TPZSYsmpMatrix<REAL> K00(nEqPres,nEqPres);
    
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //Altera range da matriz stiffness.
    // Stiffness.SetEquationRange(0,nEqPres);
    {
        // TPZSYsmpMatrix<REAL> *Stiff = dynamic_cast<TPZSYsmpMatrix<REAL> *>(Stiffness.Create());
        // K00.SetData(Stiff->IA(),Stiff->JA(),Stiff->A());
        // Stiffness.Assemble(K00,rhsFull,guiInterface);

        // std::cout << "IA K00 = " << Stiff->IA() << std::endl;
        // std::cout << "JA K00 = " << Stiff->JA() << std::endl;
        // std::cout << "A K00 = " << Stiff->A() << std::endl;

        Stiffness.EquationFilter().Reset();
        // Stiffness.SetEquationRange(nEqPres,nEqFull);
        TPZSYsmpMatrix<REAL> *StiffK11 = dynamic_cast<TPZSYsmpMatrix<REAL> *>(Stiffness.Create());

        //Fazer uma rotina para separar IA, JA e A de K01, K10 e K11;
        TPZVec<int64_t> IA_K00(nEqPres+1,0), IA_K01(nEqPres+1,0), IA_K10(nEqFlux+1,0), IA_K11(nEqFlux+1,0);
        
        std::vector<int64_t> auxK00, auxK01, auxK10, auxK11;

        IA_K00[0] = 0;
        IA_K01[0] = 0;
        IA_K10[0] = 0;
        IA_K11[0] = 0;
        for (int i = 0; i < nEqPres; i++){
            for (int j = StiffK11->IA()[i]; j < StiffK11->IA()[i+1]; j++){
                if (StiffK11->JA()[j] < nEqPres){
                    // Faz parte da matriz K00
                    // std::cout << "Termo IA = " << StiffK11->IA()[i] << ", JA = " << StiffK11->JA()[j] << ", K00\n";
                    auxK00.push_back(StiffK11->JA()[j]);
                } else {
                    // Faz parte da matriz K01
                    // std::cout << "Termo IA = " << StiffK11->IA()[i] << ", JA = " << StiffK11->JA()[j] << ", K01\n";
                    auxK01.push_back(StiffK11->JA()[j]-nEqPres);
                }
            }
            IA_K00[i+1] = auxK00.size();
            IA_K01[i+1] = auxK01.size();
        }
        for (int i = nEqPres; i < nEqFull; i++){
            for (int j = StiffK11->IA()[i]; j < StiffK11->IA()[i+1]; j++){
                if (StiffK11->JA()[j] >= nEqPres){
                    // Faz parte da matriz K11
                    // std::cout << "Termo IA = " << StiffK11->IA()[i] << ", JA = " << StiffK11->JA()[j] << ", K01\n";
                    auxK11.push_back(StiffK11->JA()[j]-nEqPres);
                }
            }
            // IA_K00[i+1] = auxK00.size();
            IA_K11[i-nEqPres+1] = auxK11.size();
        }
        
        //Do the transpose - Matriz K10
        IA_K10[0]=0;
        for (int i = 0  ; i < nEqFlux; i++){
            int nNonZeros = std::count(auxK01.begin(),auxK01.end(),i);
            IA_K10[i+1] = IA_K10[i] + nNonZeros;
        }
        // std::cout << "IA_K10 " << IA_K10 << std::endl;
        auxK10.resize(auxK01.size());

        // //Simetrização da matriz K11;
        // std::vector<int64_t> auxK11_sym;
        // for (int i = 0; i < nEqFlux; i++){
        //     //Adiciona os termos da matriz triangular inferior
        //     if (i>0){
        //         for (int k = 0; k < i; k++){
        //             for (int j = IA_K11[k]; j < IA_K11[k+1]; j++)
        //             {
        //                 if (auxK11[j] == i){
        //                     auxK11_sym.push_back(auxK11[k]);
        //                 }
        //             }   
        //         } 
        //     }

        //     //Adiciona os termos da matriz triangular superior
        //     for (int j = IA_K11[i]; j < IA_K11[i+1]; j++)
        //     {
        //         auxK11_sym.push_back(auxK11[j]);
        //     }
        // }
        // //Atualiza IA_K11
        // for (int i = 0  ; i < nEqFlux; i++){
        //     int nNonZeros = std::count(auxK11_sym.begin(),auxK11_sym.end(),i);
        //     IA_K11[i+1] = IA_K11[i] + nNonZeros;
        // }
        



        TPZVec<int64_t> JA_K00(auxK00.size(),0), JA_K01(auxK01.size(),0), JA_K10(auxK01.size(),0), JA_K11(auxK11.size(),0);
        TPZVec<double> A_K00(auxK00.size(),0.), A_K01(auxK01.size(),0.), A_K10(auxK01.size(),0.), A_K11(auxK11.size(),0.);
        
        for (int i = 0; i < JA_K00.size(); i++) JA_K00[i] = auxK00[i];
        for (int i = 0; i < JA_K01.size(); i++) JA_K01[i] = auxK01[i];
        // for (int i = 0; i < JA_K10.size(); i++) JA_K10[i] = auxK10[i];
        for (int i = 0; i < JA_K11.size(); i++) JA_K11[i] = auxK11[i];
        
        // std::cout << "IA_K01 " << IA_K01 << std::endl;
        // std::cout << "JA_K01 " << JA_K01 << std::endl;

        // std::cout << "IA = " << StiffK11->IA() << std::endl;
        // std::cout << "JA = " << StiffK11->JA() << std::endl;
        // std::cout << "A = " << StiffK11->A() << std::endl;

        //Aloca estrutura das matrizes esparsas
        K00.SetData(IA_K00,JA_K00,A_K00);
        matRed->K01().SetData(IA_K01,JA_K01,A_K01);
        matRed->K01().SetData(IA_K01,JA_K01,A_K01);
        matRed->K10().SetData(IA_K10,JA_K10,A_K10);
        matRed->K11().SetData(IA_K11,JA_K11,A_K11);

       
    }
    
    //Transfere as submatrizes da matriz auxiliar para a matriz correta.
    step.SetMatrix(&K00);
    matRed->SetSolver(&step);
  
    Stiffness.EquationFilter().Reset();
    an.SetStructuralMatrix(Stiffness);

    //Monta a matriz auxiliar
    rhsAux.Zero();
    TPZTimer clock;
    clock.start();
    Stiffness.Assemble(*matRed,rhsAux,guiInterface);
    clock.stop();
    std::cout << "Time Assemble " << clock << std::endl;

    rhsFull=rhsAux;

    TPZBlockDiagonalStructMatrix<STATE> BDFmatrix(cmesh);
    BDFmatrix.SetEquationRange(nEqPres,nEqFull);
    TPZBlockDiagonal<REAL> KBD;
    BDFmatrix.AssembleBlockDiagonal(KBD);
    
    TPZFMatrix<REAL> *K11Red = new TPZFMatrix<REAL>(nEqFlux,nEqFlux,0.);
    TPZSYsmpMatrix<REAL> *K11Sparse = new TPZSYsmpMatrix<REAL>(nEqFlux,nEqFlux);
    // TPZVerySparseMatrix<REAL> *K11Red2 = dynamic_cast<TPZVerySparseMatrix<REAL> *>(K11Red);

    matRed->SetF(rhsFull);
    matRed->K00()->Print("K00=",out,EMathematicaInput);

    matRed->K11Reduced(*K11Red,rhsFlux);
    K11Red->Print("K11Red=",out,EMathematicaInput);
    // matRed->F1Red(rhsFlux);
    matRed->SetF(rhsFull);
        matRed->K11().Print("K11=",out,EMathematicaInput);

    KBD.Print("KBD=",out,EMathematicaInput);
    matRed->F1().Print("F1=",out,EMathematicaInput);
    matRed->K00()->Print("K00=",out,EMathematicaInput);
    matRed->K01().Print("K01=",out,EMathematicaInput);
    matRed->K10().Print("K10=",out,EMathematicaInput);
    rhsFlux.Print("RHSFlux = ",out,EMathematicaInput);

//     rhsFlux = { 0.33333331813141026, -2.5037841787550497e-16 , 2.0354088869182159e-16, -0.33333331813141037 , -0.33333331813141026 ,
//  0.33333331813141048 ,
//  4.3165702128649358e-16 ,
//  3.7007434998048114e-17 ,
//  -0.33333331813141037 ,
//  0.33333331813141037 ,
//  -0.33333331813141043 ,
//  0.33333331813141004  };

    //Creates the preconditioner 
    // TPZCopySolve<STATE> *precond = new TPZCopySolve<STATE>( &matRed->K11() );
    // TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>( &matRed->K11() );
    TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>( &KBD );
    precond->SetDirect(ELU);
    // precond->SetJacobi(1,1.e-6,0);
    // precond->SetCopySolve();
    // precond->Solve(matRed->F1(),rhsAux);
    // std::cout << "matRed->F1() " << matRed->F1() << std::endl;
    // std::cout << "rhsAux " << rhsAux << std::endl;
    int64_t nMaxIter = 30;
    TPZVec<REAL> errors(nMaxIter);
    errors.Fill(0.);

    // for (int64_t iter = 1; iter < nMaxIter; iter++){
        
        // std::cout << "ITER = " << iter << std::endl;
        TPZFMatrix<STATE> residual(nMaxIter,1,0.);
        REAL tol = 1.e-10;
        TPZFMatrix<STATE> solution(nEqFlux,1,0.);

        clock.start();
        K11Red->SolveCG(nMaxIter,*precond,rhsFlux,solution,&residual,tol);
        // K11RedSparse->SolveCG(nMaxIter,*precond,rhsFlux,solution,&residual,tol);
        clock.stop();
        std::cout << "Time CG " << clock << std::endl;

        REAL norm = 0.;
        std::cout << "Number of CG iterations = " << nMaxIter << " , residual = " << tol << std::endl;
        TPZFMatrix<STATE> press(nEqPres,1,0.);
        matRed->UGlobal(solution,press);

        //Update solution in Analysis
        for (int i = 0; i < nEqFlux; i++){
            rhsFull(nEqPres+i,0) = solution(i,0);
        }
        for (int i = 0; i < nEqPres; i++){
            rhsFull(i,0) = press(i,0);
        }

        an.Solution()=rhsFull;
        an.LoadSolution();
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


template class TPZKernelHdivUtils<STATE>;
// template class TPZKernelHdivUtils<CSTATE>;


