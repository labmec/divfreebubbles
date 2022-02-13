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
void TPZKernelHdivUtils<TVar>::PrintCompMesh(TPZCompMesh *cmesh,std::string &file_name){

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
    constexpr int nThreads{12};
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
    constexpr int nThreads{10};
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
    // TPZMatrixSolver<STATE> * precond = an.BuildPreconditioner(TPZAnalysis::EBlockJacobi , true);
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
    TPZManVector<std::string,10> scalnames(0), vecnames(2);


    // scalnames[0] = "Pressure";
    // scalnames[1] = "ExactPressure";
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


template class TPZKernelHdivUtils<STATE>;
// template class TPZKernelHdivUtils<CSTATE>;


