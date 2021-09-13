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
void TPZKernelHdivUtils<TVar>::SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    //sets number of threads to be used by the solver
    constexpr int nThreads{0};
    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(nThreads);

    //-----------------------
    TPZVec<int64_t> activeEqs;
  
    if(FilterEdgeEquations(cmesh, activeEqs)){
        return;
    }
    const int neqs = activeEqs.size();
    
    matskl.EquationFilter().SetActiveEquations(activeEqs);

    //----------------------

    an.SetStructuralMatrix(matskl);

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);

    // std::string fluxFile = "FluxCMesh";
    // // this->PrintCompMesh(cmesh,fluxFile);
    // std::filebuf fb;
    // fb.open ("test.txt",std::ios::out);
    // std::ostream os(&fb);
    // an.Print(fluxFile,os);

    //assembles the system
    an.Assemble();

    ///solves the system
    an.Solve();

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
    TPZSkylineStructMatrix<STATE> matskl(cmesh);
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
    step.SetGMRES(2000,2000, *precond, tol, 0);
    // step.SetCG(2000, *precond, tol, 0);
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

    int div = 0;
    std::string plotfile = "solutionMDFB.vtk";
    an.DefineGraphMesh(cmesh->Dimension(),scalnames,vecnames,plotfile);
    an.PostProcess(div,cmesh->Dimension());
    // Print mesh properties
    // std::ofstream out("mesh.txt");
    // an.Print("nothing",out);

}

// An Util to compute the error on Kernel Hdiv simulations
template <class TVar>
void TPZKernelHdivUtils<TVar>::ComputeError(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile)
{
    ///Calculating approximation error  
    TPZManVector<REAL,5> error;

    auto cmeshNew = an.Mesh();
    int64_t nelem = cmeshNew->NElements();
    cmeshNew->LoadSolution(cmeshNew->Solution());
    cmeshNew->ExpandSolution();
    cmeshNew->ElementSolution().Redim(nelem, 5);

    an.PostProcessError(error);
        
    std::cout << "\nApproximation error:\n";
    std::cout << "H1 Norm = " << std::scientific << std::setprecision(15) << error[0]<<'\n';
    std::cout << "L1 Norm = " << std::scientific << std::setprecision(15) << error[1]<<'\n'; 
    std::cout << "H1 Seminorm = " << std::scientific << std::setprecision(15) << error[2]<<'\n'; 
    // std::cout << "error 4 = " << error[3]<<'\n'; 
    // std::cout << "error 5 = " << error[4] << "\n\n";
}



template <class TVar>
bool TPZKernelHdivUtils<TVar>::FilterEdgeEquations(TPZAutoPointer<TPZCompMesh> cmesh,
                    TPZVec<int64_t> &activeEquations)
{

    /**TODO:
     1: for each vertex: at least one edge must be removed
     2: test for p = 1.
     3: create the class for filtering the functions in the master el.
     3.1: phi_f^{e,n}: remove one for each face.
     4: test for p = 2.
     5: think about the remaining functions.
   */
    const auto gmesh = cmesh->Reference();
    
    const auto nnodes = gmesh->NNodes();

    /**
         The i-th position contains the indices of all the 
        connects associated with the edges adjacent to the i-th node
    */
    TPZVec<std::set<int>> vertex_edge_connects(nnodes);

    /**
         The i-th is true if the i-th node has already been dealt with.
    */
    TPZVec<bool> done_vertices(nnodes, false);


    /*
    we need to keep ONE edge, globally speaking*/
    done_vertices[0] = true;

    
    /**
         Contains all the connects marked for removal. It is expected
        that all the connects are associated with edges.
    */
    std::set<int> removed_connects;

    /**
         Contains coplanar connects (to removed_connects) that should not be removed.
    */
    std::set<int> no_remove;
    constexpr int edgeDim{1};

    for(auto gel : gmesh->ElementVec()){
        const auto nEdges = gel->NSides(edgeDim);
        const auto firstEdge = gel->FirstSide(edgeDim);
        const auto firstFace = firstEdge + nEdges;

        //Auxiliary variable to avoid the selection of three connects of the same face
        int aux = 0;

        for(auto ie = firstEdge; ie < firstFace; ie++){
            TPZGeoElSide edge(gel,ie);
        
            const auto con = edge.Element()->Reference()->ConnectIndex(ie-firstEdge);

            /* checks if edge has been treated already */
            if (removed_connects.find(con) != removed_connects.end()) {
                continue;
            }
            
            const auto v1 = edge.SideNodeIndex(0);
            const auto v2 = edge.SideNodeIndex(1);
            
            vertex_edge_connects[v1].insert(con);
            vertex_edge_connects[v2].insert(con);
            //both vertices have been treated
            if(done_vertices[v1] && done_vertices[v2]){
                continue;
            } else {
                /* If aux = 2, then we can't take the next connect from this side.
                  We also store the connect that can't be removed amd check it in the next elements.*/
                if (aux == 2) {
                    no_remove.insert(con);
                    continue;
                }
                if (no_remove.find(con) != no_remove.end()) {
                    continue;
                }
                /*either one or none vertex has been treated,
                so we need to remove the edge*/
                removed_connects.insert(con);
                if(!done_vertices[v1] && !done_vertices[v2]){
                    /*if no vertex has been treated we
                        mark ONE OF THEM as treated*/
                    done_vertices[v1] = true;
                    aux++;
                } else {
                    /*one of them had been treated already,
                        instead of checking which, we mark both as true*/
                    done_vertices[v1] = true;
                    done_vertices[v2] = true;
                }
            }
        }
    }


    bool check_all_vertices{true};
    for(auto v : done_vertices){
        check_all_vertices = check_all_vertices && v;
        if(!v){
        break;
        }
    }

    // for (auto rem : removed_connects)
    // {
    //     std::cout << rem << " " ;
    // }
    // std::cout << std::endl;
    

    if(check_all_vertices == false) DebugStop();
    
    bool check_edges_left{true};
    for(auto iv = 0; iv < nnodes; iv++){
        const auto all_edges = vertex_edge_connects[iv];
        bool local_check{false};
        for(auto edge: all_edges){
        if(removed_connects.find(edge) == removed_connects.end()){
            local_check = true;
            break;
        }
        }
        check_edges_left = local_check && check_edges_left;
    }

    if (check_edges_left == false) DebugStop();
    
    activeEquations.Resize(0);
    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
        if (removed_connects.find(iCon) == removed_connects.end()) {
        auto &con = cmesh->ConnectVec()[iCon];
        if (con.HasDependency()){
            continue;
        }
        const auto seqnum = con.SequenceNumber();
        const auto pos = cmesh->Block().Position(seqnum);
        const auto blocksize = cmesh->Block().Size(seqnum);
        if (blocksize == 0){
            continue;
        }
        
        const auto vs = activeEquations.size();
        activeEquations.Resize(vs + blocksize);
        for (auto ieq = 0; ieq < blocksize; ieq++) {
            activeEquations[vs + ieq] = pos + ieq;
        }
        }
    }
    
    return 0;


}


template class TPZKernelHdivUtils<STATE>;
// template class TPZKernelHdivUtils<CSTATE>;


