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
    TPZVec<int> removed_edges(nnodes, 0);


    /* We need to keep ONE edge, globally speaking*/
    done_vertices[0] = true;

    /**
        Contains all the connects marked for removal. It is expected
        that all the connects are associated with edges.
    */
    std::set<int64_t> removed_connects;

    /**
        Contains coplanar connects (to removed_connects) that should not be removed.
    */
    std::set<int> no_remove;

    constexpr int edgeDim{1};

    auto nElements = gmesh->NElements();

    //Loop over the elements
    for (auto gel : gmesh->ElementVec()){
    
        //Only tetrahedra is considered in the algorithm
        if (gel->Dimension() < 3) continue;
        
        const auto nEdges = gel->NSides(edgeDim);
        const auto firstEdge = gel->FirstSide(edgeDim);
        const auto firstFace = firstEdge + nEdges;
        const auto index = gel->Index();
        const auto ncorner = gel->NCornerNodes();

        // Number of faces removed
        TPZManVector<int> nface_removed(4,0);
        // Edge blocked - 1 if it can't be removed
        TPZManVector<int> block_edges(6,0);

        //First - check which edges have already been treated and relates them to the element faces
        for(auto ie = firstEdge; ie < firstFace; ie++){
            
            TPZGeoElSide edge(gel,ie);
            const auto con = edge.Element()->Reference()->ConnectIndex(ie-firstEdge);
            
            //Store vertex edge connects for checking purposes - not part of the algorithm
            //------------------------------------
            const auto v1 = edge.SideNodeIndex(0);
            const auto v2 = edge.SideNodeIndex(1);
            
            vertex_edge_connects[v1].insert(con);
            vertex_edge_connects[v2].insert(con);
            //------------------------------------

            // checks if the edge has already been treated
            if (removed_connects.find(con) != removed_connects.end()) {
                
                // Gets all the element faces containing the edge ie
                TPZStack<TPZGeoElSide> gelsides;
                gel->AllHigherDimensionSides(ie,2,gelsides);
                
                // Increment the faces with the edge ie 
                for (int iside = 0; iside < gelsides.size(); iside++)
                {
                    int face = gelsides[iside].Side() - firstFace;
                    nface_removed[face]++;
                    if (nface_removed[face] > 2){
                        DebugStop();
                    }
                }//iside
            }//if connect has been treated
        }//edges

        //Second - Blocks the remaining edge in
        for(auto ie = firstEdge; ie < firstFace; ie++){
            
            TPZGeoElSide edge(gel,ie);
            const auto con = edge.Element()->Reference()->ConnectIndex(ie-firstEdge);

            // Gets all the element faces containing the edge ie
            TPZStack<TPZGeoElSide> gelsides;
            gel->AllHigherDimensionSides(ie,2,gelsides);
            
            // Blocks remaining edge in the faces with 2 edges already treated
            for (int iside = 0; iside < gelsides.size(); iside++)
            {
                int face = gelsides[iside].Side() - firstFace;
                // If 2 edges were removed, block the third edge;
                // If 3 edges were removed - forbiden! DebugStop();
                // Do nothing otherwise.
                if (nface_removed[face] == 2) {
                    block_edges[ie-firstEdge]=1;
                    //this connect should not be removed in the next elements too.
                    // no_remove.insert(con); 
                } else if (nface_removed[face] > 2){
                    DebugStop();
                }
            }//iside
        }//edge

        //Third - Removes a connect, if possible. Loops over the corner nodes
        for (int icorner = 0; icorner < ncorner; icorner++)
        {            
            // If the node has been treated, go to the next one
            auto inode = gel->NodeIndex(icorner);
            if (done_vertices[inode]) {
                continue;
            }

            // Gets the edges with the corner node
            TPZStack<TPZGeoElSide> cornerEdges;
            gel->AllHigherDimensionSides(icorner,1,cornerEdges);

            // Loops over the edges with the corner
            for (int iedge = 0; iedge < cornerEdges.size(); iedge++)
            {
                int edge = cornerEdges[iedge].Side() - ncorner;
                const auto con = gel->Reference()->ConnectIndex(edge);

                // /* Checks if the connect can be removed*/
                // if (no_remove.find(con) != no_remove.end()) {
                //     done_vertices[inode] = true;
                //     break;
                // }

                // If the edge has already been removed or it is blocked, go to the next one
                if(removed_connects.find(con) != removed_connects.end()) {
                    done_vertices[inode] = true;
                    break;
                }
                if(block_edges[edge] == 1) continue;
                
                // Else, remove the edge
                removed_connects.insert(con);
                done_vertices[inode] = true;

                // Gets the faces with the edge
                TPZStack<TPZGeoElSide> gelsides;
                gel->AllHigherDimensionSides(edge+ncorner,2,gelsides);
                
                for (int iface = 0; iface < gelsides.size(); iface++)
                {
                    //Increment the face counter with the edge
                    int face = gelsides[iface].Side() - firstFace;
                    nface_removed[face]++;
                    if (nface_removed[face] > 2){
                        DebugStop();
                    }
                    // If the counter is 2, then block the remaining edge
                    if (nface_removed[face] == 2) {
                        
                        TPZStack<int> smallsides;
                        gel->LowerDimensionSides(face+firstFace,smallsides);

                        for (size_t isize = 0; isize < smallsides.size(); isize++)
                        {
                            int side = smallsides[isize];
                            if (side < ncorner) continue;
                            
                            block_edges[side-ncorner] = 1;
                            // no_remove.insert(gel->Reference()->ConnectIndex(side-ncorner));
                        }                       
                    }   
                }//iside
            }//iedge
        }//icorner
    }//gel

    // removed_connects.insert(1);
    // removed_connects.insert(2);
    // removed_connects.insert(7);


    bool check_all_vertices{true};
    for(auto v : done_vertices){
        check_all_vertices = check_all_vertices && v;
        if(!v){
            break;
        }
    }
    
    std::cout << "Removed (" << removed_connects.size() << ") = ";
    for (auto rem : removed_connects)
    {
        std::cout << rem << " " ;
    }
    std::cout << std::endl;

    std::cout << "No_remove (" << no_remove.size() << ") = ";
    for (auto rem : no_remove)
    {
        std::cout << rem << " " ;
    }
    std::cout << std::endl;

    std::cout << "Done_vertices (" << done_vertices.size() << ") = ";
    for (auto rem : done_vertices)
    {
        std::cout << rem << " " ;
    }
    std::cout << std::endl;

    std::cout << "vertex_edge_connects \n";
    for (int i = 0; i < nnodes; i++)
    {
        std::cout << "node = " << i << "; ";
        int count = 0;
        for (auto rem : vertex_edge_connects[i])
        {
            if (removed_connects.find(rem) != removed_connects.end()){
                std::cout << "{" << rem << "} " ;
                count++;
            } else {
                std::cout << rem << " " ;
            }
            
        }
        std::cout <<";             " << count << " edge(s) removed" << std::endl;
        if (count == 0 && i != 0) std::cout << "PROBLEM! Node " << i << " did not removed an edge!\n" ;
        
    }
    
    // if(check_all_vertices == false) DebugStop();
    
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

    // if (check_edges_left == false) DebugStop();
    
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
    
    //Check if a face has all the edges removed
    for(auto gel : gmesh->ElementVec()){
        const auto nEdges = gel->NSides(edgeDim);
        const auto firstEdge = gel->FirstSide(edgeDim);
        const auto firstFace = firstEdge + nEdges;

        int aux = 0;

        std::set<int> local_removed;
        TPZVec<bool> local_index(firstFace,false);
        for(auto ie = firstEdge; ie < firstFace; ie++){
            
            TPZGeoElSide edge(gel,ie);
        
            const auto con = edge.Element()->Reference()->ConnectIndex(ie-firstEdge);

            /* checks if edge has been treated already */
            if (removed_connects.find(con) != removed_connects.end()) {
                aux++;
                local_removed.insert(con);
                local_index[ie-firstEdge] = true;
            }
        }

        if (gel->Dimension() == 3){
            std::cout << "gel = " << gel->Index() << ", removed connects = ";
            for (auto rem : local_removed)
            {
                std::cout << rem << ", " ;
            }
            // Checks of coplanar edges
            if ((local_index[0] && local_index[1] && local_index[2]) ||
                (local_index[0] && local_index[3] && local_index[4]) ||
                (local_index[2] && local_index[3] && local_index[5]) ||
                (local_index[1] && local_index[5] && local_index[4])) 
                std::cout << "Problem with this element: coplanar edges removed!" ;
            // // Checks of edges sharing the same node
            // if ((local_index[0] && local_index[2] && local_index[3]) ||
            //     (local_index[0] && local_index[1] && local_index[4]) ||
            //     (local_index[2] && local_index[1] && local_index[5]) ||
            //     (local_index[3] && local_index[5] && local_index[4])) 
            //     std::cout << "Problem with this element: at least 3 edges sharing the same node!" ;
            std::cout << std::endl;
        }

        if (aux == nEdges) {
            std::cout << "All edges have been removed in GeoEl " << gel->Index()<< std::endl;
        }

    }

    if (removed_connects.size() != nnodes-1){
        std::cout << "Removed " << removed_connects.size() << "/" << nnodes-1 << " allowed connects."  << std::endl;
        DebugStop();
    }


    return 0;
}

template class TPZKernelHdivUtils<STATE>;
// template class TPZKernelHdivUtils<CSTATE>;


