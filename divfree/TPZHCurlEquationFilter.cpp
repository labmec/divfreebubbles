//
// Created by Jeferson Fernandes on 16/09/21.
//

#include "TPZHCurlEquationFilter.h"

#include "pzcmesh.h"
#include "pzgmesh.h"
#include <TPZVTKGeoMesh.h>

template <class TVar>
void TPZHCurlEquationFilter<TVar>::CreateFilterDataStructure(TPZGeoMesh* gmesh, const int &edgeDim,
                                                         TPZVec<std::set<int>> &vertex_edge_connects)
{
    //Loop over the elements
    for (auto gel : gmesh->ElementVec()){
    
        //Only tetrahedra is considered in the algorithm
        if (gel->Dimension() < 3) continue;
        
        const auto nEdges = gel->NSides(edgeDim);
        const auto firstEdge = gel->FirstSide(edgeDim);
        const auto firstFace = firstEdge + nEdges;
        const auto index = gel->Index();
        const auto ncorner = gel->NCornerNodes();

        //First - check which edges have already been treated and relates them to the element faces
        for(auto ie = firstEdge; ie < firstFace; ie++){
            
            TPZGeoElSide edge(gel,ie);
            const auto con = edge.Element()->Reference()->ConnectIndex(ie-firstEdge);
            
            //Store the connects of a given vertex 
            const auto v1 = edge.SideNodeIndex(0);
            const auto v2 = edge.SideNodeIndex(1);
            
            vertex_edge_connects[v1].insert(con);
            vertex_edge_connects[v2].insert(con);
        }
    }


}

template <class TVar>
bool TPZHCurlEquationFilter<TVar>::FilterEdgeEquations(TPZAutoPointer<TPZCompMesh> cmesh,
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


    CreateFilterDataStructure(gmesh,edgeDim,vertex_edge_connects);

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

template class TPZHCurlEquationFilter<STATE>;
// template class TPZHCurlEquationFilter<CSTATE>;


