//
// Created by Jeferson Fernandes on 16/09/21.
//

#include "TPZHCurlEquationFilter.h"

#include "pzcmesh.h"
#include "pzgmesh.h"
#include <TPZVTKGeoMesh.h>

template <class TVar>
void TPZHCurlEquationFilter<TVar>::InitDataStructures(TPZGeoMesh* gmesh, 
                                                      TPZVec<VertexFilter> &mVertex,
                                                      std::map<int,EdgeFilter> &mEdge,
                                                      std::map<int,FaceFilter> &mFace,
                                                      std::map<int,std::set<int>> &freeEdges)
{
    const int edgeDim{1};

    //Loop over the elements
    for (auto gel : gmesh->ElementVec()){
    
        //Only tetrahedra is considered in the algorithm
        if (gel->Dimension() < 3) continue;
        
        const auto nEdges = gel->NSides(edgeDim);
        const auto firstEdge = gel->FirstSide(edgeDim);
        const auto firstFace = firstEdge + nEdges;
        const auto index = gel->Index();
        const auto ncorner = gel->NCornerNodes();
        const auto nsides = gel->NSides();

        for(auto ie = firstEdge; ie < firstFace; ie++){
            // Vertex data structure            
            TPZGeoElSide edge_(gel,ie);
            const auto con = edge_.Element()->Reference()->ConnectIndex(ie-firstEdge);
            
            //Store the connects of a given vertex 
            const auto v1 = edge_.SideNodeIndex(0);
            const auto v2 = edge_.SideNodeIndex(1);
            
            mVertex[v1].edge_connect.insert(con);
            mVertex[v2].edge_connect.insert(con);

            mEdge[con].vertex_connect.insert(v1);
            mEdge[con].vertex_connect.insert(v2);

            freeEdges[con].insert(v1);
            freeEdges[con].insert(v2);
            
            TPZStack<TPZGeoElSide> gelsides;
            gel->AllHigherDimensionSides(ie,2,gelsides);

            for (int i = 0; i < gelsides.size(); i++)
            {   
                // Edge data structure
                const auto con_face = edge_.Element()->Reference()->ConnectIndex(gelsides[i].Side()-ncorner);
                // If the edge already exists in mEdge, insert the face connect value. Create one otherwise.
                if(mEdge.find(con) != mEdge.end()){
                    mEdge[con].face_connect.insert(con_face);
                } else {
                    //create an auxiliary edge
                    EdgeFilter auxEdge;
                    mEdge.insert(std::make_pair(con,auxEdge));
                    mEdge[con].index = con;           
                    mEdge[con].face_connect.insert(con_face);
                }

                // Face data structure
                // If the face already exists in mFace, insert the edge connect value. Create one otherwise.
                if(mFace.find(con_face) != mFace.end()){
                    mFace[con_face].edge_connect.insert(con);
                } else {
                    //create an auxiliary edge
                    FaceFilter auxFace;
                    mFace.insert(std::make_pair(con_face,auxFace));
                    mFace[con_face].index = con;           
                    mFace[con_face].edge_connect.insert(con);
                }
            }
        }
    }

    //Complete the data structures   
    for (auto inode = 0; inode < mVertex.size(); inode++)
    {
        mVertex[inode].free_edges = mVertex[inode].edge_connect.size();
    }

    //Cheks if an edge has more than 2 vertices
    for(auto it = mEdge.cbegin(); it != mEdge.cend(); ++it){
        if (it->second.vertex_connect.size() != 2){
            std::cout << "Edge with more than 2 vertices!\n";
            DebugStop();
        }
    }
    for(auto it = freeEdges.cbegin(); it != freeEdges.cend(); ++it){
        if (it->second.size() != 2){
            std::cout << "Edge with more than 2 vertices!\n";
            DebugStop();
        }
    }
    
    //Cheks if a face has more than 3 edges
    for(auto it = mFace.cbegin(); it != mFace.cend(); ++it){
        if (it->second.edge_connect.size() != 3){
            std::cout << "Face with more than 3 edges!\n"; 
            DebugStop();
        }
    }

    //*************FINISH**************
    // PRINT THE DATA STRUCTURES
    // std::cout << "MEDGE.size() " << mEdge.size()<<std::endl;
    // for(auto it = mEdge.cbegin(); it != mEdge.cend(); ++it)
    // {
    //     std::cout <<  "MEDGE = " << it->first <<  ", faces = ";
    //     for (auto rem : it->second.face_connect)
    //     {
    //         std::cout << rem << " " ;
    //     }
    //     std::cout << std::endl;
    // }
    
    // std::cout << "mFace.size() " << mFace.size()<<std::endl;
    // for(auto it = mFace.cbegin(); it != mFace.cend(); ++it)
    // {
    //     std::cout <<  "mFace = " << it->first <<  ", faces = ";
    //     for (auto rem : it->second.edge_connect)
    //     {
    //         std::cout << rem << " " ;
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "mEdge.size() " << mEdge.size()<<std::endl;
    // for(auto it = mEdge.cbegin(); it != mEdge.cend(); ++it)
    // {
    //     std::cout <<  "mEdge = " << it->first <<  ", faces = ";
    //     for (auto rem : it->second.vertex_connect)
    //     {
    //         std::cout << rem << " " ;
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "freeEdges.size() " << freeEdges.size()<<std::endl;
    // for(auto it = freeEdges.cbegin(); it != freeEdges.cend(); ++it)
    // {
    //     std::cout <<  "freeEdges = " << it->first <<  ", free vertices = ";
    //     for (auto rem : it->second)
    //     {
    //         std::cout << rem << " " ;
    //     }
    //     std::cout << std::endl;
    // }
}

template <class TVar>
int64_t TPZHCurlEquationFilter<TVar>::ChooseNode(TPZVec<VertexFilter> &mVertex)
{
    int64_t nedges_min = 1e3;
    int64_t treatNode;
    for (int inode = 0; inode < mVertex.size(); inode++)
    {
        if (mVertex[inode].status) continue; // it has been already treated

        // The criterium is the node with min edges. TO BE DISCUSSED
        if (mVertex[inode].edge_connect.size() < nedges_min) {
            treatNode = inode;
            nedges_min = mVertex[inode].edge_connect.size();
        }
    }

    return treatNode;
}   

template <class TVar>
int64_t TPZHCurlEquationFilter<TVar>::ChooseEdge(int64_t &treatNode, TPZVec<VertexFilter> &mVertex, 
                                                 std::map<int,EdgeFilter> &mEdge)
{
    int64_t nfaces_min = 1e3;
    int64_t remEdge;
    bool flag = false;
    for (auto iedge : mVertex[treatNode].edge_connect)
    {
        // If the edge has been removed or is blocked, go to the next one
        if (mVertex[treatNode].status == EBlockedEdge) continue;
        if (mVertex[treatNode].status == ERemovedEdge) continue;

        // The criterium adopted is the edge with min number of faces... to be discussed
        if (mEdge[iedge].face_connect.size() < nfaces_min) {
            remEdge = iedge;
            nfaces_min = mEdge[iedge].face_connect.size();
            flag = true;
        }
    }

    if (!flag) DebugStop(); // No edge available to remove in this node

    return remEdge;
}

template <class TVar>
void TPZHCurlEquationFilter<TVar>::UpdateVertexFreeEdges(int64_t &treatNode, int64_t &remEdge, 
                                                         TPZVec<VertexFilter> &mVertex,
                                                         std::map<int,std::set<int>> &freeEdges)
{
    // Discover the two nodes connect to the chosen edge 
    TPZVec<int> twoNodes(2);
    twoNodes[0] = treatNode; //One is the treated node. Just need to find the other one
    for (int inode = 0; inode < mVertex.size(); inode++)
    {
        if (mVertex[inode].edge_connect.find(remEdge) != mVertex[inode].edge_connect.end()) {
            if (inode == treatNode) continue;
            twoNodes[1] = inode;
            break;
        }
    }
    //Reduce the number of free edges for the two nodes
    mVertex[twoNodes[0]].free_edges--;
    mVertex[twoNodes[1]].free_edges--;

    // if(mVertex[twoNodes[0]].free_edges < 0 || mVertex[twoNodes[1]].free_edges < 0) DebugStop();

    

}

template <class TVar>
void TPZHCurlEquationFilter<TVar>::UpdateEdgeAndFaceStatus(int64_t &remEdge, std::map<int,EdgeFilter> &mEdge, std::map<int,FaceFilter> &mFace)
{

    //************** MAYBE THERE IS SOMETHING MISSING HERE.**************
    for (auto iface : mEdge[remEdge].face_connect)
    {
        auto status = mFace[iface].status;
        if (status == EFreeFace){
            mFace[iface].status == ECryticalFace;
        } else if (status == ECryticalFace) {
            mFace[iface].status == EBlockedFace;
        }
    }

    // Updates blocked edges and face status
    for(auto it = mFace.begin(); it != mFace.end(); ++it)
    {
        int aux = 0;
        TPZVec<int> loc_edges(2); // the edges already removed - if >2 should give a DebugStop.
        //Loops over the edges of a face
        for (auto iedge : it->second.edge_connect)
        {
            if (mEdge[iedge].status == ERemovedEdge){
                loc_edges[aux] = iedge;
                aux++;
            }
        }
        // Now checks how many edges have been removed
        switch (aux)
        {
            case 0:
                it->second.status = EFreeFace;
                break;
            case 1:
                it->second.status = ECryticalFace;
                break;
            case 2:
                //Blocks the third edge
                for (auto iedge : it->second.edge_connect){
                    if (iedge != loc_edges[0] && iedge != loc_edges[1]){
                        mEdge[iedge].status = EBlockedEdge;
                        break;
                    }
                }
                it->second.status = EFortunateFace; // it has an edge blocked.
                break;
            default:
                DebugStop();
                break;
        }
    }
}


template <class TVar>
bool TPZHCurlEquationFilter<TVar>::FilterEdgeEquations(TPZAutoPointer<TPZCompMesh> cmesh,
                    TPZVec<int64_t> &activeEquations)
{

    const auto gmesh = cmesh->Reference();    
    const auto nnodes = gmesh->NNodes();

    // 1 - Vertex Data Structure
    TPZVec<VertexFilter> mVertex(gmesh->NNodes());
    // 2 - Edge Data Structure
    std::map<int,EdgeFilter> mEdge;
    // 3 - Face Data Structure
    std::map<int,FaceFilter> mFace;
    // 4 - The free edges and their corresponding free vertices
    std::map<int,std::set<int>> freeEdges;

    InitDataStructures(gmesh,mVertex,mEdge,mFace,freeEdges);

     /**
        The i-th is true if the i-th node has already been dealt with.
    */
    TPZVec<bool> done_vertices(nnodes, false);
    const int edgeDim{1};
    /* We need to keep ONE edge, globally speaking*/
    done_vertices[0] = true;

    /**
        Contains all the connects marked for removal. It is expected
        that all the connects are associated with edges.
    */
    // ATTENTION!! IT IS NOT BEEING CURRENTLY USED, BUT CAN BE UPDATED AT THE END OF THE ALGORITHM
    std::set<int64_t> removed_edges; 

    // while there is a vertex not treated
    auto count = std::count(done_vertices.begin(), done_vertices.end(), false);
    
    //ATENTION! Needs to perform an update before entering the loop since we want to keep node 0.
    //It is just to mark it as treated, and update the edge information

    while (count > 0)
    {    
        // Choose the node to be treated
        int64_t treatNode = ChooseNode(mVertex);

        // Choose the edge to be removed: the edge with min number of faces (but can be another criterium)
        int64_t remEdge = ChooseEdge(treatNode,mVertex,mEdge);
      
        // Update the node status
        mVertex[treatNode].status = true;
        done_vertices[treatNode] = true; // used just for the counter... can be changed/removed later
        mVertex[treatNode].removed_edge = remEdge;
        mEdge[remEdge].vertex_treated = treatNode;

        // Updates the free edges in a vertex
        UpdateVertexFreeEdges(treatNode,remEdge,mVertex,freeEdges);

        // Update the edge status 
        mEdge[remEdge].status = ERemovedEdge;

        // Updates the edges (if some needs to be blocked) and faces status.
        UpdateEdgeAndFaceStatus(remEdge,mEdge,mFace);

        
        //Update freeEdges structure
        freeEdges.erase(remEdge);
        for(auto it = freeEdges.begin(); it != freeEdges.end(); ++it)
        {
            for (auto iNode : it->second)
            {
                if (iNode == treatNode){
                    it->second.erase(iNode);
                    break;
                }
            }
        }
        
        // Update number of free edges:
        for(auto it = freeEdges.begin(); it != freeEdges.end(); ++it)
        {
            if (it->second.size() == 0)
            {
                freeEdges.erase(it->first);
            }
        }
        
        //Print for checking
        // std::cout << "freeEdges.size() " << freeEdges.size()<<std::endl;
        // for(auto it = freeEdges.cbegin(); it != freeEdges.cend(); ++it)
        // {
        //     std::cout <<  "freeEdges = " << it->first <<  ", free vertices = ";
        //     for (auto rem : it->second)
        //     {
        //         std::cout << rem << " " ;
        //     }
        //     std::cout << std::endl;
        // }

        //Count the number of nodes not treated
        count = std::count(done_vertices.begin(), done_vertices.end(), false);


    }//while
    

    for (size_t i = 0; i < nnodes; i++)
    {
        if (mVertex[i].status) removed_edges.insert(mVertex[i].removed_edge);       
    }
    

    std::cout << "Removed (" << removed_edges.size() << ") = ";
    for (auto rem : removed_edges)
    {
        std::cout << rem << " " ;
    }
    std::cout << std::endl;

    // std::cout << "No_remove (" << no_remove.size() << ") = ";
    // for (auto rem : no_remove)
    // {
    //     std::cout << rem << " " ;
    // }
    // std::cout << std::endl;

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
        for (auto rem : mVertex[i].edge_connect)
        {
            if (removed_edges.find(rem) != removed_edges.end()){
                std::cout << "{" << rem << "} " ;
                count++;
            } else {
                std::cout << rem << " " ;
            }
            
        }
        std::cout <<";             " << count << " edge(s) removed" << std::endl;
        // if (count == 0 && i != 0) std::cout << "PROBLEM! Node " << i << " did not removed an edge!\n" ;
        
    }
    
    // if(check_all_vertices == false) DebugStop();
    
    // bool check_edges_left{true};
    // for(auto iv = 0; iv < nnodes; iv++){
    //     const auto all_edges = vertex_edge_connects[iv];
    //     bool local_check{false};
    //     for(auto edge: all_edges){
    //     if(removed_edges.find(edge) == removed_edges.end()){
    //         local_check = true;
    //         break;
    //     }
    //     }
    //     check_edges_left = local_check && check_edges_left;
    // }

    // if (check_edges_left == false) DebugStop();
    
    activeEquations.Resize(0);
    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
        if (removed_edges.find(iCon) == removed_edges.end()) {
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
            if (removed_edges.find(con) != removed_edges.end()) {
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

    if (removed_edges.size() != nnodes-1){
        std::cout << "Removed " << removed_edges.size() << "/" << nnodes-1 << " allowed connects."  << std::endl;
        DebugStop();
    }
    

    return 0;
}

template class TPZHCurlEquationFilter<STATE>;
// template class TPZHCurlEquationFilter<CSTATE>;


