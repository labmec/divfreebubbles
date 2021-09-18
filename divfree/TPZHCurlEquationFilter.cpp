//
// A filter to HCurl spaces.
//

#include "TPZHCurlEquationFilter.h"

#include "pzcmesh.h"
#include "pzgmesh.h"
#include <TPZVTKGeoMesh.h>

/* 
    Initialize the data structures used to filter the equations
*/
template <class TVar>
void TPZHCurlEquationFilter<TVar>::InitDataStructures(TPZGeoMesh* gmesh)
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
        const auto nsides = gel->NSides();

        for(auto ie = firstEdge; ie < firstFace; ie++){
            // Vertex data structure            
            TPZGeoElSide edge_(gel,ie);
            const auto conIndex = edge_.Element()->Reference()->ConnectIndex(ie-firstEdge);
            
            //Store the connects of a given vertex 
            const int64_t v1 = edge_.SideNodeIndex(0);
            const int64_t v2 = edge_.SideNodeIndex(1);
            
            mVertex[v1].edge_connect.insert(conIndex);
            mVertex[v2].edge_connect.insert(conIndex);

            mEdge[conIndex].vertex_connect.insert(v1);
            mEdge[conIndex].vertex_connect.insert(v2);

            TPZStack<TPZGeoElSide> gelsides;
            gel->AllHigherDimensionSides(ie,2,gelsides);

            for (int i = 0; i < gelsides.size(); i++)
            {   
                // Edge data structure
                const auto con_face = edge_.Element()->Reference()->ConnectIndex(gelsides[i].Side()-ncorner);
                // If the edge already exists in mEdge, insert the face connect value. Create one otherwise.
                if(mEdge.find(conIndex) != mEdge.end()){
                    mEdge[conIndex].face_connect.insert(con_face);
                } else {
                    //create an auxiliary edge
                    EdgeFilter auxEdge;
                    mEdge.insert(std::make_pair(conIndex,auxEdge));
                    mEdge[conIndex].index = conIndex;           
                    mEdge[conIndex].face_connect.insert(con_face);
                }

                // Face data structure
                // If the face already exists in mFace, insert the edge connect value. Create one otherwise.
                if(mFace.find(con_face) != mFace.end()){
                    mFace[con_face].edge_connect.insert(conIndex);
                } else {
                    //create an auxiliary edge
                    FaceFilter auxFace;
                    mFace.insert(std::make_pair(con_face,auxFace));
                    mFace[con_face].index = conIndex;           
                    mFace[con_face].edge_connect.insert(conIndex);
                }
            }
        }
    }

    //Complete the data structures
    for(auto it = mVertex.begin(); it != mVertex.end(); ++it)
    // for (auto inode = 0; inode < mVertex.size(); inode++)
    {
        int size_edges = it->second.edge_connect.size();
        it->second.free_edges = size_edges;
        freeEdgesToNodes[size_edges].insert(it->first);
    }

    //Cheks if an edge has more than 2 vertices
    for(auto it = mEdge.cbegin(); it != mEdge.cend(); ++it){
        if (it->second.vertex_connect.size() != 2){
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

    // //*************FINISH**************
    // // PRINT THE DATA STRUCTURES
    // std::cout << "vertex.size() " << mVertex.size()<<std::endl;
    // for(auto it = mVertex.cbegin(); it != mVertex.cend(); ++it)
    // {
    //     std::cout <<  "mVertex = " << it->first <<  ", faces = ";
    //     for (auto rem : it->second.edge_connect)
    //     {
    //         std::cout << rem << " " ;
    //     }
    //     std::cout << std::endl;
    // }
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

    // std::cout << "freeEdges.size() " << freeEdgesToNodes.size()<<std::endl;
    // for(auto it = freeEdgesToNodes.cbegin(); it != freeEdgesToNodes.cend(); ++it)
    // {
    //     std::cout <<  "freeEdgesToNodes = " << it->first <<  ", free nodes = ";
    //     for (auto rem : it->second)
    //     {
    //         std::cout << rem << " " ;
    //     }
    //     std::cout << std::endl;
    // }
}


/* 
    Chooses a node to be treated 
*/
template <class TVar>
int64_t TPZHCurlEquationFilter<TVar>::ChooseNode()
{

    // The first one in the list, i.e., the node with the lowest number of free edges       
    int64_t treatNode = *freeEdgesToNodes.begin()->second.begin();
    freeEdgesToNodes.begin()->second.erase(treatNode);
    
    return treatNode;
}

/*  
    Chooses an edge to be removed in the node chosen to be treated
*/
template <class TVar>
int64_t TPZHCurlEquationFilter<TVar>::ChooseEdge(int64_t &treatNode)
{
    int64_t nfaces_min = 1e3;
    int64_t remEdge;
    bool flag = false;
    for (auto iedge : mVertex[treatNode].edge_connect)
    {
        // If the edge has been removed or is blocked, go to the next one
        if (mEdge[iedge].status == EBlockedEdge) continue;
        if (mEdge[iedge].status == ERemovedEdge) continue;

        // The criterium adopted is the edge with min number of faces... can be changed
        if (mEdge[iedge].face_connect.size() < nfaces_min) {
            remEdge = iedge;
            nfaces_min = mEdge[iedge].face_connect.size();
            flag = true;
        }
    }

    if (!flag) DebugStop(); // No edge available to remove in this node

    return remEdge;
}

/* 
    Insert definition
*/
template <class TVar>
void TPZHCurlEquationFilter<TVar>::UpdateVertexFreeEdges(int64_t &treatNode, int64_t &remEdge)
{
 
    // Gets the two nodes connect to the chosen edge 
    TPZVec<int> twoNodes(2);
    int64_t diffNode;
    int count = 0;
    for (auto inode : mEdge[remEdge].vertex_connect){
        twoNodes[count] = inode;
        count++;
        if (inode != treatNode) diffNode = inode;
    }

    // The Vertex different from treatNode losts a freeEdge
    for(auto it = freeEdgesToNodes.begin(); it != freeEdgesToNodes.end(); it++)
    {
        if (it->second.find(diffNode) != it->second.end())
        {
            //Found diffNode in freeEdgesToNodes, then it losts a freeEdge.
            it->second.erase(diffNode);
            freeEdgesToNodes[it->first - 1].insert(diffNode);
            if (it->first - 1 == 0){
                std::cout << "No available edge for node " << diffNode << std::endl;
                DebugStop();
            }
            break;
        }
    }

    //Reduce the number of free edges for the two nodes
    mVertex[twoNodes[0]].free_edges--;
    mVertex[twoNodes[1]].free_edges--;

    if(mVertex[twoNodes[0]].free_edges < 0 || mVertex[twoNodes[1]].free_edges < 0) DebugStop();

}

template <class TVar>
void TPZHCurlEquationFilter<TVar>::UpdateEdgeAndFaceStatus(int64_t &remEdge)
{
    std::set<int> blockedEdges;

    //Look at the faces where the removed edge is a part of and updates its status
    for (auto iface : mEdge[remEdge].face_connect)
    {
        auto status = mFace[iface].status;
        if (status == EFreeFace){
            mFace[iface].status = ECryticalFace;
        } else if (status == ECryticalFace || status == EFortunateFace) { // here there is a ?
            mFace[iface].status = EBlockedFace;

            //Loops over the edges of a face. If it is a free edge, then block it
            for (auto iedge : mFace[iface].edge_connect){
                if (mEdge[iedge].status == EFreeEdge){
                    mEdge[iedge].status = EBlockedEdge;
                    blockedEdges.insert(iedge);
                    break;
                }
            }
        } else if (status == EBlockedFace){
            DebugStop();
        } 
    }

    // Now we can update all the faces in the mesh, as we can have some blocked 
    // edges and it will produce fortunate ones.
    for (auto iedge : blockedEdges)
    {
        for (auto iface : mEdge[iedge].face_connect){
            mFace[iface].status = EFortunateFace;
        }
    }    

    // Then update freeEdgetoNodes structure. - Loops over the edges.
    // The Vertex different from treatNode losts a freeEdge
    for(auto it = freeEdgesToNodes.begin(); it != freeEdgesToNodes.end(); it++)
    {
        for (auto node = mVertex.begin(); node != mVertex.end(); node++){
            if (node->second.status) continue;
            //Count the number of free edges.
            int count = 0;
            for (auto iedges : node->second.edge_connect){
                if (mEdge[iedges].status == EFreeEdge) count++;
            }
            //If the number of free edges changed, then update freeEdgesToNodes
            if (count != it->first){
                it->second.erase(node->first);
                freeEdgesToNodes[count].insert(node->first);
                if (count == 0){
                    std::cout << "No available edge for node " << node->first << std::endl;
                    std::cout << "Edges of this node: ";
                    for (auto iedges : node->second.edge_connect){
                        std::cout << iedges << " " ;
                    }
                    DebugStop();
                }
            }
        }
    }

}


template <class TVar>
bool TPZHCurlEquationFilter<TVar>::FilterEdgeEquations(TPZAutoPointer<TPZCompMesh> cmesh,
                    TPZVec<int64_t> &activeEquations)
{

    cmesh->LoadReferences();
    const auto gmesh = cmesh->Reference();    
    const auto nnodes = gmesh->NNodes();
    
    //Initialize the data structures
    InitDataStructures(gmesh);
    
    //Needs to perform an update before entering the loop since we want to keep one node.
    //One just mark the first one as treated, and update the free edges information
    mVertex.begin()->second.status = true;
    for(auto it = freeEdgesToNodes.begin(); it != freeEdgesToNodes.end(); ++it){
        if (it->second.find(mVertex.begin()->first) != it->second.end()) {
            it->second.erase(mVertex.begin()->first);
            break;
        }
    }

    while (freeEdgesToNodes.size() > 0)
    {    
        // First - Choose the node to be treated: the first one in the list, 
        // which corresponds to the one with the lowest number of free edges.
        int64_t treatNode = ChooseNode();

        // Choose the edge to be removed: the edge with min number of faces (but can be another criterium)
        // Pegar primeiro mais livre possível, depois uma face que esteja ok. uma crítica, verificar se está bloqueado
        int64_t remEdge = ChooseEdge(treatNode);
      
        // Update the node status and the edge related to it
        mVertex[treatNode].status = true;
        mVertex[treatNode].removed_edge = remEdge;
        mEdge[remEdge].vertex_treated = treatNode;

        // Updates the free edges in a vertex
        UpdateVertexFreeEdges(treatNode,remEdge);
        
        // Update the edge status 
        mEdge[remEdge].status = ERemovedEdge;

        // Updates the edges (if some needs to be blocked) and faces status.
        UpdateEdgeAndFaceStatus(remEdge);

        //Clears the structure with 0 nodes in some free edges
        for (auto it = freeEdgesToNodes.cbegin(); it != freeEdgesToNodes.cend();) {
            if (it->second.size() == 0) {
                freeEdgesToNodes.erase(it++);    // or "it = m.erase(it)" since C++11
            } else {
                ++it;
            }
        }

        //Check face
        for(auto it = mFace.cbegin(); it != mFace.cend(); ++it){
            int aux = 0; 
            for (auto iedge : it->second.edge_connect){
                if (mEdge[iedge].status == ERemovedEdge) aux++;
            }
            if (aux == 3){
                std::cout << "Problem with face " << it->first << std::endl;
                DebugStop();
            }
        }
    
    }//while
    
    std::set<int64_t> removed_edges; 
    TPZVec<bool> done_vertices(nnodes, false);

    for (size_t i = 0; i < nnodes; i++)
    {
        if (mVertex[i].status) {
            done_vertices[i] = true;
            removed_edges.insert(mVertex[i].removed_edge);
        }
    }
        
    bool check_edges_left{true};
    for(auto iv = 0; iv < nnodes; iv++){
        const auto all_edges = mVertex[iv].edge_connect;
        bool local_check{false};
        for(auto edge: all_edges){
        if(removed_edges.find(edge) == removed_edges.end()){
            local_check = true;
            break;
        }
        }
        check_edges_left = local_check && check_edges_left;
    }

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
            if (local_removed.size() == 0) continue;
            // std::cout << "gel = " << gel->Index() << ", removed connects = ";
            // for (auto rem : local_removed)
            // {
            //     std::cout << rem << ", " ;
            // }
            // Checks of coplanar edges
            if ((local_index[0] && local_index[1] && local_index[2]) ||
                (local_index[0] && local_index[3] && local_index[4]) ||
                (local_index[2] && local_index[3] && local_index[5]) ||
                (local_index[1] && local_index[5] && local_index[4])) 
                std::cout << "Problem with this element: coplanar edges removed!" ;
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


