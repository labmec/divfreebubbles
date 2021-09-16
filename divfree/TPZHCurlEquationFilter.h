//
// Created by Jeferson Fernandes on 16/09/21.
//

#ifndef HCURL_EQUATION_FILTER_H
#define HCURL_EQUATION_FILTER_H

// TODO add doc

#include <iostream>
#include <map>
#include "pzstack.h"
#include <TPZVTKGeoMesh.h>

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZGeoMesh;
class TPZCompElSide;
class TPZInterpolatedElement;
class TPZGeoEl;

template<class T, int N>
class TPZStack;

template <class TVar>
class TPZHCurlEquationFilter {

private:
    // DATA STRUCTURES
    // 1 - VERTEX DATA STRUCTURE
    struct Vertex
    {   
        // The edges connected to a vertex
        std::set<int> edge_connect;
        // The number of free edges, i.e, not removed in a vertex
        int64_t free_edges;
    };

    // 2 - EDGE DATA STRUCTURE
    struct Edge
    {
        //The edge index
        int index;
        //The vertices connected to an edge
        std::set<int> vertex_connect;
        //The faces connected to an edge
        std::set<int> face_connect;
        //The edge status: 0 - free to be removed; 1 - removed; 2 - blocked
        int status = 0;
        //Number of faces removed from the edge
        int faces_removed = 0;
    };
    
    // 3 - FACE DATA STRUCTURE
    struct Face
    {
        //The face index
        int index;

        //The connects of a face
        // TPZVec<int> edge_connect(const int size = 3);
        std::set<int> edge_connect;
        
        //Face status:
        //  0 - free: 0 edges removed
        //  1 - crytical: 1 edge removed
        //  2 - blocked: 2 edges removed
        //  3 - fortunate: 1 edge blocked
        int status = 0;    
    };

  
public:
    /**
         @brief Removes some equations associated with edges to ensure that
        the gradient of the lowest order H1 functions cannot be represented.
        @param[in] cmesh Computational mesh.
        @param[out] indices of all remaining equations.
        @return 0 if no errors were detected, 1 if a vertex was left untreated,
        2 if a vertex had all the adjacent edges removed.
    */
    bool FilterEdgeEquations(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<int64_t> &activeEquations);
    
    void InitDataStructures(TPZGeoMesh* gmesh, TPZVec<Vertex> &mVertex, std::map<int,Edge> &mEdge,std::map<int,Face> &mFace);


};
#endif