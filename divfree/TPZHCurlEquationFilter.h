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
    // The data structures

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
    
    
    void CreateFilterDataStructure(TPZGeoMesh* gmesh, const int &edgeDim, TPZVec<std::set<int>> &vertex_edge_connects);

};
#endif