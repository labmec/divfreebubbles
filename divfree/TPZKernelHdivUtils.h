//
// Created by Jeferson Fernandes on 11/08/21.
//

#ifndef KERNELHDIV_UTILS_H
#define KERNELHDIV_UTILS_H

// TODO add doc

#include <iostream>
#include <map>
#include "pzstack.h"
#include "TPZLinearAnalysis.h"

//#include "pzgeoelrefless.h"

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZGeoMesh;
class TPZMaterial;
class TPZCompElSide;
class TPZInterpolatedElement;
class TPZMultiphysicsCompMesh;
class TPZGeoEl;

template<class T, int N>
class TPZStack;

template <class TVar>
class TPZKernelHdivUtils {

public:
    /**
     * @brief Prints the computational mesh information (specially the connects)
     * 
     * @param cmesh 
     */
    void PrintCMeshConnects(TPZCompMesh *cmesh);

    /**
     * @brief Prints geometric mesh information (specially neighbours)
     * 
     * @param geomesh 
     */
    void PrintGeoMesh(TPZGeoMesh *geomesh);

    /**
     * @brief Prints computational mesh information to a file
     * 
     * @param cmesh 
     * @param file_name 
     */
    void PrintCompMesh(TPZCompMesh *cmesh, std::string &file_name);

    /**
     * @brief Solves an algebraic system by means of an iterative method
     * 
     * @param an 
     * @param cmesh 
     */
    void SolveProblemIterative(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

    /**
     * @brief Solves an algebraic system by means of a direct method
     * 
     * @param an 
     * @param cmesh 
     */
    void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh); 

    /**
     * @brief Prints the results of a multiphysics mesh to a .vtk file
     * 
     * @param meshvector 
     * @param an 
     * @param cmesh 
     */
    void PrintResultsMultiphysics(TPZVec<TPZCompMesh *> &meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh);

    /**
     * @brief Util to compute the solution error
     * 
     * @param an 
     * @param anPostProcessFile 
     */
    void ComputeError(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);

    /**
         @brief Removes some equations associated with edges to ensure that
        the gradient of the lowest order H1 functions cannot be represented.
        @param[in] cmesh Computational mesh.
        @param[out] indices of all remaining equations.
        @return 0 if no errors were detected, 1 if a vertex was left untreated,
        2 if a vertex had all the adjacent edges removed.
    */
    bool FilterEdgeEquations(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<int64_t> &activeEquations);

};
#endif