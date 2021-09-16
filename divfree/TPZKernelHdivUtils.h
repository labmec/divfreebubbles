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
#include <TPZVTKGeoMesh.h>

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


};
#endif