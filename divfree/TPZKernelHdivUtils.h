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
    void PrintCMeshConnects(TPZCompMesh *cmesh);

    void PrintGeoMesh(TPZGeoMesh *geomesh);

    void PrintCompMesh(TPZCompMesh *cmesh, std::string &file_name);

    void SolveProblemIterative(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

    void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh); 

};
#endif