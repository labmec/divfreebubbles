//
// Created by Jeferson Fernandes on 11/08/21.
//

#ifndef KERNELHDIV_UTILS_H
#define KERNELHDIV_UTILS_H

// TODO add doc

#include <iostream>
#include <map>
#include "pzstack.h"
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

struct TPZKernelHdivUtils {

    /**
     * @brief Default constructor
     */
    TPZKernelHdivUtils() = default;


    void PrintCMeshConnects(TPZCompMesh *cmesh);

    void PrintGeoMesh(TPZGeoMesh *geomesh);



};
#endif