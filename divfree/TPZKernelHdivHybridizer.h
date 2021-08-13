//
// Created by Jeferson Fernandes on 11/08/21.
//

#ifndef TPZKERNELHDIV_HYBRIDIZER_H
#define TPZKERNELHDIV_HYBRIDIZER_H

// TODO add doc

#include <iostream>
#include <map>
#include "pzstack.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZMultiphysicsCompMesh.h"

//#include "pzgeoelrefless.h"

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZMaterial;
class TPZCompElSide;
class TPZInterpolatedElement;
class TPZMultiphysicsCompMesh;
class TPZGeoEl;

template<class T, int N>
class TPZStack;

struct TPZKernelHdivHybridizer {

    // Wrap geometric elements material id
    int fEWrap = -10;
    // Hybridized pressure geometric elements material id
    int fEPressureHyb = -9;
    // Interface elements material id
    int fEInterface = -8;
    // Point-Wrap elements material id
    int fEPont = -7;
    // Domain material id
    int fEDomain = 1;
    // number of state variables
    int fNState = 1;

    /**
     * @brief Default constructor
     */
    TPZKernelHdivHybridizer() = default;

    //Creates Wrap, Interface, Points and Hybridized pressure geometric elements 
    //for both domain and Boundary Condition Hybridizations. The set matIdBC must contains the 
    //boundary material id's to be hybridized
    void CreateWrapElements(TPZGeoMesh *geomesh, std::set<int> &matIdBC, bool domainHyb);

    void SetPeriferalMaterialIds(int Wrap, int Lagrange, int Interface, int Point, int Domain)
    {
        fEWrap = Wrap;
        fEPressureHyb = Lagrange;
        fEInterface = Interface;
        fEPont = Point;
        fEDomain = Domain;

    }

    void SemiHybridizeFlux(TPZCompMesh *cmesh, std::set<int> &matBCId);

    void CreateMultiphysicsInterfaceElements(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh, TPZVec<TPZCompMesh *> &meshvector, std::set<int> &matIdNeumann);

    void GroupAndCondenseElements(TPZMultiphysicsCompMesh *cmesh);

    void AssociateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup, TPZVec<int64_t> &elementgroup2);

};

#endif //TPZKERNELHDIVHYBRIDIZER_H
