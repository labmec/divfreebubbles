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
    // Edge Remove material id
    int fEEdgeRemove = -6;
    // number of state variables
    int fNState = 1;

    /**
     * @brief Default constructor
     */
    TPZKernelHdivHybridizer() = default;

    /**
     * @brief This method adapts a given geometric mesh to both domain and BC hybridizations.
     * It creates Wrap, Interface, Points and Hybridized pressure geometric elements 
     * for both Hybridization cases. 
     * 
     * @param geomesh = the given geometric mesh
     * @param matIdBC = a std::set containing all BC' material id's to be hybridized
     * @param domainHyb = true if the domain is also hybridized, false otherwise
     */
    void CreateWrapElements(TPZGeoMesh *geomesh, std::set<int> &matIdBC, bool domainHyb, HDivFamily &hdivfam);

    /**
     * @brief Set the material id's needed to separate the geometric elements and perform the hybridation.
     * 
     * @param Wrap 
     * @param Lagrange 
     * @param Interface 
     * @param Point 
     * @param Domain 
     */
    void SetMaterialIds(int Wrap, int Lagrange, int Interface, int Point, int Domain)
    {
        fEWrap = Wrap;
        fEPressureHyb = Lagrange;
        fEInterface = Interface;
        fEPont = Point;
        fEDomain = Domain;

    }
    void SetEdgeRemove(int edge){
        fEEdgeRemove = edge;
    }
    
    /**
     * @brief Updates the wrap element's neighbour to perform a domain semi-hybridization, i.e, 
     * an hybridization with constant pressure through elements
     * 
     * @param cmesh = the computational mesh
     * @param matBCId = BC material id's
     */
    void SemiHybridizeFlux(TPZCompMesh *cmesh, std::set<int> &matBCId);

    /**
     * @brief Updates the approximation order for the BC hybridization Lagrange multipliers.
     * As in the semi-hybrid case constant pressure through elements is adopted to the domain,
     * the BC Lagrange multipliers order needs to be updated.
     * 
     * @param cmesh = the computational mesh
     * @param pOrder = the polynomial order to be updated
     * @param matBCId = BC material id's 
     */
    void SemiHybridizePressure(TPZCompMesh *cmesh, int pOrder, std::set<int> &matBCId);
    void EdgeRemove(TPZCompMesh *cmesh);


    /**
     * @brief Create the Multiphysics Interface Elements after all wrap, interface 
     * and Lagrange multipliers are set
     * 
     * @param cmesh 
     * @param gmesh 
     * @param meshvector 
     * @param matIdNeumann = hybridized BC material is's
     */
    void CreateMultiphysicsInterfaceElements(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh *gmesh, TPZVec<TPZCompMesh *> &meshvector, std::set<int> &matIdBCHyb);

    /**
     * @brief Groups and performs static condensation to elements
     * 
     * @param cmesh 
     * @param matIdBC 
     */
    
    void GroupAndCondenseElements(TPZMultiphysicsCompMesh *cmesh, std::set<int> &matIdBC);
    void GroupAndCondenseElementsDuplicatedConnects(TPZMultiphysicsCompMesh *cmesh, std::set<int> &matIdBC);

    /**
     * @brief Associate element connects to be condensed
     * 
     * @param cmesh 
     * @param elementgroup 
     * @param matIdBC 
     */
    void AssociateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup, std::set<int> &matIdBC);
    void AssociateElementsToMatId(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup, std::set<int> &matIdBC);
    void AssociateElementsDuplicatedConnects(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup, std::set<int> &matIdBC);

};

#endif //TPZKERNELHDIVHYBRIDIZER_H
