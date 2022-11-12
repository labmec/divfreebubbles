

#ifndef MHM_COMPMESHUTILS
#define MHM_COMPMESHUTILS

#include "pzcmesh.h"
#include "pzgmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMHMGeoMeshCreator.h"
#include "TPZAnalyticSolution.h"
#include "TPZHDivApproxCreator.h"

class TPZCompMesh;
class TPZGeoMesh;
class TPZMultiphysicsCompMesh;
class TPZCompEl;
class TPZGeoElSide;

class TPZMHMHDivApproxCreator : public TPZHDivApproxCreator {

public:
    TPZMHMGeoMeshCreator fGeoMeshCreator;
    int fAvPresLevel;
    int fDistFluxLevel;
    int fPOrderSkeleton;

public:
    TPZMHMHDivApproxCreator(TPZMHMGeoMeshCreator &mhm_gcreator, TPZAutoPointer<TPZGeoMesh> &gmesh);

    ~TPZMHMHDivApproxCreator(){};
    
    TPZMultiphysicsCompMesh * BuildMultiphysicsCMesh();

    void InsertMaterialObjects(TPZAnalyticSolution &analytic);

    void PutinSubstructures(TPZCompMesh &cmesh);

    void CondenseElements(TPZCompMesh &cmesh);

    /// Set skeleton default polynomial order
    void SetPOrderSkeleton(const int ord){
        fPOrderSkeleton = ord;
    }

    /// Get skeleton default polynomial order
    int &GetPOrderSkeleton(){return fPOrderSkeleton;}

    
};

#endif