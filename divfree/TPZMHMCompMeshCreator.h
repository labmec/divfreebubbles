

#ifndef MHM_COMPMESHUTILS
#define MHM_COMPMESHUTILS

#include "pzcmesh.h"
#include "pzgmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMHMGeoMeshCreator.h"

class TPZMHMCompMeshCreator {

public:
    TPZMHMGeoMeshCreator fGeoMeshCreator;
    int fAvPresLevel;
    int fDistFluxLevel;

public:
    TPZMHMCompMeshCreator(TPZMHMGeoMeshCreator &mhm_gcreator);

    ~TPZMHMCompMeshCreator(){};
    
    TPZMultiphysicsCompMesh * BuildMultiphysicsCMesh(int pOrder_vol, int pOrder_skel, TPZGeoMesh * gmesh);

    /// generate a mesh with HDiv elements
    TPZCompMesh *GenerateFluxMesh(TPZGeoMesh* gmesh, int pOrder_vol, int pOrder_skel);

    /// generate a mesh with L2 elements
    TPZCompMesh *GeneratePressureMesh(TPZGeoMesh* gmesh, int pOrder);

    /// generate a mesh with constant elements
    TPZCompMesh *GenerateConstantMesh(TPZGeoMesh* gmesh, int level);

    void InsertMaterialObjects(TPZMultiphysicsCompMesh &cmesh);

    void PutinSubstructures(TPZCompMesh &cmesh);

    void CondenseElements(TPZCompMesh &cmesh);
};

#endif