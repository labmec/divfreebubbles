

#ifndef MHM_COMPMESHUTILS
#define MHM_COMPMESHUTILS

#include "pzcmesh.h"
#include "pzgmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMHMGeoMeshCreator.h"
#include "TPZAnalyticSolution.h"
#include "TPZHDivApproxSpaceCreator.h"

class TPZMHMCompMeshCreator {

public:
    TPZMHMGeoMeshCreator fGeoMeshCreator;
    int fAvPresLevel;
    int fDistFluxLevel;
    bool fDuplConnects = false;
    TPZHDivApproxSpaceCreator<double> *fHDivCreator;
    HDivFamily fHdivFamily;

public:
    TPZMHMCompMeshCreator(TPZMHMGeoMeshCreator &mhm_gcreator, HDivFamily hdivfam = HDivFamily::EHDivStandard);

    ~TPZMHMCompMeshCreator(){};
    
    TPZMultiphysicsCompMesh * BuildMultiphysicsCMesh(int pOrder_vol, int pOrder_skel, TPZAutoPointer<TPZGeoMesh> &gmesh, TPZAnalyticSolution &analytic);

    /// generate a mesh with HDiv elements
    TPZCompMesh *GenerateFluxMesh(TPZAutoPointer<TPZGeoMesh> &gmesh, int pOrder_vol, int pOrder_skel);

    /// generate a mesh with L2 elements
    TPZCompMesh *GeneratePressureMesh(TPZAutoPointer<TPZGeoMesh> &gmesh, int pOrder);

    /// generate a mesh with constant elements
    TPZCompMesh *GenerateConstantMesh(TPZAutoPointer<TPZGeoMesh> &gmesh, int level);

    void InsertMaterialObjects(TPZMultiphysicsCompMesh &cmesh,TPZAnalyticSolution &analytic);

    void PutinSubstructures(TPZCompMesh &cmesh);

    void CondenseElements(TPZCompMesh &cmesh);

    void DuplicateConnects(){fDuplConnects = true;}

    void CreateFluxDuplicatedConnects(TPZCompMesh *cmesh);

};

#endif