

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

class TPZMHMCompMeshCreator {

public:
    TPZMHMGeoMeshCreator fGeoMeshCreator;
    int fAvPresLevel;
    int fDistFluxLevel;
    TPZHDivApproxCreator *hDivCreator;

public:
    TPZMHMCompMeshCreator(TPZMHMGeoMeshCreator &mhm_gcreator, TPZAutoPointer<TPZGeoMesh> &gmesh);

    ~TPZMHMCompMeshCreator(){};
    
    TPZMultiphysicsCompMesh * BuildMultiphysicsCMesh(int pOrder_vol, int pOrder_skel, TPZAutoPointer<TPZGeoMesh> &gmesh, TPZAnalyticSolution &analytic);

    void InsertMaterialObjects(TPZAnalyticSolution &analytic);

    void PutinSubstructures(TPZCompMesh &cmesh);

    void CondenseElements(TPZCompMesh &cmesh);

    
};

#endif