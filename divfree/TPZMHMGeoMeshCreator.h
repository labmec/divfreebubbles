

#ifndef MHM_GEOMESHUTILS
#define MHM_GEOMESHUTILS

#include "pzgmesh.h"

class TPZMHMGeoMeshCreator {

public:
    int fSkeletonMatId;
    int fDomainMatId;
    std::set<int> fBoundMatId;
    int fNumSubGrids;

    TPZVec<int64_t> fElementPartition;

public:
    TPZMHMGeoMeshCreator() = default;

    ~TPZMHMGeoMeshCreator(){};
    
    void CreateSkeleton(TPZGeoMesh *gmesh);

    void CreateSubGrids(TPZGeoMesh *gmesh);

    void RefineSubGrids(TPZGeoMesh *gmesh);

    void RefineSkeleton(TPZGeoMesh *gmesh);
};

#endif