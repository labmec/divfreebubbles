

#ifndef MHM_GEOMESHUTILS
#define MHM_GEOMESHUTILS

#include "pzgmesh.h"

class TPZMHMGeoMeshCreator {

public:
    int fSkeletonMatId;
    int fDomainMatId;
    int fBC1;
    int fBC2;
    int fBC3;
    int fBC4;
    int fNumSubGrids;
    std::set<int> fBoundMatId;

    TPZVec<int64_t> fElementPartition;

public:
    TPZMHMGeoMeshCreator() = default;

    ~TPZMHMGeoMeshCreator(){};

    void SetBCMatId(int bcMatId){
        fBoundMatId.insert(bcMatId);
    }
    
    void CreateSkeleton(TPZAutoPointer<TPZGeoMesh> &gmesh);

    void CreateSubGrids(TPZAutoPointer<TPZGeoMesh> &gmesh);

    void RefineSubGrids(TPZAutoPointer<TPZGeoMesh> &gmesh);

    void RefineSkeleton(TPZAutoPointer<TPZGeoMesh> &gmesh);

    void AddBoundaryElements(TPZAutoPointer<TPZGeoMesh> &gmesh);

    void CreateTriangleElements(TPZAutoPointer<TPZGeoMesh> gmesh, std::map<int,int> &matmap, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenter);

    
};

#endif