
#include "TPZMHMGeoMeshCreator.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"
#include "TPZRefPatternDataBase.h"


//Creates the skeleton geometric elements
void TPZMHMGeoMeshCreator::CreateSkeleton(TPZGeoMesh *gmesh){

    for (int64_t iEl = 0; iEl < gmesh->NElements(); iEl++)
    {
        int dim = gmesh->Dimension();
        TPZGeoEl* gel = gmesh->ElementVec()[iEl];
        if (!gel || gel->Dimension() != dim) continue;

        fNumSubGrids++;
        int nSides = gel->NSides();
        int nFacets = gel->NSides(dim-1);

        for (int iFace = nSides-1-nFacets; iFace < nSides-1; iFace++)
        {
            TPZGeoElSide gelside(gel,iFace);
            TPZGeoElSide neighbour = gelside.Neighbour();

            if (neighbour.Element()->Dimension() == dim-1) {
                // neighbour.Element()->ResetReference();
                continue;
            }

            TPZGeoElBC gelSkeleton(gelside, fSkeletonMatId);

        }
    }  
};

void TPZMHMGeoMeshCreator::CreateSubGrids(TPZGeoMesh *gmesh){

    //This method should always be called after CreateSkeleton
    fElementPartition.resize(gmesh->NElements());
    fElementPartition.Fill(-1);

    for (int64_t iEl = 0; iEl < gmesh->NElements(); iEl++)
    {
        int dim = gmesh->Dimension();
        TPZGeoEl* gel = gmesh->ElementVec()[iEl];
        if (gel->MaterialId() == fSkeletonMatId) continue;
    
        fElementPartition[iEl] = gel->Index();
    }  

}

void TPZMHMGeoMeshCreator::RefineSubGrids(TPZGeoMesh *gmesh){


    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t> extpartition(fElementPartition);
    extpartition.Expand(nel*3);
    static int once = 1;
    int dim = gmesh->Dimension();
    if(once) {
        switch (dim)
        {
        case 2:
            gRefDBase.InitializeUniformRefPattern(EOned);    
            gRefDBase.InitializeUniformRefPattern(ETriangle);
            gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
            break;
        
        case 3:
            gRefDBase.InitializeUniformRefPattern(ETriangle);
            gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
            gRefDBase.InitializeUniformRefPattern(ETetraedro);
            gRefDBase.InitializeUniformRefPattern(ECube);
            break;

        default:
            DebugStop();
        }
        once = 0;
    }
    
    for (int64_t el = 0; el < nel; el++) {

        TPZGeoEl *gel = gmesh->Element(el);
        int gelmat = gel->MaterialId();
        bool boundmat = false;

        if (fBoundMatId.find(gelmat) != fBoundMatId.end()) boundmat = true;

        if(boundmat)
        {
            TPZGeoElSide gelside(gel);
            for(TPZGeoElSide neigh = gelside.Neighbour(); neigh != gelside; neigh++)
            {
                if(neigh.Element()->Dimension() == dim)
                {
                    fElementPartition[el] = fElementPartition[neigh.Element()->Index()];
                    break;
                }
            }
            auto skel = gelside.HasNeighbour(fSkeletonMatId);
            if(skel)
            {
                skel.Element()->RemoveConnectivities();
                delete skel.Element();
            }
        }

        // if(gel->Type() != ETriangle && !boundmat) continue;
        if(gel->Dimension() != dim && !boundmat) continue;
        
        auto partition = fElementPartition[el];
        
        if(partition < 0) DebugStop();

        TPZManVector<TPZGeoEl *> subels;
        gel->Divide(subels);
       
        
        int nsubels = subels.size();
        for (int sub=0; sub<nsubels; sub++) {
            int64_t index = subels[sub]->Index();
            if(extpartition.size() < index+1) extpartition.resize(index+1);
            extpartition[index] = partition;
        }
    }
    fElementPartition = extpartition;

}

void TPZMHMGeoMeshCreator::RefineSkeleton(TPZGeoMesh *gmesh){

    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t> extpartition(fElementPartition);
    extpartition.Expand(nel*3);
    
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;
        int gelmat = gel->MaterialId();
        if(gelmat != fSkeletonMatId) continue;
        auto partition = fElementPartition[el];
        if(partition >= 0) DebugStop();
        TPZManVector<TPZGeoEl *> subels(4);
        gel->Divide(subels);
        int nsubels = subels.size();
        for (int sub=0; sub<nsubels; sub++) {
            int64_t index = subels[sub]->Index();
            if(extpartition.size() < index+1) extpartition.resize(index+1);
            extpartition[index] = partition;
        }
    }
    fElementPartition = extpartition;

}