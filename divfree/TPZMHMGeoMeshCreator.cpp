
#include "TPZMHMGeoMeshCreator.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"
#include "TPZRefPatternDataBase.h"


//Creates the skeleton geometric elements
void TPZMHMGeoMeshCreator::CreateSkeleton(TPZAutoPointer<TPZGeoMesh> &gmesh){

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

void TPZMHMGeoMeshCreator::CreateSubGrids(TPZAutoPointer<TPZGeoMesh> &gmesh){

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

void TPZMHMGeoMeshCreator::RefineSubGrids(TPZAutoPointer<TPZGeoMesh> &gmesh){


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

void TPZMHMGeoMeshCreator::RefineSkeleton(TPZAutoPointer<TPZGeoMesh> &gmesh){

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


void TPZMHMGeoMeshCreator::AddBoundaryElements(TPZAutoPointer<TPZGeoMesh> &gmesh)
{
    std::set<int64_t> setbottom,setright,settop,setleft;
    int64_t nnodes = gmesh->NNodes();
    int dim = gmesh->Dimension();
    for (int64_t in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> xco(3);
        gmesh->NodeVec()[in].GetCoordinates(xco);
        if (fabs(xco[1]+1.) < 1.e-3)
        {
            setbottom.insert(in);
        }
        if (fabs(xco[0]-1.) < 1.e-3)
        {
            setright.insert(in);
        }
        if (fabs(xco[1]-1.) < 1.e-3)
        {
            settop.insert(in);
        }
        if (fabs(xco[0]+1.) < 1.e-3)
        {
            setleft.insert(in);
        }
    }
    int64_t nelem = gmesh->NElements();
    for (int64_t el=0; el<nelem; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            int nsidenodes = gel->NSideNodes(is);
            int nfoundbottom = 0;
            int nfoundright = 0;
            int nfoundtop = 0;
            int nfoundleft = 0;
            for (int in=0; in<nsidenodes; in++) {
                int64_t nodeindex = gel->SideNodeIndex(is, in);
                if (setbottom.find(nodeindex) != setbottom.end()) {
                    nfoundbottom++;
                }
                if (setright.find(nodeindex) != setright.end()) {
                    nfoundright++;
                }
                if (settop.find(nodeindex) != settop.end()) {
                    nfoundtop++;
                }
                if (setleft.find(nodeindex) != setleft.end()) {
                    nfoundleft++;
                }
            }
            if (nfoundbottom == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,fBC1);
            }
            if (nfoundright == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,fBC2);
            }
            if (nfoundtop == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,fBC3);
            }
            if (nfoundleft == nsidenodes) {
                TPZGeoElBC gelbc(gel,is,fBC4);
            }
            else
            {
                TPZGeoElSide gelside(gel,is);
                TPZGeoElSide neighbour = gelside.Neighbour();
                if (neighbour == gelside) {
                    int EPoly = 100;
                    TPZGeoElBC(gelside,EPoly);
                }
            }
        }
    }
}

void TPZMHMGeoMeshCreator::CreateTriangleElements(TPZAutoPointer<TPZGeoMesh> gmesh, std::map<int,int> &matmap, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenter)
{
    int64_t nel = gmesh->NElements();
    TPZManVector<int64_t> extpartition(elpartition);
    extpartition.Expand(nel*3);
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        int gelmat = gel->MaterialId();
        if(gel->Type() != EOned) DebugStop();
        auto partition = elpartition[el];
        if(partition >= 0)
        {
            if(matmap.find(gelmat) == matmap.end()) DebugStop();
            int64_t center = scalingcenter[partition];
            TPZManVector<int64_t,3> nodes(3);
            nodes[0] = gel->NodeIndex(0);
            nodes[1] = gel->NodeIndex(1);
            nodes[2] = center;
            int64_t index;
            gmesh->CreateGeoElement(ETriangle, nodes, matmap[gelmat], index);
            if(extpartition.size() < index+1) extpartition.resize(index+1);
            extpartition[index] = partition;
            // set the oned element partition as -1. It will serve as skeleton element
            extpartition[el] = -1;
        }
    }
    gmesh->BuildConnectivity();
    elpartition = extpartition;
}