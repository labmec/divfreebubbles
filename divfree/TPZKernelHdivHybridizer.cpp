//
// Created by Jeferson Fernandes on 11/08/21.
//

#include "TPZKernelHdivHybridizer.h"

#include <pzgmesh.h> 
#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"
#include "pznoderep.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"
#include "pzcompel.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

void TPZKernelHdivHybridizer::CreateWrapElements(TPZGeoMesh *gmesh, std::set<int> &matIdBC, bool domainHyb){

    int dim = gmesh->Dimension();

    //Disconnect elements
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) DebugStop();
        auto type = gel -> Type();
        int64_t index;
        auto matid = gel->MaterialId();

        using namespace pzgeom;
        using namespace pzshape;

        if (gel->Dimension() != dim) continue;
        
        int nsides = gel->NSides();
        // gel->ResetReference();
        for (int side = 0; side < nsides; side++) {

            if(gel->SideDimension(side) != dim-1) continue;
            TPZGeoElSide geoside(gel,side);
            TPZGeoElSide neighbour = geoside.Neighbour();

            if (domainHyb)
            {
                if (neighbour.Element()->MaterialId() == fEDomain)
                {
                    TPZGeoElBC gelbcWrap(geoside, fEWrap);
                                
                    TPZGeoElSide gelWrapSide(gelbcWrap.CreatedElement(),2);
                    // // CRIA OS ELEMENTOS GEOMÉTRICOS WRAP E INTERFACE
                    TPZGeoElBC gelbc(gelWrapSide, fEInterface); // AQUI CRIA O ELEMENTO DE INTERFACE GEOMÉTRICO
                    
                    bool flag = false;
                    TPZGeoElSide sidePresHyb;
                    auto ncorner = gel->NCornerNodes();

                    //Checks if the geometric element for the hybridized pressure already exists, otherwise creates it
                    for (int k = ncorner; k < nsides-1; k++)
                    {
                        TPZGeoElSide geosideNeig(neighbour.Element(),k);
                        if (geosideNeig.Element()->Neighbour(k).Element()->MaterialId()==fEWrap) {
                            flag = true;
                            sidePresHyb = geosideNeig.Element()->Neighbour(2).Element()->Neighbour(2).Element()->Neighbour(2);
                        }
                    }
                    if (flag==false){
                        TPZGeoElSide gelPresHSide(gelbc.CreatedElement(),2);
                        TPZGeoElBC gelPHyb(gelPresHSide, fEPressureHyb);
                    }
                } 
            }

            //Creates interface and wrap geometric elements for hybridized BC
            if (matIdBC.find(neighbour.Element()->MaterialId()) != matIdBC.end()){
                TPZGeoElBC gelbcWrap(geoside, fEWrap);
                TPZGeoElSide gelWrapSide(gelbcWrap.CreatedElement(),2);
                TPZGeoElBC gelbc(gelWrapSide, fEInterface);
                // gelbcWrap.CreatedElement()->ResetReference();
                // gelbc.CreatedElement()->ResetReference();
            }
        }

        //Creates a point for each hybrizides volumetric finite element
        if (domainHyb){
            TPZGeoElSide geosidePoint(gel,0);
            TPZGeoElBC gelbcPoint(geosidePoint, fEPont);
        // }
            for (int side = 0; side < nsides; side++) {
                TPZGeoElSide gelside(gel,side);
                gelside.Element()->ResetReference();
                TPZGeoElSide neighbour = gelside.Neighbour();
                neighbour.Element()->ResetReference();
            }   
            gel->ResetReference();
        }
    }
}




void TPZKernelHdivHybridizer::SemiHybridizeFlux(TPZCompMesh *cmesh, std::set<int> &matBCId)
{
    int dim = cmesh->Dimension();
    int nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int64_t index;

        if (cel->Dimension() != dim) continue;

        auto nsides = gel->NSides();
        for (int side = 0; side < nsides-1; side++)
        {       
            if (gel->SideDimension(side) != dim-1) continue;
            TPZGeoElSide gelside(gel,side);
            if (gelside.HasNeighbour(matBCId)) continue;
            int64_t cIndex = cel->ConnectIndex(side);
            TPZGeoElSide neighbour = gelside.Neighbour();

            while (neighbour != gelside)
            {
                TPZCompEl * neighcel = neighbour.Element()->Reference();
                if (neighcel)
                {
                    // TPZGeoElBC gelbc(gelWrapSide, fEInterface);
                    // cel->SetConnectIndex(neighbour.Side(),cIndex);
                    neighcel->SetConnectIndex(neighbour.Side(),cIndex);
                }
                neighbour=neighbour.Neighbour();
            }
        }
    }
    
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes();

}

void TPZKernelHdivHybridizer::CreateMultiphysicsInterfaceElements(TPZMultiphysicsCompMesh *cmesh, TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> &meshvector, std::set<int> &matIdBCHyb){

    cmesh->LoadReferences();
    
    for (auto gel : gmesh->ElementVec())
    {
        if (!gel || gel->MaterialId() != fEWrap) continue;
        auto nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        // here I generalized - an interface is created whenever a wrap element exists
        auto gelsidepr = gelside.HasNeighbour(matIdBCHyb);
        if (!gelsidepr)
        {
            DebugStop();
        }

        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide celneigh = gelsidepr.Reference();
        if (!celside || !celneigh) {
            DebugStop();
        }

        TPZGeoEl *gelIntface = gel->Neighbour(2).Element();
        if (gelIntface->MaterialId() != fEInterface) DebugStop();
        
        // Creates Multiphysics Interface element
        int64_t index;
        TPZMultiphysicsInterfaceElement *intf = new TPZMultiphysicsInterfaceElement(*cmesh,gelIntface,index,celneigh,celside);
    }

}

void TPZKernelHdivHybridizer::AssociateElements(TPZCompMesh *cmesh, TPZVec<int64_t> &elementgroup, std::set<int> &matIdBC)
{
    int64_t nel = cmesh->NElements();
    elementgroup.Resize(nel, -1);
    elementgroup.Fill(-1);
    // elementgroup2.Resize(nel, -1);
    // elementgroup2.Fill(-1);
    int64_t nconnects = cmesh->NConnects();
    TPZVec<int64_t> groupindex(nconnects, -1);
    TPZVec<int64_t> groupindex2(nconnects, -1);
    int dim = cmesh->Dimension();
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        cel->LoadElementReference();
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != dim) {
            continue;
        }
        elementgroup[cel->Index()] = cel->Index();
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);

        // std::cout << "Connect List = " << connectlist << std::endl;
        int k = -1;
        for (auto cindex : connectlist) {
            // k++;
            // auto gel = cel->Reference();
            // TPZGeoElSide geoside(gel,k);
            // if (geoside.Dimension() == dim-1){
            //     TPZGeoElSide neighbour = geoside.Neighbour();
            //     auto Nneighbour = neighbour.Element()->Neighbour(2);
            //     std::cout << "Nmaterial = " << neighbour.Element()->MaterialId() << " " << Nneighbour.Element()->MaterialId()<< " " << Nneighbour.Element()->Neighbour(2).Element()->MaterialId() << std::endl;
            //     if (matIdBC.find(Nneighbour.Element()->Neighbour(2).Element()->MaterialId()) != matIdBC.end()) {
            //         continue;
            //     }
            // }

            if (groupindex[cindex] != -1) {
                groupindex2[cindex] = cel->Index();
            } else {
                groupindex[cindex] = cel->Index();
            }
        }
    }
    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
//        std::cout << "Analysing element " << cel->Index();
        int64_t groupfound = -1;
        int k = -1;
        for (auto cindex : connectlist) {
            // k++;
            // auto gel = cel->Reference();
            // TPZGeoElSide geoside(gel,k);
            // if (geoside.Dimension() == dim-1){
            //     TPZGeoElSide neighbour = geoside.Neighbour();
            //     auto Nneighbour = neighbour.Element()->Neighbour(2);
            //     std::cout << "Nmaterial = " << neighbour.Element()->MaterialId() << " " << Nneighbour.Element()->MaterialId()<< " " << Nneighbour.Element()->Neighbour(2).Element()->MaterialId() << std::endl;
            //     if (matIdBC.find(Nneighbour.Element()->Neighbour(2).Element()->MaterialId()) != matIdBC.end()) {
            //         continue;
            //     }
            // }
            if (groupindex[cindex] != -1) {
                // assign the element to the group
                if(groupfound != -1 && groupfound != groupindex[cindex])
                {
                    //Do nothing
                    groupfound = groupindex2[cindex];
                }else{
                    elementgroup[cel->Index()] = groupindex[cindex];
                    groupfound = groupindex[cindex];
                }
            }
        }
//        std::cout << std::endl;
    }
}

void TPZKernelHdivHybridizer::GroupAndCondenseElements(TPZMultiphysicsCompMesh *cmesh, std::set<int> &matIdBC){

    int64_t nel = cmesh->NElements();
    TPZVec<int64_t> groupnumber(nel,-1);
    /// compute a groupnumber associated with each element
    AssociateElements(cmesh,groupnumber,matIdBC);

    std::map<int64_t, TPZElementGroup *> groupmap;
    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (int64_t el = 0; el<nel; el++) {
        int64_t groupnum = groupnumber[el];

        if(groupnum == -1) continue;
        auto iter = groupmap.find(groupnum);
        if (groupmap.find(groupnum) == groupmap.end()) {
            int64_t index;
            TPZElementGroup *elgr = new TPZElementGroup(*cmesh,index);
            groupmap[groupnum] = elgr;
            elgr->AddElement(cmesh->Element(el));
        }
        else
        {
            iter->second->AddElement(cmesh->Element(el));
        }
    }

    cmesh->ComputeNodElCon();
    nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (elgr) {
            TPZCondensedCompEl *cond = new TPZCondensedCompEl(elgr);
            cond->SetKeepMatrix(false);
        }
    }

    cmesh->InitializeBlock();
    cmesh->ComputeNodElCon();
}


void TPZKernelHdivHybridizer::SemiHybridizePressure(TPZCompMesh *cmesh, int pOrder, std::set<int> &matBCId)
{
 
    int64_t nel = cmesh->NElements();
    cmesh->Reference()->ResetReference();
    for(int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel ||cel->Reference()->Dimension() != cmesh->Dimension()) continue;
        auto gel = cel->Reference();
        int matid = gel->MaterialId();
        auto nconnects = cel->NConnects();
        
        if (matBCId.find(matid) == matBCId.end()) continue;

        auto *intel = dynamic_cast<TPZCompElDisc *>(cel);

        if(!intel) DebugStop();
        intel->PRefine(pOrder-1);
    }
    
    cmesh->ExpandSolution();

}