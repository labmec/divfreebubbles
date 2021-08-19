//
//  TPZApprxSpaceKernelHdiv.cpp
//  DivFreeBubbles project
//
//  Created by Jeferson Fernandes on 18/08/21.
//
//  Based on TPZCreateMultiphysicsSpace from ErrorEstimate project
//

#include "TPZApproxSpaceKernelHdiv.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoelbc.h"
#include <TPZNullMaterial.h>
#include "TPZCompElKernelHdiv.h"
#include "TPZCompElKernelHdivBC.h"
#include "TPZKernelHdivHybridizer.h"
#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"

TPZApproxSpaceKernelHdiv::TPZApproxSpaceKernelHdiv(TPZGeoMesh *gmesh, MSpaceType spacetype) :
            fSpaceType(spacetype), fGeoMesh(gmesh) {
    fDimension = gmesh->Dimension();
}

/// copy constructor
TPZApproxSpaceKernelHdiv::TPZApproxSpaceKernelHdiv(const TPZApproxSpaceKernelHdiv &copy)
{
    
}

/// = operator
TPZApproxSpaceKernelHdiv & TPZApproxSpaceKernelHdiv::operator=(const TPZApproxSpaceKernelHdiv &copy)
{
    return *this;
}


//Creates the flux mesh
TPZCompMesh * TPZApproxSpaceKernelHdiv::CreateFluxCMesh(std::set<int> &matIdVec)
{    
    fGeoMesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(fDimension);
    std::set<int> allmat;  

    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        allmat.insert(*it);
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        cmesh->InsertMaterialObject(mat);
        mat->SetDimension(fDimension);
        mat->SetBigNumber(1.e10);
    } 

    //Creates computational elements
    int64_t nel = fGeoMesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel) DebugStop();
        auto type = gel -> Type();
        int64_t index;
        auto matid = gel->MaterialId();
        if(allmat.find(matid) == allmat.end()) continue;

        using namespace pzgeom;
        using namespace pzshape;

        if (type == EPoint){
            if (fSpaceType != ENormal) continue;
            new TPZIntelGen<TPZShapePoint>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(0);
        } else if (type == EOned){
            if (fSpaceType != ENormal) continue;
            new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(1);
        } else if (type == EQuadrilateral){
            gel->ResetReference();   
            new TPZCompElKernelHDiv<TPZShapeQuad>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(2);
            
            if (fSpaceType == ENormal) continue;
            //Creates point
            TPZGeoElSide gelside(gel,0);
            TPZGeoElSide neighbour = gelside.Neighbour();

            if (neighbour.Element()->MaterialId() == fConfig.fPoint){
                new TPZIntelGen<TPZShapePoint>(*cmesh,neighbour.Element(),index);
                TPZMaterial *mat = cmesh->FindMaterial(fConfig.fPoint);
                TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            }

            for (int side=4; side<gel->NSides()-1; side++){
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();

                if (neighbour.Element()->MaterialId() == fConfig.fWrap){
                    new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,neighbour.Element(),index);
                    TPZMaterial *mat = cmesh->FindMaterial(fConfig.fWrap);
                    TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                    nullmat->SetDimension(1);
                    neighbour.Element()->ResetReference();
                } else {
                    if (allmat.find(neighbour.Element()->MaterialId()) != allmat.end()){
                        new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,neighbour.Element(),index);
                        TPZMaterial *mat = cmesh->FindMaterial(neighbour.Element()->MaterialId());
                        TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                        nullmat->SetDimension(1);
                        neighbour.Element()->ResetReference();
                    }
                }
            }
            for (int side = 0; side < gel->NSides(); side++)
            {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                neighbour.Element()->ResetReference();
            }
            // gel->ResetReference();    
        } else if(type == ETriangle) {
            new TPZCompElKernelHDiv<TPZShapeTriang>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(2);
        }
    }    

    cmesh->ExpandSolution();
    cmesh->InitializeBlock();
    cmesh->ComputeNodElCon();
    cmesh->LoadReferences();

    if (fSpaceType == EDomainSemiHybrid){
        TPZKernelHdivHybridizer hybridizer;
        hybridizer.SemiHybridizeFlux(cmesh,fBCMaterialIds);
    }

    return cmesh;
}

TPZCompMesh * TPZApproxSpaceKernelHdiv::CreatePressureCMesh(std::set<int> &matIdVec, std::set<int> &matBC){

    fGeoMesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    std::set<int> matIdNeumann;

    // Sets matid to BC geometric elements
    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        mat->SetDimension(1);
        cmesh->InsertMaterialObject(mat);
        mat->SetBigNumber(1.e10);
        matIdNeumann.insert(*it); 
    }

    cmesh->InitializeBlock();

    if(fDefaultPOrder == 0 || fSpaceType == EDomainSemiHybrid){
        cmesh->SetDefaultOrder(0);
        cmesh->SetDimModel(fDimension-1);
        cmesh->SetAllCreateFunctionsDiscontinuous();
    } else {
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
        cmesh->SetDefaultOrder(fDefaultPOrder);
        cmesh->SetDimModel(fDimension-1);
    }

    cmesh->AutoBuild(matIdNeumann);
    
    for(auto &newnod : cmesh->ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }

    if (fSpaceType == EDomainSemiHybrid){
        TPZKernelHdivHybridizer hybridizer;
        hybridizer.SemiHybridizePressure(cmesh,fDefaultPOrder,matBC);
    }

    return cmesh;
}

