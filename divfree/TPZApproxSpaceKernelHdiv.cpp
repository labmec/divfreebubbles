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
#include "TPZL2ProjectionCS.h"
#include "TPZLagrangeMultiplierCS.h"
#include <TPZNullMaterialCS.h>
#include "TPZMixedDarcyFlowHybrid.h"

#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"

template<class TVar>
TPZApproxSpaceKernelHdiv<TVar>::TPZApproxSpaceKernelHdiv(TPZGeoMesh *gmesh, MSpaceType spacetype) :
            fSpaceType(spacetype), fGeoMesh(gmesh) {
    fDimension = gmesh->Dimension();
}

/// copy constructor
template<class TVar>
TPZApproxSpaceKernelHdiv<TVar>::TPZApproxSpaceKernelHdiv(const TPZApproxSpaceKernelHdiv<TVar> &copy)
{
    
}

/// = operator
template<class TVar>
TPZApproxSpaceKernelHdiv<TVar> & TPZApproxSpaceKernelHdiv<TVar>::operator=(const TPZApproxSpaceKernelHdiv<TVar> &copy)
{
    return *this;
}

template<class TVar>
void TPZApproxSpaceKernelHdiv<TVar>::Initialize()
{
    if (fSpaceType != ENone)
    {
        hybridizer.SetPeriferalMaterialIds(fConfig.fWrap,fConfig.fLagrange,fConfig.fInterface,fConfig.fPoint,fConfig.fDomain);
        hybridizer.CreateWrapElements(fGeoMesh,fConfig.fBCHybridMatId,true);
    } else {
        hybridizer.SetPeriferalMaterialIds(fConfig.fWrap,fConfig.fLagrange,fConfig.fInterface,fConfig.fPoint,fConfig.fDomain);
        hybridizer.CreateWrapElements(fGeoMesh,fConfig.fBCHybridMatId,false);
    }
    
}


template<class TVar>
void TPZApproxSpaceKernelHdiv<TVar>::Solve(TPZLinearAnalysis &an, TPZMultiphysicsCompMesh * cmesh, bool direct)
{
    if (direct)
    {
        util->SolveProblemDirect(an,cmesh);
    } else {
        util->SolveProblemIterative(an,cmesh);
    }
}



//Creates the flux mesh
template<class TVar>
TPZCompMesh * TPZApproxSpaceKernelHdiv<TVar>::CreateFluxCMesh()
{   
    fGeoMesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(fDimension);

    //Inserts Null Materials
    std::set<int> allMat={};
    // if (fSpaceType == ENone) allMat = fConfig.fBCMatId;
    //Just the BC's not hybridized
    set_symmetric_difference(fConfig.fBCMatId.begin(), fConfig.fBCMatId.end(), fConfig.fBCHybridMatId.begin(), fConfig.fBCHybridMatId.end(), inserter(allMat, allMat.begin()));
    allMat.insert(fConfig.fDomain);
    allMat.insert(fConfig.fPoint);
    allMat.insert(fConfig.fWrap);

    for (std::set<int>::iterator it=allMat.begin(); it!=allMat.end(); ++it)
    {
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

        using namespace pzgeom;
        using namespace pzshape;

        if (type == EPoint){
            if (fSpaceType != ENone) continue;
            new TPZIntelGen<TPZShapePoint>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            // nullmat->SetDimension(0);
        } else if (type == EOned){
            if (fSpaceType != ENone) continue;
            if (allMat.find(matid) == allMat.end()) continue;
            
            new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(1);
            if (matid == fConfig.fWrap){
                gel->ResetReference();
            }

        } else if (type == EQuadrilateral){
            new TPZCompElKernelHDiv<TPZShapeQuad>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(2);
        } else if(type == ETriangle) {
            new TPZCompElKernelHDiv<TPZShapeTriang>(*cmesh,gel,index);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(2);
        }
        
        if (fSpaceType == ENone) continue;
        if (gel->Dimension() != fDimension) continue;
        //Creates point
        TPZGeoElSide gelside(gel,0);
        TPZGeoElSide neighbour = gelside.Neighbour();

        if (neighbour.Element()->MaterialId() == fConfig.fPoint){
            new TPZIntelGen<TPZShapePoint>(*cmesh,neighbour.Element(),index);
            TPZMaterial *mat = cmesh->FindMaterial(fConfig.fPoint);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
        }

        for (int side=gel->NCornerNodes(); side<gel->NSides()-1; side++){
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();

            if (neighbour.Element()->MaterialId() == fConfig.fWrap){
                new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,neighbour.Element(),index);
                TPZMaterial *mat = cmesh->FindMaterial(fConfig.fWrap);
                TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                nullmat->SetDimension(1);
                // neighbour.Element()->ResetReference();
            } else {
                if (allMat.find(neighbour.Element()->MaterialId()) != allMat.end()){
                    new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,neighbour.Element(),index);
                    TPZMaterial *mat = cmesh->FindMaterial(neighbour.Element()->MaterialId());
                    TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                    nullmat->SetDimension(1);
                    // neighbour.Element()->ResetReference();
                }
            }
        }
        for (int side = 0; side < gel->NSides(); side++)
        {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            neighbour.Element()->ResetReference();
        }
    }    

    cmesh->ExpandSolution();
    cmesh->InitializeBlock();
    cmesh->ComputeNodElCon();
    cmesh->LoadReferences();

    if (fSpaceType == ESemiHybrid){
        hybridizer.SemiHybridizeFlux(cmesh,fConfig.fBCHybridMatId);
    }

    return cmesh;
}

template<class TVar>
TPZCompMesh * TPZApproxSpaceKernelHdiv<TVar>::CreatePressureCMesh()
{
    std::set<int> matIdVec = fConfig.fBCHybridMatId;
    matIdVec.insert(fConfig.fLagrange);
    
    fGeoMesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);

    // Sets matid to BC geometric elements
    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        mat->SetDimension(1);
        cmesh->InsertMaterialObject(mat);
        mat->SetBigNumber(1.e10);
    }

    cmesh->InitializeBlock();

    if(fDefaultPOrder == 0 || fSpaceType == ESemiHybrid){
        cmesh->SetDefaultOrder(0);
        cmesh->SetDimModel(fDimension-1);
        cmesh->SetAllCreateFunctionsDiscontinuous();
    } else {
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
        cmesh->SetDefaultOrder(fDefaultPOrder-1);
        cmesh->SetDimModel(fDimension-1);
    }
    
    cmesh->AutoBuild(matIdVec);
    
    if (fSpaceType == ESemiHybrid){
        hybridizer.SemiHybridizePressure(cmesh,fDefaultPOrder,fConfig.fBCHybridMatId);
    }

    for(auto &newnod : cmesh->ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }

    

    return cmesh;
}

template<class TVar>
TPZMultiphysicsCompMesh * TPZApproxSpaceKernelHdiv<TVar>::CreateMultiphysicsCMesh(TPZVec<TPZCompMesh *> &meshvector, ForcingFunctionBCType<TVar> exactSol, std::set<int> &BCNeumann, std::set<int> &BCDirichlet)
{
    // gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(fDimension);
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();

    // eh preciso criar materiais para todos os valores referenciados no enum
    auto mat = new TPZMixedDarcyFlowHybrid(fConfig.fDomain,fDimension);
    mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
        
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);

    //Dirichlet Boundary Conditions
    for (auto matId : BCDirichlet)
    {
        TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, matId, 0, val1, val2);
        BCond->SetForcingFunctionBC(exactSol);
        cmesh->InsertMaterialObject(BCond);
    }

    //Neumann Boundary Conditions
    for (auto matId : BCNeumann)
    {
        TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, matId, 1, val1, val2);
        BCond->SetForcingFunctionBC(exactSol);
        cmesh->InsertMaterialObject(BCond);
    }


    auto *matL2 = new TPZL2ProjectionCS<>(fConfig.fPoint,0,1);
    cmesh->InsertMaterialObject(matL2);

    auto * nullmat2 = new TPZNullMaterialCS<>(fConfig.fWrap,fDimension-1,1);
    cmesh->InsertMaterialObject(nullmat2);

    auto * nullmat3 = new TPZNullMaterialCS<>(fConfig.fLagrange,fDimension-1,1);
    cmesh->InsertMaterialObject(nullmat3);

    TPZManVector<int> active(2,1);
    active[0]=1;
    active[1]=1;
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes(); 

    auto mat3 = new TPZLagrangeMultiplierCS<STATE>(fConfig.fInterface, fDimension-1);
    cmesh->InsertMaterialObject(mat3);

    auto matIdBCHyb = fConfig.fBCHybridMatId;
    matIdBCHyb.insert(fConfig.fLagrange);
    hybridizer.CreateMultiphysicsInterfaceElements(cmesh,fGeoMesh,meshvector,matIdBCHyb);

    return cmesh;
}

template class TPZApproxSpaceKernelHdiv<STATE>;
// template class TPZApproxSpaceKernelHdiv<CSTATE>;
