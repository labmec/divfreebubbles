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
#include "TPZCompElKernelHDiv.h"
#include "TPZCompElKernelHDiv3D.h"
#include "TPZCompElKernelHDivBC.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZLagrangeMultiplierCS.h"
#include <TPZNullMaterialCS.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZMixedDarcyFlowHybrid.h"
#include "TPZCompElDisc.h"
#include "TPZCompElH1.h"
#include <pzbuildmultiphysicsmesh.h>
#include "Projection/TPZHCurlProjection.h"
#include "pzelchdiv.h"
#include "pzelchdivbound2.h"


#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"

auto forcefunction = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u){
    const auto &x=loc[0];
    const auto &y=loc[1];

    // //Nabla u = 1
    u[0] = 0.;
};


template<class TVar>
TPZApproxSpaceKernelHdiv<TVar>::TPZApproxSpaceKernelHdiv(TPZGeoMesh *gmesh, MSpaceType spacetype, MShapeType shapetype) :
            fSpaceType(spacetype), fShapeType(shapetype), fGeoMesh(gmesh) {
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
        hybridizer.SetEdgeRemove(fConfig.fEdgeRemove);
        hybridizer.CreateWrapElements(fGeoMesh,fConfig.fBCHybridMatId,true);
    } else {
        hybridizer.SetPeriferalMaterialIds(fConfig.fWrap,fConfig.fLagrange,fConfig.fInterface,fConfig.fPoint,fConfig.fDomain);
        hybridizer.CreateWrapElements(fGeoMesh,fConfig.fBCHybridMatId,false);
    }

    if (fDimension == 3) CreateOrientedBoundaryElements();
    
}


template<class TVar>
void TPZApproxSpaceKernelHdiv<TVar>::Solve(TPZLinearAnalysis &an, TPZCompMesh * cmesh, bool direct, bool filterEquations)
{
    if (direct)
    {
        bool domainHybridization = false;
        if (fSpaceType != ENone) domainHybridization = true;
        util->SolveProblemDirect(an,cmesh,filterEquations,domainHybridization);
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
    //Just the BC's not hybridized
    
    allMat.insert(fConfig.fDomain);
    allMat.insert(fConfig.fPoint);
    allMat.insert(fConfig.fWrap);
    allMat.insert(fConfig.fEdgeRemove);



    for (std::set<int>::iterator it=allMat.begin(); it!=allMat.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        cmesh->InsertMaterialObject(mat);
        mat->SetDimension(fDimension);
        if (*it == fConfig.fEdgeRemove) mat->SetDimension(1);
        mat->SetBigNumber(1.e10);
    } 

    for (std::set<int>::iterator it=fConfig.fBCMatId.begin(); it!=fConfig.fBCMatId.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it,fDimension-1);
        cmesh->InsertMaterialObject(mat);
        mat->SetBigNumber(1.e10);
    } 

    set_symmetric_difference(fConfig.fBCMatId.begin(), fConfig.fBCMatId.end(), fConfig.fBCHybridMatId.begin(), fConfig.fBCHybridMatId.end(), inserter(allMat, allMat.begin()));
    
    // switch (fShapeType)
    // {
    // case EHDivKernel:
    //     cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivKernel);
    //     break;
    // case EHCurlNoGrads:
    //     cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHCurlNoGrads);
    //     break;
    // case EHDivConstant:
    //     cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
    //     break;

    // default:
    //     break;
    // }
    
    // cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(fDimension);
    // cmesh->AutoBuild();

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
            if (fDimension == 3) continue;
            new TPZCompElH1<TPZShapePoint>(*cmesh,gel);
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            // nullmat->SetDimension(0);
        } else if (type == EOned){
            if (matid==fConfig.fEdgeRemove && fConfig.fEdgeRemove > 0){
                TPZCompElDisc* disc = new TPZCompElDisc(*cmesh,gel);
                TPZMaterial *mat = cmesh->FindMaterial(matid);
                TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                nullmat->SetDimension(1);
            }
            if (fSpaceType != ENone) continue;
            if (allMat.find(matid) == allMat.end()) continue;
            
            if (fShapeType == EHDivKernel || fShapeType == EHCurlNoGrads){
                new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,gel);
            } else if (fShapeType == EHDivConstant) {
                new TPZCompElHDivBound2<TPZShapeLinear>(*cmesh,gel,HDivFamily::EHDivConstant);
            } else {
                DebugStop();
            }
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(1);
            if (matid == fConfig.fWrap){
                gel->ResetReference();
            }

        } else if (type == EQuadrilateral){
            if (fDimension == 2){
                if (fShapeType == EHDivKernel || fShapeType == EHCurlNoGrads){
                    new TPZCompElKernelHDiv<TPZShapeQuad>(*cmesh,gel);
                } else if (fShapeType == EHDivConstant) {
                    new TPZCompElHDiv<TPZShapeQuad>(*cmesh,gel,HDivFamily::EHDivConstant);
                } else {
                    DebugStop();
                }
                TPZMaterial *mat = cmesh->FindMaterial(matid);
                TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                nullmat->SetDimension(2);
            } else if (fDimension == 3){
                if (fSpaceType != ENone) continue;
                if (allMat.find(matid) == allMat.end()) continue;
                if (fShapeType == EHDivKernel || fShapeType == EHCurlNoGrads){
                    new TPZCompElKernelHDiv3D<TPZShapeQuad>(*cmesh,gel,HDivFamily::EHDivKernel);
                } else if (fShapeType == EHDivConstant) {
                    new TPZCompElHDivBound2<TPZShapeQuad>(*cmesh,gel,HDivFamily::EHDivConstant);
                } else {
                    DebugStop();
                }
                TPZMaterial *mat = cmesh->FindMaterial(matid);
                TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);                
                // TPZHCurlProjection<> *nullmat = dynamic_cast<TPZHCurlProjection<> *>(mat);                
                nullmat->SetDimension(2);
                if (matid == fConfig.fWrap){
                    gel->ResetReference();
                }
            }
        } else if(type == ETriangle) {
            if (fDimension == 2){
                if (fShapeType == EHDivKernel){
                    new TPZCompElKernelHDiv<TPZShapeTriang>(*cmesh,gel);
                } else if (fShapeType == EHDivConstant){
                    new TPZCompElHDiv<TPZShapeTriang>(*cmesh,gel,HDivFamily::EHDivConstant);
                } else {
                    DebugStop();
                }
                TPZMaterial *mat = cmesh->FindMaterial(matid);
                TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                nullmat->SetDimension(2);
            } else if (fDimension == 3){
                if (fSpaceType != ENone) continue;
                if (allMat.find(matid) == allMat.end()) continue;
                if (fShapeType == EHDivKernel || fShapeType == EHCurlNoGrads){
                    new TPZCompElKernelHDiv3D<TPZShapeTriang>(*cmesh,gel,HDivFamily::EHDivKernel);
                } else if (fShapeType == EHDivConstant) {
                    // new TPZCompElHDivConstantBC<TPZShapeTriang>(*cmesh,gel);
                } else {
                    DebugStop();
                }
                TPZMaterial *mat = cmesh->FindMaterial(matid);
                TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);                
                nullmat->SetDimension(2);
                if (matid == fConfig.fWrap){
                    gel->ResetReference();
                }
            }
        } else if(type == ETetraedro) {
            if (fShapeType == EHDivKernel || fShapeType == EHCurlNoGrads){
                new TPZCompElKernelHDiv3D<TPZShapeTetra>(*cmesh,gel,HDivFamily::EHDivKernel);
            } else if (fShapeType == EHDivConstant){
                // new TPZCompElHDivConstant<TPZShapeTetra>(*cmesh,gel);
            } else {
                DebugStop();
            }
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(3);
        } else if(type == ECube) {
            if (fShapeType == EHDivKernel || fShapeType == EHCurlNoGrads){
                new TPZCompElKernelHDiv3D<TPZShapeCube>(*cmesh,gel,HDivFamily::EHDivKernel);
            } else if (fShapeType == EHDivConstant){
                // new TPZCompElHDivConstant<TPZShapeCube>(*cmesh,gel);
            } else {
                DebugStop();
            }
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(3);
        } else if(type == EPrisma) {
            if (fShapeType == EHDivKernel || fShapeType == EHCurlNoGrads){
                new TPZCompElKernelHDiv3D<TPZShapePrism>(*cmesh,gel,HDivFamily::EHDivKernel);
            } else if (fShapeType == EHDivConstant){
                // new TPZCompElHDivConstant<TPZShapePrism>(*cmesh,gel);
            } else {
                DebugStop();
            }
            TPZMaterial *mat = cmesh->FindMaterial(matid);
            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            nullmat->SetDimension(3);
        }
        
        if (fSpaceType == ENone) continue;//No hybridization
        if (gel->Dimension() != fDimension) continue;
        //Creates point
        if (fDimension == 2){
            TPZGeoElSide gelside(gel,0);
            TPZGeoElSide neighbour = gelside.Neighbour();
            if (fShapeType == EHDivConstant) continue;
            if (neighbour.Element()->MaterialId() == fConfig.fPoint){
                new TPZCompElH1<TPZShapePoint>(*cmesh,neighbour.Element());
                TPZMaterial *mat = cmesh->FindMaterial(fConfig.fPoint);
                TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
            }
        }      

        for (int side=gel->NCornerNodes(); side<gel->NSides()-1; side++){
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            if (gelside.Dimension() != gel->Dimension()-1) continue;

            if (neighbour.Element()->MaterialId() == fConfig.fWrap){
                if (fDimension == 2){
                    if (fShapeType == EHDivKernel || fShapeType == EHCurlNoGrads){
                        new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,neighbour.Element());
                    } else if (fShapeType == EHDivConstant){
                        // new TPZCompElHDivConstantBC<TPZShapeLinear>(*cmesh,neighbour.Element());
                    } else {
                        DebugStop();
                    }
                    TPZMaterial *mat = cmesh->FindMaterial(fConfig.fWrap);
                    TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                    nullmat->SetDimension(1);
                } else if (fDimension == 3){
                    switch (type)
                    {
                    case ETetraedro:
                        new TPZCompElKernelHDiv3D<TPZShapeTriang>(*cmesh,neighbour.Element(),HDivFamily::EHDivKernel);
                        break;
                    case ECube:
                        new TPZCompElKernelHDiv3D<TPZShapeQuad>(*cmesh,neighbour.Element(),HDivFamily::EHDivKernel);
                        break;
                    
                    default:
                        DebugStop();
                        break;
                    }
                    
                    TPZMaterial *mat = cmesh->FindMaterial(fConfig.fWrap);
                    TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                    nullmat->SetDimension(2);
                }
            } else {
                if (allMat.find(neighbour.Element()->MaterialId()) != allMat.end()){
                    if (fDimension == 2){
                        if (fShapeType == EHDivKernel || fShapeType == EHCurlNoGrads){
                            new TPZCompElKernelHDivBC<TPZShapeLinear>(*cmesh,neighbour.Element());
                        } else if (fShapeType == EHDivConstant) {
                            // new TPZCompElHDivConstantBC<TPZShapeLinear>(*cmesh,neighbour.Element());
                        } else {
                            DebugStop();
                        }
                            TPZMaterial *mat = cmesh->FindMaterial(neighbour.Element()->MaterialId());
                            TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                            nullmat->SetDimension(1);
                        
                    } else if (fDimension == 3) {
                        switch (type)
                        {
                        case ETetraedro:
                            new TPZCompElKernelHDiv3D<TPZShapeTriang>(*cmesh,neighbour.Element(),HDivFamily::EHDivKernel);
                            break;
                        case ECube:
                            new TPZCompElKernelHDiv3D<TPZShapeQuad>(*cmesh,neighbour.Element(),HDivFamily::EHDivKernel);
                            break;
    
                        default:
                            DebugStop();
                            break;
                        }
                        
                        TPZMaterial *mat = cmesh->FindMaterial(neighbour.Element()->MaterialId());
                        TPZNullMaterial<> *nullmat = dynamic_cast<TPZNullMaterial<> *>(mat);
                        nullmat->SetDimension(2);
                    }
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

    // TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cmesh->Element(0));
    // intel->SetSideOrder(20,2);
    // TPZConnect &con = intel->Connect(12);
    // con.SetOrder(0);
    // con.SetNShape(0);

    cmesh->InitializeBlock(); 
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    cmesh->LoadReferences();
    

    // if (fDimension == 3){
    //     OrientFaces(cmesh);
    // }

    if (fSpaceType == ESemiHybrid){
        hybridizer.SemiHybridizeFlux(cmesh,fConfig.fBCHybridMatId);
    }

    //Set the connect as equal to the removed edge
    // cmesh->LoadReferences();
    // int count = 0;
    // for (auto cel : cmesh->ElementVec())
    // {
    //     cel->LoadElementReference();
    //     auto gel = cel->Reference();
    //     if (gel->Dimension() != 1) continue;
    //     if (gel->MaterialId() != fConfig.fEdgeRemove) continue;
        
    //     cel->SetgOrder(0);
    //     TPZGeoElSide geoside(gel,2);
    //     TPZGeoElSide neighbour = geoside.Neighbour();

    //     while (neighbour.Element()->MaterialId() != fConfig.fDomain){
    //         neighbour = neighbour.Neighbour();
    //     }
        
    //     //Descobrir qual o connect de aresta vizinha ao elemento 1D, que será removida
    //     auto s = neighbour.Side() - 4;
    //     TPZGeoElSide neigh(neighbour.Element(),s);
    //     auto neigh2 = neigh.Neighbour();
    //     auto index = neighbour.Element()->Reference()->ConnectIndex(s);
        
    //     cel->SetConnectIndex(0,index);
       
    // }
    
    // hybridizer.EdgeRemove(cmesh);


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

    if (fShapeType == EHDivConstant){
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(fConfig.fDomain);
        mat->SetDimension(fDimension);
        cmesh->InsertMaterialObject(mat);
        mat -> SetBigNumber(1.e10);

        cmesh->SetAllCreateFunctionsDiscontinuous();

        cmesh->SetDefaultOrder(0);
        cmesh->SetDimModel(fDimension);
        cmesh->AutoBuild();

        int ncon = cmesh->NConnects();
        for(int i=0; i<ncon; i++)
        {
            TPZConnect &newnod = cmesh->ConnectVec()[i]; 
            newnod.SetLagrangeMultiplier(1);
        }

        int nel = cmesh->NElements();
        for(int i=0; i<nel; i++){
            TPZCompEl *cel = cmesh->ElementVec()[i];
            TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
            if(!celdisc) continue;
            celdisc->SetConstC(1.);
            celdisc->SetTrueUseQsiEta();
            // espera-se elemento de pressao apenas para o contorno
            auto aaa = celdisc->Reference()->Dimension();
            // if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
            // {
            //     DebugStop();
            // }
        }
    } else {

        cmesh->InitializeBlock();

        if(fDefaultPOrder <= 1 || fSpaceType == ESemiHybrid){
            cmesh->SetDefaultOrder(0);
            cmesh->SetDimModel(fDimension-1);
            cmesh->SetAllCreateFunctionsDiscontinuous();
            // cmesh->ApproxSpace().CreateDisconnectedElements(true);
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
   
    // eh preciso criar materiais para todos os valores referenciados no enum
    auto mat = new TPZMixedDarcyFlow(fConfig.fDomain,fDimension);
    // auto mat = new TPZMixedDarcyFlowHybrid(fConfig.fDomain,fDimension);
    mat->SetConstantPermeability(1.);
    mat->SetForcingFunction(forcefunction,1);

    // mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
    // mat->NStateVariables(3);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);

    //Dirichlet Boundary Conditions
    for (auto matId : BCDirichlet)
    {
        TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, matId, 0, val1, val2);
        BCond->SetForcingFunctionBC(exactSol);
        // BCond->SetForcingFunctionBC(exactSol,2);
        cmesh->InsertMaterialObject(BCond);
    }

    //Neumann Boundary Conditions
    for (auto matId : BCNeumann)
    {
        TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, matId, 1, val1, val2);
        BCond->SetForcingFunctionBC(exactSol);
        cmesh->InsertMaterialObject(BCond);
    }

    if (fConfig.fEdgeRemove > 0){
        TPZL2ProjectionCS<> *matEdge = new TPZL2ProjectionCS<>(fConfig.fEdgeRemove,0,1);
        matEdge->SetScaleFactor(1.e10);
        cmesh->InsertMaterialObject(matEdge);
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
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    if (fShapeType == EHDivConstant) {
        TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);
    }
    cmesh->LoadReferences();
    cmesh->CleanUpUnconnectedNodes(); 

    auto mat3 = new TPZLagrangeMultiplierCS<STATE>(fConfig.fInterface, fDimension-1);
    cmesh->InsertMaterialObject(mat3);

    auto matIdBCHyb = fConfig.fBCHybridMatId;
    matIdBCHyb.insert(fConfig.fLagrange);
    hybridizer.CreateMultiphysicsInterfaceElements(cmesh,fGeoMesh,meshvector,matIdBCHyb);

    return cmesh;
}


template<class TVar>
void TPZApproxSpaceKernelHdiv<TVar>::CreateOrientedBoundaryElements()
{   
    for(auto gel : fGeoMesh->ElementVec())
    {
        if (gel->Dimension() < 3) continue;
        int sideStart = 0;
        int sideEnd = 0;
        if (gel->Type() == ETetraedro){
            sideStart = 10;
            sideEnd = 14;
        } else if (gel->Type() == ECube){
            sideStart = 20;
            sideEnd = 26;
        } else if (gel->Type() == EPrisma){
            sideStart = 15;
            sideEnd = 20;
        }
        //For tetrahedra only, loop over the surface sides
        for (int side = sideStart; side < sideEnd; side++){
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            //Neighbour material id
            auto Nmatid = neighbour.Element()->MaterialId();

            /*  If the boundary has BC, delete the neighbour GeoElement and  
                create another one from TPZGeoElBC with the same material id
            */
            if (fConfig.fBCMatId.find(Nmatid) == fConfig.fBCMatId.end()) continue;
            fGeoMesh->DeleteElement(neighbour.Element(),neighbour.Element()->Index()); 
            TPZGeoElBC gelbcWrap(gelside, Nmatid);
           
        }
    }
}


template class TPZApproxSpaceKernelHdiv<STATE>;
// template class TPZApproxSpaceKernelHdiv<CSTATE>;


