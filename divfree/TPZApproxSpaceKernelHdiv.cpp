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
#include "TPZCompElHDivSemiHybrid.h"
#include "TPZCompElHDivSemiHybridBound.h"

#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"
using namespace pzshape;

auto forcefunction = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u){
    const auto &x=loc[0];
    const auto &y=loc[1];

    // //Nabla u = 1
    u[0] = 0.;
};


template<class TVar>
TPZApproxSpaceKernelHdiv<TVar>::TPZApproxSpaceKernelHdiv(TPZGeoMesh *gmesh, MSpaceType spacetype, HDivFamily shapetype) :
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
        if (fSpaceType != ESemiHybrid && fShapeType != HDivFamily::EHDivKernel) {
            //In this case a different approach is used. Instead of using a geometric elements to hybridize the formu Dont create the wrap elements
        } else {
            hybridizer.CreateWrapElements(fGeoMesh,fConfig.fBCHybridMatId,true,fShapeType);
        }
    } else {
        hybridizer.SetPeriferalMaterialIds(fConfig.fWrap,fConfig.fLagrange,fConfig.fInterface,fConfig.fPoint,fConfig.fDomain);
        hybridizer.CreateWrapElements(fGeoMesh,fConfig.fBCHybridMatId,false,fShapeType);
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
    // allMat.insert(fConfig.fEdgeRemove);
    std::set<int> domainBC;// = fConfig.fBCMatId;
    domainBC.insert(fConfig.fDomain);


    for (std::set<int>::iterator it=allMat.begin(); it!=allMat.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        cmesh->InsertMaterialObject(mat);
        mat->SetDimension(fDimension);
        if (*it == fConfig.fEdgeRemove) mat->SetDimension(1);
        mat->SetBigNumber(fBigNumber);
    } 

    for (std::set<int>::iterator it=fConfig.fBCMatId.begin(); it!=fConfig.fBCMatId.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it,fDimension-1);
        cmesh->InsertMaterialObject(mat);
        mat->SetBigNumber(fBigNumber);
    } 

    set_symmetric_difference(fConfig.fBCMatId.begin(), fConfig.fBCMatId.end(), fConfig.fBCHybridMatId.begin(), fConfig.fBCHybridMatId.end(), inserter(allMat, allMat.begin()));
    
    //Creates computational elements
    if (fSpaceType == ENone && fConfig.fBCHybridMatId.size() == 0){//No hybridization (so easy...)
       
        cmesh->ApproxSpace().SetHDivFamily(fShapeType);
        cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(fDimension);
        cmesh->AutoBuild();

    } else { // Hybridization active
        if (fShapeType == HDivFamily::EHDivKernel){
            CreateFluxHybridezedHDivKernel(cmesh);
        } else if (fShapeType == HDivFamily::EHDivConstant) {
            CreateFluxHybridezedHDivConstant(cmesh);
        } else {
            std::cout << "You should hybridize your approximation space in other way." << std::endl;
            DebugStop();
        }
    }// end if hybridization

    return cmesh;
}

template<class TVar>
TPZCompMesh * TPZApproxSpaceKernelHdiv<TVar>::CreatePressureCMesh()
{
    std::set<int> matIdVec = fConfig.fBCHybridMatId;
    if (fShapeType != HDivFamily::EHDivConstant) matIdVec.insert(fConfig.fLagrange);
    
    fGeoMesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);

    // Sets matid to BC geometric elements
    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        mat->SetDimension(1);
        cmesh->InsertMaterialObject(mat);
        mat->SetBigNumber(fBigNumber);
    }

    if (fShapeType == HDivFamily::EHDivConstant || fShapeType == HDivFamily::EHDivStandard){
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(fConfig.fDomain);
        mat->SetDimension(fDimension);
        cmesh->InsertMaterialObject(mat);
        mat -> SetBigNumber(fBigNumber);

        if (fShapeType == HDivFamily::EHDivConstant) {
            cmesh->SetAllCreateFunctionsDiscontinuous();
            cmesh->SetDefaultOrder(0);
        }
        if (fShapeType == HDivFamily::EHDivStandard){
            cmesh->SetDefaultOrder(fDefaultPOrder);
            cmesh->SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        }
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

        if(fSpaceType == ESemiHybrid){
            cmesh->SetDefaultOrder(0);
            cmesh->SetDimModel(fDimension-1);
            cmesh->SetAllCreateFunctionsDiscontinuous();
            // cmesh->ApproxSpace().CreateDisconnectedElements(true);
        } else {
            cmesh->SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
            cmesh->SetDefaultOrder(fDefaultPOrder);
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
    mat -> SetBigNumber(fBigNumber);
    // mat->NStateVariables(3);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);

    //Dirichlet Boundary Conditions
    for (auto matId : BCDirichlet)
    {
        TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, matId, 0, val1, val2);
        // BCond->SetForcingFunctionBC(exactSol);
        BCond->SetForcingFunctionBC(exactSol,2);
        cmesh->InsertMaterialObject(BCond);
    }

    //Neumann Boundary Conditions
    for (auto matId : BCNeumann)
    {
        TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, matId, 1, val1, val2);
        BCond->SetForcingFunctionBC(exactSol,2);
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

    TPZManVector<int> active(meshvector.size(),1);
    // active[2]=0;
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    if (fShapeType == HDivFamily::EHDivConstant) {
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
        if (!gel || gel->Dimension() < fGeoMesh->Dimension()) continue;
        
        int nSides = gel->NSides();
        //For tetrahedra only, loop over the surface sides
        for (int side = 0; side < nSides; side++){
            if (gel->SideDimension(side) != gel->Dimension()-1) continue;

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


template<class TVar>
void TPZApproxSpaceKernelHdiv<TVar>::CreateFluxHybridezedHDivKernel(TPZCompMesh *cmesh)
{   
    int64_t nel = fGeoMesh->NElements();

    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);

        if(!gel) DebugStop();
        auto type = gel -> Type();
        auto matid = gel->MaterialId();

        // We only create points here if we are hybridizing HDivFamily::EHDivConstant
        if (gel->Dimension() == 0  && fSpaceType == ENone){
            new TPZCompElH1<TPZShapePoint>(*cmesh,gel,H1Family::EH1Standard);
        }

        if (matid != fConfig.fDomain) continue;

        using namespace pzgeom;
        using namespace pzshape;
        
        // First: create the volumetric element. Notice that there is no difference EHCurlNoGrads and EHDivKernel in 2D.
        if (fDimension == 2){
            switch (type){
            case ETriangle:
                CreateHDivKernelTriangleEl(gel,*cmesh,fShapeType);
                break;
            case EQuadrilateral:
                CreateHDivKernelQuadEl(gel,*cmesh,fShapeType);
                break;
            default:
                DebugStop();
                break;
            }
        } else if (fDimension == 3){
            switch (type){
            case ECube:
                CreateHDivKernelCubeEl(gel,*cmesh,fShapeType);
                break;
            case ETetraedro:
                CreateHDivKernelTetraEl(gel,*cmesh,fShapeType);
                break;
            default:
                DebugStop();
                break;
            }
        }
    
        // For dim = 2 we need to create a point element
        if (fDimension == 2 && fSpaceType != ENone){
            TPZGeoElSide gelside(gel,0);
            TPZGeoElSide neighbour = gelside.Neighbour();
            if (fShapeType == HDivFamily::EHDivConstant) continue;
            if (neighbour.Element()->MaterialId() == fConfig.fPoint){
                new TPZCompElH1<TPZShapePoint>(*cmesh,neighbour.Element(),H1Family::EH1Standard);
            } else {
                std::cout << "You need to provide a geometric point to each element for the hybridization in 2D\n" ;
                DebugStop();
            }
        }      

        // Go to the boundaries.
        for (int side=0; side<gel->NSides()-1; side++){
            
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            if (gelside.Dimension() != gel->Dimension()-1) continue;

            //Creates the computational elements. Both wrap, BC and interface
            if (fDimension == 2){
                if (neighbour.Element()->MaterialId() == fConfig.fDomain) continue;
                CreateHDivKernelBoundLinearEl(neighbour.Element(),*cmesh,fShapeType);                       
            } else if (fDimension == 3){
                switch (type){
                case ETetraedro:
                    CreateHDivKernelBoundTriangleEl(neighbour.Element(),*cmesh,fShapeType);
                    break;
                case ECube:
                    CreateHDivKernelBoundQuadEl(neighbour.Element(),*cmesh,fShapeType);
                    break;
                default:
                    DebugStop();
                    break;
                }
            }
        }
        if (fSpaceType != ENone){
            //Finally reset the reference of all neighbours
            for (int side = 0; side < gel->NSides(); side++)
            {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                neighbour.Element()->ResetReference();
            }
        }
    } 

    cmesh->InitializeBlock(); 
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    cmesh->LoadReferences();

    if (fSpaceType == ESemiHybrid){
        hybridizer.SemiHybridizeFlux(cmesh,fConfig.fBCHybridMatId);
    }
}

template<class TVar>
void TPZApproxSpaceKernelHdiv<TVar>::CreateFluxHybridezedHDivConstant(TPZCompMesh *cmesh)
{   
    int64_t nel = fGeoMesh->NElements();

    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);

        if(!gel) DebugStop();
        auto type = gel -> Type();
        auto matid = gel->MaterialId();
        if (matid != fConfig.fDomain) continue;

        using namespace pzgeom;
        using namespace pzshape;
        
        // First: create the volumetric element. Notice that there is no difference EHCurlNoGrads and EHDivKernel in 2D.
        if (fDimension == 2){
            switch (type){
            case ETriangle:
                if (fSpaceType == ESemiHybrid){
                    CreateHDivSemiHybridTriangleEl(gel,*cmesh,fShapeType);
                } else if (fSpaceType == EFullHybrid){
                    CreateHDivTriangleEl(gel,*cmesh,fShapeType);
                } else {
                    DebugStop();
                }
                break;
            case EQuadrilateral:
                if (fSpaceType == ESemiHybrid){
                    CreateHDivSemiHybridQuadEl(gel,*cmesh,fShapeType);
                } else if (fSpaceType == EFullHybrid){
                    CreateHDivQuadEl(gel,*cmesh,fShapeType);
                } else {
                    DebugStop();
                }
                break;
            default:
                DebugStop();
                break;
            }
        } else if (fDimension == 3){
            switch (type){
            case ECube:
                if (fSpaceType == ESemiHybrid){
                    CreateHDivSemiHybridCubeEl(gel,*cmesh,fShapeType);
                } else if (fSpaceType == EFullHybrid){
                    CreateHDivCubeEl(gel,*cmesh,fShapeType);
                } else {
                    DebugStop();
                }
                break;
            case ETetraedro:
                if (fSpaceType == ESemiHybrid){
                    CreateHDivSemiHybridTetraEl(gel,*cmesh,fShapeType);
                } else if (fSpaceType == EFullHybrid){
                    CreateHDivTetraEl(gel,*cmesh,fShapeType);
                } else {
                    DebugStop();
                }
                break;
            default:
                DebugStop();
                break;
            }
        }
    
        // Go to the boundaries.
        for (int side=0; side<gel->NSides()-1; side++){
            
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            if (gelside.Dimension() != gel->Dimension()-1) continue;

            //Creates the computational elements. Both wrap, BC and interface
            if (fDimension == 2){
                if (fSpaceType == ESemiHybrid){
                    if (fConfig.fBCMatId.find(neighbour.Element()->MaterialId()) != fConfig.fBCMatId.end()){
                        CreateHDivSemiHybridLinearEl(neighbour.Element(),*cmesh,fShapeType);
                    }
                }else {
                    CreateHDivBoundLinearEl(neighbour.Element(),*cmesh,fShapeType);   
                }
            } else if (fDimension == 3){
                switch (type){
                case ETetraedro:
                    CreateHDivBoundTriangleEl(neighbour.Element(),*cmesh,fShapeType);
                    break;
                case ECube:
                    CreateHDivBoundQuadEl(neighbour.Element(),*cmesh,fShapeType);
                    break;
                default:
                    DebugStop();
                    break;
                }
            }
        }
        //Finally reset the reference of all neighbours
        if (fSpaceType == EFullHybrid){
            for (int side = 0; side < gel->NSides(); side++)
            {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.Neighbour();
                neighbour.Element()->ResetReference();
            }
        }
    } 

    if (fSpaceType == ESemiHybrid){
        DuplicateInternalConnects(cmesh);
    }

    cmesh->InitializeBlock(); 
    cmesh->ExpandSolution();
    cmesh->ComputeNodElCon();
    cmesh->LoadReferences();

    // When hybridization in active, the side orient needs to be set as one so there is no need to vector compatibility
    // between elements. What is needed is that the flux orientation be outward the element.
    if (fSpaceType == EFullHybrid){
        for (auto cel:cmesh->ElementVec())
        {   
            if (cel->Reference()->Dimension() != fDimension) continue;
            
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            auto nsides = cel->Reference()->NSides();
            auto nfacets = cel->Reference()->NSides(fDimension-1);

            int firstside = nsides-nfacets-1;
            for (int i = firstside; i<nsides-1; i++){
                intel->SetSideOrient(i,1);
            }
        }
    }
    // if (fSpaceType == ESemiHybrid){
    //     hybridizer.SemiHybridizeFlux(cmesh,fConfig.fBCHybridMatId);
    // }

}


template<class TVar>
TPZCompMesh * TPZApproxSpaceKernelHdiv<TVar>::CreatePressureCMeshHybridizedHDivConstant()
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
        mat->SetBigNumber(fBigNumber);
    }

    if(fSpaceType == ESemiHybrid){
        cmesh->SetDefaultOrder(0);
        cmesh->SetDimModel(fDimension-1);
        cmesh->SetAllCreateFunctionsDiscontinuous();
        // cmesh->ApproxSpace().CreateDisconnectedElements(true);
    } else {
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
        cmesh->SetDefaultOrder(fDefaultPOrder);
        cmesh->SetDimModel(fDimension-1);
    }
    
    cmesh->AutoBuild(matIdVec);
    
    // if (fSpaceType == ESemiHybrid){
    //     // DebugStop();
    //     hybridizer.SemiHybridizePressure(cmesh,fDefaultPOrder,fConfig.fBCHybridMatId);
    // }

    
    for(auto &newnod : cmesh->ConnectVec())
    {
        newnod.SetLagrangeMultiplier(1);
    }

    return cmesh;
}

template<class TVar>
void TPZApproxSpaceKernelHdiv<TVar>::DuplicateInternalConnects(TPZCompMesh *cmesh)
{
    std::map<int64_t,int64_t> conn2duplConn; 
    util->PrintCMeshConnects(cmesh);
    int dim = cmesh->Dimension();

    //Loop over the computational elements
    for (auto cel:cmesh->ElementVec())
    {
        auto gel = cel->Reference();
        if (!gel) DebugStop();
        
        auto nFacets = gel->NSides(dim-1);

        //Loop over the element facets - which are the connects the be duplicated (edges in 2D and faces in 3D)
        for (int iFacet = 0; iFacet < nFacets; iFacet++)
        {
            // Algoritmo: Percorre cada facet do elemento, verifica se o connect já está no mapa, se já está,
            // apenas fala que o connect seguinte (duplicado) é igual ao valor retornado pelo std::map
            // Caso contrário, aloca um novo connect e insere o seu index no std::map.

            auto conn = cel->ConnectIndex(2*iFacet);           

            if (conn2duplConn.find(conn) == conn2duplConn.end()){
                //not found, so allocate a new connect
                auto pOrder = cmesh->GetDefaultOrder();
                int nshape = 0;
                int nstate = 1;
                int64_t newConnect = cmesh->AllocateNewConnect(nshape,nstate,pOrder);//Atualizar fConnectIndexes para ficar na ordem desejada.        
                conn2duplConn[conn] = newConnect;
                cel->SetConnectIndex(2*iFacet+1,newConnect);
            } else {
                //found, so just set the proper connect index
                cel->SetConnectIndex(2*iFacet+1,conn2duplConn[conn]);
            }
        }
        //Updates the number of shape functions and also the integration rule
        TPZCompElHDivSemiHybrid<TPZShapeQuad> *celHybrid = dynamic_cast<TPZCompElHDivSemiHybrid<TPZShapeQuad> *> (cel); 
        if (celHybrid){     
            for (int icon = 0; icon < celHybrid->NConnects(); icon++)
            {
                TPZConnect &c = celHybrid->Connect(icon);
                int nShapeF = celHybrid->NConnectShapeF(icon,c.Order());
                c.SetNShape(nShapeF);
                int64_t seqnum = c.SequenceNumber();
                int nvar = 1;
                TPZMaterial * mat = celHybrid->Material();
                if (mat) nvar = mat->NStateVariables();
                celHybrid->Mesh()->Block().Set(seqnum, nvar * nShapeF);
            }
        }
        TPZCompElHDivSemiHybridBound<TPZShapeLinear> *celHybBound = dynamic_cast<TPZCompElHDivSemiHybridBound<TPZShapeLinear> *> (cel);
        if (celHybBound){
            for (int icon = 0; icon < celHybBound->NConnects(); icon++)
            {
                TPZConnect &c = celHybBound->Connect(icon);
                int nShapeF = celHybBound->NConnectShapeF(icon,c.Order());
                c.SetNShape(nShapeF);
                int64_t seqnum = c.SequenceNumber();
                int nvar = 1;
                TPZMaterial * mat = celHybBound->Material();
                if (mat) nvar = mat->NStateVariables();
                celHybBound->Mesh()->Block().Set(seqnum, nvar * nShapeF);
            }
        
        }
        // this->Mesh()->ExpandSolution();
        // this->Mesh()->Block().Resequence();
    }
    util->PrintCMeshConnects(cmesh);
    cmesh->Block().Resequence();
    cmesh->ExpandSolution();

}

template class TPZApproxSpaceKernelHdiv<STATE>;
// template class TPZApproxSpaceKernelHdiv<CSTATE>;


