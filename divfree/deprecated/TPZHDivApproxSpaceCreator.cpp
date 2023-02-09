//
//  TPZApprxSpaceKernelHdiv.cpp
//  DivFreeBubbles project
//
//  Created by Jeferson Fernandes on 18/08/21.
//
//  Based on TPZCreateMultiphysicsSpace from ErrorEstimate project
//

#include "TPZHDivApproxSpaceCreator.h"

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
#include "TPZCompElHDivDuplConnects.h"
#include "TPZCompElHDivDuplConnectsBound.h"
#include "TPZCompElConstFluxHybrid.h"
#include "Elasticity/TPZElasticity2D.h"
#include "TPZMixedElasticityND.h"
#include "TPZAnalyticSolution.h"
#include "TPZCompelDiscScaled.h"

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
    const auto &z=loc[2];

    // //Nabla u = 1
    u[0] = 0.;
    // u[0] = 2*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
    // u[0] =  2*(x-1)*x*(y-1)*y + 2*(x-1)*x*(z-1)*z + 2*(y-1)*y*(z-1)*z;
    
    // double fx =-4144.653167389283*pow(10,
    //                                      -pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x) + 4771.70829943056*
    //    pow(10,-pow(-2*M_PI + 15.*x,2) -
    //        pow(-2*M_PI + 15.*y,2))*pow(-2*M_PI + 15.*x,3) +
    //    4771.70829943056*pow(10,
    //                         -pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x)*pow(-2*M_PI + 15.*y,2);
    
    // u[0] = -4144.653167389282*pow(2,2 - pow(-2*M_PI + 15.*x,2) -
    //                                       pow(-2*M_PI + 15.*y,2))*
    //    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x) + 4771.708299430558*
    //    pow(2,2 - pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    pow(-2*M_PI + 15.*x,3) +
    //    4771.708299430558*pow(2,
    //                          2 - pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    pow(5,-pow(-2*M_PI + 15.*x,2) - pow(-2*M_PI + 15.*y,2))*
    //    (-2*M_PI + 15.*x)*pow(-2*M_PI + 15.*y,2);
};


template<class TVar>
TPZHDivApproxSpaceCreator<TVar>::TPZHDivApproxSpaceCreator(TPZGeoMesh *gmesh, MSpaceType spacetype, HDivFamily shapetype) :
            fSpaceType(spacetype), fShapeType(shapetype), fGeoMesh(gmesh) {
    fDimension = gmesh->Dimension();
}

/// copy constructor
template<class TVar>
TPZHDivApproxSpaceCreator<TVar>::TPZHDivApproxSpaceCreator(const TPZHDivApproxSpaceCreator<TVar> &copy)
{
    
}

/// = operator
template<class TVar>
TPZHDivApproxSpaceCreator<TVar> & TPZHDivApproxSpaceCreator<TVar>::operator=(const TPZHDivApproxSpaceCreator<TVar> &copy)
{
    return *this;
}

template<class TVar>
void TPZHDivApproxSpaceCreator<TVar>::Initialize()
{
    if (fSpaceType != ENone)
    {
        hybridizer.SetMaterialIds(fConfig.fWrap,fConfig.fLagrange,fConfig.fInterface,fConfig.fPoint,fConfig.fDomain);
        // if (fSpaceType != EDuplicatedConnects) {
            hybridizer.CreateWrapElements(fGeoMesh,fConfig.fBCHybridMatId,true,fShapeType);
        // }
    } else {
        hybridizer.SetMaterialIds(fConfig.fWrap,fConfig.fLagrange,fConfig.fInterface,fConfig.fPoint,fConfig.fDomain);
        hybridizer.CreateWrapElements(fGeoMesh,fConfig.fBCHybridMatId,false,fShapeType);
    }

    if (fDimension == 3) CreateOrientedBoundaryElements();
    
}


template<class TVar>
void TPZHDivApproxSpaceCreator<TVar>::Solve(TPZLinearAnalysis &an, TPZCompMesh * cmesh, bool direct, bool filterEquations)
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
TPZCompMesh * TPZHDivApproxSpaceCreator<TVar>::CreateFluxCMesh()
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
    std::set<int> domainBC;// = fConfig.fBCMatId;
    domainBC.insert(fConfig.fDomain);


    for (std::set<int>::iterator it=allMat.begin(); it!=allMat.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        cmesh->InsertMaterialObject(mat);
        mat->SetDimension(fDimension);
        mat->SetBigNumber(fBigNumber);
        if (mixedElasticity) mat->SetNStateVariables(fDimension);
    } 

    for (std::set<int>::iterator it=fConfig.fBCMatId.begin(); it!=fConfig.fBCMatId.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it,fDimension-1);
        cmesh->InsertMaterialObject(mat);
        mat->SetBigNumber(fBigNumber);
        if (mixedElasticity) mat->SetNStateVariables(fDimension);
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
        } else if (fShapeType == HDivFamily::EHDivConstant || fShapeType == HDivFamily::EHDivStandard) {
            CreateFluxHybridezedHDivConstant(cmesh);
        } else {
            std::cout << "You should hybridize your approximation space in other way." << std::endl;
            DebugStop();
        }
    }// end if hybridization

    cmesh->InitializeBlock();

    return cmesh;
}

template<class TVar>
TPZCompMesh * TPZHDivApproxSpaceCreator<TVar>::CreatePressureCMesh()
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
        if (mixedElasticity) mat->SetNStateVariables(fDimension);
    }

    if (fShapeType == HDivFamily::EHDivConstant || fShapeType == HDivFamily::EHDivStandard){
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(fConfig.fDomain);
        mat->SetDimension(fDimension);
        cmesh->InsertMaterialObject(mat);
        mat -> SetBigNumber(fBigNumber);
        if (mixedElasticity) mat->SetNStateVariables(fDimension);

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
            newnod.SetLagrangeMultiplier(3);
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



    
    

    if(fSpaceType == EDuplicatedConnects){
        int64_t nconnects_i = cmesh->NConnects();
        // Sets matid to BC geometric elements
        std::set<int> matIdVec2 = fConfig.fBCHybridMatId;
        matIdVec2.insert(fConfig.fLagrange);
        for (std::set<int>::iterator it=matIdVec2.begin(); it!=matIdVec2.end(); ++it)
        {
            TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
            mat->SetDimension(fDimension-1);
            cmesh->InsertMaterialObject(mat);
            mat->SetBigNumber(fBigNumber);
            if(mixedElasticity) mat->SetNStateVariables(fDimension);
        }


        cmesh->SetDefaultOrder(0);
        cmesh->SetDimModel(fDimension-1);
        cmesh->SetAllCreateFunctionsDiscontinuous();
        cmesh->AutoBuild(matIdVec2);
        

        for(int64_t i = nconnects_i; i<cmesh->NConnects(); i++)
        {
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(10);
        }

    }
    
    cmesh->InitializeBlock();
    return cmesh;
}

template<class TVar>
TPZMultiphysicsCompMesh * TPZHDivApproxSpaceCreator<TVar>::CreateMultiphysicsCMesh(TPZVec<TPZCompMesh *> &meshvector, ForcingFunctionBCType<TVar> exactSol, std::set<int> &BCNeumann, std::set<int> &BCDirichlet)
{
    // gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(fDimension);
    
    // eh preciso criar materiais para todos os valores referenciados no enum
    auto mat = new TPZMixedDarcyFlow(fConfig.fDomain,fDimension);
    mat->SetConstantPermeability(1.);
    mat->SetForcingFunction(forcefunction,4);

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
        BCond->SetForcingFunctionBC(exactSol,4);
        cmesh->InsertMaterialObject(BCond);
    }

    //Neumann Boundary Conditions
    for (auto matId : BCNeumann)
    {
        TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, matId, 1, val1, val2);
        BCond->SetForcingFunctionBC(exactSol,4);
        cmesh->InsertMaterialObject(BCond);
    }

    int nstate = 1;
    if (mixedElasticity) nstate = fDimension;

    auto *matL2 = new TPZL2ProjectionCS<>(fConfig.fPoint,0,nstate);
    cmesh->InsertMaterialObject(matL2);

    auto * nullmat2 = new TPZNullMaterialCS<>(fConfig.fWrap,fDimension-1,nstate);
    cmesh->InsertMaterialObject(nullmat2);

    auto * nullmat3 = new TPZNullMaterialCS<>(fConfig.fLagrange,fDimension-1,nstate);
    cmesh->InsertMaterialObject(nullmat3);

    TPZManVector<int> active(meshvector.size(),1);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    // cmesh->AdjustBoundaryElements();
    cmesh->BuildMultiphysicsSpace(active, meshvector);
    // if (fShapeType == HDivFamily::EHDivConstant) {
    //     TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    //     TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);
    // }
    // cmesh->LoadReferences();
    // cmesh->CleanUpUnconnectedNodes(); 

    auto mat3 = new TPZLagrangeMultiplierCS<STATE>(fConfig.fInterface, fDimension-1);
    cmesh->InsertMaterialObject(mat3);

    auto matIdBCHyb = fConfig.fBCHybridMatId;
    matIdBCHyb.insert(fConfig.fLagrange);
    if (fSpaceType == EDuplicatedConnects){
        hybridizer.CreateInterfaceDuplConnects(cmesh,matIdBCHyb);
        // hybridizer.CreateMultiphysicsInterfaceElements(cmesh,fGeoMesh,meshvector,matIdBCHyb);
    } else {
        hybridizer.CreateMultiphysicsInterfaceElements(cmesh,fGeoMesh,matIdBCHyb);
    }

    return cmesh;
}

template<class TVar>
TPZMultiphysicsCompMesh * TPZHDivApproxSpaceCreator<TVar>::CreateMultiphysicsCMeshElasticity(TPZVec<TPZCompMesh *> &meshvector, TPZAnalyticSolution * gAnalytic, std::set<int> &BCNeumann, std::set<int> &BCDirichlet)
{
    // gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(fDimension);
    
    REAL E = 250.; //* @param E elasticity modulus
    REAL nu = 0.25; //* @param nu poisson coefficient

    
    TElasticity2DAnalytic * analytic2D = 0;
    TElasticity3DAnalytic * analytic3D = 0;
    if(fDimension == 2)
    {
        analytic2D = dynamic_cast<TElasticity2DAnalytic *>(gAnalytic);
        if(!analytic2D) DebugStop();
        TPZManVector<REAL,3> x(3,0.);
        // analytic2D->Elastic(x, E, nu);
    }
    else if(fDimension == 3)
    {
        analytic3D = dynamic_cast<TElasticity3DAnalytic *>(gAnalytic);
        if(!analytic3D) DebugStop();
        TPZManVector<REAL,3> x(3,0.);
        // analytic3D->Elastic(x, E, nu);
    }

    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = 0.; //* @param fx forcing function \f$ -x = fx \f$
    int plainStress = 0; //* @param plainstress = 1 \f$ indicates use of plainstress
    if(fDimension == 2)
    {
        if (analytic2D->fPlaneStress == 0) {
            plainStress = 0;
        } else {
            plainStress = 1;
        }
    }
    TPZMixedElasticityND * mat = new TPZMixedElasticityND(fConfig.fDomain, E, nu, fx, fy, plainStress, fDimension);
    mat->SetForcingFunction(gAnalytic->ForceFunc(),3);

    // // eh preciso criar materiais para todos os valores referenciados no enum
    // auto mat = new TPZMixedElasticityND(fConfig.fDomain,fDimension);
    // mat->SetConstantPermeability(1.);
    // mat->SetForcingFunction(forcefunction,4);

    // mat->SetPermeabilityFunction(1.);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(fBigNumber);
    // mat->NStateVariables(3);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(fDimension,fDimension,1.);
    TPZManVector<STATE> val2(fDimension,0.);

    //Dirichlet Boundary Conditions
    for (auto matId : BCDirichlet)
    {
        TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, matId, 0, val1, val2);
        BCond->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
        cmesh->InsertMaterialObject(BCond);
    }

    //Neumann Boundary Conditions
    for (auto matId : BCNeumann)
    {
        TPZBndCondT<STATE> * BCond = mat->CreateBC(mat, matId, 1, val1, val2);
        BCond->SetForcingFunctionBC(gAnalytic->ExactSolution(),4);
        cmesh->InsertMaterialObject(BCond);
    }

    auto *matL2 = new TPZL2ProjectionCS<>(fConfig.fPoint,0,fDimension);
    cmesh->InsertMaterialObject(matL2);

    auto * nullmat2 = new TPZNullMaterialCS<>(fConfig.fWrap,fDimension-1,fDimension);
    cmesh->InsertMaterialObject(nullmat2);

    auto * nullmat3 = new TPZNullMaterialCS<>(fConfig.fLagrange,fDimension-1,fDimension);
    cmesh->InsertMaterialObject(nullmat3);

    TPZManVector<int> active(meshvector.size(),1);
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

    auto mat3 = new TPZLagrangeMultiplierCS<STATE>(fConfig.fInterface, fDimension-1, fDimension);
    cmesh->InsertMaterialObject(mat3);

    auto matIdBCHyb = fConfig.fBCHybridMatId;
    matIdBCHyb.insert(fConfig.fLagrange);
    if (fSpaceType == EDuplicatedConnects){
        hybridizer.CreateInterfaceDuplConnects(cmesh,matIdBCHyb);
        // hybridizer.CreateMultiphysicsInterfaceElements(cmesh,fGeoMesh,meshvector,matIdBCHyb);
    } else {
        hybridizer.CreateMultiphysicsInterfaceElements(cmesh,fGeoMesh,matIdBCHyb);
    }

    return cmesh;
}


template<class TVar>
void TPZHDivApproxSpaceCreator<TVar>::CreateOrientedBoundaryElements()
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
void TPZHDivApproxSpaceCreator<TVar>::CreateFluxHybridezedHDivKernel(TPZCompMesh *cmesh)
{   
    int64_t nel = fGeoMesh->NElements();

    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);

        if(!gel) DebugStop();
        auto type = gel -> Type();
        auto matid = gel->MaterialId();
        if(gel->HasSubElement()) continue;

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
    cmesh->ComputeNodElCon();
    cmesh->LoadReferences();

    if (fSpaceType == ESemiHybrid || fSpaceType == EDuplicatedConnects){
        hybridizer.SemiHybridizeFlux(cmesh,fConfig.fBCHybridMatId);
    }
}

template<class TVar>
void TPZHDivApproxSpaceCreator<TVar>::CreateFluxHybridezedHDivConstant(TPZCompMesh *cmesh)
{   
    int64_t nel = fGeoMesh->NElements();

    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);

        if(!gel) DebugStop();
        auto type = gel -> Type();
        auto matid = gel->MaterialId();
        if (matid != fConfig.fDomain) continue;
        if(gel->HasSubElement()) continue;
        
        using namespace pzgeom;
        using namespace pzshape;
        
        // First: create the volumetric element. Notice that there is no difference EHCurlNoGrads and EHDivKernel in 2D.
        if (fDimension == 2){
            switch (type){
            case ETriangle:
                if (fSpaceType == EDuplicatedConnects){
                    CreateHDivDuplConnectsTriangleEl(gel,*cmesh,fShapeType);
                } else if (fSpaceType == EFullHybrid){
                    CreateHDivTriangleEl(gel,*cmesh,fShapeType);
                } else {
                    DebugStop();
                }
                break;
            case EQuadrilateral:
                if (fSpaceType == EDuplicatedConnects){
                    CreateHDivDuplConnectsQuadEl(gel,*cmesh,fShapeType);
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
                if (fSpaceType == EDuplicatedConnects){
                    CreateHDivDuplConnectsCubeEl(gel,*cmesh,fShapeType);
                } else if (fSpaceType == EFullHybrid){
                    CreateHDivCubeEl(gel,*cmesh,fShapeType);
                } else {
                    DebugStop();
                }
                break;
            case ETetraedro:
                if (fSpaceType == EDuplicatedConnects){
                    CreateHDivDuplConnectsTetraEl(gel,*cmesh,fShapeType);
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
                if (fSpaceType == EDuplicatedConnects){
                    // if (fConfig.fBCMatId.find(neighbour.Element()->MaterialId()) != fConfig.fBCMatId.end()){
                        CreateHDivDuplConnectsBoundLinearEl(neighbour.Element(),*cmesh,fShapeType);
                    // } else {
                        // new TPZCompElConstFluxHybrid(*cmesh,neighbour.Element());
                    // }
                }else {
                    CreateHDivBoundLinearEl(neighbour.Element(),*cmesh,fShapeType);   
                }
            } else if (fDimension == 3){
                switch (type){
                case ETetraedro:
                    if (fSpaceType == EDuplicatedConnects){
                        // if (fConfig.fBCMatId.find(neighbour.Element()->MaterialId()) != fConfig.fBCMatId.end()){
                            CreateHDivDuplConnectsBoundTriangEl(neighbour.Element(),*cmesh,fShapeType);
                        // }
                    } else {   
                        CreateHDivBoundTriangleEl(neighbour.Element(),*cmesh,fShapeType);
                    }
                    break;
                case ECube:
                    if (fSpaceType == EDuplicatedConnects){
                        // if (fConfig.fBCMatId.find(neighbour.Element()->MaterialId()) != fConfig.fBCMatId.end()){
                            CreateHDivDuplConnectsBoundQuadEl(neighbour.Element(),*cmesh,fShapeType);
                        // }
                    }else {
                        CreateHDivBoundQuadEl(neighbour.Element(),*cmesh,fShapeType);
                    }
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

    if (fSpaceType == EDuplicatedConnects){
        for (int64_t i = 0; i < cmesh->NElements(); i++)
        {
            auto celType = cmesh->Element(i)->Type();
            switch (celType)
            {
            case EOned:
                {
                    TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *> (cmesh->Element(i)); 
                    if (celb) celb->ActiveDuplConnects(fConnDuplicated);
                }
                break;

            case EQuadrilateral:
                {
                    TPZCompElHDivDuplConnects<TPZShapeQuad> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeQuad> *> (cmesh->Element(i)); 
                    if (celd) celd->ActiveDuplConnects(fConnDuplicated);
                    TPZCompElHDivDuplConnectsBound<TPZShapeQuad> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeQuad> *> (cmesh->Element(i)); 
                    if (celb) celb->ActiveDuplConnects(fConnDuplicated);
                }
                break;

            case ETriangle:
                {
                    TPZCompElHDivDuplConnects<TPZShapeTriang> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeTriang> *> (cmesh->Element(i)); 
                    if (celd) celd->ActiveDuplConnects(fConnDuplicated);
                    TPZCompElHDivDuplConnectsBound<TPZShapeTriang> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeTriang> *> (cmesh->Element(i)); 
                    if (celb) celb->ActiveDuplConnects(fConnDuplicated);
                }
                break;

            case ECube:
                {
                    TPZCompElHDivDuplConnects<TPZShapeCube> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeCube> *> (cmesh->Element(i)); 
                    if (celd) celd->ActiveDuplConnects(fConnDuplicated);
                }
                break;

            case ETetraedro:
                {
                    TPZCompElHDivDuplConnects<TPZShapeTetra> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeTetra> *> (cmesh->Element(i)); 
                    if (celd) celd->ActiveDuplConnects(fConnDuplicated);
                }
                break;

            default:
                DebugStop();
                break;
            }
        }

        PartitionDependMatrix(cmesh);
        // for (int64_t i = 0; i < cmesh->NElements(); i++)
        // {
        //     TPZCompElHDivDuplConnects<TPZShapeQuad> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeQuad> *> (cmesh->Element(i)); 
        //     TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *> (cmesh->Element(i)); 
        //     if (celd) celd->InactiveDuplConnects();
        //     if (celb) celb->InactiveDuplConnects();
        // }
        // RegroupDependMatrix(cmesh);
        
    }

    cmesh->InitializeBlock(); 

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
    if (fSpaceType == EDuplicatedConnects){
        hybridizer.SemiHybridizeDuplConnects(cmesh,fConfig.fBCHybridMatId);
    }

}


template<class TVar>
TPZCompMesh * TPZHDivApproxSpaceCreator<TVar>::CreatePressureCMeshHybridizedHDivConstant()
{
    std::set<int> matIdVec = fConfig.fBCHybridMatId;
    matIdVec.insert(fConfig.fLagrange);
    
    fGeoMesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);

    // Sets matid to BC geometric elements
    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        mat->SetDimension(fDimension-1);
        cmesh->InsertMaterialObject(mat);
        mat->SetBigNumber(fBigNumber);
        if(mixedElasticity) mat->SetNStateVariables(fDimension);
    }

    if(fSpaceType == ESemiHybrid || fSpaceType == EDuplicatedConnects){
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
        newnod.SetLagrangeMultiplier(10);
    }
    cmesh->InitializeBlock();

    return cmesh;
}

template<class TVar>
void TPZHDivApproxSpaceCreator<TVar>::DuplicateInternalConnects(TPZCompMesh *cmesh)
{
    // std::map<int64_t,int64_t> conn2duplConn; 

    //Loop over the computational elements
    for (int icel = 0; icel < cmesh->NElements(); icel++)
    {

        TPZCompEl* cel = cmesh->ElementVec()[icel];
        auto gel = cel->Reference();
        if (!gel) DebugStop();
        auto type = cel->Type();
        
        auto nFacets = gel->NSides(fDimension-1);

        // gel->ResetReference();
        // for (int side = 0; side < gel->NSides(); side++)
        // {
        //     TPZGeoElSide gelside(gel,side);
        //     TPZGeoElSide neighbour = gelside.Neighbour();
        //     neighbour.Element()->ResetReference();
        // }
        for (int i = 0; i<cel->NConnects(); i++){
            std::cout << "Connects indexes = "<< i << " " << cel->ConnectIndex(i) << "\n";
        }
        //Loop over the element facets - which are the connects the be duplicated (edges in 2D and faces in 3D)
        for (int iFacet = 0; iFacet < nFacets; iFacet++)
        {
            // Algorithm description: for each element facet, checks if the corresponding original connect is in the map fConnDuplicated.
            // If Yes, just sets the returning value from fConnDuplicated to the duplicated connect in the current element;
            // If No, allocate a new connect and inserts its index to fConnDuplicated using the original connect as key

            auto conn = cel->ConnectIndex(2*iFacet);           

            if (fConnDuplicated.find(conn) == fConnDuplicated.end()){
                //not found, so allocate a new connect
                auto pOrder = cmesh->GetDefaultOrder();
                int nshape = 0;//It is updated in the next loop
                int nstate = 1;//It can possibly change
                int64_t newConnect = cmesh->AllocateNewConnect(nshape,nstate,pOrder);
                fConnDuplicated[conn] = newConnect;
                cel->SetConnectIndex(2*iFacet+1,newConnect);
            } else {
                //found, so just set the proper index of the duplicated connect
                cel->SetConnectIndex(2*iFacet+1,fConnDuplicated[conn]);
            }
        }

        
        for (int i = 0; i<cel->NConnects(); i++){
            std::cout << "Connects indexes = "<< i << " " << cel->ConnectIndex(i) << "\n";
        }
        //Updates the number of shape functions and also the integration rule. 
        //We need different casts because the element can be volumetric or boundary
        TPZInterpolatedElement *celHybrid = dynamic_cast<TPZInterpolatedElement *> (cel); 
        if (celHybrid){
            int nConnects = celHybrid->NConnects();
            for (int icon = 0; icon < nConnects; icon++)
            {
                TPZConnect &c = celHybrid->Connect(icon);
                int nShapeF = celHybrid->NConnectShapeF(icon,c.Order());
                c.SetNShape(nShapeF);
                int64_t seqnum = c.SequenceNumber();
                int nvar = 1;
                TPZMaterial * mat = celHybrid->Material();
                if (mat) nvar = mat->NStateVariables();
                c.SetNState(nvar);
                celHybrid->Mesh()->Block().Set(seqnum, nvar * nShapeF);
            }
        }
        TPZInterpolatedElement *celHybBound = dynamic_cast<TPZInterpolatedElement *> (cel);
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
                c.SetNState(nvar);
                celHybBound->Mesh()->Block().Set(seqnum, nvar * nShapeF);
            }
        }
    }
    // util->PrintCMeshConnects(cmesh);
    cmesh->InitializeBlock();    
}

/**
 * @brief Generates the constant computational mesh
 * @param Gmesh: Geometric mesh
 * @param third_LM: Bool Third Lagrange multiplier
 * @return Constant computational mesh
 */
template<class TVar>
TPZCompMesh *TPZHDivApproxSpaceCreator<TVar>::CreateConstantCmesh(TPZGeoMesh *Gmesh, int lagLevel)
{
    TPZCompMesh *Cmesh= new TPZCompMesh (Gmesh);
    
    Cmesh->SetDimModel(Gmesh->Dimension());
    Cmesh->SetDefaultOrder(0);
    Cmesh->SetAllCreateFunctionsDiscontinuous();
    
    //Add material to the mesh
    int dimen = Gmesh->Dimension();
    int MaterialId = 1;
    
    TPZNullMaterial<> *mat =new TPZNullMaterial<>(MaterialId);
    mat->SetDimension(dimen);
    mat->SetNStateVariables(1);
    
    if (mixedElasticity){
        if(fDimension == 2){
            mat->SetNStateVariables(3);
        } else if(fDimension == 3) {
            mat->SetNStateVariables(6);
        } else {
            DebugStop();
        }
    }

    //Insert material to mesh
    Cmesh->InsertMaterialObject(mat);
    
    if (mixedElasticity){
        std::set<int> materialids;
        materialids.insert(fConfig.fDomain);
        //materialids.insert(3);
        {
            Gmesh->ResetReference();
            int64_t nel = Gmesh->NElements();
            for (int64_t el = 0; el < nel; el++) {
                TPZGeoEl *gel = Gmesh->Element(el);
                if (!gel)continue;
                if(gel->HasSubElement()) continue;
                int matid = gel->MaterialId();
                if (materialids.find(matid) == materialids.end()) {
                    continue;
                }
                TPZCompElDisc *disc = new TPZCompElDisc(*Cmesh, gel);
                disc->SetFalseUseQsiEta();
                gel->ResetReference();
            }
        }
        int ncon = Cmesh->NConnects();
        for (int i = 0; i < ncon; i++) {
            TPZConnect &newnod = Cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(lagLevel);
        }

        Cmesh->CleanUpUnconnectedNodes();
        Cmesh->ExpandSolution();

    }else {
        //Autobuild
        Cmesh->AutoBuild();
        
        int ncon = Cmesh->NConnects();
        for(int i=0; i<ncon; i++)
        {
            TPZConnect &newnod = Cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(lagLevel);
        }
        Cmesh->InitializeBlock();
    }
    return Cmesh;
}

/// change the order of the internal connect to the given order
template<class TVar>
void TPZHDivApproxSpaceCreator<TVar>::ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder) {
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;

        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != cmesh->Dimension()) {
            continue;
        }
        int nc = cel->NConnects();
        int64_t conIndex = cel->ConnectIndex(nc-1);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel) DebugStop();

        // intel->ForceSideOrder(gel->NSides() - 1, pOrder);
        TPZConnect &c = cmesh->ConnectVec()[conIndex];
        int64_t index;
        c.SetOrder(pOrder,index);
        int64_t seqnum = c.SequenceNumber();
        int nvar = 1;
        TPZMaterial * mat = cel->Material();
        if (mat) nvar = mat->NStateVariables();
        int nshape = intel->NConnectShapeF(nc-1,pOrder);
        c.SetNShape(nshape);
        c.SetNState(nvar);
        cmesh->Block().Set(seqnum, nvar * nshape);


        // int64_t seqnum = c.SequenceNumber();
        // int nvar = 1;
        // TPZMaterial * mat =this-> Material();
        // if(mat) nvar = mat->NStateVariables();
        // int nshape = NConnectShapeF(connectaux,order);
        // c.SetNShape(nshape);
        // c.SetNState(nvar);
        // this-> Mesh()->Block().Set(seqnum,nshape*nvar);
        // if(connectaux == NConnects()-1)
        // {
            cel->SetIntegrationRule(2*pOrder);
        // }


    }
    cmesh->InitializeBlock();
}

template<class TVar>
TPZCompMesh* TPZHDivApproxSpaceCreator<TVar>::CreateRotationCmesh(TPZGeoMesh *gmesh, int pOrder, REAL elementdim) {
    
    int dim = gmesh->Dimension();
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

    cmesh->SetAllCreateFunctionsDiscontinuous();

    //    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    //Criando material cujo nSTATE = 1:
    TPZNullMaterial<> *material = new TPZNullMaterial<>(fConfig.fDomain); //criando material que implementa a formulacao fraca do problema modelo
    material->SetDimension(dim);
    if(dim == 3)
    {
        material->SetNStateVariables(3);
    }

    cmesh->InsertMaterialObject(material); //Insere material na malha
    std::set<int> materialids;
    materialids.insert(fConfig.fDomain);
    //materialids.insert(3);
    {
        gmesh->ResetReference();
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el < nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel)continue;
            if(gel->HasSubElement()) continue;
            int matid = gel->MaterialId();
            if (materialids.find(matid) == materialids.end()) {
                continue;
            }
            new TPZCompElDiscScaled(*cmesh, gel);
            gel->ResetReference();
        }
    }

    //cmesh->LoadReferences();
    //    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    //    cmesh->AutoBuild();


    int ncon = cmesh->NConnects();
    for (int i = 0; i < ncon; i++) {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }

    int64_t nelem = cmesh->NElements();
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZCompElDiscScaled *disc = dynamic_cast<TPZCompElDiscScaled *> (cel);
        if (!disc) {
            continue;
        }
        disc->SetTotalOrderShape();
        disc->SetFalseUseQsiEta();
        disc->SetConstC(elementdim);
        disc->SetScale(1./elementdim);
    }
    cmesh->InitializeBlock();
    return cmesh;

}

template<class TVar>
void TPZHDivApproxSpaceCreator<TVar>::PartitionDependMatrix(TPZCompMesh *cmesh) {

    // Now that all edge connects have been duplicated, we need to expand the dependency matrix of a connect if it exists
    // For this purpose we can use the map fConnDuplicated which relates the original connect and its corresponding duplicated one
    // PS: This operation is performed outside the element duplication of connects because all connects have to be already duplicated 
    // to the dependency matrix expansion to be possible
    std::map<int64_t,bool> fDepMatTreated;
    for (int64_t iCon = 0; iCon < cmesh->NConnects(); iCon++)
    {
        TPZConnect &c = cmesh->ConnectVec()[iCon];
        if (fDepMatTreated[iCon]) continue;
        if(c.HasDependency()){
            // std::cout << "Need to partition the dependency matrix." << std::endl;
            auto *ptr = c.FirstDepend();
            
            while(ptr) {
                int64_t cIndex_old1 = iCon;
                int64_t cIndex_old2 = ptr->fDepConnectIndex;
                int64_t cIndex_new1 = fConnDuplicated[cIndex_old1];
                int64_t cIndex_new2 = fConnDuplicated[cIndex_old2];

                int rows = ptr->fDepMatrix.Rows();
                int cols = ptr->fDepMatrix.Cols();

                TPZFMatrix<REAL> DepMat00(1,1,0.), DepMat01(1,cols-1,0.), DepMat10(rows-1,1,0.), DepMat11(rows-1,cols-1,0.);
                DepMat00(0,0) = ptr->fDepMatrix(0,0);
                for (int i = 1; i < rows; i++){
                    DepMat10(i-1,0) = ptr->fDepMatrix(i,0);
                    for (int j = 1; j < cols; j++){
                        if (i == 1) DepMat01(0,j-1) = ptr->fDepMatrix(0,j);
                        DepMat11(i-1,j-1) = ptr->fDepMatrix(i,j);
                    }                  
                }
                // std::cout << "K00 " << DepMat00 << std::endl;
                // std::cout << "K01 " << DepMat01 << std::endl;
                // std::cout << "K10 " << DepMat10 << std::endl;
                // std::cout << "K11 " << DepMat11 << std::endl;
                
                c.RemoveDepend(cIndex_old1,cIndex_old2);
                
                fDepMatTreated[cIndex_old1] = true;
                fDepMatTreated[cIndex_new1] = true;
                
                //Get the duplicated connect and set the dependency matrix for all
                TPZConnect &c2 = cmesh->ConnectVec()[cIndex_new1];
                c2.RemoveDepend();

                //Dependency 00 - old1 + old2
                TPZConnect::TPZDepend *depend00 = c.AddDependency(cIndex_old1, cIndex_old2, DepMat00, 0, 0, 1, 1);

                //Dependency 01 - old1 + new2
                TPZConnect::TPZDepend *depend01 = c.AddDependency(cIndex_old1, cIndex_new2, DepMat01, 0, 0, 1, cols-1);

                //Dependency 10 - new1 + old2
                TPZConnect::TPZDepend *depend10 = c2.AddDependency(cIndex_new1, cIndex_old2, DepMat10, 0, 0, rows-1, 1);

                //Dependency 11 - new1 + new2
                TPZConnect::TPZDepend *depend11 = c2.AddDependency(cIndex_new1, cIndex_new2, DepMat11, 0, 0, rows-1, cols-1);

                ptr = ptr->fNext;
            }
        }
    }

}

template<class TVar>
void TPZHDivApproxSpaceCreator<TVar>::RegroupDependMatrix(TPZCompMesh *cmesh) {

    // Now that we have disabled the duplicated connect, we also need to restore the corresponding dependency matrix
    // We again use the map fConnDuplicated which relates the original connect and its corresponding duplicated one
    std::map<int64_t,bool> fDepMatTreated;

    for (int64_t iCon = 0; iCon < cmesh->NConnects(); iCon++){
        TPZConnect &c = cmesh->ConnectVec()[iCon];

        if (fDepMatTreated[iCon]) continue;

        if(c.HasDependency()){
            // std::cout << "Need to partition the dependency matrix." << std::endl;
            auto *ptr = c.FirstDepend();
            
            while(ptr) {

                int64_t cIndex_old1 = iCon;
                int64_t cIndex_old2;
                if (fConnDuplicated.find(ptr->fDepConnectIndex) != fConnDuplicated.end()){
                    cIndex_old2 = fConnDuplicated[ptr->fDepConnectIndex];
                } else {
                    int64_t findVal = ptr->fDepConnectIndex;    
                    auto it = find_if(fConnDuplicated.begin(), fConnDuplicated.end(), [findVal](const std::map<int64_t,int64_t>::value_type & p) {
                        return p.second == findVal;
                    });
                    cIndex_old2 = it->first;
                }

                int64_t cIndex_new1 = fConnDuplicated[cIndex_old1];
                int64_t cIndex_new2 = fConnDuplicated[cIndex_old2];

                TPZConnect &c2 = cmesh->ConnectVec()[cIndex_new1];
                auto *ptr2 = c2.FirstDepend();

                auto dep00 = ptr->HasDepend(cIndex_old2);
                auto dep01 = ptr->HasDepend(cIndex_new2);

                auto dep10 = ptr2->HasDepend(cIndex_old2);
                auto dep11 = ptr2->HasDepend(cIndex_new2);

                if (!dep00 || !dep01 || !dep10 || !dep11) DebugStop();

                int rows = dep00->fDepMatrix.Rows() + dep11->fDepMatrix.Rows();
                int cols = dep00->fDepMatrix.Cols() + dep11->fDepMatrix.Cols();
                
                TPZFMatrix<REAL> DepMat(rows,cols,0.);

                DepMat(0,0) = dep00->fDepMatrix(0,0);
                for (int i = 1; i < rows; i++){
                    DepMat(i,0) = dep10->fDepMatrix(i-1,0);
                    for (int j = 1; j < cols; j++){
                        if (i == 1) DepMat(0,j) = dep01->fDepMatrix(0,j-1);
                        DepMat(i,j) = dep11->fDepMatrix(i-1,j-1);
                    }                  
                }
                // std::cout << "K00 " << dep00->fDepMatrix << std::endl;
                // std::cout << "K01 " << dep01->fDepMatrix << std::endl;
                // std::cout << "K10 " << dep10->fDepMatrix << std::endl;
                // std::cout << "K11 " << dep11->fDepMatrix << std::endl;
                // std::cout << "K " << DepMat << std::endl;
                
                c.RemoveDepend();
                c2.RemoveDepend();

                //Only one dependency - old1 + old2
                TPZConnect::TPZDepend *depend = c.AddDependency(cIndex_old1, cIndex_old2, DepMat, 0, 0, rows, cols);
                               
                fDepMatTreated[cIndex_old1] = true;
                fDepMatTreated[cIndex_new1] = true;
                
                ptr = ptr->fNext;
                ptr = ptr->fNext;
            
            }
        }
    }


}



template class TPZHDivApproxSpaceCreator<STATE>;
// template class TPZHDivApproxSpaceCreator<CSTATE>;


