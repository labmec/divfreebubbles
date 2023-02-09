#include "TPZMixedElasticityCMeshCreator.h"
#include "TPZNullMaterial.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "TPZCompElDiscScaled.h"
#include <Elasticity/TPZMixedElasticityND.h>
#include "pzintel.h"
#include "TPZAnalyticSolution.h"

TPZCompMesh * TPZMixedElasticityCMeshCreator::CMesh_S(TPZGeoMesh *gmesh, int pOrder, HDivFamily &hdivfam) {
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    int dim = gmesh->Dimension();
    cmesh->SetDefaultOrder(pOrder); //Default polynomial order of the approximation
    cmesh->SetDimModel(dim); //Dimesion of the model

    //Definition of the approximation space:
    cmesh->ApproxSpace().SetHDivFamily(hdivfam);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);

    //Criando material cujo nSTATE = 2:
    TPZNullMaterial<STATE> * material = new TPZNullMaterial<STATE>(fDomainMatId);
    material->SetNStateVariables(dim);
    material->SetDimension(dim);
    cmesh->InsertMaterialObject(material); //Insere material na malha

    //Boundary conditions:
    //Dirichlet Boundary Conditions
    TPZFMatrix<STATE> val1(2, 2, 0.);
    TPZManVector<STATE> val2(2, 0.);
    for (auto matId : fBCDirichlet)
    {
        TPZBndCondT<STATE> * BCond = material->CreateBC(material, matId, 0, val1, val2);
        cmesh->InsertMaterialObject(BCond);
    }

    for (auto matId : fBCNeumann)
    {
        TPZBndCondT<STATE> * BCond = material->CreateBC(material, matId, 1, val1, val2);
        cmesh->InsertMaterialObject(BCond);
    }
    
    for (auto matId : fBCDirichletVar)
    {
        DebugStop();
    }

    for (auto matId : fBCMixed)
    {
        DebugStop();
    }

    for (auto matId : fBCPoint)
    {
        DebugStop();
    }

    cmesh->InsertMaterialObject(material->CreateBC(material, matLagrange, 1, val1, val2)); //Insere material na malha

    //Criando elementos computacionais que gerenciarão o espaco de aproximacao da malha:
    int ncel = cmesh->NElements();
    for (int i = 0; i < ncel; i++) {
        TPZCompEl * compEl = cmesh->ElementVec()[i];
        if (!compEl) continue;
        TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *> (compEl);
        if (facel)DebugStop();
    }

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;

}



TPZCompMesh* TPZMixedElasticityCMeshCreator::CMesh_U(TPZGeoMesh *gmesh, int pOrder) {
    
    int dim = gmesh->Dimension();
    if (pOrder == 0) {
        // Constant-per-element
        TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
        cmesh->SetDimModel(dim); //Insere dimensão do modelo

        cmesh->SetAllCreateFunctionsDiscontinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);

        
        TPZNullMaterial<> * material = new TPZNullMaterial<>(fDomainMatId);
        material->SetDimension(dim);
        material->SetNStateVariables(dim);

        cmesh->InsertMaterialObject(material); //Insere material na malha

        std::set<int> materialids;
        materialids.insert(fDomainMatId);
        {
            gmesh->ResetReference();
            int64_t nel = gmesh->NElements();
            for (int64_t el = 0; el < nel; el++) {
                TPZGeoEl *gel = gmesh->Element(el);
                if (!gel)continue;
                int matid = gel->MaterialId();
                if (materialids.find(matid) == materialids.end()) {
                    continue;
                }
                new TPZCompElDiscScaled(*cmesh, gel);
                gel->ResetReference();
            }
        }

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
            disc->SetTrueUseQsiEta();
        }

        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
        return cmesh;
    } else {
        //Criando malha computacional:
        TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
        cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
        cmesh->SetDimModel(dim); //Insere dimensão do modelo

        cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
        cmesh->ApproxSpace().CreateDisconnectedElements(true);


        TPZNullMaterial<> * material = new TPZNullMaterial<>(fDomainMatId);
        material->SetDimension(dim);
        material->SetNStateVariables(dim);

        cmesh->InsertMaterialObject(material); //Insere material na malha




        //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha

        int64_t ncel = cmesh->NElements();
        for (int64_t i = 0; i < ncel; i++) {
            TPZCompEl * compEl = cmesh->ElementVec()[i];
            if (!compEl) continue;
            TPZInterfaceElement * facel = dynamic_cast<TPZInterfaceElement *> (compEl);
            if (facel)DebugStop();

        }
        std::set<int> materialids;
        materialids.insert(fDomainMatId);
        cmesh->AutoBuild(materialids);
        cmesh->LoadReferences();
        cmesh->ApproxSpace().CreateDisconnectedElements(false);
        cmesh->AutoBuild();

        int64_t ncon = cmesh->NConnects();
        for (int64_t i = 0; i < ncon; i++) {
            TPZConnect &newnod = cmesh->ConnectVec()[i];
            newnod.SetLagrangeMultiplier(1);
        }
        for (int64_t i = 0; i < ncel; i++) {
            TPZCompEl * compEl = cmesh->ElementVec()[i];
            if (!compEl) continue;
            compEl->Connect(0).SetLagrangeMultiplier(3);

        }


        //    cmesh->AdjustBoundaryElements();
        //    cmesh->CleanUpUnconnectedNodes();

        return cmesh;
    }
}



TPZCompMesh* TPZMixedElasticityCMeshCreator::CMesh_P(TPZGeoMesh *gmesh, int pOrder, REAL elementdim) {
    
    int dim = gmesh->Dimension();
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

    cmesh->SetAllCreateFunctionsDiscontinuous();

    //    cmesh->SetAllCreateFunctionsContinuous(); //Criando funções H1
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    //Criando material cujo nSTATE = 1:
    TPZNullMaterial<> *material = new TPZNullMaterial<>(fDomainMatId); //criando material que implementa a formulacao fraca do problema modelo
    material->SetDimension(dim);
    if(dim == 3)
    {
        material->SetNStateVariables(3);
    }

    cmesh->InsertMaterialObject(material); //Insere material na malha
    std::set<int> materialids;
    materialids.insert(fDomainMatId);
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
    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    return cmesh;

}



TPZCompMesh* TPZMixedElasticityCMeshCreator::CMesh_RigidBody(TPZGeoMesh *gmesh, int lagrange) {
    //Criando malha computacional:
    int dim = gmesh->Dimension();
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(0); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

    cmesh->SetAllCreateFunctionsDiscontinuous();

    //Criando material cujo nSTATE = 1:
    TPZNullMaterial<> *material = new TPZNullMaterial<>(fDomainMatId); //criando material que implementa a formulacao fraca do problema modelo
    material->SetDimension(dim);
    if(dim == 2)
    {
        material->SetNStateVariables(3);
    }
    else if(dim == 3)
    {
        material->SetNStateVariables(6);
    }

    cmesh->InsertMaterialObject(material); //Insere material na malha
    std::set<int> materialids;
    materialids.insert(fDomainMatId);
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
            TPZCompElDisc *disc = new TPZCompElDisc(*cmesh, gel);
            disc->SetFalseUseQsiEta();
            gel->ResetReference();
        }
    }
    int ncon = cmesh->NConnects();
    for (int i = 0; i < ncon; i++) {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(lagrange);
    }

    cmesh->CleanUpUnconnectedNodes();
    cmesh->ExpandSolution();
    return cmesh;

}


TPZCompMesh* TPZMixedElasticityCMeshCreator::CMesh_Girk(TPZGeoMesh *gmesh, int pOrder) {

    int dim = gmesh->Dimension();
    //Criando malha computacional:
    int bc_inte_order = 10;
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    //    TElasticityExample1 example;

    // Criando material:


    REAL E = 20.59; //* @param E elasticity modulus
    REAL nu = 0.; //* @param nu poisson coefficient


    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = -32.69; //* @param fx forcing function \f$ -x = fx \f$
    int plain = 1.; //* @param plainstress = 1 \f$ indicates use of plainstress

    TPZMixedElasticityND * material = new TPZMixedElasticityND(1, E, nu, fx, fy, plain, dim);
    //TPZMixedElasticityMaterial * material2 = new TPZMixedElasticityMaterial(3,E,nu,fx,fy,plain,dim);
    //material1->SetAxisSymmetric();
    //material2->SetAxisSymmetric();
    //material->SetForcingFunction(example.ForcingFunction());
    // Inserindo material na malha
    //    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (f_source);
    //    TPZAutoPointer<TPZFunction<STATE> > pp = new TPZDummyFunction<STATE> (p_exact);
    //    TPZAutoPointer<TPZFunction<STATE> > vp = new TPZDummyFunction<STATE> (v_exact);

    //    material->SetForcingFunction(fp);
    //    material->SetForcingFunctionExactPressure(pp);
    //    material->SetForcingFunctionExact(vp);

    cmesh->InsertMaterialObject(material);
    //cmesh->InsertMaterialObject(material2);

    //Condições de contorno:

    TPZFMatrix<REAL> val1(2, 2, 0.);
    TPZManVector<REAL> val2(2, 0.);
    REAL x;
    val1(0, 0) = 0;
    val1(1, 1) = material->BigNumber();

    for (auto matId : fBCDirichlet)
    {
        TPZBndCondT<STATE> * BCond = material->CreateBC(material, matId, 0, val1, val2);
        cmesh->InsertMaterialObject(BCond);
    }

    for (auto matId : fBCNeumann)
    {
        TPZBndCondT<STATE> * BCond = material->CreateBC(material, matId, 1, val1, val2);
        cmesh->InsertMaterialObject(BCond);
    }
    
    for (auto matId : fBCDirichletVar)
    {
        DebugStop();
    }

    for (auto matId : fBCMixed)
    {
        DebugStop();
    }

    for (auto matId : fBCPoint)
    {
        DebugStop();
    }




    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;


}



TPZCompMesh* TPZMixedElasticityCMeshCreator::CMesh_m(TPZGeoMesh *gmesh, int pOrder, TPZAnalyticSolution * gAnalytic) {

    int dim = gmesh->Dimension();
    //Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo
    cmesh->SetAllCreateFunctionsMultiphysicElem();

    // Criando material:

    // example is initialized in the calling method
    //    example.fProblemType = TElasticityExample1::EThiago;
    //    example.fStressState   = TElasticityExample1::EPlaneStrain;

    REAL E = 250.; //* @param E elasticity modulus
    REAL nu = 0.25; //* @param nu poisson coefficient

    
    TElasticity2DAnalytic * analytic2D = 0;
    TElasticity3DAnalytic * analytic3D = 0;
    if(dim == 2)
    {
        analytic2D = dynamic_cast<TElasticity2DAnalytic *>(gAnalytic);
        if(!analytic2D) DebugStop();
        TPZManVector<REAL,3> x(3,0.);
        // analytic2D->Elastic(x, E, nu);
    }
    else if(dim == 3)
    {
        analytic3D = dynamic_cast<TElasticity3DAnalytic *>(gAnalytic);
        if(!analytic3D) DebugStop();
        TPZManVector<REAL,3> x(3,0.);
        // analytic3D->Elastic(x, E, nu);
    }

    REAL fx = 0.; //* @param fx forcing function \f$ -x = fx \f$
    REAL fy = 0.; //* @param fx forcing function \f$ -x = fx \f$
    int plainStress = 0; //* @param plainstress = 1 \f$ indicates use of plainstress
    if(dim == 2)
    {
        if (analytic2D->fPlaneStress == 0) {
            plainStress = 0;
        } else {
            plainStress = 1;
        }
    }
    TPZMixedElasticityND * material = new TPZMixedElasticityND(fDomainMatId, E, nu, fx, fy, plainStress, dim);

    // if (TElasticityExample1::fStressState == TElasticityExample1::EAxiSymmetric) {
    //     material->SetAxisSymmetric();
    // }
    material->SetForcingFunction(gAnalytic->ForceFunc(),3);
    // Inserindo material na malha
    //    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE> (f_source);
    //    TPZAutoPointer<TPZFunction<STATE> > pp = new TPZDummyFunction<STATE> (p_exact);
    //    TPZAutoPointer<TPZFunction<STATE> > vp = new TPZDummyFunction<STATE> (v_exact);

    //    material->SetForcingFunction(fp);
    //    material->SetForcingFunctionExactPressure(pp);
    //    material->SetForcingFunctionExact(vp);

    cmesh->InsertMaterialObject(material);


    //Condições de contorno:

    TPZFMatrix<REAL> val1(dim, dim, 0.);
    TPZManVector<REAL> val2(dim, 0.);
    
    for (auto matId : fBCDirichlet)
    {
        TPZBndCondT<STATE> * BCond = material->CreateBC(material, matId, 0, val1, val2);
        BCond->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
        cmesh->InsertMaterialObject(BCond);
    }

    for (auto matId : fBCNeumann)
    {
        TPZBndCondT<STATE> * BCond = material->CreateBC(material, matId, 1, val1, val2);
        BCond->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
        cmesh->InsertMaterialObject(BCond);
    }
    
    for (auto matId : fBCDirichletVar)
    {
        DebugStop();
    }

    for (auto matId : fBCMixed)
    {
        DebugStop();
    }

    for (auto matId : fBCPoint)
    {
        DebugStop();
    }

    auto * BCond4 = material->CreateBC(material, matLagrange, 1, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond4->SetForcingFunctionBC(gAnalytic->ExactSolution(),3);
    cmesh->InsertMaterialObject(BCond4); //Insere material na malha

    //Ponto

    //    TPZFMatrix<REAL> val3(1,1,0.), val4(1,1,0.);
    //    val4(0,0)=0.0;
    //
    //    TPZMaterial * BCPoint = material->CreateBC(material, matPoint, pointtype, val3, val4); //Cria material que implementa um ponto para a pressão
    //    cmesh->InsertMaterialObject(BCPoint); //Insere material na malha

    //Criando elementos computacionais que gerenciarão o espaco de aproximação da malha:

    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;

}


/// change the order of the internal connect to the given order
void TPZMixedElasticityCMeshCreator::ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder) {
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;

        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != cmesh->Dimension()) {
            continue;
        }
        int nc = cel->NConnects();
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel) DebugStop();

        intel->ForceSideOrder(gel->NSides() - 1, pOrder);
    }
    cmesh->ExpandSolution();
}
