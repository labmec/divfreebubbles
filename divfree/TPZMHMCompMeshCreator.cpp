#include "TPZMHMCompMeshCreator.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZAnalyticSolution.h"
#include "pzsubcmesh.h"
#include "pzcondensedcompel.h"
#include "Common.h"
#include "TPZHDivApproxSpaceCreator.h"
#include "TPZCompElHDivDuplConnects.h"
#include "TPZCompElHDivDuplConnectsBound.h"
#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzshapeprism.h"
#include "pzgeoelside.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzconnect.h"

using namespace pzshape;

TPZMHMCompMeshCreator::TPZMHMCompMeshCreator(TPZMHMGeoMeshCreator &mhm_gcreator, HDivFamily hdivfam) : fGeoMeshCreator(mhm_gcreator){
    fHdivFamily = hdivfam;
};


TPZMultiphysicsCompMesh * TPZMHMCompMeshCreator::BuildMultiphysicsCMesh(int pOrder_vol, int pOrder_skel, TPZAutoPointer<TPZGeoMesh> &gmesh, TPZAnalyticSolution &analytic){

    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(gmesh);
    auto flux = GenerateFluxMesh(gmesh,pOrder_vol,pOrder_skel);
    auto pressure = GeneratePressureMesh(gmesh, pOrder_vol);
    auto distflux = GenerateConstantMesh(gmesh, fDistFluxLevel);
    auto avpressure = GenerateConstantMesh(gmesh, fAvPresLevel);
    TPZCompMesh *LagMesh;
    if(fDuplConnects){
        LagMesh = GenerateLagranceCMeshDuplConnects(gmesh,10);
    }

    if(1)
    {
        std::ofstream outflux("flux.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(flux, outflux);
        std::ofstream outpressure("pressure.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(pressure, outpressure);
        std::ofstream outconstflux("constflux.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(distflux, outconstflux);
        std::ofstream outconstpress("constpressure.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(avpressure, outconstpress);
    }
    if(1)
    {
        std::ofstream outflux("flux.txt");
        flux->Print(outflux);
        std::ofstream outpressure("pressure.txt");
        pressure->Print(outpressure);
        std::ofstream outconstflux("constflux.txt");
        distflux->Print(outconstflux);
        std::ofstream outconstpress("constpressure.txt");
        avpressure->Print(outconstpress);
    }
    // TPZManVector<TPZCompMesh *> meshvec = {flux,pressure,distflux,avpressure,LagMesh};
    TPZManVector<TPZCompMesh *> meshvec = {flux,pressure,distflux,avpressure};
    // TPZManVector<TPZCompMesh *> meshvec = {flux,pressure};
    InsertMaterialObjects(*cmesh,analytic);
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    cmesh->BuildMultiphysicsSpace(meshvec);

    if (fDuplConnects){
        auto dim = cmesh->Dimension();
        auto mat3 = new TPZLagrangeMultiplierCS<STATE>(fInterface, dim-1);
        cmesh->InsertMaterialObject(mat3);

        std::set<int> matIdBCHyb = {fLagrange};
        fHybridizer.CreateInterfaceDuplConnects(cmesh,matIdBCHyb);
    }

    return cmesh;
}

/// generate a mesh with HDiv elements
TPZCompMesh * TPZMHMCompMeshCreator::GenerateFluxMesh(TPZAutoPointer<TPZGeoMesh> &gmesh, int pOrder_vol,int pOrder_skel)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = gmesh->Dimension();
    
    auto *mat = new TPZNullMaterial<>(fGeoMeshCreator.fDomainMatId,dim,1);
    cmesh->InsertMaterialObject(mat);

    for (auto matId:fGeoMeshCreator.fBoundMatId)
    {
        mat = new TPZNullMaterial<>(matId,dim-1,1);
        cmesh->InsertMaterialObject(mat);
    }
    mat = new TPZNullMaterial<>(fGeoMeshCreator.fSkeletonMatId,dim-1,1);
    cmesh->InsertMaterialObject(mat);

    // cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivStandard);
    cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    
    std::set<int> matvol = fGeoMeshCreator.fBoundMatId;
    matvol.insert(fGeoMeshCreator.fDomainMatId);
    std::set<int> matskel = {fGeoMeshCreator.fSkeletonMatId};
 
    if (fDuplConnects){
        //Configure an HDiv Approx Space Creator
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(fWrap);
        cmesh->InsertMaterialObject(mat);        
        
        fHybridizer.SetMaterialIds(fWrap,fLagrange,fInterface,fPoint,fGeoMeshCreator.fDomainMatId);
        std::set<int> BCHybrid = {};

        // std::cout << "Number of Elements = " << gmesh->NElements()<< std::endl;
        // fHybridizer.CreateWrapElements(gmesh.operator->(),BCHybrid,true,fHdivFamily);
        // std::cout << "Number of Elements2 = " << gmesh->NElements()<< std::endl;
        // // UpdateElementPartition();

        std::ofstream vtkfile("GMESH.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile);
        
        cmesh->SetDefaultOrder(pOrder_vol);
        CreateFluxDuplicatedConnects(cmesh);
        cmesh->SetDefaultOrder(pOrder_skel);
        CreateSkeletonDuplicatedConnects(cmesh);
        // cmesh->AutoBuild(matskel);

        ActivateDuplicatedConnects(cmesh);
        DisableDuplicatedConnects(cmesh);
        // fHybridizer.SemiHybridizeDuplConnects(cmesh,BCHybrid);

    } else {
        cmesh->SetDefaultOrder(pOrder_vol);
        cmesh->AutoBuild(matvol);
        cmesh->SetDefaultOrder(pOrder_skel);
        cmesh->AutoBuild(matskel);
        cmesh->InitializeBlock();
        
    }

    int64_t nc = cmesh->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        auto &c = cmesh->ConnectVec()[ic];
        c.SetLagrangeMultiplier(0);
    }

    // skeleton elements have lagrange multiplier 4
    int64_t nel = cmesh->NElements();
    for (int el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->MaterialId() == fGeoMeshCreator.fSkeletonMatId)
        {
            cel->Connect(0).SetLagrangeMultiplier(4);
        }
    }
    gmesh->ResetReference();

    
    return cmesh;
}

/// generate a mesh with L2 elements
TPZCompMesh * TPZMHMCompMeshCreator::GeneratePressureMesh(TPZAutoPointer<TPZGeoMesh> &gmesh, int pOrder)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

    switch (fHdivFamily)
    {
    case HDivFamily::EHDivConstant:
        cmesh->SetAllCreateFunctionsDiscontinuous();
        cmesh->SetDefaultOrder(0);
        break;
    case HDivFamily::EHDivStandard:
        cmesh->SetDefaultOrder(pOrder);
        cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
        break;
    
    default:
        DebugStop();
        break;
    }  

    int dim = gmesh->Dimension();
    auto *mat = new TPZNullMaterial<>(fGeoMeshCreator.fDomainMatId,dim,1);
    cmesh->InsertMaterialObject(mat);
    
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    int64_t ncon = cmesh->ConnectVec().NElements();
    for (int64_t ic = 0; ic<ncon; ic++) {
        cmesh->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }

    // if (fDuplConnects){
    //     auto *mat2 = new TPZNullMaterial<>(fLagrange,dim-1,1);
    //     cmesh->InsertMaterialObject(mat2);

    //     std::set<int> lag ={fLagrange};
    //     cmesh->AutoBuild(lag);
    //     for(auto &newnod : cmesh->ConnectVec())
    //     {
    //         newnod.SetLagrangeMultiplier(10);
    //     }
    //     cmesh->InitializeBlock();
    //     cmesh->ExpandSolution();
    // }
    
//    std::set<int64_t> partitionprocessed;
//    int64_t nel = cmesh->NElements();
//    for (int el = 0; el<nel; el++) {
//        TPZCompEl *cel = cmesh->Element(el);
//        TPZGeoEl *gel = cel->Reference();
//        int64_t partition = elpartition[gel->Index()];
//        if(partitionprocessed.find(partition) == partitionprocessed.end())
//        {
//            partitionprocessed.insert(partition);
//            cel->Connect(0).SetLagrangeMultiplier(5);
//        }
//    }
    gmesh->ResetReference();
    return cmesh;
}

/// generate a mesh with constant elements
TPZCompMesh * TPZMHMCompMeshCreator::GenerateConstantMesh(TPZAutoPointer<TPZGeoMesh> &gmesh, int level)
{
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(0);
    int dim = gmesh->Dimension();
    auto *mat = new TPZNullMaterial<>(fGeoMeshCreator.fDomainMatId,dim,1);
    cmesh->InsertMaterialObject(mat);
    cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        auto index = gel->Index();
        cel->Connect(0).SetLagrangeMultiplier(level);
    }
    gmesh->ResetReference();
    return cmesh;

}


void TPZMHMCompMeshCreator::InsertMaterialObjects(TPZMultiphysicsCompMesh &cmesh,TPZAnalyticSolution &analytic)
{

    auto gmesh = cmesh.Reference();
    int dim = gmesh->Dimension();
    auto *darcy = new TPZMixedDarcyFlow(fGeoMeshCreator.fDomainMatId,dim);
    // auto *darcy = new TPZMixedDarcyFlow(Emat2,dim);
    
    darcy->SetForcingFunction(analytic.ForceFunc(), 3);
    darcy->SetExactSol(analytic.ExactSolution(), 3);
    cmesh.InsertMaterialObject(darcy);
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZVec<STATE> val2(1,1.);

    for (auto mat:fGeoMeshCreator.fBoundMatId)
    {
        auto *bnd2 = darcy->CreateBC(darcy, mat, 0, val1, val2);
        bnd2->SetForcingFunctionBC(analytic.ExactSolution(), 3);
        cmesh.InsertMaterialObject(bnd2);
    }

    auto *nullmat = new TPZNullMaterialCS<>(fGeoMeshCreator.fSkeletonMatId,dim,1);
    cmesh.InsertMaterialObject(nullmat);

    if (fDuplConnects){
        auto * nullmat2 = new TPZNullMaterialCS<>(fWrap,dim-1,1);
        cmesh.InsertMaterialObject(nullmat2);

        auto * nullmat3 = new TPZNullMaterialCS<>(fLagrange,dim-1,1);
        cmesh.InsertMaterialObject(nullmat3);
    }

}

void TPZMHMCompMeshCreator::PutinSubstructures(TPZCompMesh &cmesh)
{
    // create subcompmeshes
    std::map<int64_t,TPZSubCompMesh *> submeshes;
    for(auto part : fGeoMeshCreator.fElementPartition) {
        if(part == -1) continue;
        if(submeshes.find(part) == submeshes.end())
        {
            auto *sub = new TPZSubCompMesh(cmesh);
            submeshes[part] = sub;
        }
    }
    int64_t nel = cmesh.NElements();
    for (int64_t el = 0; el<nel; el++) {
        auto *cel = cmesh.Element(el);
        auto *sub = dynamic_cast<TPZSubCompMesh *>(cel);
        if(sub) continue;
        auto *gel = cel->Reference();
        auto index = gel->Index();
        auto part = fGeoMeshCreator.fElementPartition[index];
        if(part >= 0) {
            auto *sub = submeshes[part];
            sub->TransferElement(&cmesh, el);
        }
    }
    cmesh.ComputeNodElCon();
    for(auto itsub : submeshes) {
        TPZCompEl *submesh = itsub.second;
        int64_t ncon = submesh->NConnects();
        for(int64_t ic = 0; ic<ncon; ic++) {
            auto &c = submesh->Connect(ic);
            if(c.LagrangeMultiplier() == fAvPresLevel) {
                c.IncrementElConnected();
                break;
            }
        }
    }
    for(auto itsub : submeshes) {
        TPZSubCompMesh *submesh = itsub.second;
        submesh->MakeAllInternal();
        submesh->CleanUpUnconnectedNodes();
    }
    cmesh.ComputeNodElCon();
    {
        int64_t ncon = cmesh.NConnects();
        for(int64_t ic = 0; ic<ncon; ic++) {
            auto &c = cmesh.ConnectVec()[ic];
            if(c.NElConnected() == 0 && c.HasDependency()) {
                c.RemoveDepend();
            }
        }
    }
    cmesh.CleanUpUnconnectedNodes();
    if(1)
    {
        std::ofstream out("cmesh_substruct.txt");
        cmesh.Print(out);
        std::ofstream out2("cmesh_substruct.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(&cmesh, out2);
    }
}

void TPZMHMCompMeshCreator::CondenseElements(TPZCompMesh &cmesh)
{
    cmesh.ComputeNodElCon();
    int64_t nel = cmesh.NElements();
    cmesh.ElementSolution().Resize(nel, 5);
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        auto *subcel = dynamic_cast<TPZSubCompMesh *>(cel);
        if(subcel) {
            CondenseElements(*subcel);
            subcel->CleanUpUnconnectedNodes();
            subcel->SetAnalysisSparse(0);
        //    subcel->SetAnalysisSkyline();
        }
        else if(cel)
        {
            TPZGeoEl *gel = cel->Reference();
            if(gel && gel->Dimension() == cmesh.Dimension()) {
                // prevent the pressure connect from being condensed
                int nc = cel->NConnects();
                for (int ic = 0; ic<nc; ic++) {
                    TPZConnect &c = cel->Connect(ic);
                    if(c.LagrangeMultiplier() == fAvPresLevel) {
                        c.IncrementElConnected();
                    }
                }
                auto *cond = new TPZCondensedCompEl(cel);
            }
        }
    }
    auto *sub = dynamic_cast<TPZSubCompMesh *>(&cmesh);
    if(0 && !sub)
    {
        std::ofstream out("CondensedCompMesh.txt");
        cmesh.Print(out);
    }
}

void TPZMHMCompMeshCreator::CreateFluxDuplicatedConnects(TPZCompMesh *cmesh){

    int64_t nel = cmesh->Reference()->NElements();
    int fDimension = cmesh->Dimension();

    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = cmesh->Reference()->Element(el);
        // TPZCompEl *cel = cmesh->Element(el);

        if(!gel) DebugStop();
        auto type = gel -> Type();
        auto matid = gel->MaterialId();
        if (matid != fGeoMeshCreator.fDomainMatId) continue;
        if(gel->HasSubElement()) continue;
        
        using namespace pzshape;
        
        // First: create the volumetric element. Notice that there is no difference EHCurlNoGrads and EHDivKernel in 2D.
        if (fDimension == 2){
            switch (type){
            case ETriangle:
                CreateHDivDuplConnectsTriangleEl(gel,*cmesh,fHdivFamily);
                break;
            case EQuadrilateral:
                CreateHDivDuplConnectsQuadEl(gel,*cmesh,fHdivFamily);
                break;
            default:
                DebugStop();
                break;
            }
        } else if (fDimension == 3){
            switch (type){
            case ECube:
                CreateHDivDuplConnectsCubeEl(gel,*cmesh,fHdivFamily);
                break;
            case ETetraedro:
                CreateHDivDuplConnectsTetraEl(gel,*cmesh,fHdivFamily);
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
            auto neighMatid=neighbour.Element()->MaterialId();

            if(fGeoMeshCreator.fBoundMatId.find(neighMatid) == fGeoMeshCreator.fBoundMatId.end() && neighMatid != fWrap) continue;
            if (gelside.Dimension() != gel->Dimension()-1) continue;

            //Creates the computational elements. Both wrap, BC and interface
            if (fDimension == 2){
                CreateHDivDuplConnectsBoundLinearEl(neighbour.Element(),*cmesh,fHdivFamily);
            } else if (fDimension == 3){
                switch (type){
                case ETetraedro:
                    CreateHDivDuplConnectsBoundTriangEl(neighbour.Element(),*cmesh,fHdivFamily);
                    break;
                case ECube:
                    CreateHDivDuplConnectsBoundQuadEl(neighbour.Element(),*cmesh,fHdivFamily);
                    break;
                default:
                    DebugStop();
                    break;
                }
            }
        }
    } 

    cmesh->InitializeBlock();

    
}


TPZCompMesh * TPZMHMCompMeshCreator::GenerateLagranceCMeshDuplConnects(TPZAutoPointer<TPZGeoMesh> &gmesh, int level)
{
    std::set<int> matIdVec = {fLagrange};
    
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    auto dim = gmesh->Dimension();

    // Sets matid to BC geometric elements
    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        mat->SetDimension(dim);
        cmesh->InsertMaterialObject(mat);
    }

    cmesh->SetDefaultOrder(0);
    cmesh->SetDimModel(dim-1);
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    cmesh->AutoBuild(matIdVec);
    
    for(auto &newnod : cmesh->ConnectVec())
    {
        newnod.SetLagrangeMultiplier(level);
    }
    cmesh->InitializeBlock();
    cmesh->ExpandSolution();

    return cmesh;
}


void TPZMHMCompMeshCreator::UpdateElementPartition(TPZAutoPointer<TPZGeoMesh> &gmesh){


    // we have to update this data structure after creating the geometrical elements for the hybridization: 

    // fGeoMeshCreator.fElementPartition;
}


void TPZMHMCompMeshCreator::ActivateDuplicatedConnects(TPZCompMesh *cmesh){
    //Prints computational mesh properties
    // std::string vtk_name = "CMESH_DELETE.vtk";
    // std::ofstream vtkfile(vtk_name.c_str());
    // TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

    for (int64_t i = 0; i < cmesh->NElements(); i++)
    {
        TPZCompElHDivDuplConnects<TPZShapeTriang> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeTriang> *> (cmesh->Element(i)); 
        TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *> (cmesh->Element(i)); 
        if (celd) celd->ActiveDuplConnects(fConnDuplicated);
        if (celb) celb->ActiveDuplConnects(fConnDuplicated);
    
    }

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

    cmesh->InitializeBlock(); 
    cmesh->ExpandSolution();

}


void TPZMHMCompMeshCreator::DisableDuplicatedConnects(TPZCompMesh *cmesh){

    for (int64_t i = 0; i < cmesh->NElements(); i++)
    {
        TPZCompElHDivDuplConnects<TPZShapeTriang> *celd = dynamic_cast<TPZCompElHDivDuplConnects<TPZShapeTriang> *> (cmesh->Element(i)); 
        TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *celb = dynamic_cast<TPZCompElHDivDuplConnectsBound<TPZShapeLinear> *> (cmesh->Element(i)); 
        if (celd) celd->InactiveDuplConnects();
        if (celb) celb->InactiveDuplConnects();
    
    }

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
                
                c.RemoveDepend(cIndex_old1,cIndex_old2);
                c.RemoveDepend(cIndex_old1,cIndex_new2);
                c2.RemoveDepend(cIndex_new1,cIndex_old2);
                c2.RemoveDepend(cIndex_new1,cIndex_new2);

                //Only one dependency - old1 + old2
                TPZConnect::TPZDepend *depend = c.AddDependency(cIndex_old1, cIndex_old2, DepMat, 0, 0, rows, cols);
                               
                fDepMatTreated[cIndex_old1] = true;
                fDepMatTreated[cIndex_new1] = true;
                
                ptr = ptr->fNext; 
            }
        }
    }

    cmesh->InitializeBlock(); 
    cmesh->ExpandSolution();
}

void TPZMHMCompMeshCreator::CreateSkeletonDuplicatedConnects(TPZCompMesh *cmesh){

    int64_t nel = cmesh->Reference()->NElements();
    int fDimension = cmesh->Dimension();

    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = cmesh->Reference()->Element(el);
        // TPZCompEl *cel = cmesh->Element(el);

        if(!gel) DebugStop();
        auto type = gel -> Type();
        auto matid = gel->MaterialId();
        if (matid != fGeoMeshCreator.fSkeletonMatId) continue;
        if(gel->HasSubElement()) continue;
        
        using namespace pzshape;
        
        // First: create the volumetric element. Notice that there is no difference EHCurlNoGrads and EHDivKernel in 2D.
        if (fDimension == 3){
            switch (type){
            case ETriangle:
                CreateHDivDuplConnectsBoundTriangEl(gel,*cmesh,fHdivFamily);
                break;
            case EQuadrilateral:
                CreateHDivDuplConnectsBoundQuadEl(gel,*cmesh,fHdivFamily);
                break;
            default:
                DebugStop();
                break;
            }
        } else if (fDimension == 2){
            switch (type){
            case EOned:
                CreateHDivDuplConnectsBoundLinearEl(gel,*cmesh,fHdivFamily);
                break;
            default:
                DebugStop();
                break;
            }
        }
    } 

    cmesh->InitializeBlock();

    
}