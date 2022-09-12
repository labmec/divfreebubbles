#include "TPZMHMCompMeshCreator.h"
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
    TPZManVector<TPZCompMesh *> meshvec = {flux,pressure,distflux,avpressure};
    // TPZManVector<TPZCompMesh *> meshvec = {flux,pressure};
    InsertMaterialObjects(*cmesh,analytic);
    cmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    cmesh->BuildMultiphysicsSpace(meshvec);

    if (fDuplConnects){
        std::set<int> matIdBCHyb;
        fHDivCreator->Hybridizer().CreateInterfaceDuplConnects(cmesh,matIdBCHyb);

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
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(15);
        cmesh->InsertMaterialObject(mat);

        fHDivCreator = new TPZHDivApproxSpaceCreator<double>(cmesh->Reference(),TPZHDivApproxSpaceCreator<double>::ENone,HDivFamily::EHDivConstant);
        std::set<int> BCNull;
        fHDivCreator->SetMaterialIds(15,16,17,18,BCNull,fGeoMeshCreator.fBoundMatId);
        fHDivCreator->fConfig.fDomain = Emat2;
        fHDivCreator->SetPOrder(pOrder_vol);
        fHDivCreator->Initialize();

        std::ofstream out2("GMESH.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fHDivCreator->GetGeoMesh(), out2);
        
        CreateFluxDuplicatedConnects(cmesh);
        // fHDivCreator->DuplicateInternalConnects(cmesh);
        // fHDivCreator->ExpandConnects(cmesh);
        cmesh->InitializeBlock();
        cmesh->ExpandSolution();

        cmesh->SetDefaultOrder(pOrder_skel);
        cmesh->AutoBuild(matskel);

    } else {
        cmesh->SetDefaultOrder(pOrder_vol);
        cmesh->AutoBuild(matvol);
        cmesh->SetDefaultOrder(pOrder_skel);
        cmesh->AutoBuild(matskel);
    }

    cmesh->InitializeBlock();
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
    cmesh->SetDefaultOrder(pOrder);

    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->SetDefaultOrder(0);

    int dim = gmesh->Dimension();
    auto *mat = new TPZNullMaterial<>(fGeoMeshCreator.fDomainMatId,dim,1);
    cmesh->InsertMaterialObject(mat);
    // cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    // cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    int64_t ncon = cmesh->ConnectVec().NElements();
    for (int64_t ic = 0; ic<ncon; ic++) {
        cmesh->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    
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
    // auto *darcy = new TPZMixedDarcyFlow(fGeoMeshCreator.fDomainMatId,dim);
    auto *darcy = new TPZMixedDarcyFlow(Emat2,dim);
    
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
        auto * nullmat2 = new TPZNullMaterialCS<>(15,dim-1,1);
        cmesh.InsertMaterialObject(nullmat2);

        auto * nullmat3 = new TPZNullMaterialCS<>(16,dim-1,1);
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
            // subcel->SetAnalysisSparse(0);
           subcel->SetAnalysisSkyline();
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

            if(fGeoMeshCreator.fBoundMatId.find(neighMatid) == fGeoMeshCreator.fBoundMatId.end()) continue;
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

    //Prints computational mesh properties
    std::string vtk_name = "CMESH_DELETE.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

    // DuplicateInternalConnects(cmesh);

    // cmesh->InitializeBlock(); 
    // cmesh->ExpandSolution();
    // cmesh->ComputeNodElCon();
    // cmesh->LoadReferences();

    // hybridizer.SemiHybridizeDuplConnects(cmesh,fConfig.fBCHybridMatId);
    
}