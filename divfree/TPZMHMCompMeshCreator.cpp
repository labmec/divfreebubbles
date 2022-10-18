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
#include "TPZHDivApproxCreator.h"

using namespace pzshape;

TPZMHMCompMeshCreator::TPZMHMCompMeshCreator(TPZMHMGeoMeshCreator &mhm_gcreator, TPZAutoPointer<TPZGeoMesh> &gmesh) : fGeoMeshCreator(mhm_gcreator){
    hDivCreator = new TPZHDivApproxCreator(gmesh.operator->());
};


TPZMultiphysicsCompMesh * TPZMHMCompMeshCreator::BuildMultiphysicsCMesh(int pOrder_vol, int pOrder_skel, TPZAutoPointer<TPZGeoMesh> &gmesh, TPZAnalyticSolution &analytic){

    hDivCreator->HdivFamily() = HDivFamily::EHDivConstant;
    hDivCreator->ProbType() = ProblemType::EDarcy;
    hDivCreator->IsRigidBodySpaces() = true;
    hDivCreator->SetDefaultOrder(pOrder_vol);
    hDivCreator->SetExtraInternalOrder(0);
    hDivCreator->SetShouldCondense(false);
    hDivCreator->HybridType() = HybridizationType::ENone;
    InsertMaterialObjects(analytic);

    int fNumMeshes = 4;
    TPZManVector<TPZCompMesh*,7> meshvec(fNumMeshes);
    int countMesh = 0;

    meshvec[countMesh++] = hDivCreator->CreateHDivSpace();

    TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(fGeoMeshCreator.fSkeletonMatId,1,1);
    meshvec[0]->InsertMaterialObject(nullmat);

    meshvec[0]->SetDefaultOrder(pOrder_skel);
    std::set<int> matSkel={fGeoMeshCreator.fSkeletonMatId};
    meshvec[0]->AutoBuild(matSkel);

    // skeleton elements have lagrange multiplier 4
    int64_t nel = meshvec[0]->NElements();
    for (int el = 0; el<nel; el++) {
        TPZCompEl *cel = meshvec[0]->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->MaterialId() == fGeoMeshCreator.fSkeletonMatId)
        {
            cel->Connect(0).SetLagrangeMultiplier(4);
        }
    }


    int lagLevelCounter = 1;
    meshvec[countMesh++] = hDivCreator->CreateL2Space(pOrder_vol,lagLevelCounter++);
    fDistFluxLevel = lagLevelCounter;
    meshvec[countMesh++] = hDivCreator->CreateConstantSpace(lagLevelCounter++);
    fAvPresLevel = lagLevelCounter;
    meshvec[countMesh++] = hDivCreator->CreateConstantSpace(lagLevelCounter++);    

    std::ofstream outflux("flux2.txt");
    meshvec[0]->Print(outflux);
    std::ofstream outpressure("pressure2.txt");
    meshvec[1]->Print(outpressure);
    std::ofstream outconstflux("constflux2.txt");
    meshvec[2]->Print(outconstflux);
    std::ofstream outconstpress("constpressure2.txt");
    meshvec[3]->Print(outconstpress);

    TPZMultiphysicsCompMesh *cmesh = hDivCreator->CreateMultiphysicsSpace(meshvec);


    return cmesh;
}


void TPZMHMCompMeshCreator::InsertMaterialObjects(TPZAnalyticSolution &analytic)
{

    int dim = 2;//cmesh.Dimension();
    auto *darcy = new TPZMixedDarcyFlow(fGeoMeshCreator.fDomainMatId,dim);
    // auto *darcy = new TPZMixedDarcyFlow(Emat2,dim);
    
    darcy->SetForcingFunction(analytic.ForceFunc(), 3);
    darcy->SetExactSol(analytic.ExactSolution(), 3);
    hDivCreator->InsertMaterialObject(darcy);
    // cmesh.InsertMaterialObject(darcy);
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZVec<STATE> val2(1,1.);

    for (auto mat:fGeoMeshCreator.fBoundMatId)
    {
        auto *bnd2 = darcy->CreateBC(darcy, mat, 0, val1, val2);
        bnd2->SetForcingFunctionBC(analytic.ExactSolution(), 3);
        hDivCreator->InsertMaterialObject(bnd2);
        // cmesh.InsertMaterialObject(bnd2);
    }

    auto *nullmat = new TPZNullMaterialCS<>(fGeoMeshCreator.fSkeletonMatId,dim-1,1);
    hDivCreator->InsertMaterialObject(nullmat);
    // cmesh.InsertMaterialObject(nullmat);

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

