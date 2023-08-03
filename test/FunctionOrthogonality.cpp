/*
  This unit test verifies if the hybridization and semi hybridization techniques are working
  for any specified polynomial order and topology.
  
*/
#include <catch2/catch.hpp>
#include "TPZShapeHDivConstant.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "TPZShapeData.h"
#include "pzgmesh.h"
#include "TPZGeoMeshTools.h"
#include "pzcmesh.h"
#include "pzelchdiv.h"
#include "TPZMaterialDataT.h"
#include "TPZLapackEigenSolver.h"
std::ofstream rprint("results_genEigval.txt",std::ios_base::app);
std::ofstream rprint2("results_allgenEigval.txt",std::ios_base::app);
template<class TSHAPE>
void getShapeFunctions(TPZGeoEl *gel, TPZCompMesh *cmesh, int &pOrder, TPZFMatrix<REAL> &stiff);

int main(){
    

    int DIMENSION = 2;
    // int pOrder = 1;
    for (int pOrder = 1; pOrder < 16; pOrder++)
    {
       
        TPZGeoMesh* gmesh;
        if (DIMENSION == 2){
            gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(MMeshType::EQuadrilateral,1,false, -1);
        } else if (DIMENSION == 3){
            gmesh = TPZGeoMeshTools::CreateGeoMeshSingleEl(MMeshType::EHexahedral,1,false, -1);
        } else {
            DebugStop();
        }

        TPZGeoEl *gel = gmesh->ElementVec()[0];

        TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

        TPZFMatrix<REAL> stiff;
        if (DIMENSION == 2){
            getShapeFunctions<pzshape::TPZShapeQuad>(gel,cmesh,pOrder,stiff);
        } else if (DIMENSION == 3) {
            getShapeFunctions<pzshape::TPZShapeCube>(gel,cmesh,pOrder,stiff);
        } else {
            DebugStop();
        }
        
        // std::cout << "STIFF = " << stiff << std::endl;

        //Assemble the generalized eigenvalue problem to be solved - separate the matrix
        int nshape = stiff.Rows();
        TPZAutoPointer<TPZFMatrix<REAL>> matA = new TPZFMatrix<REAL>(nshape,nshape,0.);
        TPZAutoPointer<TPZFMatrix<REAL>> matB = new TPZFMatrix<REAL>(nshape,nshape,0.);
        
        int nrt0;
        if (gel->Dimension() == 2) nrt0 = 4;
        if (gel->Dimension() == 3) nrt0 = 6;
        for (int i = 0; i < nrt0; i++){
            for (int j = 0; j < nrt0; j++){
                matB->PutVal(i,j,stiff(i,j));
            }
            for (int j = nrt0; j < nshape; j++){
                matA->PutVal(i,j,stiff(i,j));
                matA->PutVal(j,i,stiff(j,i));
            }
        }
        for (int i = nrt0; i < nshape; i++){
            for (int j = nrt0; j < nshape; j++){
                matB->PutVal(i,j,stiff(i,j));
            }
        }
        
        // matA->Print("A=",std::cout,EMathematicaInput);
        // matB->Print("B=",std::cout,EMathematicaInput);

        TPZFMatrix<std::complex<REAL>> eigVectors(nshape,nshape,0.);
        TPZVec<std::complex<REAL>> eigValues(nshape,0.);

        TPZLapackEigenSolver<REAL> eigSolver;

        eigSolver.SetMatrixA(matA);
        eigSolver.SetMatrixB(matB);
        int a = eigSolver.SolveGeneralisedEigenProblem(eigValues,eigVectors);

        double maxEigval = 0.;
        double minEigval = 1000.;
        rprint2 << "pOrder = " << pOrder << std::endl;
        for (int i = 0; i < nshape; i++)
        {
            rprint2 << "Eigval [" << i<< "]= " << eigValues[i].real() << std::endl;
            if (fabs(eigValues[i].real()) > maxEigval ) maxEigval = eigValues[i].real();
            if (fabs(eigValues[i].real()) < fabs(minEigval) && fabs(eigValues[i].real()) > 0) 
                minEigval = eigValues[i].real();
        }
        rprint << "pOrder = " << pOrder << ", MaxEigVal = " << maxEigval << ", minEigVal = " << minEigval << std::endl;
    }
    return 0;
}

template<class TSHAPE>
void getShapeFunctions(TPZGeoEl *gel, TPZCompMesh *cmesh, int &pOrder, TPZFMatrix<REAL> &stiff){

    int dim = gel->Dimension();
    cmesh->SetDefaultOrder(pOrder);
    if (dim == 2){
        CreateHDivQuadEl(gel,*cmesh,HDivFamily::EHDivConstant);
    } else {
        CreateHDivCubeEl(gel,*cmesh,HDivFamily::EHDivConstant);
    }
    TPZCompEl *cel = cmesh->ElementVec()[0];
    TPZCompElHDiv<TSHAPE> *celhdiv = dynamic_cast<TPZCompElHDiv<TSHAPE>*> (cel);

    TPZIntPoints *intrule = gel->CreateSideIntegrationRule(gel->NSides()-1, pOrder+2);
    int npoints = intrule->NPoints();
    TPZManVector<REAL,3> xi(dim);
    TPZFNMatrix<9,REAL> jac(dim,dim),jacinv(dim,dim),axes(dim,3);
    REAL detjac;

    int nshape = 0;
    for (size_t i = 0; i < TSHAPE::NFacets+1; i++){
        nshape += TPZShapeHDivConstant<TSHAPE>::ComputeNConnectShapeF(i,pOrder);
    }

    stiff.Resize(nshape,nshape);
    stiff.Zero();

    for (int ip =0; ip<npoints; ip++) {
        REAL weight;
        TPZMaterialDataT<REAL> data;
        celhdiv->InitMaterialData(data);
        
        intrule->Point(ip, xi, weight);
        // gel->Jacobian(xi, jac, axes, detjac, jacinv);
        celhdiv->ComputeRequiredData(data, xi);
        celhdiv->ComputeShape(xi, data);

        // std::cout << "DeformedDirections = " << data.fDeformedDirections << std::endl;
        // std::cout << "DeformedDirections = " << data.phi << std::endl;
        
        double inv_perm=1.;
        //Calculate the matrix contribution for flux. Matrix A
        for (int iq = 0; iq < nshape; iq++) {
            //ef(iq, 0) += 0.;
            int ivecind = data.fVecShapeIndex[iq].first;
            int ishapeind = data.fVecShapeIndex[iq].second;
            TPZFNMatrix<3, REAL> ivec(3, 1, 0.);
            for (int id = 0; id < 3; id++) {
                ivec(id, 0) = data.fDeformedDirections(id, ivecind);
            }

            TPZFNMatrix<3, REAL> ivecZ(3, 1, 0.);
            TPZFNMatrix<3, REAL> jvecZ(3, 1, 0.);
            for (int jq = 0; jq < nshape; jq++) {
                TPZFNMatrix<3, REAL> jvec(3, 1, 0.);
                int jvecind = data.fVecShapeIndex[jq].first;
                int jshapeind = data.fVecShapeIndex[jq].second;

                for (int id = 0; id < 3; id++) {
                    jvec(id, 0) = data.fDeformedDirections(id, jvecind);
                }

                //dot product between Kinv[u]v
                jvecZ.Zero();
                for (int id = 0; id < 3; id++) {
                    jvecZ(id, 0) += inv_perm * jvec(id, 0);
                }
                double prod1 = ivec(0, 0) * jvecZ(0, 0) + ivec(1, 0) * jvecZ(1, 0) + ivec(2, 0) * jvecZ(2, 0);
                stiff(iq, jq) += weight * data.phi(ishapeind, 0) * data.phi(jshapeind, 0) * prod1;
            }
        }


    }
    delete intrule;
    

}




























