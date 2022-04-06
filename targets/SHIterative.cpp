#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

// #include "TPZGenGrid2D.h"

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZGmshReader.h>
#include <TPZVTKGeoMesh.h>
#include "Poisson/TPZMatPoisson.h" //for TPZMatLaplacian
#include "Projection/TPZL2Projection.h" //for BC in a single point
#include "pzmultiphysicscompel.h"
#include <TPZNullMaterial.h>
#include <TPZNullMaterialCS.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"// for Hdiv problem
#include <TPZBndCond.h> //for TPZBndCond
#include "TPZLinearAnalysis.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZLagrangeMultiplier.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzstrmatrixor.h"
#include "pzlog.h"
#include "pzshapecube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapepoint.h"
#include "pzshapetriang.h"

#include "divfree_config.h"
#include "TPZMatDivFreeBubbles.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZCompElKernelHDiv.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZKernelHdivUtils.h"
#include "TPZHDivApproxSpaceCreator.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "tpzsparseblockdiagonal.h"
#include "TPZFileStream.h"
#include "TPZAnalyticSolution.h"
#include "pzsmanal.h"
#include "pzsubcmesh.h"
#include "tpzverysparsematrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZSparseMatRed.h"
#include "pzsysmp.h"

void GetNSubEquations(TPZMultiphysicsCompMesh* cmesh, int64_t &nEqPres, int64_t &nEqFlux);

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

enum EMatid  {ENone, EDomain, EBottom, ERight, ETop, ELeft, EPont, EWrap, EIntface, EPressureHyb};

int main(int argc, char* argv[])
{
    //dimension of the problem
    constexpr int dim{2};
    constexpr int pOrder{1};
      

#ifdef PZ_LOG
TPZLogger::InitializePZLOG();
#endif
    
    //read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    {
        TPZGmshReader reader;
        // essa interface permite voce mapear os nomes dos physical groups para
        // o matid que voce mesmo escolher
        TPZManVector<std::map<std::string,int>,4> stringtoint(4);
        stringtoint[2]["Surface"] = 1;
        stringtoint[1]["Bottom"] = 2;
        stringtoint[1]["Right"] = 3;
        stringtoint[1]["Top"] = 4;
        stringtoint[1]["Left"] = 5;
        stringtoint[0]["Point"] = 6;
        stringtoint[1]["Top2"] = 7;
        reader.SetDimNamePhysical(stringtoint);
        reader.GeometricGmshMesh("../mesh/1element.msh",gmesh);
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    

    //............................Div Free Bubbles............................
    TPZCompMesh * cmeshflux = 0;
    TPZCompMesh * cmeshpressure = 0;
   
    TPZKernelHdivUtils<STATE> util;

    //Insert here the BC material id's to be hybridized
    std::set<int> matBCHybrid={};
    //Insert here the type of all boundary conditions
    std::set<int> matIDNeumann{};
    std::set<int> matIDDirichlet{ERight,ETop,EBottom,ELeft};
    /// All bc's mat ID's
    std::set<int> matBC;
    std::set_union(matIDNeumann.begin(),matIDNeumann.end(),matIDDirichlet.begin(),matIDDirichlet.end(),std::inserter(matBC, matBC.begin()));

    /// Creates the approximation space - Set the type of domain hybridization
    TPZHDivApproxSpaceCreator<STATE> createSpace(gmesh,TPZHDivApproxSpaceCreator<STATE>::ESemiHybrid);

    //Setting material ids
    createSpace.fConfig.fDomain = EDomain;
    createSpace.SetMaterialIds(EWrap,EPressureHyb,EIntface,EPont,matBCHybrid,matBC);
    createSpace.SetPOrder(pOrder+1);
    createSpace.Initialize();
    // util.PrintGeoMesh(gmesh);

    //Flux mesh
    TPZCompMesh * cmeshfluxNew = createSpace.CreateFluxCMesh();
    // std::cout << "FLUX \n";
    // util.PrintCMeshConnects(cmeshfluxNew);
    std::string fluxFile = "FluxCMesh";
    util.PrintCompMesh(cmeshfluxNew,fluxFile);

    //Pressure mesh
    TPZCompMesh * cmeshpressureNew = createSpace.CreatePressureCMesh();
    // std::cout << "PRESSURE \n";
    // util.PrintCMeshConnects(cmeshpressureNew);
    std::string pressureFile = "PressureCMesh";
    util.PrintCompMesh(cmeshpressureNew,pressureFile);

    TLaplaceExample1 LaplaceExact;
    LaplaceExact.fExact = TLaplaceExample1::EHarmonic;

    //Multiphysics mesh
    TPZManVector< TPZCompMesh *, 2> meshvectorNew(2);
    meshvectorNew[0] = cmeshfluxNew;
    meshvectorNew[1] = cmeshpressureNew;      
    auto * cmeshNew = createSpace.CreateMultiphysicsCMesh(meshvectorNew,LaplaceExact.ExactSolution(),matIDNeumann,matIDDirichlet);
    auto nEqNoCondense = cmeshNew->NEquations();
    // std::cout << "MULTIPHYSICS \n";
    // util.PrintCMeshConnects(cmeshNew);
    // Group and condense the elements
    createSpace.Condense(cmeshNew);
    std::string multiphysicsFile = "MultiPhysicsMeshNew";
    util.PrintCompMesh(cmeshNew,multiphysicsFile);


    //HERE STARTS THE ITERATIVE SOLVER SET
    //sets number of threads to be used by the solver
    constexpr int nThreads{0};
    // Solve the problem
    TPZLinearAnalysis anNew(cmeshNew,false);
    auto anNewSparse = anNew;
    // createSpace.Solve(anNew, cmeshNew, true, false); 

    // Compute the number of equations in the system
    int64_t nEqFull = cmeshNew->NEquations();
    int64_t nEqPres, nEqFlux;
    GetNSubEquations(cmeshNew,nEqPres,nEqFlux);
    std::cout << "NUMBER OF EQUATIONS:\n No condense = " << nEqNoCondense << 
                 "\n Condensed = " << nEqFull << 
                 "\n Flux = " << nEqFlux << 
                 " \n Pressure = " << nEqPres << std::endl;

    // Create the RHS vectors
    TPZFMatrix<STATE> rhsFull(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsFullSparse(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsAux(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsAuxSparse(nEqFull,1,0.);
    TPZFMatrix<STATE> rhsFlux(nEqFlux,1,0.);
    TPZFMatrix<STATE> rhsFluxSparse(nEqFlux,1,0.);
    TPZAutoPointer<TPZGuiInterface> guiInterface,guiInterfaceSparse;

    //Creates the problem matrix    
    TPZSkylineStructMatrix<STATE> Stiffness(cmeshNew);
    Stiffness.SetNumThreads(nThreads);
#ifdef PZ_USING_MKL
    TPZSpStructMatrix<STATE> StiffnessSparse(cmeshNew);
    StiffnessSparse.SetNumThreads(nThreads);
#endif
    for (auto &con : cmeshNew->ConnectVec())
    {
        int64_t seqNum = con.SequenceNumber();
        if (con.IsCondensed()) continue;
        if (seqNum < 0) continue;

        con.ResetElConnected();
        con.IncrementElConnected();
        
        if (seqNum < nEqPres){
            con.SetSequenceNumber(seqNum+nEqPres);
        } else {
            con.SetSequenceNumber(seqNum-nEqPres);
        }
        //Em cada caso precisa atualizar o tamanho do bloco
        int neq = con.NDof()*con.NState();
        seqNum=con.SequenceNumber();
        cmeshNew->Block().Set(seqNum,neq);
        // std::cout << "Connect = " << con << std::endl;
    }   
    cmeshNew->ExpandSolution();
    
    //Cria duas matrizes, para inverter a ordem das matrizes em bloco
    TPZSparseMatRed<STATE, TPZFMatrix<STATE>> *matRed;
    matRed = new TPZSparseMatRed<STATE, TPZFMatrix<STATE>>(nEqFull,nEqPres);
    TPZSparseMatRed<STATE, TPZVerySparseMatrix<STATE>> *matRedSparse;
    matRedSparse = new TPZSparseMatRed<STATE, TPZVerySparseMatrix<STATE>>(nEqFull,nEqPres);

    //Primeiro cria a matriz auxiliar
    TPZFMatrix<REAL> K00(nEqPres,nEqPres,0.);
    TPZVerySparseMatrix<REAL> K00Sparse(nEqPres,nEqPres,0.);

    TPZStepSolver<STATE> step, stepSparse;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    stepSparse.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    anNew.SetSolver(step);
    anNewSparse.SetSolver(stepSparse);

    //Altera range da matriz stiffness e cria a K00 correta no matRed correto.
    Stiffness.SetEquationRange(0,nEqPres);
    StiffnessSparse.SetEquationRange(0,nEqPres);

    //Transfere as submatrizes da matriz auxiliar para a matriz correta.
    step.SetMatrix(&K00);
    stepSparse.SetMatrix(&K00Sparse);
    matRed->SetSolver(&step);
    matRedSparse->SetSolver(&stepSparse);
    
    Stiffness.EquationFilter().Reset();
    StiffnessSparse.EquationFilter().Reset();
    anNew.SetStructuralMatrix(Stiffness);
    anNewSparse.SetStructuralMatrix(StiffnessSparse);

    //Monta a matriz auxiliar
    // matRed->K00()->Zero();
    rhsAux.Zero();
    rhsAuxSparse.Zero();
    Stiffness.Assemble(*matRed,rhsAux,guiInterface);
    StiffnessSparse.Assemble(*matRedSparse,rhsAuxSparse,guiInterfaceSparse);

    rhsFull=rhsAux;
    rhsFullSparse=rhsAuxSparse;
    
    // std::ofstream outMatRed("outMatRed.txt");
    // std::ofstream outMatRedSparse("outMatRedSparse.txt");
    // matRed->Print("MatRed",outMatRed,EMathematicaInput);
    // matRedSparse->Print("MatRedSparse",outMatRedSparse,EMathematicaInput);
    // matRed->Print("MATRED "); //Deve ser a matriz de pressão, depois fluxo.
    
    //Cria precondicionador bloco diagonal
    TPZBlockDiagonalStructMatrix<STATE> BDFmatrix(cmeshNew);
    BDFmatrix.SetEquationRange(nEqPres,nEqFull);
    // BDFmatrix.SetEquationRange(nEqFlux,nEqFull);
    // TPZSparseBlockDiagonal<REAL> KBD;
    TPZBlockDiagonal<REAL> KBD;
    BDFmatrix.CreateAssemble(rhsAux,guiInterface);
    BDFmatrix.EndCreateAssemble(&KBD);
    // KBD.Print("KBD");

    TPZFMatrix<REAL> *K11Red;
    TPZFMatrix<STATE> *K11RedSparse;
    K11Red = new TPZFMatrix<REAL>;
    K11RedSparse = new TPZFMatrix<REAL>;
    K11Red->Redim(nEqFlux,nEqFlux);
    K11RedSparse->Redim(nEqFlux,nEqFlux);
    
    matRed->SetF(rhsFull);
    matRedSparse->SetF(rhsFullSparse);
    
    matRed->K11Reduced(*K11Red,rhsFlux);
    matRedSparse->K11Reduced(*K11RedSparse,rhsFluxSparse);
    
    std::ofstream outMatRed("outMatRed.txt");
    std::ofstream outMatRedSparse("outMatRedSparse.txt");
    matRed->Print("MatRed",outMatRed,EMathematicaInput);
    matRedSparse->Print("MatRedSparse",outMatRedSparse,EMathematicaInput);
    
    std::ofstream outK11("outK11.txt");
    K11Red->Print("K11Red=",outK11,EMathematicaInput);
    std::ofstream outK11Sp("outK11sp.txt");
    K11RedSparse->Print("K11RedSparse=",outK11Sp,EMathematicaInput);

    matRed->F1Red(rhsFlux);
    matRed->SetF(rhsFull);
     
    TPZStepSolver<STATE> *precond = new TPZStepSolver<STATE>( &KBD );
    precond->SetDirect(ELU);
    int nMaxIter = 30;
    TPZVec<REAL> errors(nMaxIter);
    errors.Fill(0.);

    for (int64_t iter = 1; iter < nMaxIter; iter++){
        
        std::cout << "ITER = " << iter << std::endl;
        TPZFMatrix<STATE> residual(nMaxIter,1,0.);
        REAL tol = 1.e-10;
        REAL overrelax = 1.;
        TPZFMatrix<STATE> solution(nEqFlux,1,0.);
        TPZFMatrix<STATE> solution2(nEqFlux,1,0.);
        K11Red->SolveCG(iter,*precond,rhsFlux,solution,&residual,tol);
        // matRed->SolveCG(iter,*precond,rhsFlux,solution,&residual,tol);
        auto diff = solution2-solution;
        REAL norm = 0.;
        
        TPZFMatrix<STATE> press(nEqPres,1,0.);
        matRed->UGlobal(solution,press);

        //Update solution in Analysis
        for (int i = 0; i < nEqFlux; i++){
            rhsFull(nEqPres+i,0) = solution(i,0);
        }
        for (int i = 0; i < nEqPres; i++){
            rhsFull(i,0) = press(i,0);
        }

        anNew.Solution()=rhsFull;
        anNew.LoadSolution();
        
        //Compute error and print results
        anNew.SetExact(LaplaceExact.ExactSolution());
        //Print results
        util.PrintResultsMultiphysics(meshvectorNew,anNew,cmeshNew);

        // std::ofstream out4("mesh_MDFB.txt");
        // anNew.Print("nothing",out4);
        // std::ofstream anPostProcessFileMDFB("postprocessMDFB.txt");
                
        ///Calculating approximation error  
        TPZManVector<REAL,5> error;

        auto cmeshNew = anNew.Mesh();
        int64_t nelem = cmeshNew->NElements();
        cmeshNew->LoadSolution(cmeshNew->Solution());
        cmeshNew->ExpandSolution();
        cmeshNew->ElementSolution().Redim(nelem, 5);

        anNew.PostProcessError(error);
            
        // std::cout << "\nApproximation error:\n";
        // std::cout << "H1 Norm = " << std::scientific << std::setprecision(15) << error[0]<<'\n';
        // std::cout << "L1 Norm = " << std::scientific << std::setprecision(15) << error[1]<<'\n'; 
        // std::cout << "H1 Seminorm = " << std::scientific << std::setprecision(15) << error[2]<<'\n'; 

        if (tol < 1.e-10) break;
        errors[iter-1] = error[1];
        // std::cout << "error 4 = " << error[3]<<'\n'; 
        // std::cout << "error 5 = " << error[4] << "\n\n";
    }

    std::cout << "ERRORS = " << std::endl;
    for (int i = 0; i < errors.size(); i++){
        std::cout <<std::scientific << std::setprecision(15)<< errors[i] << std::endl;
    }

    return 0;
}

void GetNSubEquations(TPZMultiphysicsCompMesh* cmesh, int64_t &nEqPres, int64_t &nEqFlux)
{
    cmesh->LoadReferences();
    std::set<int> auxConnects;
    TPZVec<int64_t> activeEquationsP;
    TPZSkylineStructMatrix<STATE> Pmatrix(cmesh);
    TPZSkylineStructMatrix<STATE> Fmatrix(cmesh);

    auto gmesh = cmesh->Reference();

    for (auto gel:gmesh->ElementVec()){

        if (gel->MaterialId() != EPressureHyb) continue;

        auto nconnects = gel->Reference()->NConnects();
        for (size_t i = 0; i < nconnects; i++)
        {
            auxConnects.insert(gel->Reference()->ConnectIndex(i));
        }
    }

    int size = auxConnects.size();
    activeEquationsP.Resize(0);

    for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
        if (auxConnects.find(iCon) != auxConnects.end()) {
            auto &con = cmesh->ConnectVec()[iCon];
            if (con.HasDependency()){
                continue;
            }
            const auto seqnum = con.SequenceNumber();
            if (seqnum<0) continue;
            const auto pos = cmesh->Block().Position(seqnum);
            const auto blocksize = cmesh->Block().Size(seqnum);
            if (blocksize == 0){
                continue;
            }
            
            const auto vs = activeEquationsP.size();
            activeEquationsP.Resize(vs + blocksize);
            for (auto ieq = 0; ieq < blocksize; ieq++) {
                activeEquationsP[vs + ieq] = pos + ieq;
            }
        }
    }
    
    // std::cout << "SIZE = " << size << std::endl;
    const int neqs = activeEquationsP.size();
    // std::cout << "Active equations = " << activeEquationsP << std::endl;
        
    Pmatrix.EquationFilter().SetActiveEquations(activeEquationsP);
    

    //Flux and BD Flux matrices:
    int64_t sizeF = activeEquationsP[0];
    TPZVec<int64_t> activeEquationsF(sizeF);
    for (int i = 0; i < sizeF; i++)
    {
        activeEquationsF[i] = i;
    }
    Fmatrix.EquationFilter().SetActiveEquations(activeEquationsF);
    // BDFmatrix.EquationFilter().SetActiveEquations(activeEquationsF);
    
    
    nEqPres = Pmatrix.NReducedEquations();
    nEqFlux = Fmatrix.NReducedEquations();

}

