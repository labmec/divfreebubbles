#ifdef HAVE_CONFIG_H
  #include <pz_config.h>
#endif

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
#include "../headers/TPZMatDivFreeBubbles.h" //THE NEW MATERIAL!
#include "Poisson/TPZMatPoisson.h" //for TPZMatLaplacian
#include "Projection/TPZL2Projection.h" //for BC in a single point
#include <TPZNullMaterial.h>
#include "DarcyFlow/TPZMixedDarcyFlow.h"// for Hdiv problem
#include <TPZBndCond.h> //for TPZBndCond
#include "TPZLinearAnalysis.h"
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include "TPZMultiphysicsCompMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzcompel.h"
#include "TPZInterfaceEl.h"
#include "pzstrmatrixor.h"
// #include "pzlog.h"

TPZCompMesh *FluxCMesh(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMesh(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *MultiphysicCMesh(int dim, int pOrder, int *matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
TPZCompMesh *CMeshH1(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *CMeshDivFreeBubbles(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh);
void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsH1(int dim, TPZLinearAnalysis &an);
void PrintResultsDivFreeBubbles(int dim, TPZLinearAnalysis &an);
void ComputeError(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);
void ComputeErrorHdiv(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile);
void CompareSolution(int dim, int pOrder, int *matIdVec, TPZCompMesh *cmeshHdiv,TPZCompMesh *cmeshDFB,TPZGeoMesh * gmesh);

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

int main(int argc, char* argv[]){
  //dimension of the problem
  constexpr int dim{2};
  constexpr int pOrder{2};
  //Materials - See the .geo file
  int matIdVec[]={1,2,3,4,5,6,7,8};
  //1 = Injection Well
  //2 = Production Well
  //3 = Bottom Line
  //4 = Top Line
  //5 = Left Line
  //6 = Right Line
  //7 = Domain
  //8 = Upper Left point

// #ifdef PZ_LOG
//   TPZLogger::InitializePZLOG();
// #endif
  
  //read mesh from gmsh
  TPZGeoMesh *gmesh;
  gmesh = new TPZGeoMesh();
  TPZGmshReader *reader;
  reader = new TPZGmshReader();
  reader -> GeometricGmshMesh4("../mesh/five-spot.msh",gmesh);
  
  //.................................Hdiv.................................
  //Flux mesh
  TPZCompMesh * cmeshflux= FluxCMesh(dim,pOrder,matIdVec,gmesh);

  //Pressure mesh
  TPZCompMesh * cmeshpressure= PressureCMesh(dim,pOrder,matIdVec,gmesh);

  //Multiphysics mesh
  TPZManVector< TPZCompMesh *, 2> meshvector(2);
  meshvector[0] = cmeshflux;
  meshvector[1] = cmeshpressure;
  TPZCompMesh * cmesh = MultiphysicCMesh(dim,pOrder,matIdVec,meshvector,gmesh);

  //Solve Multiphysics
  TPZLinearAnalysis an(cmesh,true);
  SolveProblem(an,cmesh);

  //Print results
  PrintResultsMultiphysic(dim,meshvector,an,cmesh);

  //..................................H1..................................
  //Creates H1 problem
  // TPZCompMesh * cmeshH1 = CMeshH1(dim,pOrder,matIdVec,gmesh);

  // //Solve H1
  // TPZLinearAnalysis anH1(cmeshH1,false);
  // SolveProblem(anH1,cmeshH1);

  // //Print results
  // PrintResultsH1(dim,anH1);

  //...........................Div Free Bubbles...........................
  //Creates DFB problem
  TPZCompMesh * cmeshDFB = CMeshDivFreeBubbles(dim,pOrder,matIdVec,gmesh);

  //Solve DFB
  TPZLinearAnalysis anDFB(cmeshDFB,false);
  SolveProblem(anDFB,cmeshDFB);

  //Print results
  PrintResultsDivFreeBubbles(dim,anDFB);
  
  //...........................ERROR EVALUATION...........................
  // std::ofstream anPostProcessFileHdiv("postprocessHdiv.txt");
  // ComputeErrorHdiv(an,anPostProcessFileHdiv);

  // std::ofstream anPostProcessFileDFB("postprocessDFB.txt");
  // ComputeError(anDFB,anPostProcessFileDFB);

  // std::ofstream anPostProcessFileH1("postprocessH1.txt");
  // ComputeError(anH1,anPostProcessFileH1);

  CompareSolution(dim,pOrder,matIdVec,cmesh,cmeshDFB,gmesh);


  return 0;
}

//Analytical solution
constexpr int solOrder{2};
auto exactSol = [](const TPZVec<REAL> &loc,
  TPZVec<STATE>&u,
  TPZFMatrix<STATE>&gradU){
  const auto &x=loc[0];
  const auto &y=loc[1];
  const auto &d = 1.; // distance between injection and production wells
  u[0]= log(hypot(x,y)) - log(hypot(x-d,y-d)) - log(hypot(x+d,y-d)) - log(hypot(x-d,y+d)) - log(hypot(x+d,y+d));
  gradU(0,0) = x/(x*x+y*y) - (x-d)/(pow(x-d,2)+pow(y-d,2)) - (x+d)/(pow(x+d,2)+pow(y-d,2)) - (x-d)/(pow(x-d,2)+pow(y+d,2)) - (x+d)/(pow(x+d,2)+pow(y+d,2));
  gradU(1,0) = y/(x*x+y*y) - (y-d)/(pow(x-d,2)+pow(y-d,2)) - (y-d)/(pow(x+d,2)+pow(y-d,2)) - (y+d)/(pow(x-d,2)+pow(y+d,2)) - (y+d)/(pow(x+d,2)+pow(y+d,2));

  // u[0]= x*x-y*y;
  // gradU(0,0) = 2.*x;
  // gradU(1,0) = -2.*y;

  // gradU(2,0) = 0;//optional
};

//Flux computational mesh
TPZCompMesh *FluxCMesh(int dim, int pOrder,int *matIdVec, TPZGeoMesh *gmesh) 
{
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  TPZNullMaterial<> *mat = new TPZNullMaterial<>(matIdVec[6]);
//  mat->NStateVariables();
  cmesh->InsertMaterialObject(mat);

  cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
  mat->SetDimension(dim);
  mat -> fBigNumber = 1.e10;
  //Boundary Conditions
  TPZFMatrix<STATE> val1(1,1,1.);
  TPZManVector<STATE> val2(1,-4.3820281971726605);
  auto * BCond0 = mat->CreateBC(mat, matIdVec[0], 1, val1, val2);
  TPZFMatrix<STATE> val3(1,1,1.); 
  TPZManVector<STATE> val4(1,0.9288056258294208);
  auto * BCond1 = mat->CreateBC(mat, matIdVec[1], 1, val3, val4);
  BCond0->SetForcingFunctionBC(exactSol);
  BCond1->SetForcingFunctionBC(exactSol);
  cmesh->InsertMaterialObject(BCond1);
  cmesh->InsertMaterialObject(BCond0);

  TPZFMatrix<STATE> val5(1,1,1.);
  TPZManVector<STATE> val6(1,0.);
  auto * BCond2 = mat->CreateBC(mat, matIdVec[2], 1, val5, val6);
  auto * BCond3 = mat->CreateBC(mat, matIdVec[3], 1, val5, val6);
  auto * BCond4 = mat->CreateBC(mat, matIdVec[4], 1, val5, val6);
  auto * BCond5 = mat->CreateBC(mat, matIdVec[5], 0, val5, val6);
  BCond2->SetForcingFunctionBC(exactSol);
  BCond3->SetForcingFunctionBC(exactSol);
  BCond4->SetForcingFunctionBC(exactSol);
  BCond5->SetForcingFunctionBC(exactSol);
  cmesh->InsertMaterialObject(BCond2);
  cmesh->InsertMaterialObject(BCond3);
  cmesh->InsertMaterialObject(BCond4);
  cmesh->InsertMaterialObject(BCond5);

  cmesh->SetDefaultOrder(pOrder);
  cmesh->AutoBuild();

  // Print flux mesh
  // std::ofstream myfile("FluxMesh.txt");
  // cmesh->Print(myfile);

  return cmesh;
}

// Pressure computational mesh
TPZCompMesh *PressureCMesh(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh)
{
  //Change if not triangular elements
  bool fTriang = false;
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  TPZNullMaterial<> *mat = new TPZNullMaterial<>(matIdVec[6]);
  mat->SetDimension(dim);
  cmesh->InsertMaterialObject(mat);
  mat -> fBigNumber = 1.e10;
  cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetDefaultOrder(pOrder);
  cmesh->SetDimModel(dim);
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
      celdisc->SetConstC(1.);
      celdisc->SetTrueUseQsiEta();
      if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
      {
          if(fTriang==true) celdisc->SetTotalOrderShape();
          else celdisc->SetTensorialShape();
      }
  }
  //Print pressure mesh
  // std::ofstream myfile("PressureMesh.txt");
  // cmesh->Print(myfile);

  return cmesh;
}

// Multiphysics computational mesh
TPZCompMesh *MultiphysicCMesh(int dim, int pOrder, int *matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
  gmesh->ResetReference();
  auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
  cmesh->SetDefaultOrder(pOrder);
  cmesh->SetDimModel(dim);
  auto mat = new TPZMixedDarcyFlow(matIdVec[6], dim);

  mat->SetPermeabilityFunction(1.);
  cmesh->InsertMaterialObject(mat);
  mat -> fBigNumber = 1.e10;
    
  //Boundary Conditions
  TPZFMatrix<STATE> val1(1,1,1.);
  TPZManVector<STATE> val2(1,-4.3820281971726605);
  auto * BCond0 = mat->CreateBC(mat, matIdVec[0], 1, val1, val2);
  TPZFMatrix<STATE> val3(1,1,1.); 
  TPZManVector<STATE> val4(1,0.9288056258294208);
  auto * BCond1 = mat->CreateBC(mat, matIdVec[1], 1, val3, val4);
  BCond0->SetForcingFunctionBC(exactSol);
  BCond1->SetForcingFunctionBC(exactSol);
  cmesh->InsertMaterialObject(BCond1);
  cmesh->InsertMaterialObject(BCond0);

  TPZFMatrix<STATE> val5(1,1,1.);
  TPZManVector<STATE> val6(1,0.);
  auto * BCond2 = mat->CreateBC(mat, matIdVec[2], 1, val5, val6);
  auto * BCond3 = mat->CreateBC(mat, matIdVec[3], 1, val5, val6);
  auto * BCond4 = mat->CreateBC(mat, matIdVec[4], 1, val5, val6);
  auto * BCond5 = mat->CreateBC(mat, matIdVec[5], 0, val5, val6);
  BCond2->SetForcingFunctionBC(exactSol);
  BCond3->SetForcingFunctionBC(exactSol);
  BCond4->SetForcingFunctionBC(exactSol);
  BCond5->SetForcingFunctionBC(exactSol);
  cmesh->InsertMaterialObject(BCond2);
  cmesh->InsertMaterialObject(BCond3);
  cmesh->InsertMaterialObject(BCond4);
  cmesh->InsertMaterialObject(BCond5);

  TPZManVector<int> active(2,1);
  cmesh->BuildMultiphysicsSpace(active, meshvector);
  cmesh->SetAllCreateFunctionsMultiphysicElem();
  cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();

  TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);

  // Prints Multiphysics mesh
  std::ofstream myfile("MultiPhysicsMesh.txt");
  cmesh->Print(myfile);

  return cmesh;
}

void SolveProblem(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
  //sets number of threads to be used by the solver
  constexpr int nThreads{4};
  //defines storage scheme to be used for the FEM matrices
  //in this case, a symmetric skyline matrix is used
  TPZSkylineStructMatrix<STATE> matskl(cmesh);
  matskl.SetNumThreads(nThreads);
  an.SetStructuralMatrix(matskl);

  ///Setting a direct solver
  TPZStepSolver<STATE> step;
  ///Setting an iterative solver
  // TPZMatrixSolver<STATE> * precond = an.BuildPreconditioner(TPZAnalysis::EBlockJacobi , true);
  // TPZCopySolve<STATE> * precond = new TPZCopySolve<STATE>( matskl.Create() );  step.ShareMatrix( *precond );
  // TPZStepSolver<STATE> * precond = new TPZStepSolver<STATE>( matskl.Create() ); step.ShareMatrix( *precond ); precond->SetJacobi(1, 0.0, 0);
  // TPZStepSolver<STATE> jac;
  // REAL tol = 1.e-10;
  // jac.SetSSOR(1,1.1,0.,0);
  // jac.ShareMatrix(step);
  // an.SetSolver(step);
  // step.SetGMRES(20,20, *precond, tol, 0);
  // step.SetCG(2000, *precond, tol, 0);
  step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
  an.SetSolver(step);

  //assembles the system
  an.Assemble();

  ///solves the system
  an.Solve();

  return;
}

void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{

  an.SetExact(exactSol,solOrder);

  TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
	TPZManVector<std::string,10> scalnames(2), vecnames(2);
    
	scalnames[0] = "Pressure";
	scalnames[1] = "ExactPressure";
	vecnames[0]= "Flux";
	vecnames[1]= "ExactFlux";

	int div = 0;
  std::string plotfile = "solutionHdiv.vtk";
	an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	an.PostProcess(div,dim);
  // Print mesh properties
	// std::ofstream out("mesh.txt");
	// an.Print("nothing",out);

  return;
}

TPZCompMesh *CMeshH1(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh)
{
  //Creates cmesh object
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetAllCreateFunctionsContinuous();

  //Sets materials
  TPZMatPoisson<> *mat = new TPZMatPoisson<>(matIdVec[6],dim);
//  mat->SetPermeability(1);
  cmesh->InsertMaterialObject(mat);

  //Insert boundary conditions
  TPZFMatrix<STATE> val1(1,1,1.);
  TPZManVector<STATE> val2(2,-4.3820281971726605);
  auto * BCond = mat->CreateBC(mat, matIdVec[0], 0, val1, val2);
  TPZFMatrix<STATE> val3(1,1,1.);
  TPZManVector<STATE> val4(2,0.9288056258294208);
  auto * BCond1 = mat->CreateBC(mat, matIdVec[1], 0, val3, val4);
  BCond->SetForcingFunctionBC(exactSol);
  BCond1->SetForcingFunctionBC(exactSol);
  cmesh->InsertMaterialObject(BCond);
  cmesh->InsertMaterialObject(BCond1);
  
  cmesh->SetDefaultOrder(pOrder);
  cmesh->AutoBuild();

  return cmesh;
}

void PrintResultsH1(int dim, TPZLinearAnalysis &an)
{
  an.SetExact(exactSol,solOrder);
  TPZVec<std::string> scalarVars(1), vectorVars(1);
  scalarVars[0] = "Solution";
  vectorVars[0] = "Derivative";
  an.DefineGraphMesh(dim,scalarVars,vectorVars,"SolutionH1.vtk");
  constexpr int resolution{0};
  an.PostProcess(resolution,dim);	

  return;
}

void PrintResultsDivFreeBubbles(int dim, TPZLinearAnalysis &an)
{
  an.SetExact(exactSol,solOrder);
  TPZVec<std::string> vectorVars(1), scalarVars(1);
  scalarVars[0] = "Solution";
  vectorVars[0] = "Derivative";
  an.DefineGraphMesh(dim,scalarVars,vectorVars,"SolutionDFB.vtk");
  constexpr int resolution{0};
  an.PostProcess(resolution,dim);	

  return;
}

void ComputeError(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile)
{
  an.SetExact(exactSol,solOrder);
  ///Calculating approximation error  
  TPZManVector<REAL,3> error;
  an.PostProcess(error,anPostProcessFile);
	
  std::cout << "\nApproximation error:\n";
  std::cout << "H1 Norm = " << error[0]<<'\n';
  std::cout << "L1 Norm = " << error[1]<<'\n'; 
  std::cout << "H1 Seminorm = " << error[2] << "\n\n";

  return;
}

void ComputeErrorHdiv(TPZLinearAnalysis &an, std::ofstream &anPostProcessFile)
{
  an.SetExact(exactSol,solOrder);
  ///Calculating approximation error  
  TPZManVector<REAL,5> error;
  an.PostProcess(error,anPostProcessFile);
	
  std::cout << "\nApproximation error:\n";
  std::cout << "H1 Norm = " << error[0]<<'\n';
  std::cout << "L1 Norm = " << error[1]<<'\n'; 
  std::cout << "H1 Seminorm = " << error[2] << "\n\n";

  return;
}


TPZCompMesh *CMeshDivFreeBubbles(int dim, int pOrder, int *matIdVec, TPZGeoMesh *gmesh)
{
  //Creates cmesh object
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetAllCreateFunctionsContinuous();

  //Sets materials
  auto *mat = new TPZMatDivFreeBubbles<STATE>(matIdVec[6],dim);
  cmesh->InsertMaterialObject(mat);
   
  //Insert boundary conditions
  TPZFMatrix<STATE> val1(1,1,1.);
  TPZManVector<STATE> val2(1,-4.3820281971726605);
  TPZManVector<STATE> val5(1,0.);
  constexpr int boundType{0};
  auto * BCond = mat->CreateBC(mat, matIdVec[0], 1, val1, val5);//Injection
  TPZFMatrix<STATE> val3(1,1,1.);
  TPZManVector<STATE> val4(1,0.9288056258294208);
  auto * BCond1 = mat->CreateBC(mat, matIdVec[1], 1, val3, val5);//Production
  BCond->SetForcingFunctionBC(exactSol);
  BCond1->SetForcingFunctionBC(exactSol);
  cmesh->InsertMaterialObject(BCond);
  cmesh->InsertMaterialObject(BCond1);
  mat -> fBigNumber = 1.e10;
  
  auto * BCond2 = mat->CreateBC(mat, matIdVec[2], 1, val3, val5);//Bottom
  auto * BCond3 = mat->CreateBC(mat, matIdVec[3], 1, val3, val5);//Top
  auto * BCond4 = mat->CreateBC(mat, matIdVec[4], 1, val3, val5);//Left
  auto * BCond5 = mat->CreateBC(mat, matIdVec[5], 0, val3, val5);//Right
  BCond2->SetForcingFunctionBC(exactSol);
  BCond3->SetForcingFunctionBC(exactSol);
  BCond4->SetForcingFunctionBC(exactSol);
  BCond5->SetForcingFunctionBC(exactSol);
  cmesh->InsertMaterialObject(BCond2);
  cmesh->InsertMaterialObject(BCond3);
  cmesh->InsertMaterialObject(BCond4);
  cmesh->InsertMaterialObject(BCond5);

  auto *mat2 = new TPZL2Projection<>(matIdVec[7],0,1);
  cmesh->InsertMaterialObject(mat2);

  //Prints computational mesh properties
  // std::stringstream text_name;
  // std::stringstream vtk_name;
  // text_name   << "geometry" << ".txt";
  // vtk_name    << "geometry" << ".vtk";
  // std::ofstream textfile(text_name.str().c_str());
  // gmesh->Print(textfile);
  // std::ofstream vtkfile(vtk_name.str().c_str());
  // TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);

  cmesh->SetDefaultOrder(pOrder);
  cmesh->AutoBuild();

  // Prints DFB mesh
  // std::ofstream myfile("DFBMesh.txt");
  // CMeshDFB->Print(myfile);

  return cmesh;
}

void CompareSolution(int dim, int pOrder, int *matIdVec, TPZCompMesh *cmeshHdiv,TPZCompMesh *cmeshDFB, TPZGeoMesh * gmesh)
{
  //Copied from TPZAnalysis

  TPZManVector<REAL,10> values(10,0.);
	TPZManVector<REAL,10> values2(10,0.);
  TPZManVector<REAL,10> values2DFB(10,0.);
  TPZSolutionMatrix fSolutionHdiv, fSolutionDFB;
	fSolutionHdiv = cmeshHdiv->Solution();
  fSolutionDFB = cmeshDFB->Solution();
  TPZAdmChunkVector<TPZCompEl *> &elvecHdiv = cmeshHdiv->ElementVec();
  TPZAdmChunkVector<TPZCompEl *> &elvecDFB = cmeshDFB->ElementVec();
	TPZManVector<REAL,10> errors(10);
  TPZManVector<REAL,10> errorsDFB(10);
	errors.Fill(0.0);
	int64_t nel = elvecHdiv.NElements();
	int matId0 = 0;
	for(int i=0;i<nel;i++) {
        TPZCompEl *cel = elvecHdiv[i];
        if(!cel) continue;
		matId0=cel->Material()->Id();
		if(matId0 > 0)
			break;
	}
	bool lastEl=false;
	for(int i=0;i<nel;i++) {
    TPZCompEl *el = (TPZCompEl *) elvecHdiv[i];


    int nc = el->NConnects();
    for (int ic = 0; ic < nc; ic++) {
      TPZManVector<STATE> connectsol;
      int64_t cindex = el->ConnectIndex(ic);
      TPZConnect &c = cmeshHdiv->ConnectVec()[cindex];
      cmeshHdiv->ConnectSolution<STATE>(cindex, cmeshHdiv, cmeshHdiv->Solution(), connectsol);
    }

		if(el) {
			errors.Fill(0.0);
      TPZGeoEl *gel = el->Reference();
      if(gel->Dimension() != cmeshHdiv->Dimension()) continue;
      auto *matError =
          dynamic_cast<TPZMatError<STATE>*>(el->Material());
      if (matError || matError->HasExactSol()) {
        el->EvaluateError(errors, 0);
      } else {
          PZError<<__PRETTY_FUNCTION__;
          PZError<<" the material has no associated exact solution\n";
          PZError<<"Aborting...";
          DebugStop();
      }
      if(matId0==el->Material()->Id()){
				for(int ier = 0; ier < errors.NElements(); ier++) 	values[ier] += errors[ier] * errors[ier];
				lastEl=false;
			}
			else {
				for(int ier = 0; ier < errors.NElements(); ier++)	values2[ier] += errors[ier] * errors[ier];
				lastEl=true;
			}
		}
	}
  
  nel = elvecDFB.NElements();
	matId0 = 0;
  lastEl=false;
	for(int i=0;i<nel;i++) {
    TPZCompEl *elDFB = (TPZCompEl *) elvecDFB[i];

    int nc = elDFB->NConnects();
    for (int ic = 0; ic < nc; ic++) {
      TPZManVector<STATE> connectsol;
      int64_t cindex = elDFB->ConnectIndex(ic);
      TPZConnect &c = cmeshDFB->ConnectVec()[cindex];
      cmeshDFB->ConnectSolution<STATE>(cindex, cmeshDFB, cmeshDFB->Solution(), connectsol);
    }

		if(elDFB) {
			errors.Fill(0.0);
      TPZGeoEl *gel = elDFB->Reference();
      
      auto *matError =
          dynamic_cast<TPZMatError<STATE>*>(elDFB->Material());
      if (matError || matError->HasExactSol()) {
        elDFB->EvaluateError(errorsDFB, 0);
        TPZCompMesh *mesh = elDFB->Mesh();
        TPZFMatrix<STATE> &elementSol = mesh->ElementSolution();
      } else {
          PZError<<__PRETTY_FUNCTION__;
          PZError<<" the material has no associated exact solution\n";
          PZError<<"Aborting...";
          DebugStop();
      }
      if(matId0==elDFB->Material()->Id()){
				for(int ier = 0; ier < errors.NElements(); ier++) 	values[ier] += errors[ier] * errors[ier];
				lastEl=false;
			}
			else {
        for(int ier = 0; ier < errorsDFB.NElements(); ier++)	values2DFB[ier] += errorsDFB[ier] * errorsDFB[ier];
				lastEl=true;
			}
		}
	}
  std::cout << "ERROR Hdiv" << scientific << std::cout.precision(15) << std::endl; 
  // std::cout << "ERROR[0] = " << sqrt(values2[0]) << std::endl; 
  // std::cout << "ERROR[1] = " << sqrt(values2[1]) << std::endl; 
  // std::cout << "ERROR[2] = " << sqrt(values2[2]) << std::endl; 
  std::cout << "ERROR[3] = " << sqrt(values2[3]) << std::endl; 
  // std::cout << "ERROR[4] = " << sqrt(values2[4]) << std::endl; 
  std::cout << "\nERROR DFB" << scientific << std::cout.precision(15) << std::endl; 
  // std::cout << "ERROR[0] = " << sqrt(values2DFB[0]) << std::endl; 
  // std::cout << "ERROR[1] = " << sqrt(values2DFB[1]) << std::endl; 
  std::cout << "ERROR[2] = " << sqrt(values2DFB[2]) << std::endl; 

}