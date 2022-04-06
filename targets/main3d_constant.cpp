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
#include "DarcyFlow/TPZIsotropicPermeability.h"
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
#include "pzshapetetra.h"
#include <pzfstrmatrix.h>


#include "divfree_config.h"
#include "TPZMatDivFreeBubbles.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZCompElKernelHDiv.h"
#include "DarcyFlow/TPZMixedDarcyFlow.h"
#include "TPZKernelHdivUtils.h"
#include "TPZApproxSpaceKernelHdiv.h"
#include <TPZGeoMeshTools.h>
#include "TPZCompElHCurlNoGrads.h"
#include "TPZMatCurlDotCurl.h"
#include "TPZHCurlEquationFilter.h"
#include "TPZCompMeshTools.h"


TPZCompMesh *FluxCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZCompMesh *PressureCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh);
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
TPZGeoMesh *CreateGeoMesh(const MMeshType meshType, const TPZVec<int> &nDivs, const int volId, const int bcId);
TPZGeoMesh *CreateGeoMeshTetra(const MMeshType meshType, const TPZVec<int> &nDivs, const int volId, const int bcId);
TPZCompMesh *FluxCMeshCurl(int dim, int pOrder, TPZGeoMesh *gmesh);
void PrintEdges(TPZGeoMesh* cmesh, std::set<int64_t> &rem_edges, std::map<int64_t, TPZHCurlEquationFilter<STATE>::VertexFilter> &vertexData, std::map<int64_t, TPZHCurlEquationFilter<STATE>::EdgeFilter> &edgeData);

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
using namespace std;

//Analytical solution
constexpr int solOrder{2};
auto exactSol = [](const TPZVec<REAL> &loc,
    TPZVec<STATE>&u,
    TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];

    // u[0] = x*x-y*y;//x*x*x*y*z - y*y*y*x*z;
    // gradU(0,0) = 2.*x;//(3.*x*x*y*z - y*y*y*z);
    // gradU(1,0) = -2.*y;//(x*x*x*z - 3.*y*y*x*z);
    // gradU(2,0) = 0.;//(x*x*x*y - y*y*y*x);

    // u[0] = x*x+y*y+z*z;//x*x*x*y*z - y*y*y*x*z;
    // gradU(0,0) = 2.*x;//(3.*x*x*y*z - y*y*y*z);
    // gradU(1,0) = 2.*y;//(x*x*x*z - 3.*y*y*x*z);
    // gradU(2,0) = 2.*z;//(x*x*x*y - y*y*y*x);

    // u[0] = x*x*x*y*z - y*y*y*x*z;
    // gradU(0,0) = (3.*x*x*y*z - y*y*y*z);
    // gradU(1,0) = (x*x*x*z - 3.*y*y*x*z);
    // gradU(2,0) = (x*x*x*y - y*y*y*x);

    REAL aux = 1./sinh(sqrt(2)*M_PI);
    u[0] = sin(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2)*M_PI*z)*aux;
    gradU(0,0) = -M_PI*cos(M_PI*x)*sin(M_PI*y)*sinh(sqrt(2)*M_PI*z)*aux;
    gradU(1,0) = -M_PI*cos(M_PI*y)*sin(M_PI*x)*sinh(sqrt(2)*M_PI*z)*aux;
    gradU(2,0) = -sqrt(2)*M_PI*cosh(sqrt(2)*M_PI*z)*sin(M_PI*x)*sin(M_PI*y)*aux;
};


enum EMatid {ENone, EDomain, ESurfaces, EPont, EWrap, EIntface, EPressureHyb, EEdgeRemove};

int main(int argc, char* argv[])
{
    //dimension of the problem
    constexpr int dim{3};
    constexpr int pOrder{5};
      

#ifdef PZ_LOG
TPZLogger::InitializePZLOG();
#endif
    
    // //read mesh from gmsh
    // TPZGeoMesh *gmesh;
    // gmesh = new TPZGeoMesh();
    // {
    //     TPZGmshReader reader;
    //     // essa interface permite voce mapear os nomes dos physical groups para
    //     // o matid que voce mesmo escolher
    //     TPZManVector<std::map<std::string,int>,4> stringtoint(4);
    //     stringtoint[3]["Domain"] = 1;
    //     stringtoint[2]["Surfaces"] = 2;
    //     // stringtoint[2]["Hybrid"] = 3;
    //     reader.SetDimNamePhysical(stringtoint);
    //     reader.GeometricGmshMesh("../mesh/span_tree.msh",gmesh);
    //     std::ofstream out("gmesh.vtk");
    //     TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    // }
    
    //for now this should suffice
    const int xdiv = 16;
    // const int ydiv = 1;
    // const int zdiv = 1;
    const MMeshType meshType = MMeshType::ETetrahedral;
    const TPZManVector<int,3> nDivs = {xdiv,xdiv,xdiv};

    TPZGeoMesh *gmesh = CreateGeoMesh(meshType,nDivs,EDomain,ESurfaces);
    // TPZGeoMesh *gmesh = CreateGeoMeshTetra(meshType,nDivs,EDomain,ESurfaces);


    TPZKernelHdivUtils<STATE> util;

    // ..................................HDiv..................................

    {
        TPZCompMesh * cmeshflux = 0;
        TPZCompMesh * cmeshpressure = 0;    

        std::set<int> matIdVecHdiv={EDomain,ESurfaces};
        std::set<int> matIdNeumannHdiv;
        std::set<int> matBCHybrid = {};
        std::set<int> matBC = {ESurfaces};
        
        TPZApproxSpaceKernelHdiv<STATE> createSpace(gmesh,
                                                    TPZApproxSpaceKernelHdiv<STATE>::ENone,        //Hybridization
                                                    HDivFamily::EHDivConstant); // Shape Type

    //     //Setting material ids
        createSpace.fConfig.fDomain = EDomain;
        createSpace.SetMaterialIds(EWrap,EPressureHyb,EIntface,EPont,matBCHybrid,matBC);
        createSpace.Initialize();
        // createSpace.SetEdgeRemove(EEdgeRemove);
        // util.PrintGeoMesh(gmesh);

        //Flux mesh
        cmeshflux = FluxCMesh(dim,pOrder,matIdVecHdiv,gmesh);     

        //Pressure mesh
        cmeshpressure = PressureCMesh(dim,0,matIdVecHdiv,gmesh);

        //Multiphysics mesh
        TPZManVector< TPZCompMesh *, 2> meshvector(2);
        meshvector[0] = cmeshflux;
        meshvector[1] = cmeshpressure;
        TPZMultiphysicsCompMesh * cmesh = MultiphysicCMesh(dim,pOrder,matIdVecHdiv,meshvector,gmesh);
        std::string multiphysicsFile = "MultiPhysicsMeshNew1";
        util.PrintCompMesh(cmesh,multiphysicsFile);

        // util.PrintCMeshConnects(cmesh);
        std::cout << "Equations without condense = " << cmesh->NEquations() - gmesh->NNodes()+1 << std::endl;
        // Group and condense the elements
        // createSpace.Condense(cmesh);
        TPZCompMeshTools::CreatedCondensedElements(cmesh,true,false);

        std::string multiphysicsFile1 = "MultiPhysicsMeshNew2";
        util.PrintCompMesh(cmesh,multiphysicsFile1);

        //Solve Multiphysics
        TPZLinearAnalysis an(cmesh,true);
        bool domainhyb=false;
                
        util.SolveProblemDirect(an,cmesh,false,domainhyb);
        std::cout << "Number of equations = " << cmesh->NEquations() << std::endl;
        
        // PrintEdges(gmesh,util.GetRemoveEdges(),util.GetVertexData(),util.GetEdgeData());

        //Print results
        an.SetExact(exactSol);
        util.PrintResultsMultiphysics(meshvector,an,cmesh);
        an.SetExact(exactSol);
        std::ofstream anPostProcessFileHdiv("postprocessHdiv.txt");
        util.ComputeError(an,anPostProcessFileHdiv);
    }




    // //............................Div Free Bubbles............................
    // {
    //     TPZCompMesh * cmeshflux = 0;
    //     TPZCompMesh * cmeshpressure = 0;    

    //     //Insert here the BC material id's to be hybridized 
    //     std::set<int> matBCHybrid={};
    //     //Insert here the type of all boundary conditions
    //     std::set<int> matIDNeumann{};
    //     std::set<int> matIDDirichlet{ESurfaces};
    //     /// All bc's mat ID's
    //     std::set<int> matBC;
    //     std::set_union(matIDNeumann.begin(),matIDNeumann.end(),matIDDirichlet.begin(),matIDDirichlet.end(),std::inserter(matBC, matBC.begin()));

    //     /// Creates the approximation space - Set the type of domain hybridization
    //     TPZApproxSpaceKernelHdiv<STATE> createSpace(gmesh,
    //                                                 TPZApproxSpaceKernelHdiv<STATE>::ENone,        //Hybridization
    //                                                 TPZApproxSpaceKernelHdiv<STATE>::EHDivKernel); // Shape Type

    //     //Setting material ids
    //     createSpace.fConfig.fDomain = EDomain;
    //     createSpace.SetMaterialIds(EWrap,EPressureHyb,EIntface,EPont,matBCHybrid,matBC);
    //     // createSpace.SetEdgeRemove(EEdgeRemove);
    //     createSpace.SetPOrder(pOrder);
    //     createSpace.Initialize();
    //     // util.PrintGeoMesh(gmesh);

    //     //Flux mesh
    //     TPZCompMesh * cmeshfluxNew = createSpace.CreateFluxCMesh();
    //     // std::cout << "FLUX \n";
    //     // util.PrintCMeshConnects(cmeshfluxNew);
    //     // std::string fluxFile = "FluxCMesh";
    //     // util.PrintCompMesh(cmeshfluxNew,fluxFile);
        
    //     // auto con = cmeshfluxNew->ConnectVec();
    //     // for (auto con :cmeshfluxNew->ConnectVec())
    //     // {
    //     //     std::cout << "con " << con << std::endl;
    //     // }
        


    //     //Pressure mesh
    //     TPZCompMesh * cmeshpressureNew = createSpace.CreatePressureCMesh();
    //     // std::cout << "PRESSURE \n";
    //     // util.PrintCMeshConnects(cmeshpressureNew);
    //     // std::string pressureFile = "PressureCMesh";
    //     // util.PrintCompMesh(cmeshpressureNew,pressureFile);        

    //     //Multiphysics mesh
    //     TPZManVector< TPZCompMesh *, 2> meshvectorNew(2);
    //     meshvectorNew[0] = cmeshfluxNew;
    //     meshvectorNew[1] = cmeshpressureNew;      
    //     auto * cmeshNew = createSpace.CreateMultiphysicsCMesh(meshvectorNew,exactSol,matIDNeumann,matIDDirichlet);
    //     // std::cout << "MULTIPHYSICS \n";
    //     // util.PrintCMeshConnects(cmeshNew);
    //     // Group and condense the elements
    //     // createSpace.Condense(cmeshNew);
    //     // std::string multiphysicsFile = "MultiPhysicsMeshNew";
    //     // util.PrintCompMesh(cmeshNew,multiphysicsFile);

    //     // Solve the problem
    //     TPZLinearAnalysis anNew(cmeshNew,false);
    //     // for (int i = 0; i < anNew.Solution().Rows(); i++)
    //     // {
    //         // TPZVec<int64_t> eqs = {0};
    //         // TPZManVector<std::string,10> scalnames(0), vecnames(1);
    //         // vecnames[0]= "State";

    //         // constexpr int resolution{2};
    //         // std::string plotfile = "ShapeRUIM.vtk";
    //         // anNew.DefineGraphMesh(3,scalnames,vecnames,plotfile);
    //         // anNew.ShowShape("ShapeRUIM",eqs);
    //     // }
        
        
        
    //     std::cout << "Number of equations = " << anNew.Mesh()->NEquations() << std::endl;        

    //     createSpace.Solve(anNew, cmeshNew, true, true); 
    //     std::cout << "Number of equations = " << anNew.Mesh()->NEquations() << std::endl;        
    //     anNew.SetExact(exactSol,solOrder);
    //     //Print results
        
    //     // // anNew.Solution().Zero();
    //     // const int rows = anNew.Solution().Rows();
    //     // const int cols = anNew.Solution().Cols();
    //     // TPZFMatrix<STATE> solNew(rows,cols,0.);
    //     // solNew(11,0)=1.;
        
        
    //     // anNew.Solution()=solNew;
    //     // anNew.LoadSolution();
    //     // std::ofstream outSol("Solution.txt");
    //     // anNew.Solution().Print("Sol",outSol,EMathematicaInput);

    //     // {
    //     //     TPZLinearAnalysis anFlux(cmeshfluxNew,false);
    //     //     for (int i = 0; i < rows; i++)
    //     //     {   
    //     //         solNew.Zero();
    //     //         solNew(i,0)=1.;
    //     //         anFlux.Solution()=solNew;
    //     //         anFlux.LoadSolution();
    //     //         // cmeshfluxNew->Solution().Print("Sol");
                
    //     //         TPZStack<std::string> scalnames;
    //     //         scalnames.Push("Curl");
    //     //         TPZStack<std::string> vecnames;
    //     //         vecnames.Push("Solution");
    //     //         const std::string plotfile = "shapeContornoNovo.vtk";
    //     //         std::set<int> matids = {ESurfaces};
                
    //     //         anFlux.DefineGraphMesh(2, matids, scalnames, vecnames, plotfile);
    //     //         anFlux.PostProcess(2,2);


    //     //         anNew.Solution()=solNew;
    //     //         anNew.LoadSolution();

    //     //         // TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
    //     //         TPZManVector<std::string,10> scalnames2(0), vecnames2(1);
    //     //         vecnames2[0]= "State";
    
    //     //         constexpr int resolution{2};
    //     //         std::string plotfile2 = "ShapeRUIM.vtk";
    //     //         anNew.DefineGraphMesh(3,scalnames2,vecnames2,plotfile2);
    //     //         anNew.PostProcess(resolution,3);
    //     //     }
            
           
    //     // }

    //     util.PrintResultsMultiphysics(meshvectorNew,anNew,cmeshNew);

    //     anNew.SetExact(exactSol,solOrder);
        

    //     std::ofstream out4("mesh_MDFB.txt");
    //     anNew.Print("nothing",out4);
    //     std::ofstream anPostProcessFileMDFB("postprocessMDFB.txt");
        
    //     util.ComputeError(anNew,anPostProcessFileMDFB);
    // }

  
    return 0;
}



//Flux computational mesh
TPZCompMesh *FluxCMesh(int dim, int pOrder,std::set<int> &matIdVec, TPZGeoMesh *gmesh) 
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh,false);

    for (std::set<int>::iterator it=matIdVec.begin(); it!=matIdVec.end(); ++it)
    {
        TPZNullMaterial<> *mat = new TPZNullMaterial<>(*it);
        mat->SetDimension(dim);
        mat->SetBigNumber(1.e10);
        cmesh->InsertMaterialObject(mat);
    }
    
    cmesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivConstant);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    
    // cmesh->ApproxSpace().SetHCurlFamily(HCurlFamily::EHCurlNoGrads);
    // cmesh->ApproxSpace().SetAllCreateFunctionsHCurl(dim);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->AutoBuild();

    // Print flux mesh
    std::ofstream myfile("FluxMesh.txt");
    cmesh->Print(myfile);

  return cmesh;
}

// Pressure computational mesh
TPZCompMesh *PressureCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZGeoMesh *gmesh)
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh,false);
    if(matIdVec.size() == 0) return cmesh;

    TPZNullMaterial<> *mat = new TPZNullMaterial<>(EDomain);
    mat->SetDimension(dim);
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
    // distincao de ordem zero
    if(pOrder == 0)
    {
        // os elementos H1 nao tem opcao de funcao constante
        cmesh->SetAllCreateFunctionsDiscontinuous();
    }
    else
    {
        // Disconnected = true faz com que os espaços sao descontinuos
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
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
        if(!celdisc) continue;
        celdisc->SetConstC(1.);
        celdisc->SetTrueUseQsiEta();
        // espera-se elemento de pressao apenas para o contorno
        if(celdisc && celdisc->Reference()->Dimension() == cmesh->Dimension())
        {
            // DebugStop();
        }
    }
    // // Print pressure mesh
    // std::ofstream myfile("PressureMesh.txt");
    // cmesh->Print(myfile);

    return cmesh;
}

// Multiphysics computational mesh
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, std::set<int> &matIdVec, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    TPZMixedDarcyFlow* mat = new TPZMixedDarcyFlow(EDomain, dim);
    mat->SetConstantPermeability(1.);
    
    cmesh->InsertMaterialObject(mat);
    mat -> SetBigNumber(1.e10);
        
    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,0.);
    auto * BCond0 = mat->CreateBC(mat, ESurfaces, 0, val1, val2);
    BCond0->SetForcingFunctionBC(exactSol);
    cmesh->InsertMaterialObject(BCond0);

    TPZManVector<int> active(2,1);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    cmesh->BuildMultiphysicsSpace(active, meshvector);

    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);

    // Prints Multiphysics mesh
    std::ofstream myfile("MultiPhysicsMesh.txt");
    cmesh->Print(myfile);

    // Prints computational mesh properties
    // std::stringstream vtk_name;
    // vtk_name    << "MultiPhysics" << ".vtk";
    // std::ofstream vtkfile(vtk_name.str().c_str());
    // TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

    return cmesh;
}


TPZGeoMesh *CreateGeoMesh(const MMeshType meshType, const TPZVec<int> &nDivs,
              const int volId, const int bcId)
{
  constexpr int dim{3};
  const TPZManVector<REAL,3> minX = {0,0,0};
  const TPZManVector<REAL,3> maxX = {1,1,1};
  constexpr bool createBoundEls{true};

  //all bcs share the same id
  const TPZManVector<int,7> matIds(7,bcId);
  matIds[0] = volId;
  
  TPZGeoMesh *gmesh = 
    TPZGeoMeshTools::CreateGeoMeshOnGrid(dim, minX, maxX,
            matIds, nDivs, meshType,createBoundEls);

  
  return gmesh;
}

TPZGeoMesh *CreateGeoMeshTetra(const MMeshType meshType, const TPZVec<int> &nDivs,
              const int volId, const int bcId)
{
    constexpr int dim{3};
    const TPZManVector<REAL,3> minX = {0,0,0};
    const TPZManVector<REAL,3> maxX = {1,1,1};
    constexpr bool createBoundEls{true};

    //all bcs share the same id
    const TPZManVector<int,7> matIds(7,bcId);
    matIds[0] = volId;
    
    TPZGeoMesh *gmesh = 
    TPZGeoMeshTools::CreateGeoMeshSingleEl(meshType,volId,createBoundEls);

    for (auto gel : gmesh->ElementVec())
    {
        if (gel->Dimension() == 3) continue;

        gel->SetMaterialId(bcId);

    }
    

  return gmesh;
}


//Flux computational mesh
TPZCompMesh *FluxCMeshCurl(int dim, int pOrder, TPZGeoMesh *gmesh) 
{
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

    // cmesh->SetDefaultOrder(pOrder);
    // cmesh->SetDimModel(dim);

    // //insert volumetric material
    // auto volMat = new TPZMatCurlDotCurl(EDomain);
    // cmesh->InsertMaterialObject(volMat);
    // //insert boundary material
    // const int bcType = 0;//dirichlet
    // TPZFNMatrix<1, STATE> val1(1, 1, 1);
    // TPZManVector<STATE,1> val2(1, 0.);
    // auto bcMat = volMat->CreateBC(volMat, ESurfaces, bcType, val1, val2);
    // cmesh->InsertMaterialObject(bcMat);
    
    // //Creates computational elements
    // int64_t nel = gmesh->NElements();
    // for (int64_t el = 0; el < nel; el++) {
    //     TPZGeoEl *gel = gmesh->Element(el);
    //     if(!gel) DebugStop();
    //     gel->ResetReference();   
    //     const MElementType type = gel->Type();
    //     const auto matid = gel->MaterialId();
    //     int64_t index;
    //     switch(type){
    //     case ETriangle:
    //     new TPZCompElHCurlNoGrads<pzshape::TPZShapeTriang>(*cmesh,gel,index);
    //     break;
    //     case ETetraedro:
    //     new TPZCompElHCurlNoGrads<pzshape::TPZShapeTetra>(*cmesh,gel,index);
    //     break;
    //     default:
    //     const auto elName =  MElementType_Name(type);
    //     // CAPTURE(elName);
    //     // CHECK(false);
    //     PZError<<__PRETTY_FUNCTION__
    //             <<"\n type not yet supported. Aborting..."<<std::endl;
    //     DebugStop();
    //     }
    // }
    // cmesh->SetAllCreateFunctionsHCurl();
    // cmesh->AutoBuild();
    // cmesh->CleanUpUnconnectedNodes();

    return cmesh;
}


void PrintEdges(TPZGeoMesh * gmesh, std::set<int64_t> &rem_edges, std::map<int64_t, TPZHCurlEquationFilter<STATE>::VertexFilter> &vertexData, std::map<int64_t, TPZHCurlEquationFilter<STATE>::EdgeFilter> &edgeData){

    int idRem=50;
    int idEdg=40;
    int node0=idRem;
    idRem++;
    int nnodes = vertexData.size();

    for (int i = 0; i < edgeData.size(); i++)
    {
        edgeData[i].status = TPZHCurlEquationFilter<STATE>::EdgeStatusType::EFreeEdge;
    }
    

    for (auto gel:gmesh->ElementVec())
    {
        if (gel->Dimension() < 3) continue;
        int sideStart = 0;
        int sideEnd = 0;
        if (gel->Type() == ETetraedro){
            sideStart = 4;
            sideEnd = 10;
        } else if (gel->Type() == ECube){
            sideStart = 8;
            sideEnd = 20;
        } else {
            DebugStop();
        }
        //For tetrahedra only, loop over the surface sides
        for (int side = sideStart; side < sideEnd; side++){
            TPZGeoElSide gelside(gel,side);
            auto con = gel->Reference()->ConnectIndex(side-sideStart);

            if (edgeData[con].status != TPZHCurlEquationFilter<STATE>::EdgeStatusType::EFreeEdge) continue;

            edgeData[con].status = TPZHCurlEquationFilter<STATE>::EdgeStatusType::ERemovedEdge;
            int Nmatid = 0;
            if (rem_edges.find(con) != rem_edges.end()){
                Nmatid = idRem++;
                for (int iNode = 0; iNode < nnodes; iNode++)
                {
                    if (con == vertexData[iNode].removed_edge){
                        int64_t index;
                        TPZVec<int64_t> corner{iNode};
                        if (iNode == 5){
                            gmesh->CreateGeoElement(EPoint,corner,node0,index);
                        } else {
                            gmesh->CreateGeoElement(EPoint,corner,Nmatid,index);
                        }
                    }
                }
            }else{
                Nmatid = idEdg;
            }
            TPZGeoElBC gelbcWrap(gelside, Nmatid);
            
            
        }
    }
    
    //Prints computational mesh with edges
    std::string vtk_name = "FluxEdges.vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, vtkfile, true);
    

}