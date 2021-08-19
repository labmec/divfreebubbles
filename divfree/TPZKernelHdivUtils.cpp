//
// Created by Jeferson Fernandes on 11/08/21.
//

#include "TPZKernelHdivUtils.h"

#include "pzcmesh.h"
#include "pzgmesh.h"
#include <TPZVTKGeoMesh.h>
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver


void TPZKernelHdivUtils::PrintCMeshConnects(TPZCompMesh *cmesh){
    
    for (int i = 0; i < cmesh->NElements(); i++)
    {
        TPZCompEl *cel = cmesh->Element(i);
        cel->LoadElementReference();
        int matid = cel->Reference()->MaterialId();
        auto nconnects = cel->NConnects();
        std::cout << "Element = " << i << ", dim= " << cel->Dimension() << ",mat = " << cel->Reference()->MaterialId() << ", nconnects= " << nconnects << ": ";
        for (int j = 0; j < nconnects; j++)
        {
            std::cout << cel->ConnectIndex(j) << ", ";
        }
        std::cout << std::endl;
    
        // std::cout << cel->Connect() << std::endl;

        // // loop only over volumetric elements
        // if(matid != EDomain) continue;
        // if (cel->Reference()->Dimension() != dim) {
        //     DebugStop();
        // }

        // int nsides = cel->Reference()->NSides();
        // int ncorner = cel->Reference()->NCornerNodes();
        // for (int side = 0; side < nsides; side++) {
        //     if(cel->Reference()->SideDimension(side) != dim-1) continue;
        //     TPZGeoElSide gelside(cel->Reference(),side);
        //     TPZGeoElSide neighbour = gelside.Neighbour();
            
        //     std::cout << "Element = " << i << ", side = " << side  
        //             << ", NEl = " << neighbour.Element()->Index()
        //             << ", Nmatid = " << neighbour.Element()->MaterialId()
        //             << ", NNEl = " << neighbour.Neighbour().Element()->Index()
        //             << ", NNmatid = " << neighbour.Neighbour().Element() -> MaterialId() << std::endl;
        //     std::cout << "Neigh connect : " ;
        //     nconnects = neighbour.Element()->Reference()->NConnects();
        //     for (int j = 0; j < nconnects; j++)
        //     {
        //         std::cout << neighbour.Element()->Reference()->ConnectIndex(j) << ", ";
        //     }
        //     std::cout << std::endl;
        // }
    }
}

void TPZKernelHdivUtils::PrintGeoMesh(TPZGeoMesh *gmesh){
    
    for (int i = 0; i < gmesh->NElements(); i++)
    {

        auto *gel = gmesh->Element(i);
        int matid = gel->MaterialId();
        auto nsides = gel->NSides();
        // auto nconnects = gel->Reference()->NConnects();
        std::cout << "ELGeometric = " << i << ", dim= " << gel->Dimension() << ",mat = " << gel->MaterialId() << std::endl;
  
        nsides = gel->NSides();
        int ncorner = gel->NCornerNodes();
        for (int side = 0; side < nsides; side++) {
            // if(gel->SideDimension(side) != 1) continue;
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neighbour = gelside.Neighbour();
            
            std::cout << "Element = " << i << ", side = " << side  
                    << ", NEL = " << neighbour.Element()->Index() 
                    << ", Nmatid = " << neighbour.Element()->MaterialId()
                    << ", NNEL = " << neighbour.Neighbour().Element()->Index() 
                    << ", NNmatid = " << neighbour.Neighbour().Element() -> MaterialId() << std::endl;
        }
    }
}

void TPZKernelHdivUtils::PrintCompMesh(TPZCompMesh *cmesh,std::string &file_name){

    // Print pressure mesh
    std::string txt = file_name + ".txt";
    std::ofstream myfile(txt);
    cmesh->Print(myfile);

    //Prints computational mesh properties
    std::string vtk_name = file_name + ".vtk";
    std::ofstream vtkfile(vtk_name.c_str());
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, vtkfile, true);

}

