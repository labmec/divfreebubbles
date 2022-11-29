//
// Created by Jeferson Fernandes on 11/08/21.
//

#ifndef KERNELHDIV_UTILS_H
#define KERNELHDIV_UTILS_H

// TODO add doc

#include <iostream>
#include <map>
#include "pzstack.h"
#include "TPZLinearAnalysis.h"
#include <TPZVTKGeoMesh.h>
#include "TPZHCurlEquationFilter.h"
#include "TPZMultiphysicsCompMesh.h"
//#include "pzgeoelrefless.h"

template<class T>
class TPZVec;

class TPZCompMesh;
class TPZGeoMesh;
class TPZMaterial;
class TPZCompElSide;
class TPZInterpolatedElement;
class TPZGeoEl;

template<class T, int N>
class TPZStack;

template <class TVar>
class TPZKernelHdivUtils {

public:
    /**
     * @brief Prints the computational mesh information (specially the connects)
     * 
     * @param cmesh 
     */
    void PrintCMeshConnects(TPZCompMesh *cmesh);

    /**
     * @brief Prints geometric mesh information (specially neighbours)
     * 
     * @param geomesh 
     */
    void PrintGeoMesh(TPZGeoMesh *geomesh);

    /**
     * @brief Prints computational mesh information to a file
     * 
     * @param cmesh 
     * @param file_name 
     */
    void PrintCompMesh(TPZCompMesh *cmesh, std::string file_name);

    /**
     * @brief Solves an algebraic system by means of an iterative method
     * 
     * @param an 
     * @param cmesh 
     */
    void SolveProblemIterative(TPZLinearAnalysis &an, TPZCompMesh *cmesh);

    /**
     * @brief Solves an algebraic system by means of a direct method
     * 
     * @param an 
     * @param cmesh 
     */
    void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh, bool filterEquations, bool &domanHybridization); 
    void SolveProblemCholesky(TPZLinearAnalysis &an, TPZCompMesh *cmesh, bool filterEquations, bool &domanHybridization); 
    
    /**
     * @brief Prints the results of a multiphysics mesh to a .vtk file
     * 
     * @param meshvector 
     * @param an 
     * @param cmesh 
     */
    void PrintResultsMultiphysics(TPZVec<TPZCompMesh *> &meshvector, TPZLinearAnalysis &an, TPZMultiphysicsCompMesh *cmesh);

    /**
     * @brief Util to compute the solution error
     * 
     * @param an 
     * @param anPostProcessFile 
     */
    void ComputeError(TPZLinearAnalysis &an, std::ostream &anPostProcessFile);

    std::map<int64_t, TPZHCurlEquationFilter<STATE>::VertexFilter> &GetVertexData(){return vertexData;}
    std::map<int64_t, TPZHCurlEquationFilter<STATE>::EdgeFilter> &GetEdgeData(){return edgeData;}

    /**
     @brief Reads the test mesh from gmsh
    @param[in] file_name the .msh mesh file.
    */
    template<class tshape>
    TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

    /**
     @brief Creates a geometric mesh with elements of a given type on a unit square or cube (depending on the mesh dimension).
    @param[in] meshType element type to be created.
    @param[in] nDivs Number of divisions (rows of elements) in x, y and z.
    @param[in] volId Material identifier for the volumetric region.
    @param[in] bcId Material identifier for the boundary.
    */
    template<class tshape>
    TPZGeoMesh* CreateGeoMesh(TPZVec<int> &nDivs, int volId, int bcId, bool createBoundEls);

private:
    std::map<int64_t, TPZHCurlEquationFilter<STATE>::VertexFilter> vertexData;
    std::map<int64_t, TPZHCurlEquationFilter<STATE>::EdgeFilter> edgeData;
};
#endif