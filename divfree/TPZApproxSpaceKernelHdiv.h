//
//  TPZApprxSpaceKernelHdiv.h
//  DivFreeBubbles project
//
//  Created by Jeferson Fernandes on 18/08/21.
//
//  Based on TPZCreateMultiphysicsSpace from ErrorEstimate project
//

#ifndef TPZApproxSpaceKernelHdiv_hpp
#define TPZApproxSpaceKernelHdiv_hpp

#include <stdio.h>
#include <set>
#include "pzmanvector.h"
#include "TPZMatTypes.h"
#include "TPZKernelHdivHybridizer.h"
#include "TPZKernelHdivUtils.h"
#include "TPZEnumApproxFamily.h"

class TPZCompMesh;
class TPZGeoMesh;
class TPZMultiphysicsCompMesh;
class TPZCompEl;
class TPZGeoElSide;

template <class TVar>
class TPZApproxSpaceKernelHdiv
{
public:
    /// types of spaces this class can create
    enum MSpaceType {ENone, EFullHybrid, ESemiHybrid};
    // enum MShapeType {EHDivKernel, EHDivConstant, EHCurlNoGrads};

private:
    /// the type of space this object will generate
    MSpaceType fSpaceType = ENone;

    /// the type of space this object will generate
    HDivFamily fShapeType = HDivFamily::EHDivKernel;

    /// default internal order for the H1 elements
    int fDefaultPOrder = 3;
    
    /// default Lagrange multiplier order
    int fDefaultLagrangeOrder = 1;
    
    /// the dimension of the geometric elements that will be used to generate computational elements
    int fDimension = -1;
    
    /// the geometric mesh which will generate the computational mesh
    TPZGeoMesh *fGeoMesh = 0;

    /// The hybridizer
    TPZKernelHdivHybridizer hybridizer;

    /// The util
    TPZKernelHdivUtils<TVar> * util;

    REAL fBigNumber = 1.e10;

public:
    /// default constructor
    TPZApproxSpaceKernelHdiv(TPZGeoMesh *gmesh, MSpaceType spacetype = ENone, HDivFamily shapetype = HDivFamily::EHDivKernel);
    
    /// copy constructor
    TPZApproxSpaceKernelHdiv(const TPZApproxSpaceKernelHdiv &copy);
    
    /// = operator
    TPZApproxSpaceKernelHdiv &operator=(const TPZApproxSpaceKernelHdiv &copy);

    //Sets polynomial order
    void SetPOrder(int order){
        fDefaultPOrder = order;
    }
    
    //Sets Lagrange multipliers polynomial order
    void SetLagrangeOrder(int order){
        fDefaultLagrangeOrder = order;
    }

    /**
     * @brief Sets the material id's to the hybridizer and creates the 
     * additional geometric elements needed to compute the solution
     * 
     */
    void Initialize();

    /**
     * @brief Solves the problem
     * 
     * @param an 
     * @param cmesh 
     * @param direct = true if direct solver, false if iterative solver
     */
    void Solve(TPZLinearAnalysis &an, TPZCompMesh * cmesh, bool direct, bool filterEquation);

    /**
     * @brief Performs the static condensation
     * 
     * @param cmesh 
     */
    void Condense(TPZMultiphysicsCompMesh * cmesh)
    {
        hybridizer.GroupAndCondenseElements(cmesh,fConfig.fBCHybridMatId);
    }

    /**
     * @brief Creates the flux KernelHdiv computational mesh
     * 
     * @return TPZCompMesh* 
     */
    TPZCompMesh * CreateFluxCMesh();

    /**
     * @brief Create the Pressure KernelHdiv computational mesh 
     * (only has elements when some hybridization is implemented)
     * 
     * @return TPZCompMesh* 
     */
    TPZCompMesh * CreatePressureCMesh();
    TPZCompMesh * CreatePressureCMeshHybridizedHDivConstant();

    /**
     * @brief Create the Multiphysics mesh (always, even when the pressure
     * computational mesh is empty)
     * 
     * @param meshvector 
     * @param exactSol 
     * @param BCNeumann 
     * @param BCDirichlet 
     * @return TPZMultiphysicsCompMesh* 
     */
    TPZMultiphysicsCompMesh * CreateMultiphysicsCMesh(TPZVec<TPZCompMesh *> &meshvector, ForcingFunctionBCType<TVar> exactSol, std::set<int> &BCNeumann, std::set<int> &BCDirichlet);
    
    /// Parameters needed for creating a hybrid KernelHdiv space
    struct TConfig
    {
        /// Domain material ID
        int fDomain = 1;
        
        /// Wrap material ID
        int fWrap = -1;

        /// Point material ID
        int fPoint = -1;

        /// Interface material ID
        int fInterface = -1;

        /// Lagrange Multiplier Material ID
        int fLagrange = -1;

        /// Edge to be removed Material ID
        int fEdgeRemove = -1;

        /// Hybridized Boundary conditions material ID's
        std::set<int> fBCHybridMatId = {};
        
        /// All Boundary conditions material ID's
        std::set<int> fBCMatId = {};

        /// default constructor
        TConfig(){}
        
        /// copy constructor
        TConfig(const TConfig &copy);
        
        /// copy operator
        TConfig &operator=(const TConfig &copy);

        
    };

    /// object which contains the relevant information for create a hybrid H1 mesh
    TConfig fConfig;

    /**
     * @brief Sets the additional material ids needed for the problem solution
     * 
     * @param Wrap 
     * @param Lagrange 
     * @param Interface 
     * @param Point 
     * @param matBChybrid 
     * @param matBC 
     */
    void SetPeriferalMaterialIds(int Wrap, int Lagrange, int Interface, int Point, std::set<int > &matBChybrid, std::set<int > &matBC)
    {
        fConfig.fWrap = Wrap;
        fConfig.fLagrange = Lagrange;
        fConfig.fInterface = Interface;
        fConfig.fPoint = Point;
        fConfig.fBCHybridMatId = matBChybrid;
        fConfig.fBCMatId = matBC;
    }

    void SetEdgeRemove(int edge){
        fConfig.fEdgeRemove = edge;
    }

    void CreateOrientedBoundaryElements();
    
    void CreateFluxHybridezedHDivKernel(TPZCompMesh *cmesh);
    void CreateFluxHybridezedHDivConstant(TPZCompMesh *cmesh);

};

#endif //TPZApproxSpaceKernelHdiv