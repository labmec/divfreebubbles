//
//  TPZHDivApproxSpaceCreator.h
//  DivFreeBubbles project
//
//  Created by Jeferson Fernandes on 18/08/21.
//
//

#ifndef TPZHDiv_ApproxSpace_Creator
#define TPZHDiv_ApproxSpace_Creator

#include <stdio.h>
#include <set>
#include "pzmanvector.h"
#include "TPZMatTypes.h"
#include "TPZKernelHdivHybridizer.h"
#include "TPZKernelHdivUtils.h"
#include "TPZEnumApproxFamily.h"
#include "TPZAnalyticSolution.h"

class TPZCompMesh;
class TPZGeoMesh;
class TPZMultiphysicsCompMesh;
class TPZCompEl;
class TPZGeoElSide;

template <class TVar>
class TPZHDivApproxSpaceCreator
{
public:
    /// types of spaces this class can create
    enum MSpaceType {ENone, EFullHybrid, ESemiHybrid, EDuplicatedConnects};

private:
    /// the type of space this object will generate
    MSpaceType fSpaceType = ENone;

    /// the type of space this object will generate
    HDivFamily fShapeType = HDivFamily::EHDivKernel;

    /// default internal order for the H1 elements
    int fDefaultPOrder = 1;
    
    /// default Lagrange multiplier order
    int fDefaultLagrangeOrder = 1;
    
    /// the dimension of the geometric elements that will be used to generate computational elements
    int fDimension = -1;

    bool mixedElasticity = false;
    
    /// the geometric mesh which will generate the computational mesh
    TPZGeoMesh *fGeoMesh = 0;

    /// The hybridizer
    TPZKernelHdivHybridizer hybridizer;

    /// The util
    TPZKernelHdivUtils<TVar> * util;

    REAL fBigNumber = 1.e10;

public:
    /// default constructor
    TPZHDivApproxSpaceCreator(TPZGeoMesh *gmesh, MSpaceType spacetype = ENone, HDivFamily shapetype = HDivFamily::EHDivKernel);
    
    /// copy constructor
    TPZHDivApproxSpaceCreator(const TPZHDivApproxSpaceCreator &copy);
    
    /// = operator
    TPZHDivApproxSpaceCreator &operator=(const TPZHDivApproxSpaceCreator &copy);

    //Sets polynomial order
    void SetPOrder(int order){
        fDefaultPOrder = order;
    }
    
    //Sets Lagrange multipliers polynomial order
    void SetLagrangeOrder(int order){
        fDefaultLagrangeOrder = order;
    }

    /**
     * @brief Initialize the data structure to provide an interface to create several approximation spaces using HDiv family of functionss
     * //In the case of duplicated connects a different approach is used. Instead of using new geometric (wrap) elements, the edge (in 2D) or face (in 3D) connects
        //are duplicated
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

    void CondenseDuplicatedConnects(TPZMultiphysicsCompMesh * cmesh)
    {
        hybridizer.GroupAndCondenseElementsDuplicatedConnects(cmesh,fConfig.fBCMatId);
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
    TPZMultiphysicsCompMesh * CreateMultiphysicsCMeshElasticity(TPZVec<TPZCompMesh *> &meshvector, TPZAnalyticSolution * gAnalytic, std::set<int> &BCNeumann, std::set<int> &BCDirichlet);
    
    HDivFamily GetHDivFamily() {return fShapeType;}

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
    void SetMaterialIds(int Wrap, int Lagrange, int Interface, int Point, std::set<int > &matBChybrid, std::set<int > &matBC)
    {
        fConfig.fWrap = Wrap;
        fConfig.fLagrange = Lagrange;
        fConfig.fInterface = Interface;
        fConfig.fPoint = Point;
        fConfig.fBCHybridMatId = matBChybrid;
        fConfig.fBCMatId = matBC;
    }

    void CreateOrientedBoundaryElements();
    
    void CreateFluxHybridezedHDivKernel(TPZCompMesh *cmesh);
    void CreateFluxHybridezedHDivConstant(TPZCompMesh *cmesh);

    void DuplicateInternalConnects(TPZCompMesh *cmesh);

    TPZCompMesh * CreateConstantCmesh(TPZGeoMesh *Gmesh, int lagLevel);

    void SetMixedElasticity(){mixedElasticity = true;}
    
    void ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder);
    TPZCompMesh *CreateRotationCmesh(TPZGeoMesh *gmesh, int pOrder, REAL elementdim);
};

#endif //TPZHDivApproxSpaceCreator