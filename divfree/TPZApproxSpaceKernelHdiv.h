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

private:
    /// the type of space this object will generate
    MSpaceType fSpaceType = ENone;

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

public:
    /// default constructor
    TPZApproxSpaceKernelHdiv(TPZGeoMesh *gmesh, MSpaceType spacetype = ENone);
    
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

    void Initialize();

    void Solve(TPZLinearAnalysis &an, TPZMultiphysicsCompMesh * cmesh, bool direct);

    void Condense(TPZMultiphysicsCompMesh * cmesh)
    {
        hybridizer.GroupAndCondenseElements(cmesh,fConfig.fBCHybridMatId);
    }

    TPZCompMesh * CreateFluxCMesh();

    TPZCompMesh * CreatePressureCMesh();

    TPZMultiphysicsCompMesh * CreateMultiphysicsCMesh(TPZVec<TPZCompMesh *> &meshvector, ForcingFunctionBCType<TVar> exactSol, std::set<int> &BCNeumann, std::set<int> &BCDirichlet);
    
    /// All parameters needed for creating a hybrid H1 space
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

    void SetPeriferalMaterialIds(int Wrap, int Lagrange, int Interface, int Point, std::set<int > &matBChybrid, std::set<int > &matBC)
    {
        fConfig.fWrap = Wrap;
        fConfig.fLagrange = Lagrange;
        fConfig.fInterface = Interface;
        fConfig.fPoint = Point;
        fConfig.fBCHybridMatId = matBChybrid;
        fConfig.fBCMatId = matBC;
    }

};

#endif //TPZApproxSpaceKernelHdiv