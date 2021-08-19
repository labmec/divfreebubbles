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
class TPZCompMesh;
class TPZGeoMesh;
class TPZMultiphysicsCompMesh;
class TPZCompEl;
class TPZGeoElSide;


class TPZApproxSpaceKernelHdiv
{
public:
    /// types of spaces this class can create
    enum MSpaceType {ENormal, EDomainHybrid, EBCHybrid, EAllHybrid, EDomainSemiHybrid};

private:
    /// the type of space this object will generate
    MSpaceType fSpaceType = ENormal;

    /// the materialids which will be used to create the atomic meshes
    std::set<int> fMaterialIds;
    
    /// the boundary condition material ids
    std::set<int> fBCMaterialIds;

    /// default internal order for the H1 elements
    int fDefaultPOrder = 3;
    
    /// default Lagrange multiplier order
    int fDefaultLagrangeOrder = 1;
    
    /// the dimension of the geometric elements that will be used to generate computational elements
    int fDimension = -1;
    
    /// the geometric mesh which will generate the computational mesh
    TPZGeoMesh *fGeoMesh = 0;

public:
    /// default constructor
    TPZApproxSpaceKernelHdiv(TPZGeoMesh *gmesh, MSpaceType spacetype = ENormal);
    
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

    TPZCompMesh * CreateFluxCMesh(std::set<int> &matIdVec);

    TPZCompMesh * CreatePressureCMesh(std::set<int> &matIdVec, std::set<int> &matBC);

    void CreateMultiphysicsCMesh();
    
    /// All parameters needed for creating a hybrid H1 space
    struct TConfig
    {
        /// Domain material ID
        int fDomain = 1;
        
        /// Wrap material ID
        int fWrap = -1;

        /// Point material ID
        int fPoint = -1;

        /// Boundary conditions material ID's
        std::set<int> fBCMatId;

        /// default constructor
        TConfig(){}
        
        /// copy constructor
        TConfig(const TConfig &copy);
        
        /// copy operator
        TConfig &operator=(const TConfig &copy);
    };

    /// object which contains the relevant information for create a hybrid H1 mesh
    TConfig fConfig;

};

#endif //TPZApproxSpaceKernelHdiv