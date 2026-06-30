//
// Created by GitHub Copilot based on NeoPZ patterns.
//

#ifndef PZ_TPZH1HYBRIDAPPROXCREATOR_H
#define PZ_TPZH1HYBRIDAPPROXCREATOR_H

#include "TPZH1ApproxCreator.h"
#include "pzinterpolationspace.h"
#include "pzfmatrix.h"

#include <map>

class TPZMultiphysicsCompMesh;
class TPZMultiphysicsElement;

/// @class TPZH1HybridApproxCreator
/// @brief Derived H1 approximation creator for hybrid H1 approximation spaces
class TPZH1HybridApproxCreator : public TPZH1ApproxCreator {
public:
    TPZH1HybridApproxCreator() = default;
    TPZH1HybridApproxCreator(TPZGeoMesh * gmesh);
    ~TPZH1HybridApproxCreator() = default;

    /// Build the approximation space for a hybrid H1 problem
    TPZMultiphysicsCompMesh *CreateHybridApproximationSpace();

    typedef std::map<int64_t, int64_t> CtoMFCel;
    /// @brief Compute the constraints to orthogonalize the restraints
    // this function will compute restraints for all boundary flux connects
    static void ComputeOrthogonalizingRestraints(TPZMultiphysicsCompMesh &mfmesh, CtoMFCel &geltogel, TPZApproxCreator::HybridizationData hybridData);

    /// Put elements in element groups
    virtual void GroupElements(TPZMultiphysicsCompMesh *mcmesh) override;
    
    /// Create condensed elements around group elements
    /// this method will adjust the connect count of DOF's that should not be condensed
    virtual void CondenseElements(TPZMultiphysicsCompMesh *mcmesh) override;


  /// @brief hybridize the low order fluxes
  void HybridizeLowOrderFluxes(TPZMultiphysicsCompMesh &mfmesh, CtoMFCel &geltogel);

protected:
    /// Add the hybridization tree entry for an interface
    void AddHybrid(std::map<int64_t,TPZHybrid> &tree, int64_t interface_id, TPZHybrid &hybrid);

    /// Add geometric elements to represent squared hybridization
    virtual void AddHybridSquareGeoElements() override;

    /// Compute projection directions for a geometric element
    static void ProjectionDirections(TPZGeoEl *gel, TPZFMatrix<REAL> &projdir);

    /// Compute projection values for a set of relative coordinates
    static void ProjectionValues(TPZVec<REAL> &xrelative, TPZFMatrix<REAL> &projdir, int ncorner, TPZFMatrix<REAL> &funcval);

    /// Compute the projection matrix for an interpolation space
    static void ComputeProjectionMatrix(TPZInterpolationSpace *intel, TPZFMatrix<REAL> &projection);

  /// @brief compute the restraints of a connect for a geometric element
  static void RestraintConnect(TPZInterpolationSpace *intel, TPZMultiphysicsElement *mfcel, int64_t newind1, int64_t newind2);

  /// @brief compute an equivalent B matrix
  static void ComputeBMatrix(TPZInterpolationSpace *intel, TPZFMatrix<STATE> &B);

  /// @brief  compute orthogonalizing restraint
  static void ComputeOrthogonalizingR(TPZFMatrix<REAL> &B, TPZFMatrix<REAL> &Restraint);

  /// @brief boolean indicating if the low order fluxes have already been hybridized
  bool fLowOrderFluxHybridized = false;

};

#endif
