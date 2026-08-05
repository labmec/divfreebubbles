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

  /// @brief boolean indicating if the low order fluxes have already been hybridized
  bool fLowOrderFluxHybridized = false;

};

#endif
