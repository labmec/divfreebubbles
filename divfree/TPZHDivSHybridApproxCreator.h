#ifndef TPZHDIVSHYBRIDAPPROXCREATOR_H
#define TPZHDIVSHYBRIDAPPROXCREATOR_H

#include "TPZHDivApproxCreator.h"

class TPZGeoMesh;

/// @class TPZHDivSHybridApproxCreator
/// @brief Derived HDiv approximation creator for semi-hybridized spaces
class TPZHDivSHybridApproxCreator : public TPZHDivApproxCreator {
public:
  /// Default constructor
  TPZHDivSHybridApproxCreator() = default;

  /// Constructor with a geometric mesh
  explicit TPZHDivSHybridApproxCreator(TPZGeoMesh *gmesh);

  /// Destructor
  ~TPZHDivSHybridApproxCreator() = default;

  /// Add a hybrid geometric configuration
  virtual void AddHybridizationGeoElements() override;

  /// Creates and HDiv approximation space/cmesh
  virtual TPZCompMesh *CreateHDivSpace() override;

  /// Insert interface periferal material objects related to geometric objects created during hybridization
  virtual void InsertInterfaceMaterialObjects(TPZMultiphysicsCompMesh *mphys);

  /// Create interface elements on hybridizes spaces
  /// @param mphys multiphysics compmesh
  virtual void AddInterfaceComputationalElements(TPZMultiphysicsCompMesh *mphys) override;

  /// Groups the elements in data structure to be condensed
  /// @param mcmesh multiphysics compmesh with elements to be condensed
  virtual void GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh) override;

  void SetShouldHybridizeLowOrderFluxes(bool value) { fShouldHybridizeLowOrderFluxes = value; }
  const bool &ShouldHybridizeLowOrderFluxes() const { return fShouldHybridizeLowOrderFluxes; }

  /// @brief hybridize the low order fluxes
  void HybridizeLowOrderFluxes(TPZMultiphysicsCompMesh &mfmesh);

private:
  bool fShouldHybridizeLowOrderFluxes = true;
};

#endif
