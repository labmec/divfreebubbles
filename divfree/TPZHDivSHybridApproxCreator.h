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
    ~TPZHDivSHybridApproxCreator() override = default;

    /// Add a hybrid geometric configuration
    virtual void AddHybridizationGeoElements() override;

        /// Creates and HDiv approximation space/cmesh
    virtual TPZCompMesh * CreateHDivSpace() override;

    /// Insert interface periferal material objects related to geometric objects created during hybridization
    virtual void InsertInterfaceMaterialObjects(TPZMultiphysicsCompMesh *mphys) override;

    /// Create interface elements on hybridizes spaces
    /// @param mphys multiphysics compmesh
    virtual void AddInterfaceComputationalElements(TPZMultiphysicsCompMesh *mphys) override;
  

};

#endif
