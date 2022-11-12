//
// Created by Gustavo Batistela on 5/13/21.
//

#ifndef TPZOrtotropicPermeability_H
#define TPZOrtotropicPermeability_H

#include <functional>
#include "pzreal.h"
#include "pzvec.h"
#include "pzfmatrix.h"

// Alias to improve readability of the permeability function type
using OrtoPermeabilityFunctionType = std::function<TPZManVector<STATE,3>(const TPZVec<REAL> &coord)>;

// Forward declaration of dummy BC interface class
class TPZOrtotropicPermeabilityBC;

/**
 * @brief  This class implements the interface with the methods required to
 * handle the permeability field of an isotropic material.
 */
class TPZOrtotropicPermeability : public TPZSavable {

public:
    using TInterfaceBC = TPZOrtotropicPermeabilityBC;

    /**
     * @brief Set a varying permeability field to the material
     * @param [in] perm_function function that describes the permeability field
     */
    void SetPermeabilityFunction(OrtoPermeabilityFunctionType &perm_function);

    /**
     * @brief Return the permeability value at a coordinate
     * @param [in] coord coordinate of interest
     */
    TPZManVector<REAL, 3> GetPermeability(const TPZVec<REAL> &coord);

    [[nodiscard]] int ClassId() const override;

    void Read(TPZStream &buf, void *context) override {};

    void Write(TPZStream &buf, int withclassid) const override {};

private:

    // Member variable to describe a constant permeability field
    STATE fConstantPermeability = 1.;

    // Member variable to describe a varying permeability field
    OrtoPermeabilityFunctionType fPermeabilityFunction{};
};

// Dummy BC interface class
class TPZMaterial;
class TPZOrtotropicPermeabilityBC : public TPZOrtotropicPermeability {
protected:
    // this method is your chance to verify if the material to which this
    // BC interface applies is compatible with this boundary interface
    // it is called in the method SetMaterial of class TPZBndCondBase
    static void SetMaterialImpl(TPZMaterial *) {}
};

#endif //TPZOrtotropicPermeability_H
