//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZOrtotropicPermeability.h"


void TPZOrtotropicPermeability::SetPermeabilityFunction(OrtoPermeabilityFunctionType &perm_function) {
    fPermeabilityFunction = perm_function;
}

TPZManVector<REAL,3> TPZOrtotropicPermeability::GetPermeability(const TPZVec<REAL> &coord) {
  if(fPermeabilityFunction){
    return fPermeabilityFunction(coord);
  }
  else{
    DebugStop();
    return TPZManVector<REAL,3>() = {fConstantPermeability,fConstantPermeability,fConstantPermeability};
  }
  
}

int TPZOrtotropicPermeability::ClassId() const {
    return Hash("TPZOrtotropicPermeability");
}
