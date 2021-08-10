#ifndef _TPZMATCURLDOTCURL_H_
#define _TPZMATCURLDOTCURL_H_

#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"

/**
   @brief TPZMatCurlDotCurl will perform $curl(phi_i)\cdot curl(phi_j)$.
   It is used for verifying if the curl of a given hcurl-conforming
   approximation space forms a linearly independent set.
   @note It currently only supports 3D problems.
*/
class TPZMatCurlDotCurl  : public TPZMatBase<STATE,TPZMatSingleSpaceT<STATE>>
{
  using TBase = TPZMatBase<STATE,TPZMatSingleSpaceT<STATE>>;
public:
  
  //!Default constructor
  TPZMatCurlDotCurl() = default;

  /**
	 * @brief Constructor setting a material identifier.
	 * @param id material id
	 */
  TPZMatCurlDotCurl(const int id);

  //! Material name
  std::string Name() const override { return "TPZMatCurlDotCurl"; }
  
  //! Material spatial dimension (three).
  int Dimension() const override { return this->fDim;}
  //! Number of state variables (one).
  int NStateVariables() const override { return this->fNStateVars; }
  //! Class identifier.
  int ClassId() const override;
  //! Contribution of the FEM matrices at each volumetric integration point.
  void Contribute(const TPZMaterialDataT<STATE> &data,
                  REAL weight,
                  TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
  //! Contribution of the FEM matrices at each BC integration point (nothing).
  void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight,
                    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                    TPZBndCondT<STATE> &bc) override {}
  //! To create another material of the same type
	TPZMaterial * NewMaterial() const override;
private:
  //! Problem dimension
	static constexpr int fDim{3};
	
	//! Number of state variables
	static constexpr int fNStateVars{1};
};

#endif /* _TPZMATCURLDOTCURL_H_ */
