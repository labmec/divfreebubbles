#ifndef _TPZCOMPELHCURLNOGRADS_H_
#define _TPZCOMPELHCURLNOGRADS_H_

#include <TPZCompElHCurlFull.h>


/**
   @brief TPZCompElHCurlNoGrads<TSHAPE> is quite a fancy class.
*/
template<class TSHAPE>
class TPZCompElHCurlNoGrads  : public TPZCompElHCurlFull<TSHAPE> {
public:
  //!Default constructor.
  TPZCompElHCurlNoGrads();
  //! Ctor taking mesh, geoel and returning index
  TPZCompElHCurlNoGrads(TPZCompMesh &mesh, TPZGeoEl *gel, int64_t &index);
  /**
   * @brief Number of shapefunctions of the connect associated
   * @param connect connect number
   * @return number of shape functions
   */
  int NConnectShapeF(int connect, int order) const override;

  void InitMaterialData(TPZMaterialData &data) override;

  void ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) override{
    ComputeRequiredDataT(data,qsi);
  }
  void ComputeRequiredData(TPZMaterialDataT<CSTATE> &data, TPZVec<REAL> &qsi) override{
    ComputeRequiredDataT(data,qsi);
  }

  void ReallyComputeSolution(TPZMaterialDataT<STATE>& data) override{
    ReallyComputeSolutionT(data);
  }
  void ReallyComputeSolution(TPZMaterialDataT<CSTATE>& data) override{
    ReallyComputeSolutionT(data);
  }
protected:
  void AdjustConnects();
  void ComputeShape(TPZMaterialData &data, TPZFMatrix<REAL> &phiHCurl);
  template<int DIM>
  void ComputeCurl(TPZMaterialData &data);
  template<class TVar>
  void ComputeRequiredDataT(TPZMaterialDataT<TVar> &data, TPZVec<REAL> &qsi);
  template<class TVar>
  void ReallyComputeSolutionT(TPZMaterialDataT<TVar> &data);
  
};
#endif /* _TPZCOMPELHCURLNOGRADS_H_ */
