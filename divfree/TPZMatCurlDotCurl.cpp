#include "TPZMatCurlDotCurl.h"
#include <TPZMaterialDataT.h>


TPZMatCurlDotCurl::TPZMatCurlDotCurl(const int id) : TBase(id)
{}

int TPZMatCurlDotCurl::ClassId() const
{
  return Hash("TPZMatCurlDotCurl") ^ TBase::ClassId() << 1;
}

void TPZMatCurlDotCurl::Contribute(const TPZMaterialDataT<STATE> &data,
                REAL weight,
                TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{

  const int nshape = data.curlphi.Cols();

  const auto &cphi = data.curlphi;
  for(int i = 0; i < nshape; i++){
    for(int j = 0; j < nshape; j++){
      STATE curlIcurlJ = 0;
      for(int x = 0; x < fDim; x++){
        curlIcurlJ += cphi.GetVal(x,i) * cphi.GetVal(x,j);
      }
      const int posI = i;
      const int posJ = j;
      ek(posI, posJ) += weight*curlIcurlJ;
    }//for j
  }//for i
}

TPZMaterial * TPZMatCurlDotCurl::NewMaterial() const
{
  return new TPZMatCurlDotCurl(*this);
}