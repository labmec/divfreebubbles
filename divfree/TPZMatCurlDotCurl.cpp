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

void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight,
                    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                    TPZBndCondT<STATE> &bc) 
{
    TPZFMatrix<REAL> phiQ = data.phi;
    int phrq = phiQ.Rows();

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0, 0);

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(3, 1);

        bc.ForcingFunctionBC()(data.x, res, gradu);

        if (bc.Type() == 0 || bc.Type() == 4) {
            v2 = res[0];
        } else {
            DebugStop();
        }
    } else {
        v2 = bc.Val2()[0];
    }


    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
    
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += (-1.) * v2 * phiQ(iq, 0) * weight;
            }
    
            break;

    }
}

TPZMaterial * TPZMatCurlDotCurl::NewMaterial() const
{
  return new TPZMatCurlDotCurl(*this);
}