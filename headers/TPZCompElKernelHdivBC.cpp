/**
 * @file
 * @brief Contains the implementation of the TPZCompElHDiv methods.
 */

#include "TPZCompElKernelHdivBC.h"

template<class TVar>
TPZCompElKernelHDivBC<TVar>::TPZCompElKernelHDivBC(TPZMaterial * material, int matid, int type, TPZFMatrix<TVar> &val1,TPZFMatrix<TVar> &val2) : 
TPZRegisterClassId(&TPZCompElKernelHDivBC::ClassId), TPZMatCombinedSpacesBC<TVar>(type){

}

template<class TVar>
TPZCompElKernelHDivBC<TVar>::~TPZCompElKernelHDivBC(){
    this->~TPZMatCombinedSpacesBC<TVar>();
}

template<class TVar>
void TPZCompElKernelHDivBC<TVar>::ContributeBC(TPZMaterialData &data, REAL weight,
                                               TPZFMatrix<TVar> &ek,TPZFMatrix<TVar> &ef, TPZBndCondT<TVar> &bc){
    
    const auto &phi = data.phi;
    const auto &dphi = data.dphix;
    const auto &dphi_ = data.dphi;
	const int phr = phi.Rows();

    if (dphi.Rows() != 1){
        // 
        std::cout << "MatDivFreeBubbles only works for one dimensional elements!" << std::endl;
        std::cout << bc.Id() << std::endl;
        TBase::Print(std::cout);
        DebugStop();
    }
    const int nvars = 1;
    const auto nloads = this->fNumLoadCases;
    const auto &bcNumLoads = dynamic_cast<TPZMatLoadCasesBC<TVar>&>(bc);
    const TPZFMatrix<REAL> &axes=data.axes;
    TPZFMatrix<TVar> v1 = bc.Val1();
    TPZVec<TVar> v2 = bc.Val2();
    TPZFMatrix<TVar> gradU(2,1);
    TPZVec<TVar> u(1);
    

    if(bc.HasForcingFunctionBC()){
        bc.ForcingFunctionBC()(data.x,u,gradU);
        if (bc.Type() == 0)
        {
            v2[0]=u[0];
        }else if (bc.Type()==1){
            v2[0] = gradU(0,0)*axes(0,1) - gradU(1,0)*axes(0,0);
        }
    }

	switch (bc.Type()){		
        // Dirichlet condition
    case 0 : {      
        for(auto in = 0 ; in < phr; in++) {
            ef(in,0) += v2[0] * dphi.GetVal(0,in) * (TVar)weight;
        }//in
        break;
    }
		// Neumann condition
    case 1 : {
        for(auto in = 0 ; in < phr; in++) {
            ef(in,0) += (TVar)TPZMaterial::fBigNumber * v2[0] * (TVar)dphi.GetVal(0,in) * (TVar)weight;
            for (auto jn = 0 ; jn < phr; jn++) {
                ek(in,jn) += (TVar)TPZMaterial::fBigNumber * dphi.GetVal(0,in) * dphi.GetVal(0,jn) * weight;
            }//jn
        }//in
        break;
    }
        //Robin condition
    case 2 : {
        // for(auto iv = 0; iv < nvars; iv++){
        //     for(auto in = 0 ; in < phr; in++) {
        //         for(auto l=0; l < nloads; l++){
        //             const TPZVec<TVar> &v2 = bcNumLoads.GetBCRhsVal(l);
        //             ef(nvars*in+iv,l) +=
        //                 (TVar)TPZMaterial::fBigNumber * v2[nvars*l+iv] * (TVar)phi.GetVal(in,0) * (TVar)weight;
        //         }
        //         for (auto jn = 0 ; jn < phr; jn++) {
        //             ek(nvars*in+iv,nvars*jn+iv) +=
        //                 TPZMaterial::fBigNumber * v1.GetVal(iv,0) * dphi.GetVal(0,in) * dphi.GetVal(0,jn) * weight;
        //         }//jn
        //     }//in
        // }//iv
        DebugStop();
    }		
    default:{
        std::cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " not implemented\n";
    }
	}//switch

}
