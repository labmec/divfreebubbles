#include "TPZMatDivFreeBubbles.h"
#include "TPZMaterialDataT.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"
#include <iostream>
#include <fstream>
template<class TVar>
TPZMatDivFreeBubbles<TVar>::TPZMatDivFreeBubbles(int id, int dim) :
    TPZRegisterClassId(&TPZMatDivFreeBubbles::ClassId),
    TBase(id), fDim(dim), fSol(0)
{
}

template<class TVar>
TPZMaterial * TPZMatDivFreeBubbles<TVar>::NewMaterial() const{
	return new TPZMatDivFreeBubbles(*this);
}

template<class TVar>
void TPZMatDivFreeBubbles<TVar>::Contribute(const TPZMaterialDataT<TVar> &data,
                                       REAL weight,
                                       TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef){
	
	const auto nLoads = this->fNumLoadCases;
    TPZManVector<TVar,10> force(nLoads,0.); 
    if(this->HasForcingFunction()){
        DebugStop();
    }
    const auto &phi = data.phi;
    const auto &dphi = data.dphix;
    const auto nshape = data.phi.Rows();
	for(int i = 0; i < nshape; i++){
		for(int j = 0; j < nshape; j++){
            STATE dphiIdphiJ = 0;
            for(int x = 0; x < fDim; x++){
                dphiIdphiJ += dphi.GetVal(x,i) * dphi.GetVal(x,j);
            }
            ek(i, j) += weight*dphiIdphiJ;
        }//forj
        // for(auto l = 0; l < nLoads; l++)
        //     for(int x = 0; x < fDim; x++)
        //         ef(fDim*i+x,l) += weight*fScale*phi.GetVal(i,0)*force[l];
            
    }//for i
}

template<class TVar>
void TPZMatDivFreeBubbles<TVar>::ContributeBC(const TPZMaterialDataT<TVar> &data,
                                         REAL weight,
                                         TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                                         TPZBndCondT<TVar> &bc)
{
	
	const auto &phi = data.phi;
    const auto &dphi = data.dphix;
    const auto &dphi_ = data.fDPhi;
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
            ef(in,0) += v2[0] * (TVar)dphi.GetVal(0,in) * (TVar)weight;
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

template<class TVar>
void TPZMatDivFreeBubbles<TVar>::GetSolDimensions(uint64_t &u_len,
                                             uint64_t &du_row,
                                             uint64_t &du_col) const
{
    u_len=1;
    du_row=3;
    du_col=1;
}


template<class TVar>
int TPZMatDivFreeBubbles<TVar>::VariableIndex(const std::string &name) const{
	if(!strcmp("Solution",name.c_str())) return ESolution;
    if(!strcmp("Derivative",name.c_str())) return EDerivative;
	return TPZMaterial::VariableIndex(name);
}

template<class TVar>
int TPZMatDivFreeBubbles<TVar>::NSolutionVariables(int var) const{
	if(var == ESolution) return 1;
    else if (var == EDerivative) {
        return 3;
    }
	
    return TPZMaterial::NSolutionVariables(var);
}

template<class TVar>
void TPZMatDivFreeBubbles<TVar>::Solution(const TPZMaterialDataT<TVar> &data,
                                     int var, TPZVec<TVar> &solOut)
{
    TPZFMatrix<TVar> dsolX(3,1);
    const auto &sol = data.sol[0];
    const auto &dsol = data.dsol[0];
    TPZAxesTools<TVar>::Axes2XYZ(dsol,dsolX,data.axes);

    
    // std::cout << "SOL " << sol << ", x= " << data.x << std::endl;
	if (var == ESolution){
        solOut.Resize(sol.size());
        for (int i=0; i<sol.size(); i++) {
            solOut[i] = sol[i];
        }
		return;
	}
    if (var == EDerivative) {
        solOut.Resize(3);
        // for (int i=0; i<3; i++) {
        //     solOut[i] = dsolX.GetVal(i,0);
        // }
        //Changed here to properly orient the flux vector
        solOut[0] = -dsolX.GetVal(1,0);
        solOut[1] =  dsolX.GetVal(0,0);
        solOut[2] = 0.;
        return;
    }
}

template<class TVar>
void TPZMatDivFreeBubbles<TVar>::Errors(const TPZMaterialDataT<TVar>&data,
                                 TPZVec<REAL> &values) {
    const auto &x = data.x;
    const auto &u = data.sol[this->fPostProcIndex];
    const auto &dudx = data.dsol[this->fPostProcIndex];
    const auto &axes = data.axes;

    TPZManVector<TVar,1> u_exact={0.};
    TPZFNMatrix<3,TVar> du_exact(3,1,0.);
    this->ExactSol()(x,u_exact,du_exact);
    values.Resize(this->NEvalErrors());
    values.Fill(0.0);
    TPZManVector<TVar> sol(3,0.),dsol(3,0.);
    TPZFNMatrix<3,TVar> gradu(3,1), gradu_(3,1);
    gradu_.Zero();
    TPZAxesTools<TVar>::Axes2XYZ(dudx,gradu_,axes);
    gradu(0,0) = -gradu_(1,0);
    gradu(1,0) =  gradu_(0,0);
    
    //values[0] : error in H1 norm
    //values[1] : eror in L2 norm
    //values[2] : erro in H1 semi-norm
    TVar diff = (u[0] - u_exact[0]);
    if constexpr (is_complex<TVar>::value){
        values[1]  = std::real((diff*std::conj(diff)));
    }else{
        values[1]  = diff*diff;
    }
  
    values[2] = 0.;

    for(auto id=0; id<fDim; id++) {
      diff = (gradu(id) - du_exact(id,0));
      if constexpr(is_complex<TVar>::value){
          values[2]  += std::real(diff*std::conj(diff));
      }else{
          values[2]  += diff*diff;
      }
    }
    values[0]  = values[1]+values[2];
}

template<class TVar>
int TPZMatDivFreeBubbles<TVar>::ClassId() const{
    return Hash("TPZMatDivFreeBubbles") ^ TBase::ClassId() << 1;
}


template class TPZMatDivFreeBubbles<STATE>;
template class TPZMatDivFreeBubbles<CSTATE>;