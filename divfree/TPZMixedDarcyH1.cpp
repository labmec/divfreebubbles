//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZMixedDarcyH1.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"
// #include "TPZLapack.h"

#define USEBLAS

TPZMixedDarcyH1::TPZMixedDarcyH1() : TPZRegisterClassId(&TPZMixedDarcyH1::ClassId),
                                         TBase(), fDim(-1) {}

[[maybe_unused]] TPZMixedDarcyH1::TPZMixedDarcyH1(int id, int dim) : TPZRegisterClassId(&TPZMixedDarcyH1::ClassId),
                                                        TBase(id), fDim(dim)
{
}

/**
         copy constructor
 */
TPZMixedDarcyH1::TPZMixedDarcyH1(const TPZMixedDarcyH1 &copy) : TBase(copy), fDim(copy.fDim)
{
    
}
/**
         copy constructor
 */
TPZMixedDarcyH1 &TPZMixedDarcyH1::operator=(const TPZMixedDarcyH1 &copy)
{
    TBase::operator=(copy);
    fDim = copy.fDim;
    return *this;
}


void TPZMixedDarcyH1::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef) {

    STATE force = 0;
    if (fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction(datavec[1].x, res);
        force = res[0];
    }
    const STATE perm = GetPermeability(datavec[0].x);
    const STATE inv_perm = 1 / perm;
    
    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFMatrix<REAL> &divQ = datavec[0].divphi;
    TPZFNMatrix<9, REAL> dphiPXY(3, dphiP.Cols());
    TPZFNMatrix<9, REAL> dphiQXY(3, dphiQ.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPXY, datavec[1].axes);
    TPZAxesTools<REAL>::Axes2XYZ(dphiQ, dphiQXY, datavec[0].axes);

    REAL &faceSize = datavec[0].HSize;

    // std::cout<<"phiQ = " << phiQ << std:: endl << "phip = "<< phip<<std::endl;

    int phrq, phrp;
    phrp = phip.Rows();
    phrq = phiQ.Rows();

    int nactive = 0;
    for (const auto &i : datavec) {
        if (i.fActiveApproxSpace) {
            nactive++;
        }
    }
#ifdef PZDEBUG
    // if (nactive >= 4) {
    //     if(nactive%2 != 0) DebugStop();
    //     int numavg = (nactive-2)/2;
    //     for(int iavg = 0; iavg<numavg; iavg++)
    //     {
    //         if(datavec[2+2*iavg].phi.Rows() != 1) DebugStop();
    //         if(datavec[2+2*iavg+1].phi.Rows() != 1) DebugStop();
    //     }
    //     if (phrp + phrq + 2*numavg != ek.Rows()) {
    //         DebugStop();
    //     }
    // } else {
    //     if (phrp + phrq != ek.Rows()) {
    //         DebugStop();
    //     }
    // }
#endif
    
    int dim = this->Dimension();
    //Calculate the matrix contribution for flux. Matrix A
    for (int iq = 0; iq < phrq; iq++) {
        for (int jq = 0; jq < phrq; jq++) {
            for (int k = 0; k < dim; k++) {
                for (int l = 0; l < dim; l++) {
                    double K = dphiQXY(l,iq) * dphiQXY(k,jq);
                    if (k==l) for (int m = dim; m--; ) K += dphiQXY(m,iq) * dphiQXY(m,jq);
                    ek(dim*iq+k,dim*jq+l) += K * weight;
                }
            }
        }
        for (int jp = 0; jp < phrp; jp++){
            for (int k = 0; k < dim; k++) {
                double Q = dphiQXY(k,iq) * phip[jp];
                ek(dim*phrq+jp,dim*iq+k) -= Q * weight;
                ek(dim*iq+k,dim*phrq+jp) -= Q * weight;
            }
        }
        
    }
    // for (int ip = 0; ip < phrp; ip++) {
    //     ek(dim*phrq+ip,dim*phrq+ip) += 1 * weight;
    // }
   
}

void TPZMixedDarcyH1::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    int dim = Dimension();

    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    int phrq = phiQ.Rows();
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    int phrp = phip.Rows();

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0, 0);
    REAL u_D = 0;
    REAL normflux = 0.;
    TPZManVector<STATE> res(3);
    TPZFNMatrix<9, STATE> gradu(3, 1);
    if (bc.HasForcingFunctionBC()) {
        
        bc.ForcingFunctionBC()(datavec[0].x, res, gradu);

        const STATE perm = GetPermeability(datavec[0].x);

        for (int i = 0; i < 3; i++) {
            normflux += datavec[0].normal[i] * perm * gradu(i, 0);
        }

        if (bc.Type() == 0 || bc.Type() == 4) {
            v2 = res[0];
            u_D = res[0];
            normflux *= (-1.);
        } else if (bc.Type() == 1 || bc.Type() == 2) {
            v2 = -normflux;
            if (bc.Type() == 2) {
                v2 = -res[0] + v2 / v1;
            }
        } else {
            DebugStop();
        }
    } else {
        gradu(0,0) = bc.Val2()[0];
        gradu(1,0) = bc.Val2()[1];
    }

    int nstate = this->NStateVariables();

    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
            for (int iq = 0; iq < phrq; iq++) {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                for (int k = 0; k < dim; k++){
                    ef(dim*iq+k, 0) += TPZMaterial::fBigNumber * gradu(k,0) * phiQ(iq, 0) * weight;
                    for (int jq = 0; jq < phrq; jq++) {
                        ek(dim*iq+k,dim*jq+k) += TPZMaterial::fBigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
                    }
                }
                // ef(iq, 0) += (-1.) * v2 * phiQ(iq, 0) * weight;
            }
            break;

        case 1 :            // Neumann condition
            // for(auto iv = 0; iv < dim; iv++){
            //     for(auto ip = 0 ; ip < phrq; ip++) {
            //             ef(dim*ip+iv,0) += res[0] * phiQ.GetVal(ip,0) * weight;
            //     }//in
            // }//iv
            break;

        case 2 :            // mixed condition
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += v2 * phiQ(iq, 0) * weight;
                for (int jq = 0; jq < phrq; jq++) {
                    ek(iq, jq) += weight / v1 * phiQ(iq, 0) * phiQ(jq, 0);
                }
            }
            break;

        case 4:
            //this case implemented the general Robin boundary condition
            // sigma.n = Km(u-u_D)+g
            //val1(0,0) = Km
            //val2(1,0) = g
            if (IsZero(bc.Val1()(0, 0))) {

                for (int iq = 0; iq < phrq; iq++) {
                    ef(iq, 0) += TPZMaterial::fBigNumber * normflux * phiQ(iq, 0) * weight;
                    for (int jq = 0; jq < phrq; jq++) {
                        ek(iq, jq) += TPZMaterial::fBigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
                    }
                }

            } else {

                REAL InvKm = 1. / bc.Val1()(0, 0);
                REAL g = normflux;
                for (int in = 0; in < phiQ.Rows(); in++) {
                    //<(InvKm g - u_D)*(v.n)
                    ef(in, 0) += (STATE) (InvKm * g - u_D) * phiQ(in, 0) * weight;
                    for (int jn = 0; jn < phiQ.Rows(); jn++) {
                        //InvKm(sigma.n)(v.n)
                        ek(in, jn) += (STATE) (InvKm * phiQ(in, 0) * phiQ(jn, 0) * weight);
                    }
                }
            }

            break;

    }
}

void TPZMixedDarcyH1::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut) {
    solOut.Resize(this->NSolutionVariables(var));
    solOut.Fill(0.);
    TPZManVector<STATE, 10> SolP, SolQ;
    const STATE perm = GetPermeability(datavec[0].x);
    const STATE inv_perm = 1 / perm;

    // SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];
    if(SolP.size() == 0) SolP.Resize(1,0.);

    if (var == 1) { //function (state variable Q)
        for (int i = 0; i < Dimension(); i++) {
            solOut[i] = datavec[0].sol[0][i];

        }
        return;
    }

    if (var == 2) {
        solOut[0] = SolP[0];//function (state variable p)
        return;
    }

    if (var == 3) {
        solOut[0] = datavec[0].dsol[0](0, 0);
        solOut[1] = datavec[0].dsol[0](1, 0);
        solOut[2] = datavec[0].dsol[0](2, 0);
        return;
    }

    if (var == 4) {
        solOut[0] = datavec[0].dsol[0](0, 1);
        solOut[1] = datavec[0].dsol[0](1, 1);
        solOut[2] = datavec[0].dsol[0](2, 1);
        return;
    }

    if (var == 5) {
        solOut[0] = datavec[0].dsol[0](0, 0) + datavec[0].dsol[0](1, 1);
        return;
    }

    // Exact solution
    if (var == 6) {
        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> flux(3, 1);
        if (fExactSol) {
            fExactSol(datavec[0].x, exactSol, flux);
        }
        solOut[0] = exactSol[0];
        return;
    } // var6

    if (var == 7) {

        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> gradu(3, 1);

        if (fExactSol) {
            fExactSol(datavec[0].x, exactSol, gradu);
        }

        for (int i = 0; i < 3; i++) {
            solOut[i] = -perm * gradu(i, 0);
        }

        return;
    } // var7

    if (var == 8) {
        solOut[0] = datavec[1].p;
        return;
    }

    if (var == 9) {

        if(datavec[1].fShapeType == TPZMaterialData::EEmpty) return;
        TPZFNMatrix<9, REAL> dsoldx(3, 1.,0.);
        TPZFNMatrix<9, REAL> dsoldaxes(fDim, 1,0.);

        dsoldaxes = datavec[1].dsol[0];
        TPZAxesTools<REAL>::Axes2XYZ(dsoldaxes, dsoldx, datavec[1].axes);

        for (int i = 0; i < fDim; i++) {
            solOut[i] = dsoldx(i, 0);
        }

        return;
    }

    if (var == 10) {
        solOut[0] = 0.;
        // solOut[0]=datavec[0].dsol[0](0,0)+datavec[0].dsol[0](1,1);
        for (int j = 0; j < fDim; j++) {
            solOut[0] += datavec[0].dsol[0](j, j);
        }
        return;
    }

    if (var == 11) {
        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> flux(3, 1);
        fExactSol(datavec[0].x, exactSol, flux);
        solOut[0] = flux(2, 0);
        return;
    }
    if (var == 12) {
        for (int i = 0; i < fDim; i++) {
            solOut[i] = 0.;
        }
        for (int i = 0; i < fDim; i++) {
            solOut[i] -= inv_perm * datavec[0].sol[0][i];
        }
        return;
    }
    if (var == 13) {
        solOut[0] = perm;
        return;
    }

    if (datavec.size() == 4) {
        if (var == 14) {
            solOut[0] = datavec[2].sol[0][0];
            return;
        }
        if (var == 15) {
            solOut[0] = datavec[3].sol[0][0];
            return;
        }

    }

    if (var == 16) { //ExactFluxShiftedOrigin
        // Solution EArcTan returns NAN for (x,y) == (0,0). Replacing data.x by
        // inf solves this problem.
        STATE infinitesimal = 0.0000000001;
        TPZManVector<REAL, 3> inf = {infinitesimal, infinitesimal, infinitesimal};

        TPZVec<STATE> exactSol(1);
        TPZFNMatrix<3, STATE> gradu(3, 1);

        if (fExactSol) {
            if (datavec[0].x[0] == 0. && datavec[0].x[1] == 0.) {
                fExactSol(inf, exactSol, gradu);
            } else {
                fExactSol(datavec[0].x, exactSol, gradu);
            }
        }
        for (int i = 0; i < 3; i++) {
            solOut[i] = -perm * gradu(i, 0);
        }

        return;
    }
}

void TPZMixedDarcyH1::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {

    /**
     * datavec[0]= Flux
     * datavec[1]= Pressure
     *
     * Errors:
     * [0] L2 for pressure
     * [1] L2 for flux
     * [2] L2 for div(flux)
     * [3] Grad pressure (Semi H1)
     * [4] Hdiv norm
    **/
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);

    TPZManVector<STATE, 3> fluxfem(3), pressurefem(1,0);
    fluxfem = data[0].sol[0];
    STATE divsigmafem = data[0].divsol[0][0];

    auto dsol = data[1].dsol;

    TPZManVector<STATE,1> divsigma(1,0.);

    TPZManVector<STATE,1> u_exact(1, 0);
    TPZFMatrix<STATE> du_exact(3, 1, 0);
    if (this->fExactSol) {
        this->fExactSol(data[0].x, u_exact, du_exact);
    }
    if (this->fForcingFunction) {
        this->fForcingFunction(data[0].x, divsigma);
    }

    REAL residual = (divsigma[0] - divsigmafem) * (divsigma[0] - divsigmafem);
    if(data[1].sol[0].size())
        pressurefem[0] = data[1].sol[0][0];

    const STATE perm = GetPermeability(data[0].x);
    const STATE inv_perm = 1 / perm;

    TPZManVector<STATE, 3> gradpressurefem(3, 0.);
    this->Solution(data, VariableIndex("GradPressure"), gradpressurefem);

    TPZManVector<STATE, 3> fluxexact(3, 0);
    TPZManVector<STATE, 3> gradpressure(3, 0);
    for (int i = 0; i < 3; i++) {
        gradpressure[i] = du_exact[i];
        fluxexact[i] = -perm * gradpressure[i];
    }

    REAL L2flux = 0., L2grad = 0.;
    for (int i = 0; i < 3; i++) {
        L2flux += (fluxfem[i] - fluxexact[i]) * inv_perm * (fluxfem[i] - fluxexact[i]);
        L2grad += (du_exact[i] - gradpressurefem[i]) * (du_exact[i] - gradpressurefem[i]);
    }
    errors[0] = (pressurefem[0] - u_exact[0]) * (pressurefem[0] - u_exact[0]);//L2 error for pressure
    errors[1] = L2flux;//L2 error for flux
    errors[2] = residual;//L2 for div
    errors[3] = L2grad;
    errors[4] = L2flux + residual;
}

int TPZMixedDarcyH1::VariableIndex(const std::string &name) const {
    if (!strcmp("Flux", name.c_str())) return 1;
    if (!strcmp("Pressure", name.c_str())) return 2;
    if (!strcmp("GradFluxX", name.c_str())) return 3;
    if (!strcmp("GradFluxY", name.c_str())) return 4;
    if (!strcmp("DivFlux", name.c_str())) return 5;
    if (!strcmp("ExactPressure", name.c_str())) return 6;
    if (!strcmp("ExactFlux", name.c_str())) return 7;
    if (!strcmp("POrder", name.c_str())) return 8;
    if (!strcmp("GradPressure", name.c_str())) return 9;
    if (!strcmp("Divergence", name.c_str())) return 10;
    if (!strcmp("ExactDiv", name.c_str())) return 11;
    if (!strcmp("Derivative", name.c_str())) return 12;
    if (!strcmp("Permeability", name.c_str())) return 13;
    if (!strcmp("g_average", name.c_str())) return 14;
    if (!strcmp("u_average", name.c_str())) return 15;
    if (!strcmp("ExactFluxShiftedOrigin", name.c_str())) return 16;
    DebugStop();
    return -1;
}

int TPZMixedDarcyH1::NSolutionVariables(int var) const {
    if (var == 1) return 3;
    if (var == 2) return 1;
    if (var == 3) return 3;
    if (var == 4) return 3;
    if (var == 5) return 1;
    if (var == 6) return 1;
    if (var == 7) return 3;
    if (var == 8) return 1;
    if (var == 9) return 3;
    if (var == 10 || var == 11) return 1;
    if (var == 12) return 3;
    if (var == 13) return 1;
    if (var == 14) return 1;
    if (var == 15) return 1;
    if (var == 16) return 3;
    DebugStop();
    return -1;
}

void TPZMixedDarcyH1::SetDimension(int dim) {
    if (dim > 3 || dim < 1) DebugStop();
    fDim = dim;
}

int TPZMixedDarcyH1::ClassId() const {
    return Hash("TPZMixedDarcyH1") ^ TBase::ClassId() << 1;
}

TPZMaterial *TPZMixedDarcyH1::NewMaterial() const {
    return new TPZMixedDarcyH1(*this);
}

void TPZMixedDarcyH1::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << this->Id() << "\n";
    out << "Dimension: " << this->Dimension() << "\n\n";
}

void TPZMixedDarcyH1::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    int nref = datavec.size();
    for (int i = 0; i < nref; i++) {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
        datavec[i].fNeedsHSize = false;
    }
}

void
TPZMixedDarcyH1::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
    // default is no specific data requirements
    int nref = datavec.size();
    for (int iref = 0; iref < nref; iref++) {
        datavec[iref].SetAllRequirements(false);
        datavec[iref].fNeedsSol = false;
    }
    datavec[0].fNeedsNormal = true;
    if (type == 50) {
        for (int iref = 0; iref < nref; iref++) {
            datavec[iref].fNeedsSol = false;
        }
    }
}
#undef USEBLAS
