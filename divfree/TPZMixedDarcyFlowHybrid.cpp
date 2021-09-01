//
// Created by Gustavo Batistela on 5/13/21.
//

#include "TPZMixedDarcyFlowHybrid.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"

TPZMixedDarcyFlowHybrid::TPZMixedDarcyFlowHybrid() : TPZRegisterClassId(&TPZMixedDarcyFlowHybrid::ClassId),
                                                        TBase(), fDim(-1) {}

[[maybe_unused]] TPZMixedDarcyFlowHybrid::TPZMixedDarcyFlowHybrid(int id, int dim) : TPZRegisterClassId(&TPZMixedDarcyFlowHybrid::ClassId),
                                                        TBase(id), fDim(dim) {}

void TPZMixedDarcyFlowHybrid::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef) {

    STATE force = 0;
    if (fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction(datavec[1].x, res);
        force = res[0];
    }

    TPZFNMatrix<1, STATE> K(3, 3, 0);
    TPZFNMatrix<1, STATE> InvK(3, 3, 0);

    fPermeabilityFunction(datavec[0].x, K, InvK);
    const REAL perm = K(0, 0);
    const REAL inv_perm = 1 / perm;

    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFNMatrix<9, REAL> dphiPXY(3, dphiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPXY, datavec[1].axes);

    REAL &faceSize = datavec[0].HSize;

    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();

    int nactive = 0;
    for (const auto &i : datavec) {
        if (i.fActiveApproxSpace) {
            nactive++;
        }
    }
#ifdef PZDEBUG
    if (nactive == 4) {
        int phrgb = datavec[2].phi.Rows();
        int phrub = datavec[3].phi.Rows();
        if (phrp + phrq + phrgb + phrub != ek.Rows()) {
            DebugStop();
        }
    } else {
        if (phrp + phrq != ek.Rows()) {
            DebugStop();
        }
    }
#endif

    //Calculate the matrix contribution for flux. Matrix A
    for (int iq = 0; iq < phrq; iq++) {
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3, REAL> ivec(3, 1, 0.);
        for (int id = 0; id < 3; id++) {
            ivec(id, 0) = datavec[0].fDeformedDirections(id, ivecind);
        }

        TPZFNMatrix<3, REAL> ivecZ(3, 1, 0.);
        TPZFNMatrix<3, REAL> jvecZ(3, 1, 0.);
        for (int jq = 0; jq < phrq; jq++) {
            TPZFNMatrix<3, REAL> jvec(3, 1, 0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;

            for (int id = 0; id < fDim; id++) {
                jvec(id, 0) = datavec[0].fDeformedDirections(id, jvecind);
            }

            //dot product between Kinv[u]v
            jvecZ.Zero();
            for (int id = 0; id < 3; id++) {
                jvecZ(id, 0) += inv_perm * jvec(id, 0);
            }
            REAL prod1 = ivec(0, 0) * jvecZ(0, 0) + ivec(1, 0) * jvecZ(1, 0) + ivec(2, 0) * jvecZ(2, 0);
            ek(iq, jq) += weight * phiQ(ishapeind, 0) * phiQ(jshapeind, 0) * prod1;
        }
    }


    // Coupling terms between flux and pressure. Matrix B
    for (int iq = 0; iq < phrq; iq++) {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;

        TPZFNMatrix<3, REAL> ivec(3, 1, 0.);
        for (int id = 0; id < 3; id++) {
            ivec(id, 0) = datavec[0].fDeformedDirections(id, ivecind);
            //ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
            //ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
        }
        TPZFNMatrix<3, REAL> axesvec(3, 1, 0.);
        datavec[0].axes.Multiply(ivec, axesvec);

        REAL divwq = 0.;
        for (int iloc = 0; iloc < fDim; iloc++) {
            divwq += axesvec(iloc, 0) * dphiQ(iloc, ishapeind);
        }
        for (int jp = 0; jp < phrp; jp++) {

            REAL fact = (-1.) * weight * phip(jp, 0) * divwq;
            // Matrix B
            ek(iq, phrq + jp) += fact;

            // Matrix B^T
            ek(phrq + jp, iq) += fact;

        }
    }

    // source term related to the pressure equation
    for (int ip = 0; ip < phrp; ip++) {
        ef(phrq + ip, 0) += (-1.) * weight * force * phip(ip, 0);
    }

    if (nactive == 4) {
        for (int ip = 0; ip < phrp; ip++) {
            ek(phrq + ip, phrq + phrp) += phip(ip, 0) * weight;
            ek(phrq + phrp, phrq + ip) += phip(ip, 0) * weight;
        }
        ek(phrp + phrq + 1, phrq + phrp) += -weight;
        ek(phrq + phrp, phrp + phrq + 1) += -weight;
    }
}

void TPZMixedDarcyFlowHybrid::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) {

    int dim = Dimension();

    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phiP = datavec[1].phi;
    int phrq = phiQ.Rows();
    int phrp = phiP.Rows();

    REAL v2 = bc.Val2()[0];
    REAL v1 = bc.Val1()(0, 0);
    REAL u_D = 0;
    REAL normflux = 0.;

    if (bc.HasForcingFunctionBC()) {
        TPZManVector<STATE> res(3);
        TPZFNMatrix<9, STATE> gradu(3, 1);

        TPZMaterialDataT<STATE> dataAux;
        if (phrq > 0) dataAux = datavec[0];
        if (phrp > 0) dataAux = datavec[1];

        bc.ForcingFunctionBC()(dataAux.x, res, gradu);

        TPZFNMatrix<9, STATE> K(3, 3, 0);
        TPZFNMatrix<9, STATE> InvK(3, 3, 0);
        fPermeabilityFunction(dataAux.x, K, InvK);
        REAL perm = K(0, 0);

        for (int i = 0; i < 3; i++) {
            normflux += dataAux.normal[i] * perm * gradu(i, 0);
        }

        if (bc.Type() == 0 || bc.Type() == 4) {
            v2 = res[0];
            u_D = res[0];
            normflux *= (-1.);
        } else if (bc.Type() == 1 || bc.Type() == 2) {
            v2 = -normflux;
            if (phrp > 0) v2 = gradu(0,0)*dataAux.axes(0,1) - gradu(1,0)*dataAux.axes(0,0);
            if (bc.Type() == 2) {
                v2 = -res[0] + v2 / v1;
            }
        } else {
            DebugStop();
        }
    } else {
        v2 = bc.Val2()[0];
    }

    // TPZFNMatrix<220,REAL> dphix(3,datavec[1].dphix.Cols());
    // TPZFMatrix<REAL> &dphi = datavec[1].dphix;
    // TPZAxesTools<REAL>::Axes2XYZ(dphi, dphix, datavec[1].axes);

    switch (bc.Type()) {
        case 0 :        // Dirichlet condition
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += (-1.) * v2 * phiQ(iq, 0) * weight;
            }
            //New implementation for Div-free bubbles
            for (int ip = 0; ip < phrp; ip++) {
                DebugStop();
                // ef(ip, 0) += (-1.) * v2 * phiP(ip, 0) * weight;
            }
            break;

        case 1 :            // Neumann condition
            for (int iq = 0; iq < phrq; iq++) {
                ef(iq, 0) += TPZMaterial::fBigNumber * v2 * phiQ(iq, 0) * weight;
                for (int jq = 0; jq < phrq; jq++) {
                    ek(iq, jq) += TPZMaterial::fBigNumber * phiQ(iq, 0) * phiQ(jq, 0) * weight;
                }
            }
            //New implementation for Div-free bubbles
            // std::cout << phiP<<std::endl;
            for (int ip = 0; ip < phrp; ip++) {
                ef(ip, 0) += (-1.) * v2 * phiP(ip, 0) * weight;
            }
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

void TPZMixedDarcyFlowHybrid::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &solOut) {
    solOut.Resize(this->NSolutionVariables(var));

    TPZManVector<STATE, 10> SolP, SolQ;

    TPZFNMatrix<9, STATE> K(3, 3, 0);
    TPZFNMatrix<9, STATE> InvK(3, 3, 0);

    fPermeabilityFunction(datavec[0].x, K, InvK);
    const REAL perm = K(0, 0);
    const REAL inv_perm = 1 / perm;

    // SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];

    if (var == 1) { //function (state variable Q)
        for (int i = 0; i < 3; i++) {
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

        TPZFNMatrix<9, REAL> dsoldx(fDim, 1.);
        TPZFNMatrix<9, REAL> dsoldaxes(fDim, 1);

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
        solOut[0] = K(0, 0);
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

void TPZMixedDarcyFlowHybrid::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {

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

    TPZManVector<STATE, 3> fluxfem(3), pressurefem(1);
    fluxfem = data[0].sol[0];
    STATE divsigmafem = data[0].divsol[0][0];

    auto dsol = data[1].dsol;

    TPZVec<STATE> divsigma(1);

    TPZVec<STATE> u_exact(1, 0);
    TPZFMatrix<STATE> du_exact(3, 1, 0);
    if (this->fExactSol) {
        this->fExactSol(data[0].x, u_exact, du_exact);
    }
    if (this->fForcingFunction) {
        this->fForcingFunction(data[0].x, divsigma);
    }

    REAL residual = (divsigma[0] - divsigmafem) * (divsigma[0] - divsigmafem);
    // pressurefem[0] = data[1].sol[0][0];

    TPZFNMatrix<9, STATE> K(3, 3, 0);
    TPZFNMatrix<9, STATE> InvK(3, 3, 0);

    fPermeabilityFunction(data[0].x, K, InvK);
    const REAL perm = K(0, 0);
    const REAL inv_perm = 1 / perm;

    TPZManVector<STATE, 3> gradpressurefem(3, 0.);
    // this->Solution(data, VariableIndex("GradPressure"), gradpressurefem);

    TPZManVector<STATE, 3> fluxexactneg(3, 0);
    TPZManVector<STATE, 3> gradpressure(3, 0);
    for (int i = 0; i < 3; i++) {
        gradpressure[i] = du_exact[i];
        fluxexactneg[i] = -perm * gradpressure[i];
    }

    if (data.size()>2) {
        for (int i = 0; i < 3; i++) {
            fluxexactneg[i] = -data[2].sol[0][i];
        }
    }

    REAL L2flux = 0., L2grad = 0.;
    for (int i = 0; i < 3; i++) {
        L2flux += (fluxfem[i] + fluxexactneg[i]) * inv_perm * (fluxfem[i] + fluxexactneg[i]);
        L2grad += (du_exact[i] - gradpressurefem[i]) * (du_exact[i] - gradpressurefem[i]);
    }
    errors[0] = (pressurefem[0] - u_exact[0]) * (pressurefem[0] - u_exact[0]);//L2 error for pressure
    errors[1] = L2flux;//L2 error for flux
    errors[2] = residual;//L2 for div
    errors[3] = L2grad;
    errors[4] = L2flux + residual;
}

int TPZMixedDarcyFlowHybrid::VariableIndex(const std::string &name) const {
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

int TPZMixedDarcyFlowHybrid::NSolutionVariables(int var) const {
    if (var == 1) return 3;
    if (var == 2) return 1;
    if (var == 3) return 3;
    if (var == 4) return 3;
    if (var == 5) return 1;
    if (var == 6) return 1;
    if (var == 7) return 3;
    if (var == 8) return 1;
    if (var == 9) return 3;
    if (var == 10 || var == 41) return 1;
    if (var == 12) return 3;
    if (var == 13) return 1;
    if (var == 14) return 1;
    if (var == 15) return 1;
    if (var == 16) return 3;
    DebugStop();
    return -1;
}

void TPZMixedDarcyFlowHybrid::SetDimension(int dim) {
    if (dim > 3 || dim < 1) DebugStop();
    fDim = dim;
}

int TPZMixedDarcyFlowHybrid::ClassId() const {
    return Hash("TPZMixedDarcyFlowHybrid") ^ TBase::ClassId() << 1;
}

TPZMaterial *TPZMixedDarcyFlowHybrid::NewMaterial() const {
    return new TPZMixedDarcyFlowHybrid(*this);
}

void TPZMixedDarcyFlowHybrid::Print(std::ostream &out) const {
    out << "Material Name: " << this->Name() << "\n";
    out << "Material Id: " << this->Id() << "\n";
    out << "Dimension: " << this->Dimension() << "\n\n";
}

void TPZMixedDarcyFlowHybrid::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
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
TPZMixedDarcyFlowHybrid::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const {
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
