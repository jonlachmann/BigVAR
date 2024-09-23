// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/BigVAR.h"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// VARXCons
MatrixXd VARXCons(NumericMatrix Y1, NumericMatrix X1, const int k, const int p, const int m, int s, bool oos, bool contemp);
RcppExport SEXP _BigVAR_VARXCons(SEXP Y1SEXP, SEXP X1SEXP, SEXP kSEXP, SEXP pSEXP, SEXP mSEXP, SEXP sSEXP, SEXP oosSEXP, SEXP contempSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< bool >::type oos(oosSEXP);
    Rcpp::traits::input_parameter< bool >::type contemp(contempSEXP);
    rcpp_result_gen = Rcpp::wrap(VARXCons(Y1, X1, k, p, m, s, oos, contemp));
    return rcpp_result_gen;
END_RCPP
}
// ARFitVARXR
List ARFitVARXR(NumericMatrix K21, const int k, const int p, int m, int s);
RcppExport SEXP _BigVAR_ARFitVARXR(SEXP K21SEXP, SEXP kSEXP, SEXP pSEXP, SEXP mSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type K21(K21SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(ARFitVARXR(K21, k, p, m, s));
    return rcpp_result_gen;
END_RCPP
}
// ICX
List ICX(NumericMatrix Y1, NumericMatrix X1, double k, int pmax, int smax, double m, std::string pen, int h);
RcppExport SEXP _BigVAR_ICX(SEXP Y1SEXP, SEXP X1SEXP, SEXP kSEXP, SEXP pmaxSEXP, SEXP smaxSEXP, SEXP mSEXP, SEXP penSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type pmax(pmaxSEXP);
    Rcpp::traits::input_parameter< int >::type smax(smaxSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< std::string >::type pen(penSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(ICX(Y1, X1, k, pmax, smax, m, pen, h));
    return rcpp_result_gen;
END_RCPP
}
// ST1a
double ST1a(double z, double gam);
RcppExport SEXP _BigVAR_ST1a(SEXP zSEXP, SEXP gamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    rcpp_result_gen = Rcpp::wrap(ST1a(z, gam));
    return rcpp_result_gen;
END_RCPP
}
// ST3a
colvec ST3a(colvec z, double gam);
RcppExport SEXP _BigVAR_ST3a(SEXP zSEXP, SEXP gamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< colvec >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    rcpp_result_gen = Rcpp::wrap(ST3a(z, gam));
    return rcpp_result_gen;
END_RCPP
}
// ST3ares
colvec ST3ares(colvec z, double gam, colvec restrictions);
RcppExport SEXP _BigVAR_ST3ares(SEXP zSEXP, SEXP gamSEXP, SEXP restrictionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< colvec >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< colvec >::type restrictions(restrictionsSEXP);
    rcpp_result_gen = Rcpp::wrap(ST3ares(z, gam, restrictions));
    return rcpp_result_gen;
END_RCPP
}
// gamloopFista
cube gamloopFista(NumericVector beta_, const mat& Y, const mat& Z, const mat gammgrid, const double eps, const colvec& YMean2, const colvec& ZMean2, mat& B1, int k, int p, double tk, int k1, int s, mat restrictions, bool sep_lambda);
RcppExport SEXP _BigVAR_gamloopFista(SEXP beta_SEXP, SEXP YSEXP, SEXP ZSEXP, SEXP gammgridSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP B1SEXP, SEXP kSEXP, SEXP pSEXP, SEXP tkSEXP, SEXP k1SEXP, SEXP sSEXP, SEXP restrictionsSEXP, SEXP sep_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const mat >::type gammgrid(gammgridSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< mat& >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type tk(tkSEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< mat >::type restrictions(restrictionsSEXP);
    Rcpp::traits::input_parameter< bool >::type sep_lambda(sep_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(gamloopFista(beta_, Y, Z, gammgrid, eps, YMean2, ZMean2, B1, k, p, tk, k1, s, restrictions, sep_lambda));
    return rcpp_result_gen;
END_RCPP
}
// gamloopFistaEN
cube gamloopFistaEN(NumericVector beta_, const mat& Y, const mat& Z, const mat gammgrid, vec& alpha, const double eps, const colvec& YMean2, const colvec& ZMean2, mat& B1, int k, int p, double tk, int k1, int s, mat restrictions, bool sep_lambda);
RcppExport SEXP _BigVAR_gamloopFistaEN(SEXP beta_SEXP, SEXP YSEXP, SEXP ZSEXP, SEXP gammgridSEXP, SEXP alphaSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP B1SEXP, SEXP kSEXP, SEXP pSEXP, SEXP tkSEXP, SEXP k1SEXP, SEXP sSEXP, SEXP restrictionsSEXP, SEXP sep_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const mat >::type gammgrid(gammgridSEXP);
    Rcpp::traits::input_parameter< vec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< mat& >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type tk(tkSEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< mat >::type restrictions(restrictionsSEXP);
    Rcpp::traits::input_parameter< bool >::type sep_lambda(sep_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(gamloopFistaEN(beta_, Y, Z, gammgrid, alpha, eps, YMean2, ZMean2, B1, k, p, tk, k1, s, restrictions, sep_lambda));
    return rcpp_result_gen;
END_RCPP
}
// Eigencomp
List Eigencomp(mat& Z1, List groups, int n1, int k1);
RcppExport SEXP _BigVAR_Eigencomp(SEXP Z1SEXP, SEXP groupsSEXP, SEXP n1SEXP, SEXP k1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    rcpp_result_gen = Rcpp::wrap(Eigencomp(Z1, groups, n1, k1));
    return rcpp_result_gen;
END_RCPP
}
// EigencompOO
List EigencompOO(mat& ZZ1, List groups, int n1, int k);
RcppExport SEXP _BigVAR_EigencompOO(SEXP ZZ1SEXP, SEXP groupsSEXP, SEXP n1SEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat& >::type ZZ1(ZZ1SEXP);
    Rcpp::traits::input_parameter< List >::type groups(groupsSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(EigencompOO(ZZ1, groups, n1, k));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopGL2
List GamLoopGL2(NumericVector beta_, List Activeset, NumericVector gamm, const mat& Y1, const mat& Z1, List jj, List jjfull, List jjcomp, double eps, const colvec& YMean2, const colvec& ZMean2, int k, int pk, const List M2f_, const List eigvalF_, const List eigvecF_, mat restrictions);
RcppExport SEXP _BigVAR_GamLoopGL2(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP Y1SEXP, SEXP Z1SEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigvalF_SEXP, SEXP eigvecF_SEXP, SEXP restrictionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< const List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< const List >::type eigvalF_(eigvalF_SEXP);
    Rcpp::traits::input_parameter< const List >::type eigvecF_(eigvecF_SEXP);
    Rcpp::traits::input_parameter< mat >::type restrictions(restrictionsSEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopGL2(beta_, Activeset, gamm, Y1, Z1, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M2f_, eigvalF_, eigvecF_, restrictions));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopGLOO
List GamLoopGLOO(NumericVector beta_, List Activeset, NumericVector gamm, const mat& Y, const mat& Z, List jj, List jjfull, List jjcomp, double eps, colvec& YMean2, colvec& ZMean2, int k, int pk, List M2f_, List eigvalF_, List eigvecF_, int k1);
RcppExport SEXP _BigVAR_GamLoopGLOO(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigvalF_SEXP, SEXP eigvecF_SEXP, SEXP k1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< List >::type eigvalF_(eigvalF_SEXP);
    Rcpp::traits::input_parameter< List >::type eigvecF_(eigvecF_SEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopGLOO(beta_, Activeset, gamm, Y, Z, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M2f_, eigvalF_, eigvecF_, k1));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGLOO
List GamLoopSGLOO(NumericVector beta_, const List Activeset_, const NumericVector gamm, const double alpha, const mat& Y, const mat& Z, List jj_, const List jjfull_, List jjcomp_, const double eps, const colvec& YMean2, const colvec& ZMean2, const int k1, const int pk, const List M2f_, const NumericVector eigs_, double m);
RcppExport SEXP _BigVAR_GamLoopSGLOO(SEXP beta_SEXP, SEXP Activeset_SEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP jj_SEXP, SEXP jjfull_SEXP, SEXP jjcomp_SEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP k1SEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigs_SEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const List >::type Activeset_(Activeset_SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< List >::type jj_(jj_SEXP);
    Rcpp::traits::input_parameter< const List >::type jjfull_(jjfull_SEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp_(jjcomp_SEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< const int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< const int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< const List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type eigs_(eigs_SEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGLOO(beta_, Activeset_, gamm, alpha, Y, Z, jj_, jjfull_, jjcomp_, eps, YMean2, ZMean2, k1, pk, M2f_, eigs_, m));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGLOODP
List GamLoopSGLOODP(NumericVector beta_, const List Activeset_, mat gamm, const colvec alpha, const mat& Y, const mat& Z, List jj_, const List jjfull_, List jjcomp_, const double eps, const colvec& YMean2, const colvec& ZMean2, const int k1, const int pk, const List M2f_, const NumericVector eigs_, double m);
RcppExport SEXP _BigVAR_GamLoopSGLOODP(SEXP beta_SEXP, SEXP Activeset_SEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP YSEXP, SEXP ZSEXP, SEXP jj_SEXP, SEXP jjfull_SEXP, SEXP jjcomp_SEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP k1SEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigs_SEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const List >::type Activeset_(Activeset_SEXP);
    Rcpp::traits::input_parameter< mat >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< List >::type jj_(jj_SEXP);
    Rcpp::traits::input_parameter< const List >::type jjfull_(jjfull_SEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp_(jjcomp_SEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< const int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< const int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< const List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type eigs_(eigs_SEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGLOODP(beta_, Activeset_, gamm, alpha, Y, Z, jj_, jjfull_, jjcomp_, eps, YMean2, ZMean2, k1, pk, M2f_, eigs_, m));
    return rcpp_result_gen;
END_RCPP
}
// Fistapar
mat Fistapar(const mat Y, const mat Z, const mat phi, const int L, const rowvec lambda, const double eps, const double tk, const int k, bool sep_lambda);
RcppExport SEXP _BigVAR_Fistapar(SEXP YSEXP, SEXP ZSEXP, SEXP phiSEXP, SEXP LSEXP, SEXP lambdaSEXP, SEXP epsSEXP, SEXP tkSEXP, SEXP kSEXP, SEXP sep_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const mat >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const int >::type L(LSEXP);
    Rcpp::traits::input_parameter< const rowvec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const double >::type tk(tkSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type sep_lambda(sep_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(Fistapar(Y, Z, phi, L, lambda, eps, tk, k, sep_lambda));
    return rcpp_result_gen;
END_RCPP
}
// gamloopHLAG
cube gamloopHLAG(NumericVector beta_, const mat& Y, const mat& Z, mat gammgrid, const double eps, const colvec& YMean2, const colvec& ZMean2, mat& B1, const int k, const int p, bool sep_lambda);
RcppExport SEXP _BigVAR_gamloopHLAG(SEXP beta_SEXP, SEXP YSEXP, SEXP ZSEXP, SEXP gammgridSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP B1SEXP, SEXP kSEXP, SEXP pSEXP, SEXP sep_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< mat >::type gammgrid(gammgridSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< mat& >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type sep_lambda(sep_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(gamloopHLAG(beta_, Y, Z, gammgrid, eps, YMean2, ZMean2, B1, k, p, sep_lambda));
    return rcpp_result_gen;
END_RCPP
}
// gamloopOO
cube gamloopOO(NumericVector beta_, const mat Y, const mat Z, mat gammgrid, const double eps, const colvec YMean2, const colvec ZMean2, mat B1, const int k, const int p, colvec w, List groups_, bool sep_lambda);
RcppExport SEXP _BigVAR_gamloopOO(SEXP beta_SEXP, SEXP YSEXP, SEXP ZSEXP, SEXP gammgridSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP B1SEXP, SEXP kSEXP, SEXP pSEXP, SEXP wSEXP, SEXP groups_SEXP, SEXP sep_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< mat >::type gammgrid(gammgridSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< mat >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< colvec >::type w(wSEXP);
    Rcpp::traits::input_parameter< List >::type groups_(groups_SEXP);
    Rcpp::traits::input_parameter< bool >::type sep_lambda(sep_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(gamloopOO(beta_, Y, Z, gammgrid, eps, YMean2, ZMean2, B1, k, p, w, groups_, sep_lambda));
    return rcpp_result_gen;
END_RCPP
}
// FistaElem
mat FistaElem(const mat& Y, const mat& Z, mat phi, const int p, const int k, rowvec lambda, const double eps, const double tk, bool sep_lambda);
RcppExport SEXP _BigVAR_FistaElem(SEXP YSEXP, SEXP ZSEXP, SEXP phiSEXP, SEXP pSEXP, SEXP kSEXP, SEXP lambdaSEXP, SEXP epsSEXP, SEXP tkSEXP, SEXP sep_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< mat >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< rowvec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const double >::type tk(tkSEXP);
    Rcpp::traits::input_parameter< bool >::type sep_lambda(sep_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(FistaElem(Y, Z, phi, p, k, lambda, eps, tk, sep_lambda));
    return rcpp_result_gen;
END_RCPP
}
// gamloopElem
cube gamloopElem(NumericVector beta_, const mat& Y, const mat& Z, mat gammgrid, const double eps, const colvec YMean2, const colvec ZMean2, mat B1, const int k, const int p, bool sep_lambda);
RcppExport SEXP _BigVAR_gamloopElem(SEXP beta_SEXP, SEXP YSEXP, SEXP ZSEXP, SEXP gammgridSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP B1SEXP, SEXP kSEXP, SEXP pSEXP, SEXP sep_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< mat >::type gammgrid(gammgridSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< mat >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type sep_lambda(sep_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(gamloopElem(beta_, Y, Z, gammgrid, eps, YMean2, ZMean2, B1, k, p, sep_lambda));
    return rcpp_result_gen;
END_RCPP
}
// powermethod
List powermethod(mat A, colvec x1);
RcppExport SEXP _BigVAR_powermethod(SEXP ASEXP, SEXP x1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< colvec >::type x1(x1SEXP);
    rcpp_result_gen = Rcpp::wrap(powermethod(A, x1));
    return rcpp_result_gen;
END_RCPP
}
// norm2
double norm2(NumericVector x);
RcppExport SEXP _BigVAR_norm2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(norm2(x));
    return rcpp_result_gen;
END_RCPP
}
// RelaxedLS
mat RelaxedLS(const mat K, mat B2);
RcppExport SEXP _BigVAR_RelaxedLS(SEXP KSEXP, SEXP B2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< mat >::type B2(B2SEXP);
    rcpp_result_gen = Rcpp::wrap(RelaxedLS(K, B2));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGLX
List GamLoopSGLX(NumericVector beta_, List Activeset, NumericVector gamm, double alpha, const mat& Y1, const mat& Z1, List jj, List jjfull, List jjcomp, double eps, colvec YMean2, colvec ZMean2, int k, int pk, List M2f_, NumericVector eigs, int k1);
RcppExport SEXP _BigVAR_GamLoopSGLX(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP Y1SEXP, SEXP Z1SEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigsSEXP, SEXP k1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eigs(eigsSEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGLX(beta_, Activeset, gamm, alpha, Y1, Z1, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M2f_, eigs, k1));
    return rcpp_result_gen;
END_RCPP
}
// proxvx2
colvec proxvx2(colvec v2, int L, double lambda, int m, int k, int F1);
RcppExport SEXP _BigVAR_proxvx2(SEXP v2SEXP, SEXP LSEXP, SEXP lambdaSEXP, SEXP mSEXP, SEXP kSEXP, SEXP F1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< colvec >::type v2(v2SEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type F1(F1SEXP);
    rcpp_result_gen = Rcpp::wrap(proxvx2(v2, L, lambda, m, k, F1));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGL
List GamLoopSGL(NumericVector beta_, List Activeset, const NumericVector gamm, const double alpha, const mat& Y1, const mat& Z1, List jj, const List jjfull, const List jjcomp, const double eps, const colvec YMean2, const colvec ZMean2, const int k, const int pk, const List M1f_, const List M2f_, const NumericVector eigs_);
RcppExport SEXP _BigVAR_GamLoopSGL(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP Y1SEXP, SEXP Z1SEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M1f_SEXP, SEXP M2f_SEXP, SEXP eigs_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< const List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< const List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< const List >::type M1f_(M1f_SEXP);
    Rcpp::traits::input_parameter< const List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type eigs_(eigs_SEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGL(beta_, Activeset, gamm, alpha, Y1, Z1, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M1f_, M2f_, eigs_));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGLDP
List GamLoopSGLDP(NumericVector beta_, List Activeset, const mat gamm, const colvec alpha, const mat& Y1, const mat& Z1, List jj, const List jjfull, const List jjcomp, const double eps, const colvec YMean2, const colvec ZMean2, const int k, const int pk, const List M1f_, const List M2f_, const NumericVector eigs_);
RcppExport SEXP _BigVAR_GamLoopSGLDP(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP Y1SEXP, SEXP Z1SEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M1f_SEXP, SEXP M2f_SEXP, SEXP eigs_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< const mat >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< const colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< const List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< const List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< const List >::type M1f_(M1f_SEXP);
    Rcpp::traits::input_parameter< const List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type eigs_(eigs_SEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGLDP(beta_, Activeset, gamm, alpha, Y1, Z1, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M1f_, M2f_, eigs_));
    return rcpp_result_gen;
END_RCPP
}
// GamLoopSGLXDP
List GamLoopSGLXDP(NumericVector beta_, List Activeset, mat gamm, colvec alpha, const mat& Y1, const mat& Z1, List jj, List jjfull, List jjcomp, double eps, colvec YMean2, colvec ZMean2, int k, int pk, List M2f_, NumericVector eigs, int k1);
RcppExport SEXP _BigVAR_GamLoopSGLXDP(SEXP beta_SEXP, SEXP ActivesetSEXP, SEXP gammSEXP, SEXP alphaSEXP, SEXP Y1SEXP, SEXP Z1SEXP, SEXP jjSEXP, SEXP jjfullSEXP, SEXP jjcompSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP kSEXP, SEXP pkSEXP, SEXP M2f_SEXP, SEXP eigsSEXP, SEXP k1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< List >::type Activeset(ActivesetSEXP);
    Rcpp::traits::input_parameter< mat >::type gamm(gammSEXP);
    Rcpp::traits::input_parameter< colvec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y1(Y1SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< List >::type jj(jjSEXP);
    Rcpp::traits::input_parameter< List >::type jjfull(jjfullSEXP);
    Rcpp::traits::input_parameter< List >::type jjcomp(jjcompSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< colvec >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< colvec >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type pk(pkSEXP);
    Rcpp::traits::input_parameter< List >::type M2f_(M2f_SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eigs(eigsSEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    rcpp_result_gen = Rcpp::wrap(GamLoopSGLXDP(beta_, Activeset, gamm, alpha, Y1, Z1, jj, jjfull, jjcomp, eps, YMean2, ZMean2, k, pk, M2f_, eigs, k1));
    return rcpp_result_gen;
END_RCPP
}
// mcp_loop
cube mcp_loop(mat Y, mat Z, cube B, const vec lambda, const double tol, double gamma, bool mcp);
RcppExport SEXP _BigVAR_mcp_loop(SEXP YSEXP, SEXP ZSEXP, SEXP BSEXP, SEXP lambdaSEXP, SEXP tolSEXP, SEXP gammaSEXP, SEXP mcpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< cube >::type B(BSEXP);
    Rcpp::traits::input_parameter< const vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type mcp(mcpSEXP);
    rcpp_result_gen = Rcpp::wrap(mcp_loop(Y, Z, B, lambda, tol, gamma, mcp));
    return rcpp_result_gen;
END_RCPP
}
// gamloopMCP
cube gamloopMCP(NumericVector beta_, const mat& Y, const mat& Z, vec lambda, const double eps, const colvec& YMean2, const colvec& ZMean2, double gamma, bool mcp);
RcppExport SEXP _BigVAR_gamloopMCP(SEXP beta_SEXP, SEXP YSEXP, SEXP ZSEXP, SEXP lambdaSEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP gammaSEXP, SEXP mcpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type mcp(mcpSEXP);
    rcpp_result_gen = Rcpp::wrap(gamloopMCP(beta_, Y, Z, lambda, eps, YMean2, ZMean2, gamma, mcp));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BigVAR_VARXCons", (DL_FUNC) &_BigVAR_VARXCons, 8},
    {"_BigVAR_ARFitVARXR", (DL_FUNC) &_BigVAR_ARFitVARXR, 5},
    {"_BigVAR_ICX", (DL_FUNC) &_BigVAR_ICX, 8},
    {"_BigVAR_ST1a", (DL_FUNC) &_BigVAR_ST1a, 2},
    {"_BigVAR_ST3a", (DL_FUNC) &_BigVAR_ST3a, 2},
    {"_BigVAR_ST3ares", (DL_FUNC) &_BigVAR_ST3ares, 3},
    {"_BigVAR_gamloopFista", (DL_FUNC) &_BigVAR_gamloopFista, 15},
    {"_BigVAR_gamloopFistaEN", (DL_FUNC) &_BigVAR_gamloopFistaEN, 16},
    {"_BigVAR_Eigencomp", (DL_FUNC) &_BigVAR_Eigencomp, 4},
    {"_BigVAR_EigencompOO", (DL_FUNC) &_BigVAR_EigencompOO, 4},
    {"_BigVAR_GamLoopGL2", (DL_FUNC) &_BigVAR_GamLoopGL2, 17},
    {"_BigVAR_GamLoopGLOO", (DL_FUNC) &_BigVAR_GamLoopGLOO, 17},
    {"_BigVAR_GamLoopSGLOO", (DL_FUNC) &_BigVAR_GamLoopSGLOO, 17},
    {"_BigVAR_GamLoopSGLOODP", (DL_FUNC) &_BigVAR_GamLoopSGLOODP, 17},
    {"_BigVAR_Fistapar", (DL_FUNC) &_BigVAR_Fistapar, 9},
    {"_BigVAR_gamloopHLAG", (DL_FUNC) &_BigVAR_gamloopHLAG, 11},
    {"_BigVAR_gamloopOO", (DL_FUNC) &_BigVAR_gamloopOO, 13},
    {"_BigVAR_FistaElem", (DL_FUNC) &_BigVAR_FistaElem, 9},
    {"_BigVAR_gamloopElem", (DL_FUNC) &_BigVAR_gamloopElem, 11},
    {"_BigVAR_powermethod", (DL_FUNC) &_BigVAR_powermethod, 2},
    {"_BigVAR_norm2", (DL_FUNC) &_BigVAR_norm2, 1},
    {"_BigVAR_RelaxedLS", (DL_FUNC) &_BigVAR_RelaxedLS, 2},
    {"_BigVAR_GamLoopSGLX", (DL_FUNC) &_BigVAR_GamLoopSGLX, 17},
    {"_BigVAR_proxvx2", (DL_FUNC) &_BigVAR_proxvx2, 6},
    {"_BigVAR_GamLoopSGL", (DL_FUNC) &_BigVAR_GamLoopSGL, 17},
    {"_BigVAR_GamLoopSGLDP", (DL_FUNC) &_BigVAR_GamLoopSGLDP, 17},
    {"_BigVAR_GamLoopSGLXDP", (DL_FUNC) &_BigVAR_GamLoopSGLXDP, 17},
    {"_BigVAR_mcp_loop", (DL_FUNC) &_BigVAR_mcp_loop, 7},
    {"_BigVAR_gamloopMCP", (DL_FUNC) &_BigVAR_gamloopMCP, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_BigVAR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
