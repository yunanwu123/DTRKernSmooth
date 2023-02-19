// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// obj_value_C
double obj_value_C(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd a, const Eigen::VectorXd& eta, const Eigen::VectorXd& prob);
RcppExport SEXP _DTRKernSmooth_obj_value_C(SEXP XSEXP, SEXP ySEXP, SEXP aSEXP, SEXP etaSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(obj_value_C(X, y, a, eta, prob));
    return rcpp_result_gen;
END_RCPP
}
// Smooth_C
List Smooth_C(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd a, Eigen::VectorXd initial, const Eigen::VectorXd& prob, const char* kn, double phi, const double gamma, const double err_tol, const int iter_tol);
RcppExport SEXP _DTRKernSmooth_Smooth_C(SEXP XSEXP, SEXP ySEXP, SEXP aSEXP, SEXP initialSEXP, SEXP probSEXP, SEXP knSEXP, SEXP phiSEXP, SEXP gammaSEXP, SEXP err_tolSEXP, SEXP iter_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< const char* >::type kn(knSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const double >::type err_tol(err_tolSEXP);
    Rcpp::traits::input_parameter< const int >::type iter_tol(iter_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(Smooth_C(X, y, a, initial, prob, kn, phi, gamma, err_tol, iter_tol));
    return rcpp_result_gen;
END_RCPP
}
// Boots_C
List Boots_C(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd a, Eigen::VectorXd initial, const Eigen::VectorXd& prob, const char* kn, Eigen::MatrixXd weights, const double alpha, double phi, const double gamma, const double err_tol, const int iter_tol);
RcppExport SEXP _DTRKernSmooth_Boots_C(SEXP XSEXP, SEXP ySEXP, SEXP aSEXP, SEXP initialSEXP, SEXP probSEXP, SEXP knSEXP, SEXP weightsSEXP, SEXP alphaSEXP, SEXP phiSEXP, SEXP gammaSEXP, SEXP err_tolSEXP, SEXP iter_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< const char* >::type kn(knSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const double >::type err_tol(err_tolSEXP);
    Rcpp::traits::input_parameter< const int >::type iter_tol(iter_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(Boots_C(X, y, a, initial, prob, kn, weights, alpha, phi, gamma, err_tol, iter_tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_DTRKernSmooth_obj_value_C", (DL_FUNC) &_DTRKernSmooth_obj_value_C, 5},
    {"_DTRKernSmooth_Smooth_C", (DL_FUNC) &_DTRKernSmooth_Smooth_C, 10},
    {"_DTRKernSmooth_Boots_C", (DL_FUNC) &_DTRKernSmooth_Boots_C, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_DTRKernSmooth(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}