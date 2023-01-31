// Created: 2023-01-30 15:32:18
#include <RcppArmadillo.h>
using namespace Rcpp;

// paircomps
SEXP paircomps( IntegerMatrix pair_geno, NumericMatrix LOD, RawMatrix geno1, RawMatrix geno2, bool symmo, int granulum, int granulum_loci);
// RcppExport SEXP WRAP_kinference_paircomps( SEXP pair_genoSEXP, SEXP LODSEXP, SEXP geno1SEXP, SEXP geno2SEXP, SEXP symmoSEXP, SEXP granulumSEXP, SEXP granulum_lociSEXP) {
RcppExport SEXP _paircomps( SEXP pair_genoSEXP, SEXP LODSEXP, SEXP geno1SEXP, SEXP geno2SEXP, SEXP symmoSEXP, SEXP granulumSEXP, SEXP granulum_lociSEXP) {
BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix>::type pair_geno( pair_genoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix>::type LOD( LODSEXP);
    Rcpp::traits::input_parameter< RawMatrix>::type geno1( geno1SEXP);
    Rcpp::traits::input_parameter< RawMatrix>::type geno2( geno2SEXP);
    Rcpp::traits::input_parameter< bool>::type symmo( symmoSEXP);
    Rcpp::traits::input_parameter< int>::type granulum( granulumSEXP);
    Rcpp::traits::input_parameter< int>::type granulum_loci( granulum_lociSEXP);
  __result = Rcpp::wrap( paircomps( pair_geno, LOD, geno1, geno2, symmo, granulum, granulum_loci));
  return __result;
END_RCPP
}

// HSP_paircomps_lots
SEXP HSP_paircomps_lots( IntegerMatrix pair_geno, NumericMatrix LOD, RawMatrix geno1, RawMatrix geno2, bool symmo, double eta, double min_keep_PLOD, int keep_n, double minbin, double binterval, int nbins);
// RcppExport SEXP WRAP_kinference_HSP_paircomps_lots( SEXP pair_genoSEXP, SEXP LODSEXP, SEXP geno1SEXP, SEXP geno2SEXP, SEXP symmoSEXP, SEXP etaSEXP, SEXP min_keep_PLODSEXP, SEXP keep_nSEXP, SEXP minbinSEXP, SEXP bintervalSEXP, SEXP nbinsSEXP) {
RcppExport SEXP _HSP_paircomps_lots( SEXP pair_genoSEXP, SEXP LODSEXP, SEXP geno1SEXP, SEXP geno2SEXP, SEXP symmoSEXP, SEXP etaSEXP, SEXP min_keep_PLODSEXP, SEXP keep_nSEXP, SEXP minbinSEXP, SEXP bintervalSEXP, SEXP nbinsSEXP) {
BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix>::type pair_geno( pair_genoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix>::type LOD( LODSEXP);
    Rcpp::traits::input_parameter< RawMatrix>::type geno1( geno1SEXP);
    Rcpp::traits::input_parameter< RawMatrix>::type geno2( geno2SEXP);
    Rcpp::traits::input_parameter< bool>::type symmo( symmoSEXP);
    Rcpp::traits::input_parameter< double>::type eta( etaSEXP);
    Rcpp::traits::input_parameter< double>::type min_keep_PLOD( min_keep_PLODSEXP);
    Rcpp::traits::input_parameter< int>::type keep_n( keep_nSEXP);
    Rcpp::traits::input_parameter< double>::type minbin( minbinSEXP);
    Rcpp::traits::input_parameter< double>::type binterval( bintervalSEXP);
    Rcpp::traits::input_parameter< int>::type nbins( nbinsSEXP);
  __result = Rcpp::wrap( HSP_paircomps_lots( pair_geno, LOD, geno1, geno2, symmo, eta, min_keep_PLOD, keep_n, minbin, binterval, nbins));
  return __result;
END_RCPP
}

// POP_paircomps_lots
SEXP POP_paircomps_lots( RawMatrix geno1, RawMatrix geno2, bool symmo, double eta, int max_keep_Nexclu, int keep_n, NumericVector bins, int AAO, int BBO);
// RcppExport SEXP WRAP_kinference_POP_paircomps_lots( SEXP geno1SEXP, SEXP geno2SEXP, SEXP symmoSEXP, SEXP etaSEXP, SEXP max_keep_NexcluSEXP, SEXP keep_nSEXP, SEXP binsSEXP, SEXP AAOSEXP, SEXP BBOSEXP) {
RcppExport SEXP _POP_paircomps_lots( SEXP geno1SEXP, SEXP geno2SEXP, SEXP symmoSEXP, SEXP etaSEXP, SEXP max_keep_NexcluSEXP, SEXP keep_nSEXP, SEXP binsSEXP, SEXP AAOSEXP, SEXP BBOSEXP) {
BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< RawMatrix>::type geno1( geno1SEXP);
    Rcpp::traits::input_parameter< RawMatrix>::type geno2( geno2SEXP);
    Rcpp::traits::input_parameter< bool>::type symmo( symmoSEXP);
    Rcpp::traits::input_parameter< double>::type eta( etaSEXP);
    Rcpp::traits::input_parameter< int>::type max_keep_Nexclu( max_keep_NexcluSEXP);
    Rcpp::traits::input_parameter< int>::type keep_n( keep_nSEXP);
    Rcpp::traits::input_parameter< NumericVector>::type bins( binsSEXP);
    Rcpp::traits::input_parameter< int>::type AAO( AAOSEXP);
    Rcpp::traits::input_parameter< int>::type BBO( BBOSEXP);
  __result = Rcpp::wrap( POP_paircomps_lots( geno1, geno2, symmo, eta, max_keep_Nexclu, keep_n, bins, AAO, BBO));
  return __result;
END_RCPP
}

// POP_wt_paircomps_lots
SEXP POP_wt_paircomps_lots( RawMatrix geno1, RawMatrix geno2, NumericVector w, bool symmo, double eta, double max_keep_wpsex, int keep_n, int AAO, int BBO, int nbins, double binterval);
// RcppExport SEXP WRAP_kinference_POP_wt_paircomps_lots( SEXP geno1SEXP, SEXP geno2SEXP, SEXP wSEXP, SEXP symmoSEXP, SEXP etaSEXP, SEXP max_keep_wpsexSEXP, SEXP keep_nSEXP, SEXP AAOSEXP, SEXP BBOSEXP, SEXP nbinsSEXP, SEXP bintervalSEXP) {
RcppExport SEXP _POP_wt_paircomps_lots( SEXP geno1SEXP, SEXP geno2SEXP, SEXP wSEXP, SEXP symmoSEXP, SEXP etaSEXP, SEXP max_keep_wpsexSEXP, SEXP keep_nSEXP, SEXP AAOSEXP, SEXP BBOSEXP, SEXP nbinsSEXP, SEXP bintervalSEXP) {
BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< RawMatrix>::type geno1( geno1SEXP);
    Rcpp::traits::input_parameter< RawMatrix>::type geno2( geno2SEXP);
    Rcpp::traits::input_parameter< NumericVector>::type w( wSEXP);
    Rcpp::traits::input_parameter< bool>::type symmo( symmoSEXP);
    Rcpp::traits::input_parameter< double>::type eta( etaSEXP);
    Rcpp::traits::input_parameter< double>::type max_keep_wpsex( max_keep_wpsexSEXP);
    Rcpp::traits::input_parameter< int>::type keep_n( keep_nSEXP);
    Rcpp::traits::input_parameter< int>::type AAO( AAOSEXP);
    Rcpp::traits::input_parameter< int>::type BBO( BBOSEXP);
    Rcpp::traits::input_parameter< int>::type nbins( nbinsSEXP);
    Rcpp::traits::input_parameter< double>::type binterval( bintervalSEXP);
  __result = Rcpp::wrap( POP_wt_paircomps_lots( geno1, geno2, w, symmo, eta, max_keep_wpsex, keep_n, AAO, BBO, nbins, binterval));
  return __result;
END_RCPP
}

// DUP_paircomps_lots
SEXP DUP_paircomps_lots( RawMatrix geno1, RawMatrix geno2, bool symmo, double max_diff_loci, int keep_n, int nbins, double binterval, double maxbin);
// RcppExport SEXP WRAP_kinference_DUP_paircomps_lots( SEXP geno1SEXP, SEXP geno2SEXP, SEXP symmoSEXP, SEXP max_diff_lociSEXP, SEXP keep_nSEXP, SEXP nbinsSEXP, SEXP bintervalSEXP, SEXP maxbinSEXP) {
RcppExport SEXP _DUP_paircomps_lots( SEXP geno1SEXP, SEXP geno2SEXP, SEXP symmoSEXP, SEXP max_diff_lociSEXP, SEXP keep_nSEXP, SEXP nbinsSEXP, SEXP bintervalSEXP, SEXP maxbinSEXP) {
BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< RawMatrix>::type geno1( geno1SEXP);
    Rcpp::traits::input_parameter< RawMatrix>::type geno2( geno2SEXP);
    Rcpp::traits::input_parameter< bool>::type symmo( symmoSEXP);
    Rcpp::traits::input_parameter< double>::type max_diff_loci( max_diff_lociSEXP);
    Rcpp::traits::input_parameter< int>::type keep_n( keep_nSEXP);
    Rcpp::traits::input_parameter< int>::type nbins( nbinsSEXP);
    Rcpp::traits::input_parameter< double>::type binterval( bintervalSEXP);
    Rcpp::traits::input_parameter< double>::type maxbin( maxbinSEXP);
  __result = Rcpp::wrap( DUP_paircomps_lots( geno1, geno2, symmo, max_diff_loci, keep_n, nbins, binterval, maxbin));
  return __result;
END_RCPP
}

// DUP_paircomps_incomplete_lots
SEXP DUP_paircomps_incomplete_lots( RawMatrix geno1, RawMatrix geno2, bool symmo, double max_diff_ppn, int limit);
// RcppExport SEXP WRAP_kinference_DUP_paircomps_incomplete_lots( SEXP geno1SEXP, SEXP geno2SEXP, SEXP symmoSEXP, SEXP max_diff_ppnSEXP, SEXP limitSEXP) {
RcppExport SEXP _DUP_paircomps_incomplete_lots( SEXP geno1SEXP, SEXP geno2SEXP, SEXP symmoSEXP, SEXP max_diff_ppnSEXP, SEXP limitSEXP) {
BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< RawMatrix>::type geno1( geno1SEXP);
    Rcpp::traits::input_parameter< RawMatrix>::type geno2( geno2SEXP);
    Rcpp::traits::input_parameter< bool>::type symmo( symmoSEXP);
    Rcpp::traits::input_parameter< double>::type max_diff_ppn( max_diff_ppnSEXP);
    Rcpp::traits::input_parameter< int>::type limit( limitSEXP);
  __result = Rcpp::wrap( DUP_paircomps_incomplete_lots( geno1, geno2, symmo, max_diff_ppn, limit));
  return __result;
END_RCPP
}

// indiv_lglk_geno
SEXP indiv_lglk_geno( NumericMatrix lpgeno, RawMatrix geno);
// RcppExport SEXP WRAP_kinference_indiv_lglk_geno( SEXP lpgenoSEXP, SEXP genoSEXP) {
RcppExport SEXP _indiv_lglk_geno( SEXP lpgenoSEXP, SEXP genoSEXP) {
BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix>::type lpgeno( lpgenoSEXP);
    Rcpp::traits::input_parameter< RawMatrix>::type geno( genoSEXP);
  __result = Rcpp::wrap( indiv_lglk_geno( lpgeno, geno));
  return __result;
END_RCPP
}

// K_indiv
SEXP K_indiv( NumericVector tt, RawMatrix geno, NumericVector vec_LOD, NumericMatrix Pg);
// RcppExport SEXP WRAP_kinference_K_indiv( SEXP ttSEXP, SEXP genoSEXP, SEXP vec_LODSEXP, SEXP PgSEXP) {
RcppExport SEXP _K_indiv( SEXP ttSEXP, SEXP genoSEXP, SEXP vec_LODSEXP, SEXP PgSEXP) {
BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector>::type tt( ttSEXP);
    Rcpp::traits::input_parameter< RawMatrix>::type geno( genoSEXP);
    Rcpp::traits::input_parameter< NumericVector>::type vec_LOD( vec_LODSEXP);
    Rcpp::traits::input_parameter< NumericMatrix>::type Pg( PgSEXP);
  __result = Rcpp::wrap( K_indiv( tt, geno, vec_LOD, Pg));
  return __result;
END_RCPP
}


static const R_CallMethodDef CallEntries[] = {
    {"_kinference_paircomps", (DL_FUNC) &_paircomps, 7},
    {"_kinference_HSP_paircomps_lots", (DL_FUNC) &_HSP_paircomps_lots, 11},
    {"_kinference_POP_paircomps_lots", (DL_FUNC) &_POP_paircomps_lots, 9},
    {"_kinference_POP_wt_paircomps_lots", (DL_FUNC) &_POP_wt_paircomps_lots, 11},
    {"_kinference_DUP_paircomps_lots", (DL_FUNC) &_DUP_paircomps_lots, 8},
    {"_kinference_DUP_paircomps_incomplete_lots", (DL_FUNC) &_DUP_paircomps_incomplete_lots, 5},
    {"_kinference_indiv_lglk_geno", (DL_FUNC) &_indiv_lglk_geno, 2},
    {"_kinference_K_indiv", (DL_FUNC) &_K_indiv, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_kinference(DllInfo* info) {
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols( info, TRUE);
}
