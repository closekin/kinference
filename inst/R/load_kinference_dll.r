


local({
  debug_dll <- getOption( 'kinference_debug_dll', '')

  loadenv <- if( 'kinference' %in% loadedNamespaces()) { 
      asNamespace( 'kinference') 
    } else {
      .GlobalEnv 
    }

  loadenv$source_DLL <- if( nzchar( debug_dll)) {
      dyn.load( debug_dll)
    } else {
      library.dynam( 'kinference', package='kinference', lib.loc=libname)
    }

    assign( "paircomps", Rcpp:::sourceCppFunction( function( pair_geno, LOD, geno1, geno2, symmo, granulum, granulum_loci) {}, FALSE, source_DLL, "C_paircomps"), loadenv)
    assign( "HSP_paircomps_lots", Rcpp:::sourceCppFunction( function( pair_geno, LOD, geno1, geno2, symmo, eta, min_keep_PLOD, bins) {}, FALSE, source_DLL, "C_HSP_paircomps_lots"), loadenv)
    assign( "POP_paircomps_lots", Rcpp:::sourceCppFunction( function( geno1, geno2, symmo, eta, max_keep_Nexclu, bins, AAO, BBO) {}, FALSE, source_DLL, "C_POP_paircomps_lots"), loadenv)
    assign( "POP_wt_paircomps_lots", Rcpp:::sourceCppFunction( function( geno1, geno2, w, symmo, eta, max_keep_wpsex, bins, AAO, BBO) {}, FALSE, source_DLL, "C_POP_wt_paircomps_lots"), loadenv)
    assign( "DUP_paircomps_lots", Rcpp:::sourceCppFunction( function( geno1, geno2, symmo, max_diff_genos) {}, FALSE, source_DLL, "C_DUP_paircomps_lots"), loadenv)
    assign( "DUP_paircomps_incomplete_lots", Rcpp:::sourceCppFunction( function( geno1, geno2, symmo, max_diff_ppn, limit) {}, FALSE, source_DLL, "C_DUP_paircomps_incomplete_lots"), loadenv)
    assign( "indiv_lglk_geno", Rcpp:::sourceCppFunction( function( lpgeno, geno) {}, FALSE, source_DLL, "C_indiv_lglk_geno"), loadenv)
    assign( "K_indiv", Rcpp:::sourceCppFunction( function( tt, geno, vec_LOD, Pg) {}, FALSE, source_DLL, "C_K_indiv"), loadenv)
})
