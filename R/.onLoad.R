".onLoad" <-
function( libname, pkgname) {
  oa <- system.file( sprintf( 'R/load_%s_dll.R', pkgname),
      package=pkgname, lib.loc=libname)
  if( !nzchar( oa)) {
    oa <- system.file( sprintf( 'R/load_%s_dll.r', pkgname),
        package=pkgname, lib.loc=libname)
  }

  source( oa, local=TRUE)
} ## .onLoad transferred in from MVB2
