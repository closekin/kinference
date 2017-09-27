map6to4 <- function(g6p0, g6p1){

  define_genotypes()

  # For the 4-ways, must condense g6p's

  map6to4 <- matrix( 0, 6, 4, dimnames=list( genotypes6, genotypes4_ambig))
  # AB & OO are OK; AAO should receive both AA and AO; etc
  mm <- match( genotypes6, substring( genotypes4_ambig, 1, 2), 0) # the "AA" bit of "AAO"...
  yup <- cbind( which(mm>0), mm[ mm>0])
  map6to4[ yup] <- 1
  mm <- match( genotypes6, substring( genotypes4_ambig, 2, 3), 0) # ... and the "AO" bit
  yup <- cbind( which(mm>0), mm[ mm>0])
  map6to4[ yup] <- 1

  # Really want g4p0[l,i,j] := map6to4[i,k6] %[k6]% g6p0[l,k6,m6] %[m6]% map6to4[m6,j]
  # ... but vecless can't presently handle multi-stages

  A[l,i,k] := g6p0[l,i,j] %[j]% map6to4[j,k]
  g4p0[l,i,j] := map6to4[ k,i] %[k]% A[l,k,j]

  A[l,i,k] := g6p1[l,i,j] %[j]% map6to4[j,k]
  g4p1[l,i,j] := map6to4[ k,i] %[k]% A[l,k,j]

  return(list(g4p0=g4p0, g4p1=g4p1))

}
