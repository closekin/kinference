"[.playback" <-
function( x, ...) {
  where <- x@where
  who <- x@whoami
  which <- where$counter[[ who]]
  where$counter[[ who]] <- which+1
  index <- where$subs[[ who]][[ which]]
  acall <- call( '[', quote( x)) # quote( x[]) returns length-3 with 3rd arg missing
  lpf <- length( where$prefixdims)
  for( ipref in 1 %upto% lpf) {
    acall[[ length( acall)+1]] <- alist(y=)$y # missing
  }
  acall[[ length( acall)+1]] <- as.vector( index) # should be vec anyway; this strips dim
  oldClass( x) <- NULL
  if( lpf) { # condense all trailing dimensions
    dim( x) <- c( where$prefixdims, prod( dim( x)[-(1:lpf)]))
  } # if not lpf, vector subscripting works fine for arrays
  res <- eval( acall)
  # dimnames( res) <- index@dimnames
res
}


"[.selfrecording_array" <-
function( x, ...) {
  where <- x@where
  whoami <- x@whoami
  if( is.null( odimx <- dim( x))) {
    odimx <- length( x)
  }
  si <- x <- unclass( x)
  si[] <- seq_along( si)
  mc <- match.call()
  si@where <- NULL
  mc[[1]] <- `[`
  mc[[2]] <- si

  res <- NextMethod('[')

  index <- unname( c( eval.parent( mc)))
  index@dim <- dim( res)
  if( !length( where$subs[[ whoami]])) { # record original dim(x), used by make_playback
    index@orig_full_dim <- odimx
  }

  # index@dimnames <- dimnames( res)
  where$subs[[ whoami]] <- c( where$subs[[ whoami]], list( index))
res
}


"[<-.playback" <-
function( x, ..., value) {
  where <- x@where
  who <- x@whoami
  which <- where$counter[[ who]]
  where$counter[[ who]] <- which+1
  index <- as.vector( where$subs[[ who]][[ which]])
  acall <- call( '[<-', quote( x))
  lpf <- length( where$prefixdims)
  for( ipref in 1 %upto% lpf) {
    acall[[ length( acall)+1]] <- alist(y=)$y # missing
  }
  acall[[ length( acall)+1]] <- index
  acall[[ length( acall)+1]] <- value
  oldClass( x) <- NULL
  if( lpf) { # condense all trailing dimensions
    odimx <- dim( x)
    dim( x) <- c( where$prefixdims, prod( dim( x)[-(1:lpf)]))
  } # if not lpf, vector subscripting works fine for arrays

  x[] <- eval( acall)
  if( lpf) {
    dim( x) <- odimx
  }
  oldClass( x) <- 'playback'
x
}


#' @importFrom mvbutils %is.a%
#' @importFrom atease @ @<-
"==.snpgeno" <- function(e1, e2) {
  if(e1 %is.a% 'snpgeno'){
    if(is.character(e2)){
      e2 <- as.raw(match( e2, e1@diplos, 0))
    }
    return(unclass(e1)==e2)
  }else{
    if(is.character(e1)) {
      e1 <- as.raw(match(e1, e1@diplos, 0))
    }
  }
  return( unclass(e2)==e1)
}


#' @importFrom mvbutils %without.name% ?
"add_list_defaults" <- function( l, ...) {
###### Add defaults to list 'l' if not already in 'l'
  defaults <- list(...)
  l <- c( l, defaults %without.name% names( l))
return( l)
}

#' @importFrom mvbutils named mlocal ? FOR
"add_pairprob_error" <- function( nlocal=sys.parent()) mlocal({
## Needs pp_true and snerr

  g_err <- named( genotypes_C) ? 0 # recorded because it makes perr easy to define
  perr <- 0*nchar( g_err) ? 1 # named

  g_1 <- substring( genotypes_C, 1, 1)
  g_2 <- substring( genotypes_C, 2, 2)
  AA_BB <- g_1==g_2 & (g_1==A | g_1==B)
  g_err[ AA_BB] <- g_1[ AA_BB] %&% 'O'
  perr[ AA_BB] <- snerr[ sprintf( '%1$s%1$s2%1$sO', g_1[ AA_BB])] ? 0 # XX2XO

  AO_BO <- (g_1==A | g_1==B) & g_2==O
  g_err[ AO_BO] <- g_1[ AO_BO] %&% g_1[ AO_BO]
  perr[ AO_BO] <- snerr[ sprintf( '%1$sO2%1$s%1$s', g_1[ AO_BO])] ? 0 # XO2XX

  # Original version: not sure it's selfrecordable
#  make_pp_err <- function( g1, g2) {
#      g1_err <- g_err[ g1]
#      g2_err <- g_err[ g2]
#      pp_true[ cbind( g1, g2) ] * (1-perr[ g1]) * (1-perr[ g2]) +
#      pp_true[ cbind( g1_err, g2) ] * perr[ g1_err] * (1-perr[ g2]) +
#      pp_true[ cbind( g1, g2_err) ] * (1-perr[ g1]) * perr[ g2_err] +
#      pp_true[ cbind( g1_err, g2_err) ] * perr[ g1_err] * perr[ g2_err]
#    }
#
  # pp_err[] <- outer( named( genotypes_C), named( genotypes_C), make_pp_err)

  # Need 2-stage setup here so selfrec works: "pure" assignment can't be overloaded in R
  pp_err <- matrix( -1, length( genotypes_C), length( genotypes_C), dimnames=rep( list( genotypes_C), 2)) ? 1

  gg <- list() ? 0 # this is a trick to allow with(gg,...) in playback, since the "real" gg isn't kept

  gg <- expand.grid( g1=genotypes_C, g2=genotypes_C, stringsAsFactors=FALSE)
  gg <- within( gg, {
     g1_err <- g_err[ g1]
     g2_err <- g_err[ g2]
  })
  extract.named( gg)

  pp_err[ cbind( g1, g2)] <-
      pp_true[ cbind( g1, g2) ] * (1-perr[ g1]) * (1-perr[ g2]) +
      pp_true[ cbind( g1_err, g2) ] * perr[ g1_err] * (1-perr[ g2]) +
      pp_true[ cbind( g1, g2_err) ] * (1-perr[ g1]) * perr[ g2_err] +
      pp_true[ cbind( g1_err, g2_err) ] * perr[ g1_err] * perr[ g2_err]        ? 0


  # Merge C with O. Has to be mlocal() to avoid <<- which buggers playback. Things before 'nlocal' are args; things after are temporaries
  add_up <- function( i, result, ..., nlocal=sys.parent(), dimor, other)  mlocal({
      result <- match( result, genotypes_C)
      dimor <- slice.index( pp_err, i)

      for( other in FOR( list( ...), match( ., genotypes_C))) {
        pp_err[ dimor==result] <- pp_err[ dimor==result] + pp_err[ dimor==other] ? 0
        pp_err[ dimor==other] <- 0 # avoid double-use ? 0
      }
    })

  for( i in 1:2) {
    add_up( i, OO, CC, CO)
    add_up( i, AO, AC)
    add_up( i, BO, BC)
  }

# IE:
# pp_err2[,OO] <- pp_err2[,CC] + pp_err2[,CO] + pp_err2[,OO]
# pp_err2[,CC] <- pp_err2[,CO] <- 0
# pp_err2[AO,] <- pp_err2[AC,] + pp_err2[AO,]
# pp_err2[AC,] <- 0

  pp6_err <- pp_err[ genotypes6, genotypes6] ? 0
})

#' @importFrom atease @ @<-
"as.character.snpgeno" <- function( x, ...) {
  y <- matrix( '', nrow( x), ncol( x))
  y[] <- x@diplos[ as.integer( x)]
  dimnames( y) <- list( rownames( x@info), x@locinfo$Locus)
return( y)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils FOR
"calc_g6probs" <-
function( pA, pB, pC, snerr) {
  # snerr = P( misclassifying true XX as XO, and vice versa)
  define_genotypes()
  p0really <- 1 - (pA+pB+pC)
  pO <- p0really + pC

  pAA <- pA^2
  pAB <- 2*pA*pB
  pBB <- pB^2
  pAC <- 2*pA*pC
  pBC <- 2*pB*pC

  pA0 <- 2*pA*p0really # really true null called 0
  pB0 <- 2*pB*p0really
  pOO <- pO^2

  # Calcs assume matrix (ie Xtuple loci); coerce if not, then undim later
  if( no_dim <- is.null( dim( snerr))) {
    snerr <- matrix( snerr, 1, ncol=length( snerr), dimnames=list( NULL, names( snerr)))
  }

  # Error possible iff (i) true AA or (ii) A present C not

  pAO <- pAC + pA0 * (1-snerr[,'AO2AA']) + pAA * snerr[,'AA2AO']
  pAA <- pAA * (1-snerr[,'AA2AO']) + pA0 * snerr[,'AO2AA']

  pBO <- pBC + pB0 * (1-snerr[,'BO2BB']) + pBB * snerr[,'BB2BO']
  pBB <- pBB * (1-snerr[,'BB2BO']) + pB0 * snerr[,'BO2BB']

  p6 <- do.call( 'cbind', FOR( genotypes6, get( 'p' %&% .)))
  perr <- cbind(
      AAO=pA0 * snerr[,'AO2AA'] + pAA * snerr[,'AA2AO'],
      BBO=pB0 * snerr[,'BO2BB'] + pBB * snerr[,'BB2BO'])
  colnames( p6) <- genotypes6
  if( no_dim) {
    p6 <- p6[1,]
    perr <- perr[1,]
  }

  p6@perr <- perr
return( p6)
}


#' @importFrom mvbutils cq ?
"calc_g6probs_IBD0_scalar" <-
function( P, snerr, record=FALSE) {
## SCALAR-ONLY VERSION... this is hard enough!
## Though can be called with 1-row matrix args, eg with( x@locinfo[1,], calc_g6probs_IBD1( pbonzer, snerr))

  # snerr = P( misclassifying true XX as XO, and vice versa)

  set_recording( cq( P, snerr, pp_true, pr2, pp_err, perr, pp6_err), record)
  define_genotypes()           ? 0

  P <- drop( P) # for scalar version
  P <- P              ? 0
stopifnot( my.all.equal( names( P), names( ABCO)))

  snerr <- drop( snerr) # for scalar version

  snerr <- snerr            ? 0
  pp_true <- matrix( 0, length( genotypes_C), length( genotypes_C), dimnames=list( genotypes_C, genotypes_C))         ? 1


  g_1 <- substring( genotypes_C, 1, 1)
  g_2 <- substring( genotypes_C, 2, 2)
  pr2 <-  nchar( named( genotypes_C))    ? 1 # named

  is_het <- g_1 != g_2
  pr2[] <- P[ g_1] * P[ g_2]             ? 0
  pr2[ is_het] <- 2 * pr2[ is_het]       ? 0

  pp_true <- matrix( 0, length( genotypes_C), length( genotypes_C),
      dimnames=rep( list( genotypes_C), 2))        ? 1

  extract.named( expand.grid( gp1=genotypes_C, gp2=genotypes_C, stringsAsFactors=FALSE))
  pp_true[ cbind( gp1, gp2)] <- pr2[ gp1] * pr2[ gp2]                       ? 0

  # NB that for this UP case, XX/XO errors shouldn't change the overall probs because the cutoffs are chosen to do exactly that!
  add_pairprob_error()

return( pp6_err)
}


#' @importFrom mvbutils cq ?
"calc_g6probs_IBD1_scalar" <-
function( P, snerr, record=FALSE) {
## SCALAR-ONLY VERSION... this is hard enough!
## Though can be called with 1-row matrix args, eg with( x@locinfo[1,], calc_g6probs_IBD1( pbonzer, snerr))

  set_recording( cq( P, snerr, pp_true, pp3, pp_err, perr, pp6_err), record)
  define_genotypes()           ? 0

  P <- drop( P) # for scalar version
  P <- P               ? 0
stopifnot( my.all.equal( names( P), names( ABCO)))

  snerr <- drop( snerr) # for scalar version
  snerr <- snerr        ? 0
  pp_true <- matrix( 0, length( genotypes_C), length( genotypes_C), dimnames=list( genotypes_C, genotypes_C))          ? 1

  # Note: AB/AB can share either A or B... so have to accumulate

  # always write as S,U1 and S,U2 so probs are unambig
  # then accumulate into sorted versions since *several* S,U1,U2 could contribute to one XY/ZW
  # ... or *none*, of course, if there is no shared allele in XY/ZW

  # 3-col DF of (shared, 1st unshared, 2nd unshared) alleles in the pair
  su1u2 <- as.matrix( expand.grid( ABCO, ABCO, ABCO))
  implied <- su1u2[ ,c(1,2,1,3)]
  swap12 <- implied[,2] < implied[,1]
  implied[ swap12, 1:2] <- implied[ swap12, 2:1]
  swap34 <- implied[,4] < implied[,3]
  implied[ swap34, 3:4] <- implied[ swap34, 4:3]
  cimplied <- cbind( implied[,1] %&% implied[,2], implied[,3] %&% implied[,4])

  pp3 <- P[ su1u2[,1]] * P[ su1u2[,2]] * P[ su1u2[,3]]             ? 0

  # Each row can contribute to at most 2 implied XY/ZW combinations
  # eg AB/AB <- ABB or BAA
  # AB/CD <- nothing
  # AB/AC <- ABC

  implstr <- cimplied[,1] %&% cimplied[,2]
  has_impl2 <- which( duplicated( implstr))

  first_impl <- match( implstr[ -has_impl2], implstr)
  second_impl <- match( implstr[has_impl2], implstr, 0)
  cimpl2 <- cimplied[ second_impl,]

  pp_true[ cimplied[ first_impl,]] <- pp3[ -has_impl2]            ? 0
  pp_true[ cimpl2] <- pp_true[ cimpl2] + pp3[ has_impl2]          ? 0

  # Allow for XX <-> XO errors--- hopefully the only ones! (Watch out for scaffoldy version)
  # Assumes errors INDEPENDENT even tho there could be heritability (eg weak grabbing of mutated primer)
  # CC <-> CO is ignored in this version, since C gets scored as O

  add_pairprob_error()

return( pp6_err)
}


#' @importFrom mvbutils cq extract.named named ?
"cg" <- function( P, snerr, record=FALSE) {
## SCALAR-ONLY VERSION... this is hard enough!
## Though can be called with 1-row matrix args, eg with( x@locinfo[1,], calc_g6probs_IBD1( pbonzer, snerr))

  set_recording( cq( P, snerr, pp_true, pp3, pp_err, perr, pp6_err), record)
  define_genotypes() ? 0

  # Could include the next lines in previous recordar call, since grouping boring things inside {} is OK...
  # ... but there's no need!
  ABCO <- named( cq( A, B, C, O))
  extract.named( ABCO) # A, B, C, and O

  P <- drop( P) # for scalar version
  P <- P ? 0
stopifnot( my.all.equal( names( P), names( ABCO)))

  snerr <- drop( snerr) # for scalar version
  snerr <- snerr ? 0
  pp_true <- matrix( 0, length( genotypes_C), length( genotypes_C), dimnames=list( genotypes_C, genotypes_C)) ? 1

  # Note: AB/AB can share either A or B... so have to accumulate

  # always write as S,U1 and S,U2 so probs are unambig
  # then accumulate into sorted versions since *several* S,U1,U2 could contribute to one XY/ZW
  # ... or *none*, of course, if there is no shared allele in XY/ZW

  # 3-col DF of (shared, 1st unshared, 2nd unshared) alleles in the pair
  su1u2 <- as.matrix( expand.grid( ABCO, ABCO, ABCO))
  implied <- su1u2[ ,c(1,2,1,3)]
  swap12 <- implied[,2] < implied[,1]
  implied[ swap12, 1:2] <- implied[ swap12, 2:1]
  swap34 <- implied[,4] < implied[,3]
  implied[ swap34, 3:4] <- implied[ swap34, 4:3]
  cimplied <- cbind( implied[,1] %&% implied[,2], implied[,3] %&% implied[,4])

  pp3 <- P[ su1u2[,1]] * P[ su1u2[,2]] * P[ su1u2[,3]] ? 0

  # Each row can contribute to at most 2 implied XY/ZW combinations
  # eg AB/AB <- ABB or BAA
  # AB/CD <- nothing
  # AB/AC <- ABC

  implstr <- cimplied[,1] %&% cimplied[,2]
  has_impl2 <- which( duplicated( implstr))

  first_impl <- match( implstr[ -has_impl2], implstr)
  second_impl <- match( implstr[has_impl2], implstr, 0)
  cimpl2 <- cimplied[ second_impl,]

  pp_true[ cimplied[ first_impl,]] <- pp3[ -has_impl2] ? 0
  pp_true[ cimpl2] <- pp_true[ cimpl2] + pp3[ has_impl2] ? 0

  # Allow for XX <-> XO errors--- hopefully the only ones! (Watch out for scaffoldy version)
  # Assumes errors INDEPENDENT even tho there could be heritability (eg weak grabbing of mutated primer)
  # CC <-> CO is ignored in this version, since C gets scored as O

  add_pairprob_error()

return( pp6_err)
}


#' @importFrom atease @ @<-
#' @importFrom tools md5sum
#' @importFrom mvbutils %such.that% %not.in% %except%
"check_baits" <-
function(
    pbtidy= NULL,

    # Names in default lists "irrel" since should never actually be accessed,
    # ... but says what real names should be
    args_read_in= list(
        filename= NULL,
        filtered= NA
      ),
    args_zap_contams= list(
        contam_file= 'KEEP EM ALL FOR NOW'
      ),
    args_dump_low_count_loci= list(
        min_med_het_tot= stop(),
        het_pc_limit= stop(),
        mintol_het= stop(),
        nearly_max= stop()
      ),
    args_check_polyploids= list(
        polycutoff_ppn= stop(),
        maxpolysamps= stop(),
        polyrange_to_inspect= NULL
      ),
    args_refalt=list(
        OOthresh_tc= stop(),
        OK_med_min= stop(),
        selector= stop(),
        rescale_alleles= stop(),
        plot.=FALSE
      ),
    args_geno_deambig_ABC= list(
        OOthresh_tc= stop(),
        het_cut= stop(),
        tc_hist_pars= list(), # OK to have default, for plot; no comput effects;
        # ... having default should mean non-triggering re-run
        ppnA_hist_pars= list()
      ),
    args_est_ALF_ABCO= list(
        max_useful_prob= stop() # eg 0.95
      ),
    args_genofreq_check_4way= list(
        thresh_pchisq_loci= stop(),
        test= stop(),
        seq_paxis=NULL # graphical, default OK
      ),
    args_dump_badhetz_fish= list(
        multhresh_badhetz_fish= stop(),
        badhetz_hist_pars= list(),
        method= stop()
      ),
#    args_dump_genomic_paralogs=list(
#        which_BLAST=stop()
#      ),
    args_choose_geno6_thresholds= list(
        nquants_bump= stop(),
        max_dat_quantile= stop(),
        distro= stop(),
        minpO_rethresh= stop(),
        m_for_rethresh= stop(),
        max_refits= stop(),
        minbin= NULL, # only for plots
        plots= NULL,
        comment_file= 'locus_comments.txt'
      ),
    args_check6and4= list(
        thresh_pchisq_6and4= stop(),
        manual_6and4_file= stop() # '' for "everything fine"
      ),
    args_hsp_power= list(
      ),
    stop_after= 0, # do 'em all--- only thing with default
    mtrace.= FALSE
) {
## "Pipeline" for checking adequacy of baited loci
## Generate "self-documenting" dataset, with R class 'locar'

## self-configuration
## read in
## weed out duplicate fish--- for now
## normalize counts
## weed out contam fish: 3-alleles
## Zap low-count loci
## >2-allele loci zap
## pick Ref/Alt per cluster
## rough MALF/NALF
## threshold setting
### OO
### ambig AA/AB/BB
### (ambig fish)
### zap over-ambig loci
#
## genotyping (or in paige.mix?)
## check P4 fit (paralogs) (or in paige.mix?)
## paige.mix

#### Here we go:

 og <- options( vecless.print=FALSE) # suppress renamex(...) output!
 on.exit( options( og))

 define_genotypes() # AA etc, as variables: mode raw

## Self-configuration (prepare to record steps; always updating 'pbtidy')
setup_recordo_check_baits()

## Read in
  recordo( 'read_in', {
    if( !is.null( pbtidy)) {
      if( !is.null( filename)) { # Check it's the same. Flawed in that only one file is allowed...
        MD5 <- md5sum( filename)
        if( !all( MD5 == pbtidy@info$MD5)) {
          scatn( 'File "%s" has changed--- re-reading', filename)
          pbtidy <- NULL
        }
      }
    }

    if( is.null( pbtidy)) { # force read from file
      redo <- TRUE
      pbtidy <- read_cluster_dart3( filename, filtered=filtered, use_rownames=FALSE)
    }
  })

## Weed out duplicate fish--- for now
  # need them later for error rates
  recordo( 'drop_replicates', {
    dupfish <- duplicated( pbtidy@info$Our_sample)
    pbtidy <- pbtidy[ !dupfish,]
  })

## Normalize counts fishwise
  recordo( 'norm_by_fish', {
    ft <- pbtidy@info$Fishtot  # total number of reads according to DArT (always more than they report to us)
    pbtidy[] <- unclass( pbtidy) * mean( ft) / ft
    pbtidy@mean_fish_tot <- mean( ft)
  })

## Weed out contam fish: 3-alleles ideal method (?), but needs to loop with paralog-chop
  # for now just kill the ones in the file...
  recordo( 'zap_contams', {
    if( file.exists( contam_file)) {
      badhetz_fishnames <- scan( contam_file, what='', quiet=TRUE) # SOME have too few hetz!!
      # drop comments in file to say where it's from
      badhetz_fishnames <- badhetz_fishnames %such.that% (substring( ., 1, 1) != '#')
      pbtidy <- pbtidy[ pbtidy@info$Our_sample %not.in% badhetz_fishnames,]
    }
  })

## Counts too low (per locus)
  recordo( 'dump_low_count_loci', # also adds 'med_het_tot' field to locinfo
    pbtidy <- dump_low_count_loci( pbtidy,
        min_med_het_tot=min_med_het_tot,
        het_pc_limit=het_pc_limit,
        mintol_het=mintol_het,
        nearly_max=nearly_max,
        show_progress=TRUE)
  )

## 3+-allele cluster zap
  # How many fish with 3+ biggish scores? IE "polyzygotes".
  # Not rock-solid (EG paralog with null, but no SNP at either place)
  # ... but should pick up  *most* paralogs if number-of-fish is reasonably big
  # XS hetz should pick up the rest
  # Much more sophisticated version *could* try to find multiple genuine clones within each cluster...
  # ... based on (if alleles ABCD occur) rarely seeing a total of more than 2.5M within AB and CD
  # ... but often seeing eg 4M in AC and in BD etc
  # ... But that's a lot of work (DArT probably does this), so for now just kill 'em all
  recordo( 'check_polyploids', {
    pbtidy <- check_polyploids( pbtidy,
        cutoff_ppn=polycutoff_ppn,
        cutoff_npolysamps=max_npolysamps,
        range_to_inspect=polyrange_to_inspect,
        show_progress=TRUE)
  })

## Ref & Alt chosen by locus, based on total counts: no options here, but nice to have this version saved
  recordo( 'refalt', {
    pbtidy <- pick_ref_alt( pbtidy) # counts now stored in Ref/Alt/Others F*L*3 array masquerading as 'loc.ar'
    pbtidy <- ppn_ref_alt_check2( pbtidy,
      OOthresh_tc= OOthresh_tc,
      OK_med_min= OK_med_min,
      selector= selector,
      rescale_alleles= rescale_alleles,
      plot.=plot.)
  })

## Apply thresholds and score genotypes--- maybe remove ambig fish & loci
  recordo( 'geno_deambig_ABC', {
    pbtidy <- geno_deambig_ABC( pbtidy,
        OOthresh_tc= OOthresh_tc,
        het_cut= het_cut,
        ppnA_hist_pars= ppnA_hist_pars,
        tc_hist_pars= tc_hist_pars
      )
  })

## MALF/NALF estimates--- also adds attribs for obs & pred geno4 per locus
  recordo( 'est_ALF_ABCO', {
    pbtidy <- est_ALF_ABCO( pbtidy)
    # Junk pointlessly uninformative loci
    pbtidy <- pbtidy[, pbtidy@locinfo$pambig[,1] %in.range% c( 1-max_useful_prob, max_useful_prob)]
  })

## How do they look? Paralogs and other locus weirdness...
  recordo( 'genofreq_check_4way',
    pbtidy <- genofreq_check_4way( pbtidy, thresh_pchisq_loci= thresh_pchisq_loci, seq_paxis=seq_paxis, test=test)
  )

## Final bad fish check: XS hetz
  recordo( 'dump_badhetz_fish', {
    pbtidy <- dump_badhetz_fish( pbtidy, multhresh_badhetz_fish, badhetz_hist_pars, method=method)
    # Next was in chisq_geno_freq, but that seems "heavy".
    # I guess it's true that we don't need these atts now
    # ... but we do need them for newer xshetz check
    pbtidy@gpred <- pbtidy@gobs <- NULL
    pbtidy@subset_like_loci <- pbtidy@subset_like_loci %except% cq( gpred, gobs)
  })

## Prep for 6-way genotyping. Lots of graphs
  recordo( 'choose_geno6_thresholds', {
    pbtidy <- choose_geno6_thresholds( pbtidy,
      nquants_bump= nquants_bump,
      max_dat_quantile= max_dat_quantile,
      distro= distro,
      minpO_rethresh= minpO_rethresh,
      m_for_rethresh= m_for_rethresh,
      max_refits= max_refits,
      minbin= minbin,
      plots= plots,
      show_progress=TRUE,
      comment_file=comment_file)
  })

## 6-way genotyping: this discards orig counts
  recordo( 'geno6way', {
    pbtidy <- geno6way( pbtidy)
  })


## Check 6way GoF and decide 4-or-6
  recordo( 'check6and4', {
    pvals <- check6and4( pbtidy,
      thresh_pchisq_6and4= thresh_pchisq_6and4,  # can use 2 numbers to highlight marginal ones
      return_what='just_pvals')
    pbtidy@locinfo$pval6 <- pvals$pval6
    pbtidy@locinfo$pval4 <- pvals$pval4

    # Only now do the selection/elimination, manually
    pbtidy@locinfo$use6 <- TRUE
    if( nzchar( manual_6and4_file) && file.exists( manual_6and4_file)) {
      # Had trouble with next line--- factors, and spaces :/
      # dodgy <- read.table( manual_6and4_file, header=TRUE, row=NULL, sep=':')
      dodgy_text <- scan("jcnov-geno-comments.txt", sep='\n', what='', quiet=T)

      # Needs to start Locus: Comment
stopifnot( grepl( '^ *Locus *:', dodgy_text[1]))

      dodgy <- data.frame( Locus=sub( ' +', '', sub( ' *:.*', '', dodgy_text)),
          Comment=sub( ' +$', '', sub( '.*: *', '', dodgy_text)),
          stringsAsFactors=FALSE)[ -1,]

      # Any line where first non-space is 4, gets treated as 4way
      # Any other type of comment causes locus to be dropped
      as4f <- grepl( '^ *4', dodgy$Comment)
      if( any( as4f)) {
        as4 <- match( dodgy$Locus[ as4f], pbtidy@locinfo$Locus, 0)
        if( !all( as4)) {
  stop( sprintf( "Can't find locus [%s] to use-as-4...", paste( dodgy$Locus[ which( as4f)[ !as4]], collapse=', ')))
        }
        pbtidy@locinfo$use6[ as4] <- FALSE
      }
      if( any( !as4f)) {
        discardo <- match( dodgy$Locus[ !as4f], pbtidy@locinfo$Locus, 0)
        if( !all( discardo)) {
  stop( sprintf( "Can't find locus [%s] to discard...",
      paste( dodgy$Locus[ which( !as4f)[ !discardo]], collapse=', ')))
        }
        pbtidy <- pbtidy[ , -discardo]
      }
      pbtidy@manual_6and4 <- manual_6and4_file
    } # if nzchar file
    pbtidy
  })

  recordo( 'hsp_power', {
    pbtidy <- hsp_power( pbtidy, want_LOD_table=TRUE, k=0.5)
    scatn( 'Ediff HSP: %5.3f', sum( pbtidy@locinfo$Ediff))
    scatn( 'SD(UP): %5.3f', sqrt( sum( pbtidy@locinfo$V.UP)))
    # plot( something)
    pbtidy
  })


## (Moved here from higher up, but maybe shouldn't do at all...) Paralog check from genome assembly
#  recordo( 'dump_genomic_paralogs',
#    if( !is.null( which_BLAST)) {
#      pbtidy <- dump_genomic_paralogs( pbtidy, which_BLAST=which_BLAST)
#    }
#  )
#

## DID HAVE THIS BUT NOT USING IT NOW... Rough MALF/NALF estimates. Not required for solo clones (biallelic by defn) and
#  recordo( 'rough_MALF_NALF',
#    pbtidy <- rough_MALF_NALF( pbtidy, hetrule)
#  )


    # make Ref commoner
    # threshold setting
    ## OO
    ## ambig AA/AB/BB
    ## (ambig fish)
    ## zap over-ambig loci

    # genotyping (or in paige.mix?)
    # check P4 fit (paralogs) (or in paige.mix?)
    # paige.mix

return( pbtidy)
}


#' @importFrom mvbutils cq %without.name%
#' @importFrom handy2 sqr
#' @importFrom atease @ @<-
"check_FPosity" <- function( snpg, nsim=0){
## snpg should have been thru 'prepare_PLOD_SPA' so it has @PPS
stopifnot( 'Kenv' %in% names( attributes( snpg)))

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  define_genotypes()
  for( iwhat in cq( LOD, PUP, PUPLOD, PUPLOD2)) {
    assign( 'O' %&% iwhat, snpg@Kenv[[ iwhat]])
  }
  mg <- OLOD@mg

  use6 <- snpg@locinfo$use6
  use4 <- !use6
  temp_snpg <- snpg
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  temp_snpg[ , use4] <- recode4to6temp( snpg[, use4]) # (AA,AO) -> AA; (BB,BO) -> BB

  # For 4way loci, temporarily treat XO as XX...
  # ... have already adjusted the LOD entries so that new_LOD6( XX/..) <- LOD4( XXO/..)
  # ... use the LOD that's in Kenv, where SPA is calculated

  # Genofreqs could've/should've be done at the start in hsp_power, but here will do

  # LOD and PUP are stored in compacted 2D form to save space... need to fix that
  # Can't quite do this with vecless!
  OLOD[ is.na( OLOD)] <- 0 # set to NA for 4way loci
  n_loci <- nrow( OPUP)
  NPUP <- PUP <- LOD <- array( 0, c( n_loci, 6, 6))
  for( ig in 1:6) {
    gjseq <- mg[ , ig]
    XXi[ l, gj] := OPUP[ l, gj=gjseq] # Shouldn't work with new vecless syntax... but does !?
    PUP[ l, {ig}, gj] := XXi[ l, gj]
    NPUP[,ig,] <- OPUP[ , mg[,ig]]

    XXi[ l, gj] := OLOD[ l, gj=gjseq]
    LOD[ l, {ig}, gj] := XXi[ l, gj]
  }

  Pg[ l, gi] := sqrt( PUP[ l, gi, gi])
  # PHSP[l,gj,gi] := exp( LOD[ l, gj, gi]) * Pg[ l, gj] * Pg[ l, gi]
  # Pg2_g1_H[ l, gj, gi] := PHSP[ l, gj, gi] / Pg[ l, gi]

  e_CLOD[ l, gi] := LOD[ l, gi, gj] %[gj]% Pg[ l, gj]  # since gj indept gi
  e2_CLOD[ l, gi] := sqr( LOD[ l, gi, gj]) %[gj]% Pg[ l, gj]
  # e_CLOD_HSP[ l, gi] := LOD[ l, gi, gj] %[gj]% Pg2_g1_H[ l, gj, gi]
  # returnList( e_CLOD, e2_CLOD, e_CLOD_HSP)

  # Int version needed for vecless lookups
  geno <- as.integer( c( snpg))
  dim( geno) <- dim( snpg)

  # Ugly way to do lookups
  # my_e_CLOD[ i] := SUM_ %[l]% e_CLOD[ l, geno[ i, l]]
  # my_e2_CLOD[ i]:= SUM_ %[l]% e2_CLOD[ l, geno[ i, l]]
  g1seq <- seq_along( genotypes6)
  my_e_CLOD[ l, i] := e_CLOD[ l, g1] %[g1]% (geno[i,l] == g1seq[ g1])
  my_e2_CLOD[ l, i]:= e2_CLOD[l,g1] %[g1]% (geno[i,l] == g1seq[ g1])
  my_v_CLOD[ l, i]:= my_e2_CLOD[ l, i] - sqr( my_e_CLOD[ l, i])
  my_rat_CLOD[ i]:= (SUM_ %[l]% my_e_CLOD[ l, i] ) / sqrt( SUM_ %[l]% my_v_CLOD[ l, i])

  if( nsim) {
    gsim <- matrix( 0L, nsim, n_loci)
    for( il in seq_len( n_loci)) {
      gsim[,il] <- rsample( nsim, 1:6, prob=Pg[il,], replace=TRUE)
    }

    sim_e_CLOD[ l, i] := e_CLOD[ l, g1] %[g1]% (gsim[i,l] == g1seq[ g1])
    sim_e2_CLOD[ l, i]:= e2_CLOD[l,g1] %[g1]% (gsim[i,l] == g1seq[ g1])
    sim_v_CLOD[ l, i]:= sim_e2_CLOD[ l, i] - sqr( sim_e_CLOD[ l, i])
    sim_rat_CLOD[ i]:= (SUM_ %[l]% sim_e_CLOD[ l, i] ) / sqrt( SUM_ %[l]% sim_v_CLOD[ l, i])
  }


  # Vectorized individual KGFs, for each sample (columns) and numerous t-values (row)
  K <- function( tt) {
    # need to lookupize this until vecless 2.0 is out...
    ETT[ it, i, l, g] := exp( tt[ it] * LOD[ l, geno[ i, l], g])
    log_S[ it, i, l] := log( Pg[ l, g] %[g]% ETT[ it, i, l, g])
    K[ it, i] := log_S[ it, i, +.]
  return( K)
  }

  dK <- function( tt) {
    ETT[ it, l, g12] := exp( tt[ it] * LODOK[ l, g12])
    S[ it, l] := PUP[ l, g12] %[g12]% ETT[ it, l, g12]
    SL[ it, l] := PUPLOD[ l, g12] %[g12]% ETT[ it, l, g12]
    rowSums( SL/S)
  }

  ddK <- function( tt) {
    ETT[ it, l, g12] := exp( tt[ it] * LODOK[ l, g12])
    S[ it, l] := PUP[ l, g12] %[g12]% ETT[ it, l, g12]
    SL[ it, l] := PUPLOD[ l, g12] %[g12]% ETT[ it, l, g12]
    SLL[ it, l] := PUPLOD2[ l, g12] %[g12]% ETT[ it, l, g12]
    rowSums( (SLL/S-sqr( SL/S)))
  }



stop()

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  if( is.null( bins)) {
    qq <- (2:nq-1)/nq
    bins <- inv_CDF( qq)
  }
  binprobs <- CDF( bins)

  mean_theory <- dK( 0)
  var_theory <- ddK( 0)

  # Trying special-cases here to minimize copying
  if( symmo) {
    if( !my.all.equal( subset1, 1 %upto% ncol( temp_snpg))) {
      temp_snpg <- temp_snpg[, subset1]
    }

    result <- HSP_cond_paircomps_lots(
      vec_LOD= LOD,
      geno1= temp_snpg,
      geno2= temp_snpg,
      e_CLOD= e_CLOD,
      e2_CLOD= e2_CLOD,
      e_CLOD_HSP= e_CLOD_HSP,
      e_typical_PLOD= mean_theory,
      v_typical_PLOD= var_theory,
      symmo= TRUE,
      eta= eta,
      min_keep_PLOD= keep_thresh,
      bins= bins)
  } else { # different subsets
stop( "Fix the non-symm code, bozo...")
    result <- HSP_cond_paircomps_lots( this+will+fail,
        pair_geno= temp_LOD@mg,
        LOD= t( temp_LOD),
        geno1= temp_snpg[ , subset1],
        geno2= temp_snpg[ , subset2],
        symmo= FALSE,
        eta= eta,
        min_keep_PLOD= keep_thresh,
        bins= bins
      )
  }

  result$bigs <- with( result, data.frame( PLOD=big_PLOD, i=big_i, j=big_j))
  result <- result %without.name% cq( big_PLOD, big_i, big_j)
  result$bins <- bins
  result$binprobs <- binprobs
  result$eta <- eta
  result$keep_thresh <- keep_thresh
  result$call <- sys.call()

return( result)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils cq %without.name% %except% as.cat
#' @importFrom utils capture.output
"check_polyploids" <- function(
    darta,
    cutoff_npolysamps,
    cutoff_ppn,
    range_to_inspect=NULL,
    show_progress=TRUE,
    info= darta@info,
    locinfo= darta@locinfo) {
### Be sure to drop contaminated fish...
  n_loci <- ncol( darta)
  n_fish <- nrow( darta)
  extract.named( locinfo)
  roundup <- 2 ^ (ceiling( log( n_loci, 2)) + 1)

  presence_threshold <- cutoff_ppn * 0.5 * rep( med_het_tot, n_alleles)
  tdc <- t( unclass( darta)) # for speed
  tpres <- tdc > presence_threshold

  # Ensure unique combo of fish and locus
  pres_at_lf <- tpres * outer( rep( 1 %upto% n_loci, n_alleles), (1 %upto% n_fish) * roundup, '+')
  tabbo <- table( pres_at_lf %except% 0) # each elt counts num alleles present at fish-locus

  clft <- as.numeric( names( tabbo))
  lf_tabbo <- cbind( fish= clft %/% roundup, locus= clft %% roundup)

  # Could/should also check ratio of 3rd-to-top count; only bona fide poly if "reasonable", say > 0.5
  # see geno_deambig for exact criterion used
  # Currently, just using higher cutoff for presence to handle both cases

  n_al_pres <- matrix( 0, n_fish, n_loci)
  n_al_pres[ lf_tabbo] <- c( tabbo) # c() avoids class <- table

  is_polyzyg <- n_al_pres > 2L

  n_polyzygs_by_locus <- colSums( is_polyzyg)
  n_polyzygs_by_fish <- rowSums( is_polyzyg)

  # For inspecting *loci*, you don't want contam fish which will be super-likely
  # ... to generate polyzygs for ANY locus

  # Here's useful (but I *didn't* remove dodgy fish first...):
#  D(13)> table( n_polyzygs_by_locus)
#  n_polyzygs_by_locus
#     0    1    2    3    4    5    6    7   17   23   24   28  ...  870  947  970 1004 1118 1229 1300 1419 1426
#  1618  119   17    3    2    2    2    1    2    1    1    1  ...    1    1    1    1    1    1    1    1    1
#

  # So it looks like there's a clear gap in the distro from 17 up
  # I would set eg 'range_to_inspect < c( 2, 25)' in the debugger
  # ... and execute the loop below

  if( length( range_to_inspect)) {
    # Improve "printing experience"
    adart <- darta
    adart@calls <- adart@args <- NULL
    adart@info <- adart@info %without.name% cq( Our_sample, File, MD5)
    these <- character()
    for( iloc in which( n_polyzygs_by_locus %in.range% range_to_inspect)) {
      this_one <- capture.output( file=NULL,
          print( adart[ is_polyzyg[ ,iloc], iloc], dot_for_0=TRUE))
      these <- c(
          these,
          character( 2),
          formatC( width=max( nchar( this_one)),
              sprintf( 'MEDIAN SINGLE: %4.0f', med_het_tot[ iloc])),
          this_one)
    }
    # In the debugger, do handy2::CC2( as.cat( these)) to bring the results up in your text
    # editor. If you don't have handy2 and also 'fixr' set up, you will have to do it another way
    # which will be worse
    # and take longer
    # so I recommend you do get those things sorted out!
return( as.cat( these))
  }

  darta@locinfo$npoly <- n_polyzygs_by_locus
  darta@info$npoly <- n_polyzygs_by_fish

  darta <- darta[ , n_polyzygs_by_locus <= cutoff_npolysamps]
return( darta)
}


#' @importFrom mvbutils returnList
#' @importFrom atease @ @<-
"check6and4" <- function( geno6,
    thresh_pchisq_6and4,
    return_what=c( 'just_pvals', 'all')) {
##########
  n_fish <- nrow( geno6)
  n_loci <- ncol( geno6)
  define_genotypes()
  diplos <- geno6@diplos
stopifnot( my.all.equal( sort( genotypes6), sort( diplos)))

  p6 <- with( geno6@locinfo,
      calc_g6probs( pbonzer[,'A'], pbonzer[,'B'], pbonzer[,'C'],
          snerr=snerr))[,genotypes6]
  # At some point, will also need p6_IBD... later...

  # Chi-sq check: requires gpred & gobs attr on geno6
  g6pred <- p6 * n_fish
  g6obs <- matrix( 0, n_loci, 6, dimnames=list( NULL, genotypes6))

  for( ig in genotypes6) {
    g6obs[,ig] <- colSums( geno6==match( ig, diplos))
  }

  pval6 <- chisq_genofreq_check( geno6, gobs=g6obs, gpred=g6pred, test='G',
      thresh_pchisq_loci=thresh_pchisq_6and4, trim=FALSE)@locinfo$pval

  g4obs <- gtab6to4( g6obs)
  g4pred <- gtab6to4( g6pred)
  pval4 <- chisq_genofreq_check( geno6, gobs=g4obs, gpred=g4pred, test='G',
      thresh_pchisq_loci=thresh_pchisq_6and4, trim=FALSE)@locinfo$pval

  return_what <- match.arg( return_what)
  if( return_what=='all') {
    geno6@locinfo$pval6 <- pval6
    geno6@locinfo$pval4 <- pval4
return( geno6)
  } else {
return( returnList( pval6, pval4))
  }
}


#' @importFrom atease @ @<-
#' @importFrom handy2 sqr
#' @importFrom stats pchisq
#' @importFrom graphics plot points mtext legend
"chisq_genofreq_check" <-
function( lociar,
    gpred= lociar@gpred,
    gobs= lociar@gobs,
    thresh_pchisq_loci,  # NULL to not worry; 1 value a threshold; 2 vals to inspect "iffy" ones
    test,  # 'Pearson' or 'G'
    trim, # TRUE to keep only above max thresh_pchisq_loci. Arguably better done post hoc...
    seq_paxis=0.025) {
##########
# Either 6- or 4-geno version should work
# 1 DoF in either case
# Assumed null distro of chisq(1) is pretty approximate
  n_loci <- nrow( gobs)

  # After ML, should never happen that gobs>0 & gpred==0... but we'll check
stopifnot( !any( gobs>0 & gpred==0))

  chistat <- if( test=='Pearson')
      rowSums( sqr( gobs - gpred) / gpred)
    else if( test=='G') # must handle 0log0 which is 0 but R doesn't know that (cf nlogp func in my Pascal armoury)
      2 * rowSums( gobs * log( ifelse( gobs>0, gobs/gpred, 1)))
    else
stop( 'test must be "Pearson" or "G"')

  DoF <- ncol( gpred)-3 # ... maybe ... !?
  pval <- pchisq( chistat, df=DoF, lower.tail=FALSE) ### df = ???
  lociar@locinfo$pval <- pval
  liffies <- length( thresh_pchisq_loci)
  keep_loci <- if( liffies) pval > max( thresh_pchisq_loci) else rep( TRUE, n_loci)
  iffy_loci <- !keep_loci & (pval > min( thresh_pchisq_loci))

  # Plot histogram of pvals - should be approximately uniformly distributed if the loci are behaving as we would like
  # [ should follow recordo paradigm as per geno_deambig ]
  par( mfrow=c(1,1)) # just one plot on first page
  hist(pval,
      main=sprintf( "%s GoF of %i-genotypes: pval from chisq( %i)", test, ncol( gpred), DoF),
      xlab="P-value: LOW == BAD",
      breaks=seq(0,1,seq_paxis),
      xlim=c(0,1))
  abline( v=thresh_pchisq_loci, col='red')

  # Locussy fits
  opar <- par(mfrow=c(2, ncol( gpred) %/% 2),
      mar=c( 3, 3, 0, 0)+0.1, oma=c( 2, 2, 3, 1))
  # omi=c(0,0,.6,0),mai=c(.8,.8,.2,.2)) Paige uses absolute margins
  on.exit( par( opar))

  gtypes <- colnames( gobs)
  for( g in 1:ncol( gpred)) {
    # all of 'em
    plot(gpred[,g], gobs[,g], pch=16, cex=0.4, col='lightblue', # pch='.' is too small
        xlab='', ylab='', main='',
        xlim=c( 0, max( c( gpred[,g], gobs[,g]))),
        ylim=c( 0, max( c( gpred[,g], gobs[,g]))))
    abline(0,1,col=8,lwd=2)
    points( gpred[ keep_loci,g], gobs[ keep_loci,g], pch=16, cex=0.4, col='green') # overplot to make visible against the line
    mtext( side=3, gtypes[ g], line=-1)
    # bad ones with X
    points( gpred[ !keep_loci, g], gobs[ !keep_loci, g], col='orange', pch=4, cex=0.6)
    if( any( iffy_loci)) {
      points( gpred[ iffy_loci, g], gobs[ iffy_loci, g], col='magenta', pch=16, cex=1)
    }
  }
  mtext( 'Green = "good"', side=3, cex=1.5, outer=TRUE, line=1.5)
  mtext( 'Expected', side=1, cex=1.5, outer=TRUE)
  mtext( 'Observed', side=2, cex=1.5, outer=TRUE)

  if( liffies) {
    legend( 'bottomright', pch=c( 4, rep( 16, liffies)), col=c( 'orange', if( liffies>1) 'magenta', 'green'), pt.cex=c( 0.6, if( liffies>1) 1, 0.6),
        legend=c( sprintf( '%4.1e < Pr ...', thresh_pchisq_loci[ 1]), if( liffies>1) sprintf( '... < %4.1e', thresh_pchisq_loci[2]), '... < 1') )
  }

  if( trim) {
    lociar <- lociar[ , keep_loci, ,drop=FALSE]
  }
return( lociar)
}


#' @importFrom mvbutils mlocal %except% %such.that% do.on
#' @importFrom tools md5sum
"deduce_dart_header_rows" <- function( nlocal=sys.parent()) mlocal({
  metafish <- matrix( unlist( strsplit( metafish, ',')),
      ncol=ncols, byrow=TRUE)
  # What are the header rows? Following usually gives no information...
  # metafish.type <- apply( metafish, 1, function( x)
  #     (x[ !is.fish.col] %except% '*') %such.that% nzchar(.)) # 5/5/2015: some are '' and some '*' :/

  metafish <- metafish[,is.fish.col]
  suppressWarnings( # fucking NA warning--- that's the fucking POINT
    metamean <- apply( metafish, 1, function(x) mean( as.numeric(x), na.rm=FALSE))
  )

  which.dartjob <- which( !is.na( metamean) & (metamean > 1000000000))
  which.fishtot <- which( !is.na( metamean)) %except% which.dartjob
  if( length( which.fishtot)>1) {
    # ... see gripe below
warning( 'Silly plate "names" apparently used (numerical)--- taking a guess...')
    which.fishtot <- which.fishtot %such.that% (metamean[.] > 10000) # ... guess, to weed out 1,2,etc
  }

  # Also need a check that fishtot really is total-reads, even if numeric. EG it is not in the giant ABT/NBT file... fuck knows what it is there.
  if( metamean[ which.fishtot] < 10000) {
warning( sprintf( 'WTF is line %i?? Supposed to be total reads per fish :/ Either job went wrong or datafile SNAFU', which.fishtot))
  }

  is.well <- grepl( '^[A-Z][0-9]{1,2}$', metafish)
  dim( is.well) <- dim( metafish)
  which.well <- which.max( rowSums( is.well)) # was: which( rowSums( is.well) > 10)  but mebbe not 10 fish
  is.fishname <- grepl( '^[^_]*_Bx[0-9]+_[A-Z][0-9]+', metafish)
  # is.sbtname <- grepl( 'S[bB][^_]*_Bx[0-9]+_[A-Z][0-9]+', metafish) # old SBT-specific version
  dim( is.fishname) <- dim( metafish)
  which.fishname <- which.max( rowSums( is.fishname)) # was: which( rowSums( is.fishname) > 10) see which.well

  # We SHOULD send Dart a plate specifier but eg School-shark first two plates are "1" and "2"
  # ... which is not very bloody informative :/
  # Code below worked for single-plate files, but a multi-plate analysis is different
  # which.ourplate <- which( apply( metafish, 1, function( x) all( x==x[1]))) %except% which.dartjob
  # So: by elimination I guess :/
  which.ourplate <- (1:nrow( metafish)) %except%  c( which.dartjob, which.fishtot, which.fishname, which.well)

  if( !length( which.fishname) && (length( which.ourplate)==2)) { # because SOMEONE did not supply Dart with names in the agreed format :/
    lu <- do.on( which.ourplate, length( unique( metafish[.,])))
    which.fishname <- which.ourplate[ which.max( lu)]
    which.ourplate <- which.ourplate[ which.min( lu)]
warning( sprintf( '"Fishname" not in agreed format. Guessing that row %i is "our plate name" and row %i is "fishname"',
        which.ourplate, which.fishname))
  }

stopifnot( length( which.fishname)==1 &&
    length( which.ourplate)==1 &&
    length( which.fishtot)==1 &&
    length( which.well)==1 &&
    length( which.dartjob)==1)

  metafish.type <- rep( '??', nrow( metafish))
  metafish.type[ which.well] <- 'Well'
  metafish.type[ which.ourplate] <- 'Our_plate'
  metafish.type[ which.fishtot] <- 'Fishtot'
  metafish.type[ which.dartjob] <- 'Dart_job'
  metafish.type[ which.fishname] <- 'Our_sample'

  metafish <- data.frame( t( metafish), stringsAsFactors=FALSE)
  metafish[ !is.na( metamean)] <- lapply( metafish[ !is.na( metamean)], as.numeric)
  names( metafish) <- metafish.type
  # rownames( metafish) <- wells[ is.fish.col] # not sure why needed... and wells may not be unique
  metafish$File <- filename
  metafish$MD5 <- md5sum( filename)

  nfish <- ncol( fishinfo)
})


#' @importFrom atease @ @<-
"define" <-
function( expr, expand_dim=FALSE) {
  expr <- substitute( expr)
  stopifnot( ((expr %is.a% '<-') || ((expr %is.a% 'call') && (expr[[1]]==as.name( '<<-')) ) ) && (expr[[2]] %is.a% 'name'))
  obj <- expr[[2]] # name/symbol
  name <- as.character( obj)

  eval.parent( expr)
  if( expand_dim && length( prefixdims)) { # known from environment
    prefixdims <- prefixdims # since these live in the parent, and substitute() won't find them
    dimorlen <- function( x) { if( is.null( d <- dim( x))) d <- length( x); d}
    eval.parent( substitute(
        obj <- structure( rep( obj, prod( prefixdims)), dim=c( prefixdims, dimorlen( obj)))
      ))
  }

  e <- environment( sys.function())
  eval.parent( substitute( {
    dimnames( obj) <- NULL # they just slow things down
    oldClass( obj) <- 'playback'
    obj@where <- e
    obj@whoami <- name
  }))
eval.parent( obj)
}


#' @importFrom mvbutils cq extract.named named mlocal
"define_genotypes" <- function( nlocal=sys.parent()) mlocal({
  ABCO <- named( cq( A, B, C, O))
  extract.named( ABCO) # A, B, C, and O

  genotypes <- cq( OO, AO, BO, AB, AA, BB, AAO, BBO, AC, BC, CO, CC, CCO)
  genotypes_ambig <- cq( OO, AB, AC, BC, AAO, BBO, CCO)
  genotypes4_ambig <- cq( OO, AB, AAO, BBO)
  genotypes6 <- cq( AA, AB, AO, BB, BO, OO)
  genotypes_C <- cq( AA, AB, AC, AO, BB, BC, BO, CC, CO, OO)

  for( ig in genotypes) {
    # assign( ig, structure( as.raw( match( ig, genotypes)), class='ABOSNP'))
    assign( ig, structure( ig, class='noquote')) # for nicer printing
  }
})


"defish_csv" <- function(csvfile, just_count_the_fish=FALSE, subset_fish=NULL,
                         loci=TRUE,
                         outfile=sub( '([.][^.]+)$', '-short\\1', csvfile)){

  first_line <- scan( csvfile, what='', sep='\n', nmax=1, quiet=TRUE)
  first_nonblank <- regexpr( '[^,*]', first_line)
  n_locinfo_fields <- length( strsplit( substring( first_line, 1, first_nonblank), split=',')[[1]])-1
  n_fish <- length( strsplit( first_line, split=',')[[1]]) - n_locinfo_fields

  if( just_count_the_fish) {
return( n_fish)
  }

  # Needs 'cut' utility
  # ... Actually would be much fater without it, using instead something like
  # scannee <- scan(..., what='', sep='\n') to read in all lines, then
  # line_len <- regexpr( '([^,]*,){2000}', scannee, perl=T)@match.length-1
  # cat( substring( humungo7, 1, line_len), sep='\n', file=outfile)
  # perl regex needed to handle large repeat counts eg {2000}
  # but that's untested for really large counts so I'm sticking with this for now...

stopifnot( nzchar( Sys.which( 'cut')[1]))

  subset_fish <- seq.int( n_fish)[ subset_fish]

  if( any( is.na( subset_fish))) {
stop( sprintf( "Invalid subset: only max %i available", n_fish))
  }

  if( any( duplicated( subset_fish))) {
    warn( "Duplicates in 'subset_fish' are being discarded") # but carry on
  }

  subset_fish <- sort( unique( subset_fish))
  subset_fish <- c(
      if( loci) seq.int( n_locinfo_fields), # else don't...
      subset_fish + n_locinfo_fields)

  if( !length( subset_fish)) {
stop( "Pick some fish (or at least set loci=TRUE), dammit")
  }

  # Concoct the "--fields=..." arg to 'cut', called 'subbo': contiguous ranges eg 2,3-7,9-11,12
  if( length( subset_fish)==1) {
    subbo <- as.character( subset_fish)
  } else { # Organize subset into contiguous ranges, ready for "cut"
    subgap <- which( diff( subset_fish) > 1)
    if( !length( subgap)) { # all together now
      subbo <- sprintf( '%i-%i', subset_fish[1], tail( subset_fish, 1))
    } else {
      subgroups <- matrix( 0, length( subgap)+1, 2)
      subgroups[ 1, ] <- subset_fish[ c( 1, subgap[ 1])]
      for( igap in 2 %upto% (length( subgap))) {
        subgroups[ igap, ] <- subset_fish[ c( subgap[ igap-1]+1, subgap[ igap]) ]
      }
      subgroups[ length( subgap)+1,] <- subset_fish[ c( subgap[ length( subgap)]+1, length( subset_fish))]

      subranges <- as.character( subgroups[,1])
      singletons <- subgroups[,2] == subgroups[,1]
      subranges[ !singletons] <- paste( subgroups[ !singletons,1], subgroups[ !singletons,2], sep='-')
      subbo <- paste( subranges, collapse=',')
    }
    if( max( subset_fish) == n_locinfo_fields + n_fish) { # drop final thing to avoid hunt for final comma
      subbo <- sub( '[0-9]+$', '', subbo)
    }
  }

  if( grepl( ' ', csvfile)) { # stupid stupid
    csvfile <- sprintf( '"%s"', csvfile)
  }

  OK <- system2( 'cut',
      args=c( '--delimiter=,', sprintf( '--fields=%s', subbo), csvfile),
      stdout=outfile,
      invisible=TRUE)

  if( OK!=0) {
stop( sprintf( 'Summat wrong... look at %s', outfile))
  }
return( outfile)
}


#' @importFrom atease @ @<-
"dim.NGS_count_ar" <-
function( x) c( nrow( unclass( x)), nrow( x@locinfo))


#' @importFrom atease @ @<-
#' @importFrom handy2 sqr
#' @importFrom stats approx
#' @importFrom mvbutils do.on
"dump_badhetz_fish" <-
function( lociar,  # either 'loc.ar' with 'geno_amb' field, or 'snpgeno'
    multhresh_badhetz_fish,
    badhetz_hist_pars,
    method # 'Nhetz' or, prolly better, 'hetzminoo'
  ){
##############
# multhresh defined rel to number of fish
# EG if there are 1000 fish and multhresh = 0.5, then we keep fish with no worse than 0.5/1000 prob of observed nhet
  define_genotypes() # AB etc
  n_fish <- nrow( lociar)
  gamb <- if( lociar %is.a% 'snpgeno') { # ready to roll on completed genotypes
      lociar # should some atts be stripped for speed
    } else {
      lociar@geno_amb # AB AC AAO OO etc
    }

  gamb@locinfo <- gamb@info <- NULL # harmless
  pamb <- lociar@locinfo$pbonzer # preferred
  if( is.null( pamb)) {
    pamb <- lociar@locinfo$pambig
  }

  # Code below expects
  if( my.all.equal( gamb@diplos, genotypes6)) {
    # No C in this encoding of genotypes; fudge the probs accordingly
    pamb[,'O'] <- pamb[,'O'] + pamb[,'C']
    pamb[,'C'] <- 0
  } else if( !my.all.equal( gamb@diplos, genotypes_ambig)) {
stop( "'dump_badhetz_fish' requires either 'genotypes6' or 'genotypes_ambig' encoding")
  }


  cutto <- multhresh_badhetz_fish / n_fish
  pcuts <- c( cutto, 1-cutto)

  if( method=='Nhetz') {
stop( "Haven't finished SPA version of Nhetz method yet")
    phet_loci <- 2*pamb[ ,'A'] * pamb[,'B']
    hetz_stat <- nhet <- rowSums( gamb==AB)

    # bsp <- bernoulli.sum.probs( phet_loci)
    # Should reeeeeeally do the interp backwards & forwards to avoid roundoff error---
    # ... these are fairly small probs. But...
    # patleast <- rev( cumsum( rev( bsp)))

    # bernoulli.sum.probs() is exact, but uses obsolete DLL
    # Use SPA instead--- code pinched from find_POPs
    K <- function( tt) {
        pel <- outer( phet_loci / (1-phet_loci), exp( tt))
        colSums( log1p( -phet_loci) + log1p( pel))
      }

    dK <- function( tt) {
        pel <- outer( phet_loci / (1-phet_loci), exp( tt))
        length( phet_loci) - colSums( 1/(1+pel))
      }

    ddK <- function( tt) {
        pel <- outer( phet_loci / (1-phet_loci), exp( tt))
        ipel <- 1/(1+pel)
        colSums( ipel) - colSums( sqr( ipel))
      }

    # check: dK(0) == sum( phet_loci)
    # ddK(0) == sum( p) - sum( phet_loci)^2
    CDF <- renorm_SPA_cumul( K, dK, ddK)$CDF
    patleast <- 1-CDF( ....) # or use renorm_SPA as below

    ncuts <- approx( patleast, seq_along( bsp)-1, xout=pcuts)$y
  } else if( method=='hetzminoo') {
    delta <- (gamb==AB) + (gamb==AC) + (gamb==BC) - (gamb==OO)
    pXY <- 2*(pamb[,'A'] * pamb[,'B'] + pamb[,'A']*pamb[,'C'] + pamb[,'B']*pamb[,'C'])
    pOO <- sqr( pamb[,'O'])
    comp_pq <- 1-pXY-pOO

    hetz_stat <- hetzminoo <- rowSums( delta)

    # The none-too-clear use of eval etc next, is so that the intermediate variables
    # ... are shared by the K-related functions, without needing to use '<<-'
    # I don't understand why it's necessary to add 'list( s=s)'; sometimes things
    # ... seems to work without that, but sometimes they do not. Ideally
    # ... (i) 's' would not need to be 'force'd, and
    # ... (ii) there would be an automatic way to add in all the args
    # ... presumably all this crap is why 'do.in.envir()' really does have a role!

    Kdiff <- function( s) eval( envir=environment( sys.function()), substitute({
        last_s <- s
        es <- exp( s)
        piepl <- pXY * es
        qienl <- pOO * (1/es)
        pPq <- piepl + qienl
        pNq <- piepl - qienl
        pPqC <- pPq + comp_pq
      return( sum( log( pPqC)))
      }, list( s=s)))
    dKdiff <- function( s) eval( envir=environment( sys.function()), substitute({
        if( !identical( s, last_s)) {
          Kdiff( s)
        }
        last_s <- s
        last_dK <- pNq / pPqC
      return( sum( last_dK))
      }, list(s=s)))
    ddKdiff <- function( s) eval( envir=environment( sys.function()), substitute({
        Kdiff( s)
      return( sum( ( pPqC * pPq - sqr( pNq)  )  / sqr( pPqC)))
      }, list(s=s)))
    e <- new.env( environment())
    e$last_s <- Inf
    # evalq( env=e, { last_s <- last_dK <- pPq <- pNq <- pPqC <- 0 })
    environment( Kdiff) <- environment( dKdiff) <- environment( ddKdiff) <- e
    ncuts <- do.on( pcuts, inv_CDF_SPA2( ., Kdiff, dKdiff, ddKdiff) )
    dens_SPA <- renorm_SPA( Kdiff, dKdiff, ddKdiff, 'func')
  } else {
stop( "Method==WTF?? %s")
  }

  # optional graphics and/or user-specified outputs
  switch( mode( badhetz_hist_pars),
    list = {
        badhetz_hist_pars <- add_list_defaults( badhetz_hist_pars,
            main=sprintf( '%s by FISH: multhresh=%5.2f', method, multhresh_badhetz_fish),
            xlim= range( c( hetz_stat, ncuts)), # so cutoff lines show
            xlab='', nclass=50)
        lv <- do.call( 'hist', c( list( x=hetz_stat), badhetz_hist_pars))
        if( method=='hetzminoo') with( lv, # then plot predicted density. Could do but harder for other
          lines( mids, diff( breaks) * dens_SPA( mids) * sum( counts), col='green')
        )
        abline( v=ncuts, col='red')
      },
    expression = eval( badhetz_hist_pars),
    NULL = NULL
  )

  keep_fish <- hetz_stat %in.range% ncuts
  lociar <- lociar[ keep_fish,,drop=FALSE]

return( lociar)
}


#' @importFrom atease @ @<-
#' @importFrom abind abind
"dump_low_count_loci" <- function(countar,
    min_med_het_tot, # eg 100
    het_pc_limit, # eg 1.3
    mintol_het, # eg 0.1
    nearly_max, # 1 to use max( top2) * mintol_het as crit for hetz, <1 to use that quantile of top2 instead
    MByte_comfy=50, just_top2=FALSE, show_progress=TRUE) {
# Want a "robust" criterion for "this locus (AKA cluster) is crap"
# Mean counts are stuffed by Nulls
# Max counts are vulnerable to occasional Copy-Number Variation
# Use median total count for "heterozygotes"

# Memory is an issue here, so there's a few rm() calls en route

  n_samps <- nrow( countar)
  n_loci <- ncol( countar)

  top_count <- next_count <- matrix( 0L, n_samps, n_loci)

  median_het_total <- numeric( n_loci)
  keep_locus <- rep( TRUE, n_loci)
  names( median_het_total) <- names( keep_locus) <- loci_names <- countar@locinfo$Locus

  ## Find HIGHEST count by fish and locus, across *all* alleles at that locus...
  ## ... and 2nd highest count
  ## ... and pick "likely hetz" based on ratio of highest/2nd being "close to 1"
  ## ... and filter out too-tiny "hetz" that are really double-nulls with splatter
  ## ... and take the median

  # Some loci have LOADS of alleles. Expanding to a full-size array straightaway will eat too much memory
  # Rather, loop over groups of loci sorted by n_alleles, trying to avoid massive arrays
  n_alleles <- countar@locinfo$n_alleles
  tablo <- table( n_alleles)
  nal_group <- as.integer( names( tablo))
  MBytes <- 4 * n_samps * cumsum( tablo * nal_group)
  next_first_undone <- 1
  while( next_first_undone <= length( MBytes)) {
    # Do all locus-selection first, to avoid forgetting it after the "main" calcs!
    # Exceed MBytes_comfy, but by as little as possible...
    # ... can't risk not being able to do any at all!

    first_undone <- next_first_undone
    upto <- first_undone
    for( upto in first_undone %upto% length( MBytes)) {
      group_MBytes <- 4e-6 * n_samps * nal_group[ upto] * sum( tablo[ first_undone %upto% upto])
      if( group_MBytes > MByte_comfy) {
    break
      }
    }
    next_first_undone <- upto + 1
    if( show_progress) {
      cat( sprintf( '\r  Tackling %i loci with %i:%i alleles',
          sum( tablo[ first_undone %upto% upto]),
          first_undone,
          upto))
    }

    # Make regular 3D array--- NB 3rd dimension being (max) num alleles
    these_loci <- which( n_alleles %in.range% nal_group[ c( first_undone, upto)])
    counts <- countar[ , these_loci, 1 %upto% nal_group[ upto]]
    storage.mode( counts) <- 'integer'

    # Test another method...
    # Theoretically inefficient, since it sorts EVERYTHING and we already know most of the ordering
    # and we only need 2 alleles
    # but it avoids dreaded "apply"
    ordo <- order( slice.index( counts, 1), slice.index( counts, 2), -counts)
    ocounts <- aperm( counts, 3:1)
    ocounts[] <- counts[ ordo]
    ocounts <- aperm( ocounts, 3:1)
    top_count[,these_loci] <- ocounts[,,1]
    next_count[,these_loci] <- ocounts[,,2]
  } # while (n_allele - groups)
  if( show_progress) {
    cat( '\nDone\n')
  }

  rm( ordo, counts, ocounts)

  if( just_top2) { # used by geno_deambig()
    return( abind( top_count, next_count, along=3))
  }

  # Are they similar enough to be a het, say <=30% diff? (and both > 0)
  top2_counts <- top_count + next_count
  couldbe_het <- (next_count > 0) & (top_count <= het_pc_limit * next_count)

  # Trick for speed, using max.col to avoid apply. there's no max.row, so transpose
  if( nearly_max==1) {
    t_het_counts <- t( top2_counts * couldbe_het)
    maxsamp_by_locus <- max.col( t_het_counts)
    maxcount_het <- t_het_counts[ cbind( 1:n_loci, maxsamp_by_locus)]
    rm( t_het_counts)
  } else { # slightly slower thx2 apply
    t2 <- top2_counts
    t2[ couldbe_het] <- NA
    maxcount_het <- apply( t2, 2, quantile, probs=nearly_max, na.rm=TRUE)
    rm( t2)
  }

  # Eliminate tiny "hetz"--- really OO--- assuming CNVariation is not more than say 10fold...
  # ... could do some kind of graphical check here, dunno wot tho
  couldbe_het <- couldbe_het & (top2_counts > rep( maxcount_het, each=n_samps) * mintol_het)

  # No point if no hetz! true this criterion is *slightly* strict on hetz, but even so
  keep_locus[] <- colSums( couldbe_het) > 0
  nohetz <- loci_names[ !keep_locus]

  top2_counts[ !couldbe_het] <- NA
  median_het_total[keep_locus] <- apply( top2_counts[,keep_locus], 2, median, na.rm=TRUE)
  keep_locus <- keep_locus & (median_het_total >  min_med_het_tot)
  # ... just dumped the low ones
  # EG for checks: hist( median_het_total, nc=200, xlim=c( 0, 500))
  # abline( v=min_med_het_tot, col='red')

  rm( top2_counts, couldbe_het, top_count, next_count)

  countar@locinfo$med_het_tot <- median_het_total

return( countar[,keep_locus])
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils cq
#' @importFrom handy2 NEG
"est_ALF_ABCO" <- function(lociar){
########## Taken largely from "pipeline_for_SBT_baits.r"
########## MVB: I'd like to clean this up
########## Careful "parallel Newton-Raphson" could allow vectorization and whoosh-factor, but NFN I guess

  define_genotypes() # AAO etc
  geno_amb <- lociar@geno_amb
  if( is.null( geno_amb)) {
stop( "No 'geno_amb' attribute :(")
  }
  stopifnot( is.character( geno_amb) || my.all.equal( geno_amb@diplos, genotypes_ambig))

  expected <- NULL
  lglk <- function( params, nobs, return_expected=FALSE) {
      has_C <- length( params)==3

      # Reparamed for with-C case to logit scale, to avoid probs when pC~=0
      pA <<- inv.logit( params[1])
      pB <<- (1-pA) * inv.logit( params[2])
      pC <<- if( has_C) (1-pA-pB) * inv.logit( params[ 3]) else 0
      pO <<- max( 0, 1 - pA - pB - pC) # rounding error guard

      phat <- make_pgeno( pA, pB, pC, which_genotypes=genotypes_ambig)
      expected <<- n_fish * phat
      lglk <- nobs %*% log(phat + (nobs==0))        # avoid 0log0 gotcha
      pen <- penscale * sum( log( cosh( params-start_par)))
    return( lglk - pen)
    }
  pA <- pB <- pC <- pO <- (-1) # overwritten when lglk runs

  n_fish <- nrow( geno_amb)
  n_loci <- ncol( geno_amb)
  gobs <- gpred <- matrix( 0, n_loci, length( genotypes_ambig), dimnames=list( NULL, genotypes_ambig))
  for( g in genotypes_ambig) {
    gobs[,g] <- colSums( geno_amb==g)
  }

  # MVB: Old code looks pretty slow.
  # Calc g.freq for all loci **BUT ONLY ACCEPTABLE FISH** to preserve matrix size
  # g.freq <- apply(gABO.obs[!iamb.f, ], 2, function(x) table(factor(x,levels=c("AA","AB","BB","OO"))))

  pambig_est <- matrix(NA, n_loci, 4, dimnames=list( NULL, cq( A, B, C, O)))   # nloci rows, 2 cols (pA, pB) (p0 = 1-rowSums(p.est))
  conv <- rep( NA, n_loci) # convergence diagnostic

  tiny <- 2^-12 # avoid rounding error

  scatn( 'Starting ambig-geno_amb MALF/NALF estimation on %i loci:\n', n_loci)
  evalq( # for debug speed
  for( ll in 1:n_loci)  {
    if( ll %% 50 == 0) { cat( '\r', ll); flush.console() }
    has_C <- sum( gobs[ll,cq( AC, BC, CCO)])>0

    # Rough ests based on presence: mild overflow guard
    pA <- 1 - sqrt( 1- sum( gobs[ll, cq( AAO, AB, AC)]) / (1 + n_fish))
    pB <- 1 - sqrt( 1- sum( gobs[ll, cq( BBO, AB, BC)]) / (1 + n_fish))
    pC <- if( !has_C) 0 else 1 - sqrt( 1- sum( gobs[ll, cq( CCO, AC, BC)]) / (1 + n_fish))

    # Really, loci with rubbish pB (or pA) should have been chucked by now... but just in case...
    if( pA==0) pA <- 1/(2*n_fish)
    if( pB==0) pB <- 1/(2*n_fish)


    # Hard overflow guard...
    duhhh <- pA + pB + pC
    if( duhhh > 0.99) {
      duhhh <- duhhh + 0.011 # push it over 1
      pA <- pA / duhhh
      pB <- pB / duhhh
      pC <- pC / duhhh
    }


    start_par <- c( logit( pA), logit( pB / (1-pA)), if( has_C) logit( pC / (1-pA-pB)))

    # Set reasonable penalty scale
    penscale <- 0
    testo <- numeric( 3)
    for( i in 1:3) {
      testo[ i] <- lglk( start_par+c( -1, 0, 1)[i], nobs=gobs[ll,])
    }
    penscale <- max( abs( diff( testo))) / 1e4

    fit <- nlminb( start_par, NEG( lglk), nobs=gobs[ll,])
    conv[ll] <- fit$convergence
    besto <- fit$par

    # Try refit with reduced and recentred penalty
    start_par <- besto
    penscale <- penscale / 10
    fit <- nlminb( start_par, NEG( lglk), nobs=gobs[ll,])

    # Make sure 'expected' is up-to-date
    lglk( fit$par, nobs=gobs[ll,])
    gpred[ll,] <- expected
    pambig_est[ll,] <- c( pA, pB, pC, pO)
  }) # for, evalq

  scatn( 'Convergence result table (0 is ideal)')
  print( table( conv))

  # Gene frequency for subsequent analysis (pA, pB, p0)
  lociar@locinfo$pambig <- pambig_est
  lociar@gobs <- gobs
  lociar@gpred <- gpred # for chi-sq checks
  lociar@subset_like_loci <- c( lociar@subset_like_loci, 'gobs', 'gpred')

return( lociar)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils cq %without.name% returnList
#' @importFrom handy2 sqr
"find_HSPs_cond" <- function(snpg, subset1=1 %upto% nrow(snpg), subset2=subset1,
    one_in_X_eta,
    rough_n_pairs_to_keep,
    eta= NULL,
    keep_thresh= NULL,
    nq= 50,
    bins= NULL) {
## snpg should have been thru 'prepare_PLOD_SPA' so it has @PPS
stopifnot( 'Kenv' %in% names( attributes( snpg)))

  # Sanity...
stopifnot( is.numeric( subset1) && is.numeric( subset2))
stopifnot( all( !duplicated( subset1)) && all( !duplicated( subset2)))
stopifnot( my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)))

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  # Here I'm using L-R tail approx SPA for CDF
  # ... although Kenv$inv_CDF is likely more accurate for "moderate" tails but I don't quite trust it in the extremes
  # ... actually they are pretty similar
  # ... Possibly, Kenv$inv_CDF should check if arg exceeds the range it was fitted to, and if so call
  # ... inv_CDF_SPA2() instead
  # ... but the range used in fitting is very goddamn wide (say +/- 10 SD) !

  define_genotypes()

  for( iwhat in cq( K, dK, ddK, inv_CDF)) {
    assign( iwhat, snpg@Kenv[[ iwhat]])
  }
  set_thresholds( keeping='hi')

  # For 4way loci, temporarily treat XO as XX...
  # ... have already adjusted the LOD entries so that new_LOD6( XX/..) <- LOD4( XXO/..)
  # ... use the LOD that's in Kenv, where SPA is calculated

  make_CLOD <- function( LOD, PUP) {
    # Could've/should've be done at the start in hsp_power, but here will do
    Pg[ l, gi] := sqrt( PUP[ l, gi, gi])
    PHSP[l,gj,gi] := exp( LOD[ l, gj, gi]) * Pg[ l, gj] * Pg[ l, gi]
    Pg2_g1_H[ l, gj, gi] := PHSP[ l, gj, gi] / Pg[ l, gi]

    e_CLOD[ l, gi] := LOD[ l, gi, gj] %[gj]% Pg[ l, gj]  # since gj indept gi
    e2_CLOD[ l, gi] := sqr( LOD[ l, gi, gj]) %[gj]% Pg[ l, gj]
    e_CLOD_HSP[ l, gi] := LOD[ l, gi, gj] %[gj]% Pg2_g1_H[ l, gj, gi]
  returnList( e_CLOD, e2_CLOD, e_CLOD_HSP)
  }

  extract.named( snpg@locinfo[ cq( use6, LOD6, LOD4)])
  use4 <- !use6
  temp_snpg <- snpg
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  temp_snpg[ , use4] <- recode4to6temp( snpg[, use4]) # (AA,AO) -> AA; (BB,BO) -> BB
  temp_LOD <- LOD # from Kenv; already done in prepare_PLOD_SPA

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  if( is.null( bins)) {
    qq <- (2:nq-1)/nq
    bins <- inv_CDF( qq)
  }
  binprobs <- CDF( bins)

  mean_theory <- dK( 0)
  var_theory <- ddK( 0)

  # Trying special-cases here to minimize copying
  if( symmo) {
    if( !my.all.equal( subset1, 1 %upto% ncol( temp_snpg))) {
      temp_snpg <- temp_snpg[, subset1]
    }

    result <- HSP_cond_paircomps_lots(
      vec_LOD= LOD,
      geno1= temp_snpg,
      geno2= temp_snpg,
      e_CLOD= e_CLOD,
      e2_CLOD= e2_CLOD,
      e_CLOD_HSP= e_CLOD_HSP,
      e_typical_PLOD= mean_theory,
      v_typical_PLOD= var_theory,
      symmo= TRUE,
      eta= eta,
      min_keep_PLOD= keep_thresh,
      bins= bins)
  } else { # different subsets
stop( "Fix the non-symm code, bozo...")
    result <- HSP_cond_paircomps_lots( this+will+fail,
        pair_geno= temp_LOD@mg,
        LOD= t( temp_LOD),
        geno1= temp_snpg[ , subset1],
        geno2= temp_snpg[ , subset2],
        symmo= FALSE,
        eta= eta,
        min_keep_PLOD= keep_thresh,
        bins= bins
      )
  }

  result$bigs <- with( result, data.frame( PLOD=big_PLOD, i=big_i, j=big_j))
  result <- result %without.name% cq( big_PLOD, big_i, big_j)
  result$bins <- bins
  result$binprobs <- binprobs
  result$eta <- eta
  result$keep_thresh <- keep_thresh
  result$call <- sys.call()

return( result)
}


#' @importFrom mvbutils %such.that%
"find_new_OOthresh" <-
function( pO, mulo, siglo, dat, m_for_rethresh) {
  pcut <- 2*pO/(2*(1-pO))
  while( pcut < 0.84) { # pnorm( 1) ish
    cutto <- qnorm( pcut, mean=mulo, sd=siglo)
    # but that's theory; for small pO, cutto could be so small that too few from bump...
    if( sum( dat < cutto) - sum( dat==0) > 2 * m_for_rethresh) {
  break # successfully grabbed a few points
    } else {
      pcut <- pcut * 2
    }
  }

  if( pcut >= 0.84) {
    cutto <- mulo+siglo # really shouldn't need to go this high, but if pO is not crazy-small there must be a middle bump
  }

return( widio( c( 0, dat %such.that% ((.>0) & (.<cutto))), m=m_for_rethresh))
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils returnList
"geno_deambig" <- function(lociar, OOthresh_tc, het_cut, tc_hist_pars,
                           ppnA_hist_pars) {

  # Pick OO, hetz
  # Leave possible homoz as "possible" for now

  define_genotypes() # knows OO AB etc

  mht <- lociar@locinfo$med_het_tot
  if( is.null( mht)) {
stop( "Need med_het_tot column in locinfo")
  }
  count3ar <- unclass( lociar)

  # count3ar normed by locus
  rmht <- rep( mht, each=nrow( count3ar))
  cA <- count3ar[,,1,drop=TRUE] / rmht
  cB <- count3ar[,,2,drop=TRUE] / rmht
  cX <- count3ar[,,3,drop=TRUE] / rmht

  # Turns out there's little to gain from using cX at this 4-geno stage
  # ... it comes in at the 6-geno stage

  is_OO <- (cA+cB) < OOthresh_tc
  switch( mode( tc_hist_pars),
    list = {
        # Default settings for 'hist' that *can* be overridden, but...
        tc_hist_pars <- add_list_defaults( tc_hist_pars,
            xlim=c( 0, 1.5),
            main='Tot count3ar in FISHLOCI, normed so non-OO ~= 1',
            xlab='',
            nclass=1000)
        do.call( 'hist', c( list( x=c( cA+cB)), tc_hist_pars))
        abline( v=OOthresh_tc, col='red')
        abline( v=1, col='blue', lty=2)
      },
    expression = eval( tc_hist_pars),
    NULL = NULL
  )


  # Proportion of "A" count3ar--- only wanted for non-nulls, but easier to calc for all.
  # Clusterized organization means count 1 always >= count 2 so lower limit is 0.5
  ppnA <- cA / (cA + cB)

  # optional graphics and/or user-specified outputs
  switch( mode( ppnA_hist_pars),
    list = {
        # Default settings for 'hist' that *can* be overridden, but...
        ppnA_hist_pars <- add_list_defaults( ppnA_hist_pars,
            main='Ppn Ref count3ar in non-OO FISHLOCI, normed by locus',
            xlab='',
            nclass=200)
        do.call( 'hist', c( list( x=c( ppnA[ !is_OO])), ppnA_hist_pars))
        abline( v=c( het_cut, 1-het_cut), col='lightblue')
      },
    expression = eval( ppnA_hist_pars),
    NULL = NULL
  )

# interesting to split by mht, eg
#  hist( ppnA[,mht<80][ !is_OO[,mht<80]], nclass=100, ylim=c( 0, 1000))
#  hist( ppnA[,mht>200][ !is_OO[,mht>200]], nclass=100, ylim=c( 0, 10000))
#  ... didn't seem to make much diff actually for SBT baits#2 [with lowcount3ar at tc==50]. Which is good.

  is_AB <- !is_OO & (ppnA %in.range% c( het_cut, 1-het_cut))

  # couldbe_homoz <- !is_OO & !is_AB

  # 4-way genotypes: use raw for memory-size reasons. Lower-case!
  geno4 <- matrix( AB, nrow( count3ar), ncol( count3ar),
      dimnames=list( count3ar@info$Our_sample, count3ar@locinfo$AlleleID))

  geno4[ is_OO] <- OO
  geno4[ !is_OO & (ppnA < het_cut)] <- BBO
  geno4[ !is_OO & (ppnA > 1-het_cut)] <- AAO
  class( geno4) <- c( 'noquote', oldClass( geno4)) # so shows ABO not "ABO" etc

  lociar@geno4 <- geno4 # but don't print it...
  lociar@subset_like_both <- c( lociar@subset_like_both, 'geno4')

  lociar@print <- expression({
      attr( x, 'locinfo') <- attr( x, 'seqinfo') <- attr( x, 'geno4') <- attr( x, 'print') <- NULL
      class( x) <- class( x) %except% 'specialprint'
      print( x)
      cat( '\n\n<<"locinfo" and "geno4" attributes hidden for brevity>>\n')
    })


return( lociar)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils cq
"geno_deambig_ABC" <-
function( lociar,
  mht= lociar@locinfo$med_het_tot,
  OOthresh_tc,
  het_cut,
  tc_hist_pars,
  ppnA_hist_pars,
  return_what=c( 'lociar', 'just_geno', 'polyzyg_and_geno')) {
##############
# Pick OO, hetz
# Leave possible homoz as "possible" for now

  define_genotypes() # knows OO AB etc

  c3 <- unclass( lociar) # but can be passed in directly
stopifnot( length( dim( c3))==3 && dim( c3)[3]==3) # expect 3 allele counted (3rd might always be 0)

  # Zap irrel attribs for speed
  attributes( c3) <- list( dim=dim( c3), dimnames=list( NULL, NULL, cq( A, B, C)))
  n_loci <- ncol( c3)
  n_fish <- nrow( c3)

  # Norm by locus
  c3[f,l,all] := c3[f,l,all] / mht[ l]
  cA <- c3[,,A]
  cB <- c3[,,B]
  cC <- c3[,,C]
  # guard against length-1 dims cRazies
  dim( cA) <- dim( cB) <- dim( cC) <- dim( c3)[1:2]

  # Double nulls:
  if( length( OOthresh_tc)==1) {
    OOthresh_tc <- rep( OOthresh_tc, n_loci)
  } # else already locus-specific, eg when regenotyping in geno6()

  is_OO[f,l] :=
      ((cA+cB)[f,l] < OOthresh_tc[l]) &
      ((cA+cC)[f,l] < OOthresh_tc[l]) &
      ((cB+cC)[f,l] < OOthresh_tc[l])

  switch( mode( tc_hist_pars),
    list = {
        # Default settings for 'hist' that *can* be overridden, but...
        tc_hist_pars <- add_list_defaults( tc_hist_pars,
            xlim=c( 0, 1.5),
            main='Tot c3 in FISHLOCI, normed so non-OO ~= 1',
            xlab='',
            nclass=200)
        # For some reason, directly using 'do.call( "hist"...)' is INCREDIBLY slow
        # Dunno why
        # do.call( 'hist', c( list( x=pmax( c( cA+cB), c( cA+cC), c( cB+cC))), tc_hist_pars))
        yy <- pmax( c( cA+cB), c( cA+cC), c( cB+cC))
        # do.call( 'hist', c( list( x=yy), tc_hist_pars))
        FFS_gimme_a_histo <- c( list( quote( hist), quote( yy)), tc_hist_pars)
        eval( as.call( FFS_gimme_a_histo))
        abline( v=mean( OOthresh_tc), col='red')
        abline( v=1, col='blue', lty=2, lwd=3)
      },
    expression = eval( tc_hist_pars),
    NULL = NULL
  )

  # Proportion of "A" c3--- only wanted for non-nulls, but easier to calc for all.
  ppnA <- cA / (cA + pmax( cB, cC))
  ppnB <- cB / (cB + pmax( cA, cC))
  ppnC <- cC / (cC + pmax( cA, cB))

  # optional graphics and/or user-specified outputs
  switch( mode( ppnA_hist_pars),
    list = {
        # Default settings for 'hist' that *can* be overridden, but...
        ppnA_hist_pars <- add_list_defaults( ppnA_hist_pars,
            main='Ppn Ref c3 in non-OO FISHLOCI, normed by locus',
            xlab='',
            nclass=200)
        # See previous histo comments
        #  do.call( 'hist', c( list( x=c( ppnA[ !is_OO])), ppnA_hist_pars))
        yy <- ppnA[ !is_OO]
        # do.call( 'hist', c( list( x=yy), tc_hist_pars))
        FFS_gimme_a_histo <- c( list( quote( hist), quote( yy)), ppnA_hist_pars)
        eval( as.call( FFS_gimme_a_histo))
        abline( v=c( het_cut, 1-het_cut), col='lightblue')
      },
    expression = eval( ppnA_hist_pars),
    NULL = NULL
  )

  inhetrange <- function( x) x %in.range% c( het_cut, 1-het_cut)
  is_AB <- !is_OO & inhetrange(ppnA) & inhetrange( ppnB)
  is_AC <- !is_OO & inhetrange(ppnA) & inhetrange( ppnC)
  is_BC <- !is_OO & inhetrange(ppnB) & inhetrange( ppnC)

  # Only one locus in jcslack4 actually seemed like a real problem here
  is_polyzyg <- is_AB + is_AC + is_BC > 1 # hopefully very few... in these cases, pick the biggest pair
  A_lowest <- cA < pmin( cB, cC)
  B_lowest <- cB < pmin( cA, cC)
  C_lowest <- cC < pmin( cA, cB)

  is_AB[ is_polyzyg & (A_lowest | B_lowest)] <- FALSE
  is_AC[ is_polyzyg & (A_lowest | C_lowest)] <- FALSE
  is_BC[ is_polyzyg & (B_lowest | C_lowest)] <- FALSE

  # Use 1-byte 'snpgeno' format. requires info & locinfo atts; give it minimal ones
  geno_amb <- snpgeno( nrow( c3), ncol( c3), diplos=genotypes_ambig,
      info=lociar@info[ 'Our_sample'], locinfo=lociar@locinfo[ 'Locus'])
  # matrix( AB, nrow( c3), ncol( c3), dimnames=list( c3@info$Our_sample, c3@locinfo$Locus))
  # `[<-.snpgeno` allows character assignments. Not tested yet...

  geno_amb[] <- AB

  geno_amb[ is_OO] <- OO
  eps <- .Machine$double.eps   # since a feeeeeeew cases hit rounding error...
  geno_amb[ !is_OO & (ppnB+eps >= 1-het_cut)] <- BBO # for now; hetz will be fixed below
  geno_amb[ !is_OO & (ppnA+eps >= 1-het_cut)] <- AAO
  geno_amb[ !is_OO & (ppnC+eps >= 1-het_cut)] <- CCO

  geno_amb[ is_AB] <- AB
  geno_amb[ is_AC] <- AC
  geno_amb[ is_BC] <- BC

  # class( geno_amb) <- c( 'noquote', oldClass( geno_amb)) # so shows ABO not "ABO" etc

  return_what <- match.arg( return_what)

  if( return_what=='just_geno') { # used when regenotyping during 'geno6', at least originally
return( geno_amb)
  } else if( return_what=='polyzyg_and_geno') {
return( returnList( is_polyzyg, geno_amb))
  }

  lociar@geno_amb <- geno_amb # but don't print it...
  lociar@subset_like_both <- c( lociar@subset_like_both, 'geno_amb')
  lociar@het_cut <- het_cut

  lociar@print <- expression({
      attr( x, 'locinfo') <- attr( x, 'seqinfo') <- attr( x, 'geno_amb') <- attr( x, 'print') <- NULL
      class( x) <- class( x) %except% 'specialprint'
      print( x)
      cat( '\n\n<<"locinfo" and "geno_amb" attributes hidden for brevity>>\n')
    })


return( lociar)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils cq
#' @importFrom handy2 %<-%
#' @export
"geno6way" <- function( lociar, li=NULL, het_cut=lociar@het_cut) {
  n_loci <- ncol( lociar)
  n_fish <- nrow( lociar)

  define_genotypes()

  dimnames( lociar)[[3]] <- cq( A, B, C)

  if( is.null( li)) { # use own data
    li <- lociar@locinfo
  } else if( is.null( lociar@locinfo)) { # "raw" stuff passed; just check dims
stopifnot( ncol( lociar) == nrow( li))
  } else { # check preset definitions match this OK...
    mm <- match( lociar@locinfo$consensus, li$consensus, 0)
stopifnot( all( mm>0))
    li <- li[ mm,]
stopifnot( all( lociar@locinfo$FullAltSeq==li$FullAltSeq))
stopifnot( all( lociar@locinfo$FullRefSeq==li$FullRefSeq))
stopifnot( my.all.equal( lociar@locinfo$FullThirdSeq, li$FullThirdSeq)) # NA-friendly
  }

  c3 <- unclass( lociar)

  mht <- li$med_het_tot
  pambig <- li$pambig
  if( is.null( mht)) {
stop( "Need med_het_tot column in locinfo/li")
  }

  # counts normed by locus
  c3[f,l,all] := c3[f,l,all] / mht[ l]
  cA <- c3[,,A]
  cB <- c3[,,B]
  cC <- c3[,,C]
  # guard against length-1 dims cRazies
  dim( cA) <- dim( cB) <- dim( cC) <- dim( c3)[1:2]

  # 4-way genotype first:
  { is_polyzyg; geno4} %<-% geno_deambig_ABC( lociar, # does its own norming (duplicated effort, but...)
      mht=mht,
      OOthresh_tc=li$OOthresh_tc,
      het_cut=het_cut,
      tc_hist_pars= NULL,
      ppnA_hist_pars= NULL,
      return_what= 'polyzyg_and_geno')

  g6 <- split_geno_4to6( geno4, cA, cB, info=lociar@info, li) # snpgeno object
  g6@mean_fish_tot <- lociar@mean_fish_tot # this is needed later for rescaling new fish
  g6@het_cut <- het_cut
  g6@is_polyzyg <- is_polyzyg
  g6@subset_like_both <- c( g6@subset_like_both, 'is_polyzyg')
return( g6)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils cq
"genofreq_check_4way" <-
function( lociar,
    gpred= lociar@gpred,
    gobs= lociar@gobs,
    thresh_pchisq_loci,
    test,  # 'Pearson' or 'G'
    seq_paxis=0.025) {
##########
# Either 6- or 4-geno version should work
# 1 DoF in either case
# Assumed null distro of chisq(1) is pretty approximate
  gpred <- lociar@gpred
  gobs <- lociar@gobs
  n_loci <- nrow( gobs)
  define_genotypes()

  agg_g <- function( g) {
      g[,AAO] <- g[,AAO] + g[,AC]
      g[,BBO] <- g[,BBO] + g[,BC]
      g[,OO] <- g[,OO] + g[,CCO]
    return( g[, genotypes4_ambig])
  }

  # Do the geno4 check. don't trim yet, to keep dim of gobs & gpred consistent
  lociar <- chisq_genofreq_check( lociar, gpred=agg_g( gpred), gobs=agg_g( gobs), test=test,
      thresh_pchisq_loci=thresh_pchisq_loci, trim=FALSE, seq_paxis=seq_paxis)
  pv1 <- lociar@locinfo$pval

  # Check the C-stuff
  agg_g2 <- function( g) {
      gg <- matrix( g[1,1], nrow( g), 4, dimnames=list( NULL, cq( XYXXO, XC, CCO, OO)))
      gg[,'XYXXO'] <- g[,AAO] + g[,BBO] + g[,AB]
      gg[,'XC'] <- g[,AC] + g[,BC]
      gg[,'OO'] <- g[,OO]
      gg[,'CCO'] <- g[,CCO]
    return( gg)
    }

  # pre-4/2017 version didn't have BC in next--- looks like a mistake
  has_C <- gobs[,'AC'] + gobs[,'BC'] + gobs[,'CCO'] > 0.01 * nrow( lociar)
  pv2 <- 1+0*pv1
  pv2[ has_C] <- chisq_genofreq_check( lociar[ , has_C], gpred=agg_g2( gpred[ has_C,]), gobs=agg_g2( gobs[ has_C,]),
      test=test, thresh_pchisq_loci=thresh_pchisq_loci, trim=FALSE, seq_paxis=seq_paxis)@locinfo$pval
  lociar@locinfo$pval2 <- pv2

  keep_loci <- pmin( pv1, pv2) > thresh_pchisq_loci

  lociar <- lociar[ , keep_loci, ,drop=FALSE]
return( lociar)
}


#' @importFrom mvbutils %is.not.a% %where%
"get_chain" <- function( thing, seed) {
  if( thing %is.not.a% 'data.frame') {
    thing <- thing$bigs
  }
  extract.named( thing)

  oset <- integer()
  set <- seed
  while( length( set) != length( oset)) {
    newj <- j[ i %in% set]
    newi <- i[ j %in% set]
    oset <- set
    set <- unique( c( set, newi, newj))
  }

  thing %where% (i %in% set | j %in% set)
}


"gtab4" <-
function( x) {
#### Combine 3rd-allele genos into geno4 (with ambig)
  # FUCKING R keeps the FUCKING names on FUCKING subsets of FUCKING vectors
  # for FUCK's sake
  # hence the FUCKING calls to unname()
  result <- c(
      AAO=unname( x['AAO'] + x['AC']),
      AB=unname( x['AB']),
      BBO=unname( x['BBO'] + x['BC']),
      OO=unname( x['OO'] + x['CCO']))
return( result)
}


"gtab6to4" <-
function( gt6) {
######### Condense 6way-genotype counts to 4way (AAO instead of AA, AO)
  define_genotypes()
  gt4 <- matrix( 0, nrow( gt6), 4, dimnames=list( dimnames( gt6)[[1]], genotypes4_ambig))
  gt4[,AB] <- gt6[,AB]
  gt4[,OO] <- gt6[,OO]
  gt4[,AAO] <- gt6[,AA] + gt6[,AO]
  gt4[,BBO] <- gt6[,BB] + gt6[,BO]
return( gt4)
}


#' @importFrom mvbutils cq
#' @importFrom atease @ @<-
#' @importFrom handy2 sqr
#' @importFrom vecless compile_vecless
"hetzminoo_fancy" <-
function( snpg, target=c( 'rich', 'poor'), hist_pars=list(), multhresh=1) {
###################
  define_genotypes()
  extract.named( snpg@locinfo[ cq( use6, PUP4, pbonzer)])
  p0 <- pbonzer[,'O'] + pbonzer[,'C']
  pA <- pbonzer[,'A']
  pB <- pbonzer[,'B']

  v <- 2*pA*pB + sqr( p0) - sqr( 2*pA*pB-sqr( p0))
  target <- match.arg( target)
  edash <- if( target=='rich') {
      (2*pA*p0+sqr(pA)) * (1-sqr(1-pB)) + (2*pB*p0+sqr( pB))*(1-sqr(1-pB)) + sqr( p0) * (1-sqr( p0))
    } else {
      edash <- 2*(pA*pB + pA*p0 + pB*p0) # minus sign
    }

  ww <- edash / v
  ww <- c( ww / sum( ww) ) # else get 1-col matrix
stopifnot( all( ww>0))

  use_loci <- which( ww > 0) # all of them, for now
  temp_snpg <- snpg[ , use_loci]
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  temp_snpg <- recode4to6temp( temp_snpg) # (AA,AO) -> AA; (BB,BO) -> BB

  delta <- (temp_snpg==AB) - (temp_snpg==OO)
  whmo[ f]:= delta[ f, l] %[l]% ww[ l]

  # Null distro: P[S==1] = pAB; P[S==-1] = pOO; P[S==0] = 1-pAB - pOO
  # E[ exp( tt*S)] = (e(tt)+1) * pAB + 1 + (e(-tt) -1)* pOO

  pAB <- 2*pA*pB
  pOO <- sqr( p0)
  four_pab_poo <- 4*pAB*pOO
  compaboo <- 1 - pAB - pOO

  K <- function( tt) {
      etwab[ l, j] := pAB[l] * exp( tt[j] * ww[ l])
      etwoo[ l, j] := pOO[l] * exp( -tt[j] * ww[ l])
      KK[ j]:= SUM_ %[l]% log( compaboo[l] + etwab[ l, j] + etwoo[ l, j])
    return( c( KK)) # without the c(), you get a scalar xtensor, and trouble...
    }

  dK <- function( tt) {
      etwab[ l, j] := pAB[l] * exp( tt[j] * ww[ l])
      etwoo[ l, j] := pOO[l] * exp( -tt[j] * ww[ l])
      denom[ l, j] := compaboo[l] + etwab[l,j] + etwoo[l,j]
      dKK[ j] := ww[l] %[l]% ((etwab[l,j]-etwoo[l,j])/denom[l,j])
    return( c( dKK))
    }

  ddK <- function( tt) {
      etwab[ l, j] := pAB[l] * exp( tt[j] * ww[ l])
      etwoo[ l, j] := pOO[l] * exp( -tt[j] * ww[ l])
      denom[ l, j] := compaboo[l] + etwab[l,j] + etwoo[l,j]
      ddKK[ j] := sqr( ww[l]) %[l]% ( (compaboo[l] * (etwab[ l, j] + etwoo[l,j]) + four_pab_poo[ l]) / sqr( denom[l, j]) )
    return( c( ddKK))
    }

  K <- compile_vecless( K(0))
  dK <- compile_vecless( dK(0))
  ddK <- compile_vecless( ddK(0))

  dens_SPA <- renorm_SPA( K, dK, ddK, 'func')

  # optional graphics and/or user-specified outputs
  switch( mode( hist_pars),
    list = {
        hist_pars <- add_list_defaults( hist_pars,
            main=sprintf( '%s: multhresh=%5.2f', target, multhresh),
            xlim= range( whmo), # so cutoff lines show
            xlab='', nclass=50)
        lv <- do.call( 'hist', c( list( x=whmo), hist_pars))
        with( lv, lines( mids, diff( breaks) * dens_SPA( mids) * sum( counts), col='green'))
        # abline( v=ncuts, col='red')
      },
    expression = eval( badhetz_hist_pars),
    NULL = NULL
  )

return( c( whmo))
}


#' @importFrom atease @ @<-
"hsp_power" <-
function( lociar,
    want_LOD_table, # T/F
    k # 0.5 for HSPs
){
############
  define_genotypes()
  li <- lociar@locinfo
  li1 <- li[1,]

  temp0 <- with( li1, calc_g6probs_IBD0_scalar( pbonzer, snerr, record=TRUE))
  cg6p0 <- make_playback( calc_g6probs_IBD0_scalar, temp0)

  temp1 <- with( li1, calc_g6probs_IBD1_scalar( pbonzer, snerr, record=TRUE))
  cg6p1 <- make_playback( calc_g6probs_IBD1_scalar, temp1)

  g6p0 <- with( li, cg6p0( pbonzer, snerr))
  g6p1 <- with( li, cg6p1( pbonzer, snerr))

  s6 <- predict_hsp_util( g6p0, g6p1, want_LOD_table, k=k)

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

  s4 <- predict_hsp_util( g4p0, g4p1, want_LOD_table, k=k)

  if( want_LOD_table) {
    li$LOD6 <- s6@LOD # matrix
    li$PUP6 <- s6@PUP
    li$LOD4 <- s4@LOD
    li$PUP4 <- s4@PUP
    s6@LOD <- s6@LOD <- s4@LOD <- s4@PUP <- NULL
  }

  li <- cbind( li, s6)
  li[ !li$use6, names( s4)] <- s4[ !li$use6,]

  lociar@locinfo <- li
return( lociar)
}

"inv_CDF_SPA2" <-
function( p, K, dK, ddK, tol=formals( ridder)$tol) {
######## Invert L-R SPA approx to CDF on "s-scale"
######## Avoids "double iteration" of nonlinearity
######## Should work for vector p (the target) but only if your K etc do

  isqrt_2pi <- 1/sqrt(2*pi)
  x <- sqrt_ddK_s <- Leg_trans <- u <- w <- 0*p # vectorized

  CSPA <- function( s) {
      x <<- dK(s)
      Leg_trans <<- s*x - K(s)
      w <<- sign(s) * sqrt( 2*Leg_trans)
      sqrt_ddK_s <<- sqrt( ddK( s))
      u <<- s * sqrt_ddK_s
    return( pnorm( w) + dnorm( w) * (1/w - 1/u) - p_target)
    }

  p_target <- 0
  isK2 <- 1/sqrt(ddK(0))
  seps <- 0.001 * isK2 # fraction of 1 SD; (tol/2) tends to give numeric errors...
  p0 <- (CSPA( seps) + CSPA( -seps)) / 2 # CSPA(0)==NA--- avoid!
  is_lower <- p0 > p
  # K3 <- 3 * (sqr(u)-sqr(w)) / seps^3 # yeh not bad FWIW

  hi <- lo <- 0*p

  # Shifted start...
  q0 <- qnorm( inv.logit( logit( p) - logit( p0)))
  bingo <- q0==0.5 # this case won't converge
  mean_x <- dK( 0)

  # L-R does not work at x==E[X]..!
  if( any( bingo)) { # bingo!
    # This was in the original code, and needed vectorizing...
    # I think p[bingo] corresponds exactly (??) to mean(X)
    if( all( bingo)) {
return( 0*x + mean_x)
    }

    # change p for those cases to something that will converge and carry on
    # sub back correct x (ie mean) on exit
    p[ bingo] <- p[ !bingo][1]
    q0[ bingo] <- q0[ !bingo][1]
  }

  p_target <- p # so CSPA gives 0 at solution

  # Only risk I can see, is that if true s ~= 0, start may be on wrong side... and CSPA calcs go wrong.
  s0 <- 0.9 * q0 *isK2 # "zeroth-order" approx, times 0.9 for guess at lower bound
  while( any(
      bad <- !is.finite(
        C0 <- CSPA( s0)))) {
    s0[ bad] <- s0[ bad]/2
  }

  toobig <- C0 > 0
  multor <- 2 * xor( toobig, s0>0) - 1 # -1 or +1
  bracketed <- s0 != s0
  step <- 1.2 ^ multor
  snext <- s0

  repeat{
    snext[ !bracketed] <- s0[ !bracketed] * step[ !bracketed]
    bracketed[ !bracketed] <- xor( toobig, CSPA( snext) > 0)[ !bracketed]
    if( all( bracketed))
  break
    s0[ !bracketed] <- snext[ !bracketed]
  }

  hi <- pmax( s0, snext)
  lo <- pmin( s0, snext)

  s <- ridder( CSPA, lo, hi, tol=tol, skip_bounds=TRUE) # root finder
  CSPA( s)

  x[ bingo] <- mean_x # any that hit first time
return( x)
}


#' @importFrom mvbutils FOR %**%
"lglk_shebang" <- function( pars) {
#### Combined lglk for all ABCO-ambig genos *and* for total counts
#### Must have its environment reset before use (in geno6)...
  # ... to something already containing predigested data
  # ... and with CDF and PDF defined ("normal"ly pnorm and dnorm)

  # muhi = 2*mulo must be lower than hix... surely...
  # and data nobs, count, n_hi where...
  # ... count is a list, grouped by genotype but with all Hetz together
  # Very high scores just classed as "high" for outlier-robustness, grouped by genotype

  # The "real params", and useful stuff. Note double-arrow assignment
  mulo <<- (hix/2) * inv.logit( pars[1])
  siglo <<- exp( pars[2])
  pA <<- inv.logit( pars[3])
  pB <<- (1-pA) * inv.logit( pars[4])
  pC <<- if( length( pars)==5) (1-pA-pB)*inv.logit( pars[ 5]) else 0
  pO <<- max( 0, 1 - pA - pB - pC) # rounding error guard

  # HW totals
  phat <<- make_pgeno( pA, pB, pC, which_genotypes=genotypes) # in caller
  expected <- n_fish * phat[ genotypes_ambig]
  lglk_hw <- nobs %*% log(phat[ genotypes_ambig] + (nobs==0))        # avoid 0log0 gotcha
  # ... just gonna assume (count<OOthresh) <=> OOness exactly

  pen <- penscale * sum( log( cosh( tail( pars-pstart, -2)))) # tame pA/pB/pC only

  # Count scores
  muhi <<- 2*mulo
  sighi <<- sqrt(2) * siglo

  # For high scores:
  p_hi <- c( CDF( hix, mulo, siglo, lower=FALSE), CDF( hix, muhi, sighi, lower=FALSE))

  p12 <- with( as.list( phat), list( # all things are probs
    Hetz= c( 0, 1),
    AAO= c( AO, AA) / (AO+AA),
    BBO= c( BO, BB) / (BO+BB),
    CCO= if( pC>0) c( CO, CC) / (CO+CC) else c( 0.5, 0.5) # anti-NA
  ))

  # and <OOthresh assumed negligible for non-truly-OO animals

  dens <- FOR( names( p12), cbind( PDF( counts[[.]], mulo, siglo), PDF( counts[[.]], muhi, sighi)))

  lglk_counts <- 0
  for( g in names( p12)) {
    lglk_counts <- lglk_counts + sum( log( dens[[g]] %**% p12[[g]])) + n_hi[[g]] * log( p_hi %**% p12[[g]])
  }

return( lglk_counts + lglk_hw - pen)
}


#' @importFrom mvbutils %such.that% %not.in%
"load_several_whoppers" <- function(zipfile, exclude=character(),
                                    template_gobj, cutoff_npoly,
                                    max_fish_at_once=5000) {

  stopifnot( template_gobj %is.a% 'snpgeno') # could check presence of certain atts, too

  zippies <- unzip( zipfile, list=TRUE)$Name %such.that% (basename( .) %not.in% exclude)
  template_gobj <- template_gobj[0,]

  for( ifile in zippies) {
    scatn( 'Reading %s', ifile); flush.console()
    template_gobj <- load_whopper(
        csvfile= c( zipfile, ifile),
        gobj= template_gobj,
        cutoff_npoly=cutoff_npoly,
        max_fish_at_once=max_fish_at_once,
        dropped_fish_file= 'dropped-' %&% ifile)
  }

return( template_gobj)
}


#' @importFrom atease @ @<-
#' @importFrom utils unzip
"load_whopper" <- function(
    csvfile, # or length-2 where #1 is zipfilename, #2 is subfile
    gobj,
    max_fish_at_once=2000,
    cutoff_npoly= 0.5,
    dropped_fish_file= NULL) {
## 'gobj' must have been thru 'geno6way' (or 'load_whopper') already
## To discard existing fish before loading new ones: gobj=mygobj[0,]

  define_genotypes()
  locinfo <- gobj@locinfo
  n_loci <- nrow( locinfo)
  if( length( csvfile) > 1) { # presumably zipfile, then csvname
    unzippo <- csvfile[2]
    unzip( csvfile[1], unzippo) # NB extract into getwd() by default
    on.exit( unlink( unzippo))
    csvfile <- unzippo
  }

  n_new_fish <- defish_csv( csvfile, just_count_the_fish=TRUE)

  # for speed, work with stripped-down raw genotypes
  diplos <- gobj@diplos
  if( nrow( gobj)) {
    geno <- as.raw( gobj) # drop all atts
    dim( geno) <- dim( gobj)
  } else {
    geno <- matrix( as.raw( 0), 0,  n_loci)
  }

  # Grow the matrix only once, at the start
  n_preex <- nrow( geno)
  # Use of NA_integer_ works even if geno is 0-row
  geno <- rbind( geno, geno[rep(NA_integer_,n_new_fish),])

  n_read <- 0
  n_remaining <- n_new_fish

  tf <- tempfile() # will be re-used for each cut operation
  on.exit( unlink( tf))

  # Read all the fish-info
  n_finfo_lines <- 0
  repeat{
    first_line <- scan( csvfile, what='', sep='\n', n=1, skip=n_finfo_lines, quiet=TRUE)
    if( grepl( '^[*]?,', first_line)) { # star or blank
      n_finfo_lines <- n_finfo_lines+1
    } else
  break
  }

  new_finfo <- read_cluster_dart3( csvfile, n_lines_max=n_finfo_lines+3)@info # at least 2 fish

  scatn( 'checking locus info')
  defish_csv( csvfile, just_count_the_fish=FALSE, subset_fish=1:2, loci=TRUE, outfile=tf)

  # Cannot use 'read_cluster_dart3' since it switches the order of rows without telling...
  # so would be permuted WRTO direct scan later on
  # What we really need is the "Tag" field, so we can later look up the original "FullAlleleSeq" etc
  # ... in the list of Tags in the new file

  newloci <- scan( tf, skip=n_finfo_lines, what='', sep='\n', quiet=TRUE) # keep field-name line
  fields <- strsplit( newloci[1], ',')[[1]]

  i_tag_field <- match( 'Tag', fields, 0)
  if( !i_tag_field) {
stop( 'No "Tag" field..?')
  }

  if( i_tag_field>1) {
    new_alleles <- sub( sprintf( '([^,]*,){%i}([^,]*),.*', i_tag_field-1), '\\2', newloci[-1])
  } else {
    new_alleles <- sub( ',.*', '', newloci-1)
  }

  # Now match alleles...

  mmseq1 <- match( locinfo$FullRefSeq, new_alleles, 0)
stopifnot( all( mmseq1>0))

  mmseq2 <- match( locinfo$FullAltSeq, new_alleles, 0)
stopifnot( all( mmseq2>0))

  polyallelic <- !is.na( locinfo$ThirdAllele)
  mmseq3 <- match( locinfo$FullThirdSeq[ polyallelic], new_alleles, 0)
stopifnot( all( mmseq3>0)) # some loci don't have 3rd alleles

  n_new_seqs <- length( new_alleles)

  mean_fish_tot <- gobj@mean_fish_tot
  poly_fish <- rep( 0, n_preex + n_new_fish)
  rescalor <- locinfo$rescalor

  while( n_remaining) {
    this_subset <- seq( from=n_read+1, to=min( n_new_fish, n_read + max_fish_at_once))

    # Originally tried not to use 'cut' if all-in-one read
    # But, locus fields at row-start in orginal are not expected after cut
    # ... and cause trouble. So lways cut
#    if( tail( this_subset, 1) == n_new_fish) { # no need to cut
#      tf <- csvfile
#    } else { # snippity snippity snip

    scatn( 'Reading fish %i:%i', this_subset[ 1], tail( this_subset, 1))
    defish_csv( csvfile, just_count_the_fish=FALSE, subset_fish=this_subset, loci=FALSE, outfile=tf)

    n_fish_this_time <- length( this_subset)
    n_remaining <- n_remaining - n_fish_this_time
    n_read <- n_read + n_fish_this_time

    counts <- scan( tf, skip=n_finfo_lines+1, what=0L, sep=',', quiet=TRUE) # +1 cos of field-names line
    dim( counts) <- c( n_fish_this_time, n_new_seqs)

    # Put these counts into a temp F*L*3 array, like 'pick_ref_alt'
    # I hate vectoRization...

    c3 <- array( 0, c( n_fish_this_time, n_loci, 3)) # Ref, Alt, Others
    c3[ ,,1] <- counts[ ,mmseq1] # ... at last got this right ...
    c3[ ,,2] <- counts[ ,mmseq2 ]
    c3[ , polyallelic,3] <- counts[ ,mmseq3 ]

    # Renorm by fish...
    new_fish_tot_subset <- new_finfo$Fishtot[ this_subset]
    # vecless can't handle scalar multiply yet! so
    # c3[f,l,a] := c3[f,l,a] * mean_fish_tot / new_fish_tot_subset[ f]
    c3[f,l,a] := c3[f,l,a] / new_fish_tot_subset[ f]
    c3 <- c3 * mean_fish_tot

    # Rescale alleles, as in ppn_ref_alt_check2
    c3[f,l,a]:= c3[f,l,a] * rescalor[l,a]

    # Polyploid check--- though using only 3 chosen alleles
#    n_al_pres[ f, l]:= SUM_ %[a]% (c3[f,l,a] > 0.5 * cutoff_ppn * locinfo$med_het_tot[ l])
#    is_poly[ f, l]:= n_al_pres[ f, l] == 3
#    toopoly_fish[ n_preex + this_subset] <- rowSums( is_poly) > cutoff_npoly
#
    # Just do it...
    g6 <- geno6way( c3, locinfo, het_cut=gobj@het_cut)
stopifnot( my.all.equal( g6@diplos, diplos)) # paranoia...

    poly_fish[ n_preex + this_subset] <- rowSums( g6@is_polyzyg)
    geno[ n_preex + this_subset, ] <- g6
  }

  # Tack on attributes...
  attributes( geno) <- c( attributes( geno), attributes( gobj)[ cq( locinfo, diplos, mean_fish_tot, het_cut, calls, args)])
  # May be an old "npoly" field in fish-info--- drop, since it won't exist for new data
  # Drop any other info fields that aren't in both old and new
  old_finfo <- gobj@info[ names( gobj@info) %that.are.in% names( new_finfo)]
  new_finfo <- new_finfo[ names( new_finfo) %that.are.in% names( gobj@info)]
  geno@info <- rbind( old_finfo, new_finfo)
  geno@locinfo$npoly <- NULL # also meaningless
  oldClass( geno) <- 'snpgeno'

  toopoly <- poly_fish > cutoff_npoly
  if( length( as.character( dropped_fish_file)==1)) {
    dropped_fish <- geno@info$Our_sample[ toopoly] # really should have something to ID the *tube* kinda...
    cat( sprintf( '%s: %i', dropped_fish, poly_fish[ toopoly]), file=dropped_fish_file, sep='\n')
  }

  geno <- geno[ !toopoly,]
return( geno)
}


#' @importFrom handy2 sqr
#' @importFrom mvbutils FOR
"make_pgeno" <- function( pA, pB, pC, which_genotypes) {
###
# which_genotypes eg genotypes_ambig. If no C-genos requested then pC is forced to 0 so new O includes C
  if( !any( grepl( 'C', which_genotypes))) {
    pC <- 0
  }

  pO <- pmax( 0, 1 - pA - pB - pC) # rounding error guard
  pAO <- 2*pA*pO
  pAA <- sqr( pA)
  pAAO  <- pAA + pAO
  pAB  <- 2*pA*pB
  pBO <- 2*pB*pO
  pBB <- sqr( pB)
  pBBO  <- pBB + pBO
  pAC <- 2*pA*pC
  pCO <- 2*pC*pO
  pCC <- sqr( pC)
  pCCO <- pCC + pCO
  pBC <- 2*pB*pC
  pOO  <- pO^2

  # Could be scalar or vector: c or cbind
  funco <- if( length( pO) > 1) cbind else c
  phat <- do.call( funco, FOR( which_genotypes, get( 'p' %&% .)))
  names( phat) <- sub( '[.].*', '', names( phat)) # loco R name-extrusion habit FFS
return( phat)
}

#' @export
#' @importFrom atease @ @<-
#' @importFrom mvbutils %without.name%
"make_playback" <-
function( fhard, template, record_arg_name='record') {
  # template should come from calling 'fhard' with recording ON but with "normal" args

  original_fhard_env <- environment( fhard)
  recordo <- template@where

  # Check for prefixdims
  recorded_formals <- names( formals( fhard)) %that.are.in% names( recordo$subs)
  one_recorded_arg <- recorded_formals[1]

  # These need to exist when first called in playback-mode, so that <<- below will overwrite them
  recordo$prefixdims <- integer()
  recordo$counter <- NA # will get set

  e <- new.env( parent=original_fhard_env)
  environment( define) <- recordo # where counter & subs & prefixdim etc live
  e$define <- define
  e$`[.playback` <- `[.playback`
  e$`[<-.playback` <- `[<-.playback`
  recordo$ee <- e #

  vf <- fhard
  environment( vf) <- recordo
  formals( vf) <- formals( vf) %without.name% record_arg_name
  body( vf) <- substitute({
      new_dim <- dim( get( one_recorded_arg))
      old_dim <- subs[[ one_recorded_arg]][[1]]@orig_full_dim
      prefixdims <<- head( new_dim, -length( old_dim))
      if( !my.all.equal( tail( new_dim, length( old_dim)), old_dim)) {
    stop( sprintf( "Basic dimensions of '%s' have changed--- no can do!", one_recorded_arg))
      }
      # Could/should check whether ALL recorded formals have the same prefixdims, but it does mean force()ing the args
      # ... which might not be so bad...

      counter <<- 1+0*lengths( subs) # named counters

      # Arrange for args to be found, being paranoid to avoid name-clashes with things in ee
      argenv <- new.env( parent=ee)
      thisenv <- environment()
      for( iarg in names( formals( sys.function()))) {
        getme <- as.name( iarg)
        eval( substitute( delayedAssign( iarg, getme, eval.env=thisenv, assign.env=argenv), list( getme=getme)))
      }

      res <- eval( exprs, argenv)
      if( res %is.a% 'playback') { # almost always..?
        res <- as.vector( res) # bareass nekkid

        if( length( prefixdims) || (length( last_dim) > 1)) {
          dim( res) <- c( prefixdims, last_dim)
        } # else scalar is OK

        rldn <- last_dimnames
        if( length( rldn)) {
          if( length( prefixdims) || (length( rldn) > 1)) {
            dimnames( res) <- c( rep( list( NULL), length( prefixdims)), rldn)
          } else {
            names( res) <- rldn[[1]]
          }
        }
      }
    return( res)
    })

return( vf)
}


"make_tot2" <-
function( cA, cB, cC, geno) {
  define_genotypes()
  tot2 <- pmax( cA+cB, cA+cC, cB+cC) # keep OO "data"
  tot2[ geno==AAO] <- cA[ geno==AAO]
  tot2[ geno==BBO] <- cB[ geno==BBO]
  tot2[ geno==CCO] <- cC[ geno==CCO]
  tot2[ geno==AB] <- (cA + cB)[ geno==AB]
  tot2[ geno==AC] <- (cA + cC)[ geno==AC]
  tot2[ geno==BC] <- (cB + cC)[ geno==BC]
return( tot2)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils cq %without.name% %except%
#' @importFrom handy2 multimatch do.to2
"merge_new_samples" <- function(lociar, new_fish, drop_already){
  # rescale all counts and coerce to 3-allele format
  # could: zap poly fish
  # could: dump xshetz fish
  # apply geno6 rules

  n_loci <- ncol( lociar)
  li <- lociar@locinfo

  if( drop_already) {
    already <- multimatch( new_fish@info, lociar@info[ names( new_fish@info)], 0)
    new_fish <- new_fish[ !already, ]
  }

  n_new <- nrow( new_fish)

  # Drop discarded loci
  lox <- match( lociar@locinfo$CloneID, new_fish@locinfo$CloneID, 0)
stopifnot( all( lox)) # missing loci in new
  new_fish <- new_fish[ , lox]

  # Norm new fish
  nc <- unclass( new_fish)
  new_fish[] <- nc * lociar@mean_fish_tot / new_fish@info$Fishtot #

 # Rescale alleles-within-locus
  newals <- new_fish@seqinfo$FullAlleleSeq
  new_counts <- array( 0, c( n_new, n_loci, 3))
  new_counts[,,1] <- nc[,match( li$FullRefSeq, newals)]
  new_counts[,,2] <- nc[,match( li$FullAltSeq, newals)]

  hasC <- !is.na( li$FullThirdSeq)
  new_counts[,hasC,3] <- nc[,match( li$FullThirdSeq[ hasC], newals)]

  for( i in 1:3) {
    new_counts[,,i] <- new_counts[,,i] * rep( li$rescalor[,i], each=n_new)
  }

  # Genotype...
  mht <- li$med_het_tot
  rimht <- rep( 1/mht, each=n_new)
  new_geno_amb <- geno_deambig_ABC( new_counts,
        mht=mht,
        OOthresh_tc=li$OOthresh_tc,
        het_cut=lociar@het_cut,
        tc_hist_pars= NULL,
        ppnA_hist_pars= NULL,
        return_what= 'just_geno')
  new_geno6 <- split_geno_4to6( new_geno_amb, new_counts[,,1] * rimht, new_counts[,,2] * rimht, best_cut= li$best_cut)

  counts <- unclass( lociar)
  counts <- abind( counts, new_counts, along=1)

  ninfo <- new_fish@info
  miss_cols <- names( lociar@info) %except% names( ninfo)
  if( length( miss_cols)) {
    misso <- lociar@info[ rep(1,n_new), miss_cols,drop=FALSE]
    misso[ ,] <- NA
    ninfo <- cbind( ninfo, misso)
  }
  row.names( ninfo) <- row.names( new_fish@info) # since data.frame ops reliably **** 'em up

  updated_atts <- list(
    info= rbind( lociar@info, ninfo),
    geno6= rbind( lociar@geno6, new_geno6),
    geno_amb= rbind( lociar@geno_amb, new_geno_amb)
  )

  # Update counts and info parts of lociar, leaving all else unchanged
  other_atts <- attributes( lociar) %without.name% c( names( updated_atts), cq( dim, dimnames))
  nlociar <- do.call( 'structure', c( list( counts), updated_atts, other_atts))
return( nlociar)
}


#' @importFrom atease @ @<-
#' @importFrom handy2 multimatch
#' @importFrom mvbutils %except%
"old_pick_ref_alt" <- function(lociar) {
  li <- lociar@locinfo
  n_fish <- nrow( lociar)
  n_loci <- ncol( lociar)
  count_fs <- unclass( lociar)
  count_rao <- array( 0, c( n_fish, n_loci, 3)) # Ref, Alt, Others

  seqi <- lociar@seqinfo
  first_seq <- match( unique( seqi$Locus), seqi$Locus)

  # Find alleles with highest & 2nd-highest total counts
  # Seemingly, this has already been done (sorted into total-count order with locus)
  # but it's quick so why not

  ref <- do.to2( seqi, by=list( Locus), to=list( CountSum), fun='max')
  iseq <- multimatch( ref, seqi[ c( 'Locus', 'CountSum')])
  li$FullRefSeq <- seqi$FullAlleleSeq[ iseq]
  li$RefAllele <- seqi$Allele[ iseq]
  # li$RefNum <- iseq - first_seq + 1 # always 1 AFAICS
  count_rao[,,1] <- count_fs[ , iseq]

  # Now get 2nd-placer, by -ving the max...
  seqi$CountSum[ iseq] <- -seqi$CountSum[ iseq]
  alt <- do.to2( seqi, by=list( Locus), to=list( CountSum), fun='max')
  iseq <- multimatch( alt, seqi[ c( 'Locus', 'CountSum')])
  li$FullAltSeq <- seqi$FullAlleleSeq[ iseq]
  li$AltAllele <- seqi$Allele[ iseq]
  # li$AltNum <- iseq - first_seq + 1 # always 2 AFAICS
  count_rao[,,2] <- count_fs[ , iseq]

  # Now 3rd-placer, by also -ving the 2nd-placer...
  seqi$CountSum[ iseq] <- -seqi$CountSum[ iseq]
  alt <- do.to2( seqi, by=list( Locus), to=list( CountSum), fun='max')
  iseq <- multimatch( alt, seqi[ c( 'Locus', 'CountSum')])
  li$FullThirdSeq <- seqi$FullAlleleSeq[ iseq]
  li$ThirdAllele <- seqi$Allele[ iseq]
  count_rao[,,3] <- count_fs[ , iseq]
  count_rao[, li$n_alleles<3,3] <- 0
  li$FullThirdSeq[ li$n_alleles<3] <- NA
  li$ThirdAllele[ li$n_alleles<3] <- NA

  if( FALSE) { # used to have...
    # Count for Others := total_count - ref_count - alt_count
    # Hopefully efficient: add up all allele-1 counts, then allele-2 counts, ...

    tot <- count_fs[ , first_seq]

    for( ial in 2 %upto% max( li$n_alleles)) {
      these_loci <- which( li$n_alleles >= ial)
      tot[,these_loci] <- tot[,these_loci] + count_fs[ , first_seq[ these_loci]-1+ial]
    }

    count_rao[,,3] <- pmax( 0, tot - count_rao[,,1] - count_rao[,,2]) # rounding errors...
  }

  count_rao@locinfo <- li
  count_rao@info <- lociar@info
  other_atts <- names( attributes( lociar)) %except%
      c( names( attributes( count_rao)), 'class', 'seqinfo')
  attributes( count_rao) <- c( attributes( count_rao),
      attributes( lociar)[ other_atts]) # calls, args, etc
  oldClass( count_rao) <- 'loc.ar'
return( count_rao)
}


#' @importFrom atease @ @<-
"orig_geno_deambig_ABC" <-
function( lociar,
  mht= lociar@locinfo$med_het_tot,
  OOthresh_tc,
  het_cut,
  tc_hist_pars,
  ppnA_hist_pars,
  return_what=c( 'lociar', 'just_geno')) {
##############
# Pick OO, hetz
# Leave possible homoz as "possible" for now

  define_genotypes() # knows OO AB etc

  c3 <- unclass( lociar) # but can be passed in directly
stopifnot( length( dim( c3))==3 && dim( c3)[3]==3) # expect 3 allele counted (3rd might always be 0)

  # Zap irrel attribs for speed
  attributes( c3) <- list( dim=dim( c3), dimnames=list( NULL, NULL, cq( A, B, C)))
  n_loci <- ncol( c3)
  n_fish <- nrow( c3)

  mht <- lociar@locinfo$med_het_tot

  # Norm by locus
  c3[f,l,all] := c3[f,l,all] / mht[ l]
  cA <- c3[,,A]
  cB <- c3[,,B]
  cC <- c3[,,C]
  # guard against length-1 dims cRazies
  dim( cA) <- dim( cB) <- dim( cC) <- dim( c3)[1:2]

  # Double nulls:
  if( length( OOthresh_tc)==1) {
    OOthresh_tc <- rep( OOthresh_tc, n_loci)
  } # else already locus-specific, eg when regenotyping in geno6()

  is_OO[f,l] :=
      ((cA+cB)[f,l] < OOthresh_tc[l]) &
      ((cA+cC)[f,l] < OOthresh_tc[l]) &
      ((cB+cC)[f,l] < OOthresh_tc[l])

  switch( mode( tc_hist_pars),
    list = {
        # Default settings for 'hist' that *can* be overridden, but...
        tc_hist_pars <- add_list_defaults( tc_hist_pars,
            xlim=c( 0, 1.5),
            main='Tot c3 in FISHLOCI, normed so non-OO ~= 1',
            xlab='',
            nclass=200)
        # For some reason, directly using 'do.call( "hist"...)' is INCREDIBLY slow
        # Dunno why
        # do.call( 'hist', c( list( x=pmax( c( cA+cB), c( cA+cC), c( cB+cC))), tc_hist_pars))
        yy <- pmax( c( cA+cB), c( cA+cC), c( cB+cC))
        # do.call( 'hist', c( list( x=yy), tc_hist_pars))
        FFS_gimme_a_histo <- c( list( quote( hist), quote( yy)), tc_hist_pars)
        eval( as.call( FFS_gimme_a_histo))
        abline( v=mean( OOthresh_tc), col='red')
        abline( v=1, col='blue', lty=2, lwd=3)
      },
    expression = eval( tc_hist_pars),
    NULL = NULL
  )

  # Proportion of "A" c3--- only wanted for non-nulls, but easier to calc for all.
  ppnA <- cA / (cA + pmax( cB, cC))
  ppnB <- cB / (cB + pmax( cA, cC))
  ppnC <- cC / (cC + pmax( cA, cB))

  # optional graphics and/or user-specified outputs
  switch( mode( ppnA_hist_pars),
    list = {
        # Default settings for 'hist' that *can* be overridden, but...
        ppnA_hist_pars <- add_list_defaults( ppnA_hist_pars,
            main='Ppn Ref c3 in non-OO FISHLOCI, normed by locus',
            xlab='',
            nclass=200)
        # See previous histo comments
        #  do.call( 'hist', c( list( x=c( ppnA[ !is_OO])), ppnA_hist_pars))
        yy <- ppnA[ !is_OO]
        # do.call( 'hist', c( list( x=yy), tc_hist_pars))
        FFS_gimme_a_histo <- c( list( quote( hist), quote( yy)), ppnA_hist_pars)
        eval( as.call( FFS_gimme_a_histo))
        abline( v=c( het_cut, 1-het_cut), col='lightblue')
      },
    expression = eval( ppnA_hist_pars),
    NULL = NULL
  )

  inhetrange <- function( x) x %in.range% c( het_cut, 1-het_cut)
  is_AB <- !is_OO & inhetrange(ppnA) & inhetrange( ppnB)
  is_AC <- !is_OO & inhetrange(ppnA) & inhetrange( ppnC)
  is_BC <- !is_OO & inhetrange(ppnB) & inhetrange( ppnC)

  # Only one locus in jcslack4 actually seemed like a real problem here
  is_polyzyg <- is_AB + is_AC + is_BC > 1 # hopefully very few... in these cases, pick the biggest pair
  A_lowest <- cA < pmin( cB, cC)
  B_lowest <- cB < pmin( cA, cC)
  C_lowest <- cC < pmin( cA, cB)

  is_AB[ is_polyzyg & (A_lowest | B_lowest)] <- FALSE
  is_AC[ is_polyzyg & (A_lowest | C_lowest)] <- FALSE
  is_BC[ is_polyzyg & (B_lowest | C_lowest)] <- FALSE

  geno_amb <- matrix( AB, nrow( c3), ncol( c3),
      dimnames=list( c3@info$Our_sample, c3@locinfo$Locus))

  geno_amb[ is_OO] <- OO
  eps <- .Machine$double.eps   # since a feeeeeeew cases hit rounding error...
  geno_amb[ !is_OO & (ppnB+eps >= 1-het_cut)] <- BBO # for now; hetz will be fixed below
  geno_amb[ !is_OO & (ppnA+eps >= 1-het_cut)] <- AAO
  geno_amb[ !is_OO & (ppnC+eps >= 1-het_cut)] <- CCO

  geno_amb[ is_AB] <- AB
  geno_amb[ is_AC] <- AC
  geno_amb[ is_BC] <- BC

  class( geno_amb) <- c( 'noquote', oldClass( geno_amb)) # so shows ABO not "ABO" etc

  return_what <- match.arg( return_what)

  if( return_what=='just_geno') { # used when regenotyping during 'geno6'
return( geno_amb)
  }

  lociar@geno_amb <- geno_amb # but don't print it...
  lociar@subset_like_both <- c( lociar@subset_like_both, 'geno_amb')
  lociar@het_cut <- het_cut

  lociar@print <- expression({
      attr( x, 'locinfo') <- attr( x, 'seqinfo') <- attr( x, 'geno_amb') <- attr( x, 'print') <- NULL
      class( x) <- class( x) %except% 'specialprint'
      print( x)
      cat( '\n\n<<"locinfo" and "geno_amb" attributes hidden for brevity>>\n')
    })


return( lociar)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils do.on
"pick_FSPs_from_HSPs" <- function( snpg, HSPs) {
  # For pairs already picked as HSPs, ie PLOD(HSP,UP) > eta: they might be FSPs. H
  # How would a FSP / HSP comparison look?
  # Don't need full pairwise screening for FSPs (do post hoc on a few hundred HSPs), hence all in R.

  # HSPs should be M*2 matrix of rows in snpg that are HSPs or FSPs

  # Transform to 4way genotypes
  # based on code in find_duplicates
  # careful, since "factor level" of AB and OO is different in 4way vs 6way
  sibg <- just_sibg <- snpg[ c( HSPs),]

  sibg@diplos <- genotypes4_ambig
  sibg[ just_sibg==AO] <- AAO
  sibg[ just_sibg==AA] <- AAO
  sibg[ just_sibg==BO] <- BBO
  sibg[ just_sibg==BB] <- BBO
  sibg[ just_sibg==OO] <- OO # need to do OO & AB too, since codes are different in 4way vs 6way
  sibg[ just_sibg==AB] <- AB


  extract.named( snpg@locinfo[ cq( PUP4, LOD4)])
  PHSP4 <- exp( LOD4) / PUP4 # Pr[gg|HSP] <- 0.5 * PUP4 + 0.5 * Pr[gg|kappa=1]
  P_k1 <- 2*PHSP4 - PUP4
  P_k1[ P_k1 < 0] <- 0 # rounding error

  samo <- do.on( strsplit( colnames( LOD4), '/'), .[1]==.[2])
  PFSP4 <- ...
  # Number of ID locis should be

  nsibs <- nrow( HSP)

  n_loci_id <- sibg[ 1:nsib,] == sibg[ nsib + 1:nsib,]
  }


#' @importFrom atease @ @<-
#' @importFrom mvbutils %except%
"pick_ref_alt" <- function( lociar) {
  li <- lociar@locinfo
  n_fish <- nrow( lociar)
  n_loci <- ncol( lociar)
  count_fs <- unclass( lociar)
  count_rao <- array( 0, c( n_fish, n_loci, 3)) # Ref, Alt, Others

  seqi <- lociar@seqinfo
  first_seq <- match( unique( seqi$Locus), seqi$Locus)

  # Find alleles with highest & 2nd-highest total counts
  # Seemingly, this has already been done (sorted into total-count order with locus)
  # but it's quick so why not

  uclu <- with( seqi, match( Locus, unique( Locus)))
  o <- order( uclu, -seqi$count_sum) # NB ascending/descending; "count_sum" in case Dart's "CountSum" is missing
  if( !all( o==seq_along( o))) {
    seqi <- seqi[ o,]
    count_fs <- count_fs[ ,o]
  }

  first_loc <- match( li$Locus, seqi$Locus)
  li$FullRefSeq <- seqi$FullAlleleSeq[ first_loc] # match
  li$RefAllele <- seqi$Allele[ first_loc]
  count_rao[,,1] <- count_fs[,first_loc]

  li$FullAltSeq <- seqi$FullAlleleSeq[ first_loc+1]
  li$AltAllele <- seqi$Allele[ first_loc+1]
  count_rao[,,2] <- count_fs[,first_loc+1]

  li$FullThirdSeq <- li$ThirdAllele <- NA_character_
  has2plus <- li$n_alleles > 2
  li$FullThirdSeq[ has2plus] <- seqi$FullAlleleSeq[ first_loc[ has2plus]+2]
  li$ThirdAllele[ has2plus] <- seqi$Allele[ first_loc[ has2plus]+2]
  count_rao[,has2plus,3] <- count_fs[,first_loc[has2plus]+2]

  count_rao@locinfo <- li
  count_rao@info <- lociar@info
  other_atts <- names( attributes( lociar)) %except%
      c( names( attributes( count_rao)), 'class', 'seqinfo')
  attributes( count_rao) <- c( attributes( count_rao),
      attributes( lociar)[ other_atts]) # calls, args, etc
  oldClass( count_rao) <- 'loc.ar'
return( count_rao)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils %is.not.an%
"playback" <- function(P, snerr, recordo, prefixdims=numeric()){
  if( recordo %is.not.an% 'environment') {
    recordo <- recordo@where
  }
  recordo$prefixdims <- prefixdims
  recordo$counter <- 1+0*lengths( recordo$subs) # named
  environment( define) <- recordo # stuff it needs
  # eval.parent( recordo$exprs)
  res <- eval( recordo$exprs)
  if( res %is.a% 'playback') { # almost always..?
    res <- as.vector( res) # bareass nekkid

    if( length( prefixdims) || (length( recordo$last_dim) > 1)) {
      dim( res) <- c( prefixdims, recordo$last_dim)
    } # else scalar is OK

    rldn <- recordo$last_dimnames
    if( length( rldn)) {
      if( length( prefixdims) || (length( rldn) > 1)) {
        dimnames( res) <- c( rep( list( NULL), length( prefixdims)), rldn)
      } else {
        names( res) <- rldn[[1]]
      }
    }
  }
return( res)
}


#' @importFrom atease @ @<-
#' @importFrom graphics title abline hist
#' @importFrom mvbutils %except%
"ppn_ref_alt_check2" <- function( lociar, OOthresh_tc,
    OK_med_min, # keep loci (and keep 3rd allele if ditto)
    selector, rescale_alleles, # T/F
    plot.=FALSE) {
### 3rd count is a real sequence BTW, but we might here decide not to trust it here

### Norm by locus
  n_fish <- nrow( lociar)
  n_loci <- ncol( lociar)

  mht <- lociar@locinfo$med_het_tot
  c3 <- unclass( lociar) # Should be
stopifnot( length( dim( c3))==3 && dim( c3)[3]==3) # expect 3 allele counted (3rd might always be 0)
  attributes( c3) <- list( dim=dim( c3), dimnames=list( NULL, NULL, cq( A, B, C)))

  c3[f,l,all] := c3[f,l,all] / mht[ l]

  OK_med_range <- sort( c( OK_med_min, 1-OK_med_min))
  plotto <- function( med_hetz, allele_type) {
    hist( med_hetz, nclass=70, xlim=0:1, xlab=NULL, ylab=NULL,
        main=sprintf( 'Med Ppn %s in plaus hetz', allele_type))
    abline( v=OK_med_range, col='lightgreen', lwd=2)
    abline( v=c( selector, 1-selector), col='orange')
    title( main=sprintf( '%i%% of loci OOR', as.integer( 100*mean( ! ((med_hetz %except% NA) %in.range% OK_med_range)))), line=-1)
  }

  find_med_hetz <- function( me, other) {
    # Keep only "reasonable" ones, and not OOs
      ppnal <- me / (me+other)
      is_OO <- (me+other) < OOthresh_tc
      ppnal[ (ppnal < selector) | (ppnal > 1-selector)] <- NA
      ppnal[ is_OO] <- NA
    return( apply( ppnal, 2, median, na.rm=TRUE))
    }

  med_ppnA_hetz <- find_med_hetz( c3[,,'A'], c3[,,'B']) # just wanna make sure Ref & Alt are OK
  if( plot.) {
    plotto( med_ppnA_hetz, 'Ref')
  }

  med_ppnA_hetz[ is.na( med_ppnA_hetz)] <- (-1) # to drop NAs, which can occur if almost all OO and no hetz!
  keepo <- med_ppnA_hetz %in.range% OK_med_range
  lociar@locinfo$med_ppnA_hetz <- med_ppnA_hetz

  if( rescale_alleles) { # want approx same score for Ref and Alt (and 3rd) if present--- based on medians
    # Preserve total count A+B
    totAB[loc,al] := c3[+.,loc,al=cq(A,B)]
    tot2[loc]:= totAB[loc,+.]
    rescalor <- cbind( A= 1/med_ppnA_hetz - 1, B= 1)

    # If we just Xed the count by rescalor, then med_ppnA_hetz would become 0.5
    # But total A+B count would change, so re-rescale
    temp[loc]:= totAB[loc,al] %[al]% rescalor[loc,al]
    rescalor[loc,al] := rescalor[loc,al]  * (tot2[loc] / temp[loc])
    rescalor <- cbind( rescalor, C=1) # for now, to facilitate next line
    c3[f,l,al]:= c3[f,l,al] * rescalor[l,al]
  }

  # Now 3rd allele C (after adjusting A vs B)
  pmaxAB <- pmax( c3[,,'A'], c3[,,'B'])
  dim( pmaxAB) <- dim( c3)[1:2]
  med_ppnC_hetz <- find_med_hetz( c3[,,'C'], pmaxAB) # slightly different criterion

  # Kick out 3rd allele if NBG
  bad3 <- is.na( med_ppnC_hetz) | !(med_ppnC_hetz %in.range% OK_med_range)
  lociar[,bad3,3] <- c3[,bad3,'C'] <- 0
  lociar@locinfo$FullThirdSeq[ bad3] <- NA
  lociar@locinfo$ThirdAllele[ bad3] <- NA

  if( plot.) {
    plotto( med_ppnC_hetz[ keepo], '3rd')
  }

  if( rescale_alleles) {
    # For Cs, need to use all the hetz we can find. Redo with new scores
    med_ppnC_hetz[ bad3] <- 0.5 # ensures no adjustment to A or B for these loci!

    # Preserve total count A+B
    totCX <- colSums( c3)
    tot3[loc]:= totCX[loc,+.]
    orig_rescalor <- rescalor

    rescalor <- cbind( A=1, B=1, C=1/med_ppnC_hetz-1)

    # If we just Xed the count by rescalor, then med_ppnC_hetz would become 0.5
    # But total count would change, so re-rescale
    # Questionable whether to base this on hetz only...

    temp[loc]:= totCX[loc,al] %[al]% rescalor[loc,al]
    rescalor[loc,al] := rescalor[loc,al]  * (tot3[loc] / temp[loc])
    c3[f,l,al]:= c3[f,l,al] * rescalor[l,al]

    rescalor[l,al] := rescalor[l,al] * orig_rescalor[l,al] # for 1-step application in fture

    c3[f,l,al] := c3[f,l,al] * mht[ l]
    lociar[] <- c3
  } else { # No rescaling, but still rescalor needs to exist
    rescalor <- matrix( 1, n_loci, 3)
  }

  lociar@locinfo$rescalor <- rescalor # for adding more fish later
  lociar <- lociar[ , keepo]
return( lociar)
}


#' @importFrom atease @ @<-
#' @importFrom handy2 sqr
"predict_hsp_util" <-
function( pIBD0, pIBD1, want_LOD_table=FALSE, k=0.5) {
  # This version ignores the possibility of errors involving AB or OO...
  # ... which should be pretty rare

  define_genotypes()
  nl <- nrow( pIBD1)
  Phsp <- pIBD1 * k + pIBD0 * (1-k)
  Pup <- pIBD0

  LOD <- log( Phsp / Pup)
  LOD[ Pup==0] <- 0 # if Pup=0 then p*log(p) = 0; only happens when r=0
  if( want_LOD_table) {
    # LOD is 3D: nloci * ng1 * ng2
    # gpLOD is 2D: nloci * n_genopairs
    # Need only certain "columns" of 2D-fied LOD
    mg <- make_genopairer( dimnames( pIBD0)[[2]])
    ngp <- max( mg)
    wanted <- match( 1:ngp, mg)

    gpLOD <- gpPUP <- matrix( 0, nl, ngp)
    LOD_as_2D <- matrix( LOD, nl, prod( dim( pIBD0)[-1]))
    gpLOD <- LOD_as_2D[ , wanted]
    PUP_as_2D <- matrix( Pup, nl, prod( dim( pIBD0)[-1]))
    gpPUP <- PUP_as_2D[ , wanted]

    # Off-diagnonals appear twice, and prob should be doubled...
    omg <- mg
    omg[ wanted] <- 0
    double_wanted <- 1:ngp %in% omg
    # wrong for some reason:    double_wanted <- wanted %in% mg[ duplicated( c( mg))]
    gpPUP[ ,double_wanted] <- gpPUP[,double_wanted] * 2

    dimnames( gpLOD) <- dimnames( gpPUP) <- list( dimnames( pIBD0)[[1]], mg@what)
    gpLOD@mg <- mg # why not
  }

  # EPLOD is sum( LOD * Pup) but we want to keep it by locus for now
  E.HSP[l] := LOD[l,i,j] %[i,j]% Phsp[l,i,j]
  E.UP[l] := LOD[l,i,j] %[i,j]% Pup[l,i,j]
  E2.UP[l] := (LOD*LOD)[l,i,j] %[i,j]% Pup[l,i,j]
  V.UP <- E2.UP - sqr( E.UP)
  Ediff <- E.HSP - E.UP

  # Standardized difference ie locus power: not so useful post hoc, but possibly interesting for 6 vs 4 comps
  sdiff <- (E.HSP - E.UP) / sqrt( V.UP)

#  Ediff <- unclass( Ediff)
#  V.UP <- unclass( V.UP)
#  sdiff <- unclass( sdiff)
#

  retval <- data.frame( Ediff, V.UP, sdiff)
  if( want_LOD_table) {
    retval@LOD <- gpLOD
    retval@PUP <- gpPUP
  }
return( retval)
}

#' @importFrom mvbutils cq %except% %not.in%
#' @importFrom atease @ @<-
#' @importFrom handy2 sqr
"prepare_PLOD_SPA" <- function( geno6, n_pts_SPA_renorm=201) {
# To be run after hsp_power( ..., want_LOD_table=TRUE)
# n_pts_SPA_renorm should really be as big as R can handle without running out memory
# ... but 201 should be OK I guess. If 201 and 301 give almost-identical results then all well!

stopifnot( all( cq( LOD4, LOD6, use6) %in% names( geno6@locinfo)))

  og <- options( vecless.print=FALSE)
  on.exit( options( og))

  # Combine 4way and 6way stuff into overall LOD and PUP
  extract.named( geno6@locinfo[ cq( use6, LOD6, LOD4, PUP6, PUP4)])
  use4 <- !use6

  # DO NOT change actual genos though; they will be changed on-the-fly prior to kin-finding
  # ... code WOULD be this:
  # temp_snpg <- snpg
  # recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  # temp_snpg[ , use4] <- recode4to6temp( snpg[, use4]) # (AA,AO) -> AA; (BB,BO) -> BB

  LOD <- LOD6
  PUP <- PUP6
  cn6 <- colnames( LOD6)
  cn4 <- colnames( LOD4)

  # Change only the entries with "full homz" since "single nulls" won't be accessed
  for( ichangio in grep( 'AA|BB', cn6) %except% grep( 'AO|BO', cn6)){
    was1 <- substring( cn6[ ichangio], 1, 2)
    was2 <- substring( cn6[ ichangio], 4, 5)
    now1 <- sub( '(A|B)\\1', '\\1\\1O', was1)
    now2 <- sub( '(A|B)\\1', '\\1\\1O', was2)
    iget4 <- paste( now1, now2, sep='/')
    if( iget4 %not.in% cn4) { # reverse the order
      iget4 <- paste( now2, now1, sep='/')
    }

    LOD[ use4, ichangio] <- LOD4[ use4, iget4]
    PUP[ use4, ichangio] <- PUP4[ use4, iget4]
  }
  # For safety's sake, LOD( XO,...) := NA; should never be accessed
  hasO <- grep( '(A|B)O', cn6)
  LOD[ use4, hasO] <- NA # security in case of wrong access later for real data
  PUP[ use4, grep( '(A|B)O', cn6)] <- 0
  LOD@mg <- make_genopairer( geno6@diplos)

  make_K <- function( PUP, LOD) { # ... while the sun skines

      e <- new.env( parent=asNamespace( 'vecless'))
      e$PUP <- PUP
      e$LOD <- e$LODOK <- LOD
      e$LODOK[ is.na( LOD)] <- 0 # leaving NAs in would mess up the calcs
      e$n_pts_SPA_renorm <- n_pts_SPA_renorm

      evalq( envir=e, {
        if( !nrow( PUP)) {
          K <- dK <- ddK <- function( tt) 0*tt
          inv_CDF <- CDF <- function( x) NA+x
        } else {
          PUPLOD <- PUP * LODOK
          PUPLOD2 <- PUPLOD * LODOK

          K <- function( tt) {
            ETT[ it, l, g12] := exp( tt[ it] * LODOK[ l, g12])
            S[ it, l] := PUP[ l, g12] %[g12]% ETT[ it, l, g12]
            rowSums( log( S))
          }

          dK <- function( tt) {
            ETT[ it, l, g12] := exp( tt[ it] * LODOK[ l, g12])
            S[ it, l] := PUP[ l, g12] %[g12]% ETT[ it, l, g12]
            SL[ it, l] := PUPLOD[ l, g12] %[g12]% ETT[ it, l, g12]
            rowSums( SL/S)
          }

          ddK <- function( tt) {
            ETT[ it, l, g12] := exp( tt[ it] * LODOK[ l, g12])
            S[ it, l] := PUP[ l, g12] %[g12]% ETT[ it, l, g12]
            SL[ it, l] := PUPLOD[ l, g12] %[g12]% ETT[ it, l, g12]
            SLL[ it, l] := PUPLOD2[ l, g12] %[g12]% ETT[ it, l, g12]
            rowSums( (SLL/S-sqr( SL/S)))
          }

          extract.named( renorm_SPA_cumul( K, dK, ddK, n_pts=n_pts_SPA_renorm))
          # ... CDF and inv_CDF but *not* for extreeeme values
        } # if nrow PUP
      })

      # K and co will know PUP & co thru enviro magic
    return( e)
    } # function make_K


  Kenv <- make_K( PUP, LOD)
  geno6@Kenv <- Kenv

  # PUP and LOD are in Kenv now, so don't duplicate them in locinfo
  # They are really the "workhorse" versions, and are a bit cheaty, so don't want
  # them too public

return( geno6)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils cq do.on %except%
"print.NGS_count_ar" <- function( x,
    trailing_dot=getOption( 'trailing_dot_NGS_count_ar', FALSE),
    dot_for_0=getOption( 'dot_for_zero_NGS_count_ar', FALSE), ...){

  attx <- attributes( x)
  df <- x@info
  rownames( df) <- NULL

  x <- unclass( x)
  n_fish <- nrow( x)
  extract.named( x@locinfo)
  n_loci <- length( start_col)

  # If non-integer, can either use trailing dot, or merely round

  if( is.integer( x)) {
    mi <- as.character( x)
  } else {
    mi <- as.character( round( x))
    if( trailing_dot) {
      mi <- paste0( mi, '.')
    }
  }

  dot <- rawToChar( as.raw( 183))
  if( dot_for_0) {
    zero <- if( trailing_dot) '0.' else '0'
    mi[ mi==zero] <- dot
  }
  dim( mi) <- dim( x)

  m <- do.on( 1:n_loci, {
      z <- mi[, start_col[ .]:end_col[ .], drop=FALSE]
      fw <- max( nchar( z))
      z[] <- format( z, width=fw, justify='right')
      z[] <- gsub( ' ', dot, z, fixed=TRUE)
      sprintf( '%s  ', apply( z, 1, paste, collapse='/'))
    })
  dim( m) <- c( n_fish, n_loci) # in case only 1 fish or 1 locus, whereupon do.on will strip dim
  # dimnames( m)[[2]] <- Locus # FFS R rejects this if m is N*1 ... so
  dimnames( m) <- list( NULL, Locus) # a compulsory column, albeit "made-up by Mark" sometimes

  m[ grepl( 'NA', m, fixed=TRUE)] <- 'NA' # presumably only if binding several datasets?

  # x <- cbind( df, m) # fucken cbind turns strings into fucken factors...
  x <- cbind( df, data.frame( m, stringsAsFactors=FALSE)) # a data.frame FFS
  print( x, quote=FALSE, right=TRUE, ...)

  # Extra atts will not be printed by default (and print.data.frame ignores them)
  # Set 'x@print_atts' to character to
  dont_print_atts <- names( attx) %except% c( attx$print_atts,
      cq( class, dim, dimnames, print_atts, info, locinfo, seqinfo))
  if( length( dont_print_atts)) {
    scatn( '\nOmitting these attributes: %s', paste( dont_print_atts, collapse=', '))
  }
  for( att_to_print in attx$print_atts) {
    scatn( '\nattr( ., %s):', att_to_print)
    print( attx[[ att_to_print]], ...)
  }

  if( 'print_atts' %in% names( attx)) {


  }

  invisible( x)
}


#' @importFrom atease @ @<-
"re_est_ALF" <-
function( snpg) {
## check to be called after load_whopper loads entire dataset

  define_genotypes()
  n_samp <- nrow( snpg)
  n_loci <- ncol( snpg)
  gamb <- snpgeno( n_samp, n_loci, diplos=genotypes_ambig) # includes C but won't be used
  gamb[ snpg==AB] <- AB
  gamb[ snpg==OO] <- OO
  gamb[ snpg==AO] <- AAO
  gamb[ snpg==AA] <- AAO
  gamb[ snpg==BO] <- BBO
  gamb[ snpg==BB] <- BBO

  snpg@geno_amb <- gamb # required by...
  new_ALFs <- est_ALF_ABCO( snpg)
return( new_ALFs)
}

#' @importFrom openxlsx read.xlsx
#' @importFrom gbasics loc.ar
#' @importFrom atease @ @<-
#' @importFrom mvbutils FOR
#' @importFrom utils head
"read_cluster_dart2" <- function(filename, filtered=FALSE, use_rownames=FALSE,
                                 show_progress=TRUE) {
  # CSV file with a few header lines of plate info, then (unlike
  # read_count_dart) each CLUSTER in several rows, one for each SNP sequence:
  # still needs to be reorganized (inefficiently) into pairs
  # ... first few cols are locus summaries, remainder are fishwise counts
  # Use either raw or filtered counts (before or after 2nd ClusterIdx column)
  # use_rownames default changed to FALSE

  if( grepl( '(?i)[.]xls(x)?$', filename)) {
    temp <- openxlsx::read.xlsx( filename, colNames=FALSE)
    # No doubt could get stuff directly from data.frame 'temp', but simpler
    # to use existing code
    # stuffo <- apply( temp, 1, paste, collapse=',')
    # ****ing NA is mishandled by paste(), so...
    temp <- FOR( temp, {z <- as.character(.); z[is.na(.)] <- ''; z})
    temp$sep <- ',' # easiest
    stuffo <- do.call( 'paste', temp)
  } else {
    stuffo <- scan( filename, what='', sep='\n')
  }

  # For some reason, there are completely blank rows in one spreadsheet
  stuffo <- stuffo[ grepl( '[^,]', stuffo)]

  # Check format (standard or JC Nov 2016)
  if( !length( grep( 'ClusterIdx,', stuffo, fixed=TRUE))) {
return( read_count_dart( filename, use_rownames=use_rownames, stuffo=stuffo))
  }

  well.row <- grep( '(,[A-H][0-9]{0,2})+$', stuffo)
  #  well.row <- grep( 'Extract_well', stuffo, value=TRUE) not always labelled
stopifnot( length( well.row)==1)

  wells <- strsplit( stuffo[ well.row], ',')[[1]]
  is.fish.col <- grepl( '^[A-H][0-9]{1,2}$', wells)

  # raw or filtered

  gap.col <- max( which( !is.fish.col))
  if( any( which( is.fish.col) < gap.col)) { # then we do have to decide
    nfish <- sum( is.fish.col) / 2
    nprefishcols <- sum( !is.fish.col)-1
stopifnot( gap.col == nfish + nprefishcols + 1)
    if( filtered) {
      # drop RH end
      stuffo <- sub( sprintf( '^([^,]*)((?:,[^,]*){%i}).*', nfish+nprefishcols-1), '\\1\\2', stuffo, perl=TRUE)
      dont_drop_yet <- 1 %upto% (nprefishcols + nfish)
    } else {
      # drop middle
      stuffo <- sub( sprintf( '^((?:[^,]*,){%i})([^,]*,){%i}', nprefishcols, nfish+1), '\\1', stuffo, perl=TRUE)
      dont_drop_yet <- c( 1 %upto% nprefishcols, nfish + nprefishcols + 1 + (1 %upto% nfish))
    }

    wells <- wells[ dont_drop_yet]
    is.fish.col <- is.fish.col[ dont_drop_yet]
  } else { # only one version of counts
    nfish <- sum( is.fish.col)
  }

  star.start <- substring( stuffo, 1, 1)=='*'
  droppo <- which.max( cumprod( star.start) * cumsum( star.start))
  metafish <- stuffo[ (1 %upto% droppo)] # %except% well.row] now keeping that
  stuffo <- stuffo[ -(1 %upto% droppo)]

  chardat <- strsplit( stuffo, ',', fixed=TRUE)
stopifnot( length( unique( lengths( chardat))) == 1)

  cols <- chardat[[1]]
  ncols <- length( cols)
  data <- matrix( unlist( chardat[-1]), ncol=ncols, byrow=TRUE)
  rm( chardat)

  seqinfo <- data[,!is.fish.col,drop=FALSE]
#   fishinfo <- data[,is.fish.col,drop=FALSE]
  fishinfo <- matrix( as.integer( data[,is.fish.col]), ncol=sum( is.fish.col))
  rm( data)

  deduce_dart_header_rows()

  # Loci: check numeric convertibelity by column
  # This doesn't *guarantee* numericism (eg 1..-2 is not convertible) but FFS...
  unint.char <- grepl( '[^0-9-]', seqinfo)
  unreal.char <- grepl( '[^0-9.eE+-]', seqinfo)
  dim( unint.char) <- dim( unreal.char) <- dim( seqinfo)

  intable <- colSums( unint.char)==0
  realable <- (colSums( unreal.char)==0) & !intable

  seqinfo <- data.frame( seqinfo, stringsAsFactors=FALSE)
  names( seqinfo) <- cols[ !is.fish.col]
  seqinfo[ intable] <- lapply( seqinfo[ intable], as.integer)
  seqinfo[ realable] <- lapply( seqinfo[ realable], as.numeric)
  seqinfo <- seqinfo[ !duplicated( cols[ !is.fish.col])] # ClusterIdx can appear twice (well 3 times overall);
  # ... name will be changed by names<- above
  seqinfo$count_sum <- rowSums( fishinfo)

  dupseq <- duplicated( seqinfo$Tag) # shouldn't happen but does. This format is iffy...
  # At least the entire row seems duplicated
  fishinfo <- fishinfo[ !dupseq,] # not yet transposed
  seqinfo <- seqinfo[ !dupseq,]

  # Make sure "clusters" are contiguous, with commonest "ref" sequence first
  uclu <- with( seqinfo, match( ClusterIdx, unique( ClusterIdx)))
  o <- order( uclu, -seqinfo$count_sum) # NB ascending/descending; "count_sum" in case Dart's "CountSum" is missing
  if( !all( o==seq_along( o))) {
    seqinfo <- seqinfo[ o,]
  }

  seqinfo <- within( seqinfo, {
    CloneID <- sprintf( 'L%i', ClusterIdx)
    Allele <- '' # abbreviated version, filled in later
  })
  names( seqinfo) <- sub( '^Tag$', 'FullAlleleSeq', names( seqinfo))

  tot_alleles <- nrow( seqinfo)
  n_alleles <- table( seqinfo$ClusterIdx)
  locinfo <- data.frame(
      CloneID= sprintf( 'L%s', names( n_alleles)),
      ClusterTempIndex= names( n_alleles),
      n_alleles= c( n_alleles),
      consensus='',
      var_pos='',
      row.names=NULL, stringsAsFactors=FALSE)
  locinfo$end_col <- cumsum( n_alleles)
  locinfo$start_col <- c( 1L, head( locinfo$end_col, -1)+1)


  for( iclu in seq_along( n_alleles)) {
    if( show_progress && (iclu %% 100 == 0)) {
      cat( '\r', sprintf( '%i of %i clusters', iclu, length( n_alleles)))
    }
    which_seqs <- which( seqinfo$ClusterIdx == names( n_alleles)[ iclu])
    tags <- seqinfo$FullAlleleSeq[ which_seqs]
    raws <- matrix( charToRaw( paste( tags, collapse='')), ncol=length( tags))
    is_diff <- raws != raws[,1]
    has_var <- rowSums( is_diff) > 0
    var_pos <- which( has_var)
    locinfo$var_pos[ iclu] <- paste( var_pos, collapse=',')
    conseq <- raws[,1]
    conseq[ has_var] <- charToRaw( '.')
    locinfo$consensus[ iclu] <- rawToChar( conseq)
    seqinfo$Allele[ which_seqs] <- sprintf( 'L%s_%s',
        names( n_alleles)[iclu],
        apply( raws[ var_pos,,drop=FALSE], 2, rawToChar))
  }
  if( show_progress) {
    cat( '\n')
  }

  # For compatibility with other funcs that expect ClusterTempIndex column
  names( seqinfo) <- sub( 'ClusterIdx', 'ClusterTempIndex', names( seqinfo), fixed=TRUE)
  names( locinfo) <- sub( 'ClusterIdx', 'ClusterTempIndex', names( locinfo), fixed=TRUE)

  # Construct 'NGS_count_ar' object
  dim( fishinfo) <- c( tot_alleles, nfish)
  fishinfo <- fishinfo[ o, ] # match changes in seqinfo
  fishinfo <- t( fishinfo)
  dimnames( fishinfo) <- list( Well=wells[ is.fish.col], NULL)
  class( fishinfo) <- c( 'NGS_count_ar', class( fishinfo))

  fishinfo@info <- metafish
  fishinfo@locinfo <- locinfo
  fishinfo@seqinfo <- seqinfo

return( fishinfo)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils cq %SUCH.THAT% %is.not.a% %except% FOR
#' @importFrom gbasics loc.ar
#' @importFrom openxlsx read.xlsx
#' @importFrom utils head
"read_cluster_dart3" <-
function( filename, filtered=FALSE, use_rownames=FALSE, show_progress=TRUE, n_lines_max=-1){
  # CSV file with a few header lines of plate info, then (unlike read_count_dart) each CLUSTER in several rows, one for each SNP sequence:
  # still needs to be reorganized (inefficiently) into pairs
  # ... first few cols are locus summaries, remainder are fishwise counts
  # Use either raw or filtered counts (before or after 2nd ClusterIdx column)
  # use_rownames default changed to FALSE

  if( grepl( '(?i)[.]xls(x)?$', filename)) {
    # n_lines_max is ignored here...
    temp <- openxlsx::read.xlsx( filename, colNames=FALSE)
    # No doubt could get stuff directly from data.frame 'temp', but simpler to use existing code
    # stuffo <- apply( temp, 1, paste, collapse=',')
    # ****ing NA is mishandled by paste(), so...
    temp <- FOR( temp, {z <- as.character(.); z[is.na(.)] <- ''; z})
    temp$sep <- ',' # easiest
    stuffo <- do.call( 'paste', temp)
  } else {
    stuffo <- scan( filename, what='', sep='\n', n=n_lines_max)
  }

  # For some reason, there are completely blank rows in one spreadsheet
  stuffo <- stuffo[ grepl( '[^,]', stuffo)]

  # Check format (standard or JC Nov 2016)
  if( !length( grep( 'ClusterIdx,', stuffo, fixed=TRUE))) {
#return( read_count_dart( filename, use_rownames=use_rownames, stuffo=stuffo))
stop( "Looks like older, or 'new improved' format; try 'read_count_dart'")
  }

  well.row <- grep( '(,[A-H][0-9]{0,2})+$', stuffo)
  #  well.row <- grep( 'Extract_well', stuffo, value=TRUE) not always labelled
stopifnot( length( well.row)==1)

  wells <- strsplit( stuffo[ well.row], ',')[[1]]
  is.fish.col <- grepl( '^[A-H][0-9]{1,2}$', wells)

  # raw or filtered

  gap.col <- max( which( !is.fish.col))
  if( any( which( is.fish.col) < gap.col)) { # then we do have to decide
    nfish <- sum( is.fish.col) / 2
    nprefishcols <- sum( !is.fish.col)-1
stopifnot( gap.col == nfish + nprefishcols + 1)
    if( filtered) {
      # drop RH end
      stuffo <- sub( sprintf( '^([^,]*)((?:,[^,]*){%i}).*', nfish+nprefishcols-1), '\\1\\2', stuffo, perl=TRUE)
      dont_drop_yet <- 1 %upto% (nprefishcols + nfish)
    } else {
      # drop middle
      stuffo <- sub( sprintf( '^((?:[^,]*,){%i})([^,]*,){%i}', nprefishcols, nfish+1), '\\1', stuffo, perl=TRUE)
      dont_drop_yet <- c( 1 %upto% nprefishcols, nfish + nprefishcols + 1 + (1 %upto% nfish))
    }

    wells <- wells[ dont_drop_yet]
    is.fish.col <- is.fish.col[ dont_drop_yet]
  } else { # only one version of counts
    nfish <- sum( is.fish.col)
  }

  star.start <- substring( stuffo, 1, 1)=='*'
  droppo <- which.max( cumprod( star.start) * cumsum( star.start))
  metafish <- stuffo[ (1 %upto% droppo)] # %except% well.row] now keeping that
  stuffo <- stuffo[ -(1 %upto% droppo)]

  chardat <- strsplit( stuffo, ',', fixed=TRUE)
stopifnot( length( unique( lengths( chardat))) == 1)

  cols <- chardat[[1]]
  ncols <- length( cols)
  data <- matrix( unlist( chardat[-1]), ncol=ncols, byrow=TRUE)
  rm( chardat)

  seqinfo <- data[,!is.fish.col,drop=FALSE]
#   fishinfo <- data[,is.fish.col,drop=FALSE]
  fishinfo <- matrix( as.integer( data[,is.fish.col]), ncol=sum( is.fish.col))
  rm( data)

  deduce_dart_header_rows()

  # Loci: check numeric convertibelity by column
  # This doesn't *guarantee* numericism (eg 1..-2 is not convertible) but FFS...
  unint.char <- grepl( '[^0-9-]', seqinfo)
  unreal.char <- grepl( '[^0-9.eE+-]', seqinfo)
  dim( unint.char) <- dim( unreal.char) <- dim( seqinfo)

  intable <- colSums( unint.char)==0
  realable <- (colSums( unreal.char)==0) & !intable

  seqinfo <- data.frame( seqinfo, stringsAsFactors=FALSE)
  names( seqinfo) <- cols[ !is.fish.col]
  seqinfo[ intable] <- lapply( seqinfo[ intable], as.integer)
  seqinfo[ realable] <- lapply( seqinfo[ realable], as.numeric)
  seqinfo <- seqinfo[ !duplicated( cols[ !is.fish.col])] # ClusterIdx can appear twice (well 3 times overall);
  # ... name will be changed by names<- above
  seqinfo$count_sum <- rowSums( fishinfo)

  dupseq <- duplicated( seqinfo$Tag) # shouldn't happen but does. This format is iffy...
  # At least the entire row seems duplicated
  fishinfo <- fishinfo[ !dupseq,] # not yet transposed
  seqinfo <- seqinfo[ !dupseq,]

  # Make sure "clusters" are contiguous, with commonest "ref" sequence first
  uclu <- with( seqinfo, match( ClusterIdx, unique( ClusterIdx)))
  o <- order( uclu, -seqinfo$count_sum) # NB ascending/descending; "count_sum" in case Dart's "CountSum" is missing
  if( !all( o==seq_along( o))) {
    seqinfo <- seqinfo[ o,]
  }

  seqinfo <- within( seqinfo, {
    CloneID <- sprintf( 'L%i', ClusterIdx)
    Allele <- '' # abbreviated version, filled in later
  })
  names( seqinfo) <- sub( '^Tag$', 'FullAlleleSeq', names( seqinfo))

  tot_alleles <- nrow( seqinfo)
  n_alleles <- table( seqinfo$ClusterIdx)
  locinfo <- data.frame(
      CloneID= sprintf( 'L%s', names( n_alleles)),
      ClusterTempIndex= names( n_alleles),
      n_alleles= c( n_alleles),
      consensus='',
      var_pos='',
      row.names=NULL, stringsAsFactors=FALSE)
  locinfo$end_col <- cumsum( n_alleles)
  locinfo$start_col <- c( 1L, head( locinfo$end_col, -1)+1)


  for( iclu in seq_along( n_alleles)) {
    if( show_progress && (iclu %% 100 == 0)) {
      cat( '\r', sprintf( '%i of %i clusters', iclu, length( n_alleles)))
    }
    which_seqs <- which( seqinfo$ClusterIdx == names( n_alleles)[ iclu])
    tags <- seqinfo$FullAlleleSeq[ which_seqs]
    raws <- matrix( charToRaw( paste( tags, collapse='')), ncol=length( tags))
    is_diff <- raws != raws[,1]
    has_var <- rowSums( is_diff) > 0
    var_pos <- which( has_var)
    locinfo$var_pos[ iclu] <- paste( var_pos, collapse=',')
    conseq <- raws[,1]
    conseq[ has_var] <- charToRaw( '.')
    locinfo$consensus[ iclu] <- rawToChar( conseq)
    seqinfo$Allele[ which_seqs] <- sprintf( 'L%s_%s',
        names( n_alleles)[iclu],
        apply( raws[ var_pos,,drop=FALSE], 2, rawToChar))
  }
  if( show_progress) {
    cat( '\n')
  }

  # For compatibility with other funcs that expect ClusterTempIndex column
  names( seqinfo) <- sub( 'ClusterIdx', 'ClusterTempIndex', names( seqinfo), fixed=TRUE)
  names( locinfo) <- sub( 'ClusterIdx', 'ClusterTempIndex', names( locinfo), fixed=TRUE)

  # Construct 'NGS_count_ar' object
  dim( fishinfo) <- c( tot_alleles, nfish)
  fishinfo <- fishinfo[ o, ] # match changes in seqinfo
  fishinfo <- t( fishinfo)
  dimnames( fishinfo) <- list( Well=wells[ is.fish.col], NULL)
  class( fishinfo) <- c( 'NGS_count_ar', class( fishinfo))

  # Rename...
  locinfo <- locinfo %SUCH.THAT% (. %is.not.a% 'numeric')
  locinfo$Locus <- locinfo$CloneID
  locinfo <- locinfo[ c( 'Locus', names( locinfo) %except%
      cq( Locus, CloneID, ClusterTempIndex))]
  seqinfo$Locus <- seqinfo$CloneID
  seqinfo <- seqinfo[ c( 'Locus', names( seqinfo) %except%
      cq( Locus, CloneID, ClusterTempIndex))]

  fishinfo@info <- metafish
  fishinfo@locinfo <- locinfo
  fishinfo@seqinfo <- seqinfo

return( fishinfo)
}


#' @importFrom atease @ @<-
#' @importFrom mvbutils cq do.on FOR
"refit_whopper" <- function(
    csvfiles, # or length-2 where #1 is zipfilename, #2 is subfile
    gobj,
    max_loci_at_once=100,
    cutoff_npoly= 0.5,
    dropped_fish_file= NULL) {
## 'gobj' must have been thru 'geno6way' (or 'load_whopper') already
## To discard existing fish before loading new ones: gobj=mygobj[0,]

  define_genotypes()
  locinfo <- gobj@locinfo
  n_loci <- nrow( locinfo)

  n_fish <- do.on( csvfiles, defish_csv( ., just_count_the_fish=TRUE))

  n_read <- 0
  n_loc_left <- n_loci

  tf <- tempfile() # will be re-used for each cut operation
  on.exit( unlink( tf))

  # Read all the fish-info
  n_finfo_lines <- 0
  repeat{
    first_line <- scan( csvfiles[1], what='', sep='\n', n=1, skip=n_finfo_lines, quiet=TRUE)
    if( grepl( '^[*]?,', first_line)) { # star or blank
      n_finfo_lines <- n_finfo_lines+1
    } else
  break
  }

  mean_fish_tot <- gobj@mean_fish_tot
  rescalor <- locinfo$rescalor

  new_finfo <- FOR( csvfiles, read_cluster_dart3( ., n_lines_max=n_finfo_lines+3)@info) # at least 2 fish from each

  scatn( 'checking locus info')
  defish_csv( csvfiles[1], just_count_the_fish=FALSE, subset_fish=1:2, loci=TRUE, outfile=tf)

  mean_fish_tot <- gobj@mean_fish_tot
  rescalor <- locinfo$rescalor


  # Cannot use 'read_cluster_dart3' since it switches the order of rows without telling...
  # so would be permuted WRTO direct scan later on
  # What we really need is the "Tag" field, so we can later look up the original "FullAlleleSeq" etc
  # ... in the list of Tags in the new file

  newloci <- scan( tf, skip=n_finfo_lines, what='', sep='\n', quiet=TRUE) # keep field-name line
  fields <- strsplit( newloci[1], ',')[[1]]

  i_tag_field <- match( 'Tag', fields, 0)
  if( !i_tag_field) {
stop( 'No "Tag" field..?')
  }

  clusters <- sub( ',.*', '', newloci)
  nseqs <- table( clusters)[ as.character( unique( clusters))] # ensure matching order

  if( i_tag_field>1) {
    new_alleles <- sub( sprintf( '([^,]*,){%i}([^,]*),.*', i_tag_field-1), '\\2', newloci[-1])
  } else {
    new_alleles <- sub( ',.*', '', newloci-1)
  }

  # Now match alleles...

  mmseq1 <- match( locinfo$FullRefSeq, new_alleles, 0)
stopifnot( all( mmseq1>0))

  mmseq2 <- match( locinfo$FullAltSeq, new_alleles, 0)
stopifnot( all( mmseq2>0))

  polyallelic <- !is.na( locinfo$ThirdAllele)
  mmseq3 <- match( locinfo$FullThirdSeq[ polyallelic], new_alleles, 0)
stopifnot( all( mmseq3>0)) # some loci don't have 3rd alleles


  # Keep all fish, but discard locus header fields
  defish_csv( csvfile, just_count_the_fish=FALSE, subset_fish=seq_len( n_fish), loci=FALSE, outfile=tf)

  loci_read <- 0
  while( loci_read < n_loci) {
    n_loci_this_time <- min( max_loci_at_once, n_loci-n_loci_read)
    these_loci <- loci_read + seq_len( n_loci_this_time)

    scatn( 'Fitting loci %i:%i of %i', loci_read+1, loci_read + n_loci_this_time, n_loci)

    c3 <- array( 0, c( n_fish, n_loci_this_time, 3)) # Ref, Alt, Others

    for( ifile in csvfiles) {
      these_fish <-
      counts <- scan( ifile, skip=n_finfo_lines+1+sum( nseqs[ seq_len( loci_read)]),
          nlines= sum( nseqs[ n_loci_this_time]),
          what=0L, sep=',', quiet=TRUE) # +1 cos of field-names line
      dim( counts) <- c( length( these_fish), n_new_seqs)

      # Put these counts into a temp F*L*3 array, like 'pick_ref_alt'

      c3[ these_fish,,1] <- counts[ ,mmseq1[ these_loci]]
      c3[ these_fish,,2] <- counts[ ,mmseq2[ these_loci]]
      c3[ these_fish, polyallelic[ these_loci], 3] <- counts[ ,mmseq3[ these_loci] ]
    }

    # Renorm by fish...
    c3[f,l,a] := c3[f,l,a] * (mean_fish_tot / new_finfo$Fishtot[ f])

    # Rescale alleles, as in ppn_ref_alt_check2
    c3[f,l,a]:= c3[f,l,a] * rescalor[l,a]


    n_loci_read <- n_loci_read + n_loci_this_time
  }

  # Just do it...
  new_locinfo <- choose_geno6_thresholds( ...)

  # Tack on attributes...
  attributes( geno) <- c( attributes( geno), attributes( gobj)[ cq( locinfo, diplos, mean_fish_tot, het_cut, calls, args)])
  # May be an old "npoly" field in fish-info--- drop, since it won't exist for new data
  # Drop any other info fields that aren't in both old and new
  old_finfo <- gobj@info[ names( gobj@info) %that.are.in% names( new_finfo)]
  new_finfo <- new_finfo[ names( new_finfo) %that.are.in% names( gobj@info)]
  geno@info <- rbind( old_finfo, new_finfo)
  geno@locinfo$npoly <- NULL # also meaningless
  oldClass( geno) <- 'snpgeno'

  toopoly <- poly_fish > cutoff_npoly
  if( length( as.character( dropped_fish_file)==1)) {
    dropped_fish <- geno@info$Our_sample[ toopoly] # really should have something to ID the *tube* kinda...
    cat( sprintf( '%s: %i', dropped_fish, poly_fish[ toopoly]), file=dropped_fish_file, sep='\n')
  }

  geno <- geno[ !toopoly,]
return( geno)
}

#' @importFrom handy2 integ
#' @importFrom gbasics ridder
"renorm_SPA" <- function(K, dK, ddK, return_what=c( 'func', 'mulfuncby'),
                         tol=formals( ridder)$tol
  # , ... ; should really allow extra args to K & co, and build them in...
){
  isqrt_2pi <- 1/sqrt( 2*pi)
  absmax <- 10 / sqrt( ddK( 0)) # x between +/- 10SD of mean

  # mc <- as.list( match.call( expand.dots=FALSE)$...)
  K <- Vectorize( K) # names( formals( K)) %except% names( list( ...)) kinda thing
  dK <- Vectorize( dK)
  ddK <- Vectorize( ddK)

  itotto <- 1
  sfunc <- function( s, ddK_s=ddK(s)) {
      x <- dK( s)
      itotto * exp( K( s) - s * x) * sqrt( ddK_s) * isqrt_2pi
    }

  itotto <- 1 / integ( sfunc( x), -absmax, absmax)

  return_what <- match.arg( return_what)
  if( return_what=='mulfuncby')
return( itotto)

  # Otherwise we need the x-ready version, with (vectorized) root-finding
  xfunc <- function( x) {
    # Lower & upper bounds for s

    # Aim for exact, miss, then try to get equal dist the other side
    iddK0 <- 1/ddK( 0)
    s1 <- (x-dK(0)) * iddK0

    dK1 <- dK( s1)
    # Paranoia: might be perfect!
    bingo <- dK1==x
    if( any( bingo)) { # otherwise, if K()==Vectorize(...), it fucks up and turns it all into a list FFS
      s1[ bingo] <- 1.99 * s1[ bingo]
      dK1[ bingo] <- dK( s1[ bingo])
    }

    ddK1 <- ddK( s1)
    s2 <- s1 - 2*(dK1-x) / ddK1
    dK2 <- dK( s2)
    while( any( same_sign <- (dK1-x)*(dK2-x)>1) ) { # unlikely; try a bit further
      s2 <- s2 - same_sign * (dK1-x)/ddK1
      dK2[ same_sign] <- dK( s2[ same_sign])
    }

    # 'ridder' wants a vector func with no args
    dK_min_x <- function( s) dK(s)-x
    s <- ridder( dK_min_x, pmin( s1, s2), pmax( s1, s2), tol=tol) # root finder
    ddK_s <- ddK( s)

    # Undo the reordering
  return( sfunc( s, ddK_s) / ddK_s)
  }

return( xfunc)
}

#' @importFrom mvbutils returnList
"renorm_SPA_cumul" <- function( K, dK, ddK, sd_half_range=10, n_pts=2001) {
  x <- 0
  SPA_s_dxds <- function( s) {
    x <<- dK( s)
    # SPA is (K(x) - s*x) / sqrt( 2*pi*ddK( s))
    # but dx/ds = ddK( s)
    exp( K(s) - s*x) * sqrt( ddK( s) / (2*pi))
  }

  sd <- sqrt( ddK( 0))
  # To norm, we'd integ X over say +/- 10 sd
  # but s ~= (x-mu) / (sd)^2 hence s-range below

  spoints <- seq( -sd_half_range, sd_half_range, length=n_pts) / sd
  pdf_s <- SPA_s_dxds( spoints) * diff( spoints[1:2]) # diff() gets the integral about right
  C <- 1/sum( pdf_s) # should be close to 1
  cdf <- cumsum( pdf_s) * C # nicely renormalized

  inv_CDF_bod <- splinefun( cdf, x,  method='hyman')
  inv_CDF <- function( p) inv_CDF_bod( p) # sensibler arg name, and no deriv arg
  CDF <- splinefun( x, cdf, method='hyman')
returnList( CDF, inv_CDF)
}

#' @export
"set_recording" <-
function( vars, record=TRUE) {
  if( record) {
    selfy <- new.env( parent=environment( recordar))
    selfy$subs <- structure( vector( 'list', length=length( vars)), names=vars)
    selfy$exprs <- quote( {})
    environment( recordar) <- selfy
  } else {
    recordar <- function( expr, expand_dim)  eval.parent( substitute( expr)) # nothing will be recorded!
  }
  assign( '?', recordar, parent.frame())
return( invisible( recordar))
}

#' @importFrom mvbutils cq mlocal
"set_thresholds" <-
function( keeping, nlocal=sys.parent()) mlocal({
stopifnot( keeping %in% cq( hi, lo))

  symmo <- my.all.equal( subset1, subset2)
  if( is.null( eta) || is.null( keep_thresh)) {
    probinverts <- numeric()
    if( is.null( eta)) {
      probinverts <- 1/one_in_X_eta
    }
    if( is.null( keep_thresh)) {
      probinverts <- c( probinverts,
          rough_n_pairs_to_keep / (length( subset1) * length( subset2) / (1+symmo)))
    }

    if( keeping == 'hi') {
      probinverts <- 1-probinverts
    }

    XX <- try( inv_CDF_SPA2( probinverts, K, dK, ddK))
    if( XX %is.a% 'try-error') {
  warning( "Couldn't set thresholds via Lugannini-Rice SPA (will use alternative); ' %&%
      'probably too extreme for this distro too.")
      # Use renormed sum-of-pdf:
      XX <- inv_CDF( probinverts)
    }

    if( is.null( eta)) {
      eta <- XX[1]
      XX <- XX[-1]
    }
    if( is.null( keep_thresh)) {
      keep_thresh <- ( if( keeping == 'lo') max else min)( XX[1], eta)
    }
  }
})


#' @importFrom mvbutils mlocal FOR
"setup_counts_and_nhi" <- function( nlocal=sys.parent()) mlocal({
#### Called either at start of locus-fit, or before refit
  # First two lines may be done already, but cheap to redo
  is_OO <- genoi==OO
  is_hetz <- genoi %in% c( AB, BC, AC)

  is_moderate <- !is_OO & (dat < hix)
  counts <- FOR( c( AAO, BBO, CCO), dat[ is_moderate & geno[,i]==.])
  counts$Hetz <- dat[ is_moderate & is_hetz]

  is_hi <- dat >= hix
  n_hi <- FOR( c( AAO, BBO, CCO), sum( geno[,i]==. & is_hi))
  n_hi$Hetz <- sum( is_hetz & is_hi)
})


#' @importFrom atease @ @<-
#' @importFrom mvbutils mlocal mvb.match.call %not.in%
#' @importFrom debug mtrace
"setup_recordo_check_baits" <- function( nlocal=sys.parent()) mlocal({
### Declutter parent function
  nfish <- nloci <- -1
  force_redo <- skip_step <- FALSE
  if( !is.null( pbtidy) && any( names( pbtidy@calls)==stop_after)) {
    scatn( "Re-running: request is to do FEWER steps than already completed--- so can't trust previous results")
    force_redo <- TRUE
  }

  mc <- mvb.match.call() # to check args

  if( mtrace.) {
    record <- record # move it here
    mtrace( record)
  }

  recordo <- function( name, expr) {
      if( skip_step) { # iff 'stop_after' is in force
        scatn( 'Skipping step "%s"', name)
        pbtidy@calls[[ name]] <- NULL # erase all calls after 'stop_after'; see force_redo setting earlier
        pbtidy@args[[ name]] <- NULL
      } else {
        name_args <- mc[[ sprintf( 'args_%s', name)]] # default NULL
        # Did have: get0( sprintf( 'args_%s', name)) but triggers error in new incremental version because of stop() in default parameters

        # Don't redo (unless forced) if no name_arg specifically set--- assume happy with existing args
        doing <- force_redo || (name %not.in% names( pbtidy@calls))
        if( !doing) { # staged checks to allow message
          doing <- !is.null( name_args) && !my.all.equal( name_args, pbtidy@args[[ name]])
          if( doing) {
            scatn( 'Arguments have changed: re-doing "%s"', name)
          }
        }

        if( doing) {
          mc <- match.call()
          mc$var <- 'pbtidy'
          mc$expr <- do.call( 'substitute', list( mc$expr, name_args)) # so actual numbers appear
          mc$these_args <- name_args
          mc[[1]] <- quote( record)
          eval.parent( mc)
          force_redo <<- TRUE # all subsequent steps (until stop_after, if any)
        }

        new_nloci <- length( unique( pbtidy@locinfo$Locus))
        new_nfish <- nrow( pbtidy)
        if( doing) { # show progress
          if( nloci > 0) {
            scatn( '%s: elim %i fish, %i loci',
                name, nfish-new_nfish, nloci-new_nloci)
          } else { # reloading
            scatn( '%s: starting from %i fish, %i loci',
                name, new_nfish, new_nloci)
          }
        } # if doing
        nfish <<- new_nfish
        nloci <<- new_nloci
      } # if not skipping

      if( name == stop_after) {
        skip_step <<- TRUE # next
      }

    return( NULL)
    } # recordo

})


#' @importFrom atease @ @<-
#' @importFrom handy2 rsample
#' @importFrom stats var
"simtest_Kstuff" <- function( ck, n, nq=20) {
  # ck needs locinfo$LOD
  extract.named( ck@Kenv) # 4ways pretending to be 6ways
  mg <- LOD@mg
  n_loci <- nrow( PUP)

  # Direct simulation of genotypes is tricky because of errors
  # ... though for completeness SHOULD really try that here
  # Instead, use PUP table

  g12_code <- matrix( 0L, n, n_loci)
  for( il in 1:n_loci) {
    this_g12 <- rsample( n, colnames( PUP), prob=PUP[il,], replace=TRUE)
    g1 <- substring( this_g12, 1, 2)
    g2 <- substring( this_g12, 4, 5)
    g12_code[ ,il] <- mg[ cbind( g1, g2)]
  }

  # NYI in vecless:
  # LOD_obs[ l, i] := LOD[ l, g12[ l, i]]
  # PLOD[ i] := LOD_obs[ +., i]

  LOD_obs <- matrix( 0, n, n_loci)
  LOD_obs[] <- LOD[ cbind( rep( 1:n_loci, each=n), c( g12_code)) ]
  PLOD <- rowSums( LOD_obs)

  scatn( 'Emp mean %6.4f SPA %6.4f', mean( PLOD), dK( 0))
  scatn( 'Emp var %6.4f SPA %6.4f', var( PLOD), ddK( 0))

  # %les:
  qq <- (2:nq-1)/nq
  pciles <- inv_CDF( qq)

  bin <- 1+findInterval( PLOD, pciles)
  counts <- tabulate( bin, nbins=nq)
  scatn( 'Expect %5.1f PLODs in each of %i percentile-bins; got this:', n/nq, nq)
  print( counts)

return( invisible( PLOD))
}


"split_geno_4to6" <-
function( geno, cA, cB, info, li, best_cut=li$best_cut) {
  define_genotypes()
  # best_cut should be a loci*2 matrix
  # geno, cA, cB will be sample*loci

  if( is.null( dim( best_cut))) { # then just 1 locus
    best_cut <- matrix( best_cut, ncol=length( best_cut), dimnames=list( NULL, names( best_cut)))
  }

  n_fish <- nrow( geno)
  n_loci <- ncol( geno)
  if( is.null( n_fish)) {
    n_fish <- length( geno) # not a matrix...
    n_loci <- 1
  }

  g6 <- snpgeno( n_fish, n_loci, diplos=genotypes6, info=info, locinfo=li)
  diplos <- as.raw( seq_along( genotypes6))
  names( diplos) <- genotypes6

  g6[ geno==OO] <- diplos[ OO]
  g6[ geno==AB] <- diplos[ AB]

  g6[ geno==AC] <- diplos[ AO]
  g6[ geno==BC] <- diplos[ BO]
  g6[ geno==CCO] <- diplos[ OO]

  g6[ geno==AAO] <- diplos[ AA] # default: if has_C, or big-enough A-count
  g6[ geno==AAO & (cA < rep( best_cut[,'A'], each=n_fish) ) ] <- diplos[ AO] # ... the non-default

  g6[ geno==BBO] <- diplos[ BB]
  g6[ geno==BBO & (cB < rep( best_cut[,'B'], each=n_fish) ) ] <- diplos[ BO]
return( g6)
}


#' @importFrom stats splinefun
"splug_transform" <- function( delta) {
  # Define a controlled monotone mapping from [0,1] to [0,1]: 0 -> 0, 1 -> 1
  # Mapping will be determined by delta, a vec of vals in [-inf,inf]
  # delta=rep(0,n) should return linear map
  # Thing returned is a map function that can be applied to arby
  #  input vec (all elts in [0,1])

  n <- length( delta)
  stopifnot( n>0)

  p <- rep( 0, n+2)
  for( i in 1:n) {
    p[ i+1] <- p[ i] + (1-p[ i]) * inv.logit( logit( 1 / (n+2-i)) + delta[ i])
  }
  p[ n+2] <- 1

  # Sneakily, return the INVERSE, which allows more drastic tmfns
  #return( approxfun( p, seq( 0, 1, length=n+2)))
  # Smooth version is nice
  return( splinefun( p, seq( 0, 1, length=n+2), method='hyman'))
}


"widio" <-
function( x, m) {
### Find widest interval in x that contains at least m datapoints and split

  xx <- sort( x)
  gapi <- which.max( diff( xx, m))
  xx <- xx[ gapi+1:m-1]

  # Original version started with quantiles
  #   n <- length( x)
  #   qq <- quantile( x, seq( 0, 1, length=2*n/m))
  #   gapi <- which.max( diff( qq, 2)) # ie from qq[ gapi] to qq[ gapi+2]
  # Now just a few points; sort them and look for biggest gap(s)
  #   xx <- x %such.that% (. %in.range% qq[ gapi + c( 0, 2)])
  #  xx <- sort( xx)

  gapi <- mean( which.max( diff( xx)))
  rangio <- floor( gapi):(1 + ceiling( gapi))
  res <- mean( xx[ rangio])
return( res)
}

