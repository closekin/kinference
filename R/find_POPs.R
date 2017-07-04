#' Kin-finders for loads-of-SNPs datasets
#'
#' @aliases find_POPs find_HSPs find_duplicates
#'
#' @description
#' These take a \code{snpgeno} dataset that has been processed as far as \code{check6and4} (and for HSPs, \code{prepare_PLOD_SPA}) and find various relations between the samples. Relationships include duplicates (DUPs/dupes; \code{find_duplicates}), parent-offspring pairs (POPs; \code{find_POPs_vs}) and half-sibling pairs (HSPs; \code{find_HSPs}). Unrelated pairs are referred to as UPs. One can specify the same or different subsets of the \code{snpgeno} for comparison: e.g., first subset for the adults, second for the juveniles.
#'
#' @section Kinformation:
#' The idea is that kin-finding is based on a statistic and a threshold \code{eta}, where the latter is chosen to keep false-positives down to a user-specified level. Anything "beyond" \code{eta} will be treated as a kin-pair ("beyond" depends on how the statistic is defined, i.e. whether a kin-pair should come out very low or very high). However, you're also likely to want to look post hoc at the distro of computed statistics \emph{near} \code{eta}, to see whether separation is as clean (or otherwise) as expected--- and also very unbeyond \code{eta} into the zone where UPs are entirely dominant, to check that theory is OK. So, as well as returning the "interesting" pairs that have a statistic close to or on the non-UP size of \code{eta}, the POP and HSP versions also return \emph{summaries} of the distribution of the statistic. The thing is that there will be zillions of statistics from UPs--- enough to blow out computer memory--- and they are not individually interesting. Specifically, the main things returned are:
#'
#' \itemize{
#' \item mean and variance of stats. Computation is restricted to those on the UP-side of \code{eta} (which is nearly all of them, usually) in order to avoid distortion from non-UP cases. The latter will often be so rare that distortion would be negligible--- but means and variances are not "robust", so . Almost all will be include
#' \item counts of binned stats, regardless of whether above or below \code{eta}. The bins are set based on SPAs to the theoretical distributions, and chosen so that an equal number of UP-pairs should fall into each bin.
#' \item cases where the stat is "interesting", i.e. on the non-UP side of \code{keep_thresh}, as a \code{data.frame}. See \bold{Value} for details
#' }
#'
#' The process is controlled by three numbers: \code{nq} for number of bins, \code{eta} itself, and some nearby threshold \code{keep_thresh} on the UP-side of \code{eta} (it will be automatically set to \code{eta} otherwise) to determine which pairs are explicitly retained for your inspection. There are two ways to specify \code{eta} and \code{keep_thresh}. Usually, you would start with the indirect method, where you choose the predicted-false-positive proportion of UP-pairs via the parameter \code{one_in_X_eta}, and \code{rough_n_pairs_to_keep}. The routines then use SPAs to the corresponding values of \code{eta} and \code{keep_thresh}; the returned value of \code{eta} is what you can subsequently use to make the actual kin-decisions yourself after the event (by subsetting the "interesting" pairs, comparing the statistic for each pair to \code{eta})--- assuming that observed does match expected.
#'
#' But, sometimes it doesn't. In that case, the predicted values of \code{eta} and \code{keep_thresh} may be way off the mark, and lead to retaining faaar too few or too many pairs. If so, then look at the histogram of retained statistics from an initial run, and try setting \code{eta} and/or \code{keep_thresh} manually, rather than futzing around with the indirect parameters until you get what you were after.
#'
#' @section Duplicates:
#' You have to set the retention threshold manually, via \code{max_diff_genos} (see arguments). Post-processing step needed to to get the indices to remove -- use \code{\link{drop_dups_pairwise_equiv}}, see \bold{Examples}.
#'
#' To avoid running \code{find_duplicates} on large numbers of fish at once, one can split the dataset; see \bold{Examples}. You first need to run on each subset separately (avoiding a quadratic number of comparisons) and reduce it to non-duplicates (again, see \code{\link{drop_dups_pairwise_equiv}}), then check the pair of reduced subsets (this will compare everything in the first to everything in second, as the subsets are different). Note that when the subsets are different, comparisons are made only \emph{between} subsets, not \emph{within} each subset.
#'
#' Uses 4-way genotyping only, since these should be largely error-free. (Looks like the exceptions are from samples with dodgy DNA.)
#'
#' @section Parent-offspring pairs:
#' 4-way genotypes are used to find "pseudo-exclusions" of the form AAO/BBO, which \emph{usually} means AA/BB or AO/BB or AA/BO (a true exclusion), but \emph{could} mean AO/BO (not an exclusion).
#'
# \code{find_POPs} merely counts these, for loci where \code{Pr[O]+Pr[C]<pOC_max} in order to avoid excess noise from AO/BO cases.
#' \code{find_POPs} uses all loci (by default) but weights them semi-optimally so that pseudo-exclusions from loci with high "false pseudo-exclusion probability" (i.e., high \code{Pr[AO/BO|UP]}) count much less than ones from loci with very low null rates, for which AAO/BBO almost certainly means AA/BB. We call this "Weighted PSeudo-EXclusion" ("WPSEX").
#'
#' The case AB/OO is also a (non-pseudo) exclusion, but is rarer than AAO/BBO (non-existent for loci without nulls, of course). The count of such cases is included in the output for "interesting" pairs in \code{find_POPs}; see \bold{Value}.
#'
#' POP-finding is based on 4-way genotypes (OO, AAO, BBO, AB) to avoid complications from genotyping-error-rates, and uses pseudo-exclusions rather than likelihood-ratios; the latter is very sensitive to false-negative-exclusions arising from typing-error, or even from mutation with so many SNPs. You can get round that by including estimates of typing-error-rate, but that's not necessarily easy to estimate in advance insofar as it applies per-locus to POPs.
#'
#' @section Speed:
#' These are written in C (\code{Rcpp}) for speed, but for big datasets they might still be quite slow. In the first instance, I certainly wouldn't try them on 20,000 fish at once; I'd try with say 1000 then if that's OK 5000 etc. Bear in mind that they can always be run on different subsets of the data, and the results patched back together (results will not change by doing that). If you can run jobs in parallel, that could help a lot.
#'
#' @param snpg a \code{snpgeno} object
#' @param subset1,subset2 numeric vectors of which samples to use (not logical, not negative). Defaults to all of them. Iff the two subsets are identical, only half the comparisons are done (i.e., not i with j then j with i). Some sanity checks are done.
#' @param alpha (\code{find_POPs}) Loci receive a weight which is proportional to (difference in probability of pseudo-exclusion between UP and POP) / (variance of indicator of pseudo-exclusion). But, should this be variance assuming UP or POP? \code{alpha} sets the balance; bigger values make it more UPpity, so placing more emphasis on avoiding false-positives (which is probably the Right Thing To Do). 0.999 could be completely fine... (but hopefully \code{alpha} won't affect the result much anyway.)
#' @param one_in_X_eta expected number of false-positive UPs you can tolerate. Setting this to say \code{1e6} means you'd expect 1 per million comparisons. Used to set the threshold \code{eta}, which is returned automatically.
#' @param rough_n_pairs_to_keep For checking, you can set this to trap many more high-scoring pairs than you expect there to "really" be, say a few thousand (NB the number of pairs retained won't exactly equal this). You can subsequently look at the "lucky losers" with high but sub-'eta' stats, and then filter them out yourself by applying a cutoff of \code{eta}. If you leave \code{rough_n_pairs_to_keep} at its default of NA, the trap will be set at \code{eta}, so \code{bigs} will contain exactly the pairs you want. Values above \code{eta} will always be kept, even if you specify something silly for \code{rough_n_pairs_to_keep}.
#' @param eta,keep_thresh see \bold{Description}. Can specify either or both. These override \code{one_in_X_eta} and \code{rough_n_pairs_to_keep} respectively.
#' @param nq number of bins to group the stats from the sub-'eta' pairs into. The bins will be set at quantiles of the expected distribution for UPs.
#' @param max_diff_genos (\code{find_duplicates}) max number of discrepant 4-way genotypes to tolerate in "identical" fish. Try increasing this from say 10 upwards, and hopefully nothing much will change (though at some point things will change a lot, as you get into the non-duplicate bit of the distribution). See \bold{Duplicates} for how to remove duplicates from the data.
#' @param quick whether to "compile" the functions for SPA, which use the magic \code{:=} operator. It speeds up the SPA bit but almost all the time is spent on actual POP-finding...
#' @param bins binning for PLODs (we throw away ones outside the range and bin them according to this within)
#'
#' @return A list, whose most important element is a \code{data.frame} called \code{bigs} with 3 columns: statistic (PLOD or number-of-excluding loci or \code{similar} which is number of mismatching genotypes--- though "bigs" is misleading in the duplicates case, since mismatches need to be \bold{small} to qualify), \code{i} (index in \code{subset1} of the first pair-member), \code{j} (index in \code{subset2} of the second). Note that \code{i} and \code{j} refer to the \emph{subsets}, not to the rows of the original \code{snpg}. Note that, iff you have set \code{rough_n_pairs_to_keep}, these will include pairs below the FP cutoff (which is returned as \code{eta}).
#'
#' \code{find_POPs} adds \code{bigs$nABOO}, showing the number of AB/OO exclusions for that potential POP. This is a useful additional diagnostic; it should be close to 0 for true POPs (it can only result from genotyping error or mutation, whereas AAO/BBO can result from nulls). For UPs, I was seeing values typically in the low 20s, which is pretty good separation.
#'
#' For duplicates, \code{bigs} does not record \emph{all} pairwise duplicates, unless the subsets are different--- otherwise you could have quadratic horror of enormous numbers of pairs arising from a cluster of say 100 identical controls! Since "duplication" is transitive (ie if i & j are the same, and i & k are the same, then j & k must also be the same), only the necessary ones are recorded to allow you to filter out yourself afterwards. e.g., if samples 1, 3, 5, and 6 are all duplicates, you'll get this:
#' %#
#' \item{# $bigs without "ndiff" column}{}
#' \item{  i j}{}
#' \item{  3 1}{}
#' \item{  4 3}{}
#' \item{  6 4}{}
#' but you won't see the pairings for 1/4, 1/6, 3/6. If you just want to strip out all duplicates bar one in each group (and you don't care which one is kept), then you can use the function \code{\link{drop_dups_pairwise_equiv}} --- see \bold{Examples}.
#'
#' For POPs and HSPs, the following are also returned in the list. The main point is that the "boring" below-threshold ones get put into bins, not kept individually. The names sometimes change depending on which statistic is being used.
#'
#' \item{eta}{false-positive cutoff to be applied to the statistic in \code{bigs} (automatically done if \code{rough_n_pairs_to_keep==NA}, or up to you if not). Variance of the stat will only be calculated from values to the "UP side" of \code{eta}. However, the set of retained pairs/individuals is actually controlled by...}
#' \item{keep_thresh}{the cutoff used to retain "interesting" pairs. Usually obvious from the range of stat-values in \code{bigs}.}
#' \item{mean_sub_<stat>, var_sub_<stat>}{empirical values for the statistic when it is below \code{eta} (ie nearly always).}
#' \item{mean_theory, var_theory}{of the statistic, to compare to previous.}
#' \item{n_<stat>_in_bin}{number of pairs whose stat fell within the range of each bin}
#' \item{bins}{cutpoints for the bins. These should be quantiles, according to the SPA; so if practice matches theory, the numbers-per-bin should all be similar.}
#'
#' @examples
#' ### Don't run
#' ## duplicate checking. ckmini2 has 6 fish where 1,3,4,6 are all identical (zero differing loci).
#' ## there's only 7 crappy loci and I faked the data for this anyway, so strict identical is needed
#' ## All-in-one
#' #test <- find_duplicates( ckmini2, max=0) # strict identity
#' #test$bigs
#' ##  ndiff i j
#' ##1     0 3 1
#' ##2     0 4 3
#' ##3     0 6 4
#' ## To remove them--- subtlety of keeping ONE from each group
#' #droppies <- drop_dups_pairwise_equiv( test$bigs[,2:3])
#' #droppies # 1, 4, 6
#' #ckmini2_nodups <- ckmini2[ -droppies, ]
#' ## Two-stage
#' #first_half <- 1:3
#' #second_half <- (1:nrow( ckmini2)) \%except\% first_half
#' #test1 <- find_duplicates( ckmini2, subset1=first_half, subset2=first_half, max=0)
#' #test1$bigs
#' ##  ndiff i j
#' ##1     0 3 1
#' #droppies1 <- first_half[ drop_dups_pairwise_equiv( test1$bigs[,2:3])] # NB must do lookup in subset
#' #test2 <- find_duplicates( ckmini2, subset1=second_half, subset2=second_half, max=0)
#' #droppies2 <- second_half[ drop_dups_pairwise_equiv( test2$bigs[,2:3])] # 4
#' ## Now check 2nd half vs 1st
#' #test2_1 <- find_duplicates( ckmini2,
#' #    subset1=first_half \%except\% droppies1,
#' #    subset2=second_half \%except\% droppies2,
#' #    max=0)
#' ## Simpler since no internal checks. Just remove 2nd-halfers that match something in the 1st-half
#' #droppies2_1 <- (second_half \%except\% droppies2)[ test2_1$bigs[,'j']) # 6
#' #droppies <- c( droppies1, droppies2, droppies2_1)
#' #ckmini2_nodups2 <- ckmini2[ -droppies,]
#' ## HSPs: comparing everything with itself (not sensible for real data, should take out adults first)
#' ## set threshold for 1 FP
#' #test <- find_HSPs( ckdata, one_in_X_eta=sqr( nrow( ckdata))/2 )
#' ## POPs: Ad-Ju comps; again 1 FP
#' #test <- find_POPs( ckdata, subset1=adults, subset2=juves, alpha=0.99,
#' #    one_in_X_eta=length( adults) * length( juves), rough_n_pairs_to_keep=500)
#' ### End don't run
#' @export
#' @importFrom gbasics sqr
#' @importFrom atease @
#' @importFrom vecless := compile_vecless
#' @importFrom stats runif
#' @importFrom mvbutils cq %upto% %that.are.in% my.all.equal extract.named %without.name%
"find_POPs" <-
function( snpg, subset1=1 %upto% nrow( snpg), subset2=subset1,
    alpha,
    one_in_X_eta,
    rough_n_pairs_to_keep= NA,
    eta= NULL,
    keep_thresh= NULL,
    nq,
    quick=TRUE) {
###################
  # Sanity...
stopifnot( is.numeric( subset1) && is.numeric( subset2))
stopifnot( all( !duplicated( subset1)) && all( !duplicated( subset2)))
stopifnot( my.all.equal( subset1, subset2) || !length( intersect( subset1, subset2)))

  # Decide based #apparent exclusions of AA/BB form, using 4way genos, though
  # ... it's really AAO/BBO so not a true exclu but
  # ... close among pop-loci
  # Sticking with 4way genos so that genotyping errors are low

  # Instead of Nexclu, uses a wted sum of "exclus" in 4way genos to max expected diff between POP and UP
  # This version doesn't allow for geno errors, but does realize that AO/BO could happen
  # Doesn't bother with OO-AB

  define_genotypes()
  extract.named( snpg@locinfo[ cq( use6, PUP4, pbonzer)])
  p0 <- pbonzer[,'O'] + pbonzer[,'C']
  pA <- pbonzer[,'A']
  pB <- pbonzer[,'B']

  # "Exclusion" whenever AAO & BBO, but this *could* be AO/BO
  # I'm calling them "ex" for now anyway, hence pex_BLAH
  pex <- cbind(
      POP= 2 * p0 * pA * pB,
      UP= 2 * (2*pA*p0 + pA*pA) * (2*pB*p0 + pB*pB)
    )

  # Want to keep the POP and UP means as far apart as possible on the scale of SDs...
  # ... but, which SD? Make it alpha * SD[UP] + (1-alpha) * SD[POP]

  # Optimal wt would depend on p0 and to some extent on pA
  # wt should be 1 if p0==0 and 0 if p0==1
  delta <- pex[,'UP'] - pex[,'POP'] # mathematically I think this *can't* be -ve
  SD <- sqrt( pex * (1-pex))
  SD_combo <- alpha * SD[,'UP'] + (1-alpha) * SD[,'POP'] # %*% c( alpha, 1-alpha)
  V_combo <- sqr( SD_combo)
  ww <- delta / V_combo # considerable algebra appears to show this is optimal
  ww <- c( ww / sum( ww) ) # else get 1-col matrix
stopifnot( all( ww>0))

  pop_loci <- which( ww > 0) # all of them, for now
  # pop_loci <- which( pbonzer[,'O'] + pbonzer[,'C'] < pOC_max)

  temp_snpg <- snpg[ , pop_loci]
  recode4to6temp <- function( x) { x[ x=='AO'] <- AA; x[ x=='BO'] <- BB; x}
  temp_snpg <- recode4to6temp( temp_snpg) # (AA,AO) -> AA; (BB,BO) -> BB

  # Remove extranea
  attributes( temp_snpg) <- attributes( temp_snpg)[ 'dim']
  temp_snpg <- t( temp_snpg)

  # Distro of #excl loci for UPs
  opphetz <- c( 'AAO/BBO', 'BBO/AAO') %that.are.in% colnames( PUP4)
  pex_up <- PUP4[ pop_loci, opphetz]
  rr <- pex_up / (1-pex_up)
  log1m_pexup <- log1p( -pex_up)

  K <- function( tt) {
      # vecless2 notation could 1-step this via SUM
      KK[ l, j]:= log1m_pexup[ l] + log1p( rr[ l] * exp( tt[ j] * ww[ l]))
      K[ j] := KK[ +., j]
    return( c( K)) # without the c(), you get a scalar xtensor, and trouble...
    }

  dK <- function( tt) {
      retw[l,j] := rr[ l] * exp( tt[ j] * ww[ l])
      dKK[ l, j] := ww[ l] * retw[ l, j] / (1+retw[ l, j]) # guessing this is more accurate tahn 1-1/1+x
      dK[ j] := dKK[ +., j]
    return( c( dK))
    }

  ddK <- function( tt) {
     retw[l,j] := rr[ l] * exp( tt[ j] * ww[ l])
      ddKK[ l, j] := sqr( ww[ l]) * retw[ l, j] / sqr( 1+retw[ l, j]) # guessing this is more accurate tahn 1-1/1+x
      ddK[ j] := ddKK[ +., j]
   return( c( ddK))
   }

  n_loci <- ncol( snpg)
  n_sim_check <- 1000
  Ktest <- function( tt) {
    x <- matrix( runif( n_sim_check * n_loci) < pex_up, n_loci, n_sim_check)
    ewx <- exp( x*ww*tt)
    colSums( log( ewx))
  }

  if( quick) {
    K <- compile_vecless( K(0))
    dK <- compile_vecless( dK(0))
    ddK <- compile_vecless( ddK(0))
  }

  symmo <- my.all.equal( subset1, subset2)
  inv_CDF <- renorm_SPA_cumul( K, dK, ddK)$inv_CDF

  set_thresholds( keeping='lo') # cf 'hi' for HSPs

  # Prepare for diagnostics of #excl
  # prolly not needed but HSP C code already does it
  qq <- (2:nq-1)/nq
  pciles <- inv_CDF( qq) # inv_CDF_SPA2 may struggle with LOWER tail...

  # Trying special-cases here to minimize copying
  if( symmo){
    if( !my.all.equal( subset1, 1 %upto% ncol( temp_snpg))) {
      temp_snpg <- temp_snpg[, subset1]
    }
    result <- POP_wt_paircomps_lots(
        geno1= temp_snpg,
        geno2= temp_snpg,
        w= ww,
        symmo= TRUE,
        eta= eta,
        max_keep_wpsex= keep_thresh,
        bins= pciles,
        AAO= match( 'AA', snpg@diplos), # NB NB: AO has been recoded to AA
        BBO= match( 'BB', snpg@diplos)
      )
  } else { # different subsets
    result <- POP_wt_paircomps_lots(
        geno1= temp_snpg[ ,subset1],
        geno2= temp_snpg[ ,subset2],
        w= ww,
        symmo= FALSE,
        eta= eta,
        max_keep_wpsex= keep_thresh,
        bins= pciles,
        AAO= match( 'AA', snpg@diplos), # NB NB: AO has been recoded to AA
        BBO= match( 'BB', snpg@diplos)
      )
  }

  # bigs is a misnomer for POPs; smalls is more like it (since we want few exclusions)
  result$bigs <- with( result, data.frame( wpsex=big_wpsex, i=big_i, j=big_j))
  # Check nABOO, only for interesting pairs
  snpg_i <- snpg[ subset1[ result$bigs$i], pop_loci]
  snpg_j <- snpg[ subset2[ result$bigs$j], pop_loci]
  isABOO <- ((snpg_i==OO) & (snpg_j==AB)) + ((snpg_i==AB) & (snpg_j==OO))
  result$bigs$nABOO <- rowSums( isABOO)

  result <- result %without.name% cq( big_wpsex, big_i, big_j)
  result <- within( result, {
    bins <- pciles
    eta <- eta
    keep_thresh <- keep_thresh
    n_loci <- length( pop_loci)
    mean_theory <- dK( 0)
    var_theory <- ddK( 0)
  })
  result$call <- sys.call() # iffy within within

return( result)
}

