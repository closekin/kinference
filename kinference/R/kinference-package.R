#' Info about this package
#'
#' Package \pkg{kinference} contains pairwise kin-finding and various QC
#' routines, aimed ultimately for close-kin mark-recapture. This document ought
#' to explain the different categories of Stuff within the package (eg QC,
#' kin-finding, kin-splitting, ...). Perhaps there should be separate
#' help-files for each category, that contain links to the actual function
#' helps.
#'
#' The most important doco is the vignette (but it's broken at time of
#' writing!) which shows the operations in more detail.
#'
#' @aliases kinference kinference-package
#' @author Mark V Bravington, David L Miller, Shane M Baylis
#' @keywords misc
#' @useDynLib kinference, .registration=TRUE
#' @import Rcpp
#' @import atease
#' @import mvbutils
#' @import gbasics
#' @import vecless
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @importFrom atease "@"
#' @importFrom atease "@<-"
#' @importFrom gbasics inv.logit
#' @importFrom gbasics logit
#' @importFrom gbasics make_genopairer
#' @importFrom gbasics rsample
#' @importFrom gbasics snpgeno
#' @importFrom gbasics sqr
#' @importFrom grDevices rgb
#' @importFrom graphics legend
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom mvbutils "%&%"
#' @importFrom mvbutils "%except%"
#' @importFrom mvbutils "%is.a%"
#' @importFrom mvbutils "%is.not.a%"
#' @importFrom mvbutils "%not.in%"
#' @importFrom mvbutils "%that.are.in%"
#' @importFrom mvbutils "%upto%"
#' @importFrom mvbutils "%where%"
#' @importFrom mvbutils "%without.name%"
#' @importFrom mvbutils "?"
#' @importFrom mvbutils FOR
#' @importFrom mvbutils cq
#' @importFrom mvbutils do.on
#' @importFrom mvbutils extract.named
#' @importFrom mvbutils mlocal
#' @importFrom mvbutils my.all.equal
#' @importFrom mvbutils named
#' @importFrom mvbutils scatn
#' @importFrom stats dnorm
#' @importFrom stats na.omit
#' @importFrom stats pchisq
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom vecless ":="
#' @importFrom vecless compile_vecless
#' @importFrom vecless make_playback
#' @importFrom vecless set_recording
NULL




#' split FSPs from POPs
#'
#' Split FSPs from POPS; experimental, obsolete, and private!
#'
#' This version is not finished and probably abandoned; it will be superceded
#' by newer continuum-likelihood method. The reason for its existence, is that
#' (to my surprise) the 4way-based optimal-weight splitter for FSP/POP doesn't
#' work particularly well.
#'
#' Tries to split Full-Sibling Pairs from Parent-Offspring Pairs, using 6-way
#' genotypes. Takes advantage of the fact that POP-detection will also pick up
#' FSPs (as both relationship classes share at least one allele at every
#' locus), so FSP-detection can be run only within the subset of animals
#' identified as potential POPs, hence low computational demand.
#'
#'
#' @param snpg a \code{snpgeno} object
#' @param candiPOPs a 2-column matrix of of row-numbers in \code{snpg}, for
#' pairs known to be either POPs or FSPs (eg from \code{find_POPs}).
#' @param SDwt_POP scalar for determining "prior weighting" of POPs and FSPs.
#' 0.5 (the default) is probably OK.
#' @seealso simcheck_FSP_POP
#' @keywords internal
NULL





#' Don't call these yourself
#'
#' Some functions that need to be exported so that other related packages can
#' find them, but that you should not be messing about with.
#'
#'
#' @aliases Internals make_pgeno
#' @param pA an arg
#' @param pB an arg
#' @param pC an arg
#' @param which_genotypes arguments
#' @return Ooooh yes. Great value.
#' @keywords misc
NULL



