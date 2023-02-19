

#' Southern Bluefin Tuna data
#' 
#' An anonymised \code{snpgeno} dataset of Southern Bluefin Tuna genotypes for
#' 1038 individuals at 1510 loci
#' 
#' 
#' @name bluefin
#' @aliases bluefin data
#' @docType data
#' @format An object of class \code{gbasics::snpgeno()}
#' @keywords data
NULL





#' Red-rumped Dropbear data
#' 
#' An anonymised \code{snpgeno} dataset of Red-rumped Dropbear
#' (\emph{Thylarctos plummetus}, ssp. \emph{haemorhous}) genotypes for 480
#' individuals at 2000 loci
#' 
#' 
#' @name dropbears
#' @docType data
#' @format An object of class \code{gbasics::snpgeno()}
#' @author Shane M Baylis <email: shane.baylis@@csiro.au>
#' @keywords data
NULL





#' The kinference package: data prep for close-kin mark-recapture
#' 
#' Preparatory functions for subsequent application of close-kin
#' mark-recapture, specifically:
#' 
#' \itemize{ \item finding close-kin pairs and duplicate samples amongst large
#' (i.e., many samples) multilocus-genotype datasets; \item QC of samples and
#' loci, ultimately for the same purpose. }
#' 
#' The genetic data currently handled is diploid biallelic1 SNP genotypes.
#' Error rates should be low (so, no "3X coverage" etc!). Null alleles3 are
#' allowed for, but "missing/unknown" genotypes are not tolerated; every sample
#' and locus must be typed, even if that leads to some errors (noting that
#' double-null is a legitamte call for \code{kinference}). The types of
#' close-kin considered2 are POP, FSP, and 2nd-order kin (HSP, GGP, FTP), which
#' is the limit of resolution in the absence of genome-assembly data. Version
#' 1.x of kinference does not use the latter.
#' 
#' The kin-finding process entails several steps, each of which needs to be
#' examined by a human being to make sure it has worked properly, before moving
#' to the next. The process should not be treated as "automatic", and there is
#' deliberately no \code{kinference::shut_up_and_find_the_pairs()} function!
#' 
#' The starting point must always be a \code{snpgeno} object (see package
#' \pkg{gbasics}) containing already-called genotypes for each sample and
#' locus, plus the crucial sample-specific information ("metadata" to
#' geneticists, "covariates" to population dynamicists and statisticians!) such
#' as sampling-year, , sex, etc depending on the dataset. The original object
#' gets augmented with extra data (e.g., allele and genotype frequency
#' estimates) as the steps proceed.
#' 
#' SEE:.LOCUS.QC
#' 
#' Sample QC
#' 
#' Allele frequency estimation
#' 
#' Pairwise kin-finding statistics
#' 
#' Categorizing kinship of specific pairs
#' 
#' Predicting kin-finding power
#' 
#' @name kinference-package
#' @aliases kinference kinference-package kinference-package kinference
#' @docType package
#' @author Mark V Bravington, David L Miller, Shane M Baylis
#' @keywords misc
NULL



