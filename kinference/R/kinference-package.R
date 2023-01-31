

#' The kinference package: data prep for close-kin mark-recapture
#' 
#' Preparatory functions for subsequent application of close-kin
#' mark-recapture, specifically:
#' 
#' \itemize{ \item finding close-kin pairs and duplicate samples amongst large
#' (i.e., many samples) multilocus-genotype datasets; \item QC of samples and
#' loci, ultimately for the same purpose. }
#' 
#' The genetic data currently handled is diploid biallelic[1] SNP genotypes.
#' Error rates should be low (so, no "3X coverage" etc!). Null alleles[3] are
#' allowed for, but "missing/unknown" genotypes are not tolerated; every sample
#' and locus must be typed, even if that leads to some errors (noting that
#' double-null is a legitamte call for \code{kinference}). The types of
#' close-kin considered[2] are POP, FSP, and 2nd-order kin (HSP, GGP, FTP),
#' which is the limit of resolution in the absence of genome-assembly data.
#' Version 1.x of kinference does not use the latter.
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
#' 
#' @name kinference-package
#' @aliases kinference kinference-package kinference
#' @docType package
#' @section See: Locus QC
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
#' @keywords misc
NULL





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
#' @keywords misc
NULL



