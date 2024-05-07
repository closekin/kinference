

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
#' Preparatory functions for subsequent application of close-kin mark-recapture
#' to a dataset from many individuals with already-called multilocus genotypes.
#' Specifically, things covered are
#' 
#' QC of samples and loci, with the specific goal of then being able to...
#' 
#' find close-kin pairs and duplicate samples.
#' 
#' These tasks break down into finer categories, and the various functions
#' under each are mentioned under \bold{See also} below.
#' 
#' The genetic data currently handled is diploid biallelic SNP genotypes. Error
#' rates should be low (so, no "3X coverage" etc!). Null alleles1 are allowed
#' for, but "missing/unknown" genotypes are not tolerated; every sample and
#' locus must be called, even if that leads to some errors (noting that
#' double-null is a legitimate call for \code{kinference}). The types of
#' close-kin considered are POP, FSP, and 2nd-order kin (HSP, GGP, FTP), which
#' is the limit of resolution in the absence of genome-assembly data. Version
#' 1.x of kinference does not use the latter.
#' 
#' The kin-finding process entails several steps, each of which needs to be
#' examined by a human being to make sure it has worked properly, before moving
#' to the next. Sometimes one has to go back to an earlier step, eg to
#' tighten/loosen quality control. The whole process should not be treated as
#' "automatic", and there is deliberately no
#' \code{kinference::shut_up_and_find_the_pairs()} function! For the same
#' reason, most of the functions do not have default parameter values (except
#' perhaps for visual purposes).
#' 
#' The starting point must always be a \code{snpgeno} object (see package
#' \pkg{gbasics}) containing already-called genotypes for each sample and
#' locus, plus the crucial sample-specific information ("metadata" to
#' geneticists, "covariates" to population dynamicists and statisticians!) such
#' as sampling-year, age, sex, etc depending on the dataset. The original
#' object gets augmented with extra data (eg allele and genotype frequency
#' estimates) as the steps proceed.
#' 
#' 1 "Null" here means: heritable and repeatable undetectable allele, eg due to
#' mutations at restriction site or indels nearby. It specifically excludes
#' dropout due to low coverage!
#' 
#' @name kinference-package
#' @aliases kinference kinference-package kinference-package kinference
#' @docType package
#' @author Mark V Bravington, David L Miller, Shane M Baylis
#' @seealso The vignette "kinference-vignette".
#' 
#' The following categories of "things the package does"; this ought to be
#' indexed better for R's help system, and the large number of semi-documented
#' functions that we don't actually want people to use nowadays should be
#' trimmed down! Anyway, these are probably the functions you'll need:
#' 
#' .LOCUS.QC'CHECK_6AND4'.(QV)
#' 
#' .SAMPLE.QC'ILGLK_GENO'.(QV),.'HETZMINOO_FANCY'.(QV),.'FIND_DUPLICATES'.(QV),.'FIND_DUPS_WITH_MISSING'.(QV),.'DROP_DUPS_PAIRWISE_EQUIV'.(QV)
#' 
#' .ALLELE.FREQUENCY.ESTIMATION'EST_ALF_ABO_QUICK'.(QV),.OR.LESS.LIKELY.'EST_ALF_6WAY'.(QV),.'EST_ALF_ABCO'.(QV),.'RE_EST_ALF'.(QV)
#' 
#' .PAIRWISE.KIN.FINDING'FIND_POPS'.(QV),.'FIND_POPS_LGLK'.(QV).WHEN.IT'S.WORKING.PROPERLY.(NOT.YET),.'FIND_HSPS'.(QV),.'HISTOPLOD'.FOR.GRAPHICS,.'KINPALETTE'.DITTO,.'SPLIT_FSPS_FROM_HSPS'.(QV).AND.OTHER.'SPLIT_X_FROM_Y'.FUNCTIONS,.'AUTOPICK_THRESHOLD'.(QV),.'VAR_PLOD_KIN'.(QV)
#' 
#' .PREDICTING.LOCUS.POWER.FOR.KIN.FINDING'KIN_POWER'.(QV)
#' 
#' .EXAMPLE.DATASETS'DROPBEARS',.'BLUEFIN'
#' @keywords misc
NULL



