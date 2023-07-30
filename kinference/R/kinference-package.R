

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





#' Kin-finding power for microhaplotyped loci
#' 
#' This is a short-term fudge for checking HSP-finding power of a bunch of loci
#' that (i) can have as many haplotypes as you like, but (ii) have no errors or
#' nulls. See \emph{Examples} for how you might use it.
#' 
#' If you want to explore the impact of missing genotypes (so that e.g. only
#' 90\ reasonable and very easy option is to multiply \code{Ediff} and
#' \code{V.UP} both by 0.9, then go thru the steps. If you choose the 0.9
#' conservatively- ie it's highly likely that >0.9 of loci get co-scored- then
#' the above calc avoids any need to do much more complicated stuff (which I
#' leave to you...).
#' 
#' At some point in future, \code{\link{kinference}} might be changed so that
#' it can handle >2 non-null alleles gracefully (ie microhaplotypes). But not
#' yet. So for now this version does some ghastly "live-hacking" of existing
#' code for \code{\link{kin_power}} to implement no-errors no-nulls
#' multi-allelic case. It will be hard to follow, so use \code{mtrace} if you
#' really want to see what's going on. The guts of the code is in
#' \code{\link{kin_power}} and \code{predict_hsp_util}.
#' 
#' @param lociar Usually, a matrix of allele frequencies (Locus * Alleles).
#' Locus names are set from the rownames, or "L1", "L2" etc if there are no
#' rownames. Allele names will be set to "A", "B", "C", etc, regardless of
#' colnames; you do not have a choice there. Will be renormalized so rows sum
#' to unity. NB \code{lociar} can also be a \code{snpgeno} object, as expected
#' for \code{kin_power}. If so, then the allele freqs are assumed to live in
#' \code{lociar$locinfo$pbonzer}, and \emph{no} nulls or genotyping errors are
#' allowed for; hence, for a DartCap "ABCO"-style dataset, \code{kin_power} and
#' \code{kin_power2} will give \emph{different} answers.
#' @return If \code{lociar} is an allele-frequency matrix, then you get a
#' dataframe with one row per locus and columns "Ediff", "V.UP", and "sdiff".
#' "Ediff" is "ELOD|HSP - ELOD|UP"; "V.UP" is "VLOD|UP"; "sdiff" is
#' \code{sqrt(V.UP)/Ediff}, useful for ranking locus power. See \emph{Examples}
#' for use.
#' @seealso \code{\link{kin_power}}
#' @keywords misc
#' @examples
#' 
#' ALF <- matrix( runif( 15), 3, 5) # 3 loci; 5 alleles
#' POW <- kin_power2( ALF)
#' # look at the contents of each...
#' # Now do it for lots of loci. NB the allele freqs above are *insanely* good; you won't
#' # find anything like that in practice for lots'n'lots of loci
#' lots <- 500
#' ALF <- matrix( runif( lots*5), lots, 5)
#' POW <- kin_power2( ALF)
#' # Now say we plan 1e6 pairwise comps, and might expect 100 HSPs
#' # Work relative to E[LOD|UP] which is not returned explicitly; treat that as "origin" ie 0
#' V <- sum( POW$V.UP)  # V[PLOD|UP]
#' E <- sum( POW$Ediff) # E[PLOD|HSP] - E[PLOD|UP]
#' E / sqrt( V) # 10.5 SDs--- pretty good.  Mean of HSPs is 10.5 UP-SDs above mean of UPs,
#' # ... so v. unlikely an UP will get as far as _typical_ HSP. But we need to be a bit
#' # ... more stringent than "typical"--- and, NB weaker kin
#' bigUP <- qnorm( 1e-6, mean=0, sd=sqrt( V), lower=FALSE) # most-kinlike UP
#' smallHSP <- qnorm( 1e-2, mean=E, sd=sqrt( 4*V))      # least-kinlike HSP
#' # ... so that's probably OK...
#' 
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



