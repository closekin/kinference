\name{kinference-package}
\alias{kinference}
\alias{kinference-package}
\alias{kinference}
\docType{package}
\title{The kinference package: data prep for close-kin mark-recapture}
\description{
Preparatory functions for subsequent application of close-kin mark-recapture to a dataset from many individuals with already-called multilocus genotypes. Specifically, things covered are

\itemize{
\item QC of samples and loci, with the specific goal of then being able to...
\item find close-kin pairs and duplicate samples.
}

These tasks break down into finer categories, and the various functions under each are mentioned under \bold{See also} below.

The genetic data currently handled is diploid biallelic SNP genotypes. Error rates should be low (so, no "3X coverage" etc!). Null alleles[1] are allowed for, but "missing/unknown" genotypes are not tolerated; every sample and locus must be called, even if that leads to some errors (noting that double-null is a legitimate call for \code{kinference}). The types of close-kin considered are POP, FSP, and 2nd-order kin (HSP, GGP, FTP), which is the limit of resolution in the absence of genome-assembly data. Version 1.x of kinference does not use the latter.

The kin-finding process entails several steps, each of which needs to be examined by a human being to make sure it has worked properly, before moving to the next. Sometimes one has to go back to an earlier step, eg to tighten/loosen quality control. The whole process should not be treated as "automatic", and there is deliberately no \code{kinference::shut_up_and_find_the_pairs()} function! For the same reason, most of the functions do not have default parameter values (except perhaps for visual purposes).

The starting point must always be a \code{snpgeno} object (see package \pkg{gbasics}) containing already-called genotypes for each sample and locus, plus the crucial sample-specific information ("metadata" to geneticists, "covariates" to population dynamicists and statisticians!) such as sampling-year, age, sex, etc depending on the dataset. The original object gets augmented with extra data (eg allele and genotype frequency estimates) as the steps proceed.

[1] "Null" here means: heritable and repeatable undetectable allele, eg due to mutations at restriction site or indels nearby. It specifically excludes dropout due to low coverage!
}
\seealso{
The vignette "kinference-vignette".

The following categories of "things the package does"; this ought to be indexed better for R{}'s help system, and the large number of semi-documented functions that we don't actually want people to use nowadays should be trimmed down! Anyway, these are probably the functions you'll need:
\subsection{Locus qc}{\code{check_6and4} (qv)
\subsection{Sample qc}{\code{\link{ilglk_geno}} (qv), \code{\link{hetzminoo_fancy}} (qv), \code{\link{find_duplicates}} (qv), \code{\link{find_dups_with_missing}} (qv), \code{\link{drop_dups_pairwise_equiv}} (qv)
}
\subsection{Allele frequency estimation}{\code{\link{est_ALF_ABO_quick}} (qv), or less likely \code{\link{est_ALF_6way}} (qv), \code{\link{est_ALF_ABCO}} (qv), \code{\link{re_est_ALF}} (qv)
}
\subsection{Pairwise kin finding}{\code{\link{find_POPs}} (qv), \code{\link{find_POPs_lglk}} (qv) when it's working properly (not yet), \code{\link{find_HSPs}} (qv), \code{\link{histoPLOD}} for graphics, \code{\link{kinPalette}} ditto, \code{\link{split_FSPs_from_HSPs}} (qv) and other \code{split_X_from_Y} functions, \code{\link{autopick_threshold}} (qv), \code{\link{var_PLOD_kin}} (qv)
}
\subsection{Predicting locus power for kin finding}{\code{\link{kin_power}} (qv)
}
\subsection{Example datasets}{\code{\link{dropbears}}, \code{\link{bluefin}}
}
}
}
\author{Mark V Bravington, David L Miller, Shane M Baylis}
\keyword{misc}
