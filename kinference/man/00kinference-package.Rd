\name{kinference-package}
\alias{kinference}
\alias{kinference-package}
\alias{kinference}
\docType{package}
\title{The kinference package: data prep for close-kin mark-recapture}
\description{
Preparatory functions for subsequent application of close-kin mark-recapture, specifically:

\itemize{
\item finding close-kin pairs and duplicate samples amongst large (i.e., many samples) multilocus-genotype datasets;
\item QC of samples and loci, ultimately for the same purpose.
}

The genetic data currently handled is diploid biallelic[1] SNP genotypes. Error rates should be low (so, no "3X coverage" etc!). Null alleles[3] are allowed for, but "missing/unknown" genotypes are not tolerated; every sample and locus must be typed, even if that leads to some errors (noting that double-null is a legitamte call for \code{kinference}). The types of close-kin considered[2] are POP, FSP, and 2nd-order kin (HSP, GGP, FTP), which is the limit of resolution in the absence of genome-assembly data. Version 1.x of kinference does not use the latter.

The kin-finding process entails several steps, each of which needs to be examined by a human being to make sure it has worked properly, before moving to the next. The process should not be treated as "automatic", and there is deliberately no \code{kinference::shut_up_and_find_the_pairs()} function!

The starting point must always be a \code{snpgeno} object (see package \pkg{gbasics}) containing already-called genotypes for each sample and locus, plus the crucial sample-specific information ("metadata" to geneticists, "covariates" to population dynamicists and statisticians!) such as sampling-year, , sex, etc depending on the dataset. The original object gets augmented with extra data (e.g., allele and genotype frequency estimates) as the steps proceed.
}
\section{See}{
Locus QC

Sample QC

Allele frequency estimation

Pairwise kin-finding statistics

Categorizing kinship of specific pairs

Predicting kin-finding power
}
\section{Vignettes}{
Intro

Full monty
}
\section{Types of genotype}{
Almost all functions in \code{kinference} v1.x expect biallelic (A,B) loci, but null alleles (O) are allowed for throughout. "Null allele" here denotes a variant/allele that is repeatable, heritable, and invisible; with sequencer data, it might for example arise from a mutation at the primer site, or a nearby indel that changes the fragment length. Individuals with a single copy of a null "look like" homozygotes (with one other non-null allele present). The presence of single-nulls inflates the number of \emph{apparent} "homozygotes" for a locus, so that a "null-naive"Hardy-Weinberg statistic will appear to fail. Individuals with two null copies (which will be rare) should be specifically typed as double-null (homozygote null). Null frequency can be estimated from the numbers of calls in the four categories [heterozygote; apparent ref homoz; apparent alt homoz; double-nulls], and a "null-robust" HW statistic can then be calculated; many loci that fail the "null-naive" HW test will pass a "null-robust" one, and can still be useful for kin-finding (although e.g. paralogous loci will still tend to fail the test, and should not be used). Double-nulls are conceptually quite different from a "missing/not sure" call which might result from, say, low read depths in a poor-quality sample. For \code{kinference}, all loci need to be called in all samples; there are good statistical and practical reasons for demanding that.

Genotypes can be encoded in different ways, and the \code{\link{snpgeno}} object in package \pkg{gbasics} allows the user to control for all loci in the dataset that via the \code{diplos} attribute (qv). It is entirely possible to define and use your own encoding (e.g. if there are more than two alleles), but in practice the only encodings that \code{kinference} supports directly are those listed in the code of \code{gbasics::define_genotypes} as \code{genotypes6}, \code{genotypes4_ambig}, and \code{genotypes3}. Most \code{kinference} functions in fact \emph{require} \code{genotypes6} encoding, even though it may not make sense for your dataset to separate single-null from homozygote calls. You still have to make an arbitrary split (see XXX) but then you can tell the \code{kinference} functions to ignore the split analytically!

The genotype for a locus can be handled at various degrees of resolution, and the \pkg{kinference} package allows this via a special locus-specific field \code{useN} described below. For example: with sequence data having high-enough read-depths, it becomes possible to distinguish single-null genotypes (AO) from homozygotes (AA) with \emph{some} reliability (not 100\%). However, at lower read-depths or with other non-sequence data, that might not be possible (or the misclassification rate might be high enough to render it futile). And although double-nulls can be informative about kinship (a double-null in one animal and a heterozygote in the other cannot be co-inherited), trusting double-null calls is not robust when null frequency is low, precisely because the kinship PLODs become highly sensitive to rare events that might actually be genotyping errors.

Encodings... NB all categories in an encoding must be mutually exclusive! Expect "genotypes6" encoding, but auto-conversion should be possible in an early step...

In \code{kinference} v1.x, most code is designed to work with the maximum "6-way" categories of genotype, whereby the called genotypes are (AB, AA, AO, BB, BO, OO), i.e. with single-null BO distinct from homozygote BB. This is really a legacy from early SNP work at CSIRO, where fewer loci were targeted but with higher read-depths, so that nulls became a 3rd allele and thus quite informative for kin-finding (even though there was an appreciable and locus-specific error-rate for homozygote-vs-single-nulls). Nowadays the economics (and the statistics) favours instead using more loci but at lower read-depths, so that is only worth using 4 categories of genotype: AB, AAO, BBO, OO where AAO means "either AA or AO, i.e. A is present but B is not". This works well unless the null frequency for a locus is very small, in which case occasional genotyping errors that result in OO calls can mess up the kinship statistics appreciably. For those loci, we therefore only use 3 categories: AB, AAO, BBOO where the latter \emph{merges} BB, BO and OO; this sacrifices a small amount of statistical power but improves robustness. The number of categories can vary between loci in a dataset, so the \code{useN} field (which needs to be added manually to \code{x$locinfo}, as explained under \code{\link{sngpeno}} in package \pkg{gbasics}) should be set to the values 6, 4, or 3 accordingly. As a general recommendation for new datasets, \code{useN=4} would be sensible unless the estimated null frequency is below 0.05, at which point \code{useN=3} is safer.
}
\section{Notes}{
[1] One or two functions can potentially handle multiallelic loci, but as yet not the core PLOD-calculation routines.

[2] POP = Parent-Offspring Pair; FSP = Full-Sibling Pair; HSP = Half-Sibling Pair; GGP = Grandparent-Grandchild Pair; HTP = Half-Thiatic Pair (e.g. half-uncle/niece, where one animal is the half-sibling of a parent of the other).

[3]
}
\keyword{misc}
