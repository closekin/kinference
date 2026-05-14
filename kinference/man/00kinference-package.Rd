\name{kinference-package}
\alias{kinference}
\alias{kinference-package}
\alias{kinference-package}
\alias{kinference}
\docType{package}
\title{Kin-finding and QC for Close-Kin Mark-Recapture}
\description{
\code{kinference} covers the preparatory steps for Close-Kin Mark-Recapture, for a multilocus SNP dataset from many individuals that have already been genotyped. Specifically, it covers:

 allele frequency estimation, followed by...

 ... QC of samples and loci, with the goal of then being able to...

 ... find close-kin pairs and duplicated samples.

These tasks break down into finer categories, and the various functions under each are mentioned under \bold{See also} below. There is a vignette, \code{kinference-vignette}, which demonstrates many but not all features. For detailed explanations of the algorithms, see the manuscript under \bold{References}.

The kin-finding process entails several steps, each of which needs to be examined by a human being to make sure it has worked properly, before moving to the next. Sometimes one has to go back to an earlier step, eg to tighten/loosen quality control. The whole process should not be treated as automatic, and there is deliberately no \code{kinference::shut_up_and_find_the_pairs()} function! For the same reason, most of the functions do not have default parameter values, except perhaps for visual display; you are expected to think about appropriate choices for your dataset.

Note that \code{kinference} does not build or fit CKMR models per se. It "just" deals with the key preparatory step of finding kin-pairs, but the model per se is up to you!
\subsection{Genetic inputs}{
\code{kinference} handles diploid biallelic SNP genotypes for thousands of loci. Genotyping errors are tolerated, but error rate should be low: no "3X coverage" etc! Loci are assumed to be in HWE, which can be checked. Genuine null alleles (see below) are allowed for, but uncalled/unknown genotypes are \emph{not} permitted in most functions; every sample and locus must be called, even if some calls are wrong. The types of close-kin considered are POP, FSP, and 2KPs (2nd-order kin, ie HSP, GGP, FTP), which is the limit of resolution in the absence of genome-assembly data (outside the scope of \code{kinference} version 1).

\code{kinference} does \emph{not} call/score the genotypes- you need to do that beforehand. The starting point for \code{kinference} must always be a \code{snpgeno} object (see package \pkg{gbasics}) containing already-called genotypes for \emph{every} sample and locus, plus the crucial sample-specific information ("metadata" to geneticists, "covariates" to population dynamicists and statisticians!) such as sampling-year, age, sex, etc depending on the dataset. The original object gets augmented with extra data (eg allele and genotype frequency estimates) as the steps proceed.

There's no requirement for one particular genotyping method (\code{kinference} has been used successfully with sequencer data and with microarrays) but there \emph{are} limitations and expectations on the properties of the data. Some of these are mentioned here, and you may find more information in \code{kinference-vignette} (qv).
}
\subsection{Missing genotypes}{
You basically can't have any. Almost all functions in \code{kinference} expect a definite (albeit possibly wrong) genotype for every sample at every locus. Sometimes the call might be that neither the reference nor the alternative SNP are present; that's still a definite call, but it's treated as a double-null (see below) not as a "missing". The two exceptions are \code{\link{find_dups_with_missing}} (qv) and \code{\link{est_ALF_nonulls}}. The former is not intended for use with CKMR data, and the latter is a first step that might be useful if you need to "impute" (ie make up) missing genotypes.

The rationale for insisting on non-missingness, is that many of the statistics calculated by \code{kinference} have null distributions (i.e., under some null hypothesis about the true situation) that can be readily calculated \emph{iff} all samples are scored at all loci. Knowing the null distribution has been invaluable to us on many occasions, for seeing whether things have gone badly wrong. Trying to allow for missingness would put the calculations somewhere on the scale between unreliable, horribly complicated, or impossible.
}
\subsection{Genotyping errors}{The functions in \code{kinference} can tolerate some level of genotyping error. However, if your data has high levels of allelic dropout due to low coverage, then it's not suitable for \code{kinference}. And samples with dodgy or contaminated DNA won't work either; hopefully the QC routines will weed them out.
}
\subsection{Null alleles}{
Some genotyping methods, when applied to some species, produce a substantial proportion of null alleles. "Null" here means a \emph{heritable} and \emph{repeatable} but undetectable allele, eg due to mutations at restriction site or indels nearby. It specifically \emph{excludes} dropout due to low coverage, which is neither repeatable nor heritable.

Note that a sample might well have a \emph{single} copy of a null allele, in which case it will "look like" a homozygote for the other (non-null) allele. Only if a locus has a high null-allele frequency, will there be many double-nulls (i.e. where neither the "major/reference" nor the "minor/alternate" SNP is present). Thus, the evidence of nulls is not mostly from entirely-missing double-null calls, but rather from inflated "homozygote" proportions, and (if POPs are present) from apparent Mendelian exclusions within obvious POPs. The interpretation and handling of nulls is very important for kin-finding (perhaps less so for other genetic applications), so you really have to understand your genotyping method: some software is reluctant to call a genotype unless the result is extremely clear, leading to "missing" genotypes that are \emph{not} the same as double-nulls! In that case, force your software to damn well make the damn calls...

Of course, \code{kinference} can perfectly well cope with datasets \emph{without} nulls, too, as the vignette shows. Deciding what to do about nulls requires you to understand the properties of your genotyping method (software pipeline as well as technological underpinnings). There are several options for handling nulls in \code{kinference}, determined by how you "encode" your genotypes; see \code{gbasics::snpgeno}.
}
\subsection{Linkage and linkage disequilibrium}{
\code{kinference} is designed for datasets of a few thousand SNPs, where Linkage Disequilibrium (LD) is unlikely to be substantial; almost inevitably, a few loci \emph{will} be in LD just because they happen to be close together genomically, but that won't noticeably affect the results. If you have zillions of SNPs, eg from Whole-Genome-Sequencing (WGS), then LD \emph{will} be a problem and you won't be able to get \code{kinference} to work properly unless you thin your dataset a lot. However, a more significant problem with WGS is likely to be low coverage, leading to unacceptable levels of dropout...

Given its intended use, \code{kinference} quite reasonably ignores LD. However, \emph{linkage} per se- the fact that long consecutive stretches from the same strand of DNA are inherited at meiosis, rather than each locus being inherited independently- is absolutely \emph{not} ignorable for finding 2KPs (or FSPs). \code{kinference} has quite a bit of code that deals practically with linkage, as inferred from pairs that clearly are 2KP.
}
}
\references{Bravington, Mark and Baylis, Shane (2025): "Finding kin-pairs for Close-Kin Mark-Recapture". Journal of Joy, v1 p1.
\subsection{Historical note}{\code{kinference} has evolved considerably since its birth at CSIRO in ~2018. It was originally developed for ddRAD datasets back in the "Bronze Age", and a couple of those datasets are still in practical use within CSIRO. Thus, \code{kinference} has to support some legacy types of data that would best be avoided in new projects- most notably, the challenging (albeit successful) "6-way genotyping" that tried to explicitly distinguish single-nulls from homozygotes, based on read-depth. That, and a few other quirky cases, lead to some pretty complicated code in parts of \code{kinference}, and to several options (especially for null alleles) which might be confusing for new users who "just want to analyse their data". Such is life- sorry!
}
}
\seealso{
The vignette \code{kinference-vignette} (qv).

And for the categories of task, and the functions within each:
\subsection{Locus qc}{\code{\link{check6and4}}
}
\subsection{Sample qc}{\code{\link{ilglk_geno}}, \code{\link{hetzminoo_fancy}}, \code{\link{find_duplicates}}, \code{\link{find_dups_with_missing}}, \code{\link{drop_dups_pairwise_equiv}} (and, after kin-finding, \code{\link{chain_pairwise}}).
}
\subsection{Allele frequency estimation}{\code{\link{est_ALF_ABO_quick}}, \code{\link{est_ALF_nonulls}}, or less likely \code{\link{est_ALF_6way}}, \code{\link{est_ALF_ABCO}}, \code{\link{re_est_ALF}}
}
\subsection{Pairwise kin finding}{
\code{\link{find_POPs}} and \code{\link{find_POPs_lglk}}; \code{\link{find_HSPs}}; \code{\link{histoPLOD}} and 'k

inPalette' for graphics; \code{\link{split_FSPs_from_HSPs}} and other \code{split_x_from_y} functions; \code{\link{autopick_threshold}} and \code{\link{var_PLOD_kin}} for separating 2KP from more-distant kin. Also \code{\link{chain_pairwise}} and \code{\link{get_chain}} for post-hoc consistency checks on likely HSPs.
}
\subsection{Predicting locus power for kin finding}{\code{\link{kin_power}} (qv)
}
\subsection{Example datasets}{\code{\link{dropbears}}, \code{\link{bluefin}}
}
}
\keyword{misc}
