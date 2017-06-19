\name{choose_geno6_thresholds}
\alias{choose_geno6_thresholds}
\alias{geno6way}
\title{Genotyping pre-processed DartCap}
\description{
\code{choose_geno6_thresholds} is to be run on a generous subset of the data, to figure out the thresholds for single nulls etc etc. Then, actually \bold{doing} the genotyping is handled by \code{geno6way}. You can call the latter separately, so that you can apply thresholds etc derived from one subset of the data to a new dataset, so that your genotyping is consistent.

Classifies each fish/locus (a biallelic SNP perhaps with nulls) into one of six categories: OO, AO, BO, AA, BB, AB. Input should come from previous steps of \code{check_baits}, either "stand-alone" to compute thresholds and error-rates, or applying preset thresholds-and-error-rates to fresh count data that's had some processing. \code{geno6} \emph{discards} the raw counts and returns instead a \code{snpgeno} object with associated locus information.

If in stand-alone mode and assuming \code{plots=TRUE}, then fits to each locus are shown graphically, and you get the chance to make a comment; any non-null comments are recorded to a file. The idea is that you can look at the file afterwards to see which loci you really didn't like.

Can also be run to "shut up and genotype" a dataset using previously-determined thresholds (eg if you have new plates of fish).
}
\usage{
geno6(lociar, OOthresh_tc = lociar@args$geno_deambig$OOthresh_tc,  nquants_bump, max_dat_quantile,
    nhistbreaks = 101, distro,  minpO_rethresh, m_for_rethresh, max_refits, plots = NULL,
    minbin, show_progress = interactive())
}
\arguments{
\item{ lociar}{a \code{loc.ar} from previous calls to \code{check_baits}, in p'tic with attributes}
\item{ OOthresh_tc:}{}
\item{ nquants_bump:}{}
\item{ max_dat_quantile:}{}
\item{ nhistbreaks:}{}
\item{ distro:}{}
\item{ minpO_rethresh:}{}
\item{ m_for_rethresh:}{}
\item{ max_refits:}{}
\item{ plots:}{}
\item{ minbin:}{}
\item{ show_progress:}{}
}
\value{
A \code{snpgeno} object with the original the \code{locinfo} attribute augmented by the :
\item{ snerr:}{}
\item{ best_cut:}{}
\item{ pbonzer:}{}
\item{ perr:}{}
}
\details{Not compulsory. Other section headings, e.g. AUTHOR, should also go here. Use \bold{single} quotes around object names and code fragments, e.g. \code{bit.of.code()}. Use \bold{double} quotes for "text" or "file.name". See \code{\link{doc2Rd}} for full details of format.
}
\seealso{\code{doc2Rd}, \code{flatdoc}
}
\examples{
# Not compulsory to have an EXAMPLES -- you can put examples into other sections.
# Here's how to make a "don't run" example:
\dontrun{
reformat.my.hard.drive()
}
}
\keyword{misc
}
