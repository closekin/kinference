\name{get_pair_covars}
\alias{get_pair_covars}
\title{Get pair covariate data}
\description{Gets sample covariate data for samples \code{i} and \code{j} in a pair data.frame returned by \code{\link{find_HSPs}}, \code{\link{find_POPs}}, etc. Sample covariate data are returned as one data.frame for each sample: pairs@i_covars contains the covariate data for sample \emph{i}, and pairs@j_covars for sample \emph{j}.
}
\usage{
get_pair_covars(
  snpg,
  pairs,
  fields = NULL,
  subset1 = 1 \%upto\% nrow(snpg),
  subset2 = subset1
)
}
\arguments{
\item{  snpg}{the \code{snpgeno} dataset from which the pairs were called}
\item{  pairs}{the output from a call to a \code{find_} or \code{split_} function from package \pkg{kinference}, or a row-wise subset of such an output.}
\item{  fields}{if NULL (the default), will return all covariate fields for each sample. Otherwise, should be a vector whose values match column names in \code{snpg@info}, and will return only the named covariates.}
\item{  subset1}{should be left as the default unless subset-comparisons were specified in the \code{find_} call that generated \code{pairs}, in which case should be set matching the subset-comparisons specified in that call}
\item{  subset2}{see subset1}
}
\value{the data.frame given in \code{pairs}, with additional data.frames pairs@i_covars and pairs@j_covars containing covariate data for sample \code{i} and \code{j}, respectively.
}
\keyword{misc}
