\name{kinPalette}
\alias{kinPalette}
\title{Set the kin palette}
\description{
Sets the \code{palette} to recommended values for different kinships, so that eg \code{point( ..., col=2)} will be consistent with colours already on graphs produced by other \code{\link{kinference}} package functions such as \code{\link{histoPLOD}}. Hex values are returned in a vector, which you can use later use like so:

\preformatted{
  kincol <- kinPalette()
  # plot something...
  abline( v=17, col=kincol[ 'UP']) # for whatever colour corresponds to UPs
}

See \bold{Value} for the predefined kinships.
}
\usage{kinPalette()
}
\value{
A named character vector of hex codes, with names "POP" etc. The current numerical ordering (not that you should ever need to know it- just use the names) is:
\item{  c(:}{}
\item{  POP= 1,:}{}
\item{  GGP= 2,:}{}
\item{  HCP= 3,:}{}
\item{  FCP= 4,:}{}
\item{  UP= 5,:}{}
\item{  HTP= 6,:}{}
\item{  FTP= 7,:}{}
\item{  HSP= 8,:}{}
\item{  FSP= 9:}{}
\item{  ):}{}
}
\keyword{misc}
