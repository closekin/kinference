\name{upairid}
\alias{upairid}
\title{Label kin-pairs}
\description{\code{\link{find_HSPs}} etc return \emph{pairs} of kin, labelled by \code{i} and \code{j}. If you want to see whether one or more pairs also occur in another set of kin, eg from calling \code{\link{find_POPs_lglk}}, it's tempting to \code{paste} together each set's \code{i} and \code{j}, and then use \code{match}. But that will go fail if the order has switched, so that \code{i} in the first set corresponds to \code{j} in the second, and so on- which can definitely happen with \code{find_XXX} (for good reasons). To avoid that, \code{upairid} sorts each (i,j) pair into order before pasting them, so that you can use \code{match( upairid( first_set), upairid( second_set))} with no need to worry about order-switching.
}
\usage{upairid(i, j, sep='/')
}
\arguments{
\item{  i, j}{If \code{i} is a list (presumably from a call to \code{find_XXX}) and \code{j} is missing, then i. Otherwise, \code{i} and \code{j} can be vectors of numbers or characters.}
\item{  sep}{character (or string) to put between each \code{i} and \code{j}.}
}
\value{A character vector of pairwise-sorted, pasted, IDs.
}
\examples{
upairid( 1:5, 5:1)
d1 <- list( i=1:3, j=4:6)
d2 <- list( i=4:6, j=1:3) # swapped
match( paste( d1$i, d1$j), paste( d2$i, d2$j)) # nope
match( upairid( d1), upairid( d2)) # yep
}
\keyword{misc}
