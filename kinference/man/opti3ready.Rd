\name{opti3ready}
\alias{opti3ready}
\title{Prepare for autopick_threshold}
\description{When choosing threshold to eliminate 3rd-order kin, \code{opti3ready} can be used first to recompute an optimal PLOD for HSP::HTP, rather than just using the HSP::UP PLOD. NOT WORKING YET!
}
\usage{opti3ready(x2, kin)
}
\arguments{
\item{ x2}{a \code{snpgeno}}
\item{ kin}{result of a previous call to \code{find_HSPs} (qv) on \code{x2}}
}
\details{Some tricky manoeuvering to calculate the needful...
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
