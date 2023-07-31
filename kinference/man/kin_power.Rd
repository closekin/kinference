\name{kin_power}
\alias{kin_power}
\title{Locus selection for kin-finding}
\description{\code{kin_power} can be used to predict how well a set of loci will work for HSP-finding, and to prepare for some QC and kinference steps on serious data. It returns the input \code{snpgeno} object with extra columns added to the \code{locinfo} attribute, related to the per-locus mean and variance of LOD (presumably an HSP/UP PLOD, though not inevitably) for different true kin-types. It respects the per-locus decision about how precisely to genotype (\code{useN=6/4/3}).
}
\usage{kin_power( lociar, want_LOD_table=TRUE, k, hack_LOD=NULL)
}
\arguments{
\item{  lociar}{\code{snpgeno} objects with the necessary ingredients}
\item{  want_LOD_table}{can't think why you'd set this to FALSE}
\item{  k}{target average kinship for LOD; 0.5 for HSPs, 0.25 for HTPs, etc.}
\item{  hack_LOD}{Don't mess around with this; it's for internal black magic}
}
\details{
\itemize{
\item E_UP, V_UP mean & variance for UPs
\item E_HSP, E_POP,E_FSP as you would expect
\item Ediff E_HSP - E_POP ie the "absolute" power of that locus
\item sdiff (E_HSP-E_POP)/sqrt(V_UP) which is arguably better than \code{Ediff} for ranking loci
}

It also attaches \code{LOD}, \code{PUP}, and \code{ev01} elements (each a matrix) to the \code{locinfo}. They have been made dull (see \code{make_dull}) to improve your viewing experience, but they work fine for all normal purposes (and you can always \code{unclass} them to remove the S3 class \code{dull}).
\subsection{Notes}{\code{kin_power} (and downstream) should get a refactor. It's daft to store LODs for only one specific kin; it'd be better to always calculate P1share and P2share as well as P0share (which is PUP), and then compute whatever-is-needed later on-the-fly. As-is, we are re-computing P1 and P2 based on LOD and PUP OTF instead (which is also unsafe, because LOD could have been calculated with k != 0.5).
}
}
\value{\code{snpgeno} object with augmented columns in "locinfo" attr.
}
\examples{
## Need some examples!
# See vignette for now
}
\keyword{misc}
