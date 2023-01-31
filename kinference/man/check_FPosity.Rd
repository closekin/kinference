\name{check_FPosity}
\alias{check_FPosity}
\title{QC for kin-finding; private for now}
\description{Incomplete! Suppose to return predicted mean & variance of CLODs for each sample, ie how prone is that sample's particular genotype to yielding unusually high/low PLODs when compared with a random unrelated sample. Then you can turn this into a prediction of each sample's per-comp chance of yielding a False-Positive with another Unrelated sample. This looks like quite a powerful diagnostic, but is not fully explored yet. (Update: we used this when "chasing Weird Johnny and the Shoulder of Doom". It worked... but entirely \emph{not} in the way we expected, revealing a rather different cause..!) The document "d:/docs/genetics/Dart/sbt-baits-v3/too-many-plods.lyx" has more info in section 4.1 on "rat CLODs".
}
\usage{check_FPosity(snpg, nsim = 0)
}
\arguments{
\item{  snpg}{a \code{snpgeno} object.}
\item{  nsim}{currently inactive. A simulation option exists in the code to check the null distro (not much use for far tails, of course).}
}
\details{There is a bunch of code in the function connected with simulations and more elaborate calculations, currently commented out. So I'm "internalling" this for now. Let's hope you're not able to see this in the documentation...
}
\value{Dataframe with columns "ECLOD" and "VCLOD". See examples format.
}
\examples{
## Rough chance of yielding a PLOD>5, say
# cloddo <- check_FPosity( snpg = snpg)
# Pr_FPos_5 <- pnorm( 5, mean=cloddo$ECLOD, sd=sqrt( cloddo$CLOD), lower=FALSE)
# hist( Pr_Fpos_5, nc=50)
## highlight some known suspects
# abline( v=Pr_Fpos_5[ suspects], col='red')
}
\keyword{misc}
