\name{kinPalette}
\alias{kinPalette}
\title{Colors for different kinships}
\description{
\code{kinPalette} allows consistent colours for different kinships when plotting; it's used by various \code{\link{kinference}} plotting functions, and you can also use it yourself to add lines, points, etc. The colours are taken from \code{viridisLite::viridis} (see \bold{References}); see the code for the hex values.

\code{kinPalette} returns a named vector of hex values, either for all kinships or just for those you specify. So you can do something like this:

\preformatted{
  # plot something...
  abline( v=17, col=kinPalette('UP') # use colour for UPs
}
\subsection{Good manners}{
Optionally, \code{kinPalette} will also call \code{grdevices::palette} for you, overwriting the existing numeric color definitions. Afterwards, \code{plot(...,col=2)} will show a different color. This isn't normally a good idea (you can always refer to the kin-colors by name; numbers are flaky) so the default is not to (from \code{\link{kinference}} v1.1 onwards). If you really want to do it yourself, then a perhaps-better approach is:

\preformatted{
  kp <- kinPalette()
  old_palette <- palette( kp)
  # on.exit( old_palette) # if inside a funciotn
}

so that any palette changes are temporary and can be undone.
}
}
\usage{kinPalette( kinships= names( kincolours), setPalette= FALSE)
}
\arguments{
\item{ kinships}{which kinships to return colors for. Default is all of them.}
\item{ setPalette}{set to TRUE if you really want to set graphics palette.}
}
\value{A character vector of hex codes, with names "POP" etc. If the \code{kinships} argument is specified, then just the corresponding elements will be returned.
}
\references{
S. Garnier, N. Ross, A. Camargo, B. Rudis, K. Woo, & M. Sciaini. (2023). sjmgarnier/viridisLite: CRAN release v0.4.2 (v0.4.2CRAN). Zenodo. https://doi.org/10.5281/zenodo.7890875

(NB I had to abbreviate the authors' first names to avoid stupid unnecessary warning from R{} CMD CHECK about non-ASCII characters... in a raw string... in the days of mandatory UTF-8... sigh)
}
\keyword{misc}
