\name{kinPalette}
\alias{kinPalette}
\title{Set the kin palette}
\description{sets the palette to the recommended palette for the different kin types. Creates objects \code{POPcol}, \code{GGPcol}, etc., which are the hex codes for the recommended colours for POPs and GGPs, and so on for HCPs, FCPs, UPs, HTPs, FTPs, HSPs, and FSPs. Also sets the palette such that:
}
\usage{kinPalette(nlocal = sys.parent())
}
\arguments{\item{  nlocal}{the environment in which to define the kin colours}
}
\details{
\itemize{
\item 1 = POP
\item 2 = GGP
\item 3 = HCP
\item 4 = FCP
\item 5 = UP
\item 6 = HTP
\item 7 = FTP
\item 8 = HSP
\item 9 = FSP
}
}
\value{the named kin colours; see text.
}
\keyword{misc}
