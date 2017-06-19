\name{record}
\alias{record}
\title{Self-documenting objects}
\description{
You might have a dataset created by a sequence of steps. If you wrap each step in a call to \code{record}, your final object will carry a record of how it came to be, in its \code{calls} attribute. Useful for self-documenting objects and "sequential processing" and "pipelines" and stuff.

More elaborate use, you can arrange to skip parts of the sequence if they've already been, and to bail out at specific points.
}
\usage{record(expr, var, name_of_call = NULL)
}
\arguments{
\item{ expr}{R{} expression to be evaluated. Use braces to wrap several expressions together. "Return-value" (i.e. result of last subexpression) is what gets assigned to \code{var}.}
\item{ var}{name of the object to be returned--- presumably, something that's created or modified inside \code{expr}}
\item{ name_of_call}{what name to give this particular "step" within the \code{calls} attribute}
}
\value{The object created, with a \code{calls} attribute (a list) to which \code{expr} is appended.
}
\examples{
# Not compulsory to have an EXAMPLES -- you can put examples into other sections.
# Here's how to make a "don't run" example:
\dontrun{
\{ mydata <- mydata[1:3,1:4] \}
\{ mydata[1:3, 1:4] \}
\{ something_else <- mydata[1:3,1:4] \} # OKreturns the value of something_else
\{ mydata <- mydata[ 1:3, 1:4]; cat( 'OK\\n') \} # NOT OK! doesn't return mydata
\{ mydata <- mydata[ 1:3, 1:4]; cat( 'OK\\n'); mydata \} # OK
reformat.my.hard.drive()
}
}
\keyword{misc
}
