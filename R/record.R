#' Self-documenting objects
#'
#' You might have a dataset created by a sequence of steps. If you wrap each step in a call to \code{record}, your final object will carry a record of how it came to be, in its \code{calls} attribute. Useful for self-documenting objects and "sequential processing" and "pipelines" and stuff.
#'
#' More elaborate use, you can arrange to skip parts of the sequence if they've already been, and to bail out at specific points.
#'
#' @param expr R{} expression to be evaluated. Use braces to wrap several expressions together. "Return-value" (i.e. result of last subexpression) is what gets assigned to \code{var}.
#' @param var name of the object to be returned--- presumably, something that's created or modified inside \code{expr}
#' @param name_of_call what name to give this particular "step" within the \code{calls} attribute
#'
#' @return The object created, with a \code{calls} attribute (a \code{list}) to which \code{expr} is appended.
#' @importFrom atease @
"record" <- function( expr, var, name_of_call=NULL, these_args=NULL) {
  expr <- substitute( expr)
  pf <- parent.frame()
  calls <- get( var, pf)@calls
  calls <- c( calls, expr)
  if( !is.null( name_of_call)) {
    names( calls)[ length( calls)] <- name_of_call
  }

  # expr <- substitute( { expr; get( var)})
  eval( expr, pf) # no side-effects
  res <- pf[[ var]]
  res@calls <- calls
  if( !is.null( these_args) && !is.null( name_of_call)) {
    res@args[[ name_of_call]] <- these_args
  }

  assign( var, res, envir=pf)
  res
}
