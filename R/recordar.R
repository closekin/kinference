#' Quasi-automatic vectorization
#'
#' @description
#' Suppose you have a cooked up a complicated numerical function \code{fhard} that works for a "scalar" case (say, one locus at a time), using lots of intricate subset-via-lookups. And say it's vectorizable in principle, but the task is beyond you. (The real crunch comes if you want to use matrix-subsetting of matrix, or vector-subsetting of matrix; you cannot "just" vectorize those cases by adding dimensions at the start.) Then, after a few slight tweaks, you can apply the \code{recordar} tools to generate a Quasi-Auto-Vectorized version, that will run \bold{quickly} ie using proper R{} vectorization.
#'
#' The sequence of ops is:
#'
#' \itemize{
#' \item make some tweaks to \code{fhard}, specifically:
#' }
#'
#'  -- add a \code{\link{record}} arg to \code{fhard}, default FALSE
#'
#'  -- say which arrays/matrices/vectors/scalars in \code{fhard} need to have their ops recorded (ie could have extra dims at the start), using \code{set_recording}
#'
#'  -- postfix key numerical and definitional statements of \code{fhard} with \code{ ? 0} or sometimes \code{ ? 1}
#'
#' \itemize{
#' \item run \code{template <- fhard( <example>, record=TRUE)}
#' \item run \code{vfhard <- make_playback( fhard, template)}
#' \item then you are good to go with \code{vfhard( <realvectorargs>)}
#' }
#'
#' Any arguments of \code{fhard} that themselves need recording, must all be pre-vectorized in the same way during playback, I think... but it's OK to have other args that aren't.
#'
#' One trick lies in the redefinition of \code{?} to record an operation, rather than summon up help! This is the only way to get the tweak to work, because of operator precedence rules; if you hate it, you can instead wrap your statements in \code{recordar( <blah>)} but it makes the code harder to read.
# Easy tweaking via eg x[ cbind( a1, a2)] <- y[ lu] ? 0
#' @aliases set_recording make_playback
#' @usage recordar(assig, expand_dim = FALSE) # It's easier to tweak with '?' but you can instead wrap your statement with 'recordar(...)'
#' @usage set_recording( vars, record=TRUE) # put eg 'set_recording( c( "ar1", "vec2"))' in the first line in your 'fhard'
#' @usage make_playback( fhard, template, record_arg_name='record')
#'
#' @param assig (recordar) Statement to be recorded
#' @param expand_dim (recordar) ?In the vectorized version, should the prefixdims be explicitly prepended? Thanks to recording magic, this is is normally not needed (and wrong) if the RHS of the assignment includes recorded variables; expansion is needed only if the "scalar" version creates a fixed number of dimensions, eg \code{myrecar <- matrix( 1:6, 3, 2)}. Hence the usual tweak is to add \code{? 0} but occasionally you need \code{? 1} or \code{? TRUE} (whatever you find clearer).
#' @param vars (set_recording) Names of objects (numeric arrays, matrices, and/or vectors) that need to be tracked
#' @param record (set_recording) If FALSE, nothing is recorded, and your function executes as normal. So you would set this to be the \code{record} argument of \code{fhard}, if there is one, or leave blank if you always want to record.
#' @param fhard (make_playback) Your \code{fhard}, post tweaks.
#' @param template (make_playback) Result of one call to tweaked \code{fhard} in its recording mode
#' @param record_arg_name (make_playback) Your original function probably should have an argument named \emph{something} like \code{record}, which will be removed from the QAV version. If your argument isn't exactly named \code{record}, you can set \code{record_arg_name} instead. If you don't have any such argument (ie you always record), then set this to ''.
#'
#' @return \code{set_recording}: returns a function which is a tuned version of \code{recordar}, ready for your specific tasks. See \bold{Examples}. \code{recordar}: result of the expression. \code{make_playback}: a function, the quasi-auto-vectorized version of your
#' @importFrom atease @
#' @importFrom mvbutils %is.not.a%
"recordar" <- function( assig, expand_dim=FALSE) {
  assig <- substitute( assig)
  envo <- environment( sys.function())
  x <- eval.parent( assig)  # always needed
  if( is.null( dim( x))) {
    envo$last_dimnames <- names( x)
    envo$last_dim <- length( x)
  } else {
    envo$last_dimnames <- dimnames( x)
    envo$last_dim <- dim( x)
  }

  is_assign <- (assig %is.a% '<-') ||
                  ((assig %is.a% 'call') && (assig[[1]]==as.name( '<<-')))
  if( is_assign && (assig[[2]] %is.not.a% 'name')) { # eg x[] <- y
    # eval.parent( assig) will have created all indices for RHS
    subassig <- assig[[2]]
    eval.parent( subassig) # create index for LHS
  } else if( is_assign) { # x is being created
    name <- as.character( assig[[2]])
    if( name %in% names( envo$subs)) { # make x self-recording
      x@where <- envo
      x@whoami <- name
      oldClass( x) <- c( 'selfrecording_array', oldClass( x))
      assign( name, x, parent.frame())
      # to allow adding prefix-dims
      assig <- substitute( define( assig, expand_dim))
    }
  } # else do nothing special, eg nonreturning function call

  # How to grow a {}-object...
  envo$exprs[[ length( envo$exprs)+1]] <- assig

  return( x)
}
