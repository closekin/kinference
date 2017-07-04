#' Find chains in HSPs
#'
#' Find chains of relatives of fish \code{seed} 
#'
#' @param thing output from \code{\link{find_HSPs}}, or some subset thereof
#' @param seed a fish ID
#' @importFrom mvbutils %is.not.a% %where%
#' @export
"get_chain" <-
function( thing, seed) {

  # if we didn't give the function a $bigs already
  # this needs to check other things are there?
  if( thing %is.not.a% 'data.frame') {
    thing <- thing$bigs
  }
  extract.named( thing)

  oset <- integer()
  set <- seed
  while( length( set) != length( oset)) {
    newj <- j[ i %in% set]
    newi <- i[ j %in% set]
    oset <- set
    set <- unique( c( set, newi, newj))
  }

  thing %where% (i %in% set | j %in% set)
}
