#' @importFrom mvbutils %is.not.a% %where%
#' @export
"get_chain" <- function( thing, seed) {
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

