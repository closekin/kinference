#' Support routines for storage ngs genotypes
#'
#' Support routines for the s3 class \code{ngs_count_ar} to store ngs genotypes, assumed to consist of lots of loci measured at lots of samples, with each locus having its own number of alleles (which might vary across loci), and a count of sequences recorded for each sample-locus-allele combination.
#'
#' \code{NGS_count_ar} generalizes the \code{loc.ar} family, which are good but basically only OK for two-alleles-with-nulls loci. To subset by animals and loci, you can write \code{x[animals,loci]} or \code{x[animals,loci,]}. If you want counts for specific alleles (and a fixed number of alleles across all loci, even for loci that don't \emph{have} that many alleles), then use \code{x[animals,loci,k]}, but see \bold{Value}.
#'
#' I suppose I should add a "constructor"; at the moment, the only place that an \code{NGS_count_ar} object is actually built, is inside \code{read_cluster_dart2}.
#'
#' This is all pre-calling-the-genotypes, keeping the actual counts. Presumably at some point there will be a called-genotype version with far more compact storage.
#'
#' Internal storage is actually sample*allele, but subsetting generally works at the locus level. The complication is that the sample*locus*alleles array is actually "ragged" (the number of alleles varies per locus), which means that \code{unclass} may not be that useful. Other info about each sample is in the attribute \code{x@info}, a \code{data.frame}, which is displayed along with the counts by the \code{print} method. Info about the loci themselves, and the sequences (alleles) themselves, is in \code{x@locinfo} and \code{x@seqinfo}, which are not shown by \code{print} (ie you have to extract them yourself). Note that I'm using \code{@} as in the \pkg{atease} package; these things are plain ol' attributes.
#'
#' @param x an \code{NGS_count_ar} object
#' @param i indices of sample, locus, and allele--- some may be empty, and \code{k} can be missing altogether.Integer, logical, and character subscripts are allowed, except that \code{k} cannot be character, and \code{j} cannot duplicate a locus. [I \emph{could} "fix" that if it was needed for simulation, I guess.]
#' @param j see above
#' @param k see above
#' @param trailing_dot default \code{FALSE}; set \code{TRUE} if you want to show that count data is non-integer (eg after norming by sample-total-reads) so that all counts end with a period. The post-decimal-place digits are not normally important, but if you really need to see them, you can do so by setting a \code{k} subscript.
#' @param dot_for_0 default \code{FALSE}; set \code{TRUE} if you want zero-counts replaced by a central dot.
#'
#' @return The subset method returns another \code{NGS_count_ar} object, \emph{unless} \code{k} is set. If \code{k} is looking for just one allele per locus, then the result is a 2D array of sample*locus; if \code{k} is looking for 2 or more alleles, the result is a 3D array of sample*locus*allele, with NA counts added when \code{k} refers to allele "number" that doesn't exist for that locus. These "pure array" results don't preserve the detailed sample or locus information (except via the \code{dimnames}) but can be useful for subsequent manipulation. \code{print} displays a data.frame-esque output, with the sample-info columns preceding a column-per-locus. The central-dot character (Latin-1 and Unicode 0xb7) is used to pad the formatting, for readability. In non-English locales, this might not display properly; I'd need expert advice to fix that, though.
#'
#' @alias print.NGS_count_ar dim.NGS_count_ar
#'
#' @examples
#' # x[ i, j, k] # S3 method for NGS_count_ar
#' # print( x, trailing_dot=getOption( 'trailing_dot_NGS_count_ar', FALSE),
#' #     dot_for_0=getOption( 'dot_for_zero_NGS_count_ar', FALSE), ...) # S3 method for NGS_count_ar
#' # dim( x) # S3 method for NGS_count_ar
#' @importFrom mvbutils cq named %without.name% FOR
#' @importFrom atease @ @<-
#' @importFrom utils head
"[.NGS_count_ar" <- function( x, i, j, k){
###########################################
# drop not allowed
  # Allow missing k (ie only 1 comma), but nothing else
stopifnot( nargs() >= 3)
  no_k <- nargs() == 3

  info <- x@info
  locinfo <- x@locinfo
  max_n_alleles <- max( locinfo$n_alleles)
  seqinfo <- x@seqinfo
  samps <- dimnames( x)[[1]]
  loci <- named( dimnames( x)[[2]])

  # Other stuff, eg: print for specialprint; genotypes if x is really read-counts
  other_atts <- attributes( x) %without.name% cq( dim, dimnames, info, seqinfo, locinfo)

  x <- unclass( x)
  if( missing( i)) {
    i <- 1 %upto% nrow( x)
  } else {
    i <- structure( seq_along( samps), names=samps)[ i] # ensure integer
  }

  if( !missing( j)) {
    if( !is.logical( j) && any( duplicated( j))) {
      stop( 'Duplicated loci not allowed')
    }

    j <- structure( 1 %upto% nrow( locinfo), names=locinfo$Locus)[ j]

    # Need to allow for changes in locus ordering (sigh...)
    msa <- match( seqinfo$Locus, locinfo$Locus[ j], 0) # 222200111
    o <- order( msa)
    seq_allele <- o[ msa[ o] > 0]

    if( !my.all.equal( seq_allele, seq_along( seqinfo))) {
      x <- x[ , seq_allele, drop=FALSE]
      seqinfo <- seqinfo[ seq_allele,,drop=FALSE]
      locinfo <- locinfo[j,,drop=FALSE]
      locinfo$end_col <- cumsum( locinfo$n_alleles)
      locinfo$start_col <- c( 1L, head( locinfo$end_col, -1)+1L)
    }
  }

  # do we REALLY need to do an i-subset?
  if( !my.all.equal( i, seq_along( samps))){
    x <- x[ i,,drop=FALSE]
    info <- info[ i,,drop=FALSE]
    samps <- samps[ i]
  }

  # Under certain circumstances, we *might* subscript by k too
  # Returns a pure 3D array
  # k >= n_alleles for that locus => NA
  if( !no_k || !missing( k)) {
    k <- (1 %upto% max_n_alleles)[ k] # logical to integer
# stopifnot( all( k) <= min( locinfo$n_alleles))
    xx <- do.call( 'cbind', FOR( 1 %upto% nrow( locinfo),
        x[ , (locinfo$start_col[.] : locinfo$end_col[.])[k],drop=FALSE]))
    if( length( k) > 1) {
      dim( xx) <- c( nrow( x), length( k), nrow( locinfo))
      xx <- aperm( xx, c( 1, 3, 2))
      dimnames( xx) <- list( samps, locinfo$Locus, NULL) # cannot assign names to the alleles, since they'll differ by locus
    } else { # save some work in the common case that length( k)==1
      dim( xx) <- c( nrow( x), nrow( locinfo))
      dimnames( xx) <- list( samps, locinfo$Locus)
    }
return( xx)
  }

  x@locinfo <- locinfo
  x@seqinfo <- seqinfo
  x@info <- info
  attributes( x) <- c( attributes( x), other_atts)
return( x)
}

