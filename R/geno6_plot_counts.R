## TODO
# - take png() call out of this?

#' @importFrom mvbutils cq do.on
#' @importFrom stats pbinom pchisq
#' @importFrom graphics axis hist lines abline layout title plot legend text strwidth mtext
#' @importFrom atease @
#' @importFrom grDevices png dev.off
"geno6_plot_counts" <- function( nlocal=sys.parent()) mlocal({
    # Only do the binning once (slow); may *plot* twice
    # All counts
    h <- hist(
      tot2[,i],
      breaks=seq(0, max( tot2[,i])*1.05, length=nhistbreaks),
      plot=FALSE)
    # Just Ref & Alt--- a bit misleading in case of 3rd-allele-splatter
    # when really AC etc
    h2 <- hist(
      (cA+cB)[,i],
      breaks=seq(0, max( c( tot2[,i], (cA+cB)[,i]))*1.05, length=nhistbreaks),
      plot=FALSE) # bizarre breaks since cA+cB *can* exceed tot2... if..?

  if( !is.null( comment_file) && !file.exists( comment_file)) {
    cat( 'Locus: Comment\n', file=comment_file) # it shall be thus!
    # Any comment where the first non-space is 4, will be treated as
    # a 4-way locus
    # Other comments lead to discarding the locus
  }

  # Display logic is probably nonsensical now
  on_screen <- show_progress # TRUE # might want to save to file

  repeat{ # in case we want to re-do & save the plot to file
    # Add hists of pA/(pA+pB); two on top row, then big count hist under
    layout( matrix( c(
      1, 2,
      3, 3), 2, 2, byrow=T),
      heights=c( 1, 3))

    for( isubbo in c( -1, +1)) {
      pick <- which( isubbo * (cA[,i] + cB[,i]) > isubbo) # first <1, then >1
      if( length( pick)) {
        hist( cA[ pick, i] / (cA[ pick, i] + cB[ pick, i]), xlim=0:1, nc=40,
             xlab=NULL, ylab=NULL, main=NULL)
        title( sprintf( 'Ppn A/(A+B): %s', if( isubbo>0) 'HI' else 'LO'),
              line=-1, cex=0.5)
      } else { # skip the plot; nothing to see here
        plot( 0, 0, type='n', xlab='', ylab='', main='', axes=FALSE)
        # frame()
      }
    }

  # plot.histogram is fugly and funcontrollable, so plot separately

    plot( h, freq=FALSE,
        col='black', border='black',
        main= sub( ':[^:]*$', '', alid[i]),
        ylim=c( 0, max( h$density[-1])*1.1), # avoid domination by OO, which
                                             # is very concentrated but
                                             # just one bar
        axes=FALSE,  xlab='', ylab='', cex.main=1)
    axis( side=1)

    # Vertical lines to show A+B only
    with( h2, lines( breaks[-1]-mean(diff(breaks))/2, density, type='h',
         col='magenta'))

    x <- seq( 0, max( dat), length=101)
    lines(x, (1-phat[OO]) * dnorm(x, res$mu[1], res$sigma[1])*res$lambda[1],
          col='red')
    lines(x, (1-phat[OO]) * dnorm(x, res$mu[2], res$sigma[2])*res$lambda[2],
          col='green')
    # OK to have red & green since they don't need to be distinguished!
    abline(v=res$mu,col=c( 'red', 'green'), lty=2, lwd=2)
    # abline(v=mu.het,col='orange',lwd=2)
    # abline( v=OOthresh_tc, lwd=2, col='purple') # cutoff for OO
    abline( v=new_OOthresh[i], col='purple', lwd=2)
    abline( v=best_cut, col='cyan', lwd=2) # two of them now
    abline( v=hix, col='yellow2', lwd=2, lty=2)

    # Obs & pred geno4
    g4_pred <- gtab4( n_fish * phat)
    gamb_obs <- do.on( names( phat), sum( geno[,i]==.))
    g4_obs <- gtab4( gamb_obs)

    # ... overall chisq for title
    gofstat4 <- 2 * g4_obs %*% log( ifelse( g4_obs>0, g4_obs/g4_pred, 1))
    this_pval4 <- pchisq( gofstat4, df=1, lower.tail=FALSE)

    # phat itself includes C
    pb <- pbinom( g4_obs, prob=g4_pred[ names( g4_obs)]/n_fish, size=n_fish)
    # Don't flag excessively GOOD fits when 0 observed & predicted!
    orright4 <- (pb %in.range% c( minbin, 1-minbin)) | (g4_obs==0 & pb > minbin)

    # Hideous trick to right-align legends; from ?legend
    legend_items <- sprintf( '%3s %4i  %4.0f.', names( g4_obs), g4_obs, g4_pred)

    # Might as well do geno6 while we're at it... though this duplicates
    # post-hoc code in geno6
    # paranoid post-subscript for order
    p6 <- calc_g6probs( pA, pB, pC, snerr=snerr)[ genotypes6]
    g6_pred <- n_fish * p6
    g6 <- split_geno_4to6( geno[,i], cA[,i], cB[,i], info=lociar@info, li=li,
        best_cut= best_cut)
    g6_obs <- do.on( as.raw( match( genotypes6, g6@diplos)), sum( g6==.))

    # ... overall chisq for title
    gofstat6 <- 2 * g6_obs %*% log( ifelse( g6_obs>0, g6_obs/g6_pred, 1))
    this_pval6 <- pchisq( gofstat6, df=3, lower.tail=FALSE)

    legend_items <- c( legend_items,
        '', '',
        sprintf( '%3s %4i  %4.0f.', genotypes6, g6_obs, g6_pred))
    pb <- pbinom( g6_obs, prob=p6, size=n_fish)
    # Don't flag excessively GOOD fits when 0 observed & predicted!
    orright6 <- (pb %in.range% c( minbin, 1-minbin)) | (g6_obs==0 & pb > minbin)

    temp <- legend('topright',
                   legend=rep(' ', length( legend_items)), # presumably 4 but...
                   xjust=1, yjust=1,
                   text.width=1.05*max( strwidth( legend_items, family='mono')))
    text( temp$rect$left + temp$rect$w, temp$text$y, pos=2, family='mono',
        col= ifelse( c( orright4, TRUE, TRUE, orright6), 'black', 'orange'),
        legend_items)

    title(sprintf(
            '%s:   MHT= %i   PVAL4= %s  PVAL6= %s   Pr[C]= %4.2f NPOLYFISH= %i',
            locus_name[ i], as.integer( mht[i]), formatC( this_pval4, 2),
            formatC( this_pval6, 2), pC, li$npoly[i]), line=-2)
    if( nzchar( comments[i])) { # only on 2nd pass, if comment was made
      mtext( comments[i], side=3, outer=TRUE, col='orange2')
    }

    if( on_screen) {
      cat( 'Whaddya reckon? (Silence => assent) :')
      comments[ i] <- sub( '^ *$', '', readLines( n=1)) # ignore all-space lines
      if( nzchar( comments[i])) {
        # write it to a file NOW, fool!
        if( !is.null( comment_file)) {
          cat(locus_name[ i], ': ', comments[ i], '\n', file=comment_file,
              append=TRUE)
        }
        on_screen <- FALSE
        png(filename=sprintf( '%s-summary.png', locus_name[ i]),
            width=1024, height=768)
        par( oma=c( 0, 0, 3, 0))
next # re-do and save plot
      } # if something is spotted
    } else { # have done & saved 2nd round of plots
        dev.off() # finish saved plot. Possibly a mistake if in
                  # "just-make-the-graphs" non-interactive mode since it
                  # generates alternate empty plots
    }

  break # normally don't save, so don't repeat
  } # ... possible repeat to allow saving graph

})
