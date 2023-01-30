## tinytest

# Prune the 'geno2017.rda' dataset down to sensible samples that can be use
# This is NOT anything like a real (tiny) test yet !!!

require( tinytest) # for interactive debugging
require( gbasics)
require( kinference)

stopifnot(
    exists( 'geno2017', environment()),
    geno2017 %is.a% 'snpgeno'
  )

if( 'useN' %not.in% names( geno2017)){
  # ... ie if from before 3-way genotyping
  geno2017$locinfo$useN <- with( geno2017$locinfo,
     ifelse( use6,
       6L,
     ifelse( pbonzer[,'O'] + pbonzer[,'C'] < 0.05,
       3L,
       4L)))
}

# Zap controls
x <- geno2017[ !grepl( 'CR', geno2017$info$Our_sample),]
x <- x[ !duplicated( x$info$Our_sample),]
# Many fewer now :) 16063 * 1541


# find_duplicates() can be slow
# First, weed out DREADFUL samples

y <- ilglk_geno( x, hist_pars=list( xlim=c( -1800, -1200), nclass=180))
x <- x[ y > (-1600),] # 15855 * 1541

y <- hetzminoo_fancy( x, target='poor')
y <- hetzminoo_fancy( x, target='rich')

# I'd probably trim real data a bit harder than this, but they're not too bad

# find_duplicates with bin calcs is much slower than without
# After some experimentation...
# ... most "dups" have well under 100 pairs...

y <- find_duplicates( x, max_diff_loci=100, nbins=0, maxbin=200, limit_pairs=1000)
dd <- drop_dups_pairwise_equiv( y)
x <- x[ -dd,] # 15811 * 1541

# Better check
system.time( yy <- find_duplicates( x[1:1000,], max_diff_loci=100, nbins=10, maxbin=100, limit_pairs=1000))
system.time( yy <- find_duplicates( x[1:5000,], max_diff_loci=100, nbins=10, maxbin=100, limit_pairs=1000))
# 27 seconds on my laptop, and quadratically slower as we'd expect...
if( FALSE) yy <- find_duplicates( x, max_diff_loci=100, nbins=10, maxbin=100, limit_pairs=1000) # There aren't any.

# Unhappy memories of trying to decode Pete's random changes in nomenclature :/
ad <- grep( '(?i)^..IN', x$info$Our_sample)
ju <- grep( '^..PL', x$info$Our_sample)      # careful about case and position! "dupl"
length( ad) + length( ju) - nrow( x) # 0 phew

system.time( fp <- find_POPs( x, ad, ju[ 1:5000], keep_thresh=0.07)) # 16 POPs obvious below wpsex=0.04

# Let's extract a small subset with all parents all offs and a few more...

testad <- c( ad[ fp$i], rsample( 100, ad[-fp$i], replace=FALSE))
testju <- c( rsample( 200, ju[-fp$j], replace=FALSE), ju[ fp$j]) # so j-col should be 200 more than i-col

# New code...
fp1 <- kinference:::find_POPs_lglk( x, testad, testju, gerr=0.01, keep_thresh= -50, nbins=100)

# Let's look for AA POPs. There's likely a few AA HSP/GGPs too. Keeps dataset a bit smaller...

system.time( fp2 <- find_POPs_lglk( x, ad, ad, gerr=0.01, keep_thresh= -300, nbins=300))

xx <- prepare_PLOD_SPA( hsp_power( x, k=0.5))
system.time( hp1 <- find_HSPs( xx, ad, ad, keep_thresh=5, nbins=300))

# Running find_HSPs() on this makes it pretty clear which are HSPs and which are POPs (or possible FSPs, in extreme theory)
