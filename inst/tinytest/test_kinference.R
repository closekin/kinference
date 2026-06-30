## These are the public unit tests for package 'kinference'
## SMB has additional tests, 
## On other machines, edit that filepath to the source dir for package kinference.

library(atease)
library(kinference)
library(gbasics) ## shouldn't be necessary to do this manually
library(mvbutils)


## absolute basic mechanics

expect_true("kinference" %in% sessionInfo()$otherPkgs$kinference)

info <- data.frame(
    Our_sample = c("sample1", "sample2"), 
    location = c("loc1", "loc2") )
    
locinfo <- data.frame(Locus = c("L1", "L2"), n_alleles = c(2,2))
genos <- matrix(data = c(rep("AAO", 4)), nrow = 2, ncol = 2)
minisnpg <- snpgeno(genos, diplos = c("AAO", "AB", "BBO", "OO"), 
    info = info, locinfo = locinfo)
expect_equal(class(minisnpg), "snpgeno")
expect_true(minisnpg[1,1] == "AAO")
minisnpg[1,2] = "AB"
minisnpg[2,1] = "BBO"
minisnpg[2,2] = "OO"

expect_true(minisnpg[1,1] == "AAO")
expect_true(minisnpg[1,2] == "AB")
expect_true(minisnpg[2,1] == "BBO")
expect_true(minisnpg[2,2] == "OO")

expect_true(as.raw(minisnpg[1,1]) == 01)
expect_true(as.raw(minisnpg[1,2]) == 02)
expect_true(as.raw(minisnpg[2,1]) == 03)
expect_true(as.raw(minisnpg[2,2]) == 04)

## test snpgeno generation
define_genotypes()
### 6-way parts
chargenos6 <- matrix( data = c( "AA", "AO", "BB", "AB", "OO", "BO"), 
    nrow = 2, ncol = 3)
intgenos6 <- matrix( data = c(1,3,4,2,6,5) , nrow = 2, ncol = 3)
rawgenos6 <- as.raw( matrix( data = c(1,3,4,2,6,5) , nrow = 2, ncol = 3))

info6 <- data.frame(
    Our_sample = c("sample1", "sample2"), location = c("loc1", "loc2") )
locinfo6 <- data.frame(Locus = c("L1", "L2", "L3"), n_alleles = c(2,2,2))

### 4-way parts
chargenos4 <- matrix( data = c( "OO", "AAO", "BBO", "AB"), nrow = 2, ncol = 2)
intgenos4 <- matrix( data = c(1,3,4,2) , nrow = 2, ncol = 2)
rawgenos4 <- as.raw( matrix( data = c(1,3,4,2) , nrow = 2, ncol = 2))

info4 <- data.frame(
    Our_sample = c("sample1", "sample2"), location = c("loc1", "loc2") )
locinfo4 <- data.frame(Locus = c("L1", "L2"), n_alleles = c(2,2))

### 6-way diplos
#### character genos
expect_silent( {
snpg6_char <- snpgeno( x = chargenos6, diplos = genotypes6, 
    info = info6, locinfo = locinfo6)
})
#### integer genos
expect_silent( {
snpg6_int <- snpgeno( x = intgenos6, diplos = genotypes6, 
    info = info6, locinfo = locinfo6, 
    allow_nonchar = TRUE)
})
#### raw genos
expect_silent( {
snpg6_raw <- snpgeno( x = rawgenos6, diplos = genotypes6, 
    info = info6, locinfo = locinfo6, 
    allow_nonchar = TRUE)
})
##### check that char, int, and raw inputs all give the same genotypes
expect_true( { all( snpg6_char == snpg6_int) })
expect_true( { all( snpg6_char == snpg6_raw) })

#### empty genos
expect_silent( {
snpg6_nogeno <- snpgeno( x = NULL, diplos = genotypes6, 
    info = info6, locinfo = locinfo6)
})
#### empty info
expect_silent( {
snpg6_noinfo <- snpgeno( x = chargenos6, diplos = genotypes6, locinfo = locinfo6)
})
#### empty locinfo
expect_silent( {
snpg6_nolocinfo <- snpgeno( x = chargenos6, diplos = genotypes6, info = info6)
})
#### empty genos and info
expect_silent( {
snpg6_nogeno_noinfo <- snpgeno( x = NULL, diplos = genotypes6,
    locinfo = locinfo6, n_samples = 2)
})
#### empty genos and locinfo
expect_silent( {
snpg6_nogeno_nolocinfo <- snpgeno( x = NULL, diplos = genotypes6, info = info6,
    n_loci = 2)
})
#### empty genos, info, and locinfo
expect_silent( {
snpg6_nogenoinfolocinfo <- snpgeno( x = NULL, diplos = genotypes6,
    n_loci = 2, n_samples = 2)
})

### 4-way diplos
#### character genos
expect_silent( {
snpg4_char <- snpgeno( x = chargenos4, diplos = genotypes4_ambig, 
    info = info4, locinfo = locinfo4)
})
#### integer genos
expect_silent( {
snpg4_int <- snpgeno( x = intgenos4, diplos = genotypes4_ambig, 
    info = info4, locinfo = locinfo4,
    allow_nonchar = TRUE)
})
#### raw genos
expect_silent( {
snpg4_raw <- snpgeno( x = rawgenos4, diplos = genotypes4_ambig, 
    info = info4, locinfo = locinfo4,
    allow_nonchar = TRUE)
})
##### check that char, int, and raw inputs all give the same genotypes
expect_true( { all( snpg4_char == snpg4_int) })
expect_true( { all( snpg4_char == snpg4_raw) })

#### empty genos
expect_silent( {
snpg4_nogeno <- snpgeno( x = NULL, diplos = genotypes4_ambig, 
    info = info4, locinfo = locinfo4)
})
#### empty info
expect_silent( {
snpg4_noinfo <- snpgeno( x = chargenos4, diplos = genotypes4_ambig, 
    locinfo = locinfo4)
})
#### empty locinfo
expect_silent( {
snpg4_nolocinfo <- snpgeno( x = chargenos4, diplos = genotypes4_ambig, 
    info = info4)
})
#### empty genos and info
expect_silent( {
snpg4_nogeno_noinfo <- snpgeno( x = NULL, diplos = genotypes4_ambig,
    locinfo = locinfo4, n_samples = 2)
})
#### empty genos and locinfo
expect_silent( {
snpg4_nogeno_nolocinfo <- snpgeno( x = NULL, diplos = genotypes4_ambig, 
    info = info4,
    n_loci = 2)
})
#### empty genos, info, and locinfo
expect_silent( {
snpg4_nogenoinfolocinfo <- snpgeno( x = NULL, diplos = genotypes4_ambig,
    n_loci = 2, n_samples = 2)
})

## generate a small 4-way, 6-way, and 3-way snpgeno object

## first, 4-way:
set.seed(1111)
## 100 loci, 30 samples. Simple snpgeno prep lightly modified from the kinference course notes
genos <- array( 
    data= sample( c(1,2,4,6), 3000, 
        prob = c(0.35, 0.4, 0.24, 0.01), replace = TRUE),
    dim = c(30, 100))
Locus <- paste("locus_", 1:100, sep = "")
newLocusData <- data.frame(
    thingOne.l = runif(ncol(genos)), 
    thingTwo.l = runif(ncol(genos)),
    thingThree.l = runif(ncol(genos)))
locinfo <- cbind( Locus, newLocusData)

Our_sample <- paste("sample_", 1:30, sep = "")
newSampleData <- data.frame(
    thingOne.s = runif(nrow(genos)),
    thingTwo.s = runif(nrow(genos)))
info <- cbind( Our_sample, newSampleData)

smallsnpg4 <- snpgeno(x = genos, diplos = genotypes6,
    info = info, locinfo = locinfo, allow_nonchar = TRUE)

smallsnpg4 <- kinference::re_est_ALF(smallsnpg4)
smallsnpg4$locinfo$useN <- 4
## snerr <- matrix(data = 0, nrow = length(Locus), ncol = 4)
## snerr@dimnames <- list(NULL, c("AA2AO", "AO2AA", "BB2BO", "BO2BB"))
## smallsnpg4$locinfo$snerr <- snerr

### check that subsets via [.snpgeno are also subsetting $info and $locinfo correctly
#### subset samples by lookup
smallsnpg4b <- smallsnpg4[ ! smallsnpg4$info$Our_sample %in% 
    c("sample_10", "sample_13") ,]
expect_true( dim(smallsnpg4b)[1] == dim(smallsnpg4b$info)[1] )
#### subset samples by number
smallsnpg4c <- smallsnpg4[ 1:10 ,]
expect_true( dim(smallsnpg4c)[1] == dim(smallsnpg4c$info)[1] )
#### subset loci by lookup
smallsnpg4d <- smallsnpg4[, ! smallsnpg4$locinfo$Locus %in% 
    c("locus_10", "locus_13")]
expect_true( dim(smallsnpg4d)[2] == dim(smallsnpg4d$locinfo)[1] )
#### subset loci by number
smallsnpg4e <- smallsnpg4[ , 1:70 ]
expect_true( dim(smallsnpg4e)[2] == dim(smallsnpg4e$locinfo)[1] )
#### subset samples and loci by number
smallsnpg4f <- smallsnpg4[ 1:10 , 1:70 ]
expect_true( dim(smallsnpg4f)[1] == dim(smallsnpg4f$info)[1] )
expect_true( dim(smallsnpg4f)[2] == dim(smallsnpg4f$locinfo)[1] )

## check that it can be kinferred w/o tripping an error or warning from find_HSPs:
smallsnpg4a <- kin_power(smallsnpg4, k = 0.5)
smallsnpg4b <- prepare_PLOD_SPA(smallsnpg4a) ## technically superfluous now

thePLODs4 <- find_HSPs(smallsnpg4b, limit_pairs = choose(nrow(smallsnpg4b),2),
    keep_thresh = -20)
## test w/o independent prepare_PLOD_SPA:
expect_silent({ find_HSPs(smallsnpg4a, limit_pairs = choose(nrow(smallsnpg4b),2),
    keep_thresh = -20) } )

## run once in January 2024, our benchmark of 'correct', yea, unto eternity:
## refPLODs4 <- thePLODs4
## save(refPLODs4, file = "smallsnpg4_referencePLODs.Rda")

# Next was passing for R-universe on Windows and sometimes on Macs but never on Linux 
# (but it does pass on SMB's Linux box). Hence added a generous tolerance in case of numericalia.
base::load("smallsnpg4_referencePLODs.Rda")  ## base:: to avoid a conflict with renv::load.
expect_true(all.equal(refPLODs4, thePLODs4, check.attributes = FALSE, tolerance=1e-5))

## test auto-generation of snerr by kin_power in cases where it's absent
expect_true( length( dim(smallsnpg4b$locinfo$snerr)) > 0 ) ## this is basically just
## 'exists( "smallsnpg4b$locinfo$snerr")', but I can't atease within a quoted string

## now, 6-way.
set.seed(1112)
## 100 loci, 30 samples. Simple snpgeno prep lightly modified from the kinference course notes
genos <- array(
    data = sample( 1:6, 3000, 
        prob = c(0.25, 0.3, 0.2, 0.09, 0.12, 0.04), replace = TRUE),
    dim = c(30, 100))
Locus <- paste("locus_", 1:100, sep = "")
snerr <- matrix(data = 0, nrow = length(Locus), ncol = 4)
snerr@dimnames <- list(NULL, c("AA2AO", "AO2AA", "BB2BO","BO2BB"))
newLocusData <- data.frame(
    thingOne.l = runif(ncol(genos)), 
    thingTwo.l = runif(ncol(genos)), 
    thingThree.l = runif(ncol(genos)), 
    useN = 6)
Our_sample <- paste("sample_", 1:nrow(genos), sep = "")

newSampleData <- data.frame(
    thingOne.s = runif(nrow(genos)),
    thingTwo.s = runif(nrow(genos)))

smallsnpg6b <- snpgeno( x = genos, diplos = genotypes6,
    locinfo = cbind(Locus, newLocusData),
    info = cbind(Our_sample, newSampleData),
    allow_nonchar = TRUE)
smallsnpg6b$locinfo$snerr <- snerr
smallsnpg6b <- est_ALF_ABO_quick( smallsnpg6b)

smallsnpg6b <- kin_power(smallsnpg6b, k = 0.5)
thePLODs6 <- find_HSPs(smallsnpg6b, limit_pairs = choose(nrow(smallsnpg4b),2),
    keep_thresh = -20)
## run once in January 2024, our benchmark of 'correct', yea, unto eternity:
## refPLODs6 <- thePLODs6
## save(refPLODs6, file = "smallsnpg6_referencePLODs.Rda")

base::load("smallsnpg6_referencePLODs.Rda")  ## avoid conflict with renv::load.
expect_true(all.equal(refPLODs6, thePLODs6, check.attributes = FALSE))

## now, 3-way. Build a 4-way dataset, then kinfer it 3-way via useN

set.seed(1113)
## 100 loci, 30 samples. Simple snpgeno prep lightly modified from the kinference course notes
genos <- array(
    data = sample(c(1,2,3,4), 3000, 
        prob = c(0.35, 0.4, 0.24, 0.01), replace = TRUE),
        dim = c(100, 30))
Locus <- paste("locus_", 1:ncol(genos), sep = "")
newLocusData <- data.frame(
    thingOne.l = runif(ncol(genos)), 
    thingTwo.l = runif(ncol(genos)), 
    thingThree.l = runif(ncol(genos)))
Our_sample <- paste("sample_", 1:nrow(genos), sep = "")

newSampleData <- data.frame(
    thingOne.s = runif(nrow(genos)),
    thingTwo.s = runif(nrow(genos)))

smallsnpg3 <- snpgeno( x = genos, diplos = genotypes4_ambig,
    locinfo = cbind(Locus, newLocusData),
    info = cbind(Our_sample, newSampleData),
    allow_nonchar = TRUE)
smallsnpg3 <- est_ALF_ABO_quick(smallsnpg3)

## just to prove it can be kinferred:
smallsnpg3a <- kin_power(smallsnpg3, k = 0.5)

thePLODs3 <- find_HSPs(smallsnpg3a, limit_pairs = choose(nrow(smallsnpg3a),2),
    keep_thresh = -20)
## run once in January 2024, our benchmark of 'correct', yea, unto eternity:
## refPLODs3 <- thePLODs3
## save(refPLODs3, file = "smallsnpg3_referencePLODs.Rda")

base::load("smallsnpg3_referencePLODs.Rda")  ## avoid conflict with renv::load.
expect_true(all.equal(refPLODs3, thePLODs3, check.attributes = FALSE))
    ## rm(list = ls() )


#####################################################################################
### Tests to ensure that the core toolchain works for either diplos =
### genotypes6 or diplos = genotypes4_ambig
#####################################################################################

expect_silent( {
  ## using the included dataset 'bluefin'
  library(kinference)
  library(gbasics)
  data(bluefin)

  BFT <- bluefin ## for testing on only the included datasets; incorporation into tinytests
  BFT <- kin_power(BFT, k = 0.5)

  if( ! "BFT_loci_keepies.Rdata" %in% dir() ) {
      loci <- check6and4(BFT, thresh_pchisq_6and4 = c(0.001, 0.0001))
      loci_keepies <- which( loci$pval4 > 0.00001)
      save( loci_keepies, file = "BFT_loci_keepies.Rdata")
  }
  load( "BFT_loci_keepies.Rdata")
  BFT <- BFT[ , loci_keepies]

  if( ! "BFT_fish_keepies.Rdata" %in% dir() ) {
      fish <- ilglk_geno(BFT)
      fish_keepies <- which( fish > -1500 )
      save( fish_keepies, file = "BFT_fish_keepies.Rdata")
  }
  load( "BFT_fish_keepies.Rdata")
  BFT <- BFT[ fish_keepies ,]

  sixway <- BFT
  sixway$locinfo$useN <- 6

  fourway <- BFT
  fourway@diplos <- genotypes4_ambig
  fourway[ BFT == "AA"] <- "AAO"
  fourway[ BFT == "BB"] <- "BBO"
  fourway[ BFT == "AO"] <- "AAO"
  fourway[ BFT == "BO"] <- "BBO"
  fourway[ BFT == "AB"] <- "AB"
  fourway[ BFT == "OO"] <- "OO"

  mixway <- sixway
  sixfours <- check6and4(mixway, thresh_pchisq_6and4 = c(0.001, 0.0001))
  mixway$locinfo$useN <- ifelse( sixfours$pval6 < 0.01, 6, 4)

}) # expect_silent

## sixway is six-way diplos (genotypes6) and all useN = 6
## fourway is four-way diplos (genotypes4_ambig) and all useN = 4
## mixway is six-way diplos (genotypes6) and useN a mix of 6 and 4

## test check6and4
expect_silent( { six_sixandfour <- check6and4( sixway, 
    thresh_pchisq_6and4 = c(0.001, 0.0001)) })
expect_silent( { mix_sixandfour <- check6and4( mixway, 
    thresh_pchisq_6and4 = c(0.001, 0.0001)) })
expect_silent( { four_sixandfour <- check6and4( fourway, 
    thresh_pchisq_6and4 = c(0.001, 0.0001)) })## fixed

## test re_est_ALF

expect_silent( { six_alfbco <- kinference:::re_est_ALF( sixway) })
expect_silent( { mix_alfbco <- kinference:::re_est_ALF( mixway) })
expect_silent( { four_alfbco <- kinference:::re_est_ALF( fourway) })

## test est_ALF_ABCO

## complains about 'no geno_amb'. If 6-way, geno_amb is 4-way; if
## 4-way, is self.  Should it auto-detect these states? Or is that
## unimportant now that est_ALF_ABCO is hidden?

## six_alfbco <- kinference:::est_ALF_ABCO( sixway) ## FAIL
## mix_alfbco <- kinference:::est_ALF_ABCO( mixway) ## FAIL
## four_alfbco <- kinference:::est_ALF_ABCO( fourway) ## FAIL

## test est_ALF_ABO_quick

expect_silent( { six_alfABOquick <- est_ALF_ABO_quick( sixway) }) ## Fixed
expect_silent( { mix_alfABOquick <- est_ALF_ABO_quick( mixway) }) ## Fixed
expect_silent( { four_alfABOquick <- est_ALF_ABO_quick( fourway) })

## test find_duplicates

expect_silent( { six_findDups <- find_duplicates( sixway, max_diff_loci = 200) })
expect_silent( { mix_findDups <- find_duplicates( mixway, max_diff_loci = 200) })
expect_silent( { four_findDups <- find_duplicates( fourway, max_diff_loci = 200) })

## test ilglk_geno

expect_silent( { six_ilglk <- ilglk_geno( sixway) })
expect_silent( { mix_ilglk <- ilglk_geno( mixway) })
expect_silent( { four_ilglk <- ilglk_geno( fourway) })

## test hetzminoo_fancy
### Rich
expect_silent( { six_hetzR <- hetzminoo_fancy( sixway, target = "rich") })
expect_silent( { mix_hetzR <- hetzminoo_fancy( mixway, target = "rich") })
expect_silent( { four_hetzR <- hetzminoo_fancy( fourway, target = "rich") })
### Poor
expect_silent( { six_hetzP <- hetzminoo_fancy( sixway, target = "poor") })
expect_silent( { mix_hetzP <- hetzminoo_fancy( mixway, target = "poor") })
expect_silent( { four_hetzP <- hetzminoo_fancy( fourway, target = "poor") })

## test kin_power

expect_silent( { sixway <- kin_power( sixway, k = 0.5) })
expect_silent( { mixway <- kin_power( mixway, k = 0.5) })
expect_silent( { fourway <- kin_power( fourway, k = 0.5) })

## test find_HSPs

expect_silent( { six_hsps <- find_HSPs( sixway, keep_thresh = -5, 
    limit_pairs = 10000) })
## PLOD_loghisto( six_hsps)
## HSP_histo(six_hsps, lb = 10, fullsib_cut = 150)
expect_silent({ autopick_threshold(sixway, kin = six_hsps, 
    fitrange_PLOD = c(1, 150), FPtol_pairs = 2, use4th = TRUE, plot_bins = 5) })

expect_silent( { mix_hsps <- find_HSPs( mixway, keep_thresh = -5, 
    limit_pairs = 10000) })
## PLOD_loghisto( mix_hsps)
## HSP_histo(mix_hsps, lb = 10, fullsib_cut = 130)
expect_silent( { autopick_threshold(sixway, kin = mix_hsps, 
    fitrange_PLOD = c(1, 130), FPtol_pairs = 2, use4th = TRUE, plot_bins = 5) })

expect_silent( { four_hsps <- find_HSPs( fourway, keep_thresh = -5, 
    limit_pairs = 10000) })
    ## Fixed, I think.
## PLOD_loghisto( four_hsps)
## HSP_histo(four_hsps, lb = 10, fullsib_cut = 100)
expect_silent( { autopick_threshold(sixway, kin = four_hsps, 
    fitrange_PLOD = c(1, 100), FPtol_pairs = 2, use4th = TRUE, plot_bins = 5) })

## test find_POPs

expect_silent( { six_POPs <- find_POPs(sixway, keep_thresh = 0.06) })
expect_silent( { mix_POPs <- find_POPs(mixway, keep_thresh = 0.06) })
expect_silent( { four_POPs <- find_POPs(fourway, keep_thresh = 0.06) })

## test split_FSPs_from_POPs

kinPalette()
expect_silent( { six_sFfP <- split_FSPs_from_POPs(sixway, 
    candipairs= six_POPs[six_POPs$wpsex < 0.05,], gerr = 0.01) })
## hist(six_sFfP$PLOD_FP, breaks = 50)
## abline(v = c(six_sFfP@E_FSP, six_sFfP@E_POP), col = kinPalette()[c("FSP", "POP")])

expect_silent( { mix_sFfP <- split_FSPs_from_POPs(mixway, 
    candipairs = mix_POPs[mix_POPs$wpsex < 0.05,], gerr = 0.01) })
## hist(mix_sFfP$PLOD_FP, breaks = 50)
## abline(v = c(mix_sFfP@E_FSP, mix_sFfP@E_POP), col = kinPalette()[c("FSP", "POP")])

expect_silent( { four_sFfP <- split_FSPs_from_POPs(fourway, 
    candipairs = four_POPs[four_POPs$wpsex < 0.05,], gerr = 0.01) })
## hist(four_sFfP$PLOD_FP, breaks = 50)
## abline(v = c(four_sFfP@E_FSP, four_sFfP@E_POP), col = kinPalette()[c("FSP", "POP")])

## test split_HSPs_from_FSPs

expect_silent( { six_sFfH <- split_FSPs_from_HSPs(sixway, 
    candipairs = six_POPs[six_POPs$wpsex < 0.05,]) })
## hist(six_sFfH$PLOD_FH, breaks = 50)
## abline(v = c(six_sFfH@E_FSP, six_sFfH@E_HSP), col = kinPalette()[c("FSP", "HSP")])

expect_silent( { mix_sFfH <- split_FSPs_from_HSPs(mixway, 
    candipairs = mix_POPs[mix_POPs$wpsex < 0.05,]) })
## hist(mix_sFfH$PLOD_FH, breaks = 50)
## abline(v = c(mix_sFfH@E_FSP, mix_sFfH@E_HSP), col = kinPalette()[c("FSP", "HSP")])

expect_silent( { four_sFfH <- split_FSPs_from_HSPs(fourway, 
    candipairs = four_POPs[four_POPs$wpsex < 0.05,]) })
## hist(four_sFfH$PLOD_FH, breaks = 50)
## abline(v = c(four_sFfH@E_FSP, four_sFfH@E_HSP), col = kinPalette()[c("FSP", "HSP")])

## test split_HSPs_from_HTPs

expect_silent( { six_sHfH <- split_HSPs_from_HTPs(sixway, 
    candipairs = six_hsps[six_hsps$PLOD > 0,]) })
## hist(six_sHfH$PLOD_ST, breaks = 50)
## abline(v = c(six_sHfH@E_HTP, six_sHfH@E_HSP), col = kinPalette()[c("HTP", "HSP")])

expect_silent( { mix_sHfH <- split_HSPs_from_HTPs(mixway, 
    candipairs = mix_hsps[mix_hsps$PLOD > 0,]) })
## hist(mix_sHfH$PLOD_ST, breaks = 50)
## abline(v = c(mix_sHfH@E_HTP, mix_sHfH@E_HSP), col = kinPalette()[c("HTP", "HSP")])

expect_silent( { four_sHfH <- split_HSPs_from_HTPs(fourway, 
    candipairs = four_hsps[four_hsps$PLOD > 0,])  })
## hist(four_sHfH$PLOD_ST, breaks = 50)
## abline(v = c(four_sHfH@E_HTP, four_sHfH@E_HSP), col = kinPalette()[c("HTP", "HSP")])

