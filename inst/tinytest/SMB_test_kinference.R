## These are the unit tests for package 'kinference'
## To run (on Shane's machine), use tinytest::test_all("~/R/localPackages/kinference/kinference").
## On other machines, edit that filepath to the source dir for package kinference.

library(atease)
library(kinference)
library(gbasics) ## shouldn't be necessary to do this manually
library(mvbutils)

shanesComp <- Sys.info()["nodename"] == "pacific-hf"  ## because !Shane doesn't have the CSIRO datasets

## absolute basic mechanics

    expect_true("kinference" %in% sessionInfo()$otherPkgs$kinference)

    info <- data.frame(Our_sample = c("sample1", "sample2"), location = c("loc1", "loc2") )
locinfo <- data.frame(Locus = c("L1", "L2"), n_alleles = c(2,2))
genos <- matrix(data = c(rep("AAO", 4)), nrow = 2, ncol = 2)
    minisnpg <- snpgeno(genos, diplos = c("AAO", "AB", "BBO", "OO"), info = info, locinfo = locinfo)
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
chargenos6 <- matrix( data = c( "AA", "AO", "BB", "AB", "OO", "BO"), nrow = 2, ncol = 3)
intgenos6 <- matrix( data = c(1,3,4,2,6,5) , nrow = 2, ncol = 3)
rawgenos6 <- as.raw( matrix( data = c(1,3,4,2,6,5) , nrow = 2, ncol = 3))

info6 <- data.frame(Our_sample = c("sample1", "sample2"), location = c("loc1", "loc2") )
locinfo6 <- data.frame(Locus = c("L1", "L2", "L3"), n_alleles = c(2,2,2))

### 4-way parts
chargenos4 <- matrix( data = c( "OO", "AAO", "BBO", "AB"), nrow = 2, ncol = 2)
intgenos4 <- matrix( data = c(1,3,4,2) , nrow = 2, ncol = 2)
rawgenos4 <- as.raw( matrix( data = c(1,3,4,2) , nrow = 2, ncol = 2))

info4 <- data.frame(Our_sample = c("sample1", "sample2"), location = c("loc1", "loc2") )
locinfo4 <- data.frame(Locus = c("L1", "L2"), n_alleles = c(2,2))

### 6-way diplos
#### character genos
expect_silent( {
snpg6_char <- snpgeno( x = chargenos6, diplos = genotypes6, info = info6, locinfo = locinfo6)
})
#### integer genos
expect_silent( {
snpg6_int <- snpgeno( x = intgenos6, diplos = genotypes6, info = info6, locinfo = locinfo6,
                     allow_nonchar = TRUE)
})
#### raw genos
expect_silent( {
snpg6_raw <- snpgeno( x = rawgenos6, diplos = genotypes6, info = info6, locinfo = locinfo6,
                     allow_nonchar = TRUE)
})
##### check that char, int, and raw inputs all give the same genotypes
expect_true( { all( snpg6_char == snpg6_int) })
expect_true( { all( snpg6_char == snpg6_raw) })

#### empty genos
expect_silent( {
snpg6_nogeno <- snpgeno( x = NULL, diplos = genotypes6, info = info6, locinfo = locinfo6)
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
snpg4_char <- snpgeno( x = chargenos4, diplos = genotypes4_ambig, info = info4, locinfo = locinfo4)
})
#### integer genos
expect_silent( {
snpg4_int <- snpgeno( x = intgenos4, diplos = genotypes4_ambig, info = info4, locinfo = locinfo4,
                     allow_nonchar = TRUE)
})
#### raw genos
expect_silent( {
snpg4_raw <- snpgeno( x = rawgenos4, diplos = genotypes4_ambig, info = info4, locinfo = locinfo4,
                     allow_nonchar = TRUE)
})
##### check that char, int, and raw inputs all give the same genotypes
expect_true( { all( snpg4_char == snpg4_int) })
expect_true( { all( snpg4_char == snpg4_raw) })

#### empty genos
expect_silent( {
snpg4_nogeno <- snpgeno( x = NULL, diplos = genotypes4_ambig, info = info4, locinfo = locinfo4)
})
#### empty info
expect_silent( {
snpg4_noinfo <- snpgeno( x = chargenos4, diplos = genotypes4_ambig, locinfo = locinfo4)
})
#### empty locinfo
expect_silent( {
snpg4_nolocinfo <- snpgeno( x = chargenos4, diplos = genotypes4_ambig, info = info4)
})
#### empty genos and info
expect_silent( {
snpg4_nogeno_noinfo <- snpgeno( x = NULL, diplos = genotypes4_ambig,
                               locinfo = locinfo4, n_samples = 2)
})
#### empty genos and locinfo
expect_silent( {
snpg4_nogeno_nolocinfo <- snpgeno( x = NULL, diplos = genotypes4_ambig, info = info4,
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
genos <- array(data = sample(c(1,2,4,6), 3000, prob = c(0.35, 0.4, 0.24, 0.01), replace = TRUE),
               dim = c(30, 100))
Locus <- paste("locus_", 1:100, sep = "")
newLocusData <- data.frame(thingOne.l = runif(ncol(genos)), thingTwo.l = runif(ncol(genos)),
                           thingThree.l = runif(ncol(genos)))
locinfo <- cbind( Locus, newLocusData)

Our_sample <- paste("sample_", 1:30, sep = "")
newSampleData <- data.frame(thingOne.s = runif(nrow(genos)),thingTwo.s = runif(nrow(genos)))
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
smallsnpg4b <- smallsnpg4[ ! smallsnpg4$info$Our_sample %in% c("sample_10", "sample_13") ,]
expect_true( dim(smallsnpg4b)[1] == dim(smallsnpg4b$info)[1] )
#### subset samples by number
smallsnpg4c <- smallsnpg4[ 1:10 ,]
expect_true( dim(smallsnpg4c)[1] == dim(smallsnpg4c$info)[1] )
#### subset loci by lookup
smallsnpg4d <- smallsnpg4[, ! smallsnpg4$locinfo$Locus %in% c("locus_10", "locus_13")]
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

thePLODs4 <- find_HSPs(smallsnpg4b, limit_pairs = choose(nrow(smallsnpg4b),2), keep_thresh = -20)
## test w/o independent prepare_PLOD_SPA:
expect_silent({ find_HSPs(smallsnpg4a, limit_pairs = choose(nrow(smallsnpg4b),2), keep_thresh = -20) } )

## run once in January 2024, our benchmark of 'correct', yea, unto eternity:
## refPLODs4 <- thePLODs4
## save(refPLODs4, file = "smallsnpg4_referencePLODs.Rda")

base::load("smallsnpg4_referencePLODs.Rda")  ## base:: to avoid a conflict with renv::load.
expect_true(all.equal(refPLODs4, thePLODs4, check.attributes = FALSE))

## test auto-generation of snerr by kin_power in cases where it's absent
expect_true( length( dim(smallsnpg4b$locinfo$snerr)) > 0 ) ## this is basically just
## 'exists( "smallsnpg4b$locinfo$snerr")', but I can't atease within a quoted string

## now, 6-way.
set.seed(1112)
## 100 loci, 30 samples. Simple snpgeno prep lightly modified from the kinference course notes
genos <- array(data = sample(c(1:6), 3000, prob = c(0.25, 0.3, 0.2, 0.09, 0.12, 0.04),
                             replace = TRUE),dim = c(30, 100))
Locus <- paste("locus_", 1:100, sep = "")
snerr <- matrix(data = 0, nrow = length(Locus), ncol = 4)
snerr@dimnames <- list(NULL, c("AA2AO", "AO2AA", "BB2BO","BO2BB"))
newLocusData <- data.frame(thingOne.l = runif(ncol(genos)), thingTwo.l = runif(ncol(genos)), thingThree.l = runif(ncol(genos)), useN = 6)
Our_sample <- paste("sample_", 1:nrow(genos), sep = "")

newSampleData <- data.frame(thingOne.s = runif(nrow(genos)),thingTwo.s = runif(nrow(genos)))

smallsnpg6b <- snpgeno( x = genos, diplos = genotypes6,
                       locinfo = cbind(Locus, newLocusData),
                       info = cbind(Our_sample, newSampleData),
                       allow_nonchar = TRUE)
smallsnpg6b$locinfo$snerr <- snerr
smallsnpg6b <- est_ALF_ABO_quick( smallsnpg6b)

smallsnpg6b <- kin_power(smallsnpg6b, k = 0.5)
thePLODs6 <- find_HSPs(smallsnpg6b, limit_pairs = choose(nrow(smallsnpg4b),2), keep_thresh = -20)
## run once in January 2024, our benchmark of 'correct', yea, unto eternity:
## refPLODs6 <- thePLODs6
## save(refPLODs6, file = "smallsnpg6_referencePLODs.Rda")

    base::load("smallsnpg6_referencePLODs.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refPLODs6, thePLODs6, check.attributes = FALSE))

## now, 3-way. Build a 4-way dataset, then kinfer it 3-way via useN

set.seed(1113)
## 100 loci, 30 samples. Simple snpgeno prep lightly modified from the kinference course notes
genos <- array(data = sample(c(1,2,3,4), 3000, prob = c(0.35, 0.4, 0.24, 0.01), replace = TRUE),dim = c(100, 30))
Locus <- paste("locus_", 1:ncol(genos), sep = "")
newLocusData <- data.frame(thingOne.l = runif(ncol(genos)), thingTwo.l = runif(ncol(genos)), thingThree.l = runif(ncol(genos)))
Our_sample <- paste("sample_", 1:nrow(genos), sep = "")

newSampleData <- data.frame(thingOne.s = runif(nrow(genos)),thingTwo.s = runif(nrow(genos)))

smallsnpg3 <- snpgeno( x = genos, diplos = genotypes4_ambig,
                       locinfo = cbind(Locus, newLocusData),
                       info = cbind(Our_sample, newSampleData),
                      allow_nonchar = TRUE)
smallsnpg3 <- est_ALF_ABO_quick(smallsnpg3)

## just to prove it can be kinferred:
smallsnpg3a <- kin_power(smallsnpg3, k = 0.5)

thePLODs3 <- find_HSPs(smallsnpg3a, limit_pairs = choose(nrow(smallsnpg3a),2), keep_thresh = -20)
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

    })

## sixway is six-way diplos (genotypes6) and all useN = 6
## fourway is four-way diplos (genotypes4_ambig) and all useN = 4
## mixway is six-way diplos (genotypes6) and useN a mix of 6 and 4

## test check6and4
expect_silent( { six_sixandfour <- check6and4( sixway, thresh_pchisq_6and4 = c(0.001, 0.0001)) })
expect_silent( { mix_sixandfour <- check6and4( mixway, thresh_pchisq_6and4 = c(0.001, 0.0001)) })
expect_silent( { four_sixandfour <- check6and4( fourway, thresh_pchisq_6and4 = c(0.001, 0.0001)) })## fixed

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

expect_silent( { six_hsps <- find_HSPs( sixway, keep_thresh = -5, limit_pairs = 10000) })
## PLOD_loghisto( six_hsps)
## HSP_histo(six_hsps, lb = 10, fullsib_cut = 150)
expect_silent( { autopick_threshold(sixway, kin = six_hsps, fitrange_PLOD = c(1, 150), FPtol_pairs = 2, use4th = TRUE, plot_bins = 5) })

expect_silent( { mix_hsps <- find_HSPs( mixway, keep_thresh = -5, limit_pairs = 10000) })
## PLOD_loghisto( mix_hsps)
## HSP_histo(mix_hsps, lb = 10, fullsib_cut = 130)
expect_silent( { autopick_threshold(sixway, kin = mix_hsps, fitrange_PLOD = c(1, 130), FPtol_pairs = 2, use4th = TRUE, plot_bins = 5) })

    expect_silent( { four_hsps <- find_HSPs( fourway, keep_thresh = -5, limit_pairs = 10000) })
    ## Fixed, I think.
## PLOD_loghisto( four_hsps)
## HSP_histo(four_hsps, lb = 10, fullsib_cut = 100)
expect_silent( { autopick_threshold(sixway, kin = four_hsps, fitrange_PLOD = c(1, 100), FPtol_pairs = 2, use4th = TRUE, plot_bins = 5) })

## test find_POPs

expect_silent( { six_POPs <- find_POPs(sixway, keep_thresh = 0.06) })
expect_silent( { mix_POPs <- find_POPs(mixway, keep_thresh = 0.06) })
expect_silent( { four_POPs <- find_POPs(fourway, keep_thresh = 0.06) })

## test split_FSPs_from_POPs

kinPalette()
expect_silent( { six_sFfP <- split_FSPs_from_POPs(sixway, candipairs = six_POPs[six_POPs$wpsex < 0.05,], gerr = 0.01) })
## hist(six_sFfP$PLOD_FP, breaks = 50)
## abline(v = c(six_sFfP@E_FSP, six_sFfP@E_POP), col = kinPalette()[c("FSP", "POP")])

expect_silent( { mix_sFfP <- split_FSPs_from_POPs(mixway, candipairs = mix_POPs[mix_POPs$wpsex < 0.05,], gerr = 0.01) })
## hist(mix_sFfP$PLOD_FP, breaks = 50)
## abline(v = c(mix_sFfP@E_FSP, mix_sFfP@E_POP), col = kinPalette()[c("FSP", "POP")])

expect_silent( { four_sFfP <- split_FSPs_from_POPs(fourway, candipairs = four_POPs[four_POPs$wpsex < 0.05,], gerr = 0.01) })
## hist(four_sFfP$PLOD_FP, breaks = 50)
## abline(v = c(four_sFfP@E_FSP, four_sFfP@E_POP), col = kinPalette()[c("FSP", "POP")])

## test split_HSPs_from_FSPs

expect_silent( { six_sFfH <- split_FSPs_from_HSPs(sixway, candipairs = six_POPs[six_POPs$wpsex < 0.05,]) })
## hist(six_sFfH$PLOD_FH, breaks = 50)
## abline(v = c(six_sFfH@E_FSP, six_sFfH@E_HSP), col = kinPalette()[c("FSP", "HSP")])

expect_silent( { mix_sFfH <- split_FSPs_from_HSPs(mixway, candipairs = mix_POPs[mix_POPs$wpsex < 0.05,]) })
## hist(mix_sFfH$PLOD_FH, breaks = 50)
## abline(v = c(mix_sFfH@E_FSP, mix_sFfH@E_HSP), col = kinPalette()[c("FSP", "HSP")])

expect_silent( { four_sFfH <- split_FSPs_from_HSPs(fourway, candipairs = four_POPs[four_POPs$wpsex < 0.05,]) })
## hist(four_sFfH$PLOD_FH, breaks = 50)
## abline(v = c(four_sFfH@E_FSP, four_sFfH@E_HSP), col = kinPalette()[c("FSP", "HSP")])

## test split_HSPs_from_HTPs

expect_silent( { six_sHfH <- split_HSPs_from_HTPs(sixway, candipairs = six_hsps[six_hsps$PLOD > 0,]) })
## hist(six_sHfH$PLOD_ST, breaks = 50)
## abline(v = c(six_sHfH@E_HTP, six_sHfH@E_HSP), col = kinPalette()[c("HTP", "HSP")])

expect_silent( { mix_sHfH <- split_HSPs_from_HTPs(mixway, candipairs = mix_hsps[mix_hsps$PLOD > 0,]) })
## hist(mix_sHfH$PLOD_ST, breaks = 50)
## abline(v = c(mix_sHfH@E_HTP, mix_sHfH@E_HSP), col = kinPalette()[c("HTP", "HSP")])

expect_silent( { four_sHfH <- split_HSPs_from_HTPs(fourway, candipairs = four_hsps[four_hsps$PLOD > 0,])  })
## hist(four_sHfH$PLOD_ST, breaks = 50)
## abline(v = c(four_sHfH@E_HTP, four_sHfH@E_HSP), col = kinPalette()[c("HTP", "HSP")])


#############################################################################################
#### Individual Species output-consistency tests - Shane's computer only  ###################
#############################################################################################

shanesComp <- Sys.info()["nodename"] == "pacific-hf" ## because !Shane doesn't have the CSIRO datasets

##############################################################################################
## SBT testing  ##############################################################################
##############################################################################################

if(shanesComp) {

    base::load("~/Data/geno2019.rda") ## only on shane's computer -- might be worth copying this across, but risky IP-wise
    theSBT6and4s <- check6and4(geno2019, thresh_pchisq_6and4 = c(0.0001, 0.00001))
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSBT6and4s <- theSBT6and4s
    ## save(refSBT6and4s, file = "SBT_reference6and4s.Rda")

    base::load("SBT_reference6and4s.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSBT6and4s, theSBT6and4s, check.attributes = FALSE))

    geno2019$locinfo$useN <- ifelse( geno2019$locinfo$use6, 6, 4)
    theSBTilglks <- ilglk_geno(geno2019)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSBTilglks <- theSBTilglks
    ## save(refSBTilglks, file = "SBT_referenceilglks.Rda")

    base::load("SBT_referenceilglks.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSBTilglks, theSBTilglks, check.attributes = FALSE))

    geno2019 <- geno2019[theSBTilglks > -1680 & theSBTilglks < -1280 ,]

    theSBThetzpoors <- hetzminoo_fancy(geno2019, 'poor')
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSBThetzpoors <- theSBThetzpoors
    ## save(refSBThetzpoors, file = "SBT_referencehetzpoors.Rda")

    base::load("SBT_referencehetzpoors.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSBThetzpoors, theSBThetzpoors, check.attributes = FALSE))

    geno2019 <- geno2019[theSBThetzpoors > 0.18 & theSBThetzpoors < 0.27 ,]

    theSBThetzriches <- hetzminoo_fancy(geno2019, 'rich')
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSBThetzriches <- theSBThetzriches
    ## save(refSBThetzriches, file = "SBT_referencehetzriches.Rda")

    base::load("SBT_referencehetzriches.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSBThetzriches, theSBThetzriches, check.attributes = FALSE))

    geno2019 <- geno2019[theSBThetzriches > 0.19 & theSBThetzriches < 0.28 ,]

    geno2019 <- kin_power(geno2019, k = 0.5)
    geno2019 <- prepare_PLOD_SPA(geno2019)

    ## for the purposes of keeping the test quick-to-run, here we import a list of animals involved in high-PLOD pairs
    ## that aren't duplicates,
    ## and trim the dataset to just those animals. That way, everything runs _much_ quicker, and the test for
    ## matching kinship statistics is just the same amongst the pairs where it matters.
    ## Here's the code for generating that list, for posterity's sake:

    ## theSBTdups <- find_duplicates(geno2019, max_diff_loci = 50)
    ## geno2019 <- geno2019[-c(drop_dups_pairwise_equiv(theSBTdups[, 2:3])), ]
    ## theSBT_PLODs <- find_HSPs(geno2019, keep_thresh = 25)
    ## ijs <- unique(c(theSBT_PLODs$i, theSBT_PLODs$j))
    ## save(ijs, file = "SBT_highPLODs.Rda")

    ## Now, let's check the consistency of our results, only for those animals
    ## involved in a high PLOD pair
    base::load("SBT_highPLODs.Rda")
    geno2019 <- geno2019[c(ijs),]

    theSBTdups <- find_duplicates(geno2019, max_diff_loci = 200)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSBTdups <- theSBTdups
    ## save(refSBTdups, file = "SBT_referencedups.Rda")

    geno2019 <- geno2019[-c(drop_dups_pairwise_equiv(theSBTdups[, 2:3])), ]

    theSBTPLODs <- find_HSPs(geno2019, keep_thresh = 0)
    ## HSP_histo(theSBTPLODs, lb = 0, fullsib_cut = 25)
    ## PLOD_loghisto(theSBTPLODs)

    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSBTPLODs <- theSBTPLODs
    ## save(refSBTPLODs, file = "SBT_referencePLODs.Rda")

    base::load("SBT_referencePLODs.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSBTPLODs, theSBTPLODs, check.attributes = FALSE))

    theSBTPOPs <- find_POPs(geno2019, keep_thresh = 0.91)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSBTPOPs <- theSBTPOPs
    ## save(refSBTPOPs, file = "SBT_referencePOPs.Rda")

    ## theSBTlglkPOPs <- find_POPs_lglk(geno2019, gerr = 0.01, keep_thresh = 50)

    base::load("SBT_referencePOPs.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSBTPOPs, theSBTPOPs, check.attributes = FALSE))

    theSBTwtsame <- split_FSPs_from_POPs(geno2019, candipairs = refSBTPOPs, gerr = 0.01)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSBTwtsame <- theSBTwtsame
    ## save(refSBTwtsame, file = "SBT_referencewtsame.Rda")

    base::load("SBT_referencewtsame.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSBTwtsame, theSBTwtsame, check.attributes = FALSE))

    theSBTwpsex <- split_FSPs_from_HSPs(geno2019, candipairs = refSBTPOPs)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSBTwpsex <- theSBTwpsex
    ## save(refSBTwpsex, file = "SBT_referencewpsex.Rda")

    base::load("SBT_referencewpsex.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSBTwpsex, theSBTwpsex, check.attributes = FALSE))

    theSBTvarhtp <- var_PLOD_kin(geno2019$locinfo, emp_V_HSP = 900, n_meio = 3)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSBTvarhtp <- theSBTvarhtp
    ## save(refSBTvarhtp, file = "SBT_referencevarhtp.Rda")

    base::load("SBT_referencevarhtp.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSBTvarhtp, theSBTvarhtp, check.attributes = FALSE))

    theSBTPLODST <- split_HSPs_from_HTPs(geno2019, candipairs = theSBTPLODs)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSBTPLODST <- theSBTPLODST
    ## save(refSBTPLODST, file = "SBT_referencePLODST.Rda")

    base::load("SBT_referencePLODST.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSBTPLODST, theSBTPLODST, check.attributes = FALSE))
}



##############################################################################################
## GGar testing  #############################################################################
##############################################################################################

if(shanesComp) {

    base::load("~/Data/Glyphis garricki/ggarnear2.rda") ## only on shane's computer -- might be worth copying this across, but risky IP-wise
    theGGar6and4s <- check6and4(ggarnear2, thresh_pchisq_6and4 = c(0.0001, 0.00001))
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGGar6and4s <- theGGar6and4s
    ## save(refGGar6and4s, file = "GGar_reference6and4s.Rda")

    base::load("GGar_reference6and4s.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGGar6and4s, theGGar6and4s, check.attributes = FALSE))

    ggarnear2$locinfo$useN <- ifelse( ggarnear2$locinfo$use6 == 6, 6, 4)

    theGGarilglks <- ilglk_geno(ggarnear2)
    ## ## run once in January 2024, our benchmark of 'correct', yea, unto eternity:
    ## refGGarilglks <- theGGarilglks
    ## save(refGGarilglks, file = "GGar_referenceilglks.Rda")

    base::load("GGar_referenceilglks.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGGarilglks, theGGarilglks, check.attributes = FALSE))

    ggarnear2 <- ggarnear2[theGGarilglks > -1100 & theGGarilglks < -1000 ,]

    theGGarhetzpoors <- hetzminoo_fancy(ggarnear2, 'poor')
    ## ## run once in January 2024, our benchmark of 'correct', yea, unto eternity:
    ## refGGarhetzpoors <- theGGarhetzpoors
    ## save(refGGarhetzpoors, file = "GGar_referencehetzpoors.Rda")

    base::load("GGar_referencehetzpoors.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGGarhetzpoors, theGGarhetzpoors, check.attributes = FALSE))

    ggarnear2 <- ggarnear2[theGGarhetzpoors > 0.30 & theGGarhetzpoors < 0.36 ,]

    theGGarhetzriches <- hetzminoo_fancy(ggarnear2, 'rich')
    ## run once in January 2024, our benchmark of 'correct', yea, unto eternity:
    ## refGGarhetzriches <- theGGarhetzriches
    ## save(refGGarhetzriches, file = "GGar_referencehetzriches.Rda")

    base::load("GGar_referencehetzriches.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGGarhetzriches, theGGarhetzriches, check.attributes = FALSE))

    ggarnear2 <- ggarnear2[theGGarhetzriches > 0.32 & theGGarhetzriches < 0.37 ,]

    ggarnear2 <- kin_power(ggarnear2, k = 0.5)
    ggarnear2 <- prepare_PLOD_SPA(ggarnear2)

    ## for the purposes of keeping the test quick-to-run, here we import a list of animals involved in high-PLOD pairs
    ## and trim the dataset to just those animals. That way, everything runs _much_ quicker, and the test for
    ## matching kinship statistics is just the same amongst the pairs where it matters.
    ## Here's the code for generating that list, for posterity's sake:

    ## theGGar_PLODs <- find_HSPs(ggarnear2, keep_thresh = 25)
    ## ijs <- unique(c(theGGar_PLODs$i, theGGar_PLODs$j))
    ## save(ijs, file = "GGar_highPLODs.Rda")

    ## Now, let's check the consistency of our results, only for those animals
    ## involved in a high PLOD pair
    base::load("GGar_highPLODs.Rda")
    ## ggarnear2 <- ggarnear2[ijs,]

    theGGardups <- find_duplicates(ggarnear2, max_diff_loci = 200)
    ## run once in January 2024, our benchmark of 'correct', yea, unto eternity:
    ## refGGardups <- theGGardups
    ## save(refGGardups, file = "GGar_referencedups.Rda")

    if(nrow(theGGardups) > 0) {
        ggarnear2 <- ggarnear2[-c(drop_dups_pairwise_equiv(theGGardups[, 2:3])), ]
    }

    theGGarPLODs <- find_HSPs(ggarnear2, keep_thresh = 0)
    ## HSP_histo(theGGarPLODs, lb = 0, fullsib_cut = 25)
    ## PLOD_loghisto(theGGarPLODs)

    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGGarPLODs <- theGGarPLODs
    ## save(refGGarPLODs, file = "GGar_referencePLODs.Rda")

    base::load("GGar_referencePLODs.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGGarPLODs, theGGarPLODs, check.attributes = FALSE))

    theGGarPOPs <- find_POPs(ggarnear2, keep_thresh = 0.91)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGGarPOPs <- theGGarPOPs
    ## save(refGGarPOPs, file = "GGar_referencePOPs.Rda")

    base::load("GGar_referencePOPs.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGGarPOPs, theGGarPOPs, check.attributes = FALSE))

    theGGarwtsame <- split_FSPs_from_POPs(ggarnear2, candipairs = theGGarPOPs, gerr = 0.01)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGGarwtsame <- theGGarwtsame
    ## save(refGGarwtsame, file = "GGar_referencewtsame.Rda")

    base::load("GGar_referencewtsame.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGGarwtsame, theGGarwtsame, check.attributes = FALSE))

    theGGarwpsex <- split_FSPs_from_HSPs(ggarnear2, candipairs = theGGarPOPs)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGGarwpsex <- theGGarwpsex
    ## save(refGGarwpsex, file = "GGar_referencewpsex.Rda")

    base::load("GGar_referencewpsex.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGGarwpsex, theGGarwpsex, check.attributes = FALSE))

    theGGarvarhtp <- var_PLOD_kin(ggarnear2$locinfo, emp_V_HSP = 900, n_meio = 3)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGGarvarhtp <- theGGarvarhtp
    ## save(refGGarvarhtp, file = "GGar_referencevarhtp.Rda")

    base::load("GGar_referencevarhtp.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGGarvarhtp, theGGarvarhtp, check.attributes = FALSE))

    theGGarPLODST <- split_HSPs_from_HTPs(ggarnear2, candipairs = theGGarPLODs)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGGarPLODST <- theGGarPLODST
    ## save(refGGarPLODST, file = "GGar_referencePLODST.Rda")

    base::load("GGar_referencePLODST.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGGarPLODST, theGGarPLODST, check.attributes = FALSE))
}

##############################################################################################
## GNS testing  ##############################################################################
##############################################################################################

if(shanesComp) {

    base::load("~/Data/Grey Nurse Sharks/Kin_stuff.rda") ## only on shane's computer -- might be worth copying this across, but risky IP-wise
    theGNS6and4s <- check6and4(gns14, thresh_pchisq_6and4 = c(0.0001, 0.00001))
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGNS6and4s <- theGNS6and4s
    ## save(refGNS6and4s, file = "GNS_reference6and4s.Rda")

    base::load("GNS_reference6and4s.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGNS6and4s, theGNS6and4s, check.attributes = FALSE))

    gns14$locinfo$useN <- ifelse( gns14$locinfo$use6, 6, 4)

    theGNSilglks <- ilglk_geno(gns14)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGNSilglks <- theGNSilglks
    ## save(refGNSilglks, file = "GNS_referenceilglks.Rda")

    base::load("GNS_referenceilglks.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGNSilglks, theGNSilglks, check.attributes = FALSE))

    gns14 <- gns14[theGNSilglks > -1740 & theGNSilglks < -1560 ,]

    theGNShetzpoors <- hetzminoo_fancy(gns14, 'poor')
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGNShetzpoors <- theGNShetzpoors
    ## save(refGNShetzpoors, file = "GNS_referencehetzpoors.Rda")

    base::load("GNS_referencehetzpoors.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGNShetzpoors, theGNShetzpoors, check.attributes = FALSE))

    gns14 <- gns14[theGNShetzpoors > 0.24 & theGNShetzpoors < 0.3 ,]

    theGNShetzriches <- hetzminoo_fancy(gns14, 'rich')
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGNShetzriches <- theGNShetzriches
    ## save(refGNShetzriches, file = "GNS_referencehetzriches.Rda")

    base::load("GNS_referencehetzriches.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGNShetzriches, theGNShetzriches, check.attributes = FALSE))

    gns14 <- gns14[theGNShetzriches > 0.28 & theGNShetzriches < 0.35 ,]

    gns14 <- kin_power(gns14, k = 0.5)
    gns14 <- prepare_PLOD_SPA(gns14)

    ## for the purposes of keeping the test quick-to-run, here we import a list of animals involved in high-PLOD pairs
    ## and trim the dataset to just those animals. That way, everything runs _much_ quicker, and the test for
    ## matching kinship statistics is just the same amongst the pairs where it matters.
    ## Here's the code for generating that list, for posterity's sake:

    ## theGNS_PLODs <- find_HSPs(gns14, keep_thresh = 25)
    ## ijs <- unique(c(theGNS_PLODs$i, theGNS_PLODs$j))
    ## save(ijs, file = "GNS_highPLODs.Rda")

    ## Now, let's check the consistency of our results, only for those animals
    ## involved in a high PLOD pair
    base::load("GNS_highPLODs.Rda")
    gns14 <- gns14[ijs,]

    theGNSdups <- find_duplicates(gns14, max_diff_loci = 200)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGNSdups <- theGNSdups
    ## save(refGNSdups, file = "GNS_referencedups.Rda")

    if(nrow(theGNSdups) > 0) {
        gns14 <- gns14[-c(drop_dups_pairwise_equiv(theGNSdups[, 2:3])), ]
    }

    theGNSPLODs <- find_HSPs(gns14, keep_thresh = 0)
    ## HSP_histo(theGNSPLODs, lb = 0, fullsib_cut = 25)
    ## PLOD_loghisto(theGNSPLODs)

    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGNSPLODs <- theGNSPLODs
    ## save(refGNSPLODs, file = "GNS_referencePLODs.Rda")

    base::load("GNS_referencePLODs.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGNSPLODs, theGNSPLODs, check.attributes = FALSE))

    theGNSPOPs <- find_POPs(gns14, keep_thresh = 0.91)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGNSPOPs <- theGNSPOPs
    ## save(refGNSPOPs, file = "GNS_referencePOPs.Rda")

    base::load("GNS_referencePOPs.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGNSPOPs, theGNSPOPs, check.attributes = FALSE))

    theGNSwtsame <- split_FSPs_from_POPs(gns14, candipairs = theGNSPOPs, gerr = 0.01)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGNSwtsame <- theGNSwtsame
    ## save(refGNSwtsame, file = "GNS_referencewtsame.Rda")

    base::load("GNS_referencewtsame.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGNSwtsame, theGNSwtsame, check.attributes = FALSE))

    theGNSwpsex <- split_FSPs_from_HSPs(gns14, candipairs = theGNSPOPs)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGNSwpsex <- theGNSwpsex
    ## save(refGNSwpsex, file = "GNS_referencewpsex.Rda")

    base::load("GNS_referencewpsex.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGNSwpsex, theGNSwpsex, check.attributes = FALSE))

    theGNSvarhtp <- var_PLOD_kin(gns14$locinfo, emp_V_HSP = 900, n_meio = 3)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGNSvarhtp <- theGNSvarhtp
    ## save(refGNSvarhtp, file = "GNS_referencevarhtp.Rda")

    base::load("GNS_referencevarhtp.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGNSvarhtp, theGNSvarhtp, check.attributes = FALSE))

    theGNSPLODST <- split_HSPs_from_HTPs(gns14, candipairs = theGNSPLODs)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refGNSPLODST <- theGNSPLODST
    ## save(refGNSPLODST, file = "GNS_referencePLODST.Rda")

    base::load("GNS_referencePLODST.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refGNSPLODST, theGNSPLODST, check.attributes = FALSE))
}

##############################################################################################
## ABT testing  ##############################################################################
##############################################################################################

if(shanesComp) {

    library(plyr)
    base::load("~/Data/allFish.Rdata") ## only on shane's computer -- might be worth copying this across, but risky IP-wise
    strata <- read.csv("~/R/localPackages/kinference/kinference/inst/tinytest/GOML_FAKErecall_2800markers_8Kcluster_recallSNP_PG.csv")
    strata$TargetID <- strata$TARGET_ID; strata$TARGET_ID <- NULL
    allFish <- allFish[allFish$info$TargetID %in% strata$TargetID,]
    ## allFish$info <- merge(allFish$info, strata, sort = FALSE) ## introduces the misalignment bug
    allFish$info <- join(allFish$info, strata, type = "left", by = "TargetID")
    allFish$locinfo$pbonzer <- kinference:::re_est_ALF(allFish)$locinfo$pambig

    which_not <- which(allFish$info$STRATA %in% cq(ASH2016, ASH2017, CAN2016, CAN2017, GOL2016,
                                                      GOL2017, GOL2018, AzT, BLF2020, GoM, iTY2020, MED) )
    ## adultIDs <- allFish$info$TargetID[which_adults]
    ABT <- allFish[ -c(which_not),]; rm(allFish)
    ## technical replicates
    dups <- duplicated(ABT$info$Our_sample)
    ABT <- ABT[!dups,]

    theABT6and4s <- check6and4(ABT, thresh_pchisq_6and4 = c(0.0001, 0.00001))
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refABT6and4s <- theABT6and4s
    ## save(refABT6and4s, file = "ABT_reference6and4s.Rda")

    base::load("ABT_reference6and4s.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refABT6and4s, theABT6and4s, check.attributes = FALSE))

    ## At <- ABT[ABT$info$STRATA %in% cq(At6, At7, At8),]
    ## pew <- ABT[ABT$info$STRATA %in% cq(PEW2016, PEW2018, PEW2019),]
    ## can <- ABT[ABT$info$STRATA %in% cq(CAN2018, CAN2019),]
    ## MiA <- ABT[ABT$info$STRATA %in% cq(MiA),] ## the misalignment bug was noticed here: something in MiA breaks ilglk_geno,
    ## causing snpg4 to have no dim this appears to be linked to MiA$info having one fewer rows than there
    ## are samples in MiA.
    ## Ultimately this traced back to load_whopper reading in one blank line of all '?' genotypes and all-NA info,
    ## which caused the strata file to merge incorrectly (out-by-one) with allFish$info

    theABTilglks <- ilglk_geno(ABT, showPlot = FALSE)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refABTilglks <- theABTilglks
    ## save(refABTilglks, file = "ABT_referenceilglks.Rda")

    base::load("ABT_referenceilglks.Rda")  ## avoid conflict with renv::load.
    ## Does fixing the misalignment also fix the 'tiny troublemakers' bug?
    ## tiny_troublemakers <- is.nan(theABTilglks) | is.nan(refABTilglks)
    ## expect_true(sum(tiny_troublemakers) < 10)
    ## ABT <- ABT[! tiny_troublemakers ]
    ## theABTilglks <- theABTilglks[!tiny_troublemakers]
    ## refABTilglks <- refABTilglks[!tiny_troublemakers]
    expect_true(all.equal(refABTilglks, theABTilglks, check.attributes = FALSE))

    ABT <- ABT[refABTilglks > -2560 & refABTilglks < -1990 ,]

    theABThetzpoors <- hetzminoo_fancy(ABT, 'poor')
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refABThetzpoors <- theABThetzpoors
    ## save(refABThetzpoors, file = "ABT_referencehetzpoors.Rda")

    base::load("ABT_referencehetzpoors.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refABThetzpoors, theABThetzpoors, check.attributes = FALSE))

    ABT <- ABT[theABThetzpoors > 0.105 & theABThetzpoors < 0.180 ,]

    theABThetzriches <- hetzminoo_fancy(ABT, 'rich')
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refABThetzriches <- theABThetzriches
    ## save(refABThetzriches, file = "ABT_referencehetzriches.Rda")

    base::load("ABT_referencehetzriches.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refABThetzriches, theABThetzriches, check.attributes = FALSE))

    ABT <- ABT[theABThetzriches > 0.143 & theABThetzriches < 0.210 ,]

    ABT <- kin_power(ABT, k = 0.5)
    ABT <- prepare_PLOD_SPA(ABT)

    ## for the purposes of keeping the test quick-to-run, here we import a list of animals involved in high-PLOD pairs
    ## and trim the dataset to just those animals. That way, everything runs _much_ quicker, and the test for
    ## matching kinship statistics is just the same amongst the pairs where it matters.
    ## Here's the code for generating that list, for posterity's sake:

    ## theABTdups <- find_duplicates(ABT, max_diff_loci = 200)
    ## ABT <- ABT[-c(drop_dups_pairwise_equiv(theABTdups[, 2:3])), ]
    ## theABT_PLODs <- find_HSPs(ABT, keep_thresh = -10)
    ## ijs <- unique(c(theABT_PLODs$i, theABT_PLODs$j))
    ## save(ijs, file = "ABT_highPLODs.Rda")

    ## Now, let's check the consistency of our results, only for those animals
    ## involved in a high PLOD pair
    base::load("ABT_highPLODs.Rda")
    ABT <- ABT[ijs,]

    theABTdups <- find_duplicates(ABT, max_diff_loci = 300)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refABTdups <- theABTdups
    ## save(refABTdups, file = "ABT_referencedups.Rda")

    if(nrow(theABTdups) > 0) {
        ABT <- ABT[-c(drop_dups_pairwise_equiv(theABTdups[, 2:3])), ]
    }

    theABTPLODs <- find_HSPs(ABT, keep_thresh = -10, limit_pairs = 2000)
    ## HSP_histo(theABTPLODs, lb = 0, fullsib_cut = 25)
    ## PLOD_loghisto(theABTPLODs)

    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refABTPLODs <- theABTPLODs
    ## save(refABTPLODs, file = "ABT_referencePLODs.Rda")

    base::load("ABT_referencePLODs.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refABTPLODs, theABTPLODs, check.attributes = FALSE))

    ## investigate the effect of the misalignment bug on the kin sets used in the billycart: #####################################
    ## split out HSPs and FSPs/POPs
    FSPsOrPOPs <- theABTPLODs[theABTPLODs$PLOD > 150,]
    FSPsOrPOPs$ij <- paste(FSPsOrPOPs$i, FSPsOrPOPs$j)
    FSPsOrPOPs$iTargetID <- ABT$info$TargetID[FSPsOrPOPs$i]
    FSPsOrPOPs$jTargetID <- ABT$info$TargetID[FSPsOrPOPs$j]
    FSPsOrPOPs$pairID <- paste(FSPsOrPOPs$iTargetID, FSPsOrPOPs$jTargetID)

    theABTPLODs$ij <- paste(theABTPLODs$i, theABTPLODs$j)
    HSPs_2 <- theABTPLODs[!theABTPLODs$ij %in% FSPsOrPOPs$ij,]
    ## plot the empirical PLOD scores with the expected distribution
    ## of the PLOD for half-sibling pairs
    thresh <- autopick_threshold(ABT, kin = HSPs_2, fitrange_PLOD = c(10, 160),
                                 FPtol_pairs = 1, use4th = TRUE, selecto = "ML", plot_bins = 5)
    HSPs <- HSPs_2[HSPs_2$PLOD > thresh,]
    HSPs$iTargetID <- ABT$info$TargetID[HSPs$i]
    HSPs$jTargetID <- ABT$info$TargetID[HSPs$j]
    HSPs$pairID <- paste(HSPs$iTargetID, HSPs$jTargetID)

    originalHSPs <- read.csv("~/Dropbox/CSIRO/Reports/ABT all-ages kinference/CapCluster/halfsibs_allfish.csv")
    originalHSPs$pairID <- paste(originalHSPs$iTargetID, originalHSPs$jTargetID)

    table(originalHSPs$pairID %in% HSPs$pairID)

    originalPOPs <- read.csv("~/Dropbox/CSIRO/Reports/ABT all-ages kinference/CapCluster/pops_allfish.csv")
    originalFSPs <- read.csv("~/Dropbox/CSIRO/Reports/ABT all-ages kinference/CapCluster/fullsibs_allfish.csv")
    originalPOPs$pairID <- paste(originalPOPs$iTargetID, originalPOPs$jTargetID)
    originalFSPs$pairID <- paste(originalFSPs$iTargetID, originalFSPs$jTargetID)
    table(FSPsOrPOPs$pairID %in% c(originalPOPs$pairID, originalFSPs$pairID) )

    ## continue with unit testing ################################################################################################
    theABTPOPs <- find_POPs(ABT, limit_pairs = 400, keep_thresh = 0.91)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refABTPOPs <- theABTPOPs
    ## save(refABTPOPs, file = "ABT_referencePOPs.Rda")

    base::load("ABT_referencePOPs.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refABTPOPs, theABTPOPs, check.attributes = FALSE))

    theABTwtsame <- split_FSPs_from_POPs(ABT, candipairs = theABTPLODs[theABTPLODs$PLOD > 150,], gerr = 0.005)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refABTwtsame <- theABTwtsame
    ## save(refABTwtsame, file = "ABT_referencewtsame.Rda")

    base::load("ABT_referencewtsame.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refABTwtsame, theABTwtsame, check.attributes = FALSE))
    ## There's a 'call' attribute for both wtsames, and all.equal tests for equality of attributes (!)
    ## the ref was called manually, but the 'theABTwtsame' was called by tinytest::test_all, so these were in conflict.
    ## Fixed by NULLifying the 'call' attribute for both.

    ## ## one-off testing of split_FSPs_from_POPs
    ## FSPsOrPOPs <- theABTPLODs[theABTPLODs$PLOD > 150,]
    ## FSPsOrPOPs$ij <- paste(FSPsOrPOPs$i, FSPsOrPOPs$j)
    ## theABTwtsame$ij <- paste(theABTwtsame$i, theABTwtsame$j)
    ## theABTwtsame$istratum <- ABT$info$STRATA[theABTwtsame$i]
    ## theABTwtsame$jstratum <- ABT$info$STRATA[theABTwtsame$i]
    ## hist(theABTwtsame$PLOD_FP, breaks = 50)
    ## abline(v = theABTwtsame@E_FSP, col = FSPcol, lwd = 2)
    ## abline(v = theABTwtsame@E_POP, col = POPcol, lwd = 2)
    ## POPs <- theABTwtsame[theABTwtsame$PLOD_FP < 0,]
    ## FSPs <- theABTwtsame[theABTwtsame$PLOD_FP > 0,] ## yeah, these basically make sense. More within-year At pairs in the
    ## ## FSPs than I would expect, though.

    theABTwpsex <- split_FSPs_from_HSPs(ABT, candipairs = theABTPOPs)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refABTwpsex <- theABTwpsex
    ## save(refABTwpsex, file = "ABT_referencewpsex.Rda")

    base::load("ABT_referencewpsex.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refABTwpsex, theABTwpsex, check.attributes = FALSE))

    theABTvarhtp <- var_PLOD_kin(ABT$locinfo, emp_V_HSP = 900, n_meio = 3)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refABTvarhtp <- theABTvarhtp
    ## save(refABTvarhtp, file = "ABT_referencevarhtp.Rda")

    base::load("ABT_referencevarhtp.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refABTvarhtp, theABTvarhtp, check.attributes = FALSE))

    theABTPLODST <- split_HSPs_from_HTPs(ABT, candipairs = theABTPLODs)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refABTPLODST <- theABTPLODST
    ## save(refABTPLODST, file = "ABT_referencePLODST.Rda")

    base::load("ABT_referencePLODST.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refABTPLODST, theABTPLODST, check.attributes = FALSE))

    ## ## one-off split_HSPs_from_HTPs testing:
    ## prollyHSPs <- theABTPLODST[theABTPLODST$PLOD_ST > -5 & theABTPLODST$PLOD_ST < 50,]
    ## prollyHTPs <- theABTPLODST[theABTPLODST$PLOD_ST < -5,]
    ## prollyHSPs$istratum <- ABT$info$STRATA[prollyHSPs$i]
    ## prollyHSPs$jstratum <- ABT$info$STRATA[prollyHSPs$j]
    ## prollyHTPs$istratum <- ABT$info$STRATA[prollyHTPs$i]
    ## prollyHTPs$jstratum <- ABT$info$STRATA[prollyHTPs$j]
    ## stratpair <- function(prolly) { table(paste(prolly$istratum, prolly$jstratum)) }
    ## stratpair(prollyHSPs)
    ## stratpair(prollyHTPs) ## these are consistent with a larger age-gap between HTPs than HSPs ID'd through
    ## ## split_HSPs_from_HTPs
    ## theABTPLODST$ij <- paste(theABTPLODST$i, theABTPLODST$j)
    ## theABTPLODs$ij <- paste(theABTPLODs$i, theABTPLODs$j)

    ## bothplods <- join(theABTPLODST, theABTPLODs, by = "ij", type = "left")
    ## with(bothplods, plot.default(PLOD_ST ~ PLOD))
    ## with(bothplods, cor(PLOD_ST,  PLOD)) ## 0.99919 correlation. Good enough for me
}

##############################################################################################
## SHS testing  ##############################################################################
##############################################################################################

if(shanesComp) {

    base::load("~/Data/SHS_input.rda") ## only on shane's computer -- might be worth copying this across, but risky IP-wise

    dups <- find_duplicates( SHS, max_diff_loci=200, show_plot = TRUE)
    prs <- with( SHS$info, cbind(Our_sample[ dups$i], Our_sample[ dups$j], dups, stringsAsFactors=FALSE))
    tech.reps <- prs[,1] == prs[,2]
    oneofeach <- drop_dups_pairwise_equiv( dups[,2:3])
    SHS <- SHS[ -c( oneofeach),]

    SHS$locinfo$pbonzer <- kinference:::re_est_ALF(SHS)$locinfo$pambig
    theSHS6and4s <- check6and4(SHS, thresh_pchisq_6and4 = c(0.0001, 0.00001))
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSHS6and4s <- theSHS6and4s
    ## save(refSHS6and4s, file = "SHS_reference6and4s.Rda")

    base::load("SHS_reference6and4s.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSHS6and4s, theSHS6and4s, check.attributes = FALSE))

    theSHSilglks <- ilglk_geno(SHS)
    ## run once in January 2024, our benchmark of 'correct', yea, unto eternity:
    ## refSHSilglks <- theSHSilglks
    ## save(refSHSilglks, file = "SHS_referenceilglks.Rda")

    base::load("SHS_referenceilglks.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSHSilglks, theSHSilglks, check.attributes = FALSE))

    SHS <- SHS[theSHSilglks > -1450 & theSHSilglks < -1260 ,]
    SHS <- kin_power(SHS, k = 0.5)
    SHS <- prepare_PLOD_SPA(SHS)

    theSHShetzpoors <- hetzminoo_fancy(SHS, 'poor')

    ## run once in January 2024, our benchmark of 'correct', yea, unto eternity:
    refSHShetzpoors <- theSHShetzpoors
    save(refSHShetzpoors, file = "SHS_referencehetzpoors.Rda")

    base::load("SHS_referencehetzpoors.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSHShetzpoors, theSHShetzpoors, check.attributes = FALSE))

    SHS <- SHS[theSHShetzpoors > 0.14 & theSHShetzpoors < 0.19 ,]

    theSHShetzriches <- hetzminoo_fancy(SHS, 'rich')
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSHShetzriches <- theSHShetzriches
    ## save(refSHShetzriches, file = "SHS_referencehetzriches.Rda")

    base::load("SHS_referencehetzriches.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSHShetzriches, theSHShetzriches, check.attributes = FALSE))

    SHS <- SHS[theSHShetzriches > 0.212 & theSHShetzriches < 0.275 ,]
    SHS <- kin_power(SHS, k = 0.5)
    SHS <- prepare_PLOD_SPA(SHS)

    ## for the purposes of keeping the test quick-to-run, here we import a list of animals involved in high-PLOD pairs
    ## and trim the dataset to just those animals. That way, everything runs _much_ quicker, and the test for
    ## matching kinship statistics is just the same amongst the pairs where it matters.
    ## Here's the code for generating that list, for posterity's sake:

    ## theSHS_PLODs <- find_HSPs(SHS, keep_thresh = 25)
    ## ijs <- unique(c(theSHS_PLODs$i, theSHS_PLODs$j))
    ## save(ijs, file = "SHS_highPLODs.Rda")

    ## Now, let's check the consistency of our results, only for those animals
    ## involved in a high PLOD pair
    base::load("SHS_highPLODs.Rda")
    SHS <- SHS[ijs,]

    theSHSdups <- find_duplicates(SHS, max_diff_loci = 200)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSHSdups <- theSHSdups
    ## save(refSHSdups, file = "SHS_referencedups.Rda")

    if(nrow(theSHSdups) > 0) {
        SHS <- SHS[-c(drop_dups_pairwise_equiv(theSHSdups[, 2:3])), ]
    }

    theSHSPLODs <- find_HSPs(SHS, keep_thresh = 0)
    ## HSP_histo(theSHSPLODs, lb = 0, fullsib_cut = 25)
    ## PLOD_loghisto(theSHSPLODs)

    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSHSPLODs <- theSHSPLODs
    ## save(refSHSPLODs, file = "SHS_referencePLODs.Rda")

    base::load("SHS_referencePLODs.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSHSPLODs, theSHSPLODs, check.attributes = FALSE))

    theSHSPOPs <- find_POPs(SHS, keep_thresh = 0.91)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSHSPOPs <- theSHSPOPs
    ## save(refSHSPOPs, file = "SHS_referencePOPs.Rda")

    base::load("SHS_referencePOPs.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSHSPOPs, theSHSPOPs, check.attributes = FALSE))

    theSHSwtsame <- split_FSPs_from_POPs(SHS, candipairs = theSHSPOPs, gerr = 0.01)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSHSwtsame <- theSHSwtsame
    ## save(refSHSwtsame, file = "SHS_referencewtsame.Rda")

    base::load("SHS_referencewtsame.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSHSwtsame, theSHSwtsame, check.attributes = FALSE))

    theSHSwpsex <- split_FSPs_from_HSPs(SHS, candipairs = theSHSPOPs)
    ## run once in November 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSHSwpsex <- theSHSwpsex
    ## save(refSHSwpsex, file = "SHS_referencewpsex.Rda")

    base::load("SHS_referencewpsex.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSHSwpsex, theSHSwpsex, check.attributes = FALSE))

    theSHSvarhtp <- var_PLOD_kin(SHS$locinfo, emp_V_HSP = 900, n_meio = 3)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSHSvarhtp <- theSHSvarhtp
    ## save(refSHSvarhtp, file = "SHS_referencevarhtp.Rda")

    base::load("SHS_referencevarhtp.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSHSvarhtp, theSHSvarhtp, check.attributes = FALSE))

    theSHSPLODST <- split_HSPs_from_HTPs(SHS, candipairs = theSHSPLODs)
    ## run once in January 2023, our benchmark of 'correct', yea, unto eternity:
    ## refSHSPLODST <- theSHSPLODST
    ## save(refSHSPLODST, file = "SHS_referencePLODST.Rda")

    base::load("SHS_referencePLODST.Rda")  ## avoid conflict with renv::load.
    expect_true(all.equal(refSHSPLODST, theSHSPLODST, check.attributes = FALSE))
}
