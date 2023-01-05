## These are the unit tests for package 'kinference'

library(atease)
library(kinference)
library(gbasics) ## shouldn't be necessary to do this manually; hack for now to get in snpgeno.default
library(testthat) # pronounced 'testy hat'
library(renv)

## absolute basic mechanics

test_that("check that kinference has loaded", {
    expect_true("kinference" %in% sessionInfo()$otherPkgs$kinference)
})

test_that("objects of class snpgeno can be generated, subset replacement works, and diplos functions as intended", {
    info <- data.frame(Our_sample = c("sample1", "sample2"), location = c("loc1", "loc2") )
    locinfo <- data.frame(Locus = c("L1", "L2"), n_alleles = c(2,2))
    minisnpg <- snpgeno(2, 2, diplos = c("AAO", "AB", "BBO", "OO"), info = info, locinfo = locinfo)
    expect_equal(class(minisnpg), "snpgeno")
    expect_true(minisnpg[1,1] == "AAO")
    minisnpg[1,2] = "AB"
    minisnpg[2,1] = "BBO"
    minisnpg[2,2] = "OO"
    as.raw(minisnpg[1,2])

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
})


## PLOD calculations

test_that("find_HSPs gets the right PLODs for a really small fake dataset", {
    set.seed(1111)
    info <- data.frame(Our_sample = c(paste("sample", 1:10)), location = c(rep("loc1", 5), rep("loc2",5)) )
    locinfo <- data.frame(Locus = c(paste("L", 1:10)), n_alleles = c(rep(2,10)))
    smallsnpg <- snpgeno(10, 10, diplos = c("AAO", "AB", "BBO", "OO"), info = info, locinfo = locinfo)
    smallsnpg[] <- sample(c("AAO", "AB", "BBO", "OO"), 100, replace = TRUE, prob = c(0.45, 0.3, 0.21, 0.04))

    smallsnpg$locinfo$pbonzer <- re_est_ALF(smallsnpg)$locinfo$pambig
    smallsnpg$locinfo$snerr <- 0
    smallsnpg$locinfo$useN <- 4

    test1 <- hsp_power(smallsnpg, k = 0.5)
    test2 <- prepare_PLOD_SPA(test1)  ## error here. Probably related to atease removal
})


test_that("we get the right number of HSPs, POPs, and FSPs for SBT", {

})


