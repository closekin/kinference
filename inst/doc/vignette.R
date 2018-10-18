## ------------------------------------------------------------------------
library( kinference) ## don't overwrite it with the one from the Rmvb repo!`

## ------------------------------------------------------------------------
    library( mvbutils)
    library( atease)
    library( xtable)
    library( mixtools)
    library( Oarray)
    library( Oarray.support)
    library( vecless)
    library( tools)  # need or else get error "could not find function "md5sum""
    library( gbasics)  #has sqr, logit and inv.logit, ridder
    library( genocalldart)
    ## library(handy2) ## not requred under the current build of kinference

## ------------------------------------------------------------------------
load('~/Data/kinference/geno2017.rda')
geno_opt <- "2017"

## ------------------------------------------------------------------------
#g1<- load_whopper('Report-DTu17-2599_part1.csv', which_rows = c(
#         well = 1,
#         dartjob = 2,
#         ourplate = 3,
#         fishname = 4,
#         fishtot = 5),
#     gobj=g0[0,],
#     cutoff_npoly=0.5, dropped_fish_file='DTu17-2599_part1_dropped_fish.txt',max_fish=2500)
#g2<- load_whopper('Report-DTu17-2599_part2.csv', which_rows = c(
#         well = 1,
#         dartjob = 2,
#         ourplate = 3,
#         fishname = 4,
#         fishtot = 5),
#     gobj=g1,
#     cutoff_npoly=0.5, dropped_fish_file='DTu17-2599_part2_dropped_fish.txt',max_fish=2500)
#g3 <- load_whopper('Report-DTu17-2711.csv', which_rows = c(
#         well = 1,
#         dartjob = 2,
#         ourplate = 3,
#         fishname = 4,
#         fishtot = 5),
#     gobj=g2,
#     cutoff_npoly=0.5, dropped_fish_file='DTu17-2711_dropped_fish.txt',max_fish=2500)
#geno2017<- g3

## ------------------------------------------------------------------------
multhresh <- 0.001

geno_okhetz <- dump_badhetz_fish(geno2017, multhresh_badhetz_fish = nrow(geno2017)*multhresh,badhetz_hist_pars=list(nclass=150,xlim=c(0,600)), 'hetzminoo') 
dim(geno_okhetz)

## ------------------------------------------------------------------------
dups <- find_duplicates(geno_okhetz,max_diff_genos=150)
rev(sort(dups$ndiff))[1:20]
droppies <- drop_dups_pairwise_equiv(cbind(dups$i,dups$j))
# Mark doesn't get rid of the 3 remaining controls, but I think we probably should
controls <- unique(c(grep('CR',geno_okhetz@info$Our_sample,ignore.case=T),grep('control',geno_okhetz@info$Our_sample,ignore.case=T)))
controls[controls%not.in%droppies]
droppies <- unique(c(droppies,controls))
geno_nodups <- geno_okhetz[ -droppies, ]
dim(geno_nodups)

## ------------------------------------------------------------------------
ilglk <- ilglk_geno(geno_nodups)
ilglk_lb <- (-1550)  ## SB: arbitrary lower bound here?
abline(v=ilglk_lb,col=6,lwd=2)
geno_nodups_goodfish <- geno_nodups[ilglk>ilglk_lb,]
dim(geno_nodups_goodfish)

## ------------------------------------------------------------------------
hetz.rich<- hetzminoo_fancy(geno_nodups_goodfish, 'rich')
hetz.poor<- hetzminoo_fancy(geno_nodups_goodfish, 'poor',hist_pars=list(nclass=100))
poor_lb <- 0.196
abline(v=poor_lb,col=6,lwd=2)
geno_nodups_gooder <- geno_nodups_goodfish[hetz.poor>poor_lb,]
# Note that "gooder" was Mark's choice of terminology!
dim(geno_nodups_gooder)

## ---- eval = FALSE-------------------------------------------------------
#  pdf(paste('check6and4.pdf'))
#  check6and4(geno_nodups_gooder,c(1e-4,1e-3))  ## not shown
#  dev.off()

## ---- eval = FALSE-------------------------------------------------------
#  adults<- grep('SbIN',geno_nodups_goodfish@info$Our_sample,ignore.case=T)
#  juvs <- grep('SbPL',geno_nodups_goodfish@info$Our_sample,ignore.case=T)
#  
#  pops <- find_POPs(geno_nodups_goodfish, subset1=adults, subset2=juvs, one_in_X_eta=1e6,
#                    keep_n=1000,
#                    nbins=50) ## find_POPs has a non-zero runtime, so is not evaluated in the vignette
#  ## save(pops,file=paste('pops2018_',geno_opt,'.rda',sep=''))
#  
#  pdf(paste('pops_histo_',geno_opt,'.pdf',sep=''))
#  hist(pops$wpsex,nc=40)
#  abline(v=0.045,col=2,lwd=2)
#  dev.off()

## ---- eval = FALSE-------------------------------------------------------
#  sum(pops$wpsex<0.045)

## ------------------------------------------------------------------------
fish_opt <- "good"
#fish_opt <- "gooder"
if(fish_opt=="good") geno_for_hsps<- geno_nodups_goodfish
if(fish_opt=="gooder") geno_for_hsps<- geno_nodups_gooder

## ---- eval = FALSE-------------------------------------------------------
#  is.juv <- grepl('SbPL',geno_for_hsps@info$Our_sample,ignore.case=T)
#  geno_juves<- geno_for_hsps[is.juv,]
#  dim(geno_juves)
#  
#  fname=paste('hsps2018_',geno_opt,'_',fish_opt,'.rda',sep='')
#  if(file.exists(fname)) load(fname) ## if already generated, don't repeat next slow step.
#  if(!file.exists(fname)) {
#    # run find_HSPs (need to run prepare_PLOD_SPA first, or else get error that Kenv is missing)
#    geno_juves <- prepare_PLOD_SPA(geno_juves)
#      hsps <- find_HSPs(geno_juves, one_in_X_eta=1e6, keep_n= 1000,
#                        bins=seq(-150,250,5)) ## find_HSPs is also not run in vignette owing to runtime
#  ##    hsps_dep <- find_HSPs_deprecated(geno_juves, one_in_X_eta=1e6, rough_n_pairs_to_keep= 200,
#  ##                                  bins=seq(-150,250,5))
#  
#  ##  save(hsps,file=fname)
#  }
#  
#  hist(hsps$PLOD[hsps$PLOD > 0], nc = 40)

