###

library(gbasics)
#library(kinference)
library(devtools)
load_all("~/current/kinference")

# required to inspect @ elements
library(atease)


# need to impress on people here that putting the RIGHT data.frame in the right function is essential or BAD THINGS will happen, also keeping track of what stage a given data.frame is at

load("dave-practice.rda")
load("fsps_for_dave.rda")
dave <- rbind( dave, more_for_dave)

# note on adding in extra data and the subset arguments

# what are we **doing** here? Removing totally unlikely fish?
# with other data sets need to fiddle with the xlims?
ill <- ilglk_geno(dave)

# need a hetzminoo_fancy thing here? Need to setup prepare_PLOD_SPA too?

# get duplicates
# from manual max_diff_genos "Try increasing this from say 10 upwards"
# nice little story here where if you set to 40 you get a spurious result
# x axis is "number of loci different" -- say something about QC
dupes <- find_duplicates(dave, max_diff_genos=100)
#hist(dupes$ndiff)

# get the indices to drop
## CHECK THIS function
droppies <- drop_dups_pairwise_equiv(dupes[,2:3])

dave_nodupes <- dave[-droppies,]
# no duplicates any more
nodupes <- find_duplicates(dave_nodupes, max_diff_genos=100)

# want to look at "non-technical duplicates" here too?
# drop_dups_pairwise_equiv(..., want_groups=TRUE)
# look at $Our_sample

# find POPs
pops <- find_POPs(dave_nodupes, one_in_X_eta=1e3, nbins=50,
                     rough_n_pairs_to_keep=5000)
#hist(pops$bigs$wpsex)

# something about nABOO here?

#all_pops <- 0.051

# find HSPs
dave_nodupes <- prepare_PLOD_SPA(dave_nodupes)

is_juve <- toupper(substring(dave_nodupes@info$Our_sample, 3, 4))=="PL"
dave_juves <- dave_nodupes[is_juve, ]

# bins found using predict_hsp_util???
dd <- find_HSPs(dave_juves, bins = seq(-120, 120, by = 5),
                keep_thresh = -5, eta = 10, keep_n=1000)

# usually bump up number of classes
#hist(dd$PLOD, nc=30)
#abline(dd$mean thing)
# could remove i and j s for POP and FSP pairs

# histogram
#plot(dd@bins, dd@n_PLODs_in_bin)
# log-histogram
#plot(dd@bins, log(dd@n_PLODs_in_bin), type="b")
# HTPs around 0 in this?

# for get chain
#table(dd[dd$PLOD >0 & dd$PLOD <30,][,c("i","j")])

# getting FSPs

# HSP pairs in a matrix 2 cols
HSPs <- as.matrix(dd[, 2:3])

pp <- find_FSPs_from_HSPs(dave_juves, HSPs)

# remove FSPs from HSPs

# calculate the "upper variance" not including the FSPs
# variance sum((PLOD - EPLOD_HSP)^2)
# HTP/HCP
# UP -- find_HSPs set one_in..
# HCP/HTP -- 

# false negative rate

# get_chain stuff here?

# example table of numbers of HSPs per year
# number of comparisons


