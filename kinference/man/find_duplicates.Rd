\name{find_duplicates}
\alias{find_duplicates}
\alias{find_HSPs}
\alias{find_POPs}
\alias{\%upto\%}
\alias{find_POPs_lglk}
\title{Kin-finders for loads-of-SNPs datasets}
\description{
These take a \code{snpgeno} dataset that has been processed as far as \code{\link{check6and4}} (and for HSPs, \code{\link{kin_power}}) and find various relations between the samples. Relationships include duplicates (DUPs/dupes/dups), parent-offspring pairs (POPs) and half-sibling pairs (HSPs) or other 2nd-order kin, plus of course unrelated pairs (UPs). You can specify the same or different subsets of the \code{snpgeno} for comparison: e.g., first subset for the adults, second for the juveniles.

There are two versions aimed at POPs currently called \code{find_POPs} and \code{find_POPs_lglk}. The former uses a "weighted pseudo-exclusion" ("wpsex") statistic that allows for null alleles and is robust to genotyping errors. The latter uses a likelihood-based statistic (again allowing for nulls), but you do have to provide a guesstimate of genotyping error rate (to robustify the calculation- otherwise, a single genotyping error in a true POP could give a log-likelihood of -Inf). \code{find_POPs_lglk} is newer, easier to explain, and perhaps less "arbitrary", but we have used the "wpsex" version on all our real CKMR datasets (>10). Time will tell whether one is better/easier than the other; finding POPs ought to be pretty easy, so the results really should be the same.

\code{find_HSPs} should really be called \code{find_2OPs} because it cannot discriminate amongst second-order kin types; there is no way to distinguish genetically between HSPs, Grandparent-Grandchild Pairs, and Full-Thiatic Pairs (eg aunt/nephew) with \code{snpgeno} data alone. But, for historical reasons, it's still called \code{find_HSPs}. Note that \code{find_HSPs} can also be tricked into targeting [some] other types of kin, such as 3; see \bold{Details}, but watch out.
}
\usage{
find_duplicates(
  snpg,
  subset1 = 1 \%upto\% nrow(snpg),
  subset2 = subset1,
  max_diff_loci,
  limit_pairs = 0.5 * nrow(snpg),
  nbins = 50,
  maxbin = ncol(snpg)/2,
  show_plot = TRUE
)
find_HSPs(
  snpg,
  subset1 = 1 \%upto\% nrow(snpg),
  subset2 = subset1,
  limit_pairs = 0.5 * nrow(snpg),
  keep_thresh,
  eta = NULL,
  nbins = 50,
  minbin = NULL,
  maxbin = NULL)
find_POPs(
  snpg,
  subset1 = 1 \%upto\% nrow(snpg),
  subset2 = subset1,
  limit_pairs = 0.5 * nrow(snpg),
  keep_thresh,
  eta = NULL,
  nbins = 50,
  maxbin = NULL,
  WPSEX_UP_POP_balance = 0.99
)
find_POPs_lglk(
  snpg,
  subset1 = 1 \%upto\% nrow(snpg),
  subset2 = subset1,
  gerr,
  limit_pairs = 0.5 * nrow(snpg),
  keep_thresh,
  eta = NULL,
  nbins = 50,
  minbin = NULL,
  maxbin = NULL
)
}
\arguments{
\item{  snpg}{a \code{snpgeno} object}
\item{  subset1, subset2}{numeric vectors of which samples to use (not logical, not negative). Defaults to all of them. Iff \code{subset1} and \code{subset2} are identical, only half the comparisons are done (i.e., not \emph{i} with \emph{j} then \emph{j} with \emph{i}). Some sanity checks are made.}
\item{  max_diff_loci}{(\code{find_duplicates}) max number of discrepant 4-way genotypes to tolerate in "identical" fish. Only the pairs with fewer than \code{max_diff_loci} discrepancies will be retained. Try increasing this from say 10 upwards, and hopefully nothing much will change (though at some point things will change a lot, as you get into the non-duplicate bit of the distribution). See \emph{Duplicates} for how to remove duplicates from the data.}
\item{  limit_pairs}{Integer. Defines the \emph{maximum} number of candidate pairs to keep. Will provide a warning if the number of identified pairs equals \code{limit_pairs}.}
\item{  nbins, minbin, maxbin}{\code{find_XXX} functions summarise their pairwise comparison statistics into bins (in the part of the range where exact values are uninteresting), as well as returning specific pairs that pass the "interesting" threshold. \code{nbins} sets the number of bins, \code{minbin} sets the top value of the lowest bin (so that bin stretches from -Inf to \code{minbin} for HSPs); \code{maxbin} sets the highest. For HSPs, the minimum is 3 bins (-Inf:minbin),[minbin:maxbin),[maxbin:Inf). \code{minbin} is not used for duplicates or (at present) for POPs, since the statistics there are defined so that the lowest possible value is 0. The defaults for \code{minbin} and/or \code{maxbin} may not be what you need in all cases, so be prepared to select manually and then re-run. For duplicates, where calculations can be slow for big datasets, you can set \code{nbins=0} to disable binning and focus instead on just finding the pairs with fewer than \code{max_diff_loci} discrepancies. Each pairwise calculation normally loops over all the loci, but is aborted when the running total of discrepant loci reaches \code{maxbin} (or, if \code{nbins=0}, when it reaches \code{max_diff_loci}), thus saving considerable time. It is therefore not sensible to have \code{maxbin<max_diff_loci} (think about it!).}
\item{  show_plot}{whether to plot log histogram. Regardless, plot will not be shown if other arguments would lead to stupid result (e.g. no bins...).}
\item{  keep_thresh}{(\code{find_HSPs} and \code{find_POPs}) is the analog of \code{max_diff_loci} for \code{find_duplicates}. It determines which pairs to retain for individual inspection. For \code{find_HSPs}, this is the lowest retained PLOD; for \code{find_POPs} (currently), it's the highest retained \code{wpsex}. Set it with the aim of including "anything interesting" (ie not \emph{missing} any interesting pairs) and expect false positives - i.e. be willing to have some weaker kin in there, and to subsequently filter those out yourself "manually", as per vignette. For HSPs, values like 0 (near the HTP mean) or -5 are a good start. For POPs, experiment (at a total guess, try 0.1). You may have to re-run the function a coupla times if you have been too brutal or too generous here - though "too generous" can be fixed post hoc just by filtering the result, as long as you haven't generated tooooo many pairs (see parameter \code{limit_pairs}).}
\item{  eta}{(\code{find_HSPs} and \code{find_POPs}) Not essential; limit for calculating empirical mean and var PLOD, to compare with theoretical \code{mean_UP} and \code{var_UP}. If you care about this (and you might not, since for \code{find_HSPs} the observed/expected binwise comparison is perhaps clearest), then set it to somewhere above 0 that should include almost all UPs and exclude most strong kin; that's an \emph{upper} limit for HSPs, and a \emph{lower} limit (currently) for POPs. The defaults are guesses, and might not do that adequately, so an apparent mismatch may not actually matter; be prepared to look at the histograms and think. The general idea is that the number of UPs should dominate any other kin-type in large sparsely-sampled datasets, so "contaminating" the empirical UP statistics with a few weakish kin at the top end shouldn't mess up mean and variance too much.}
\item{  WPSEX_UP_POP_balance}{(\code{find_POPs}) loci receive a weight which is proportional to (difference in probability of pseudo-exclusion between UP and POP) / (variance of indicator of pseudo-exclusion). But, should this be variance assuming UP or POP? \code{WPSEX_UP_POP_balance} sets the balance; bigger values make it more UPpity, so placing more emphasis on avoiding false-positives (which is probably the Right Thing To Do). 0.99 could be completely fine... (but hopefully \code{WPSEX_UP_POP_balance} won't affect the result much anyway.)}
\item{  gerr}{(\code{find_POPs_lglk}) Genotyping error rate (apart from any AA/AO-type errors)- which had better be a small number. You have to pick it yourself, but it is only used to "robustify" the lglk-based (PO)PLOD for testing POPs vs UPs, and thus can be a rough guesstimate. FWIW we have used 0.01 (i.e. 1\\ replicate analysis suggests, but is "safe" while still being small enough not to muck up overall statistical performance. You should really do the same thing yourself, and if you are very paranoid then try sensitivity analyses; but in practice, the results of \code{find_POPs_lglk} are liable to be so clear-cut that you may not feel it necessary to try more than one small value...)}
}
\details{
Some categories will "catch" others (e.g. \code{find_HSPs} will certainly include any POPs too), so you may need the splitter routines such as \code{split_POPs_from_HSPs} afterwards. The safest general-purpose strategy - but often \emph{not} the most sensible, if your data is nicely organized and you know what you want - would be:

\code{find_duplicates} and then get rid of them

\code{find_HSPs} to get \emph{all} kin (though you will usually have to sacrifice some HSPs to false-neg because you'll need a threshold)

\code{split_POPs_from_HSPs} to split HSPs from POPs/FSPs

\code{split_POPs_from_FSPs} to split the latter.

The non-splitter functions, i.e. \code{find_XXX}, might be run on huge numbers of samples, entailing a \code{choose(huge, 2)} number of comparisons. You don't want all those individual comparison results, and your computer certainly wouldn't enjoy trying to keep them! So the general idea is to set a threshold for what constitutes "maybe worth keeping individually" (that you expect will be generous enough to contain everything you \emph{do} want, plus some dross), and then to retain just binned counts of the relevant comp statistic for all comps (usually, the vast majority) which don't make your threshold.

In addition, the \code{limit_pairs} argument is there to prevent your computer locking out with bazillions of unwanted pairs (in case you guess the bin limit inapproriately); the comparisons will be stopped if \code{limit_pairs} is hit, with a warning. In that case, you probably need to change a threshold, or re-run with larger \code{limit_pairs}. The default isn't meant to correspond to any biomathematical logic, it's just to stop blue smoke coming out your USB ports.

For \code{find_duplicates}, there are at least two different use-cases. First, you might want an initial run on a non-too-large subset of your data, to check that dups \emph{can} be clearly distinguished and to look at typical extent of genotyping errors (based on clear duplicates that don't match at every locus). For that, you can set \code{nbins} and choose some reasonable guess as to \code{max_diff_loci} (say, 5\\ loci). Because you set \code{nbins>0}, \emph{every} pair (almost...) gets checked at \emph{all} loci, so it can be slow. Thus, if you have done this before and have a good sense of "how bad can a real duplicate be?", then set \code{nbins=0} (and \code{max_diff_geno} to a small but safe value that won't miss any realistic duplicate-with-genotyping-error) so it will abort a comparison early as soon as it reaches \code{max_diff_geno} differing loci. That saves a \emph{lot} of time on big datasets! You won't get a histo of number-of-diffs, but you don't need one for that use-case. The "almost" is that \code{find_duplicates} uses "transitivity" (if A is a dup of B and of C, then we don't need to check B vs C), so it only counts differences for not-yet-known duplicates \emph{based on} \code{max_diff_loci}. To discard duplicates and to find entire equivalence-classes of duplicates, e.g. from a control specimen included in numerous plates, see \code{drop_dup_pairwise_equiv}.

\code{find_HSPs} relies on pre-computed values of "LOD" and "PUP" that have been set by \code{\link{kin_power}}. Normally you would call the latter with \code{k=0.5}, since that's what HSPs are. However, the devious user can try \emph{different} values of \code{k}- which is how \code{find_POPs_lglk} works- and then the target of \code{find_HSPs} will become "kin with that value of \code{k}". Be very careful!!!
\subsection{Kinformation}{
The idea is that kin-finding is based on a statistic and a threshold \code{eta}, where the latter is chosen to keep false-positives down to a user-specified level. Anything "beyond" \code{eta} will be treated as a kin-pair ("beyond" depends on how the statistic is defined, i.e. whether a kin-pair should come out very low or very high). However, you're also likely to want to look post hoc at the distro of computed statistics \emph{near} \code{eta}, to see whether separation is as clean (or otherwise) as expected - and also very unbeyond \code{eta} into the zone where UPs are entirely dominant, to check that theory is OK. So, as well as returning the "interesting" pairs that have a statistic close to or on the non-UP size of \code{eta}, the POP and HSP versions also return \emph{summaries} of the distribution of the statistic. The thing is that there will be zillions of statistics from UPs - enough to blow out computer memory - and they are not individually interesting. Specifically, the main things returned are:

mean and variance of stats. Computation is restricted to those on the UP-side of \code{eta} (which is nearly all of them, usually) in order to avoid distortion from non-UP cases. The latter will often be so rare that distortion would be negligible - but means and variances are not "robust".

counts of binned stats, regardless of whether above or below \code{eta}. The bins are set based on SPAs to the theoretical distributions, and chosen so that an equal number of UP-pairs should fall into each bin.

cases where the stat is "interesting", i.e. on the non-UP side of \code{keep_thresh}, as a \code{data.frame}. See \emph{Value} for details

The process is controlled by three numbers: \code{nbins} for number of bins, \code{eta} itself, and some nearby threshold \code{keep_thresh} on the UP-side of \code{eta} (it will be automatically set to \code{eta} otherwise) to determine which pairs are explicitly retained for your inspection. There are two ways to specify \code{eta} and \code{keep_thresh}. Usually, you would start with the indirect method, where you choose the predicted-false-positive proportion of UP-pairs via the parameter \code{one_in_X_eta}, and \code{rough_n_pairs_to_keep}. The routines then use SPAs to the corresponding values of \code{eta} and \code{keep_thresh}; the returned value of \code{eta} is what you can subsequently use to make the actual kin-decisions yourself after the event (by subsetting the "interesting" pairs, comparing the statistic for each pair to \code{eta})- assuming that observed does match expected.

But, sometimes it doesn't. In that case, the predicted values of \code{eta} and \code{keep_thresh} may be way off the mark, and lead to retaining faaar too few or too many pairs. If so, then look at the histogram of retained statistics from an initial run, and try setting \code{eta} and/or \code{keep_thresh} manually, rather than futzing around with the indirect parameters until you get what you were after.
}
}
\value{
A \code{data.frame} with extra attributes (see below) and at least 3 columns: statistic \code{PLOD} or \code{wpsex} or \code{ndiff} (number of mismatching genotypes), then \code{i} and \code{j} (index in \code{subset1} and \code{subset2} respectively of the first pair-member). Note that \code{i} and \code{j} refer to the \emph{subsets}, not to the rows of the original \code{snpg}- so if you are comparing subsets within the same overall \code{snpg}, then you have to adjust accordingly. The attributes in all cases include \code{bins} (upper boundaries), some kind of count statistic for number of comparisons in each bin (names vary), \code{binprobs} (theoretical CDF for UPs in \code{find_HSPs}; should exist for \emph{POPs} (not UPs) in \code{find_POPs} (the \code{wpsex} version) but currently doesn't), some of the input parameters, and the \code{call} that invoked the function. \code{find_POPs} adds a column named \code{nABOO}, showing the number of AB/OO exclusions for that potential POP. This is a useful additional diagnostic; it should be close to 0 for true POPs (it can only result from genotyping error or mutation, whereas AAO/BBO can result from nulls). For UPs, I was seeing values typically in the low 20s, which is pretty good separation. \code{find_HSPs} and \code{find_POPs} have a bunch of extra attributes which should be reeeeasonably clear. For \code{find_HSPs}, \code{mean_sub_PLOD} and \code{var_sub_PLOD} are the empirical means & var below \code{eta}, andy they should be close to \code{mean_UP} and \code{var_UP} \emph{iff} \code{eta} has been chosen sensibly. For \code{find_POPs}, the same goes for \code{mean_wpsex_hi} and \code{var_wpsex_hi}. For duplicates, not \emph{all} pairwise duplicates are recorded, unless the subsets are different - otherwise you could have quadratic horror of enormous numbers of pairs arising from a cluster of say 100 identical controls! Since "duplication" is transitive (ie if i & j are the same, and i & k are the same, then j & k must also be the same), only the necessary ones are recorded to allow you to filter out yourself afterwards. e.g., if samples 1, 3, 5, and 6 are all duplicates, you'll get this (without "ndiff" column) [FORMATTING SCREWED UP PROBABLY ROXYGEN'S FAULT]:
\item{  i j:::}{}
\item{  3 1:::}{}
\item{  4 3:::}{}
\item{  6 4}{but you won't see the pairings for 1/4, 1/6, 3/6. If you just want to strip out all duplicates bar one in each group (and you don't care which one is kept), then you can use the function \code{drop_dups_pairwise_equiv} - see \emph{Examples}. For POPs and HSPs, the following are also returned as attributes (that can be accessed by \code{@} if \code{atease} is loaded). The main point is that the "boring" below-threshold pairs get put into bins and are not kept individually. The names sometimes change depending on which statistic is being used.}
\item{  - eta}{false-positive cutoff to be applied to the statistic in question (automatically done if \code{rough_n_pairs_to_keep==NA}, or up to you if not). Variance of the stat will only be calculated from values to the "UP side" of \code{eta}. However, the set of retained pairs/individuals is actually controlled by...}
\item{  - keep_thresh}{the cutoff used to retain "interesting" pairs. Usually obvious from the range of statistic values.}
\item{  - mean_sub_stat, var_sub_stat (where stat is PLOD, wpsex, or ndiff)}{empirical values for the statistic when it is below \code{eta} (ie nearly always).}
\item{  - mean_theory, var_theory}{of the statistic, to compare to previous.}
\item{  - n_stat_in_bin (where stat is PLOD, wpsex, or ndiff)}{number of pairs whose statstic fell within the range of each bin}
\item{  - bins}{cutpoints for the bins. These should be quantiles, according to the SPA; so if practice matches theory, the numbers-per-bin should all be similar.}
}
\examples{
## find_duplicates
library( mvbutils)
library( atease)
gbasics::define_genotypes()
x <- matrix(sample(c("AAO", "AB", "BBO", "OO"), 10000, TRUE, prob = c(0.2, 0.45, 0.2, 0.15)), nrow = 100, ncol = 100)
minisnpg <- gbasics::snpgeno(x = x, diplos = genotypes4_ambig)
## seed some duplicates in. Sample 2 will be a copy of sample 1,
## and samples 21, 22, and 23 will be an exactly-matching group.
minisnpg[2,] <- minisnpg[1,]
minisnpg[ 22,] <- minisnpg[ 23,] <- minisnpg[ 21,]
## seed some inexact duplicates in. Make sample 3 nearly match
## sample 4, and the same for samples 13 and 14
minisnpg[4,1:80] <- minisnpg[3,1:80]
minisnpg[14,1:80] <- minisnpg[13,1:80]
## find exact duplicates
exact <- find_duplicates( minisnpg, max_diff_loci = 0)
exact ## finds ij pairs \{1, 2\} and the group \{21, 22, 23\}
## as pairs \{21, 22\} and \{22, 23\}
## find (in-)exact duplicates
inexact <- find_duplicates( minisnpg, max_diff_loci = 25)
inexact ## finds the exactly-matching pairs as before, plus
## the inexactly-matching row pairs \{3, 4\} and \{13, 14\} with
## >0 differences
## to remove duplicates, keeping only one member of each
## group, use drop_dups_pairwise_equiv
droppers <- drop_dups_pairwise_equiv( inexact[,2:3])
droppers ## note that _all but one_ of each _group_ of
## (near-)duplicates is included in 'droppers'
minisnpg_nodups <- minisnpg[ - droppers,] ## drop all-but-one
## of each set of duplicates
## find_HSPs (PLOD_HU)
data( bluefin)
bluefin <- bluefin[,1:500] ## make dataset smaller for speed in R CMD check
## stripped-down data-cleaning for example - see the
## vignette for approach for real data!
pvals <- check6and4( bluefin, thresh_pchisq_6and4 = c( 0.001, 0.0001))
bluefin_1 <- bluefin[ , pvals$pval4 > 0.01]
ilglks <- ilglk_geno( bluefin_1)
bluefin_2 <- bluefin_1[ ilglks > -1050,]
bluefin_3 <- est_ALF_ABO_quick( bluefin_2)
bluefin_4 <- bluefin_3[ , bluefin_3@locinfo$pbonzer[,"B"] > 0.02]
bluefin_5 <- kin_power( bluefin_4, k = 0.5)
dups <- find_duplicates( bluefin_5, max_diff_loci = 20)
bluefin_6 <- bluefin_5[ -c(drop_dups_pairwise_equiv( dups)) ,]
hsps <- find_HSPs( bluefin_6, keep_thresh = 0)
HSP_histo( hsps, lb = 0, fullsib_cut = 75)
head( hsps) ## PLODs and row numbers for each pair member in bluefin
hsps@mean_HSP ## mean expected PLOD for true HSPs. mean_UP (unrelated),
## mean_POP (parent-offspring), mean_FSP (full-sibling) follow the
## same format
## subset comparisons
## limit comparisons to those between animals on plate 3 and animals on
## the other two plates. Useful when, e.g., looking for kinships between
## adults and juveniles, but not kinships between adults and other
## adults. For demonstration, treat plate 3 as adults and plates 1 and 2
## as juveniles (they're not really).
adults <- which( bluefin_6@info$Our_plate == "plate3")
juvs <- which( bluefin_6@info$Our_plate == "plate1" | bluefin_6@info$Our_plate == "plate2" )
hsps_subset <- find_HSPs( bluefin_6, subset1 = adults, subset2 = juvs, keep_thresh = 0)
hsps_subset$iplate <- bluefin_6@info$Our_plate[ hsps_subset$i]
hsps_subset$jplate <- bluefin_6@info$Our_plate[ hsps_subset$j]
## all pairs are between a plate 3 fish and a fish on plate 1 or plate 2
head(hsps_subset)
## find_POPs and find_POPs_lglk
pops <- find_POPs(bluefin_6, keep_thresh = 0.1)
## keep_thresh = 0.1 is not magic, but often OK. If you
## don't see a big spike out to the right of your plot,
## you should generally set keep_thresh higher
pops@mean_UP ## expected mean 'wpsex' for unrelated pairs
pops@var_UP ## expected variance 'wpsex' for unrelated pairs
hist(pops$wpsex) ## true POPs should have a wpsex near zero, and
## a nABOO of exactly zero, unless sequencing error has occurred
with( pops, plot.default( wpsex ~ nABOO))
## plot the full distribution from binned records
plot( pops@bins, pops@n_wpsex_in_bin, type = "s")
## add expected mean wpsex for unrelated pairs
abline( v = pops@mean_UP, col = kinPalette()["UP"], lwd = 2)
pops_lglk <- find_POPs_lglk(bluefin_6, keep_thresh = 0, gerr = 0.01)
pops_lglk@mean_UP ## as before
pops_lglk@mean_POP ## both mean_POP and mean_UP, because E(stat | k = UP) is
## no longer zero
## plot the full distribution from binned records
PLOD_loghisto( hsps = pops_lglk, HSP = FALSE, FSP = FALSE)
}
\keyword{misc}
