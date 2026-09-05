library(kinference)
library(gbasics)

set.seed(17)
enc <- gbasics::get_genotype_encoding()
snpg <- gbasics::snpgeno(
  matrix(sample(c("AAO", "AB", "BBO"), 8000, replace = TRUE,
                prob = c(.25, .5, .25)), nrow = 80),
  diplos = enc$genotypes4_ambig,
  info = data.frame(Our_sample = paste0("sample_", seq_len(80)))
)
snpg <- gbasics::with_rowid_field(snpg, "Our_sample")
snpg <- kinference::kin_power(kinference::est_ALF_nonulls(snpg), k = .5)
ids <- snpg$info$Our_sample

check_identifiers <- function(fun, args) {
  numeric_pairs <- suppressMessages(do.call(fun,
    c(list(snpg = snpg, ij_numeric = TRUE), args)))
  named_pairs <- suppressMessages(do.call(fun,
    c(list(snpg = snpg, ij_numeric = FALSE), args)))
  tinytest::expect_true(nrow(numeric_pairs) > 0L)
  tinytest::expect_identical(named_pairs$i, ids[numeric_pairs$i])
  tinytest::expect_identical(named_pairs$j, ids[numeric_pairs$j])
  # Includes nABOO for parent-offspring pairs.
  for (column in setdiff(names(numeric_pairs), c("i", "j"))) {
    tinytest::expect_identical(named_pairs[[column]], numeric_pairs[[column]])
  }
}

check_identifiers(kinference::find_duplicates,
  list(max_diff_loci = 100, limit_pairs = 100,
       maxbin = 100, nbins = 0, show_plot = FALSE))
check_identifiers(kinference::find_dups_with_missing,
  list(max_diff_ppn = 1, limit = 10000))
check_identifiers(kinference::find_HSPs,
  list(keep_thresh = -100, limit_pairs = 100))
check_identifiers(kinference::find_POPs,
  list(keep_thresh = 100, limit_pairs = 100))

# Subsets must map back to original sample identifiers, not subset positions.
subsets <- list(subset1 = c(3, 8, 11), subset2 = c(20, 30, 40))
check_identifiers(kinference::find_HSPs,
  c(subsets, list(keep_thresh = -100, limit_pairs = 100)))
check_identifiers(kinference::find_POPs,
  c(subsets, list(keep_thresh = 100, limit_pairs = 100)))
