library(kinference)
library(gbasics)

enc <- gbasics::get_genotype_encoding()
snpg <- gbasics::snpgeno(
  matrix(c("AB", "AAO", "BBOO", "AB", "AAO", "BBOO"),
         nrow = 3, ncol = 2),
  diplos = enc$genotypes3,
  info = data.frame(Our_sample = c("sample_A", "sample_B", "sample_C")),
  locinfo = data.frame(Locus = c("locus_1", "locus_2"))
)
snpg$locinfo$pbonzer <- cbind(A = c(.5, .5), B = c(.5, .5), C = 0, O = 0)
snpg$locinfo$useN <- rep(3L, 2)

# No snpg4 object is supplied: the function must initialize it from snpg.
observed <- kinference::ilglk_geno(snpg, showPlot = FALSE)
tinytest::expect_equal(length(observed), 3L)
tinytest::expect_true(all(is.finite(observed)))
# Two loci, no null alleles: P(AB)=0.5, P(AAO)=P(BBOO)=0.25.
tinytest::expect_equal(as.numeric(observed), 2 * log(c(.5, .25, .25)))
