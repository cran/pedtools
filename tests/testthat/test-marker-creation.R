
test_that("marker() catches invalid genotype assignments", {
  x = nuclearPed(1)
  expect_error(marker(x, 1), "Genotypes must be named")
  expect_error(marker(x, "4" = 1), "Unknown ID label")
  expect_error(marker(x, "1" = ""), "Unclear genotype")
  expect_error(marker(x, "1" = "1/2/3"), "Unclear genotype")
  expect_error(marker(x, "1" = 1:3), "Unclear genotype")
  expect_error(marker(x, "1" = 1, "1" = 2), "Multiple genotypes given for individual")
  expect_error(marker(x, "1" = "3/4", alleles = 1:2), "Invalid allele")
  expect_error(marker(x, "1" = 1, allelematrix = 0), "When specifying genotypes, `allelematrix` must be NULL")

  expect_error(marker(x, geno = 1), "Genotypes must be named")
  expect_error(marker(x, geno = c("4" = 1)), "Unknown ID label")
  expect_error(marker(x, geno = c("1" = "")), "Unclear genotype")
  expect_error(marker(x, geno = c("1" = "1/2/3")), "Unclear genotype")
  expect_error(marker(x, geno = c("1" = 1, "1" = 2)), "Multiple genotypes given for individual")
  expect_error(marker(x, geno = c("1" = "3/4"), alleles = 1:2), "Invalid allele")
  expect_error(marker(x, geno = c("1" = 1), allelematrix = 0), "When specifying genotypes, `allelematrix` must be NULL")
})

test_that("marker() assigns genotypes with `geno`", {
  s = singleton(1)
  expect_equal(marker(s, geno = "a/b"), marker(s, '1' = "a/b"))
  expect_equal(marker(s, geno = "0/1"), marker(s, '1'= c(NA, 1)))

  x = nuclearPed(1)
  expect_equal(marker(x, geno = c(NA, "a/b", 0)), marker(x, '2' = "a/b"))
  expect_equal(marker(x, geno = c('3'=3, '2'=2, '1'=1)), marker(x, '1'=1, '2'=2, '3'=3))
})

test_that("empty markers are created correctly", {
  x = nuclearPed(2)

  # Adding an empty SNP (all genotypes are missing):
  expect_is(marker(x), "marker")
  expect_is(marker(x, '1'=0), "marker")
})

test_that("alleles and freqs are sorted together", {
  x = nuclearPed(1)
  p = .8; q = 1-p

  # empty
  m1 = marker(x, alleles=2:1, afreq=c(q, p))
  expect_equal(alleles(m1), c('1','2'))
  expect_equal(unname(afreq(m1)), c(p, q))

  # the "2" allele observed
  m2 = marker(x, '1'=2, alleles=1:2, afreq=c(p, q))
  expect_equal(alleles(m2), c('1','2'))
  expect_equal(unname(afreq(m2)), c(p, q))

})

test_that("marker() catches various errors", {
  x = nuclearPed(1)
  expect_error(marker(x, alleles = c(1,1)),
               "Duplicated allele label: 1")
  expect_error(marker(x, afreq = c(a=.1, a=.9)),
               "Duplicated allele label: a")
  expect_error(marker(x, afreq = c(.2,.3,.5)),
               "When `alleles` is NULL, `afreq` must be named")
  expect_error(marker(x, alleles = 1, afreq = c(a = 1)),
               "Argument `alleles` should not be used when `afreq` has names")
})
