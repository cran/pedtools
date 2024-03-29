test_that("addMarker() works with `locusAttr` list", {
  attr = list(name = "M1", afreq = c(a = .2, b= .3, c = .5), chrom = "X")
  x = nuclearPed() |> addMarker(locusAttr = attr)
  expect_equal(afreq(x, 1), attr$afreq)
  expect_equal(name(x, 1), "M1")
  expect_equal(chrom(x, 1), "X")

  # ped list
  y = list(nuclearPed(), singleton(1)) |> addMarker(locusAttr = attr)
  expect_equal(afreq(y, 1), attr$afreq)
  expect_equal(name(y, 1), "M1")
  expect_equal(chrom(y, 1), "X")
})

test_that("addMarker() works on ped list", {
  x = list(nuclearPed(), singleton(4))
  y = addMarker(x, geno = c("1" = "1/2", "4" = 1), alleles = 1:2)
  expect_equal(y, list(addMarker(nuclearPed(), "1" = "1/2"), addMarker(singleton(4), 1, alleles = 1:2)))

  z = addMarker(x, "1" = "1/2", alleles = 1:3)
  expect_equal(alleles(z, 1), as.character(1:3))

  expect_equal(addMarker(x, alleles = 1:2), list(addMarker(x[[1]], alleles = 1:2),
                                                 addMarker(x[[2]], alleles = 1:2)))
})


test_that("setMarkers() creates empty markers with locusAttributes", {
  x = nuclearPed(1)
  ann = list(list(alleles=1:2, name="snp1"), list(alleles=c("a", "b")))
  x = setMarkers(x, locusAttributes=ann)
  expect_equal(nMarkers(x), 2)
  expect_true(all(do.call(cbind, x$MARKERS) == 0))
  expect_equal(name(x, 1), "snp1")
  expect_equal(alleles(x, 2), c("a", "b"))
})

test_that("addMarkers() adds empty markers with locusAttributes", {
  x = nuclearPed(1) |> addMarker(name = "M0")

  ann = list(list(alleles=1:2, name="snp1"), list(alleles=c("a", "b")))
  x = addMarkers(x, locusAttributes=ann)
  expect_equal(nMarkers(x), 3)
  expect_equal(name(x, 2), "snp1")
  expect_equal(alleles(x, 3), c("a", "b"))
})


test_that("setMarkers() gives same result with and without alleleMarkers", {
  x = nuclearPed(fa="fa", mo="mo", children="boy")
  m = marker(x, mo='b', boy=c('b','c'), alleles=letters[1:3], name="snp")
  y = setMarkers(x, m)

  amat = cbind(c(0,'b','b'),c(0,'b','c'))
  ann = list(alleles=letters[1:3], name="snp")
  z = setMarkers(x, alleleMatrix = amat, locusAttributes=ann)

  expect_identical(y,z)
})


test_that("addMarkers() gives same result with and without alleleMarkers", {
  x = nuclearPed(fa="fa", mo="mo", children="boy") |>
    addMarker(name="M0")
  y = addMarker(x, mo='b/b', boy='b/c', alleles=letters[1:3], name="snp")

  amat = cbind(c(0,'b','b'),c(0,'b','c'))
  ann = list(alleles=letters[1:3], name="snp")
  z = addMarkers(x, alleleMatrix = amat, locusAttributes=ann)

  expect_identical(y,z)
})

test_that("setMarkers() deals with marker names and ordering", {
  ### singleton
  s = singleton("a")
  m1 = marker(s, a = 1:2, name = "M1")
  m2 = marker(s, a = 2:3, name = "M2")
  s12 = setMarkers(s, list(m1, m2))

  am = getAlleles(s12)
  locAttr = getLocusAttributes(s12)

  expect_identical(s12, setMarkers(s, alleleMatrix = am, locus = locAttr))
  expect_identical(s12, setMarkers(s, alleleMatrix = am, locus = locAttr[2:1]))

  am2 = am[, 3:4, drop = F]
  expect_identical(setMarkers(s, m2),
                   setMarkers(s, alleleMatrix = am2, locusAttributes = locAttr))
  expect_error(setMarkers(s, alleleMatrix = am2, locusAttributes = locAttr[1]),
               "Marker name found in `allelematrix`, but not in `locusAttributes`: M2")

  ### trio
  trio = nuclearPed(child = "a")
  m1 = marker(trio, a = "1/2", name = "M1")
  m2 = marker(trio, a = "2/3", name = "M2")
  trio12 = setMarkers(trio, list(m1, m2))

  am_all = getAlleles(trio12)
  am_child = getAlleles(trio12, ids = "a")
  locAttr = getLocusAttributes(trio12)

  expect_identical(trio12, setMarkers(trio, alleleMatrix = am_all, locus = locAttr))
  expect_identical(trio12, setMarkers(trio, alleleMatrix = am_child, locus = locAttr))

  expect_identical(trio12, setMarkers(trio, alleleMatrix = am_all, locus = locAttr[2:1]))
  expect_identical(trio12, setMarkers(trio, alleleMatrix = am_child, locus = locAttr[2:1]))

  am2_all = am_all[, 3:4, drop = F]
  am2_child = am_child[, 3:4, drop = F]
  expect_identical(setMarkers(trio, m2),
                   setMarkers(trio, alleleMatrix = am2_all, locusAttributes = locAttr))
  expect_identical(setMarkers(trio, m2),
                   setMarkers(trio, alleleMatrix = am2_child, locusAttributes = locAttr))

  expect_error(setMarkers(trio, alleleMatrix = am2_all, locusAttributes = locAttr[1]),
               "Marker name found in `allelematrix`, but not in `locusAttributes`: M2")
})

test_that("setMarkers() gives correct errors with duplicated marker names", {
  x = singleton(1)
  m1 = marker(x, name = "M")
  expect_error(setMarkers(x, list(m1,m1)), "Duplicated marker name: M")

  m2 = marker(x)
  expect_silent(setMarkers(x, list(m2,m2)))

  x = setMarkers(x, list(m1, m2))
  expect_error({name(x, 2) <- "M" }, "Duplicated marker name: M")

})
