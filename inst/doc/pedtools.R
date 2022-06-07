## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center"
)
library(kableExtra)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("pedtools")

## ---- eval = FALSE------------------------------------------------------------
#  # install.packages("devtools") # install devtools if needed
#  devtools::install_github("magnusdv/pedtools")

## ----message=FALSE------------------------------------------------------------
library(pedtools)

## -----------------------------------------------------------------------------
ped(id = 1:3, fid = c(0,0,1), mid = c(0,0,2), sex = c(1,2,2))

## -----------------------------------------------------------------------------
trio = ped(id = c("fa", "mo", "girl"), fid = c("","","fa"), mid = c("","","mo"), sex = c(1,2,2))
trio

## -----------------------------------------------------------------------------
unclass(trio)

## ----eval=FALSE---------------------------------------------------------------
#  plot(trio)

## ----echo=FALSE, fig.dim = c(2,2)---------------------------------------------
plot(trio, margins = c(1,1,1,1))

## ----eval=FALSE---------------------------------------------------------------
#  plot(trio, deceased = "fa", starred = "mo", hatched = "girl",
#       col = c("green", "red", "blue"), title = "Trio 1")

## ----echo=FALSE, fig.dim = c(2,2)---------------------------------------------
plot(trio, deceased = "fa", starred = "mo", hatched = "girl", 
     col = c("green", "red", "blue"), title = "Trio 1", margins = c(1,1,1.5,1))

## -----------------------------------------------------------------------------
trio2 = nuclearPed(nch = 1)
trio2 = swapSex(trio2, ids = 3)
trio2 = relabel(trio2, new = c("fa", "mo", "girl"))

## -----------------------------------------------------------------------------
trio3 = nuclearPed(father = "fa", mother = "mo", children = "girl", sex = 2)

## -----------------------------------------------------------------------------
trio4 = singleton("fa")
trio4 = addDaughter(trio4, parent = "fa", id = "girl")
trio4 = relabel(trio4, old = "NN_1", new = "mo")

## ----echo=FALSE---------------------------------------------------------------
x1 = halfSibPed(nch1 = 1, nch2 = 2, sex1 = 1, sex2 = 2:1)
x1 = addChildren(x1, father = 4, mother = 5, nch = 1)
plot(x1, margin = c(1,1,1,1))

## -----------------------------------------------------------------------------
x1 = halfSibPed(nch1 = 1, nch2 = 2, sex1 = 1, sex2 = 2:1)
x1 = addSon(x1, parents = 4:5)

## -----------------------------------------------------------------------------
x2 = halfCousinPed(0, child = T)
x2 = addSon(x2, parents = 2:3)
x2 = relabel(x2, old = c(7,6), new = c(6,7))

## -----------------------------------------------------------------------------
identical(x1, x2)

## -----------------------------------------------------------------------------
x2

## -----------------------------------------------------------------------------
x2 = reorderPed(x2)
identical(x1, x2)

## ----merge-example, echo = FALSE, message=FALSE-------------------------------
# Top part
x = ancestralPed(g = 2) # bottom person is "7"

# Bottom part
y = halfCousinPed(degree = 1) # top person is "2"
y = swapSex(y, 4)

# Merge
z = mergePed(x, y, by = c("7" = "2"), relabel = TRUE)

## ----merge-plot, echo = FALSE, fig.width = 3.5, fig.height = 3.7--------------
plot(z, margins = c(1,1,1,1))

## ---- label = "merge-example"-------------------------------------------------
# Top part
x = ancestralPed(g = 2) # bottom person is "7"

# Bottom part
y = halfCousinPed(degree = 1) # top person is "2"
y = swapSex(y, 4)

# Merge
z = mergePed(x, y, by = c("7" = "2"), relabel = TRUE)

## ----merge-parts, fig.width = 9, fig.height = 3.5-----------------------------
plotPedList(list(x, y, z))

## -----------------------------------------------------------------------------
marker(trio)

## -----------------------------------------------------------------------------
m1 = marker(trio, fa = "A/A", mo = "A/B", name = "snp1")

## -----------------------------------------------------------------------------
marker(trio, fa = "A/A", mo = "A/B", afreq = c(A = .2, B = .3, C = .5))

## -----------------------------------------------------------------------------
m2 = marker(trio, fa = "A/A", mo = "A/B", chrom = "X", name = "snpX")
m2

## ---- eval=FALSE--------------------------------------------------------------
#  plot(trio, marker = m1)

## ---- echo=FALSE, fig.dim=c(2,2)----------------------------------------------
plot(trio, marker = m1, margins = c(1,1,1,1))

## ---- eval=FALSE--------------------------------------------------------------
#  plot(trio, marker = m1, showEmpty = T, missing = "?", sep = "")

## ---- echo=FALSE, fig.dim=c(2,2)----------------------------------------------
plot(trio, marker = m1, sep = "", showEmpty = T, missing = "?", margins = c(1,1,1,1))

## -----------------------------------------------------------------------------
trio = setMarkers(trio, list(m1, m2))
trio

## -----------------------------------------------------------------------------
nuclearPed(1) |> 
  addMarker(name = "myMarker", alleles = c("a", "b", "c")) |>
  setGenotype(marker = 1, id = 3, geno = "a/c")

## -----------------------------------------------------------------------------
whichMarkers(trio, chrom = "X")
selectMarkers(trio, markers = "snp1")

## -----------------------------------------------------------------------------
afreq(m1)
afreq(trio, marker = "snp1")

## -----------------------------------------------------------------------------
afreq(trio, marker = "snp1") = c(A = 0.9, B = 0.1)

## -----------------------------------------------------------------------------
# Girl is not genotyped
genotype(trio, "snpX", id = "girl")

# Set genotype
genotype(trio, "snpX", id = "girl") = "A/A"

# Inspect the result
trio

## ----getset, echo = FALSE-----------------------------------------------------
getters.df = rbind(
  c("`getAlleles(x)`", 
    "extract all alleles as a matrix.", 
    "do summary stats on the marker alleles"),
  c("`getFreqDatabase(x)`", 
    "extract allele frequencies as a data.frame in *allelic ladder* format.", 
    "transfer to other objects, or write the database to a file"),
  c("`getMarkers(x)`", 
    "extract list of marker objects. Each marker is a `N * 2` allele matrix (`N = pedsize(x)`) with locus annotations as attributes", 
    "do computations")
)

setters.df = rbind(
  c("`setAlleles(x, ...)`", 
    "replace the genotypes of `x` without changing the locus attributes.", 
    "erase all genotypes"),
  c("`setFreqDatabase(x, db)`", 
    "replace all allele frequencies without changing the genotype data. The input is a data.frame in *allelic ladder* format. Conceptually equivalent to `setMarkers(x, alleleMatrix = getAlleles(x), locusAnnotations = db)`.", 
    "change the frequency database"),
  c("`setMarkers(x, ...)`", 
    "attach marker objects with or without genotype data. Locus attributes are indicated as a list; genotypes as a matrix or data.frame.", 
    "prepare joint manipulation of a pedigree and marker data")
)

conversions.df = rbind(
  c("`as.data.frame(x)`", 
    "convert `x` to a data.frame, with pedigree columns in standard format followed by genotype columns. One column per marker, with genotype format `a/b` and missing alleles indicated as `-`.", 
    "pretty-print ped objects"),
  c("`as.matrix(x)`", 
    "convert `x` to a *numerical* matrix, with additional info attached as attributes.", 
    "modify a pedigree with marker data")
)

other.df = rbind(
  c("`transferMarkers(from, to)`", 
    "transfer genotypes and attributes between pedigree objects (or lists of such).", 
    "transfer simulated marker data")
)

getset.df = rbind(getters.df, setters.df, conversions.df, other.df)
tbl.getset = kable(getset.df, 
            col.names = c("Use ...", "When you want to ...", "For example to ..."))
tbl.getset = column_spec(tbl.getset, 1, width = "4.5cm")
tbl.getset = column_spec(tbl.getset, 2, width = "8.5cm")
tbl.getset = pack_rows(tbl.getset, "Get", 1, 3, indent = F)
tbl.getset = pack_rows(tbl.getset, "Set", 4, 6, indent = F)
tbl.getset = pack_rows(tbl.getset, "Convert", 7, 8, indent = F)
tbl.getset = pack_rows(tbl.getset, "Transfer", 9, 9, indent = F)
tbl.getset

