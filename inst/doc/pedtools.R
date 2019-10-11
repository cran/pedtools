## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(kableExtra)

## ---- eval = FALSE-------------------------------------------------------
#   # First install devtools if needed
#  if(!require(devtools)) install.packages("devtools")
#  
#  # Install pedtools from github
#  devtools::install_github("magnusdv/pedtools")

## ----message=FALSE-------------------------------------------------------
library(pedtools)

## ------------------------------------------------------------------------
ped(id = 1:3, fid = c(0,0,1), mid = c(0,0,2), sex = c(1,2,2))

## ------------------------------------------------------------------------
trio = ped(id = c("fa", "mo", "girl"), fid = c("","","fa"), mid = c("","","mo"), sex = c(1,2,2))
trio

## ------------------------------------------------------------------------
unclass(trio)

## ----eval=FALSE----------------------------------------------------------
#  plot(trio)

## ----echo=FALSE, fig.dim = c(2,2)----------------------------------------
plot(trio, margins = c(1,1,1,1))

## ----eval=FALSE----------------------------------------------------------
#  plot(trio, deceased = "fa", starred = "mo", shaded = "girl",
#       col = c("green", "red", "blue"), title = "Trio 1")

## ----echo=FALSE, fig.dim = c(2,2)----------------------------------------
plot(trio, deceased = "fa", starred = "mo", shaded = "girl", 
     col = c("green", "red", "blue"), title = "Trio 1", margins = c(1,1,1.5,1))

## ------------------------------------------------------------------------
trio2 = nuclearPed(nch = 1)
trio2 = swapSex(trio2, id = 3)
trio2 = relabel(trio2, new = c("fa", "mo", "girl"))

## ------------------------------------------------------------------------
trio3 = nuclearPed(father = "fa", mother = "mo", children = "girl", sex = 2)

## ------------------------------------------------------------------------
trio4 = singleton("fa")
trio4 = addDaughter(trio4, parent = "fa", id = "girl")
trio4 = relabel(trio4, old = "NN_1", new = "mo")

## ----echo=FALSE----------------------------------------------------------
x1 = halfSibPed(nch1 = 1, nch2 = 2, sex1 = 1, sex2 = 2:1)
x1 = addChildren(x1, father = 4, mother = 5, nch = 1)
plot(x1, margin = c(1,1,1,1))

## ------------------------------------------------------------------------
x1 = halfSibPed(nch1 = 1, nch2 = 2, sex1 = 1, sex2 = 2:1)
x1 = addChildren(x1, father = 4, mother = 5, nch = 1)

## ------------------------------------------------------------------------
x2 = halfCousinPed(0, child = T)
x2 = addChildren(x2, father = 1, mother = 4, nch = 1)
x2 = relabel(x2, old = c(4,3,7,6), new = c(3,4,6,7))

## ------------------------------------------------------------------------
identical(x1, x2)

## ------------------------------------------------------------------------
x2

## ------------------------------------------------------------------------
x2 = reorderPed(x2)
identical(x1, x2)

## ----echo = FALSE, fig.width = 2.5, fig.width = 4------------------------
y1 = linearPed(3, sex = 1) # male line 
y2 = linearPed(2, sex = 2) # female line

y2 = relabel(y2, c(11:14, 6)) # relabel to match only for ID = 6

y = mergePed(y1, y2)
plot(y, id=NULL, margins=c(0,1,1,1), symbolsize=0.9)

## ----eval=FALSE----------------------------------------------------------
#  y1 = linearPed(3, sex = 1) # male line
#  y2 = linearPed(2, sex = 2) # female line
#  
#  y2 = relabel(y2, c(11:14, 6)) # relabel to match at ID = 6
#  
#  y = mergePed(y1, y2)
#  plot(y)

## ------------------------------------------------------------------------
marker(trio)

## ------------------------------------------------------------------------
m1 = marker(trio, fa = "A", mo = c("A","B"), name = "snp1")

## ------------------------------------------------------------------------
marker(trio, fa = "A", mo = c("A","B"), alleles = c("A","B","C"), afreq = c(.2,.3,.5))

## ------------------------------------------------------------------------
m2 = marker(trio, fa = "A", mo = c("A","B"), chrom = 23, name = "snpX")

## ---- eval=FALSE---------------------------------------------------------
#  plot(trio, marker = m1)

## ---- echo=FALSE, fig.dim=c(2,2)-----------------------------------------
plot(trio, marker = m1, margins = c(1,1,1,1))

## ---- eval=FALSE---------------------------------------------------------
#  plot(trio, marker = m1, sep = "", skip.empty.genotypes = T)

## ---- echo=FALSE, fig.dim=c(2,2)-----------------------------------------
plot(trio, marker = m1, sep = "", skip.empty.genotypes = T, margins = c(1,1,1,1))

## ------------------------------------------------------------------------
trio = setMarkers(trio, list(m1, m2))
trio

## ------------------------------------------------------------------------
whichMarkers(trio, chrom = 23)
selectMarkers(trio, markers = "snp1")

## ------------------------------------------------------------------------
afreq(m1)
afreq(trio, marker = "snp1")

## ------------------------------------------------------------------------
afreq(trio, marker = "snp1") = c(A = 0.9, B = 0.1)

## ------------------------------------------------------------------------
genotype(trio, "snpX", id = "girl")
genotype(trio, "snpX", id = "girl") = "A"
trio

## ----getset, echo = FALSE------------------------------------------------
getters.df = rbind(
  c("`getAlleles(x)`", 
    "extract all alleles as a matrix.", 
    "do summary stats on the marker alleles"),
  c("`getFrequencyDatabase(x)`", 
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
  c("`setFrequencyDatabase(x, db)`", 
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
tbl.getset = column_spec(tbl.getset, 1, width = "5cm")
tbl.getset = column_spec(tbl.getset, 2, width = "8cm")
tbl.getset = pack_rows(tbl.getset, "Get", 1, 3, indent = F)
tbl.getset = pack_rows(tbl.getset, "Set", 4, 6, indent = F)
tbl.getset = pack_rows(tbl.getset, "Convert", 7, 8, indent = F)
tbl.getset = pack_rows(tbl.getset, "Transfer", 9, 9, indent = F)
tbl.getset

