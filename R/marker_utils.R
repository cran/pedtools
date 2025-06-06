

# Check if a mutation matrix - or a list of such - is diagonal (= trivial)
trivialMut = function(mut) {
  if(is.list(mut))
    return(all(vapply(mut, trivialMut, logical(1))))

  all(diag(mut) == 1)
}

#' Test if something is a marker
#'
#' Functions for testing if something is a `marker` object, or a list of such
#' objects.
#'
#' @param x Any object
#'
#' @return A logical
#' @export
is.marker = function(x) {
  inherits(x, "marker")
}

#' @export
#' @rdname is.marker
is.markerList = function(x) {
  inherits(x, "markerList") || (is.list(x) && all(sapply(x, is.marker)))
}


#' The number of markers attached to a pedigree
#'
#' @param x A `ped` object or a list of such.
#' @param compwise A logical, only relevant if `x` is a ped list. Default FALSE.
#'
#' @return `nMarkers()` by default returns a single number; the number of marker
#'   objects attached to `x`. If `x` is a ped list, an error is raised if the
#'   components have different numbers of markers. This check can be skipped by
#'   setting `compwise = TRUE`, in which case the function returns a vector of
#'   the component-wise marker numbers.
#'
#'   The function `hasMarkers(x)` returns TRUE if (at least component of) `x`
#'   has attached markers, otherwise FALSE. If `compwise = TRUE`, a logical
#'   vector of the same length as `x`.
#'
#' @examples
#' x = nuclearPed() |> addMarker()
#' nMarkers(x) # = 1
#'
#' y = list(x, singleton(1))
#' nMarkers(y, compwise = TRUE) # c(1,0)
#'
#' hasMarkers(y) # TRUE
#' hasMarkers(y, compwise = TRUE) # c(TRUE, FALSE)
#'
#' @export
nMarkers = function(x, compwise = FALSE) {
  if(is.ped(x))
    return(length(x$MARKERS))
  else if(is.pedList(x)) {
    nvec = vapply(x, function(comp) length(comp$MARKERS), 1L)

    if(compwise)
      return(nvec)

    if(!listIdentical(nvec))
      stop2("The pedigree components have different number of markers attached")
    return(nvec[[1]])
  }
  stop2("Input to `nMarkers()` must be a `ped` object or a list of such")
}

#' @export
#' @rdname nMarkers
hasMarkers = function(x, compwise = FALSE) {
  nm = nMarkers(x, compwise = TRUE)
  if(compwise) nm > 0 else any(nm > 0)
}

checkDupNames = function(x) {
  if(is.ped(x))
    mlist = x$MARKERS
  else if(is.pedList(x))
    mlist = x[[1]]$MARKERS
  else
    mlist = x

  mnames = unlist(lapply(mlist, attr, "name"))
  if(dup <- anyDuplicated.default(mnames, incomparables = NA))
    stop2("Duplicated marker name: ", mnames[dup])
}

checkConsistency = function(x, mlist) {
  checkDupNames(mlist)

  wrongSize = unlist(lapply(mlist, nrow) != pedsize(x))
  if(any(wrongSize)) {
    erri = which(wrongSize)[1]
    errsize = nrow(mlist[[erri]])
    stop("Incompatible input: Pedigree has size ", pedsize(x),
         " but marker ", erri, " has ", errsize, " rows", call. = FALSE)
  }
  #TODO: check loop breakers, ped labels, sex
  return(TRUE)
}


#' Add allele
#'
#' Extends the allele set of a marker attached to a pedigree, by adding a single
#' allele.
#'
#' @param x A `ped` object or a list of such, or a frequency database (list of
#'   numeric vectors).
#' @param marker The name or index of a marker attached to `x`.
#' @param allele The name of the new allele.
#' @param freq The frequency of the new allele, by default 0.001.
#' @param adjust Either "previous" or "all", indicating how the frequencies
#'   should be adjusted so that they sum to 1. If "previous" (default), the
#'   frequencies of the original alleles are multiplied with `1 - freq`. If
#'   "all", scaling is performed after adding the new allele, i.e., dividing all
#'   frequencies by `1 + freq`.
#'
#' @return A copy of `x` with modified marker attributes.
#'
#' @examples
#'
#' ## Ped input
#' x = nuclearPed() |>
#'   addMarker(geno = c(NA, NA, "b/c"), afreq = c(b = 0.5, c = 0.5))
#'
#' y = addAllele(x, marker = 1, allele = "a")
#' afreq(y, 1)
#'
#' z = addAllele(y, marker = 1, allele = "d", freq = 0.1, adjust = "all")
#' afreq(z, 1)
#'
#'
#' ## Database input
#' db = list(M1 = c(a = .2, b = .3, c = .5),
#'           M2 = c("7" = .9, "8.3" = .1))
#' addAllele(db, marker = "M2", allele = "8")
#'
#' @export
addAllele = function(x, marker, allele, freq = 0.001, adjust = c("previous", "all")) {

  if(missing(marker) || length(marker) == 0)
    stop2("Argument `marker` cannot be empty")
  if(length(marker) > 1)
    stop2("Allele modifications must be done for a single marker at a time")
  if(missing(allele) || length(allele) == 0)
    stop2("Argument `allele` cannot be empty")
  if(length(allele) > 1)
    stop2("Please add one allele at a time")

  # Input type: ped or db?
  type = if(is.ped(x) || is.pedList(x)) "ped" else "db"

  if(type == "ped") {
    old = afreq(x, marker)
  }
  else {
    old = x[[marker]]
    if(is.null(old) || is.null(names(old)) || !is.numeric(old))
      stop2("Something wrong with the existing frequency vector: Not a named numeric.")
    if(round(sum(old),3) != 1)
      stop2("Something wrong with the existing frequency vector: Does not sum to 1:\n", old)
  }

  if(allele %in% names(old))
    stop2("This allele is already present")

  names(freq) = allele
  s = sum(old) # should be 1, but may deviate slightly

  # Add new allele and adjust frequencies
  switch(match.arg(adjust),
    all = {
      newfr = c(old, freq)/(s + freq)
    },
    previous = {
      newfr = c(old/s * (1 - freq), freq) # note normalisation of `old`
    }
  )

  if(type == "ped") {
    x = setAfreq(x, marker, afreq = newfr, strict = FALSE)
  }
  else {
    als = names(newfr)
    isNum = !anyNA(suppressWarnings(as.numeric(als)))
    ord = if(isNum) order(as.numeric(als)) else order(als)
    x[[marker]] = newfr[ord]
  }

  x
}


#' Swap genotypes between individuals
#'
#' @param x A `ped` object or a list of such.
#' @param ids A vector of 2 members of `x`.
#'
#' @return An object identical to `x`, except that the genotypes of the `ids`
#'   pair have been swapped.
#'
#' @seealso [transferMarkers()]
#'
#' @examples
#' x = nuclearPed() |>
#'   addMarker(geno = c("1/1", "2/2", "3/3"))
#'
#' swapGenotypes(x, ids = 1:2)
#'
#' @export
swapGenotypes = function(x, ids = NULL) {
  if(length(ids) < 2)
    stop2("Argument `ids` must have length at least 2")

  ord = order(match(ids, labels(x, unlist = TRUE), nomatch = 0L))

  # Main case: swap a pair
  if(length(ids) == 2 && ord[1] < ord[2]) {
    ids = ids[2:1]
    ord = ord[2:1]
  }

  als = getAlleles(x, ids)

  # Swap/reorder rownames
  rownames(als) = rownames(als)[ord]

  setAlleles(x, ids = ids, alleles = als)
}

#' Harmonise markers across components in a ped list
#'
#' This function ensures that all components of a ped list have the same
#' markers, in the same order. Missing markers are added with empty genotypes.
#' Note that _unnamed_ markers cannot be harmonised and will be removed by this
#' function.
#'
#' @param x A list of `ped` objects. An error is raised if any component
#'   contains an unnamed marker.
#'
#' @returns A copy of `x` where all components have the same markers attached,
#'   and in the same order. Unnamed markers are removed.
#'
#' @examples
#' x = list(
#'   singleton(1) |> addMarker(),  # unnamed marker will be removed
#'   singleton(2) |> addMarker(name = "M1"),
#'   singleton(3) |> addMarker(geno = "3/3", alleles = 1:3, name = "M2")
#' )
#' harmoniseMarkers(x)
#'
#' @export
harmoniseMarkers = function(x) {
  # If connected pedigree, return as is
  if(is.ped(x) || (length(x) == 1 && is.ped(x[[1]])))
    return(x)

  mList0 = lapply(x, name)
  mList = lapply(mList0, function(v) v[!is.na(v)])
  if(listIdentical(mList)) {
    if(identical(mList0, mList))
      return(x)
    else
      return(selectMarkers(x, markers = mList[[1]]))
  }

  firstCmp = .firstComp(mList) # first component of each
  allMarkers = names(firstCmp)

  # List of locus attributes for all markers
  allAttribs = list()
  for(i in 1:max(firstCmp)) {
    new = getLocusAttributes(x[[i]], markers = allMarkers[firstCmp == i])
    allAttribs = c(allAttribs, new)
  }
  allAttribs

  # Add missing markers to each component
  for(i in seq_along(x)) {
    missingMs = setdiff(allMarkers, mList[[i]])
    if(length(missingMs))
      x[[i]] = addMarkers(x[[i]], locusAttributes = allAttribs[missingMs])
  }

  # Ensure same order
  selectMarkers(x, markers = allMarkers)
}

# Utility: In a list of overlapping characters, find first component of each element
.firstComp = function(lst) {
  allElements = unlist(lst, recursive = FALSE, use.names = FALSE)
  u = !duplicated(allElements)
  lstIdx = rep(seq_along(lst), lengths(lst))
  .setnames(lstIdx[u], allElements[u])
}
