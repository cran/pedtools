#' Add/remove pedigree members
#'
#' Functions for adding or removing individuals in a 'ped' object.
#'
#' In `addChildren()` and `addParents()`, labels of added individuals are
#' created automatically if they are not specified by the user. In the automatic
#' case, the labelling depends on whether the existing labels are integer-like
#' or not (i.e. if `labels(x)` equals `as.character(as.integer(labels(x)))`.) If
#' so, the new labels are integers subsequent to the largest of the existing
#' labels. If not, the new labels are "NN_1", "NN_2", ... If any such label
#' already exists, the numbers are adjusted accordingly.
#'
#' `addSon()` and `addDaughter()` are wrappers for the most common use of
#' `addChildren()`, namely adding a single child to a pedigree. Note that the
#' parents can be given in any order. If only one parent is supplied, the other
#' is created as a new individual.
#'
#' `removeIndividuals()` removes the individuals indicated with `ids` along with
#' all of their ancestors OR descendants, depending on the `remove` argument.
#' Leftover spouses disconnected to the remaining pedigree are also removed. An
#' error is raised if result is a disconnected pedigree.
#'
#' The `branch()` function extracts the sub-pedigree formed by `id` and all
#' his/her spouses and descendants.
#'
#' Finally, `subset()` can be used to extract any connected sub-pedigree. (Note
#' that in the current implementation, the function does not actually check that
#' the indicated subset forms a connected pedigree; failing to comply with this
#' may lead to obscure errors.)
#'
#' @param x A `ped` object.
#' @param ids A character vector (or coercible to such) with ID labels. In
#'   `addChildren` the (optional) `ids` argument is used to specify labels for
#'   the created children. If given, its length must equal `nch`. If not given,
#'   labels are assigned automatically as explained in Details.
#' @param remove Either "ancestors" or "descendants" (default), dictating the
#'   method of removing pedigree members. Abbreviations are allowed.
#' @param returnLabs A logical, by default FALSE. If TRUE, `removeIndividuals()`
#'   returns only the labels of all members to be removed, instead of actually
#'   removing them.
#' @param id The ID label of a pedigree member.
#' @param father,mother Single ID labels. At least one of these must belong to
#'   an existing pedigree member. The other label may either: 1) belong to an
#'   existing member, 2) not belong to any existing member, or 3) be missing
#'   (i.e. not included in the function call). In cases 2 and 3 a new founder is
#'   added to the pedigree. In case 2 its label is the one given, while in case
#'   3 a suitable label is created by the program (see Details).
#' @param nch A positive integer indicating the number of children to be
#'   created. Default: 1.
#' @param sex Gender codes of the created children (recycled if needed).
#' @param verbose A logical: Verbose output or not.
#' @param parents A vector of 1 or 2 ID labels, of which at least one must be an
#'   existing member of `x`.
#'
#' @return The modified `ped` object.
#' @seealso [ped()], [relabel()], [swapSex()]
#'
#' @examples
#'
#' x = nuclearPed(1) |>
#'   addSon(3) |>
#'   addParents(4, father = 6, mother = 7) |>
#'   addChildren(father = 6, mother = 7, nch = 3, sex = c(2,1,2))
#'
#' # Remove 6 and 7 and their descendants
#' y1 = removeIndividuals(x, 6:7)
#'
#' # Remove 8-10 and their parents
#' y2 = removeIndividuals(x, 8:10, remove = "ancestors")
#'
#' @name ped_modify
NULL

#' @rdname ped_modify
#' @export
addChildren = function(x, father = NULL, mother = NULL, nch = NULL, sex = 1, ids = NULL, verbose = TRUE) {
  if(!is.ped(x) && !is.pedList(x))
    stop2("Input is not a `ped` object or a list of such")

  nch = nch %||% if(!is.null(ids)) length(ids) else length(sex)
  if(!isCount(nch))
    stop2("Argument `nch` must be a positive integer: ", nch)

  # This variable will change as new members are created
  labs = labels(x)

  # Check input
  father_exists = isTRUE(father %in% labs)
  mother_exists = isTRUE(mother %in% labs)

  if(!father_exists && !mother_exists)
    stop2("At least one parent must be an existing pedigree member")

  if(father_exists && getSex(x, father) == 2)
    stop2("Assigned father is female: ", father)
  if(mother_exists && getSex(x, mother) == 1)
    stop2("Assigned mother is male: ", mother)

  if(!is.null(ids) && length(ids) != nch)
    stop2("Length of 'ids' must equal the number of children")
  if(any(ids %in% labs))
    stop2("Individuals already exist: ", intersect(ids, labs))

  # Recycle `sex` if needed
  sex = rep_len(sex, nch)

  if(is.pedList(x)) {
    comp = getComponent(x, c(father, mother), checkUnique = TRUE)
    compSet = unique.default(comp[!is.na(comp)])
    if(length(compSet) == 2)
      y = .addChildrenAcrossComps(x, father, mother, nch = nch, sex = sex, ids = ids, verbose = verbose)
    else if(length(compSet) == 1) {
      co = x[[compSet]]
      y = x
      y[[compSet]] = addChildren(co, father, mother, nch = nch, sex = sex, ids = ids, verbose = verbose)
    }
    return(y)
  }

  n = pedsize(x)

  # Prepare manipulation of pedigree matrix
  p = as.matrix(x)
  attrs = attributes(p)
  nmark = nMarkers(x)

  # Utility for creating new labels, depending on existing labels being numeric
  nextlabs = function(labs, len) {
    if(hasNumLabs(x)) {
      mx = max(as.numeric(labs))
      seq.int(mx + 1, length.out = len)
    }
    else {
      res = character(0)
      for(i in seq_len(len)) res = c(res, nextNN(c(labs, res)))
      res
    }
  }

  if(!father_exists) {
    if(is.null(father))
      father = nextlabs(labs, 1)
    if (verbose)
      message("Father: Creating new individual with ID = ", father)

    labs = c(labs, father)
    father_int = n = n + 1
    p = rbind(p, c(father_int, 0, 0, 1, rep(0, 2*nmark)))
  }
  else {
    father_int = internalID(x, father)
  }
  if(!mother_exists) {
    if(is.null(mother))
      mother = nextlabs(labs, 1)
    if (verbose)
      message("Mother: Creating new individual with ID = ", mother)

    labs = c(labs, mother)
    mother_int = n = n + 1
    p = rbind(p, c(mother_int, 0, 0, 2, rep(0, 2*nmark)))
  }
  else {
    mother_int = internalID(x, mother)
  }
  # Children
  if(is.null(ids)) ids = nextlabs(labs, len = nch)
  labs = c(labs, ids)

  if (anyDuplicated.default(labs))
    stop2("Duplicated ID label: ", labs[duplicated(labs)])

  children_pedcols = cbind(nrow(p) + (1:nch), father_int, mother_int, sex)
  children_markers = matrix(0, nrow = nch, ncol = nmark*2)
  p = rbind(p, cbind(children_pedcols, children_markers))

  attrs$LABELS = as.character(labs)
  restorePed(p, attrs = attrs)
}

.addChildrenAcrossComps = function(x, father, mother, nch, sex, ids, verbose) {
  stop2("Adding children across components is not implemented yet")
}

#' @rdname ped_modify
#' @export
addSon = function(x, parents, id = NULL, verbose = TRUE) {

  npar = length(parents)
  if(npar == 0 || npar > 2)
    stop2("Argument `parents` must have length 1 or 2: ", parents)

  parents = as.character(parents) # remove potential names etc
  existing = parents %in% unlist(labels(x))
  if(!any(existing))
    stop2("At least one parent must be an existing pedigree member: ", parents)

  par1 = if(existing[1]) parents[1] else parents[2]
  par2 = if(npar == 2) setdiff(parents, par1) else NULL

  sex1 = getSex(x, par1)
  if(sex1 == 1)
    addChildren(x, father = par1, mother = par2, nch = 1, sex = 1, ids = id, verbose = verbose)
  else if(sex1 == 2)
    addChildren(x, father = par2, mother = par1, nch = 1, sex = 1, ids = id, verbose = verbose)
  else
    stop2("Not implemented for parents of unknown sex: ", par1)
}

#' @rdname ped_modify
#' @export
addDaughter = function(x, parents, id = NULL, verbose = TRUE) {

  npar = length(parents)
  if(npar == 0 || npar > 2)
    stop2("Argument `parents` must have length 1 or 2: ", parents)

  parents = as.character(parents) # remove potential names etc
  existing = parents %in% unlist(labels(x))
  if(!any(existing))
    stop2("At least one parent must be an existing pedigree member: ", parents)

  par1 = if(existing[1]) parents[1] else parents[2]
  par2 = if(npar == 2) setdiff(parents, par1) else NULL

  sex1 = getSex(x, par1)
  if(sex1 == 1)
    addChildren(x, father = par1, mother = par2, nch = 1, sex = 2, ids = id, verbose = verbose)
  else if(sex1 == 2)
    addChildren(x, father = par2, mother = par1, nch = 1, sex = 2, ids = id, verbose = verbose)
  else
    stop2("Not implemented for parents of unknown sex: ", par1)
}

#' @rdname ped_modify
#' @export
addParents = function(x, id, father = NULL, mother = NULL, verbose = TRUE) {
  if (length(id) > 1)
      stop2("Parents cannot be added to multiple individuals at once: ", id)
  if (id %in% nonfounders(x))
    stop2("Individual ", id, " already has parents in the pedigree")

  id_int = internalID(x, id)
  labs = labels(x)

  # Check that assigned parents are OK
  desc = descendants(x, id)
  if (!is.null(father)) {
    if (father == id) stop2("Father and child are equal")
    if (father %in% desc) stop2("Assigned father is a descendant of ", id)
    if (father %in% labs && getSex(x, father) == 2) stop2("Assigned father is female")
  }
  if (!is.null(mother)) {
    if (mother == id) stop2("Mother and child are equal")
    if (mother %in% desc) stop2("Assigned mother is a descendant of ", id)
    if (mother %in% labs && getSex(x, mother) == 1) stop2("Assigned mother is male")
  }

  # If no labels are given, create them
  if (is.null(father) || is.null(mother)) {
    if(hasNumLabs(x)) {
      mx = max(as.numeric(labs))
      if (is.null(father)) {father = mx + 1; mx = mx + 1}
      if (is.null(mother)) {mother = mx + 1}
    }
    else {
      if (is.null(father)) {father = nextNN(labs)}
      if (is.null(mother)) {mother = nextNN(c(labs, father))}
    }
  }

  p = as.matrix(x)
  attrs = attributes(p)
  newlabs = attrs$LABELS
  nmark = nMarkers(x)

  new.father = !father %in% labs
  new.mother = !mother %in% labs

  if (new.father) {
    if(verbose) message("Father: Creating new individual with ID = ", father)
    fath_int = nrow(p) + 1
    p = rbind(p, c(fath_int, 0, 0, 1, rep(0, 2*nmark)))
    newlabs = c(newlabs, father)
  }
  else fath_int = internalID(x, father)

  if (new.mother) {
    if(verbose) message("Mother: Creating new individual with ID = ", mother)
    moth_int = nrow(p) + 1
    p = rbind(p, c(moth_int, 0, 0, 2, rep(0, 2*nmark)))
    newlabs = c(newlabs, mother)
  }
  else moth_int = internalID(x, mother)

  p[id_int, 2:3] = c(fath_int, moth_int)
  attrs$LABELS = newlabs

  y = restorePed(p, attrs = attrs)
  neworder = c(seq_len(id_int - 1),
               if(new.father) fath_int,
               if(new.mother) moth_int,
               id_int:pedsize(x))
  reorderPed(y, neworder)
}



#' @rdname ped_modify
#' @export
removeIndividuals = function(x, ids, remove = c("descendants", "ancestors"),
                             returnLabs = FALSE, verbose = TRUE) {
  if(!is.ped(x))
    stop2("Input is not a `ped` object")

  if(!length(ids))
    return(x)

  if(anyDuplicated.default(ids))
    ids = unique.default(ids)

  ids_int = internalID(x, ids, errorIfUnknown = TRUE)
  labs = x$ID

  remove = match.arg(remove)

  if(verbose)
    message(sprintf("Removing individual(s) and their %s: %s", remove, toString(ids)))

  # Descendants OR ancestors.
  remov = switch(remove,
                ancestors = ancestors(x, ids_int, inclusive = TRUE, internal = TRUE),
                descendants = descendants(x, ids_int, inclusive = TRUE, internal = TRUE))

  makeFounderIdx = integer(0)

  if(remove == "descendants") {
    # Remove founders that are NOT parents of any remaining
    FOU = founders(x, internal = TRUE)
    leftovers = .mysetdiff(FOU, c(remov, x$FIDX[-remov], x$MIDX[-remov]))

    if(length(leftovers))
      remov = unique.default(c(remov, leftovers))
  }
  else if(remove == "ancestors") {
    while(TRUE) {
      # Remaining parents (idx) - includes zeroes
      remainFidx = x$FIDX[-remov]
      remainMidx = x$MIDX[-remov]

      # Remaining individuals for which at least 1 parent is removed
      makeFounder = remainFidx %in% remov | remainMidx %in% remov
      makeFounderIdx = seq_along(x$ID)[-remov][makeFounder]

      # Remove spouses who become unattached
      sp = .mysetdiff(c(remainFidx[makeFounder], remainMidx[makeFounder]), remov)
      isParent = sp %in% c(remainFidx[!makeFounder], remainMidx[!makeFounder])
      isChild = !x$FIDX[sp] %in% c(0, remov) & !x$MIDX[sp] %in% c(0, remov)
      leftovers = sp[!isParent & !isChild]

      if(length(leftovers))
        remov = unique.default(c(remov, leftovers))
      else
        break
    }
  }
  if(returnLabs)
    return(labs[.mysortInt(remov)])

  if(verbose) {
    message(sprintf("Total number to be removed: %d. (Remaining: %d)",
                    length(remov), length(labs) - length(remov)))
    if(length(makeFounderIdx))
      message("Converting to founder: ", toString(labs[makeFounderIdx]))
  }

  # The actual reduction. Using the as.matrix trick anticipating marker data a.s.o.
  xmatr = as.matrix(x)
  new = xmatr[-remov, , drop = FALSE]

  if(nrow(new) == 0) {
    if(verbose) message("Remaining pedigree is empty!")
    return(NULL)
  }

  if(length(makeFounderIdx))
    new[match(makeFounderIdx, new[,1]), 2:3] = 0

  attrs = attributes(xmatr)

  # Remove labels
  attrs$LABELS = attrs$LABELS[-remov]

  # Remove founder inbreeding
  finb = attrs$FOUNDER_INBREEDING
  if (!is.null(finb) && length(intersect(remov, FOU))) {
    aut = finb$autosomal
    xchr = finb$x
    attrs$FOUNDER_INBREEDING = list(autosomal = aut[!names(aut) %in% labs[remov]],
                                            x = xchr[!names(xchr) %in% labs[remov]])
  }

  restorePed(new, attrs = attrs)
}

#' @rdname ped_modify
#' @export
branch = function(x, id) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(id) == 0)
    stop2("`id` cannot be empty")
  if(length(id) > 1)
    stop2("`id` must contain a single ID label: ", id)
  desc = descendants(x, id)
  spous = unlist(lapply(c(id, desc), spouses, x = x))
  ids = unique.default(c(id, desc, spous))

  # sort (since subset() does not sort)
  ids = ids[order(internalID(x, ids))]

  subset(x, subset = ids)
}


#' @param subset A character vector (or coercible to such) with ID labels
#'   forming a connected sub-pedigree.
#' @param ... Not used.
#'
#' @rdname ped_modify
#' @export
subset.ped = function(x, subset, ...) {
  if(!is.null(x$LOOP_BREAKERS))
    stop2("`subset()` is not yet implemented for pedigrees with broken loops")

  sub_idx = internalID(x, subset)

  if(anyDuplicated.default(subset))
    stop2("Duplicated ID label: ", unique(subset[duplicated(subset)]))

  pedm = as.matrix(x)
  subped = pedm[sub_idx, , drop = FALSE]

  # set FID = 0 if father is not in subset
  subped[!(subped[, 2] %in% sub_idx), 2] = 0L

  # set MID = 0 if mother is not in subset
  subped[!(subped[, 3] %in% sub_idx), 3] = 0L

  # Fix labels
  attrs = attributes(pedm)
  attrs$LABELS = attrs$LABELS[sub_idx]

  restorePed(subped, attrs, ...)
}

