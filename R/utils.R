# Preferred version of stop()
stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

# Preferred version of stopifnot()
stopifnot2 = function(...) {
  exprs = list(...)

  for (i in seq_along(exprs)) {
    expri = .subset2(exprs, i)
    if (length(expri) != 1L || is.na(expri) || !expri) {
      full_call = match.call()
      call = deparse(full_call[[i + 1]])
      stop(sQuote(call), " is not TRUE", call. = FALSE, domain = NA)
    }
  }
}

# Test that input is a single integer.
isCount = function(x, minimum = 1, maximum = NA) {
  length(x) == 1 &&
    (is.integer(x) || (is.numeric(x) && x == as.integer(x))) &&
    x >= minimum && (is.na(maximum) || x <= maximum)
}

# Test that input is a single number, with optional range constraints
isNumber = function(x, minimum = NA, maximum = NA) {
  length(x) == 1 &&
    is.numeric(x) &&
    (is.na(minimum) || x >= minimum) &&
    (is.na(maximum) || x <= maximum)
}

`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

seq_ped = function(x)
  seq_len(pedsize(x))

seq_markers = function(x)
  seq_len(nMarkers(x))

# A safer version of base::sample
safe_sample <- function(x, ...) x[sample.int(length(x), ...)]

# Faster (specially for vectors of size 1 and 2) version of sort.int()
.mysortInt = function(v) {
  L = length(v)
  if(L == 1)
    return(v)
  if(L == 2) {
    if(v[1] > v[2])
      return(v[c(2, 1)])
    else return(v)
  }
  sort.int(v, method = "shell")
}

# Fast setdiff
.mysetdiff = function(x, y) unique.default(x[match(x, y, 0L) == 0L])

# Fast intersection. NB: assumes no duplicates!
.myintersect = function(x, y) y[match(x, y, 0L)]

#
.setnames = function (x, nms = x) {
    names(x) = nms
    x
}

# Fast version of t(combn(n,2))
.comb2 = function(n) {
  switch(max(min(n,5),1),
    `dim<-`(integer(0), c(0L, 2L)),
    `dim<-`(c(1L, 2L), c(1L, 2L)),
    `dim<-`(c(1L, 1L, 2L, 2L, 3L, 3L), c(3L, 2L)),
    `dim<-`(c(1L, 1L, 1L, 2L, 2L, 3L, 2L, 3L, 4L, 3L, 4L, 4L), c(6L, 2L)),
    {
      v1 = rep.int(seq_len(n - 1), (n - 1):1)
      v2 = sequence.default((n - 1):1, 2:n)
      cbind(v1, v2, deparse.level = 0)
    })
}

# Random 0/1 vector of length n
.rand01 = function(n) sample.int(2, size = n, replace = TRUE) - 1


stopifnotSimpleVector = function(x, argname = "x") {
  if(is.null(x))
    return()

  if(!is.vector(x)) {
    errmess = sprintf("argument `%s` must be a vector", argname)

    cl = class(x)[1]
    if(!cl %in% c("numeric", "integer", "character", "logical", "double"))
      errmess = sprintf("%s; received an object of class '%s'", errmess, cl)

    stop2(errmess)
  }
}

# Stripped version of expand.grid
fast.grid = function(argslist, as.list = FALSE) {
  nargs = length(argslist)
  orep = nr = prod(lengths(argslist))
  if (nargs == 0L || nr == 0L)
    return(if(as.list) list() else matrix(ncol = 0, nrow = 0))

  rep.fac = 1L
  res = NULL
  for (x in argslist) {
    nx = length(x)
    orep = orep/nx
    res = c(res, x[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep)])  #this is res[, i]
    rep.fac = rep.fac * nx
  }
  dim(res) = c(nr, nargs)
  if (as.list)
    res = lapply(seq_len(nr), function(r) res[r, ])
  res
}

# Add string to certain data.frame entries without disrupting the alignment
# df = data.frame; i = column; pred = logical(nrow(df)); comment = string
commentAndRealign = function(df, i, pred, comment) {
  stopifnot2(is.logical(pred), length(pred) == nrow(df))
  padding = strrep(" ", nchar(comment))

  if(padding == "" || !any(pred))
    return(df)

  df[[i]] = paste0(df[[i]], ifelse(pred, comment, padding))

  if(!is.null(names(df)))
    names(df)[i] = paste0(names(df)[i], padding)

  df
}

# Check that all elements (typically vectors) of a list are identical
listIdentical = function(x) {
  if(length(x) <= 1)
    return(TRUE)
  all(vapply(x[-1], identical, y = x[[1]], logical(1)))
}


# Fold a single string at roughly the given width; try to break at nice places
# Not optimised, and probably reinventing the wheel here.
# Perhaps better idea: Start with strsplit(s, "")
smartfold = function(s, width = 10, breakAt = c(' ', '-', '.', ':', ')', ']')) {
  width = as.integer(max(width, 2))
  nch = nchar(s)
  if(nch < width + 3)
    return(s)

  res = character(0)
  remaining = s
  while(nch >= width + 3) {
    b = width # default next break, if no better
    for (ch in breakAt) {
      pos = gregexpr(ch, remaining, fixed = TRUE)[[1]] |> as.integer() # all positions
      goodpos = pos[pos >= width - 2 & pos <= width + 3 & pos <= nch - 3]
      if(length(goodpos)) {
        b = max(goodpos)
        break
      }
    }
    res = c(res, substr(remaining, 1, b))
    remaining = substr(remaining, b + 1, nch)
    nch = nchar(remaining)
  }
  res = c(res, remaining)
  paste0(res, collapse = "\n")
}

trunc = function(x, printMax = 8, head = 3) {
  if(length(x) <= printMax)
    return(toString(x))
  maxHead = max(1, min(head, printMax - 2))
  y = c(head(x, maxHead), "...", x[length(x)])
  toString(y)
}
