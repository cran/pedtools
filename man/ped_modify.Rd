% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ped_modify.R
\name{ped_modify}
\alias{ped_modify}
\alias{addChildren}
\alias{addChild}
\alias{addSon}
\alias{addDaughter}
\alias{addParents}
\alias{removeIndividuals}
\alias{trim}
\alias{branch}
\alias{subset.ped}
\title{Add/remove pedigree members}
\usage{
addChildren(
  x,
  father = NULL,
  mother = NULL,
  nch = NULL,
  sex = 1L,
  ids = NULL,
  verbose = TRUE
)

addChild(x, parents, id = NULL, sex = 1, verbose = TRUE)

addSon(x, parents, id = NULL, verbose = TRUE)

addDaughter(x, parents, id = NULL, verbose = TRUE)

addParents(x, id, father = NULL, mother = NULL, verbose = TRUE)

removeIndividuals(
  x,
  ids,
  remove = c("descendants", "ancestors"),
  returnLabs = FALSE,
  verbose = TRUE
)

trim(x, uninformative, verbose = TRUE)

branch(x, id)

\method{subset}{ped}(x, subset, ...)
}
\arguments{
\item{x}{A \code{ped} object, or a list of such.}

\item{father, mother}{Single ID labels. At least one of these must be an
existing member of \code{x}. The other may be (i) another existing member, (ii)
a new founder to be created, or (iii) missing (i.e., NULL), in which case
the other parent is created and given a suitable name.}

\item{nch}{A positive integer indicating the number of children to be
created. Default: 1.}

\item{sex}{Gender codes of the created children (recycled if needed).}

\item{ids}{A vector of ID labels. In \code{addChildren()} these are the children
to be created. If NULL (default) given, automatic labels are generated.}

\item{verbose}{A logical: Verbose output or not.}

\item{parents}{A vector of 1 or 2 ID labels, of which at least one must be an
existing member of \code{x}.}

\item{id}{The ID label of a pedigree member.}

\item{remove}{Either "ancestors" or "descendants" (default), dictating the
method of removing pedigree members. Abbreviations are allowed.}

\item{returnLabs}{A logical, by default FALSE. If TRUE, \code{removeIndividuals()}
returns only the labels of all members to be removed, instead of actually
removing them.}

\item{uninformative}{A vector naming individuals considered "uninformative",
or a function (typically a helper function like \code{\link[=untypedMembers]{untypedMembers()}}).
Uninformative leaves are removed iteratively until no more can be found.}

\item{subset}{A character vector (or coercible to such) with ID labels
forming a connected sub-pedigree.}

\item{...}{Not used.}
}
\value{
The modified \code{ped} object.
}
\description{
Functions for adding or removing individuals in a 'ped' object.
}
\details{
In \code{addChildren()} and \code{addParents()}, labels of added individuals are
generated automatically if they are not specified by the user. The automatic
labelling uses the smallest integers not already in use.

\code{addChild()}, \code{addSon()} and \code{addDaughter()} are convenient wrappers for the
most common use of \code{addChildren()}, namely adding a single child to a
pedigree. Note that the parents can be given in any order. If only one parent
is supplied, the other is created as a new individual.

\code{removeIndividuals()} removes the individuals indicated with \code{ids} along with
all of their ancestors OR descendants, depending on the \code{remove} argument.
Leftover spouses disconnected from the remaining pedigree are also removed.

The \code{branch()} function extracts the sub-pedigree formed by \code{id} and all
his/her spouses and descendants.

The \code{trim()} function iteratively removes uninformative leaves (i.e., members
without children) from the pedigree. Note that the definition of
"uninformative" is entirely user-defined. For example, \code{trim(x, untypedMembers)}, will remove untyped individuals from the bottom until the
process stops.

Finally, \code{subset()} can be used to extract any connected sub-pedigree. (Note
that in the current implementation, the function does not actually check that
the indicated subset forms a connected pedigree; failing to comply with this
may lead to obscure errors.)
}
\examples{

x = nuclearPed(1) |>
  addSon(3) |>
  addParents(4, father = 6, mother = 7) |>
  addChildren(father = 6, mother = 7, nch = 3, sex = c(2,1,2))

# Remove 6 and 7 and their descendants
y1 = removeIndividuals(x, 6:7)

# Remove 8-10 and their parents
y2 = removeIndividuals(x, 8:10, remove = "ancestors")

# Adding a child across components
z = singletons(1:2, sex = 1:2) |> addDaughter(1:2)


}
\seealso{
\code{\link[=ped]{ped()}}, \code{\link[=relabel]{relabel()}}, \code{\link[=swapSex]{swapSex()}}
}
