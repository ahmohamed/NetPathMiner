\name{assignEdgeWeights}
\alias{assignEdgeWeights}
\title{Assigning weights to network edges}
\usage{
assignEdgeWeights(microarray, graph, use.attr, y, weight.method = "compCor",
  complex.method = "max", missing.method = "median",
  same.gene.penalty = "median", bootstrap = 100, verbose = TRUE)
}
\arguments{
  \item{microarray}{Microarray should be a Dataframe or a
  matrix, with genes as rownames, and samples as columns.}

  \item{graph}{An annotated igraph object.}

  \item{use.attr}{An attribute name to map
  \code{microarray} rows (genes) to graph vertices. The
  attribute must be annotated in \code{graph}, and the
  values correspond to \code{rownames} of
  \code{microarray}. You can check the coverage and if
  there are complex vertices using
  \code{\link{getAttrStatus}}. You can eliminate complexes
  using \code{\link{expandComplexes}}.}

  \item{y}{Sample labels, given as a factor or a character
  vector. This must be the same size as the columns of
  \code{microarray}}

  \item{weight.method}{A function, or a string indicating
  the name of the function to be used to compute the edge
  weights. The function is provided with 2 numerical
  verctors (2 rows from \code{microarray}), and it should
  return a single numerical value (or \code{NA}). The
  default computes Pearson's correlation.}

  \item{complex.method}{A function, or a string indicating
  the name of the function to be used in weighting edges
  connecting complexes. If a vertex has >1 attribute value,
  all possible pairwise weights are first computed, and
  given to \code{complex.method}. The default function is
  \code{\link[base]{max}}.}

  \item{missing.method}{A function, or a string indicating
  the name of the function to be used in weighting edges
  when one of the vertices lack expression data. The
  function is passed all edge weights on the graph. Default
  is \code{\link[stats]{median}}.}

  \item{same.gene.penalty}{A numerical value to be assigned
  when 2 adjacent vertices have the same attribute value,
  since correlation and similarity measure will give
  perfect scores. Alternatively, \code{same.gene.penalty}
  can be a function, computing the penalty from all edge
  weights on the graph (excluding same-gene and missing
  values). The default is to take the
  \code{\link[stats]{median}}}

  \item{bootstrap}{An integer \code{n}, where the
  \code{weight.method} is perfomed on \code{n} permutations
  of the gene profiles, and taking the median value. Set it
  to \code{NA} to disable bootstrapping.}

  \item{verbose}{Print the progress of the function.}
}
\value{
The input graph with \code{edge.weight} as an edge
attribute. The attribute can be a list of weights if
\code{y} labels were provided.
}
\description{
This function computes edge weights based on a gene
expression profile.
}
\examples{
## Conver a metabolic network to a reaction network.
 data(ex_sbml) # bipartite metabolic network of Carbohydrate metabolism.
 rgraph <- makeReactionNetwork(ex_sbml, simplify=TRUE)

	## Assign edge weights based on Affymetrix attributes and microarray dataset.
 # Calculate Pearson's correlation.
	data(ex_microarray)	# Part of ALL dataset.
	rgraph <- assignEdgeWeights(microarray = ex_microarray, graph = rgraph,
		weight.method = "cor", use.attr="miriam.uniprot",
		y=factor(colnames(ex_microarray)), bootstrap = FALSE)

 # Using Spearman correlation, assigning missing edges to -1
 \dontrun{
        assignEdgeWeights(microarray, graph, use.attr="miriam.affy.probeset",
                y=factor(colnames(microarray)),
                weight.method = function(x1,x2) cor(x1,x2, method="spearman),
                missing.method = -1)
 }
}
\author{
Ahmed Mohamed
}

