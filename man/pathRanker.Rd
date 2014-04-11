\name{pathRanker}
\alias{pathRanker}
\title{Extracting and ranking paths from a network}
\usage{
pathRanker(graph, method = c("prob.shortest.path", "pvalue"), start, end,
  verbose = TRUE, ...)
}
\arguments{
  \item{graph}{A weighted igraph object. Weights must be in
  \code{edge.weights} or \code{weight} edge attributes.}

  \item{method}{Which path ranking method to use.}

  \item{start}{A list of start vertices, given by their
  vertex id.}

  \item{end}{A list of terminal vertices, given by their
  vertex id.}

  \item{verbose}{Whether to display the progress of the
  function.}

  \item{...}{Method-specific parameters. See Details
  section.}
}
\value{
A list of paths where each path has the following items:
\item{gene}{The ordered sequence of genes visited along the
path.} \item{compounds}{The ordered sequence of compounds
visited along the path.} \item{weights}{The ordered
sequence of the log(ECDF edge weights) along the path.}
\item{distance}{The sum of the log(ECDF edge weights) along
each path.  (a sum of logs is a product)}
}
\description{
Given a weighted igraph object, path ranking finds a set of
node/edge sequences (paths) to maximize the sum of edge
weights. \code{pathRanker(method="prob.shortest.path")}
extracts the K most probable paths within a weighted
network. \code{pathRanker(method="pvalue")} extracts a list
of paths whose sum of edge weights are significantly higher
than random paths of the same length.
}
\details{
The input here is \code{graph}. A weight must be assigned
to each edge. Bootstrapped Pearson correlation edge weights
can be assigned to each edge by
\code{\link{assignEdgeWeights}}. However the specification
of the edge weight is flexible with the condition that
increasing values indicate stronger relationships between
vertices.

\subsection{Probabilistic Shortest Paths}{
\code{pathRanker(method="prob.shortest.path")} finds the K
most probable loopless paths given a weighted network.
Before the paths are ranked the edge weights are converted
into probabilistic edge weights using the Empirical
Cumulative Distribution (ECDF) over all edge weights.  This
is called ECDF edge weight. The ECDF edge weight serves as
a probabilistic rank of the most important gene-gene
interactions.  The probabilistic nature of the ECDF edge
weights allow for a significance test to determine if a
path contains any functional structure or is simply a
random walk.  The probability of a path is simily the
product of all ECDF weights along the path. This is
computed as a sum of the logs of the ECDF edge weights.

The follwing arguments can be passed to
\code{pathRanker(method="prob.shortest.path")}: \describe{
\item{\code{K}}{Maximum number of paths to extract.
Defaults to 10.} \item{\code{minPathSize}}{The minimum
number of edges for each extracted path. Defualts to 1.}
\item{\code{normalize}}{Specify if you want to normalize
the probabilistic edge weights (across different labels)
before extracting the paths. Defaults to TRUE.} } }

\subsection{P-value method}{
\code{pathRanker(method="pvalue")} searches all paths
between the specified start and end vertices, and if a
significant path is found it returns it. However, It
doesn't search for the best path between the start and
terminal vertices, as there could be many paths which lead
to the same terminal vertex, and searching through all of
them is time comsuming.  We just stop when the first
significant path is found.

All provided edge weights are recaled from 0-1. Path
significance is calculated based on the empirical
distribution of random paths of the same length. This can
be estimated using \code{\link{samplePaths}} and passed as
an argument.

The follwing arguments can be passed to
\code{pathRanker(method="pvalue")}: \describe{
\item{\code{sampledpaths}}{The emripical results from
\code{\link{samplePaths}}.} \item{\code{alpha}}{The P value
cut-off. Defualts to 0.01} } }
}
\examples{
## Prepare a weighted reaction network.
	## Conver a metabolic network to a reaction network.
 data(ex_sbml) # bipartite metabolic network of Carbohydrate metabolism.
 rgraph <- makeReactionNetwork(ex_sbml, simplify=TRUE)

	## Assign edge weights based on Affymetrix attributes and microarray dataset.
 # Calculate Pearson's correlation.
	data(ex_microarray)	# Part of ALL dataset.
	rgraph <- assignEdgeWeights(microarray = ex_microarray, graph = rgraph,
		weight.method = "cor", use.attr="miriam.uniprot",
		y=factor(colnames(ex_microarray)), bootstrap = FALSE)

	## Get ranked paths using probabilistic shortest paths.
 ranked.p <- pathRanker(rgraph, method="prob.shortest.path",
					K=20, minPathSize=6)

	## Get significantly correlated paths using "p-valvue" method.
	##   First, establish path score distribution by calling "samplePaths"
 pathsample <- samplePaths(rgraph, max.path.length=10,
                        num.samples=100, num.warmup=10)

	##   Get all significant paths with p<0.1
	significant.p <- pathRanker(rgraph, method = "pvalue",
                sampledpaths = pathsample ,alpha=0.1)
}
\author{
Timothy Hancock, Ichigaku Takigawa, Nicolas Wicker and
Ahmed Mohamed
}
\seealso{
getPathsAsEIDs, extractPathNetwork

Other Path ranking methods:
\code{\link{extractPathNetwork}};
\code{\link{getPathsAsEIDs}}; \code{\link{samplePaths}}
}

