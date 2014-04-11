\name{pathsToBinary}
\alias{pathsToBinary}
\title{Converts the result from pathRanker into something suitable for pathClassifier or pathCluster.}
\usage{
pathsToBinary(ypaths)
}
\arguments{
  \item{ypaths}{The result of \code{\link{pathRanker}}.}
}
\value{
A list with the following elements. \item{paths}{All paths
within ypaths converted to a binary string and concatenated
into the one matrix.} \item{y}{The response variable.}
\item{pidx}{An matrix where each row specifies the location
of that path within the \code{ypaths} object.}
}
\description{
Converts the result from pathRanker into something suitable
for pathClassifier or pathCluster.
}
\details{
Converts a set of pathways from \code{\link{pathRanker}}
into a list of binary pathway matrices. If the pathways are
grouped by a response label then the \emph{pathsToBinary}
returns a list labeled by response class where each element
is the binary pathway matrix for each class. If the
pathways are from \code{\link{pathRanker}} then a list wiht
a single element containing the binary pathway matrix is
returned. To look up the structure of a specific binary
path in the corresponding \code{ypaths} object simply use
matrix index by calling
\code{ypaths[[ybinpaths\$pidx[i,]]]}, where \code{i} is the
row in the binary paths object you wish to reference.
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

	## Convert paths to binary matrix.
	ybinpaths <- pathsToBinary(ranked.p)
	p.cluster <- pathCluster(ybinpaths, M=3)
	plotClusters(ybinpaths, p.cluster, col=c("red", "green", "blue") )
}
\author{
Timothy Hancock and Ichigaku Takigawa
}
\seealso{
Other Path clustering & classification methods:
\code{\link{pathClassifier}}; \code{\link{pathCluster}};
\code{\link{plotClassifierROC}};
\code{\link{plotClusterMatrix}},
\code{\link{plotClusterProbs}}, \code{\link{plotClusters}};
\code{\link{plotPathClassifier}};
\code{\link{plotPathCluster}};
\code{\link{predictPathClassifier}};
\code{\link{predictPathCluster}}
}

