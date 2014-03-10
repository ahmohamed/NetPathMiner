\name{pathCluster}
\alias{pathCluster}
\title{3M Markov mixture model for clustering pathways}
\usage{
pathCluster(ybinpaths, M, iter = 1000)
}
\arguments{
  \item{ybinpaths}{The training paths computed by
  \code{\link{pathsToBinary}}.}

  \item{M}{The number of clusters.}

  \item{iter}{The maximum number of EM iterations.}
}
\value{
A list with the following items: \item{h}{The posterior
probabilities that each path belongs to each cluster.}
\item{labels}{The cluster membership labels.}
\item{theta}{The probabilities of each gene for each
cluster.} \item{proportions}{The mixing proportions of each
path.} \item{likelihood}{The likelihood convergence
history.} \item{params}{The specific parameters used.}
}
\description{
3M Markov mixture model for clustering pathways
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
		weight.method = "cor", use.attr="miriam.uniprot", bootstrap = FALSE)

	## Get ranked paths using probabilistic shortest paths.
 ranked.p <- pathRanker(rgraph, method="prob.shortest.path",
					K=20, minPathSize=6)

	## Convert paths to binary matrix.
	ybinpaths <- pathsToBinary(ranked.p)
	p.cluster <- pathCluster(ybinpaths, M=3)
	plotClusters(ybinpaths, p.cluster)
}
\author{
Ichigaku Takigawa

Timothy Hancock
}
\references{
Mamitsuka, H., Okuno, Y., and Yamaguchi, A. 2003. Mining
biologically active patterns in metabolic pathways using
microarray expression profiles. SIGKDD Explor. News l. 5, 2
(Dec. 2003), 113-121.
}
\seealso{
Other Path clustering & classification methods:
\code{\link{pathClassifier}}; \code{\link{pathsToBinary}};
\code{\link{plotClassifierROC}};
\code{\link{plotClusterMatrix}},
\code{\link{plotClusterProbs}}, \code{\link{plotClusters}};
\code{\link{plotPathClassifier}};
\code{\link{plotPathCluster}};
\code{\link{predictPathClassifier}};
\code{\link{predictPathCluster}}
}

