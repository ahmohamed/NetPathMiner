% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netWeight.R
\name{getGeneSetNetworks}
\alias{getGeneSetNetworks}
\title{Generate geneset networks from an annotated network.}
\usage{
getGeneSetNetworks(
  graph,
  use.attr = "pathway",
  format = c("list", "pathway-class")
)
}
\arguments{
\item{graph}{An annotated igraph object..}

\item{use.attr}{The attribute by which vertices are grouped (tepically pathway, or GO)}

\item{format}{The output format. If "list" is specified, a list of subgraphs are returned (default).
If "pathway-class" is specified, a list of pathway-class objects are returned. Pathway-class
is used by graphite package to run several methods of topology-based enrichment analyses.}
}
\value{
A list of geneset networks as igraph or Pathway-class objects.
}
\description{
This function generates geneset networks based on a given netowrk, by grouping vertices sharing
common attributes (in the same pathway or compartment).
}
\examples{
 data(ex_kgml_sig)	# Ras and chemokine signaling pathways in human
 genesetnets <- getGeneSetNetworks(ex_kgml_sig, use.attr="pathway")

 # Integration with graphite package
 \dontrun{
 if(requireNamespace("graphite") & requireNamespace("clipper") & requireNamespace("ALL")){
	genesetnets <- getGeneSetNetworks(ex_kgml_sig,
						use.attr="pathway", format="pathway-class")
	path <- convertIdentifiers(genesetnets$`Chemokine signaling pathway`,
						"entrez")
	genes <- nodes(path)
	data(ALL)
	all <- as.matrix(exprs(ALL[1:length(genes),1:20]))
	classes <- c(rep(1,10), rep(2,10))
	rownames(all) <- genes

	runClipper(path, all, classes, "mean", pathThr=0.1)
 }
 }

}
\seealso{
\code{\link{getGeneSets}}
}
\author{
Ahmed Mohamed
}
