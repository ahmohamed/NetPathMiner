\name{getGeneSetNetworks}
\alias{getGeneSetNetworks}
\title{Generate geneset networks from an annotated network.}
\usage{
getGeneSetNetworks(graph, use.attr = "pathway")
}
\arguments{
  \item{graph}{An annotated igraph object..}

  \item{use.attr}{The attribute by which vertices are
  grouped (tepically pathway, or GO)}
}
\value{
A list of geneset networks.
}
\description{
This function generates geneset networks based on a given
netowrk, by grouping vertices sharing common attributes (in
the same pathway or compartment).
}
\examples{
data(ex_kgml_sig)	# Ras and chemokine signaling pathways in human
 genesetnets <- getGeneSetNetworks(ex_kgml_sig, use.attr="pathway")
}
\author{
Ahmed Mohamed
}
\seealso{
\code{\link{getGeneSets}}
}

