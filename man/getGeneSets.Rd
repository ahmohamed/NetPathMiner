\name{getGeneSets}
\alias{getGeneSets}
\title{Generate genesets from an annotated network.}
\usage{
getGeneSets(graph, use.attr = "pathway", gene.attr = "genes")
}
\arguments{
  \item{graph}{An annotated igraph object..}

  \item{use.attr}{The attribute by which vertices are
  grouped (tepically pathway, or GO)}

  \item{gene.attr}{The attribute listing genes annotated
  with each vertex (ex: miriam.ncbigene, miriam.uniprot,
  ...)}
}
\value{
A list of genesets.
}
\description{
This function generates genesets based on a given netowrk,
by grouping vertices sharing common attributes (in the same
pathway or compartment). Genes associated with each vertex
can be specified through \code{gene.attr} argument.
}
\examples{
data(ex_kgml_sig)	# Ras and chemokine signaling pathways in human
 genesets <- getGeneSets(ex_kgml_sig, use.attr="pathway", gene.attr="miriam.ncbigene")
}
\author{
Ahmed Mohamed
}
\seealso{
\code{\link{getGeneSetNetworks}}
}

