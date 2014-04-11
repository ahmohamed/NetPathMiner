\name{simplifyReactionNetwork}
\alias{simplifyReactionNetwork}
\title{Removes reactions with no gene annotations}
\usage{
simplifyReactionNetwork(reaction.graph, gene.attr = "genes",
  remove.missing.genes = TRUE, reconnect.threshold = vcount(reaction.graph))
}
\arguments{
  \item{reaction.graph}{A reaction network.}

  \item{gene.attr}{The attribute to be considered as
  "genes". Reactions missing this annotation, will be
  removed.}

  \item{remove.missing.genes}{If \code{FALSE}, only
  tranlocation and spontaneous reactions are removed,
  otherwise all rections with no gene annotations are
  removed.}

  \item{reconnect.threshold}{An argument passed to
  \code{\link{vertexDeleteReconnect}}}
}
\value{
A simplified reaction network.
}
\description{
This function removes reaction vertices with no gene
annotations as indicated by the parameter \code{gene.attr},
and connect their neighbour vertices to preserve graph
connectivity. This is particularly meaningful when
reactions are translocation or spontaneous reactions, which
are not catalysed by genes.
}
\examples{
data(ex_sbml)
 rgraph <- makeReactionNetwork(ex_sbml, simplify=FALSE)

 ## Removes all reaction nodes with no annotated genes.
 rgraph <- simplifyReactionNetwork(rgraph, remove.missing.genes=TRUE)
}
\author{
Ahmed Mohamed
}
\seealso{
Other Network processing methods:
\code{\link{expandComplexes}},
\code{\link{makeGeneNetwork}};
\code{\link{makeReactionNetwork}};
\code{\link{rmSmallCompounds}};
\code{\link{vertexDeleteReconnect}}
}

