\name{vertexDeleteReconnect}
\alias{vertexDeleteReconnect}
\title{Network editing: removing vertices and connecting their neighbours}
\usage{
vertexDeleteReconnect(graph, vids, reconnect.threshold = vcount(graph),
  copy.attr = NULL)
}
\arguments{
  \item{graph}{A reaction network.}

  \item{vids}{Vertex ids to be removed.}

  \item{reconnect.threshold}{If the shortest path between
  vertices is larger than this threshold, they are not
  reconnected.}

  \item{copy.attr}{A function, or a list of functions,
  combine edge attributes. Edge attributes of new edges
  (between reconnected neighbours) are obtained by
  combining original edges attributes along the shortest
  path between reconnected neighbors.}
}
\value{
A modified graph.
}
\description{
This function removes vertices given as \code{vids} and
connects their neighbours as long as the shortest path
beween the neighbours are below the
\code{reconnect.threshold}.
}
\examples{
## Remove all reaction vertices from a bipartite metabolic network
	##  keeping only metabolite vertices.
 data(ex_sbml)
 graph <- vertexDeleteReconnect(ex_sbml, vids=which(V(ex_sbml)$reactions))
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
\code{\link{simplifyReactionNetwork}}
}

