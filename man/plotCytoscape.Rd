\name{plotCytoscape}
\alias{plotCytoscape}
\title{Plots an annotated igraph object in Cytoscape.}
\usage{
plotCytoscape(graph, title, layout = layout.auto, vertex.size, vertex.label,
  vertex.shape, vertex.color, edge.color)
}
\arguments{
  \item{graph}{An annotated igraph object.}

  \item{title}{Will be set as a window title in Cytoscape.}

  \item{layout}{Either a graph layout function, or a
  two-column matrix specifiying vertex coordinates.}

  \item{vertex.size}{Vertex size. If missing, the vertex
  attribute "size" (\preformatted{V(g)$size)}) will be
  used.}

  \item{vertex.label}{Vertex labels. If missing, the vertex
  attribute "label" (\preformatted{V(g)$label)}) will be
  used. If missing, vertices are labeled by their name.}

  \item{vertex.shape}{Vertex shape in one of igraph shapes.
  If missing, the vertex attribute "shape"
  (\preformatted{V(g)$shape)}) will be used. Shapes are
  converted from igraph convention to Cytoscape convention.
  "square","rectangle" and "vrectangle" are converted to
  "RECT", "csquare" and "crectangle" are converted to
  "ROUND_RECT", all other shapes are considered "ELLIPSE"}

  \item{vertex.color}{A color or a list of colors for
  vertices. Vetices with multiple colors are not supported.
  If missing, the vertex attribute "color"
  (\preformatted{V(g)$color)}) will be used.}

  \item{edge.color}{A color or a list of colors for edges.
  If missing, the edge attribute "color"
  (\preformatted{E(g)$color)}) will be used.}
}
\value{
A CytoscapeWindow object constructed by
\code{\link[RCytoscape]{new.CytoscapeWindow}}.
}
\description{
This function uses RCytoscape interface to plot igraph
object in Cytoscape, enabling interactive investigation of
the network. The function requires an open Cytoscape
window, with CytoscapeRPC plugin activated.
}
\examples{
\dontrun{
 data("ex_kgml_sig")
 v.layout <- layoutVertexByAttr(ex_kgml_sig, "pathway")

	cw<-plotCytoscape(ex_kgml_sig, title="example", layout = v.layout,
				vertex.size = 5, vertex.color = v.color)
}
}
\author{
Ahmed Mohamed
}
\seealso{
Other Plotting methods: \code{\link{colorVertexByAttr}};
\code{\link{layoutVertexByAttr}};
\code{\link{plotAllNetworks}};
\code{\link{plotClassifierROC}};
\code{\link{plotClusterMatrix}},
\code{\link{plotClusterProbs}}, \code{\link{plotClusters}};
\code{\link{plotNetwork}};
\code{\link{plotPathClassifier}}; \code{\link{plotPaths}}
}

