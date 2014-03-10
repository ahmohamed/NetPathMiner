\name{toGraphNEL}
\alias{toGraphNEL}
\title{Converts an annotated igraph object to graphNEL}
\usage{
toGraphNEL(graph, export.attr = "")
}
\arguments{
  \item{graph}{An annotated igraph object..}

  \item{export.attr}{A \code{\link{regex}} experssion
  representing vertex attributes to be exported to the new
  graphNEL object. Supplying an empty string "" (default)
  will export all attributes.}
}
\value{
A graphNEL object.
}
\description{
Converts an annotated igraph object to graphNEL
}
\examples{
data(ex_kgml_sig)	# Ras and chemokine signaling pathways in human
 graphNEL <- toGraphNEL(ex_kgml_sig, export.attr="^miriam.")
}
\author{
Ahmed Mohamed
}

