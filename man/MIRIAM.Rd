\name{stdAttrNames}
\alias{fetchAttribute}
\alias{stdAttrNames}
\title{MIRIAM annotation attributes}
\usage{
stdAttrNames(graph, return.value = c("matches", "graph"))

fetchAttribute(graph, organism = "Homo sapiens", target.attr, source.attr,
  bridge.web = NPMdefaults("bridge.web"))
}
\arguments{
  \item{graph}{An annotated igraph object.}

  \item{return.value}{Specify whether to return the names
  of matched standard annotations, or modify the graph
  attribute names to match the standards.}

  \item{organism}{The latin name of the organism
  (Case-sensitive).}

  \item{target.attr}{The target annotation, given as MIRIAM
  standard in the format \code{miriam.xxx}}

  \item{source.attr}{The source annotation attribute from
  \code{graph}}

  \item{bridge.web}{The base URL for Brigde Database
  webservices.}
}
\value{
For \code{stdAttrNames}, \code{matches} gives the original
attribute names and their MIRIAM version. Since this is
done by simple text matching, mismatches may occur for
ambiguous annotations (such as GO, EC number). \code{graph}
returns the input graph with attribute names standardized.

For \code{fetchAttribute}, the input \code{graph} with the
fetched attribute mapped to vertices.
}
\description{
These functions deals with conforming with MIRIAM
annotation guidelines, conversion and mapping between
MIRIAM identifiers.
}
\examples{
data(ex_kgml_sig)	# Ras and chemokine signaling pathways in human
 ## Modify attribute names to match MIRIAM standard annotations.
 graph <- stdAttrNames(ex_kgml_sig, "graph")
# Use Attribute fetcher to get affymetrix probeset IDs for network vertices.
 \dontrun{
   graph <- fetchAttribute(graph, organism="Homo sapiens",
                   target.attr="miriam.affy.probeset")
 }
}
\author{
Ahmed Mohamed
}
\seealso{
Other Attribute handling methods:
\code{\link{getAttrNames}}, \code{\link{getAttrStatus}},
\code{\link{getAttribute}}, \code{\link{rmAttribute}},
\code{\link{setAttribute}}
}

