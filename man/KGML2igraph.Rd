\name{KGML2igraph}
\alias{KGML2igraph}
\title{Processes KGML files into igraph objects}
\usage{
KGML2igraph(filename, parse.as = c("metabolic", "signaling"),
  expand.complexes = FALSE, verbose = TRUE)
}
\arguments{
  \item{filename}{A character vector containing the KGML
  files to be processed. If a directory path is provided,
  all *.xml files in it and its subdirectories are
  included.}

  \item{parse.as}{Whether to process file into a metabolic
  or a signaling network.}

  \item{expand.complexes}{Split protein complexes into
  individual gene nodes. This argument is ignored if
  \code{parse.as="metabolic"}}

  \item{verbose}{Whether to display the progress of the
  function.}
}
\value{
An igraph object, representing a metbolic or a signaling
network.
}
\description{
This function takes KGML files as input, and returns either
a metabolic or a signaling network as output.
}
\details{
Users can specify whether files are processes as metabolic
or signaling networks.

Metabolic networks are given as bipartite graphs, where
metabolites and reactions represent vertex types. This is
constructed from <reaction> xml node in KGML file,
connecting them to their corresponding substrates and
products. Each reaction vertex has \code{genes} attribute,
listing all genes associated with the reaction. As a
general rule, reactions inherit all annotation attributes
of its catalyzig genes.

Signaling network have genes as vertices and edges
represent interactions, such as activiation / inhibition.
Genes participating in successive reactions are also
connected. Signaling parsing method processes <ECrel>,
<PPrel> and <PCrel> interactions from KGML files.

To generate a genome scale network, simply provide a list
of files to be parsed, or put all file in a directory, as
pass the directory path as \code{filename}
}
\examples{
if(is.loaded("readkgmlfile")){ # This is false if libxml2 wasn't available at installation.
    filename <- system.file("extdata", "hsa00860.xml", package="NetPathMiner")

    # Process KGML file as a metabolic network
    g <- KGML2igraph(filename)
    plotNetwork(g)

    # Process KGML file as a signaling network
    g <- KGML2igraph(filename, parse.as="signaling", expand.complexes=TRUE)
    plotNetwork(g)
}
}
\author{
Ahmed Mohamed
}
\seealso{
Other Database extraction methods:
\code{\link{SBML2igraph}}; \code{\link{biopax2igraph}}
}

