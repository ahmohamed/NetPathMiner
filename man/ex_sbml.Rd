\docType{data}
\name{ex_sbml}
\alias{ex_sbml}
\title{Metabolic network from SBML example}
\description{
An example igraph object representing bipartite metabolic
network of Carbohydrate metabolism extracted from SBML file
from Reactome database.
}
\examples{
data(ex_sbml)
plotNetwork(ex_sbml, vertex.color="compartment.name")
}

