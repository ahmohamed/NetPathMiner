###############################################################################
#
# NPM-package.R: This file contains the documentation for the package itself, and .onLoad function.
# author: Ahmed Mohamed <mohamed@kuicr.kyoto-u.ac.jp>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################

#' General framework for network extraction, path mining.
#'
#' NetPathMiner implements a flexible module-based process flow for network path mining and visualization,
#' which can be fully inte-grated with user-customized functions. 
#' NetPathMiner supports construction of various types of genome scale networks from KGML, SBML and BioPAX
#' formats, enabling its utility to most common pathway databases. 
#' NetPathMiner also provides different visualization techniques to facilitate the analysis of even 
#' thousands of output paths. 
#' 
#' @import igraph
#' @useDynLib NetPathMiner
#' @author Ahmed Mohamed \email{mohamed@@kuicr.kyoto-u.ac.jp}
#' @name NetPathMiner-package
#' @aliases NetPathMiner NPM
#' @docType package
#' 
NULL

#' Biopax example data
#' 
#' A dataset containing Porphyrin metabolism pathway in Biopax Level 3 and parsed with
#' \code{\link[rBiopaxParser]{readBiopax}}.
#' 
#' @docType data
#' @name ex_biopax
#' @examples
#' data(ex_biopax)
#' ex_biopax
#' 
NULL

#' Singaling network from KGML example
#' 
#' An example igraph object representing Ras and chemokine signaling pathways in human 
#' extracted from KGML files.
#' 
#' @docType data
#' @name ex_kgml_sig
#' @examples
#' data(ex_kgml_sig)
#' plotNetwork(ex_kgml_sig, vertex.color="pathway")
#' 
NULL

#' Metabolic network from SBML example
#' 
#' An example igraph object representing bipartite metabolic network of Carbohydrate
#' metabolism extracted from SBML file from Reactome database.
#' 
#' @docType data
#' @name ex_sbml
#' @examples
#' data(ex_sbml)
#' plotNetwork(ex_sbml, vertex.color="compartment.name")
#' 
NULL

#' An microarray data example.
#' 
#' An microarray data example. This is part of the ALL dataset, for demonstration purposes.
#' 
#' @docType data
#' @name ex_microarray
#' @examples
#' data(ex_microarray)
#' 
NULL

.onLoad<- function(lib, pkg){
    env <- new.env()
    load(system.file("extdata", "env_data.RData", package="NetPathMiner"), envir=env)
    env$memory.err <- character()
    options(NPM_ENV = env)
}

#' Internal method to register memery errors.
#' 
#' Internal method to register memery errors, caused by compiled code. This method
#' is used only by the package, and should not be invoked by users. 
#' 
#' @param method The mathod which generated the error. 
#' 
#' @return NULL
#' 
#' @author Ahmed Mohamed
#' @export
#'  
registerMemoryErr <- function(method){
    env <-options("NPM_ENV")[[1]]
    env$memory.err <- c(env$memory.err, method)
}
