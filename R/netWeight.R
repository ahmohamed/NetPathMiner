###############################################################################
#
# netWeight.R:     This file contains all functions assigning edge weights to 
# networks, as well as setting and getting annotation attributes and generating
# gene sets. 
# 
# author: Ahmed Mohamed <mohamed@kuicr.kyoto-u.ac.jp>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################

#' Get / Set vertex attribute names and coverage
#' 
#' These functions report the annotation status of the vertices of a given network, modify
#' or remove certain annotations. 
#' 
#' NetPathMiner stores all its vertex annotation attributes in a list, and stores them collectively as
#' a single \code{attr}. This is not to interfer with \code{\link[igraph]{attributes}} from \code{igraph} package.
#' All functions here target NetPathMiner annotations only.  
#' 
#' @param graph An annotated igraph object. 
#' @param pattern A \code{\link{regex}} experssion representing attribute name pattern.
#' 
#' @return For \code{getAttrStatus}, a dataframe summarizing the number of vertices with no (\code{missing}), one (\code{single}) 
#' or more than one (\code{complex}) attribute value. The coverage% is also reported to each attribute.
#' 
#' @author Ahmed Mohamed
#' @family Attribute handling methods 
#' @export
#' @rdname getAttr
#' @examples
#'  data(ex_kgml_sig)	# Ras and chemokine signaling pathways in human
#' 
#'  # Get status of attribute "pathway" only
#'  getAttrStatus(ex_kgml_sig, "^pathway$")
#'  
#'  # Get status of all attributes  starting with "pathway" and "miriam" keywords
#'  getAttrStatus(ex_kgml_sig, "(^miriam)|(^pathway)")
#' 
getAttrStatus <- function(graph, pattern="^miriam."){
    attr.names <- getAttrNames(graph, pattern)
    if(length(attr.names)==0)
        stop("No attributes matching the pattern")
    
    attr.info <- sapply(attr.names,function(x)
                            lapply(getAttribute(graph, x), length))
    
    attr.info <- data.frame(t(apply(attr.info, 2,function(x) c(sum(x==0), sum(x==1), sum(x>1)))))
    names(attr.info) <- c("missing", "single", "complex")
    attr.info <- within(attr.info,coverage.pct<- round((single+complex)/(missing+single+complex)*100) )
    
    return(attr.info)
}

#' @return For \code{getAttrNames}, a character vector of attribute names matching the pattern.
#' 
#' @export
#' @rdname getAttr
#' @examples
#'  # Get all attribute names containing "miriam"
#'  getAttrNames(ex_kgml_sig, "miriam")
getAttrNames <- function(graph, pattern=""){
    if(is.null(V(graph)$attr))
        stop("The igraph object is not annotated.")
    attr.names <- unique(names(unlist(V(graph)$attr, recursive=FALSE)))
    attr.names <- attr.names[grep(pattern, attr.names, ignore.case=TRUE)]
    return(attr.names)
}

#' @param attr.name The attribute name
#' @return For \code{getAttribute}, a list of vertex annotation values for the query attribute.
#' 
#' @export
#' @rdname getAttr
#' @examples
#'  # Get all attribute names containing "miriam"
#'  getAttribute(ex_kgml_sig, "miriam.ncbigene")
#' 
getAttribute <- function(graph, attr.name){
    if(is.null(V(graph)$attr))
        stop("The igraph object is not annotated.")
    return( lapply(V(graph)$attr, "[[", attr.name) )
}

#' @param attr.value A list of attribute values. This must be the same size as the number of vertices. 
#' @return For \code{setAttribute}, a graph with the new attribute set.
#' 
#' @export
#' @rdname getAttr
setAttribute <- function(graph, attr.name, attr.value){
    if(length(attr.value) != vcount(graph))
        stop("Number of provided attribute values doesn't match number of vertices on the graph.")
    
    attr <- V(graph)$attr
    if(is.null(attr))
        attr<-rep(list(list()), vcount(graph))    #initialize the attr as lists.
    
    attr <- mapply(function(attr_, attr_val){
                        attr_[[attr.name]] <- attr_val; return(attr_);
                    }, attr, attr.value, SIMPLIFY=FALSE)
    V(graph)$attr <- attr
    return(graph)
}

#' @return For \code{rmAttrNames}, a new igraph object with the attibute removed.
#' 
#' @export
#' @rdname getAttr
#' @examples
#'  # Remove an attribute from graph
#'  graph <- rmAttribute(ex_kgml_sig, "miriam.ncbigene")
rmAttribute <- function(graph, attr.name){
    if(is.null(V(graph)$attr))
        stop("Graph is not annotated.")
    
    V(graph)$attr <- lapply(V(graph)$attr, function(x) {x[[attr.name]] <- NULL; return(x);})
    return(graph)
}


#' MIRIAM annotation attributes
#' 
#' These functions deals with conforming with MIRIAM annotation guidelines, conversion and mapping between MIRIAM identifiers. 
#'  
#' @param graph An annotated igraph object. 
#' @param return.value Specify whether to return the names of matched standard annotations, or modify the 
#' graph attribute names to match the standards.
#' 
#' @return For \code{stdAttrNames}, \code{matches} gives the original attribute names and their MIRIAM version.
#' Since this is done by simple text matching, mismatches may occur for ambiguous annotations (such as GO, EC number). 
#' \code{graph} returns the input graph with attribute names standardized.
#' 
#' @author Ahmed Mohamed
#' @family Attribute handling methods 
#' @export
#' @rdname MIRIAM
#' @examples
#'  data(ex_kgml_sig)	# Ras and chemokine signaling pathways in human	
#'  ## Modify attribute names to match MIRIAM standard annotations.
#'  graph <- stdAttrNames(ex_kgml_sig, "graph")
#'  
stdAttrNames <- function(graph, return.value=c("matches", "graph")){
    attr.names <- getAttrNames(graph, "miriam")
    suffix <- sub("miriam.(.+)$", "", attr.names)
    db.name <- sub("^(.*)miriam.(.+)$", "miriam.\\2", attr.names)
    db.name <- sub("miriam.obo.", "miriam.",db.name)
    bridge <- NPMdefaults("bridge")
    
    miriam.matches <- sapply(db.name,  function(x)bridge$miriam[grep(x, bridge$miriam, ignore.case=TRUE)])
    sh.name.matches <- sapply(sub("miriam.", "", db.name),  function(x)bridge$miriam[agrep(x, bridge$short.name, ignore.case=TRUE)])
    name.matches <- sapply(sub("miriam.", "", db.name),  function(x)bridge$miriam[grep(x, bridge$name, ignore.case=TRUE)])
    
    matches <- mapply(function(...)head(c(...), n=1L), miriam.matches, name.matches, sh.name.matches, SIMPLIFY=FALSE)
    names(matches) <- attr.names
    matches <- unlist(matches)
    matches <- data.frame(cbind(standard=matches, type=split(bridge$type, bridge$miriam)[matches]))
    matches$standard <- paste(sub("miriam.(.+)$", "", rownames(matches)),matches$standard, sep="")
    
    if("graph" %in% return.value){
        V(graph)$attr <- lapply(1:vcount(graph), function(i)
                    setNames(V(graph)[i]$attr[rownames(matches)], matches$standard)
        )
        if("matches" %in% return.value)
            return(list(matches=matches, graph=graph))
        else return(graph)
    }
    return(matches)
}

# TODO: Check pattern of identifier.
#' @param organism The latin name of the organism (Case-sensitive). 
#' @param target.attr The target annotation, given as MIRIAM standard in the format \code{miriam.xxx}
#' @param source.attr The source annotation attribute from \code{graph}
#' @param bridge.web The base URL for Brigde Database webservices.
#' 
#' @return For \code{fetchAttribute}, the input \code{graph} with the fetched attribute mapped to vertices.  
#' 
#' @export
#' @rdname MIRIAM
#' @examples
#'  # Use Attribute fetcher to get affymetrix probeset IDs for network vertices.
#'  \dontrun{
#'    graph <- fetchAttribute(graph, organism="Homo sapiens", 
#'                    target.attr="miriam.affy.probeset")
#'  }
#'
fetchAttribute <- function(graph, organism="Homo sapiens", target.attr, source.attr, bridge.web=NPMdefaults("bridge.web")){
    if(!require(RCurl))
        stop("This function uses RCurl package. Required package not installed.")
    if(!RCurl::url.exists(bridge.web))
        stop("Couldn't access BridgeDB webservice.\nThere may be a internet connection problem, or the server is down.")    
    if(!organism %in% NPMdefaults("bridge.organisms"))
        stop(organism, " is not supported. Here are supported organisms:\n",
                toString(NPMdefaults("bridge.organisms")))
    if(missing(target.attr))
        stop("Target attribute not specified.")
    
    bridge <- NPMdefaults("bridge") 
    t.code <- bridge$short.name[match(target.attr, bridge$miriam)]
    t.code <- na.omit(t.code)
    
    if(length(t.code)==0)
        stop("Target attribute not supported. Type NPMdefaults(\"bridge\")$miriam for supported attributes.")
    
    std.ann <- stdAttrNames(graph, "matches")
    if(missing(source.attr)){    
        source.attr <- grep("^miriam", rownames(std.ann), value=TRUE)
    }else{
        if(!source.attr %in% rownames(std.ann))
            stop("Couldn't find source attribute in the network.")
    }
    s.code <- bridge$code[match(std.ann[source.attr,"standard"], bridge$miriam)]
    names(s.code) <- source.attr
    s.code <- na.omit(s.code)        

    if(length(s.code)==0)
        stop("Source attribute(s) not supported. Type NPMdefaults(\"bridge\")$miriam for supported attributes.")
    
    base.urls <- paste(NPMdefaults("bridge.web"), URLencode(organism), "xrefs",s.code,sep="/")
    
    all.res <- list()
    for (i in 1:length(s.code)){
        s.attr <- setNames(getAttribute(graph, names(s.code)[i]), V(graph)$name)
        s.attr <- do.call("rbind",lapply(names(s.attr),
                        function(x)
                            if(!is.null(s.attr[[x]]))cbind(x, s.attr[[x]])
                        ))
        
        uris = paste(base.urls[i],"/", unique(s.attr[,2]), "?dataSource=", t.code,sep="")
        query = lapply(RCurl::getURI(uris, async=TRUE, verbose=TRUE), 
                    function(x) if(x!="") as.character(read.table(text=x)$V1) )
        
        s.attr[,2] <- match(s.attr[,2], unique(s.attr[,2]))
        res <- lapply(split(s.attr[,2], s.attr[,1]), 
                    function(x) unlist(query[as.numeric(x)], use.names=FALSE))
        
        all.res[[i]] <- res[V(graph)$name]
    }
    
    final <- do.call(function(...)mapply(c, ...), all.res)
    graph <- setAttribute(graph, target.attr, final)
    return(graph)
}


#' Assigning weights to network edges
#' 
#' This function computes edge weights based on a gene expression profile. 
#' 
#' @param microarray Microarray should be a Dataframe or a matrix, with genes as rownames, and samples as columns.
#' @param graph An annotated igraph object.
#' @param use.attr An attribute name to map \code{microarray} rows (genes) to graph vertices. The attribute must 
#' be annotated in \code{graph}, and the values correspond to \code{rownames} of \code{microarray}. You can check the coverage and 
#' if there are complex vertices using \code{\link{getAttrStatus}}. You can eliminate complexes using \code{\link{expandComplexes}}.
#' @param y Sample labels, given as a factor or a character vector. This must be the same size as the columns of \code{microarray} 
#' @param weight.method A function, or a string indicating the name of the function to be used to compute the edge weights. 
#' The function is provided with 2 numerical verctors (2 rows from \code{microarray}), and it should return a single numerical 
#' value (or \code{NA}). The default computes Pearson's correlation.
#' @param complex.method A function, or a string indicating the name of the function to be used in weighting edges connecting complexes.
#' If a vertex has >1 attribute value, all possible pairwise weights are first computed, and given to \code{complex.method}. The default 
#' function is \code{\link[base]{max}}.
#' @param missing.method A function, or a string indicating the name of the function to be used in weighting edges when one of the vertices
#' lack expression data. The function is passed all edge weights on the graph. Default is \code{\link[stats]{median}}.
#' @param same.gene.penalty A numerical value to be assigned when 2 adjacent vertices have the same attribute value, since correlation and 
#' similarity measure will give perfect scores. Alternatively, \code{same.gene.penalty} can be a function, computing the penalty from all 
#' edge weights on the graph (excluding same-gene and missing values). The default is to take the \code{\link[stats]{median}}  
#' @param bootstrap An integer \code{n}, where the \code{weight.method} is perfomed on \code{n} permutations of the gene profiles, and taking
#' the median value. Set it to \code{NA} to disable bootstrapping. 
#' @param verbose Print the progress of the function.
#' 
#' @return The input graph with \code{edge.weight} as an edge attribute. The attribute can be a list of weights if \code{y} labels 
#' were provided.
#' 
#' @author Ahmed Mohamed
#' @export
#' @examples
#' 	## Convert a metabolic network to a reaction network.
#'  data(ex_sbml) # bipartite metabolic network of Carbohydrate metabolism.
#'  rgraph <- makeReactionNetwork(ex_sbml, simplify=TRUE)
#' 
#' 	## Assign edge weights based on Affymetrix attributes and microarray dataset.
#'  # Calculate Pearson's correlation.
#' 	data(ex_microarray)	# Part of ALL dataset.
#' 	rgraph <- assignEdgeWeights(microarray = ex_microarray, graph = rgraph,
#' 		weight.method = "cor", use.attr="miriam.uniprot", 
#' 		y=factor(colnames(ex_microarray)), bootstrap = FALSE)
#'  
#'  # Using Spearman correlation, assigning missing edges to -1 
#'  \dontrun{
#'    assignEdgeWeights(microarray, graph, use.attr="miriam.affy.probeset", 
#'        y=factor(colnames(microarray)), 
#'        weight.method = function(x1,x2) cor(x1,x2, method="spearman"),
#'        missing.method = -1)
#'  }
#' 
assignEdgeWeights <- function(microarray, graph, use.attr, y, weight.method="cor", 
                                complex.method="max", missing.method="median", same.gene.penalty ="median", 
                                bootstrap = 100, verbose=TRUE) 
{
    # correlation function
    compCor <- function(MA,EL,SAMEG,WEIGHT, BOOTSTRAP) {
        all.cors <- .C("corEdgeWeights",
                as.double(t(MA)),
                as.integer(EL-1),
                as.integer(SAMEG),
                weight = as.double(WEIGHT),
                as.integer(length(EL)/2),
                as.integer(ncol(MA)),
                as.integer(BOOTSTRAP))
        return(all.cors$weight)(all.cors$weight)
    }
    
    if(is.null(V(graph)$attr))
        stop("The igraph object is not annotated.")
    ## Process the weight method
    if(is.null(weight.method)) weight.method <- NA 
    if(!is.function(weight.method)){
        if(is.na(weight.method) || is.numeric(weight.method))
            wt.func <- function(x1,x2) return(weight.method)
        else{
            if(weight.method=="compCor"){
                wt.func <- "compCor"
            }else{
                if(bootstrap==FALSE || bootstrap < 2)
                    wt.func <- get(weight.method)
                else{
                    wt.func <- function(x1,x2){
                        samples <- lapply(1:bootstrap,function(x)sample(length(x1), replace=TRUE))
                        vals <- sapply(samples, function(s) get(weight.method)(x1[s],x2[s])) 
                        return(median(vals))
                    }
                }
            }
        }
    }else{wt.func <- weight.method}
    
    ## Process the complex method
    if(is.null(complex.method)) complex.method <- NA
    if(!is.function(complex.method)){
        if(is.na(complex.method) || is.numeric(complex.method))
            cp.func <- function(...) return(complex.method)
        else{
            cp.func <- get(complex.method)
        }
    }else{cp.func <- complex.method}
    
    ## Process the missing method
    if(is.null(missing.method)) missing.method <- NA
    if(!is.function(missing.method)){
        if(is.na(missing.method) || is.numeric(missing.method))
            ms.func <- function(...) return(missing.method)
        else{
            ms.func <- get(missing.method)
        }
    }else{ms.func <- missing.method}
    
    ## Process the same gene method
    if(is.null(same.gene.penalty)) same.gene.penalty <- NA
    if(!is.function(same.gene.penalty)){
        if(is.na(same.gene.penalty) || is.numeric(same.gene.penalty))
            smgene.func <- function(...) return(same.gene.penalty)
        else{
            smgene.func <- get(same.gene.penalty)
        }
    }else{smgene.func <- same.gene.penalty}
    
    
    
    microarray = as.matrix(microarray) 
    if(!missing(y) && length(y) != ncol(microarray))
        stop("Number of Y labels doesn't match samples in the microarray.")
    
    
    ### Processing the graph.
    genes <- getAttribute(graph, use.attr)
    network.genes <- unique(unlist(genes)) # Get unique list of genes participating in netowrk reactions.    
    intersect.genes <- intersect(network.genes, rownames(microarray))
    if(length(intersect.genes)==0)
        stop("No genes from the microarray are represented in the network.")
    if(verbose){
        cat(sum(!rownames(microarray) %in% intersect.genes), "genes were present in the microarray, but not represented in the network.\n")
        cat(sum(!network.genes %in% intersect.genes), "genes were couldn't be found in microarray.\n")        
    }
    
    
    # Get gene edges 
    edges <- get.edgelist(graph, names=FALSE)
    gene.connections <- do.call("rbind" ,
            lapply(1:dim(edges)[1], function(x) 
                        expand.grid(genes[[edges[x,1]]], genes[[edges[x,2]]], x, stringsAsFactors=FALSE)
            ))
    colnames(gene.connections)    <- c("from", "to", "id")
    
    gene.connections$from <- match( gene.connections$from, rownames(microarray) )
    gene.connections$to <- match( gene.connections$to, rownames(microarray) )
    gene.connections <- gene.connections[complete.cases(gene.connections),]
    if(nrow(gene.connections)==0)
        stop("Couldn't map enough vertices to weight any edge.")
    samegenes <- gene.connections$from == gene.connections$to 
    missing <- which(!(1:nrow(edges) %in% gene.connections$id))        #edges with no gene connections
    
    ### Computing Weights for each y label.
    edge.weights <- NULL
    if (missing(y) || is.null(y)) {
        if(verbose) cat("Assigning edge weights.\n")      
        if(!is.function(wt.func)){
            edge.weights <- compCor(microarray,unlist(gene.connections[,1:2]),
                                samegenes,double(nrow(gene.connections)), bootstrap)
        }else{
            edge.weights <- apply(gene.connections,1, function(x)
                        wt.func(microarray[x[[1]],],microarray[x[[2]],]))
        }
        
        edge.weights <- data.frame(edge.weights, id=gene.connections$id)

        # Add same gene penalties.
        if(sum(samegenes) >0){        
                edge.weights[samegenes,1] <- smgene.func(na.omit(edge.weights[!samegenes,1]))                    
        }
        
        # Add missing values.
        if(length(missing)>0){
            missing.val <- ms.func(na.omit(edge.weights[!samegenes,1]))
            edge.weights <- rbind(edge.weights,cbind(edge.weights=missing.val, id=missing))            
        }
        
        # Complexes
        edge.weights <- suppressWarnings(
                sapply(split(edge.weights[,1], as.numeric(edge.weights$id)), function(x) cp.func( na.omit(x) ) )) 
        
        names(edge.weights) <- "weight"        
    } else {
        y <- as.factor(y)
        y.labels <- c()
        for (yl in 1:nlevels(y)) {  
            if(sum(y == levels(y)[yl])<2){
                warning(levels(y)[yl]," has less than 2 samples. Skipping it.")
                next
            }
            data <- microarray[,y == levels(y)[yl]]
            
            y.labels <- c(y.labels, levels(y)[yl])
            if(verbose) cat("Assigning edge weights for label",levels(y)[yl], "\n") 
            if(!is.function(wt.func)){
                ed <- compCor(data,unlist(gene.connections[,1:2]),
                        samegenes,double(nrow(gene.connections)), bootstrap)
            }else{
                ed <- apply(gene.connections,1, function(x)
                            wt.func(data[x[[1]],],data[x[[2]],]))
            }            
            ed <- data.frame(ed, id=gene.connections$id)
            
            # Add same gene penalties.
            if(sum(samegenes) >0){        
                ed[samegenes,1] <- smgene.func(na.omit(ed[!samegenes,1]))                    
            }
            
            # Add missing values.
            if(length(missing)>0){
                missing.val <- ms.func(na.omit(ed[!samegenes,1]))
                ed <- rbind(ed,cbind(ed=missing.val, id=missing))            
            }

            # Complexes
            ed <- suppressWarnings(
                    sapply(split(ed[,1], as.numeric(ed$id)), function(x) cp.func( na.omit(x) ) ))
            
            if (is.null(edge.weights)) edge.weights <- data.frame(ed)
            else edge.weights <- cbind(edge.weights,data.frame(ed))
        }
        
        names(edge.weights) <- paste("weight:",y.labels,sep = "")
        
        # Edges with same genes are marked with a cor of -1. Set them to the median value.
        #edge.weights <- apply(edge.weights, 2, function(x){ x[x== -1] <- median(x, na.rm=TRUE); return(x) })
        
    }    
    
    # Convert edge.weights to a list of rows.
    E(graph)$edge.weights = as.list(as.data.frame(t(edge.weights)))    
    graph$y.labels = if(missing(y) || is.null(y)) "" else y.labels
    return(graph)
}

#' Generate genesets from an annotated network.
#' 
#' This function generates genesets based on a given netowrk, by grouping vertices sharing
#' common attributes (in the same pathway or compartment). Genes associated with each vertex
#' can be specified through \code{gene.attr} argument.
#' 
#' @param graph An annotated igraph object..
#' @param use.attr The attribute by which vertices are grouped (tepically pathway, or GO)
#' @param gene.attr The attribute listing genes annotated with each vertex (ex: miriam.ncbigene, miriam.uniprot, ...)
#' @param gmt.file Optinal. If provided, Results are exported to a GMT file. GMT files are readily used
#' by most gene set analysis packages. 
#' 
#' @return A list of genesets or written to gmt file if provided.
#' 
#' @author Ahmed Mohamed
#' @seealso \code{\link{getGeneSetNetworks}}
#' @export
#' @examples
#'  data(ex_kgml_sig)	# Ras and chemokine signaling pathways in human
#'  genesets <- getGeneSets(ex_kgml_sig, use.attr="pathway", gene.attr="miriam.ncbigene")
#' 
#' 
#' 	# Write the genesets in a GMT file, and read it using GSEABase package.
#'  getGeneSets(ex_kgml_sig, use.attr="pathway", gene.attr="miriam.ncbigene", gmt.file="kgml.gmt")
#'  \dontrun{
#' 	if(require(GSEABase))
#' 		toGmt("kgml.gmt")
#' 	}
#'
#' 	# Create genesets using compartment information
#' 	data(ex_sbml) # bipartite metabolic network of Carbohydrate metabolism.
#' 	genesets <- getGeneSets(ex_sbml, use.attr="compartment.name", gene.attr="miriam.uniprot")
#'
getGeneSets <- function(graph, use.attr="pathway", gene.attr="genes", gmt.file){    
    attr.names <- getAttrNames(graph)
    if(!use.attr %in% getAttrNames(graph))
        stop(use.attr, ": attribute not found in graph.")
    
    if(!gene.attr %in% getAttrNames(graph))
        stop(gene.attr, ": attribute not found in graph.")
    
    attr <- getAttribute(graph, use.attr)
    attr <- do.call("rbind", lapply(1:length(attr), 
                    function(i) cbind(id=i, val=if(length(attr[[i]])==0) NA else attr[[i]] )))
    
    
    genes <- getAttribute(graph, gene.attr)
    sets <- split(as.numeric(attr[,1]), attr[,2])
    genesets <- lapply(sets, function(x) unlist(genes[x]))
    
    if(!missing(gmt.file)){
        pos <- regexpr("\\.([[:alnum:]]+)$", gmt.file)
        ext <- ifelse(pos > -1L, substring(gmt.file, pos + 1L), "")
        ext <- tolower(ext)
        
        if(!ext == "gmt")
            stop("File format not suuported. Please rename the file as *.gmt")
        
        gmt <- paste(names(genesets), paste(graph$source, use.attr), 
                    lapply(genesets, paste, sep="", collapse="\t"), 
                    sep="\t", collapse="\n")
        
        write(gmt, file=gmt.file)
        
    }else{
        return(genesets)
    }
}

#' Generate geneset networks from an annotated network.
#' 
#' This function generates geneset networks based on a given netowrk, by grouping vertices sharing
#' common attributes (in the same pathway or compartment). 
#' 
#' @param graph An annotated igraph object..
#' @param use.attr The attribute by which vertices are grouped (tepically pathway, or GO)
#' @param format The output format. If "list" is specified, a list of subgraphs are returned (default). 
#' If "pathway-class" is specified, a list of pathway-class objects are returned. Pathway-class
#' is used by graphite package to run several methods of topology-based enrichment analyses.  
#' 
#' @return A list of geneset networks as igraph or Pathway-class objects.
#' 
#' @author Ahmed Mohamed
#' @seealso \code{\link{getGeneSets}}
#' @export
#' @examples
#'  data(ex_kgml_sig)	# Ras and chemokine signaling pathways in human
#'  genesetnets <- getGeneSetNetworks(ex_kgml_sig, use.attr="pathway")
#' 
#'  # Integration with graphite package
#'  \dontrun{
#'  if(require(graphite) & require(clipper) & require(ALL)){
#'		genesetnets <- getGeneSetNetworks(ex_kgml_sig, 
#' 						use.attr="pathway", format="pathway-class")
#'		path <- convertIdentifiers(genesetnets$`Chemokine signaling pathway`, 
#' 						"entrez")
#'		genes <- nodes(path)
#'		data(ALL)
#'		all <- as.matrix(exprs(ALL[1:length(genes),1:20]))
#'		classes <- c(rep(1,10), rep(2,10))
#'		rownames(all) <- genes
#'		
#'		runClipper(path, all, classes, "mean", pathThr=0.1)
#'  }
#'  }
#'
getGeneSetNetworks <- function(graph, use.attr="pathway", format=c("list", "pathway-class")){    
    attr.names <- getAttrNames(graph)
    if(!use.attr %in% getAttrNames(graph))
        stop(use.attr, ": attribute not found in graph.")
    
    attr <- getAttribute(graph, use.attr)
    attr <- do.call("rbind", lapply(1:length(attr), 
                    function(i) cbind(id=i, val=if(length(attr[[i]])==0) NA else attr[[i]] )))
    
    
    sets <- split(as.numeric(attr[,1]), attr[,2])
    genesetnet <- lapply(sets, function(x) induced.subgraph(graph, x))
    
    if(!missing(format) && format=="pathway-class"){
        pathway <- setClass("pathway",
                representation(title="vector",
                        nodes="vector",
                        edges="data.frame",
                        ident="vector",
                        database="vector",
                        timestamp="Date"))
        
        stds <- stdAttrNames(graph, "matches")
        if("miriam.ncbigene" %in% stds$standard){
            annotation <- "miriam.ncbigene" 
        }else if("miriam.uniprot" %in% stds$standard){
            annotation <- "miriam.uniprot"
        }else{
            stop("Provided graph doesn't contain Entrez IDs nor UniProt IDs needed for graphite package")
        }
        
        genesetnet <- lapply(genesetnet, function(g)
                            try(
                                expandComplexes(g, rownames(stds[stds$standard==annotation,]), missing.method="remove"),
                            silent=TRUE)
                )
        
        genesetnet <- genesetnet[ !sapply(genesetnet, class) == "try-error" ]
                
        genesetnet <- genesetnet[ !sapply(genesetnet, ecount) == 0 ]
        ann.prefix <- ifelse(annotation == "miriam.ncbigene", "EntrezGene", "UniProt") 
        genesetnet <- lapply(genesetnet, function(g) 
                                        set.vertex.attribute(g, "name", 
                                                value=paste(ann.prefix, V(g)$name, sep=":") ) 
                            )
                            
        prepareEdges <- function(g){
            ret <- data.frame( get.edgelist(g), 
                    direction="directed", 
                    type=as.character(E(g)$attr)
            )
            names(ret)[1:2] <- c("src", "dest")
            ret$src <- as.character(ret$src)
            ret$dest <- as.character(ret$dest)
            return(ret)
        }
                
                edges <- lapply(genesetnet, function(g){
                                                    
                } )
        
        pathwayset <- mapply(function(name, g){
                    pathway(title=name, 
                            nodes=V(g)$name,
                            edges=prepareEdges(g),
                            ident="native", 
                            database=paste("NetPathMiner(",graph$source,")", sep=""),
                            timestamp=Sys.Date()
                    )                    
                }, names(genesetnet), genesetnet
            )
        
        return(pathwayset)
    }# End pathway-class return
    
    return(genesetnet)
}

#' Converts an annotated igraph object to graphNEL
#' 
#' Converts an annotated igraph object to graphNEL 
#' 
#' @param graph An annotated igraph object..
#' @param export.attr  A \code{\link{regex}} experssion representing vertex attributes to be 
#' exported to the new graphNEL object. Supplying an empty string "" (default) will export
#' all attributes.  
#' 
#' @return A graphNEL object.
#' 
#' @author Ahmed Mohamed
#' @export
#' @examples
#' 	data(ex_kgml_sig)	# Ras and chemokine signaling pathways in human
#'  graphNEL <- toGraphNEL(ex_kgml_sig, export.attr="^miriam.")
#'
toGraphNEL<- function(graph, export.attr=""){
  if(!require(graph))
      stop("This function uses graph package. Required package not installed.")
    attr.names  <- getAttrNames(graph)
    attr.names  <- grep(export.attr, attr.names, value=TRUE)
    
    new.graph <- remove.vertex.attribute(graph, "attr")
    
    for(i in attr.names){
        new.graph <- set.vertex.attribute(new.graph, i, value=getAttribute(graph, i))
    }
    
    return(as_graphnel(new.graph))
}


#' Default values for NetPathMiner
#' 
#' This function gets a NetPathMiner default value for a variable. 
#' 
#' NetPathMiner defines the following defaults:
#' \itemize{
#'   \item small.comp.ls Dataframe of ubiquitous metabolites. Used by \code{\link{rmSmallCompounds}}.    
#'   \item bridge Dataframe of attributes supported by Brigde Database. Used by \code{\link{fetchAttribute}}.
#'   \item bridge.organisms A list of bridge supported organisms. Used by \code{\link{fetchAttribute}}.
#'   \item bridge.web The base URL for Brigde Database webservices. Used by \code{\link{fetchAttribute}}.
#' }
#' @param value a character string indicating the variable name. 
#' 
#' @return The defuult value for the given variable.
#' 
#' @author Ahmed Mohamed
#' @export
#' @examples
#'  # Get the default list of small compounds (uniquitous metabolites).
#'  NPMdefaults("small.comp.ls")
#'  
NPMdefaults <- function(value){
    env = options("NPM_ENV")[[1]]
    tryCatch(return(get(value, env)), error=function(e){
                warning("NetPathMiner doesn't have a default value for ", value,
                        "\n Available defaults are:\n ", toString(ls(env)))
                return(NULL)
            })
} 

