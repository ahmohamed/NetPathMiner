###############################################################################
#
# netProcess.R:     This file contains all functions for network processing.
# Conversion between different network representations, network editing are included
# here.
#
# author: Ahmed Mohamed <mohamed@kuicr.kyoto-u.ac.jp>
#
# This is released under GPL-2.
#
# Documentation was created using roxygen
#
###############################################################################

# Neccessary non-sense to pass R CMD check
utils::globalVariables(c("nei", "to", "from", "delete"))

#' Remove uniquitous compounds from a metabolic network
#'
#' This function removes uniquitous compounds (metabolites connected to numerous reactions)
#' from a metabolic network.These compounds are reaction cofactors and currency compounds,
#' such as ATP, CO2, etc. A path through these metabolites may not be bioloigcally meaningful.
#' The defualt small compound list is derived from Reactome, containing keeg.compound, pubchem.compound,
#' ChEBI and CAS identifiers.
#'
#'
#' @param graph A metabolic network.
#' @param method How to handle small compounds. Either simply delete these vertices \code{"remove"} (default),
#' or make a separate vertex for each reaction they participate in \code{"duplicate"}.
#' @param small.comp.ls A list of small compounds to be used.
#'
#' @return A modified graph, with the small compounds removed or duplicated.
#'
#' @author Ahmed Mohamed
#' @family Network processing methods
#' @export
#' @examples
#'     data(ex_sbml)
#'     \dontshow{ex_sbml}
#'
#'     sbml.removed <- rmSmallCompounds(ex_sbml, method="remove")
#'     \dontshow{sbml.removed}
#'
rmSmallCompounds <- function(graph, method=c("remove","duplicate"), small.comp.ls=NPMdefaults("small.comp.ls")){
    if(is.null(V(graph)$reactions))
        stop("graph is not a metaboite-reaction netowrk.")

    if(is.null(V(graph)$attr))
        stop("graph is not annotated.")
    term.list <- c("miriam.pubchem.compound", "miriam.kegg.compound", "miriam.chebi", "miriam.cas")

    std.names <- stdAttrNames(graph, "matches")
    terms <- std.names$standard %in% c("miriam.pubchem.compound", "miriam.kegg.compound", "miriam.chebi")
    if(sum(terms)==0)
        stop("No stuitable compound annotation found.\n
            Current suuported annotations are ", toString(term.list))

    attr.ls <- lapply(rownames(std.names)[terms],
                        function(x) getAttribute(graph, x)[!V(graph)$reactions]
                )


    res <- sapply(1:length(attr.ls), function(i){
                sapply(attr.ls[[i]],
                        function(x) sum(x%in%small.comp.ls[,std.names[terms,1][[i]]])>0
                )
            })

    vids <- which(apply(res,1,any)) # Rows with at least one match.
    vids <- as.integer(V(graph)[!V(graph)$reactions])[vids] #vertex ids on the graph

    if(method=="duplicate"){
        v.names <- V(graph)[vids]$name
        el <- get.edgelist(graph)
        from <- el[,1] %in% v.names
        el[from,1] <- paste(el[from,1], el[from,2], sep="##")
        to <- el[,2] %in% v.names
        el[to,2] <- paste(el[to,2], el[to,1], sep="##")

        graph <- graph + vertices(c(el[from,1],el[to,2])) +igraph::edges(el[from|to,])
    }
    return(graph - vertices(vids))
}

#' Convert metabolic network to reaction network.
#'
#' This function removes metabolite nodes keeping them as edge attributes. The resulting
#' network contains reaction nodes only, where edges indicate that a metabolite produced
#' by one reaction is consumed by the other.
#'
#'
#' @param graph A metabolic network.
#' @param simplify An option to remove translocation and spontaneous reactions that require
#' no catalyzing genes. Translocation reactions are detected from reaction name (SBML, BioPAX), or
#' by having identical substrates and products.
#'
#' @return A reaction network.
#'
#' @author Ahmed Mohamed
#' @family Network processing methods
#' @export
#' @examples
#' 	## Conver a metabolic network to a reaction network.
#'  data(ex_sbml) # bipartite metabolic network of Carbohydrate metabolism.
#'  rgraph <- makeReactionNetwork(ex_sbml, simplify=TRUE)
#'
makeReactionNetwork <- function(graph, simplify=FALSE){
    reactions <- V(graph)$reactions
    if(is.null(reactions) || class(reactions)!="logical" || sum(reactions)==0)
        stop("The graph contains no reactions.")

    edges <- get.edgelist(graph)

    input <- edges[edges[,2]%in% V(graph)[reactions]$name,]  #Reactant-reaction edges
    output <- edges[edges[,1]%in% V(graph)[reactions]$name,]   #Reaction-product edges

    match1 <- match(output[,2],input[,1])            #Find reactions that produces reactants in input
    match2 <- match(input[,1],output[,2])            #Find reactions that consumes products in output

    #Connect matches
    r.connections <- rbind(
            cbind(output[,],input[match1,2]),
            cbind(output[match2,], input[,2])
    )

    # remove missing and duplicate values
    r.connections <- r.connections[complete.cases(r.connections),]
    r.connections <- unique(r.connections)
    r.connections <- data.frame(r.connections[,c(1,3,2)],
            I(V(graph)[r.connections[,2]]$attr))
    names(r.connections) <- c("from", "to", "compound", "attr")

    #### Making a reaction network graph from edges#####
    reaction.graph = simplify(graph.data.frame(r.connections),edge.attr.comb="c")

    V(reaction.graph)$attr <- V(graph)[V(reaction.graph)$name]$attr

    if(simplify)
        reaction.graph <- simplifyReactionNetwork(reaction.graph, remove.missing.genes=FALSE)

    reaction.graph$source <- graph$source
    reaction.graph$type <- "R.graph"
    return(reaction.graph)
}

#' Removes reactions with no gene annotations
#'
#' This function removes reaction vertices with no gene annotations as indicated by the parameter
#' \code{gene.attr}, and connect their neighbour vertices to preserve graph connectivity. This is
#' particularly meaningful when reactions are translocation or spontaneous reactions,
#' which are not catalysed by genes.
#'
#' @param reaction.graph A reaction network.
#' @param gene.attr The attribute to be considered as "genes". Reactions missing this annotation, will be removed.
#' @param remove.missing.genes If \code{FALSE}, only tranlocation and spontaneous reactions are removed, otherwise
#' all rections with no gene annotations are removed.
#' @param reconnect.threshold An argument passed to \code{\link{vertexDeleteReconnect}}
#'
#' @return A simplified reaction network.
#'
#' @author Ahmed Mohamed
#' @family Network processing methods
#' @export
#' @examples
#'  data(ex_sbml)
#'  rgraph <- makeReactionNetwork(ex_sbml, simplify=FALSE)
#'
#'  ## Removes all reaction nodes with no annotated genes.
#'  rgraph <- simplifyReactionNetwork(rgraph, remove.missing.genes=TRUE)
#'
simplifyReactionNetwork <- function(reaction.graph, gene.attr="genes",
                            remove.missing.genes=TRUE,
                            reconnect.threshold=vcount(reaction.graph))
{
    #List of all missing-gene reactions.
    no.gene <- lapply( lapply(V(reaction.graph)$attr, "[[", gene.attr), length) == 0
    if(sum(no.gene)==0)
        return(reaction.graph)

    attr.ls <- if(sum(no.gene)==1) list(V(reaction.graph)[no.gene]$attr)
                else V(reaction.graph)[no.gene]$attr

    translocation = sapply( attr.ls,
            function(x) ( grepl("transloc", x[["name"]], ignore.case=TRUE) ||
                identical(x[["reactants"]],x[["products"]]) )
    )

    spontaneous = sapply( attr.ls,
            function(x) grepl("spontaneous", x[["name"]], ignore.case=TRUE)
    )

    new.reaction.graph = reaction.graph

    if(sum(translocation | spontaneous) > 0){
      new.reaction.graph = vertexDeleteReconnect(reaction.graph,
              V(reaction.graph)[no.gene][translocation | spontaneous]
              ,copy.attr=list(compound=function(...)paste(..., sep="->"), attr="c") )
    }

    if(remove.missing.genes & length(which(no.gene)[!translocation & !spontaneous] >0)){
        new.reaction.graph = vertexDeleteReconnect(new.reaction.graph,
                V(new.reaction.graph)[no.gene], reconnect.threshold)
    }
    return(new.reaction.graph)
}

#' Network editing: removing vertices and connecting their neighbours
#'
#' This function removes vertices given as \code{vids} and connects their neighbours as
#' long as the shortest path beween the neighbours are below the \code{reconnect.threshold}.
#'
#' @param graph A reaction network.
#' @param vids Vertex ids to be removed.
#' @param reconnect.threshold If the shortest path between vertices is larger than this threshold,
#' they are not reconnected.
#' @param copy.attr A function, or a list of functions, combine edge attributes. Edge attributes
#' of new edges (between reconnected neighbours) are obtained by combining original edges attributes
#' along the shortest path between reconnected neighbors.
#'
#' @return A modified graph.
#'
#' @author Ahmed Mohamed
#' @family Network processing methods
#' @export
#' @examples
#'  ## Remove all reaction vertices from a bipartite metabolic network
#' 	##  keeping only metabolite vertices.
#'  data(ex_sbml)
#'  graph <- vertexDeleteReconnect(ex_sbml, vids=which(V(ex_sbml)$reactions))
#'
vertexDeleteReconnect <- function(graph, vids, reconnect.threshold=vcount(graph), copy.attr=NULL){
    if(length(vids)==0){return(graph)}
    V(graph)$delete <- FALSE
    V(graph)[vids]$delete <- TRUE

    #A subgraph only including vertices to be deleted and their neighbours.
    graph.sub <- subgraph.edges(graph,
            E(graph)[V(graph)[delete] %--% V(graph)[nei(V(graph)[delete])|delete]])

    #Calculate the shortest path length in this network between "neighbours" (retained) vertices.
    #In this graph, the shortest path will have to go through deleted vertices only.
    #Path lengths indicate the number of deleted vertices lying between 2 retained vertices.
    short.paths <- shortest.paths(graph.sub, V(graph.sub)[!delete],V(graph.sub)[!delete], "out")
    new.edges <- which(short.paths!=0 & short.paths!=Inf & short.paths< reconnect.threshold+1, arr.ind=TRUE)
    new.edges <- rbind(rownames(new.edges),colnames(short.paths)[new.edges[,2]])  #edges to be added

    #Get Actual shortest paths for vertices that passed the threshold, to copy edge attributes.
    attr <- NULL
    if(!is.null(copy.attr)){
        paths <- mapply(function(from, to) igraph::get.shortest.paths(graph.sub, from, to, output="both",mode="out")$vpath,
                        new.edges[1,], new.edges[2,])

        if(!is.list(copy.attr)) copy.attr <- sapply(list.edge.attributes(graph.sub), function(x)copy.attr)

        attr <- sapply(names(copy.attr), function(y) lapply(paths,
                            function(x) do.call( copy.attr[[y]], as.list(get.edge.attribute(graph.sub, y,E(graph.sub, path=x))) )
                            )
                    ,simplify=FALSE)
    }

    new.graph <- add.edges(graph, new.edges, attr=attr)
    new.graph <- delete.vertices(new.graph, V(graph)[delete]$name)
    new.graph <- remove.vertex.attribute(new.graph, "delete")

    return(new.graph)
}

#' Expand reactions / complexes into their gene constituents.
#'
#' These are general functions to expand vertices by their attributes, i.e. create a separate
#' vertex for each attribute value.
#'
#' These functions can be very useful when merging networks constructed from different databases.
#' For example, to match a network created from Reactome to a KEGG network, you can expand metabolite
#' vertices by "miriam.kegg.compound" attribute.
#'
#' @param graph An annotated igraph object.
#' @param v.attr Name of the attribute which vertices are expanded to.
#' @param keep.parent.attr A (List of) \code{\link{regex}} experssions representing attributes to be
#' inherited by daughter vertices. If \code{"all"} is passed, all parent attributes are inherited.
#' @param expansion.method If \code{"duplicate"}, attribute values sharing more than one parent vertex
#' are duplicated for each vertex they participate in. For exmaple, if one gene G1 catalyzes reactions
#' R1, R2; then G1##R1, and G1##R2 vertices are created. If \code{"normal"} only one vertex (G1) is created,
#' and inherit all R1 and R2 connections and attributes.
#' @param missing.method How to deal with vertices with no attribute values. \code{"keep"} retains the parent
#' node, \code{"remove"} simply deletes the vertex, and \code{"reconnect"} removes the vertex and connect its
#' neighbours to each other (to prevent graph cuts).
#'
#' @return A new graph with vertices expanded.
#'
#' @author Ahmed Mohamed
#' @family Network processing methods
#' @rdname makeGeneNetwork
#' @export
#' @examples
#'  ## Make a gene network from a reaction network.
#'  data(ex_sbml)	# A bipartite metbaolic network.
#'  rgraph <- makeReactionNetwork(ex_sbml, simplify=TRUE)
#'  ggraph <- makeGeneNetwork(rgraph)
#'
#'  ## Expand vertices into their contituent genes.
#'  data(ex_kgml_sig)	# Ras and chemokine signaling pathways in human
#' 	ggraph <- expandComplexes(ex_kgml_sig, v.attr = "miriam.ncbigene",
#' 						keep.parent.attr= c("^pathway", "^compartment"))
#'
#'  ## Create a separate vertex for each compartment. This is useful in duplicating
#' 	##  metabolite vertices in a network.
#' \dontrun{
#'  graph <- expandComplexes(graph, v.attr = "compartment",
#'         keep.parent.attr = "all",
#'         expansion.method = "duplicate",
#'         missing.method = "keep")
#' }
#'
expandComplexes <- function(graph, v.attr,
        keep.parent.attr= "^pathway",
        expansion.method=c("normal", "duplicate"),
        missing.method=c("keep", "remove", "reconnect"))
{
    if(missing(v.attr))
        stop("v.attr: Vertex attribute to be used in expansion is not specified.")

    attr.names <- getAttrNames(graph)
    if(!v.attr %in% attr.names)
        stop(v.attr,": Attribute not found in graph.")

    attr.ls = lapply(V(graph)$attr, "[[", v.attr)

    if(missing(expansion.method)) expansion.method <- "normal"
    else if(!expansion.method %in% c("normal", "duplicate"))
        stop("Unknown expansion method provided: ", expansion.method)

    if(missing(missing.method)) missing.method <- "remove"
    else if(!missing.method %in% c("keep", "remove", "reconnect"))
        stop("Unknown missin method provided: ", missing.method)


    if(missing.method!="remove"){
        no.attr <- lapply(attr.ls, length) ==0
        attr.ls[no.attr] <- V(graph)[no.attr]$name
    }

    attr_func <- function(...){
        l = mapply(function(...) unique(c(...)),...,SIMPLIFY=FALSE)
        l[!is.na(names(l))]
    }
    z = .Call("expand_complexes", ATTR_LS=attr.ls,
            EL=as.integer(t(get.edgelist(graph, names=FALSE))-1),
            V = V(graph)$name,
            EXPAND=expansion.method,
            MISSING=missing.method)

    gout = graph.empty() + vertices(z$vertices)
    gout = gout +igraph::edges(z$edges)

    if(missing.method=="reconnect"){
        no.attr.vids <- which(z$parents %in% which(no.attr))
        if(length(no.attr.vids)>0)
            gout <- vertexDeleteReconnect(gout, no.attr.vids)
    }

    if(!is.null(keep.parent.attr)){
        if(length(keep.parent.attr)==1 && keep.parent.attr=="all"){
            attr_terms = lapply(V(graph)$attr, "[", keep.parent.attr)
            V(gout)$attr <- lapply(z$parents, function(x) do.call("attr_func", attr_terms[x]))
        }else{
            if(length(keep.parent.attr) > 1)
                keep.parent.attr <- do.call("paste",as.list(c(keep.parent.attr, sep="|")))

            keep.parent.attr <- attr.names[grep(keep.parent.attr, attr.names)]

            attr_terms = lapply(V(graph)$attr, "[", keep.parent.attr)
            V(gout)$attr <- lapply(z$parents, function(x) do.call("attr_func", attr_terms[x]))
        }
    }

    gout <- setAttribute(gout, v.attr, V(gout)$name)
    for(i in list.edge.attributes(graph)){
        gout <- set.edge.attribute(gout, i, value=get.edge.attribute(graph, i)[z$e.parents] )
    }

    gout$source <- graph$source
    return(gout);
}

#' @return \code{makeGeneNetwork} returns a graph, where nodes are genes, and edges represent
#' participation in succesive reactions.
#'
#' @rdname makeGeneNetwork
#' @export
makeGeneNetwork <- function(graph, v.attr="genes",
        keep.parent.attr= "^pathway",
        expansion.method="duplicate",
        missing.method="remove"){

    gene.graph <- expandComplexes(graph, v.attr, keep.parent.attr,
            expansion.method, missing.method)

    V(gene.graph)$color <- "blue"
    gene.graph$source <- graph$source
    gene.graph$type = "G.graph"

    return(gene.graph)
}
