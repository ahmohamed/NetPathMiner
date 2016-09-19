###############################################################################
#
# plotPath.R:     This file contains all functions related to plotting of 
# networks and paths.
# author: Ahmed Mohamed <mohamed@kuicr.kyoto-u.ac.jp>
#
# This is released under GPL-2.
# 
# Documentation was created using roxygen
#
###############################################################################

#' Plots an annotated igraph object.
#' 
#' This function is a wrapper function for \code{\link[igraph]{plot.igraph}}, with 2 main additions. 
#' 1. Add the ability to color vertices by their attributes (see examples), accompanied by an inofrmative
#' legend. 2. Resize vertex.size, edge.arrow.size, label.cex according to the plot size and the size of the
#' network.
#' 
#' @param graph An annotated igraph object.
#' @param vertex.color A list of colors for vertices, or an attribute names (ex: "pathway") by which vertices
#' will be colored. Complex attributes, where a vertex belongs to more than one group, are supported. This can 
#' also be the output of \code{\link{colorVertexByAttr}}. 
#' @param col.palette A color palette, or a palette generating function (ex: \preformatted{col.palette=rainbow}).
#' @param layout Either a graph layout function, or a two-column matrix specifiying vertex coordinates.
#' @param legend Wheter to plot a legend. The legend is only plotted if vertices are colored by attribute values.
#' @param ... Additional arguments passed to \code{\link[igraph]{plot.igraph}}.
#' 
#' @return 
#' Produces a plot of the network.
#' 
#' @author Ahmed Mohamed
#' @family Plotting methods
#' @export
#' @examples 
#'  data("ex_kgml_sig")
#'  plotNetwork(ex_kgml_sig, vertex.color="pathway")
#'  plotNetwork(ex_kgml_sig, vertex.color="pathway", col.palette=heat.colors)
#'  plotNetwork(ex_kgml_sig, vertex.color="pathway", 
#'              col.palette=c("red", "green","blue","grey"))
#' 
plotNetwork <- function(graph, vertex.color, col.palette = palette(), layout = layout.auto, legend=TRUE,...){
    opar <- par()
    opar[c("cin", "cra", "csi", "cxy", "din", "page")] <- NULL
    on.exit(par(opar))
    
    par(mar=c(0,0,2,0))
    
    # Leave a space in the outer margin for legend.
    if(legend)
        par(omi=c(0,0,0, par("pin")[1]*0.2), new=!par("page"))
    
    plotNetwork_internal(graph, vertex.color, col.palette, layout, legend,...)
}

#' Plots an annotated igraph object higlighting ranked paths.
#' 
#' This function plots a network highlighting ranked paths. If \code{path.clusters} are provided,
#' paths in the same cluster are assigned similar colors.
#' 
#' @param paths The result of \code{\link{pathRanker}}.
#' @param graph An annotated igraph object.
#' @param path.clusters The result from \code{\link{pathCluster}} or \code{\link{pathClassifier}}.
#' @param col.palette A color palette, or a palette generating function (ex: \preformatted{col.palette=rainbow}).
#' @param layout Either a graph layout function, or a two-column matrix specifiying vertex coordinates.
#' @param ... Additional arguments passed to \code{\link{plotNetwork}}.
#' 
#' @return Produces a plot of the network with paths highlighted. If paths are computed for several
#' labels (sample categories), a plot is created for each label.
#' 
#' @author Ahmed Mohamed
#' @family Plotting methods
#' @export
#' @examples
#' 	## Prepare a weighted reaction network.
#' 	## Conver a metabolic network to a reaction network.
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
#' 	## Get ranked paths using probabilistic shortest paths.
#'  ranked.p <- pathRanker(rgraph, method="prob.shortest.path", 
#' 					K=20, minPathSize=6)
#' 
#' 	## Plot paths.
#' 	plotPaths(ranked.p, rgraph)
#' 
#' 	## Convert paths to binary matrix, build a classifier. 
#' 	ybinpaths <- pathsToBinary(ranked.p)
#' 	p.class <- pathClassifier(ybinpaths, target.class = "BCR/ABL", M = 3)
#' 	
#'  ## Plotting with clusters, on a metabolic graph.
#' 	plotPaths(ranked.p, ex_sbml, path.clusters=p.class)
#' 
plotPaths <- function(paths, graph, path.clusters=NULL, col.palette=palette(), layout=layout.auto, ...){
    opar <- par()
    opar[c("cin", "cra", "csi", "cxy", "din", "page")] <- NULL
    on.exit({par(opar);graphics::layout(1)})
    
    process.layout(paths, 1, path.clusters, plot.clusters=FALSE)    
    col.process <- process.color(col.palette, paths, path.clusters, same.plot=FALSE )
    path.col <- col.process[[1]]
    legend.paths <- col.process[[2]]
    
    eids <- getPathsAsEIDs(paths, graph)
    
    
    if(is.null(names(path.col))){
        if( !is.null(names(eids)) ) eids <- unlist(eids, recursive=FALSE)
        legend.v <- highlightPaths(eids, graph, layout=layout, path.col = path.col, ...)
    }else{
        legend.v <- lapply(names(eids), function(x)
                    highlightPaths(eids[[x]], graph, layout=layout, path.col= path.col[[x]],
                            main=paste("y =",x), ...))
        
        legend.v <- legend.v[[1]]
    }
    
    if(length(legend.paths)>0){
        if(length(legend.v)>0){
            drawLegend(vertices=legend.v, paths=legend.paths)
        }else{
            drawLegend(paths=legend.paths)
        }
    }else if (length(legend.v)>0){
        drawLegend(vertices=legend.v)
    }
}

#' Higlighting ranked paths over multiple network representations.
#' 
#' This function highlighting ranked paths over different network representations, metabolic, reaction and
#' gene networks. The functions finds equivalent paths across different networks and marks them.
#' 
#' @param paths The result of \code{\link{pathRanker}}.
#' @param metabolic.net A bipartite metabolic network.
#' @param reaction.net A reaction network, resulting from \code{\link{makeReactionNetwork}}.
#' @param gene.net A gene network, resulting from \code{\link{makeGeneNetwork}}.
#' @param path.clusters The result from \code{\link{pathCluster}} or \code{\link{pathClassifier}}.
#' @param plot.clusters Whether to plot clustering information, as generated by \code{\link{plotClusters}}
#' @param col.palette A color palette, or a palette generating function (ex: \preformatted{col.palette=rainbow}).
#' @param layout Either a graph layout function, or a two-column matrix specifiying vertex coordinates.
#' @param ... Additional arguments passed to \code{\link{plotNetwork}}.
#' 
#' @return Highlights the path list over all provided networks.
#' 
#' @author Ahmed Mohamed
#' @family Plotting methods
#' @export
#' @examples
#' 	## Prepare a weighted reaction network.
#' 	## Conver a metabolic network to a reaction network.
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
#' 	## Get ranked paths using probabilistic shortest paths.
#'  ranked.p <- pathRanker(rgraph, method="prob.shortest.path", 
#' 					K=20, minPathSize=6)
#' 
#' plotAllNetworks(ranked.p, metabolic.net = ex_sbml, reaction.net = rgraph,
#' 					vertex.label = "", vertex.size = 4)
#' 
plotAllNetworks <- function(paths, metabolic.net=NULL, reaction.net=NULL, gene.net=NULL, 
        path.clusters=NULL, plot.clusters=TRUE, col.palette=palette(), layout=layout.auto,...){
    opar <- par()
    opar[c("cin", "cra", "csi", "cxy", "din", "page")] <- NULL
    on.exit({par(opar);graphics::layout(1)})
    
    # Counting the number of networks to be plotted
    numPlots = length(Filter(Negate(is.null), list(metabolic.net, reaction.net, gene.net)))
    if(numPlots < 1) stop("No networks provided for plotting")
    
    process.layout(paths,numPlots, path.clusters, plot.clusters)    
    col.process <- process.color(col.palette, paths, path.clusters, !is.null(path.clusters$y) && plot.clusters )
    path.col <- col.process[[1]]
    legend.paths <- col.process[[2]]
    
    params <- list(paths=paths, metabolic.net=metabolic.net, reaction.net=reaction.net, 
            gene.net=gene.net, clusters=path.clusters, layout=layout, path.col=path.col, ...)
    
    # Plotting networks    
    do.call("plotAllNetworks_internal", params)
    
    par(mar = c(5,2,0,2)+0.1, xpd=FALSE)
    # Plotting cluster information
    if(!is.null(path.clusters) && plot.clusters){        
        ybinpaths <- pathsToBinary(paths)
        
        if(!is.function(col.palette))
            col.palette <- colorRampPalette(col.palette)
        
        plotClusterMatrix(ybinpaths, path.clusters, col= col.palette(path.clusters$params$M))
        plotClusterProbs(path.clusters, col= col.palette(path.clusters$params$M))        
    }
    
    if(length(legend.paths)>0)
        drawLegend(paths=legend.paths)
}

#' A graph layout function, which groups vertices by attribute.
#' 
#' This function generates a layout for igraph objects, keeping vertices with the same attribute
#' (ex: in the same pathway, etc) close to each other. 
#' 
#' @param graph An annotated igraph object.
#' @param attr.name The attribute name by which vertices are laid out.
#' @param cluster.strength A number indicating tie strengths between vertices with the same attribute. 
#' The larger it is, the closer the vertices will be.
#' @param layout A layout function, ideally a force-directed layout fuction, such as 
#' \code{\link[igraph]{layout.fruchterman.reingold}} and \code{\link[igraph]{layout.kamada.kawai}}. 
#' 
#' @return A two-column matrix indicating the x and y postions of vertices.
#' 
#' @author Ahmed Mohamed
#' @family Plotting methods
#' @export
#' @examples 
#'   data("ex_kgml_sig")
#'   v.layout <- layoutVertexByAttr(ex_kgml_sig, "pathway") 
#'   plotNetwork(ex_kgml_sig, vertex.color="pathway", layout=v.layout)
#' 
#'   v.layout <- layoutVertexByAttr(ex_kgml_sig, "pathway", cluster.strength=5)
#'   plotNetwork(ex_kgml_sig, vertex.color="pathway", layout=v.layout)
#' 
layoutVertexByAttr <- function(graph, attr.name, cluster.strength = 1,layout=layout.auto){    
    if(!is.function(layout))
        stop("'layout' is not a function.")
    
    g <- graph.edgelist(get.edgelist(graph)) # create a lightweight graph w/o the attributes.
    V(g)$name <- as.character(1:vcount(g))
    E(g)$weight <- 1
    
    attr <- getAttribute(graph, attr.name)
    attr <- do.call("rbind", lapply(1:length(attr), 
                    function(i) cbind(id=i, val=if(length(attr[[i]])==0) NA else attr[[i]] )))
    
    attr <- attr[!is.na(attr[,2]),] #Remove vertices with no attribute values        
    
    g <- g + vertices(unique(attr[,2])) + igraph::edges(unlist(t(attr)), weight=cluster.strength)
    l <- layout(g, weights=E(g)$weight)[1:vcount(graph),]
    return(l)    
}

#' Computes colors for vertices according to their attributes.
#' 
#' This function returns a list of colors for vertices, assigned similar colors if they share a common attribute
#' (ex: in the same pathway, etc).
#' 
#' @param graph An annotated igraph object.
#' @param attr.name The attribute name (ex: "pathway") by which vertices will be colored. 
#' Complex attributes, where a vertex belongs to more than one group, are supported. 
#' @param col.palette A color palette, or a palette generating function (ex: \preformatted{col.palette=rainbow}).
#' 
#' @return A list of colors (in HEX format) for vertices. 
#' 
#' @author Ahmed Mohamed
#' @family Plotting methods
#' @export
#' @examples 
#'   data("ex_kgml_sig")
#'   v.colors <- colorVertexByAttr(ex_kgml_sig, "pathway") 
#'   plotNetwork(ex_kgml_sig, vertex.color=v.colors)
#'
colorVertexByAttr <- function(graph, attr.name, col.palette = palette()){
    attr <- getAttribute(graph, attr.name)
    attr <- do.call("rbind", lapply(1:length(attr), 
                    function(i) data.frame(id=i, val=if(length(attr[[i]])==0) NA else attr[[i]] )))
    
    attr$val <- as.factor(attr$val)
    
    if(!is.function(col.palette))
        col.palette <- colorRampPalette(col.palette)
    
    col <- col.palette(nlevels(attr$val) + 1)
    names(col) <- c(levels(attr$val), "NA_value_color")

    col.vector <- col[attr$val]
    col.vector[is.na(col.vector)] <- col[["NA_value_color"]] 
    
    return( split(col.vector, as.numeric(attr$id)) )
}

#' Plots an annotated igraph object in Cytoscape.
#'
#' \code{plotCytoscape} function has been removed because RCytoscape is no longer prensent in Bioconductor.
#' Future plans will use RCy3 for Cytoscape plotting, once RCy3 is supported on MacOS and Windows.    
#' \link{plotCytoscapeGML} exports the network plot in GML format, that can be later imported into Cytoscape
#' (using "import network from file" option). This fuction is compatible with all Cytoscape versions.
#'
#' @param graph An annotated igraph object.
#' @param file Output GML file name to which the network plot is exported.
#' @param layout Either a graph layout function, or a two-column matrix specifiying vertex coordinates.
#' @param vertex.size Vertex size. If missing, the vertex attribute "size" (\preformatted{V(g)$size)}) will be used.
#' @param vertex.label Vertex labels. If missing, the vertex attribute "label" (\preformatted{V(g)$label)}) will be used.
#' If missing, vertices are labeled by their name.
#' @param vertex.shape Vertex shape in one of igraph shapes. If missing, the vertex attribute "shape" (\preformatted{V(g)$shape)})
#' will be used. Shapes are converted from igraph convention to Cytoscape convention. "square","rectangle" and "vrectangle" are
#' converted to "RECT",  "csquare" and "crectangle" are converted to "ROUND_RECT", all other shapes are considered "ELLIPSE"
#' @param vertex.color A color or a list of colors for vertices. Vetices with multiple colors are not
#' supported. If missing, the vertex attribute "color" (\preformatted{V(g)$color)}) will be used.
#' @param edge.color A color or a list of colors for edges. If missing, the edge attribute "color"
#' (\preformatted{E(g)$color)}) will be used.
#'
#' @return For \code{plotCytoscapeGML}, results are written to file.
#' 
#' @export
#' @author Ahmed Mohamed
#' @family Plotting methods
#' @rdname plotCytoscape
#' @examples
#'  data("ex_sbml")
#' 	rgraph <- makeReactionNetwork(ex_sbml, simplify=TRUE)
#'  v.layout <- layoutVertexByAttr(rgraph, "compartment") 
#' 	v.color <- colorVertexByAttr(rgraph, "compartment")
#'  
#'  # Export network plot to GML file
#'  plotCytoscapeGML(rgraph, file="example.gml", layout=v.layout, 
#' 				vertex.color=v.color, vertex.size=10)
#' 
plotCytoscapeGML <- function(graph, file, layout=layout.auto, 
        vertex.size, vertex.label, vertex.shape, vertex.color, edge.color){
    
    #v.size
    if(!missing(vertex.size)){
        if(length(vertex.size)==1)
            vertex.size <- rep(vertex.size, vcount(graph))
        
        vertex.size <- as.integer(vertex.size)
        if(length(vertex.size) != vcount(graph)){
            warning("Vertex sizes length and number of vertices don't match")
        }else{
            V(graph)$size <- vertex.size
        }
    }else if(is.null(V(graph)$size)){
        V(graph)$size <- 15
    }
    
    #v.label
    if(!missing(vertex.label)){
        if(length(vertex.label) == vcount(graph)){
            V(graph)$label <- as.character(vertex.label)
        }else
            warning("Vertex lebels length and number of vertices don't match")        
    }
    
    #v.shape
    igraphShape2Cyto <- function(x){
        return(ifelse( x %in% c("square","rectangle","vrectangle"), "rect",
                        ifelse(x %in% c("csquare","crectangle"), "round_rect", "ellipse" ))
        )
    }
    if(!missing(vertex.shape)){
        vertex.shape <- as.character(vertex.shape)
        if(length(vertex.shape)==1)
            vertex.shape <- rep(vertex.shape, vcount(graph))
        
        if(length(vertex.shape) == vcount(graph)){
            V(graph)$type <- vertex.shape
        }else
            warning("Vertex shapes length and number of vertices don't match")
        
    }else if(!is.null(V(graph)$shape)){
        V(graph)$type <- igraphShape2Cyto(V(graph)$shape)
    }else{
        V(graph)$type <- "ellipse"
    }
    
    #v.color
    col2hex <- function(x) return( rgb(t(col2rgb( x )/255)) )
    if(!missing(vertex.color)){
        vertex.color <- as.character(vertex.color)
        if(length(vertex.color)==1)
            vertex.color <- rep(vertex.color, vcount(graph))        
        
        if(length(vertex.color) == vcount(graph)){
            V(graph)$fill <- col2hex(vertex.color)
        }else{
            warning("Vertex colors length and number of vertices don't match.\n
                            Multi-colored vertices are not supported in Cytoscape plots")
        }
    }else if(!is.null(V(graph)$color)){
        V(graph)$fill <- col2hex(V(graph)$color)
    }else{
        V(graph)$fill <- col2hex("skyblue")		
    }
    
    #e.color
    if(!missing(edge.color)){
        edge.color <- col2hex(edge.color)
        if(length(edge.color)==1)
            edge.color <- rep(edge.color, ecount(graph))
        
        if(length(edge.color) == ecount(graph)){
            E(graph)$fill <- col2hex(edge.color)
        }else
            warning("Edge colors length and number of vertices don't match")
        
    }else if(!is.null(E(graph)$color)){
        E(graph)$fill <- col2hex(E(graph)$color)
    }else{
        E(graph)$fill <- col2hex("grey")
    }
    
    if(!is.null(layout)){
        if(is.function(layout)){
            l <- layout(graph) #* max(V(graph)$size)
            V(graph)$x <- l[,1]
            V(graph)$y <- l[,2]
        }else if(ncol(layout)==2 && nrow(layout)==vcount(graph)){
            V(graph)$x <- layout[,1]
            V(graph)$y <- layout[,2]
        }else{
            warning("Incompatible layout dimensions. It should be 2 columns and rows matching number of vertices")
        }		
    }
    
    attr <- get.data.frame(graph, what="both")
    node.attr <- attr$vertices[!sapply(attr$vertices, is.list)]
    node.graphics <- node.attr[ ,colnames(node.attr) %in% c("x", "y", "size", "type", "fill"), drop=FALSE]
    
    node.attr <- node.attr[ ,!colnames(node.attr) %in% c("x", "y", "size", "type", "fill", "shape", "color"), drop=FALSE]
    node.attr <- cbind(id=as.numeric(V(graph)), node.attr, graphics=toGML(node.graphics, "graphics"))
    nodes.gml <- toGML(node.attr, "node", "\n")
    
    
    edge.attr <- attr$edges[!sapply(attr$edges, is.list)]
    names(edge.attr)[1:2] <- c("source", "target")
    edge.attr[1:2] <- get.edgelist(graph,names=FALSE)
    edge.graphics <- edge.attr[ ,colnames(edge.attr) %in% c("size", "type", "fill"), drop=FALSE]
    
    edge.attr <- edge.attr[ ,!colnames(edge.attr) %in% c("size", "type", "fill", "shape", "color"), drop=FALSE]
    edge.attr <- cbind(edge.attr, graphics=toGML(edge.graphics, "graphics"))
    edges.gml <- toGML(edge.attr, "edge", "\n")
    
    gml.ret <- paste('Creator "NetPathMiner ', packageVersion("NetPathMiner"), 
            ' ', date(),'" \n',
            'graph [\n', nodes.gml, '\n', edges.gml, '\n]',
            sep='')
    
    write(gml.ret, file=file)
}
####################### Internal functions for plotting #######################

plotNetwork_internal <- function(graph, vertex.color, col.palette = palette(), layout = layout.auto, legend,...){
    params <- list(x=graph, ...)
    legend.info <- list()
    gsizes <- graph.sizes(vcount(graph))
    
    if(is.null(params$vertex.size))
        if(is.null(V(graph)$size))
            params$vertex.size = gsizes$vsize            
    
    if(is.null(params$edge.arrow.size))
        if(is.null(E(graph)$arrow.size))
            params$edge.arrow.size = gsizes$earrow
    
    if(is.null(params$vertex.label.cex))
        if(is.null(V(graph)$label.cex))
            params$vertex.label.cex = gsizes$label
    
    if(is.null(params$vertex.frame.color))
        if(is.null(V(graph)$frame.color))
            params$vertex.frame.color = NA
    
    if(!missing(vertex.color))
        params$vertex.color = vertex.color
    else if(!is.null(V(graph)$color))
        params$vertex.color = V(graph)$color
    
    if(!is.null(params$vertex.color)){
        if(length(params$vertex.color)==1 && params$vertex.color %in% getAttrNames(graph)){ # If vertex color is an attribute.
            attr.name <- params$vertex.color
            params$vertex.color <- colorVertexByAttr(graph, attr.name, col.palette)
            
            #Store lenged info
            legend.info <- legend.info <- unlist(setNames(params$vertex.color, NULL))
            legend.info <- legend.info[!duplicated(legend.info)]
            naidx <- is.na(names(legend.info))
            if(sum(naidx)>0)
                legend.info <- c(setNames(NA,attr.name) ,legend.info[!naidx], setNames(legend.info[naidx], "N/A"))
            else legend.info <- c(setNames(NA,attr.name) ,legend.info)
            #layout vertices by attribute
            if(is.function(layout))
                params$layout <- layoutVertexByAttr(graph, attr.name, layout = layout)
        }
        
        # Vertices with multiple colors
        if(is.list(params$vertex.color) && any(sapply(params$vertex.color, length) > 1) ){
            params$vertex.pie.color <- params$vertex.color
            params$vertex.shape <- "pie"
            pie <- lapply(params$vertex.color, function(x) rep(1, length(x)))
            params$vertex.pie <- pie
            params$vertex.color <- NULL
        }else{
            params$vertex.color <- unlist(params$vertex.color)
        }
    }
    
    if(is.null(params$layout))    
        params$layout <- layout
    
    do.call("plot.igraph", params)
    
    if(legend && length(legend.info)>0){
        drawLegend(vertices = legend.info)
    }else if(length(legend.info)>0) {
        invisible(legend.info)
    }
}

highlightPaths <- function(paths.eids, graph, path.clusters,layout=layout.auto, path.col, ...){
    if(length(paths.eids)==0){
        plotNetwork_internal(graph, layout=layout, legend=FALSE, ...)
        return()
    }
    
    # Prepare the color pallette, and mark groups (paths)
    if(missing(path.col))
        path.col <- process.color(path.col, list(paths=paths.eids, y.labels=""), path.clusters, same.plot=FALSE)    
    
    path.col <- unlist(sapply(1:length(paths.eids), function(x) rep(path.col[[x]], length(paths.eids[[x]])) ))
    
    
    mark.groups <- do.call("rbind", lapply(paths.eids, function(x) get.edges(graph, x) ))
    mark.groups <- split(mark.groups,1:nrow(mark.groups))
    
    legend.vertices <- plotNetwork_internal(graph, layout=layout, legend=FALSE, 
            mark.groups=mark.groups, mark.col=path.col, mark.border=NA, mark.shape=0, mark.expand=1,
            ...)    
    
    return(legend.vertices)
}

plotAllNetworks_internal <- function(paths, metabolic.net, reaction.net, gene.net, clusters, layout, path.col, ...){
    # Plotting metabolic network
    if(!is.null(metabolic.net)){
        mr.eids <- getPathsAsEIDs(paths, metabolic.net)
        mr.layout = layout(metabolic.net)
        
        if(is.null(names(path.col))){
            if( !is.null(names(mr.eids)) ) mr.eids <- unlist(mr.eids, recursive=FALSE)
            highlightPaths(mr.eids, metabolic.net, layout=mr.layout, path.col = path.col, main="Metabolic Network", ...)
        }else
            lapply(names(mr.eids), function(x)
                highlightPaths(mr.eids[[x]], metabolic.net, layout=mr.layout, path.col= path.col[[x]],
                        main=paste("Metabolic Network y=",x), ...))
    }
    
    # Plotting reaction network
    if(!is.null(reaction.net)){        
        r.eids <- getPathsAsEIDs(paths, reaction.net)
        if(!is.null(metabolic.net)){
            reaction.coor = match(V(reaction.net)$name,V(metabolic.net)$name)
            r.layout <-    mr.layout[reaction.coor,]        #use the same layout as metabolic net.
        }else    r.layout <- layout(reaction.net)    
        
        if(is.null(names(path.col))){
            if( !is.null(names(r.eids)) ) r.eids <- unlist(r.eids, recursive=FALSE)    
            highlightPaths(r.eids, reaction.net, layout=r.layout, path.col= path.col, main="Reaction Network", ...)
        }else
            lapply(names(r.eids), function(x)
                        highlightPaths(r.eids[[x]], reaction.net, layout=r.layout, path.col= path.col[[x]],
                                main=paste("Reaction Network y=",x), ...))        
    }
    
    # Plotting gene network
    if(!is.null(gene.net)){
        g.eids <- getPathsAsEIDs(paths, gene.net)
        
        g.layout <- NULL
        if(!is.null(metabolic.net)){    #use metabolic.net layout as reference
            gene.coor <- match( sub("^(.*)##", "",V(gene.net)$name), V(metabolic.net)$name)
            g.layout <- mr.layout
        }else if(!is.null(reaction.net)){    #use reaction.net layout as reference
            gene.coor <- match( sub("^(.*)##", "",V(gene.net)$name), V(reaction.net)$name)
            g.layout <- r.layout
        }
        if(!is.null(g.layout)){
            minx = rep(-100, vcount(gene.net));    maxx = rep(100, vcount(gene.net))    
            miny = rep(-100, vcount(gene.net));    maxy = rep(100, vcount(gene.net))
            minx = g.layout[gene.coor,1]; maxx = g.layout[gene.coor,1]
            miny = g.layout[gene.coor,2]; maxy = g.layout[gene.coor,2]
            g.layout = layout.fruchterman.reingold(gene.net, minx = minx, maxx = maxx, miny = miny, maxy = maxy)            
        }else    g.layout <- layout(gene.net)    
        
        if(is.null(names(path.col))){
            if( !is.null(names(g.eids)) ) g.eids <- unlist(g.eids, recursive=FALSE)
            highlightPaths(g.eids, gene.net, layout=g.layout, path.col= path.col, main="Gene Network", ...)
        }else
            lapply(names(g.eids), function(x)
                        highlightPaths(g.eids[[x]], gene.net, layout=g.layout, path.col= path.col[[x]],
                                main=paste("Gene Network y=",x), ...))        
    }    
}

process.layout <- function(paths, numPlots, clusters, plot.clusters){
    # Setting the layout of the plot
    par(omi=c(0,0,0, par("pin")[1]*0.2), new=!par("page"))
    
    if(length(paths$y.labels)==1){
        if(!is.null(clusters) && plot.clusters){
            if(numPlots >1){layout.mat = c(1:numPlots, rep(numPlots+1, numPlots-1), numPlots+2)}
            else{layout.mat = c(1,1,2,3)}        
            graphics::layout(matrix(layout.mat, 2, byrow=TRUE))
            par(mar=c(0,0,2,0)+0.1)
        }else{ par(mfrow=c(1,numPlots), mar=c(0,0,2,0)+0.1)}    
    }else{
        if(!is.null(clusters) && plot.clusters){
            if(numPlots >1){layout.mat = c(1:numPlots, rep(numPlots+1, numPlots-1), numPlots+2)}
            else{layout.mat = c(1,1,2,3)}        
            graphics::layout(matrix(layout.mat, 2, byrow=TRUE))
            par(mar=c(0,0,2,0)+0.1)
        }else{ par(mfrow=c(length(paths$y.labels),numPlots), mar=c(0,0,2,0)+0.1)}        
    }
}

process.color <- function(col.palette, paths, clusters, same.plot){
    #cat("Entering color"); flush.console()
    numPaths <- ifelse(length(paths$y.labels)>1, length(unlist(paths$paths, recursive=FALSE)), length(paths$paths))
    
    #cat("no problem"); flush.console()
    if(!is.function(col.palette))
        col.palette <- colorRampPalette(col.palette)
        
    legend.info <- list()
    if(same.plot){        # Plotting paths belonging to all y labels on the same net.
        path.col <- col.palette(2)[clusters$y+1]
        path.col <- alpha(path.col, clusters$posterior.probs)
        legend.info <- setNames(c(NA, col.palette(2)), c("Y lables", "y = 0", "y = 1"))
    }else if(is.null(clusters)){        
        if(length(paths$y.labels)>1)
            path.col <- lapply(paths$paths, function(x) alpha(col.palette(length(x)), 0.3))
        else path.col <- alpha(col.palette(length(paths$paths)), 0.3)        
    }else{ 
        # Clusters are provided. Assign path colors by cluster.
        if(!is.null(clusters$y)){
            path.col <- col.palette(clusters$param$M)[clusters$component]
            
            #set transparency as the probability that a path belongs to its corresponding cluster.
            path.col <- alpha(path.col, clusters$path.probabilities[matrix(c(1:length(clusters$component),clusters$component), ncol=2)])
            
            path.label <- sapply(1:length(paths$paths), function(x) rep(names(paths$paths)[[x]], length(paths$paths[[x]])))
            path.col <- split(path.col, unlist(path.label))            
        }else{
            path.col <- col.palette(clusters$param$M)[clusters$labels]
            
            #set transparency as the probability that a path belongs to its corresponding cluster.
            path.col <- alpha(path.col, clusters$h[matrix(c(1:length(clusters$labels),clusters$labels), ncol=2)])
            
            if(length(paths$y.labels)>1){
                path.label <- sapply(1:length(paths$paths), function(x) rep(names(paths$paths)[[x]], length(paths$paths[[x]])))
                path.col <- split(path.col, unlist(path.label))
            }
        }
        # Legend info showing cluster colors.
        legend.info <- setNames(c(NA, col.palette(clusters$param$M)), 
                c( "Path Clusters", paste("M",1:clusters$params$M, sep="") ) )
    }
    
    return(list(path.col, legend.info))
}

graph.sizes <- function(v, devsize = min(par("pin"))){
    vsize <- (devsize*200)/v
    if(vsize < 1) vsize <- 1
    if(vsize > 15) vsize <- 15
    
    earrow <- vsize*0.08
    if(earrow < 0.2) earrow <- 0.2
    if(earrow > 1) earrow <- 1
    
    label <- earrow*1.5
    if(label>1) label <- 1
    return(list(vsize=vsize, earrow=earrow, label=label))
}

drawLegend <- function(vertices, paths){
    target.w <- par()$omi[4]
    
    par(omi=c(0,0,0,0), mar= c(0,0,0,0), new=TRUE)
    par(fig=c(0,1,0,1), new=TRUE)
    plot.new()
    
    legend.info <- character()
    lty <- numeric()
    lwd <- numeric()
    pch <- numeric()
    font <- numeric()
    
    if(!missing(vertices) && length(vertices)>0){
        legend.info <- c(legend.info, vertices)
        lty <- c(lty , rep(NA, length(vertices)) )
        lwd <- c(lwd , rep(NA, length(vertices)) )
        pch <- c(pch , rep(16, length(vertices)) )
        font <- pch <- c(font , 2 ,rep(16, length(vertices)-1) ) # Boldface for legend section.
    }
    if(!missing(paths) && length(paths)>0){
        legend.info <- c(legend.info, paths)
        lty <- c(lty , rep(1, length(paths)) )
        lwd <- c(lwd , rep(5, length(paths)) )
        pch <- c(pch , rep(NA, length(paths)) )
        font <- c(font , 2 ,rep(16, length(paths)-1) ) # Boldface for legend section.
    }
    l.names <- names(legend.info)
    l.cex <- 1
    
    usr2inches <- strwidth("exmpleexmpleexmple", units="inches")/strwidth("exmpleexmpleexmple", units="user")
    lw <-legend(x="topright", legend=l.names, lty=1, cex=l.cex, plot=FALSE)$rect$w    * usr2inches
    
    if(lw > target.w){
        l.cex <- l.cex * target.w/lw 
        if(l.cex < 0.5){
            l.cex <- 0.5
            lw <-legend(x="topright", legend=l.names, lty=1, cex=l.cex, plot=FALSE)$rect$w    * usr2inches
            
            names.w <- floor(max(nchar(l.names)) * target.w/lw)
            l.names <- strWrap(l.names, width=names.w)
            par(lheight=0.7)
        }
    }
    
    legend(x="topright", legend=l.names, lty=lty, lwd=lwd,pch=pch, col=legend.info, 
            title="Legend", text.font=font, cex=l.cex,xpd=NA)
    
} 

toGML <- function(attr, element, collapse=NULL){
    l <- length(attr)
    num <- sapply(attr,is.numeric)
    
    # Formats create sprintf formulas for GML formatting.
    # In general, GML format follow the pattern: attr.name attr.val \n
    # Attribute string values are quoted.
    #
    # Format1 deals with attribute names.
    # Format2 deals with attribute values.
    # Format interwines format1 and format2, and supply it to sprintf method.
    
    format1 <- paste('%',1:l,'$s ', sep='') # attribute names exported as strings.
    format1[names(attr)== "graphics"] <- ''
    format2 <- paste('%', (1:l) + l,'$', ifelse(num, 'g', 's'), sep='')
    format2 <- paste(ifelse(num, '', '"'), 
            format2, 
            ifelse(num, '\n', '"\n'), 
            sep='')
    
    format2[names(attr)== "graphics"] <- paste('%', which(names(attr)== "graphics") + l, '$s\n', sep='')
    
    format <- paste(t(cbind(format1, format2)), collapse='')
    format <- paste(element, '[\n', format, ']')
    
    
    ret <- paste( do.call(function(...) sprintf(format, ...), c(names(attr), attr)) ,collapse=collapse) 
    return(ret)
}

strWrap <- function(s, width=20){
    for(i in 1:(floor( max(nchar(s))/ width )+1) ){
        s <- gsub( paste( "(^| |\n)(.{", width-2, ",", width, "}) ", sep=""),"\\1\\2\n", s , perl=TRUE)
        s <- gsub( paste( "(((?!\n).){", width, "})((?!\n).{1,})", sep=""),
                "\\1\n\\3",s, perl=TRUE) #Check for words longers
    }
    return(s)
}

# Adapted from "scales" package
alpha <- function (colour, alpha = NA){
    col <- col2rgb(colour, TRUE)/255
    if (length(colour) != length(alpha)) {
        if (length(colour) > 1 && length(alpha) > 1) {
            stop("Only one of colour and alpha can be vectorised")
        }
        if (length(colour) > 1) {
            alpha <- rep(alpha, length.out = length(colour))
        }
        else if (length(alpha) > 1) {
            col <- col[, rep(1, length(alpha)), drop = FALSE]
        }
    }
    alpha[is.na(alpha)] <- col[4, ][is.na(alpha)]
    new_col <- rgb(col[1, ], col[2, ], col[3, ], alpha)
    new_col[is.na(colour)] <- NA
    new_col
}
