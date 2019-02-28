###############################################################################
#
# dbExtract.R:     This file contains all functions related to pathway file
# processing into igraph objects.
# author: Ahmed Mohamed <mohamed@kuicr.kyoto-u.ac.jp>
#
# This is released under GPL-2.
#
# Documentation was created using roxygen
#
###############################################################################


check.file <- function(filename, ext=".xml"){
    if(is.null(filename)) stop("No file provided!")
    if(length(filename)==1){
        if (!file.exists(filename)) stop("Cannot find file:", filename)
        if(length(list.files(filename))==0)
            return(filename)

        fileList = file.path(filename,list.files(filename, ext))
        return(check.file(fileList))
    }else{
        return(unlist(sapply(filename, check.file, USE.NAMES=FALSE)));
    }
}

#' Processes KGML files into igraph objects
#'
#' This function takes KGML files as input, and returns either a metabolic or a signaling
#' network as output.
#'
#' Users can specify whether files are processes as metabolic or signaling networks.
#'
#' Metabolic networks are given as bipartite graphs, where metabolites and reactions represent
#' vertex types. This is constructed from <reaction> xml node in KGML file, connecting them
#' to their corresponding substrates and products. Each reaction vertex has \code{genes} attribute,
#' listing all genes associated with the reaction. As a general rule, reactions inherit all annotation
#' attributes of its catalyzig genes.
#'
#' Signaling network have genes as vertices and edges represent interactions, such as activiation / inhibition.
#' Genes participating in successive reactions are also connected. Signaling parsing method processes <ECrel>, <PPrel>
#' and <PCrel> interactions from KGML files.
#'
#' To generate a genome scale network, simply provide a list of files to be parsed, or put all
#' file in a directory, as pass the directory path as \code{filename}
#'
#' @param filename A character vector containing the KGML files to be processed.
#' If a directory path is provided, all *.xml files in it and its subdirectories are included.
#' @param parse.as Whether to process file into a metabolic or a signaling network.
#' @param expand.complexes Split protein complexes into individual gene nodes. This argument is
#' ignored if \code{parse.as="metabolic"}
#' @param verbose Whether to display the progress of the function.
#'
#' @return An igraph object, representing a metbolic or a signaling network.
#' @author Ahmed Mohamed
#' @family Database extraction methods
#' @export
#' @examples
#' if(is.loaded("readkgmlfile")){ # This is false if libxml2 wasn't available at installation.
#'     filename <- system.file("extdata", "hsa00860.xml", package="NetPathMiner")
#'
#'     # Process KGML file as a metabolic network
#'     g <- KGML2igraph(filename)
#'     plotNetwork(g)
#'
#'     # Process KGML file as a signaling network
#'     g <- KGML2igraph(filename, parse.as="signaling", expand.complexes=TRUE)
#'     plotNetwork(g)
#' }
#'
KGML2igraph <- function(filename, parse.as=c("metabolic","signaling"), expand.complexes=FALSE, verbose=TRUE){
    if(!is.loaded("readkgmlfile"))
        stop("KGML2igraph requires libxml2 to be present. Please reinstall NetPathMiner after installing libxml2.")
    if("KGML2igraph" %in% options("NPM_ENV")[[1]]$memory.err)
        stop("KGML2igraph has previously caused a critical memory error. For safety reasons, please restart R.")

    fileList = check.file(filename)
    if(length(fileList)==0) stop("No files to process!")
    # Check the parsing method.
    if(!missing(parse.as)){
        if(parse.as=="signaling"){
            return(KGML_signal(fileList, expand.complexes, verbose))
        }else{
            if(parse.as != "metabolic")
                stop("Unknown parsing method:", parse.as)
        }
    }

    if(verbose) message("Parsing KGML files as metabolic networks")
    # If a directory is provided, all xml files are processed.
    if(length(fileList)==1){
        zkgml <- .Call("readkgmlfile", FILENAME = fileList, VERBOSE=verbose)
    }else{
        zkgml <- unlist(sapply(fileList,
                        function(x) .Call("readkgmlfile", FILENAME = x, VERBOSE=verbose)
            , USE.NAMES=FALSE), recursive=FALSE)

        dup.zkgml <- duplicated(names(zkgml))
        dup.rns <- sapply(zkgml[dup.zkgml], "[[", "miriam.kegg.pathway")
        zkgml <- zkgml[!dup.zkgml] # Remove duplicated reactions.


        if(length(dup.rns)>0){
            dup.rns <- split( unname(dup.rns), names(dup.rns))	# contains pathway info for removed rns.
            zkgml[names(dup.rns)] <- mapply(
                                    function(x, dup){
                                        x$miriam.kegg.pathway <- c(x$miriam.kegg.pathway, dup)
                                        return(x)
                                    }, zkgml[names(dup.rns)], dup.rns,
                                    SIMPLIFY=FALSE)
        }
    }

    if(verbose) message("Files processed succefully. Building the igraph object.")
    ######### Make the reaction substrate network ############
    edges <-do.call("rbind",lapply(zkgml,
                    function(x)expand.grid(x$reactants, x$name, x$reactant.stoichiometry)))

    edges <- rbind(
            do.call("rbind",lapply(zkgml,
                            function(x)expand.grid(x$name,x$products, x$product.stoichiometry)))
            , edges)

    names(edges) <- c("from","to", "stoichiometry")

    ########## Making the iGraph object ###################
    # contructing a bipartite graph and removing loop edges
    graph = simplify(graph.data.frame(edges), edge.attr.comb="first")

    reactions <- V(graph)$name %in% names(zkgml)
    V(graph)$reactions <- reactions
    V(graph)$shape<- ifelse(V(graph)$type==TRUE, "square", "circle")
    V(graph)$color <- ifelse(V(graph)$type==TRUE,"red", "skyblue")

    V(graph)[reactions]$attr <- zkgml[V(graph)[reactions]$name]
    V(graph)[!reactions]$attr = lapply(V(graph)[!reactions]$name, function(x) list(miriam.kegg.compound=x))

    graph$source = "KGML"
    graph$type = "MR.graph"

    return(graph)
}


KGML_signal <- function(fileList, expand.complexes, verbose){
    if(verbose) message("Parsing KGML files as signaling networks")
    zkgml <- .Call("readkgml_sign", FILENAME = fileList,
                EXPAND_COMPLEXES = expand.complexes, VERBOSE=verbose)

    if(verbose) message("Files processed succefully. Building the igraph object.")

    graph = graph.empty() + vertices(names(zkgml[[1]]), attr=zkgml[[1]])
    graph = graph + igraph::edges(zkgml[[2]], attr=zkgml[[3]])

    V(graph)$shape<- "circle"
    V(graph)$color <- "blue"

    graph$source = "KGML"
    graph$type = "S.graph"

    return(graph)
}

#' Processes SBML files into igraph objects
#'
#' This function takes SBML files as input, and returns either a metabolic or a signaling
#' network as output.
#'
#' Users can specify whether files are processes as metabolic or signaling networks.
#'
#' Metabolic networks are given as bipartite graphs, where metabolites and reactions represent
#' vertex types. This is constructed from \code{ListOfReactions} in SBML file, connecting them
#' to their corresponding substrates and products (\code{ListOfSpecies}). Each reaction vertex has \code{genes} attribute,
#' listing all \code{modifiers} of this reaction. As a general rule, reactions inherit all annotation
#' attributes of its catalyzig genes.
#'
#' Signaling network have genes as vertices and edges represent interactions. Since SBML format may
#' represent singling events as \code{reaction}, all species are assumed to be genes (rather than small
#' molecules). For a simple path \code{S0 -> R1 -> S1}, in signaling network, the path will be
#' \code{S0 -> M(R1) -> S1} where \code{M(R1)} is R1 modifier(s). To ditiguish gene species from small
#' molecules, user can provide \code{gene.attr} (for example: \code{miriam.uniprot} or \code{miriam.ncbigene})
#' where only annotated species are considered genes.
#'
#' All annotation attributes written according to MIRIAM guidlines (either \code{urn:miriam:xxx:xxx} or
#' \code{http://identifiers.org/xxx/xxx}) are etxracted by default. Non-conforming attributes can be extracted
#' by specifying \code{miriam.attr}.
#'
#' To generate a genome scale network, simply provide a list of files to be parsed, or put all
#' file in a directory, as pass the directory path as \code{filename}
#'
#' Note: This function requires libSBML installed (Please see the installation instructions in the Vignette).
#' Some SBML level-3 files may requires additional libraries also (An infomative error will be displayed when
#' parsing such files). Please visit \url{http://sbml.org/Documents/Specifications/SBML_Level_3/Packages} for
#' more information.
#'
#' @param filename A character vector containing the SBML files to be processed. If a directory path
#' is provided, all *.xml and *.sbml files in it and its subdirectories are included.
#' @param parse.as Whether to process file into a metabolic or a signaling network.
#' @param miriam.attr A list of annotation attributes to be extracted. If \code{"all"}, then all attibutes
#' written in MIRIAM guidelines (see Details) are extracted (Default). If \code{"none"}, then no attributes
#' are extracted. Otherwise, only attributes matching those specified are extracted.
#' @param gene.attr An attribute to distinguish \code{species} representing genes from those
#' representing small molecules (see Details). Ignored if \code{parse.as="metabolic"}.
#' @param expand.complexes Split protein complexes into individual gene nodes. Ignored if
#' \code{parse.as="metabolic"}, or when \code{gene.attr} is not provided.
#' @param verbose Whether to display the progress of the function.
#'
#' @return An igraph object, representing a metbolic or a signaling network.
#' @author Ahmed Mohamed
#' @family Database extraction methods
#' @export
#' @examples
#' if(is.loaded("readsbmlfile")){ # This is false if libSBML wasn't available at installation.
#'     filename <- system.file("extdata", "porphyrin.sbml", package="NetPathMiner")
#'
#'     # Process SBML file as a metabolic network
#'     g <- SBML2igraph(filename)
#'     plotNetwork(g)
#'
#'     # Process SBML file as a signaling network
#'     g <- SBML2igraph(filename, parse.as="signaling",
#'                     gene.attr="miriam.uniprot",expand.complexes=TRUE)
#'     dev.new()
#'     plotNetwork(g)
#' }
SBML2igraph <- function(filename, parse.as=c("metabolic","signaling"), miriam.attr="all",
                    gene.attr, expand.complexes, verbose=TRUE){
    if(!is.loaded("readsbmlfile"))
        stop("SBML2igraph requires libSBML to be present. Please reinstall NetPathMiner after installing libSBML.")

    if("SBML2igraph" %in% options("NPM_ENV")[[1]]$memory.err)
        stop("SBML2igraph has previously caused a critical memory error. For safety reasons, please restart R.")

    fileList <- check.file(filename, ext=".sbml|.xml")
    if(length(fileList)==0) stop("No files to process!")

    if(!missing(parse.as)){
        if(parse.as=="signaling"){
            return(SBML_signal(fileList, miriam.attr,gene.attr,expand.complexes, verbose))
        }else{
            if(parse.as != "metabolic")
                stop("Unknown parsing method:", parse.as)
        }
    }
    if(verbose) message("Parsing SBML files as metabolic networks")

    if(length(fileList)==1){
        zsbml <- .Call("readsbmlfile", FILENAME = fileList, ATTR_TERMS = miriam.attr, VERBOSE=verbose)
        names(zsbml) <- c("reactions", "species")
    }else{
        zsbml <- sapply(fileList, function(x)
                    .Call("readsbmlfile", FILENAME = x, ATTR_TERMS = miriam.attr, VERBOSE=verbose)
                    , USE.NAMES=FALSE)

        zsbml <- apply(zsbml, 1, function(x) c(unlist(x, recursive=FALSE)))
        names(zsbml) <- c("reactions", "species")

        #Resolve reactions and species particiapting in multiple pathways
        dup.zsbml <- duplicated(names(zsbml$reactions))
        dup.rns <- sapply(zsbml$reactions[dup.zsbml], "[[", "pathway")
        zsbml$reactions <- zsbml$reactions[!dup.zsbml] # Remove duplicated reactions.


        if(length(dup.rns)>0){
            dup.rns <- split( unname(dup.rns), names(dup.rns))	# contains pathway info for removed rns.
            zsbml$reactions[names(dup.rns)] <- mapply(
                    function(x, dup){
                        x$pathway <- c(x$pathway, dup)
                        return(x)
                    }, zsbml$reactions[names(dup.rns)], dup.rns,
                    SIMPLIFY=FALSE)
        }

        zsbml$species <- zsbml$species[!duplicated(names(zsbml$species))]
    }
    if(verbose) message("SBML files processed successfully")

    if(verbose) message("Constructing Metabolic Network")
    ######### Make the reaction substrate network ############
    edges <-do.call("rbind",lapply(1:length(zsbml$reactions),
                    function(x)expand.grid(zsbml$reactions[[x]]$reactants,
                                names(zsbml$reactions[x]),
                                zsbml$reactions[[x]]$reactant.stoichiometry)))

    edges <- rbind(
            do.call("rbind",lapply(1:length(zsbml$reactions),
                            function(x)expand.grid(names(zsbml$reactions[x]),
                                        zsbml$reactions[[x]]$products,
                                        zsbml$reactions[[x]]$product.stoichiometry)))
            , edges)

    names(edges) <- c("from","to", "stoichiometry")

    ########## Making the iGraph object ###################
    # contructing a bipartite graph and removing multiple edges
    graph <- graph.empty() + vertices(names(zsbml$reactions))
    graph <- graph + vertices(names(zsbml$species))
    graph <- graph + igraph::edges(c(t(edges[,1:2])), stoichiometry=edges$stoichiometry)

    graph <- simplify(graph, edge.attr.comb="first")
    V(graph)$attr <- unlist(zsbml, recursive=FALSE, use.names=FALSE)

    # Set graphical attributes
    reactions <- V(graph)$name %in% names(zsbml$reactions)
    V(graph)$reactions <- reactions
    V(graph)$shape<- ifelse(V(graph)$reactions==TRUE, "square", "circle")
    V(graph)$color <- ifelse(V(graph)$reactions==TRUE,"red", "skyblue")

    graph$source = "SBML"
    graph$type = "MR.graph"

    return(graph)
}

SBML_signal <- function(fileList, miriam.attr="all", gene.attr, expand.complexes, verbose){
    if(verbose) message("Parsing SBML files as signaling networks")
    zsbml <- .Call("readsbml_sign", FILENAME = fileList, ATTR_TERMS = miriam.attr, VERBOSE=verbose)

    if(verbose) message("SBML files processed successfully")

    # DeleteReconnect non-gene reactions (translocation, spontaneous)
    # DeleteReconnect all non-gene, max 3 edges (R->s->RN->s->R will be disconnected).
    graph = graph.empty() + vertices(zsbml$vertices, attr=zsbml$attr)
    graph = graph + igraph::edges(zsbml$edges)

    if(missing(gene.attr)){
        if(verbose) message("No attributes distiuishing genes were provided.")
        if(verbose) message("All species from the SBML file are included in the network.")
    }else{
        no.gene = which(lapply(lapply(zsbml$attr, "[[", gene.attr), length)==0)
        if(length(no.gene>0))
            graph = vertexDeleteReconnect(graph, no.gene, reconnect.threshold=3)

        if(!missing(expand.complexes)&& expand.complexes)
            graph = expandComplexes(graph, gene.attr, keep.parent.attr=miriam.attr)
    }

    V(graph)$shape<- "circle"
    V(graph)$color <- "blue"

    graph$source = "SBML"
    graph$type = "S.graph"

    return(graph)
}

# TODO: biopax signaling L2, Biopax metabolic L2 small molecules.
#' Processes BioPAX objects into igraph objects
#'
#' This function takes BioPAX objects (level 2 or 3) as input, and returns either a metabolic or a signaling
#' network as output.
#'
#' This function requires \code{rBiopaxParser} installed.
#'
#' Users can specify whether files are processes as metabolic or signaling networks.
#'
#' Metabolic networks are given as bipartite graphs, where metabolites and reactions represent
#' vertex types. Reactions are constructed from \code{Conversion} classes, connecting them
#' to their corresponding \code{Left}s and \code{Right}s. Each reaction vertex has \code{genes} attribute,
#' listing all \code{Catalysis} relationships of this reaction. As a general rule, reactions inherit all annotation
#' attributes of its catalyzig genes.
#'
#' Signaling network have genes as vertices and edges represent interactions, such as activiation / inhibition.
#' Genes participating in successive reactions are also connected. Signaling interactions are constructed from
#' \code{Control} classes, where edges are drawn from \code{controller} to \code{controlled}.
#'
#' All annotation attributes are exracted from \code{XRefs} associated with the vertices, and are stored according to
#' MIRIAM guidelines (\code{miraim.db}, where db is the database name).
#'
#'
#' @param biopax BioPAX object generated by \code{\link[rBiopaxParser]{readBiopax}}.
#' @param parse.as Whether to process file into a metabolic or a signaling network.
#' @param expand.complexes Split protein complexes into individual gene nodes. Ignored if
#' \code{parse.as="metabolic"}.
#' @param inc.sm.molecules Include small molecules that are participating in signaling events. Ignored if
#' \code{parse.as="metabolic"}.
#' @param verbose Whether to display the progress of the function.
#'
#' @return An igraph object, representing a metbolic or a signaling network.
#' @author Ahmed Mohamed
#' @family Database extraction methods
#' @export
#' @examples
#' if(require(rBiopaxParser)){
#'     data(ex_biopax)
#'     # Process biopax as a metabolic network
#'     g <- biopax2igraph(ex_biopax)
#'     plotNetwork(g)
#'
#'     # Process SBML file as a signaling network
#'     g <- biopax2igraph(ex_biopax, parse.as="signaling", expand.complexes=TRUE)
#' }
biopax2igraph <- function(biopax, parse.as=c("metabolic","signaling"),
                    expand.complexes=FALSE, inc.sm.molecules=FALSE, verbose=TRUE){
    if (!require(rBiopaxParser))
        stop("This functions needs the rBiopaxParser library installed. Check out the installation instructions!")
    if (!("biopax" %in% class(biopax)))
        stop("Error: biopax2igraph: parameter biopax has to be of class biopax.")

    biopax$df <- as.data.frame(biopax$dt)
    if(missing(parse.as) || parse.as=="metabolic"){
        if(biopax$biopaxlevel == 3)
            return(bpMetabolicL3(biopax, verbose))
        else
            return(bpMetabolicL2(biopax, verbose))
    }else if(parse.as=="signaling"){
        if(biopax$biopaxlevel == 3)
            return(bpSignalingL3(biopax, expand.complexes, inc.sm.molecules, verbose))
        else
            stop("Parsing singaling networks from BioPAX level 2 files is not currently supported.")
    }
}

bpMetabolicL3 <- function(biopax, verbose){
    if(verbose) message("Processing BioPAX (level 3) object as a metabolic network", appendLF=FALSE)
    to.df <- function(dt) return(as.data.frame(dt))

    df <- biopax$df
    df$property = tolower(df$property)

    # All conversion events
    lefts <- df[df$property=="left",c(from="property_attr_value", to="id")]
    lefts[,1] <- striph(lefts[,1])
    rights <- df[df$property=="right",c("id","property_attr_value")]
    rights[,2] <- striph(rights[,2])

    # Getting distiction between Metabolic and signaling reactions.
    # Metabolic reactions takes only small molecules as substrates/products.
    # The following lines removes any Biochemical reactions with non-small molecules
    #   participants.
    classes <- df[ rBiopaxParser::selectInstances(biopax, lefts[,1], returnValues=FALSE), c("class", "id")]
    sig.reactions <- lefts[match(classes$id, lefts[,1])[classes$class !="SmallMolecule"],2]
    classes <- df[ rBiopaxParser::selectInstances(biopax, rights[,2], returnValues=FALSE), c("class", "id")]
    sig.reactions <- unlist(list(sig.reactions,
                    rights[ match(classes$id, rights[,2])[classes$class !="SmallMolecule"], 1]))

    lefts <- lefts[ !lefts[,2] %in% sig.reactions, ]
    rights <- rights[ !rights[,1] %in% sig.reactions, ]

    # Putting the edges into an igraph object.
    graph <- graph.data.frame( rbind(lefts, setNames(rights, names(lefts))) )
    reactions <- V(graph)$name %in% unique(as.character(lefts[,2], rights[,1]))

    if(sum(reactions)==0)
        stop("No metabolic reactions found. Try parsing the file as a signaling netowrk.")
    if(verbose) message(": ", sum(reactions), " reactions found.")

    XRefs <- bpGetReferences(biopax, V(graph)$name)
    names(XRefs) <- V(graph)$name
    cat.r <- striph(to.df(rBiopaxParser::selectInstances(biopax, class="Catalysis", property="controlled"))$property_attr_value)
    cat.gene <- striph(to.df(rBiopaxParser::selectInstances(biopax, class="Catalysis", property="controller"))$property_attr_value)

    gene.xref <- bpGetReferences(biopax, cat.gene)
    names(gene.xref) <- cat.r
    XRefs[names(gene.xref)] <- mapply(c, XRefs[names(gene.xref)], gene.xref, SIMPLIFY=FALSE)

    ##############################################################
    # Getting attribute lists (name, compartment, pathway, MIRIAM)
    # For reactions (name, reactants, reactant.stoichiometry, products, product.stoichiomentry,
    # reversible, kinetics, genes, pathway)

    # MIRIAM attributes
    attr <- bpGetAnnFromXRef(df, XRefs[V(graph)$name])

    # Name attributes
    v.names <- to.df(rBiopaxParser::listInstances(biopax, id=V(graph)$name))
    v.names <- v.names[match(V(graph)$name, v.names$id), "name"]

    # Pathway attributes##########################
    ## Get Pathway name and annotations
    pw <- to.df(rBiopaxParser::listInstances(biopax, class="pathway"))
    pwXRef <- bpGetReferences(biopax, pw$id)
    #pw.ann <- bpGetAnnFromXRef(df, pwXRef)

    ## Get pathway components (only reactions are returned)
    pwcomp <- lapply(pw$id, function(x)rBiopaxParser::listPathwayComponents(biopax,x, returnIDonly=TRUE, includeSubPathways = FALSE))
    pwcomp <- do.call("rbind", lapply(1:length(pwcomp),
                    function(x)data.frame(id=x, comp=pwcomp[[x]])))

    ## Map pathways to reaction vertices.
    pwcomp <- pwcomp[pwcomp$comp %in% V(graph)[reactions]$name, ]
    pwcomp <- split(pwcomp$id, pwcomp$comp, drop=TRUE)

    ## For metabolites, pathway attributes are inherited from reactions they participate in.
    leftright = rbind(lefts, rights)  ## Edges connecting metaboites to reactions
    leftright$id = match(leftright$id, names(pwcomp)) ## Use numerical ids.

    pwcomp <-c(pwcomp,
            lapply(split(pwcomp[leftright$id], leftright$property_attr_value),
                    function(x) unique(unlist(x)) )
    ) # Do the magic :)

    pwcomp <- pwcomp[match(V(graph)$name, names(pwcomp))]  # Reorder to match V(graph)$name

    pw.ann <- bpGetAnnFromXRef(df,lapply(pwcomp, function(x) unlist(pwXRef[x], use.names=FALSE)))
    ###########################################

    # Compartment attributes #######################
    ## Metabolites have compaertment attributes,
    ## while reactions inherit them from their catalysts.

    ## Get Compartment name and annotations
    comp <- rBiopaxParser::listInstances(biopax, class="cellularLocationvocabulary", returnIDonly=TRUE)
    comp.terms <- as.character(bpGetAttrbyID(df, comp, "term")$property_value)
    #compXRef <- bpGetReferences(biopax, comp)
    comp.ann <- bpGetAnnFromXRef(df, bpGetReferences(biopax, comp) )

    ## Get Metbolite compartments
    met.loc <- bpGetAttrbyID(df, V(graph)$name, "cellularlocation", "property_attr_value")
    met.loc$property_attr_value <- match(striph(met.loc$property_attr_value), comp)
    met.loc <- split(met.loc$property_attr_value, met.loc$id, drop=TRUE)

    ## Get Catalyst compartments
    cat.loc <- bpGetAttrbyID(df, cat.gene, "cellularlocation", "property_attr_value")
    cat.loc$property_attr_value <- match(striph(cat.loc$property_attr_value), comp)
    cat.loc$id <- cat.r[match( cat.loc$id, cat.gene )]
    cat.loc <- split(cat.loc$property_attr_value, cat.loc$id, drop=TRUE)

    ## Combine reaction and metbolite comparments, and reorder
    loc <- c(met.loc,cat.loc)
    loc <- loc[ match(V(graph)$name, names(loc)) ]
    ###############################################

    # Reaction-Specific attributes#################
    ##Reactants
    reactants <- split(lefts[,1], lefts[,2], drop=TRUE)
    ##Products
    products <- split(rights[,2], rights[,1], drop=TRUE)
    ##Genes
    cat.gene.name <- to.df(rBiopaxParser::listInstances(biopax,cat.gene))
    cat.gene.name <- cat.gene.name[match(cat.gene,cat.gene.name$id),"name"]
    genes <- split(cat.gene.name, cat.r, drop=TRUE)

    ## Stoichiometry
    st <- bpGetAttrbyID(df, V(graph)[reactions]$name, "participantstoichiometry", "property_attr_value")
    st.met <- bpGetAttrbyID(df, st$property_attr_value, "physicalentity", "property_attr_value")
    st$met <- striph(st.met[match(striph(st$property_attr_value),st.met$id),3])
    st.cf <- bpGetAttrbyID(df, st$property_attr_value, "stoichiometriccoefficient")
    st$cf <- st.cf[match(striph(st$property_attr_value),st.cf$id),3]

    edges <- rbind(lefts,rights)
    edges<- cbind(edges, st=NA)
    for(i in 1:nrow(st)){
        edges[ edges[,1]==st[i,4] & edges[,2]==st[i,1], "st"] <- as.character(st[i,5])
    }
    edges$st <- as.numeric(edges$st)

    r.stoic <- split(edges[1:nrow(lefts),3], edges[1:nrow(lefts),2], drop=TRUE)
    p.stoic <- split(edges[-c(1:nrow(lefts)),3], edges[-c(1:nrow(lefts)),2], drop=TRUE)

    # Reorder our attributes
    reactants <- reactants[ match(V(graph)[reactions]$name, names(reactants)) ]
    products <- products[match(V(graph)[reactions]$name, names(products))]
    genes <- genes[match(V(graph)[reactions]$name, names(genes))]
    r.stoic <- r.stoic[match(V(graph)[reactions]$name, names(r.stoic))]
    p.stoic <- p.stoic[match(V(graph)[reactions]$name, names(p.stoic))]
    #######################################################
    # Putting it together

    ## A small hack to distinguish compartment annotations from vertex annotation.
    names(comp.ann) <- rep("compartment", length(comp.ann))

    ## For reactions
    V(graph)[reactions]$attr <-
            mapply(function(...){
                        args=list(...)
                        c(name=args[[1]],
                            reversible=FALSE,
                            reactants = list(args[[2]]),
                            reactant.stoichiometry= list(args[[3]]),
                            products = list(args[[4]]),
                            product.stoichiometry= list(args[[5]]),
                            kinetics=NULL,
                            genes = list(args[[6]]),
                            compartment=list(comp.terms[args[[7]] ]),
                            unlist(comp.ann[args[[7]] ], recursive=FALSE),
                            pathway=list(pw$name[args[[8]] ]),
                            unlist(args[9],recursive=FALSE),
                            args[[10]]
                        )
                    }, v.names[reactions],
                    reactants, r.stoic, products, p.stoic, genes,
                    loc[reactions], pwcomp[reactions], pathway = pw.ann[reactions],attr[reactions],
                    SIMPLIFY=FALSE)

    V(graph)[!reactions]$attr <-
            mapply(function(...){
                        args=list(...)
                        c(name=args[[1]],
                            compartment=list(comp.terms[args[[2]] ]),
                            unlist(comp.ann[args[[2]] ], recursive=FALSE),
                            pathway=list(pw$name[args[[3]] ]),
                            unlist(args[4],recursive=FALSE),
                            args[[5]]
                        )
                    }, v.names[!reactions],
                    loc[!reactions], pwcomp[!reactions], pathway= pw.ann[!reactions],attr[!reactions],
                    SIMPLIFY=FALSE)

    E(graph)$stoichiometry = edges$st

    V(graph)$reactions <- reactions
    V(graph)$shape<- ifelse(V(graph)$reactions==TRUE, "square", "circle")
    V(graph)$color <- ifelse(V(graph)$reactions==TRUE,"red", "skyblue")

    graph$source = "BioPAX_L3"
    graph$type = "MR.graph"

    return(graph)
}

bpMetabolicL2 <- function(biopax, verbose){
    if(verbose) message("Processing BioPAX (level 2) object as a metabolic network", appendLF=FALSE)
    to.df <- function(dt) return(as.data.frame(dt))

    df <- biopax$df
    df$property = tolower(df$property)

    lefts <- df[df$property=="left",c(from="property_attr_value", to="id")]
    left.mol = bpGetAttrbyID(df,lefts[,1], "physical-entity", "property_attr_value")#
    lefts[,1] = left.mol[match(striph(lefts[,1]), left.mol$id), "property_attr_value"]#
    lefts[,1] <- striph(lefts[,1])

    rights <- df[df$property=="right",c("id","property_attr_value")]
    right.mol = bpGetAttrbyID(df,rights[,2], "physical-entity", "property_attr_value")#
    rights[,2] = right.mol[match(striph(rights[,2]), right.mol$id), "property_attr_value"]#
    rights[,2] <- striph(rights[,2])

    graph <- graph.data.frame( rbind(lefts, setNames(rights, names(lefts))) )
    reactions <- V(graph)$name %in% unique(as.character(lefts[,2], rights[,1]))

    if(verbose) message(": ", sum(reactions), " reactions found.")

    XRefs <- bpGetReferences(biopax, V(graph)$name)
    names(XRefs) <- V(graph)$name
    cat.r <- striph(to.df(rBiopaxParser::selectInstances(biopax, class="Control", property="controlled"))$property_attr_value)
    cat.gene <- striph(to.df(rBiopaxParser::selectInstances(biopax, class="Control", property="controller"))$property_attr_value)

    gene.xref <- bpGetReferences(biopax, cat.gene, followProperties="physical-entity")#
    names(gene.xref) <- cat.r
    XRefs[names(gene.xref)] <- mapply(c, XRefs[names(gene.xref)], gene.xref, SIMPLIFY=FALSE)

    ##############################################################
    # Getting attribute lists (name, compartment, pathway, MIRIAM)
    # For reactions (name, reactants, reactant.stoichiometry, products, product.stoichiomentry,
    # reversible, kinetics, genes, pathway)

    # MIRIAM attributes
    attr <- bpGetAnnFromXRef(df, XRefs[V(graph)$name])

    # Name attributes
    v.names <- to.df(rBiopaxParser::listInstances(biopax, id=V(graph)$name))
    v.names <- v.names[match(V(graph)$name, v.names$id), "name"]

    # Pathway attributes##########################
    ## Get Pathway name and annotations
    pw <- rBiopaxParser::listInstances(biopax, class="pathway")
    pwXRef <- bpGetReferences(biopax, pw$id)
    #pw.ann <- bpGetAnnFromXRef(df, pwXRef)

    ## Get pathway components (only reactions are returned)
    pwcomp <- lapply(pw$id, function(x)rBiopaxParser::listPathwayComponents(biopax,x, returnIDonly=TRUE))
    pwcomp <- do.call("rbind", lapply(1:length(pwcomp),
                    function(x)data.frame(id=x, comp=pwcomp[[x]])))

    ## Map pathways to reaction vertices.
    pwcomp <- pwcomp[pwcomp$comp %in% V(graph)[reactions]$name, ]
    pwcomp <- split(pwcomp$id, pwcomp$comp, drop=TRUE)

    ## For metabolites, pathway attributes are inherited from reactions they participate in.
    leftright = rbind(lefts, rights)  ## Edges connecting metaboites to reactions
    leftright$id = match(leftright$id, names(pwcomp)) ## Use numerical ids.

    pwcomp <-c(pwcomp,
            lapply(split(pwcomp[leftright$id], leftright$property_attr_value),
                    function(x) unique(unlist(x)) )
    ) # Do the magic :)

    pwcomp <- pwcomp[match(V(graph)$name, names(pwcomp))]  # Reorder to match V(graph)$name

    pw.ann <- bpGetAnnFromXRef(df,lapply(pwcomp, function(x) unlist(pwXRef[x], use.names=FALSE)))
    ###########################################

    # Compartment attributes #######################
    ## Metabolites have compaertment attributes,
    ## while reactions inherit them from their catalysts.

    ## Get Pathway name and annotations
    comp <- rBiopaxParser::listInstances(biopax, class="openControlledVocabulary", returnIDonly=TRUE)#
    comp.terms <- as.character(bpGetAttrbyID(df, comp, "term")$property_value)
    compXRef <- bpGetReferences(biopax, comp)#
    #comp.ann <- bpGetAnnFromXRef(df, bpGetReferences(biopax, comp) )

    ## Get Metbolite compartments
    left.loc <- bpGetAttrbyID(df, left.mol$id, "cellular-location", "property_attr_value")#
    left.loc$id <- left.mol[match(left.loc$id, left.mol$id),3]#
    right.loc <- bpGetAttrbyID(df, right.mol$id, "cellular-location", "property_attr_value")#
    right.loc$id <- right.mol[match(right.loc$id, right.mol$id),3]#

    met.loc <- rbind(left.loc, right.loc)#
    met.loc$property_attr_value <- match(striph(met.loc$property_attr_value), comp)
    met.loc <- split(met.loc$property_attr_value, striph(met.loc$id), drop=TRUE)#

    ## Get Catalyst compartments
    cat.loc <- bpGetAttrbyID(df, cat.gene, "cellular-location", "property_attr_value")#
    cat.loc$property_attr_value <- match(striph(cat.loc$property_attr_value), comp)
    cat.loc$id <- cat.r[match( cat.loc$id, cat.gene )]
    cat.loc <- split(cat.loc$property_attr_value, cat.loc$id, drop=TRUE)

    ## Combine reaction and metbolite comparments, and reorder
    loc <- c(met.loc,cat.loc)
    loc <- loc[ match(V(graph)$name, names(loc)) ]
    comp.names <- lapply(loc, function(x)unique(unlist(comp.terms[ x ], use.names=FALSE)))

    locXRef <- lapply(loc, function(x)unlist(compXRef[ x ], use.names=FALSE))
    comp.ann <- bpGetAnnFromXRef(df, locXRef)
    ###############################################

    # Reaction-Specific attributes#################
    ##Reactants
    reactants <- split(lefts[,1], lefts[,2], drop=TRUE)
    ##Products
    products <- split(rights[,2], rights[,1], drop=TRUE)
    ##Genes ##totally different from level3 function
    cat.gene.ref <- bpGetAttrbyID(df,cat.gene, "physical-entity", "property_attr_value")
    cat.gene.name <- to.df(rBiopaxParser::listInstances(biopax, cat.gene.ref[,3]))
    cat.gene.ref$name <- cat.gene.name$name[match(striph(cat.gene.ref[,3]), cat.gene.name$id)]

    cat.gene.ref <- cat.gene.ref[match(cat.gene,cat.gene.ref$id),"name"]
    genes <- split(cat.gene.ref, cat.r, drop=TRUE)

    ## Stoichiometry
    lefts <- df[df$property=="left",c(from="property_attr_value", to="id")]
    rights <- df[df$property=="right",c("id","property_attr_value")]
    edges <- rbind(lefts,rights)
    edges<- cbind(edges, st=NA)

    st <- bpGetAttrbyID(df, edges$property_attr_value, "stoichiometric-coefficient")
    edges$st = st$property_value[match(striph(edges$property_attr_value), st$id)]
    edges$st <- as.numeric(as.character(edges$st))

    r.stoic <- split(edges[1:nrow(lefts),3], edges[1:nrow(lefts),2], drop=TRUE)
    p.stoic <- split(edges[-c(1:nrow(lefts)),3], edges[-c(1:nrow(lefts)),2], drop=TRUE)

    # Reorder our attributes
    reactants <- reactants[ match(V(graph)[reactions]$name, names(reactants)) ]
    products <- products[match(V(graph)[reactions]$name, names(products))]
    genes <- genes[match(V(graph)[reactions]$name, names(genes))]
    r.stoic <- r.stoic[match(V(graph)[reactions]$name, names(r.stoic))]
    p.stoic <- p.stoic[match(V(graph)[reactions]$name, names(p.stoic))]
    #######################################################
    # Putting it together

    ## A small hack to distinguish compartment annotations from vertex annotation.
    #names(comp.ann) <- rep("compartment", length(comp.ann))

    ## For reactions
    V(graph)[reactions]$attr <-
            mapply(function(...){
                        args=list(...)
                        c(name=args[[1]],
                            reversible=FALSE,
                            reactants = list(args[[2]]),
                            reactant.stoichiometry= list(args[[3]]),
                            products = list(args[[4]]),
                            product.stoichiometry= list(args[[5]]),
                            kinetics=NULL,
                            genes = list(args[[6]]),
                            compartment=list(args[[7]]),
                            unlist(args[8], recursive=FALSE),
                            pathway=list(pw$name[args[[9]] ]),
                            unlist(args[10],recursive=FALSE),
                            args[[11]]
                        )
                    }, v.names[reactions],
                    reactants, r.stoic, products, p.stoic, genes,
                    comp.names[reactions], compartment = comp.ann[reactions],
                    pwcomp[reactions], pathway = pw.ann[reactions],attr[reactions],
                    SIMPLIFY=FALSE)

    V(graph)[!reactions]$attr <-
            mapply(function(...){
                    args=list(...)
                    c(name=args[[1]],
                        compartment=list(args[[2]]),
                        unlist(args[3], recursive=FALSE),
                        pathway=list(pw$name[args[[4]] ]),
                        unlist(args[5],recursive=FALSE),
                        args[[6]]
                    )
                }, v.names[!reactions],
                comp.names[!reactions], compartment = comp.ann[!reactions],
                pwcomp[!reactions], pathway = pw.ann[!reactions],attr[!reactions],
                SIMPLIFY=FALSE)

    E(graph)$stoichiometry = edges$st

    V(graph)$reactions <- reactions
    V(graph)$shape<- ifelse(V(graph)$reactions==TRUE, "square", "circle")
    V(graph)$color <- ifelse(V(graph)$reactions==TRUE,"red", "skyblue")

    graph$source = "BioPAX_L2"
    graph$type = "MR.graph"

    return(graph)
}

# TODO: edge attributes in signaling networks
bpSignalingL3 <- function(biopax, expand.complexes=FALSE, inc.sm.molecules=FALSE, verbose=TRUE){
    if(verbose) message("Processing BioPAX (level 3) object as a signaling network", appendLF=FALSE)
    to.df <- function(dt) return(as.data.frame(dt))
    df <- biopax$df
    df$property = tolower(df$property)

    ############### Constructing the signaling graph ##########################
    # All conversion events
    lefts <- df[df$property=="left",c(from="property_attr_value", to="id")]
    lefts[,1] <- striph(lefts[,1])
    rights <- df[df$property=="right",c("id","property_attr_value")]
    rights[,2] <- striph(rights[,2])

    # Identifying Metabolic and signaling reactions.
    # Metabolic reactions takes only small molecules as substrates/products.
    # Signaling reactions have at least one non-small-molecule participant.
    classes.left <- df[ rBiopaxParser::selectInstances(biopax, lefts[,1], returnValues=FALSE), c("class", "id")]
    sig.reactions <- lefts[match(classes.left$id, lefts[,1])[classes.left$class !="SmallMolecule"],2]
    classes.right <- df[ rBiopaxParser::selectInstances(biopax, rights[,2], returnValues=FALSE), c("class", "id")]
    sig.reactions <- unlist(list(sig.reactions,
                    rights[ match(classes.right$id, rights[,2])[classes.right$class !="SmallMolecule"], 1]))


    controlled <- striph(to.df(rBiopaxParser::selectInstances(biopax, class="Catalysis", property="controlled"))$property_attr_value)
    controller <- striph(to.df(rBiopaxParser::selectInstances(biopax, class="Catalysis", property="controller"))$property_attr_value)


    # Metabolic reactions are represented in signaling networks by their "Controller"
    #   where controllers of successive reactions are connected. If no "Controller"
    #   is assigned, no edges are drawn.
    met.left <- lefts[ !lefts[,2] %in% sig.reactions,]
    met.right <- rights[ !rights[,1] %in% sig.reactions,]
    met.left <- do.call("rbind",
            mapply(bpExpand.grid, split(met.left[,1], met.left[,2], drop=TRUE)[controlled], controller,SIMPLIFY=FALSE)
    )
    met.right <- do.call("rbind",
            mapply(bpExpand.grid, controller, split(met.right[,2], met.right[,1], drop=TRUE)[controlled], SIMPLIFY=FALSE)
    )

    graph <- graph.data.frame( na.omit( rbind(met.left, met.right) ) )
    reactions <- V(graph)$name %in% unique(as.character(met.left[,2], met.right[,1]))
    if(sum(reactions)>0){
        V(graph)$reactions <- reactions
        V(graph)$attr <- list(list())
        graph <- makeReactionNetwork(graph)
    }

    # Signaling reactions are represeted as Left -> Controller -> Right
    # If no "Controller" is assigned, it's represented as Left -> Right
    if(inc.sm.molecules){
        sig.left <- lefts[ lefts[,2] %in% sig.reactions, ]
        sig.right <- rights[ rights[,1] %in% sig.reactions, ]
    }else{
        sig.left <- lefts[ lefts[,1] %in% classes.left[classes.left$class!="SmallMolecule","id"] , ]
        sig.right <- rights[ rights[,2] %in% classes.right[classes.right$class!="SmallMolecule","id"], ]
    }

    # No Controllers (Left -> Right)
    sig.no.ctrl <- as.character(unique( sig.reactions[ !sig.reactions %in% controlled ] ))
    sig.no.ctrl <- do.call("rbind",
            mapply(bpExpand.grid, split(sig.left[,1], sig.left[,2], drop=TRUE)[sig.no.ctrl],
                                  split(sig.right[,2], sig.right[,1], drop=TRUE)[sig.no.ctrl]
                    ,SIMPLIFY=FALSE)
    )

    # Controllers present (Left -> Controller -> Right)
    sig.left <- do.call("rbind",
            mapply(bpExpand.grid, split(sig.left[,1], sig.left[,2], drop=TRUE)[controlled], controller,SIMPLIFY=FALSE)
    )
    sig.right <- do.call("rbind",
            mapply(bpExpand.grid, controller, split(sig.right[,2], sig.right[,1], drop=TRUE)[controlled], SIMPLIFY=FALSE)
    )

    # Putting it all signaling reactions together
    all.sig <- c(t( na.omit(rbind(sig.no.ctrl, sig.left, sig.right)) ))
    sig.vertices <- unique(all.sig)

    # Combining the signaling reactions with the metabolic one.
    graph <- graph + vertices(sig.vertices[!sig.vertices %in% V(graph)$name], attr=list(list())) + igraph::edges(all.sig)

    if(expand.complexes){
        comsp <- bpSplitComplex(biopax,V(graph)$name, inc.sm.molecules)
        graph <- setAttribute(graph,"complex", comsp)
        graph <- igraph::simplify( expandComplexes(graph, "complex", NULL, "normal", "keep") )
    }

    if(ecount(graph)==0)
        stop("No Singling interaction were found.")
    if(verbose) message(": ", ecount(graph), " interaction found.")

    ##################################################################
    ############ Getting attributes for graph vertices ###############
    # Getting attribute lists (name, compartment, pathway, MIRIAM)
    XRefs <- bpGetReferences(biopax, V(graph)$name)
    names(XRefs) <- V(graph)$name

    # MIRIAM attributes
    attr <- bpGetAnnFromXRef(df, XRefs[V(graph)$name])

    # Name attributes
    v.names <- to.df(rBiopaxParser::listInstances(biopax, id=V(graph)$name))
    v.names <- v.names[match(V(graph)$name, v.names$id), "name"]

    ## Get Pathway name and annotations
    pw <- to.df(rBiopaxParser::listInstances(biopax, class="pathway"))
    pwXRef <- bpGetReferences(biopax, pw$id)

    ## Get pathway components (only reactions are returned)
    pwcomp <- lapply(pw$id, function(x) rBiopaxParser::listPathwayComponents(biopax,x, returnIDonly=TRUE))
    pwcomp <- do.call("rbind", lapply(1:length(pwcomp),
                    function(x)data.frame(id=x, comp=pwcomp[[x]])))

    ## Restrict Pathway components to Conversion reactions (Removing PathwayStep, Control).
    pwcomp <- pwcomp[pwcomp$comp %in% rBiopaxParser::listInstances(biopax, class="Conversion", includeSubClasses=TRUE, returnIDonly=TRUE),]

    pwcomp <- split(pwcomp$id, pwcomp$comp, drop=TRUE)

    ## Since Conversions (reactions) are no longer present in the network
    #   (either their controllers or substrate/products are retained), vertices inherit
    #   Pathway annotations from reactions they participate in.
    leftright <- rbind(lefts, rights, cbind(property_attr_value=controller,id=controlled)) # Vertices and reactions they participate in.
    leftright <- leftright[ leftright[,1] %in% V(graph)$name,]
    leftright$id = match(leftright$id, names(pwcomp)) ## Use numerical ids.

    pwcomp <-c(pwcomp,
            lapply(split(pwcomp[leftright$id], leftright$property_attr_value),
                    function(x) unique(unlist(x)) )
    ) # Do the magic :)

    pwcomp <- pwcomp[V(graph)$name]  # Reorder to match V(graph)$name
    pw.ann <- bpGetAnnFromXRef(df,lapply(pwcomp, function(x) unlist(pwXRef[x], use.names=FALSE)))
    ###########################################

    ## Compartment attributes #######################
    ### Get Compartment name and annotations
    comp <- rBiopaxParser::listInstances(biopax, class="cellularLocationvocabulary", returnIDonly=TRUE)
    comp.terms <- as.character(bpGetAttrbyID(df, comp, "term")$property_value)
    comp.ann <- bpGetAnnFromXRef(df, bpGetReferences(biopax, comp) )

    ## Get Vertices' compartments
    loc <- bpGetAttrbyID(df, V(graph)$name, "cellularlocation", "property_attr_value")
    loc$property_attr_value <- match(striph(loc$property_attr_value), comp)
    loc <- split(loc$property_attr_value, loc$id, drop=TRUE)[ V(graph)$name]
    ###############################################
    ###############################################
    # Putting it together

    ## A small hack to distinguish compartment annotations from vertex annotation.
    names(comp.ann) <- rep("compartment", length(comp.ann))

    V(graph)$attr <-
            mapply(function(...){
                        args=list(...)
                        c(name=args[[1]],
                              compartment=list(comp.terms[ args[[2]] ]),
                            unlist(comp.ann[ args[[2]] ], recursive=FALSE),
                            pathway=list(pw$name[ args[[3]] ]),
                            unlist(args[4],recursive=FALSE),
                            args[[5]]
                        )
                    }, v.names, loc, pwcomp, pathway= pw.ann ,attr,
                    SIMPLIFY=FALSE)

    V(graph)$shape<- "circle"
    V(graph)$color <- "blue"

    graph$source = "BioPAX_L3"
    graph$type = "S.graph"

    return(graph)
}
############Helper functions to process biopax objects################
bpGetReferences <- function (biopax, id,
                        followProperties = c("entityreference","component", "memberPhysicalEntity"),
                        getProperties="xref")
{
    id.u = unique(striph(id))

    isrdfresource = biopax$df$property_attr == "rdf:resource"
    propertysel = tolower(biopax$df$property) %in% tolower(followProperties)
    newIDs = biopax$df[biopax$df$id %in% id.u & isrdfresource &
                        propertysel, c("id","property_attr_value")]

    propertyget = tolower(biopax$df$property) %in% tolower(getProperties)
    attr = biopax$df[biopax$df$id %in% id.u & isrdfresource &
                    propertyget, c("id","property_attr_value")]

    attr.ls = split(as.character(attr[,2]), attr$id, drop=TRUE)
    if(nrow(newIDs)>0){
        recur.ls = bpGetReferences(biopax, newIDs$property_attr_value, followProperties, getProperties)
        recur.sp = split(recur.ls,newIDs$id, drop=TRUE)
        keys <- unique(c(names(attr.ls), names(recur.sp)))
        attr.ls = setNames(mapply(function(...)unique(unlist(c(...))),
                        attr.ls[keys], recur.sp[keys], SIMPLIFY=FALSE),
                keys)
    }

    pos = match(striph(id), names(attr.ls))
    return(attr.ls[pos])
}

bpSplitComplex<-function (biopax, complexid, inc.sm.molecules)
{
    compname = c("COMPONENTS", "PHYSICAL-ENTITY")
    if (biopax$biopaxlevel == 3) {
        compname = c("memberPhysicalEntity","component")
    }
    if(inc.sm.molecules)
        classes <- c("dna", "rna", "protein", "smallmolecule")
    else
        classes <- c("dna", "rna", "protein")

    ref = bpGetReferences(biopax, complexid, followProperties=compname, getProperties=compname)

    if (is.null(ref))
        return(NULL)

    ref = do.call("rbind", lapply( na.omit(names(ref)), function(x) cbind( x, striph(ref[[x]]) ) ) )
    referenced = biopax$df[rBiopaxParser::selectInstances(biopax, id = unique(ref[,2]), returnValues=FALSE), c("class", "id")]
    sel = referenced[ tolower(referenced$class) %in% classes, "id"]
    if (length(sel)==0)
        return(NULL)

    ref = ref[ ref[,2] %in% sel, ]
    return( split(ref[,2], ref[,1])[ complexid ] )
}
bpGetAttrbyID <- function(df, id, attr, value="property_value"){
    id = striph(id)
    df[df$id %in% id & df$property %in% attr, c("id","property",value) ]
}
bpGetAnnFromXRef <- function(df, id.ls, attr){
    ann <- do.call("rbind",lapply(na.omit(names(id.ls)),
                            function(x) cbind(x, if(is.null(id.ls[[x]])) NA else id.ls[[x]])
                ))


    id = unique(ann[,2])
    if(length(na.omit(id))==0)
        return(list()[match(names(id.ls), NULL)])

    db <- bpGetAttrbyID(df, id, "db", "property_value")
    dbid <- bpGetAttrbyID(df, id, "id", "property_value")
    dbid <- cbind(db, value=dbid[,3])

    dbid$property_value <- paste("miriam.", dbid$property_value, sep="")

    pos <- match(striph(ann[,2]), dbid$id)
    # Split the dataframe into annotations belonging to each input id (ann[,1])
    # Then split each of the daughter dataframe by the attribute name.
    ann.ls <- lapply(split(dbid[pos,c(3,4)], ann[,1], drop=TRUE),
                    function(x) split(as.character(x[,2]), x[,1], drop=TRUE) )


    return(ann.ls[match(names(id.ls), names(ann.ls))])
}
bpGetAnnotation <- function(df, id, attr){
    ann <- bpGetAttrbyID(df, id, "xref", "property_attr_value")
    ref <- bpGetAttrbyID(df, id, "entityreference", "property_attr_value")
    ref.ann <-  bpGetAttrbyID(df, ref$property_attr_value, "xref", "property_attr_value")

    pos <- match(striph(ref$property_attr_value), ref.ann$id)
    ann <- rbind(ann[,c(1,3)],data.frame(id=ref$id, property_attr_value=ref.ann$property_attr_value[pos],
                    stringsAsFactors = FALSE))

    db <- bpGetAttrbyID(df, ann$property_attr_value, "db", "property_value")
    dbid <- bpGetAttrbyID(df, ann$property_attr_value, "id", "property_value")

    ann.ls <- as.list(as.character(dbid$property_value))
    names(ann.ls) <- paste("miriam.", db$property_value, sep="")
    pos <- match(striph(ann$property_attr_value), dbid$id)
    ann.ls <- split(ann.ls[pos], ann$id, drop=TRUE)


    return(ann.ls[match(striph(id), names(ann.ls))])
}

bpExpand.grid <- function(x1,x2){
    if(is.null(x1))
        x1=NA
    if(is.null(x2))
        x2=NA
    return(expand.grid(x1,x2))
}

striph <- function(s){
    return(sub("^#", "", s))
}
