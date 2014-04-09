### R code from vignette source 'NPMVignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: no.nonsense
###################################################
rm(list=ls())


###################################################
### code chunk number 2: Load_package
###################################################
library(NetPathMiner)


###################################################
### code chunk number 3: NPMVignette.Rnw:235-237 (eval = FALSE)
###################################################
## graph <- KGML2igraph(filename = file)
## graph <- SBML2igraph(filename = file)


###################################################
### code chunk number 4: NPMVignette.Rnw:243-246 (eval = FALSE)
###################################################
## require(rBiopaxParser)
## biopax = readBiopax(file)
## graph <- BioPAX2igraph(biopax = biopax)


###################################################
### code chunk number 5: NPMVignette.Rnw:252-253 (eval = FALSE)
###################################################
## graph <- KGML2igraph(filename = c(file1, file2))


###################################################
### code chunk number 6: NPMVignette.Rnw:257-258 (eval = FALSE)
###################################################
## graph <- KGML2igraph(filename = ".")


###################################################
### code chunk number 7: NPMVignette.Rnw:263-268 (eval = FALSE)
###################################################
## # Extract all MIRIAM identifiers from an SBML file.
## graph <- SBML2igraph(filename = file, miriam = "all")
## 
## # Extract all MIRIAM identifiers from an SBML file.
## graph <- BioPAX2igraph(biopax = biopax, miriam = "go")


###################################################
### code chunk number 8: NPMVignette.Rnw:275-276
###################################################
file <- file.path(find.package("NetPathMiner"), "extdata", "hsa00860.xml")


###################################################
### code chunk number 9: NPMVignette.Rnw:278-282 (eval = FALSE)
###################################################
## graph <- KGML2igraph(filename = file, parse.as = "signaling")
## 
## graph <- KGML2igraph(filename = file, parse.as = "signaling", 
## 	expand.complexes = TRUE)


###################################################
### code chunk number 10: NPMVignette.Rnw:288-291
###################################################
data("ex_sbml")
graph <- ex_sbml
graph


###################################################
### code chunk number 11: NPMVignette.Rnw:300-301
###################################################
head( V(graph) )


###################################################
### code chunk number 12: NPMVignette.Rnw:304-305
###################################################
head( E(graph) )


###################################################
### code chunk number 13: NPMVignette.Rnw:308-309
###################################################
head( V(graph)[ reactions ] )


###################################################
### code chunk number 14: NPMVignette.Rnw:314-315
###################################################
V(graph)[ "reaction_71850" ]$attr


###################################################
### code chunk number 15: NPMVignette.Rnw:322-323
###################################################
getAttrNames(graph)


###################################################
### code chunk number 16: NPMVignette.Rnw:330-331
###################################################
getAttrStatus(graph, pattern = "^miriam.")


###################################################
### code chunk number 17: NPMVignette.Rnw:337-345 (eval = FALSE)
###################################################
## require("RCurl")
## # Fetch uniprot annotation
## graph <- fetchAttribute(graph, organism = "Homo sapiens", 
## target.attr = "miriam.ncbigene" , source.attr = "miriam.uniprot")
## 
## # Fetch ChEBI annotation. 
## graph <- fetchAttribute(graph, target.attr = "miriam.chebi", 
## source.attr = "miriam.kegg.compound")


###################################################
### code chunk number 18: NPMVignette.Rnw:355-357
###################################################
rgraph <- makeReactionNetwork(graph, simplify=FALSE)
rgraph


###################################################
### code chunk number 19: NPMVignette.Rnw:362-364 (eval = FALSE)
###################################################
## rgraph <- simplifyReactionNetwork(rgraph)
## rgraph <- makeReactionNetwork(graph, simplify=TRUE)


###################################################
### code chunk number 20: NPMVignette.Rnw:369-375
###################################################
# Expand complexes of gene network.
ggraph <- expandComplexes(rgraph, v.attr = "miriam.uniprot", 
		keep.parent.attr= c("^pathway", "^compartment"))

# Convert reaction network to gene network.
ggraph <- makeGeneNetwork(rgraph)


###################################################
### code chunk number 21: NPMVignette.Rnw:387-389
###################################################
data(ex_microarray)



###################################################
### code chunk number 22: NPMVignette.Rnw:390-395 (eval = FALSE)
###################################################
## # Assign weights to edges.
## if(require("RCurl") && url.exists( NPMdefaults("bridge.web") ))
## 	rgraph <- fetchAttribute(rgraph, organism = "Homo sapiens", 
## 						target.attr = "miriam.affy.probeset", 
## 						source.attr = "miriam.uniprot")


###################################################
### code chunk number 23: NPMVignette.Rnw:405-409 (eval = FALSE)
###################################################
## library(ALL)
## data(ALL)
## rgraph <- assignEdgeWeights(microarray = exprs(ALL), graph = rgraph,
## weight.method = "cor", use.attr="miriam.affy.probeset", y=ALL$mol.bio, bootstrap = FALSE)


###################################################
### code chunk number 24: NPMVignette.Rnw:413-416
###################################################
data(ex_microarray)
rgraph <- assignEdgeWeights(microarray = ex_microarray, graph = rgraph,
weight.method = "cor", use.attr="miriam.uniprot", y=colnames(ex_microarray), bootstrap = FALSE)


###################################################
### code chunk number 25: NPMVignette.Rnw:420-422
###################################################
rgraph$y.labels
head( E(rgraph)$edge.weights )


###################################################
### code chunk number 26: NPMVignette.Rnw:431-433
###################################################
ranked.p <- pathRanker(rgraph, method = "prob.shortest.path",
	K = 25, minPathSize = 6)


###################################################
### code chunk number 27: NPMVignette.Rnw:438-443 (eval = FALSE)
###################################################
## pathsample <- samplePaths(rgraph, max.path.length = vcount(rgraph),
## num.samples = 1000, num.warmup = 10)
## 
## ranked.p <- pathRanker(rgraph, method = "pvalue", 
## sampledpaths = pathsample ,alpha=0.1)


###################################################
### code chunk number 28: NPMVignette.Rnw:448-450
###################################################
# Get paths as edge IDs.
eids <- getPathsAsEIDs(paths = ranked.p, graph = rgraph)


###################################################
### code chunk number 29: NPMVignette.Rnw:455-457
###################################################
# Convert paths to other networks. 
eids <- getPathsAsEIDs(paths = ranked.p, graph = ggraph)


###################################################
### code chunk number 30: NPMVignette.Rnw:464-467
###################################################
# Clustering.
ybinpaths <- pathsToBinary(ranked.p)
p.cluster <- pathCluster(ybinpaths, M = 2)


###################################################
### code chunk number 31: NPMVignette.Rnw:469-470
###################################################
plotClusters(ybinpaths, p.cluster)


###################################################
### code chunk number 32: NPMVignette.Rnw:475-476
###################################################
p.class <- pathClassifier(ybinpaths, target.class = "BCR/ABL", M = 2)


###################################################
### code chunk number 33: NPMVignette.Rnw:478-479 (eval = FALSE)
###################################################
## plotClassifierROC(p.class)


###################################################
### code chunk number 34: NPMVignette.Rnw:487-488
###################################################
plotClusters(ybinpaths, p.class)


###################################################
### code chunk number 35: NPMVignette.Rnw:496-497
###################################################
plotNetwork(rgraph, vertex.color="compartment.name")


###################################################
### code chunk number 36: NPMVignette.Rnw:502-506 (eval = FALSE)
###################################################
## plotPaths(ranked.p, rgraph)
## 
## # With clusters
## plotPaths(ranked.p, graph, path.clusters=p.class)


###################################################
### code chunk number 37: NPMVignette.Rnw:511-513
###################################################
plotAllNetworks(ranked.p, metabolic.net = graph, reaction.net = rgraph,
		path.clusters=p.class, vertex.label = "", vertex.size = 4)


###################################################
### code chunk number 38: NPMVignette.Rnw:518-522 (eval = FALSE)
###################################################
## layout.c <- clusterVertexByAttr(rgraph, "pathway", cluster.strength = 3)
## v.color <- colorVertexByAttr(rgraph, "pathway")
## plotPaths(ranked.p , rgraph, clusters=p.class, 
## 	layout = layout.c, vertex.color = v.color)


###################################################
### code chunk number 39: NPMVignette.Rnw:527-530 (eval = FALSE)
###################################################
## require(RCytoscape)
## cw <- plotCytoscape(graph, "example", layout = layout.c,
## 				vertex.size = 5, vertex.color = v.color)


###################################################
### code chunk number 40: NPMVignette.Rnw:537-538
###################################################
getGeneSets(graph, use.attr="compartment", gene.attr="miriam.uniprot")


###################################################
### code chunk number 41: NPMVignette.Rnw:543-544
###################################################
getGeneSetNetworks(graph, use.attr="compartment")


###################################################
### code chunk number 42: NPMVignette.Rnw:550-551 (eval = FALSE)
###################################################
## graphNEL <- toGraphNEL(graph, export.attr="^miriam.")


