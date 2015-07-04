# TODO: Add comment
# 
# Author: mohamedahmed
###############################################################################
library("NetPathMiner")

# Create the network
filename <- file.path(find.package("NetPathMiner"), "extdata", "hsa00860.xml")
data("ex_sbml")
graph <- ex_sbml
graph

head( V(graph) )
head( E(graph) )
head( V(graph)[ reactions ] ) 
V(graph)[ "reaction_71850" ]$attr

# List annotation attributes.
getAttrNames(graph)

# Attribute status for all MIRIAM annotations.
getAttrStatus(graph, pattern = "^miriam.")

# Fetch uniprot annotation  
if(require("RCurl") && url.exists(NPMdefaults("bridge.web")))
graph <- fetchAttribute(graph, organism = "Homo sapiens", target.attr = "miriam.uniprot", source.attr = "miriam.ncbigene")

# Fetch ChEBI annotation. 
if(require("RCurl") && url.exists(NPMdefaults("bridge.web")))
graph <- fetchAttribute(graph, target.attr = "miriam.chebi", source.attr = "miriam.kegg.compound")

# Convert metablite-reaction network to reaction network.
rgraph <- makeReactionNetwork(graph, simplify=FALSE)
head( V(rgraph) )

# Remove transport and spontaneous reactions.
rgraph <- simplifyReactionNetwork(rgraph)
rgraph <- makeReactionNetwork(graph, simplify=TRUE)

# Convert reaction network to gene network.
ggraph <- makeGeneNetwork(rgraph)

# Exapnd complexes of gene network.
ggraph <- expandComplexes(ggraph, v.attr = "miriam.ncbigene", 
						keep.parent.attr= c("^pathway", "^compartment"))

# Load data from "ALL" package

# Assign weights to edges.
if(require("RCurl") && url.exists(NPMdefaults("bridge.web")))
rgraph <- fetchAttribute(rgraph, organism = "Homo sapiens", 
						target.attr = "miriam.affy.probeset", 
						source.attr = "miriam.uniprot")

getAttrStatus(rgraph, pattern = "miriam.affy.probeset")
rgraph <- assignEdgeWeights(microarray = exprs(ALL), graph = rgraph, 
							use.attr="miriam.affy.probeset", y=ALL$mol.bio)

					
data(ex_microarray)
rgraph <- assignEdgeWeights(microarray = ex_microarray, graph = rgraph,
		weight.method = "cor", use.attr="miriam.uniprot", y=factor(colnames(ex_microarray)), bootstrap = FALSE)

rgraph$y.labels
head( E(rgraph)$edge.weights )

# Rank paths by probabilistic shortest path method
ranked.p <- pathRanker(rgraph, method = "prob.shortest.path", K = 100, minPathSize = 4)

# Rank paths by p-value method
pathsample <- samplePaths(rgraph, max.path.length = vcount(rgraph),
						num.samples = 1000, num.warmup = 10)
ranked.p <- pathRanker(rgraph, method = "pvalue", 
				sampledpaths = pathsample ,alpha=0.1)

# Clustering.
ybinpaths <- pathsToBinary(ranked.p)
p.cluster <- pathCluster(ybinpaths, M = 2)
plotClusters(ybinpaths, p.cluster)

# Classification
p.class <- pathClassifier(ybinpaths, target.class = "BCR/ABL", M = 5)
plotClassifierROC(p.class)
plotClusters(ybinpaths, p.class)

# Get paths as edge IDs.
eids <- getPathsAsEIDs(paths = ranked.p, graph = rgraph)

# Convert paths to other networks. 
eids <- getPathsAsEIDs(paths = ranked.p, graph = ggraph)

plotNetwork(rgraph, vertex.color="compartment.name")

# Plot paths.
plotPaths(ranked.p, rgraph)

# With clusters
plotPaths(ranked.p, graph, path.clusters=p.class)
dev.off()
# Plot All networks.
plotAllNetworks(ranked.p, metabolic.net = graph, reaction.net = rgraph,
			vertex.label = "", vertex.size = 4)
dev.off()

# Color and cluster by attribute
layout.c <- layoutVertexByAttr(rgraph, "pathway", cluster.strength = 2)
v.color <- colorVertexByAttr(rgraph, "pathway")
plotPaths(graph, rgraph, clusters=p.class, layout = layout.c, vertex.color = v.color)

# Plot in Cytoscape
require(RCytoscape)
cw<-plotCytoscape(rgraph, "example", layout = layout.c,
				vertex.size = 5, vertex.color = v.color)
