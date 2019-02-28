###############################################################################
#
# pathCluster.R: This file contains all functions related to clustering paths
# as well as plotting clutserig results.
#
# author: Ahmed Mohamed <mohamed@kuicr.kyoto-u.ac.jp>
#
# This is released under GPL-2.
#
# Documentation was created using roxygen
#
###############################################################################

#' 3M Markov mixture model for clustering pathways
#'
#' 3M Markov mixture model for clustering pathways
#'
#' @param ybinpaths The training paths computed by \code{\link{pathsToBinary}}.
#' @param M The number of clusters.
#' @param iter The maximum number of EM iterations.
#'
#' @return A list with the following items:
#' \item{h}{The posterior probabilities that each path belongs to each cluster.}
#' \item{labels}{The cluster membership labels.}
#' \item{theta}{The probabilities of each gene for each cluster.}
#' \item{proportions}{The mixing proportions of each path.}
#' \item{likelihood}{The likelihood convergence history.}
#' \item{params}{The specific parameters used.}
#'
#' @references Mamitsuka, H., Okuno, Y., and Yamaguchi, A. 2003. Mining biologically active patterns in
#' metabolic pathways using microarray expression profiles. SIGKDD Explor. News l. 5, 2 (Dec. 2003), 113-121.
#'
#' @author Ichigaku Takigawa
#' @author Timothy Hancock
#' @family Path clustering & classification methods
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
#' 		weight.method = "cor", use.attr="miriam.uniprot", bootstrap = FALSE)
#'
#' 	## Get ranked paths using probabilistic shortest paths.
#'  ranked.p <- pathRanker(rgraph, method="prob.shortest.path",
#' 					K=20, minPathSize=8)
#'
#' 	## Convert paths to binary matrix.
#' 	ybinpaths <- pathsToBinary(ranked.p)
#' 	p.cluster <- pathCluster(ybinpaths, M=2)
#' 	plotClusters(ybinpaths, p.cluster)
#'
pathCluster <- function(ybinpaths, M, iter=1000) {
  x <- ybinpaths$paths

  # remove constant columns
  varying.cols <- which(sapply(x,sd) != 0)
  tr.x <- x[varying.cols]
  if(ncol(tr.x) <= M)
	  stop("Specified number of clusters ", M ,"is larger than varaible genes.",ncol(tr.x),"\n Choose a smaller M")

  zcluster <- sample(1:M,nrow(tr.x),replace = TRUE)
  ptheta <- as.matrix(aggregate(tr.x,list(zcluster),mean)[-1])

  pk <- rep(1/M,M)

  pmx <- matrix(0,nrow(tr.x),M)
  for (k in 1:M) {
    res <- as.matrix(tr.x) %*% diag(as.double(ptheta[k,]))
	  res[tr.x == 0] <- NA
	  pmx[,k] <- apply(res,1,prod,na.rm = TRUE)
  }

  pkm <- matrix(pk,nrow=nrow(tr.x),ncol = M,byrow = TRUE)
  hij <- pkm*pmx/rowSums(pkm*pmx)

  fit <- .C("pathMix",
    X = as.integer(as.matrix(tr.x)),
    M = as.integer(M),
    NOBS = as.integer(nrow(tr.x)),
    NX = as.integer(ncol(tr.x)),
    ITER = as.integer(iter),
    H = as.double(hij),
    THETA = as.double(t(ptheta)),
    PROPORTIONS = as.double(pk),
    LIKELIHOOD = double(iter))

  posterior.probs = data.frame(matrix(fit$H,ncol = M))
  names(posterior.probs) <- paste("M",1:M,sep = "")

  theta <- matrix(NA,nrow = M,ncol = ncol(x))
  t.complete <- matrix(as.double(fit$THETA),nrow = M,ncol = ncol(tr.x),byrow = TRUE)
  theta[,varying.cols] <- t.complete
  theta[,c(1,ncol(theta))] <- 1
  theta <- data.frame(theta)
  names(theta) <- names(x)

  ll <- fit$LIKELIHOOD[1:fit$ITER]

  # cluster labels
  zm <- apply(posterior.probs,1,max)
  cl <- apply(posterior.probs == zm,1,which)
  if (is.list(cl)) {
        cat("\nMultiple cluster labels for some paths, some clusters might not be valid\n")
        clusters <-  paste("M",sapply(cl,"[[",1),sep = "")
  } else clusters <- cl

  return(list(h = posterior.probs,
              labels = clusters,
              theta = theta,
              proportions = fit$PROPORTIONS,
              likelihood = ll,
              params = list(M = M)))
}

#' Predicts new paths given a pathCluster model
#'
#' Predicts new paths given a pathCluster model.
#'
#' @param pfit The pathway cluster model trained by \code{\link{pathCluster}} or \code{\link{pathClassifier}}.
#' @param newdata The binary pathway dataset to be assigned a cluster label.
#'
#' @return A list with the following elements:
#' \tabular{ll}{
#' \code{labels} \tab a vector indicating the 3M cluster membership. \cr
#' \code{posterior.probs} \tab a matrix of posterior probabilities for each path belonging to each cluster.
#' }
#'
#' @author Ichigaku Takigawa
#' @author Timothy Hancock
#' @family Path clustering & classification methods
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
#' 		weight.method = "cor", use.attr="miriam.uniprot", bootstrap = FALSE)
#'
#' 	## Get ranked paths using probabilistic shortest paths.
#'  ranked.p <- pathRanker(rgraph, method="prob.shortest.path",
#' 					K=20, minPathSize=8)
#'
#' 	## Convert paths to binary matrix.
#' 	ybinpaths <- pathsToBinary(ranked.p)
#' 	p.cluster <- pathCluster(ybinpaths, M=2)
#'
#' 	## just an example of how to predict cluster membership.
#' 	pclust.pred <- predictPathCluster(p.cluster,ybinpaths$paths)
#'
predictPathCluster <- function(pfit,newdata) {
  pmx <- matrix(0,nrow=nrow(newdata),ncol = nrow(pfit$theta))
  tt <- pfit$theta
  tt[is.na(tt)] <- 1

  for (k in 1:nrow(pfit$theta)) {
	  res <- as.matrix(newdata) %*% diag(as.numeric(tt[k,]))
	  res[newdata == 0] <- 1
	  pmx[,k] <- pfit$proportions[k] * apply(res,1,prod)
  }
  pmx <- data.frame(pmx/rowSums(pmx))
  names(pmx) <- paste("M",1:ncol(pmx),sep = "")

  # cluster labels
  zm <- apply(pmx,1,max)
  cl <- apply(pmx == zm,1,which)
  if (is.list(cl)) {
        cat("\nMultiple cluster labels for some paths: some clusters might not be valid\n")
        clusters <-  sapply(cl,"[[",1) # just take the first because they are all just as likely!
  } else clusters <- cl

  return(list(labels = clusters,h = pmx))
}


#' Plots the structure of specified path cluster
#'
#' Plots the structure of specified path found by pathCluster.
#'
#' @param ybinpaths The training paths computed by \code{\link{pathsToBinary}}.
#' @param clusters The pathway cluster model trained by \code{\link{pathCluster}} or \code{\link{pathClassifier}}.
#' @param m The path cluster to view.
#' @param tol A tolerance for 3M parameter \code{theta} which is the probability for
#' each edge within each cluster.  If the tolerance is set all edges with a \code{theta}
#' below that tolerance will be removed from the plot.
#'
#' @return Produces a plot of the paths with the path probabilities and cluster membership probabilities.
#' \item{Center Plot}{An image of all paths the training dataset. Rows are the paths and columns are the genes
#'  (features) included within each path.}
#' \item{Right}{The training set posterior probabilities for each path belonging to the current 3M component.}
#' \item{Top Bar Plots}{\code{Theta}, The 3M component probabilities - indicates the importance of each edge to a pathway.}
#'
#' @author Timothy Hancock and Ichigaku Takigawa
#' @family Path clustering & classification methods
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
#' 		weight.method = "cor", use.attr="miriam.uniprot", bootstrap = FALSE)
#'
#' 	## Get ranked paths using probabilistic shortest paths.
#'  ranked.p <- pathRanker(rgraph, method="prob.shortest.path",
#' 					K=20, minPathSize=8)
#'
#' 	## Convert paths to binary matrix.
#' 	ybinpaths <- pathsToBinary(ranked.p)
#' 	p.cluster <- pathCluster(ybinpaths, M=2)
#' 	plotPathCluster(ybinpaths, p.cluster, m=2, tol=0.05)
#'
plotPathCluster <- function(ybinpaths, clusters, m, tol = NULL) {
	if(m > clusters$params$M)
		stop("There are only ", clusters$params$M, " clusters!")
	pp <- clusters$theta[m,]
	fidx <- 1:length(pp)
	if (!is.null(tol)) fidx <- which(pp >= tol | is.na(pp))
	pp <- pp[fidx]

	g <- names(pp)
	x <- ybinpaths$paths[fidx]

	mpar <- par()$mar

	h <- clusters$h[,m]
	sh <- sort(h,index.return = TRUE)

	nl <- layout(matrix(c(1,4,3,2),2,2,byrow=TRUE),
			heights = c(.3,.7,.8,.2),
			widths =  c(.75,.25,.75,.25),TRUE)

	par(mar = c(0,8,4,0),pty = "m")
	xpos <- barplot(height = as.numeric(pp),space = 0,xlim = c(0,ncol(pp)),xaxs = "i",ylim = c(0,1),
			ylab = "Gene\nProbabilities", col = c(1:length(pp))%% 5 + 1)
	abline(h = 0)

	par(mar = c(8,0,0,4))
	plot(x = sh$x,y = 1:nrow(x),xlim = c(0,1),ylim = c(1,nrow(x)),
			xlab = paste("P(M = ",m,"|X)",sep =""),ylab = "",yaxs = "i",type = "l",axes = FALSE)
	abline(v = 0.5,col = 2,lwd = 2)
	abline(v = 0,col = 1,lwd = 2)
	axis(side = 1,at = seq(0,1,1/10),las = 1)
	abline(h = nrow(x)+1)

	hp <- as.factor(as.numeric(h > 0.5))
	ystats <- as.numeric(summary(hp))
	ylab <- c("","","","Not a\nMember","Member")
	ytick <- c(1,nrow(x),ystats[1],ystats[1]/2,ystats[2]/2 + ystats[1])

	par(mar = c(8,8,0,0))
	x[x == 0] <- NA
	x[h < 0.5,] <- x[h < 0.5,] + 1
	px <- x[sh$ix,]
	image(x = xpos,y = 1:nrow(x),as.matrix(t(px)),axes = FALSE,ylab = "",xlab = "",xaxs = "i",yaxs = "i",col = cm.colors(2))
	axis(side = 2,at = ytick,labels = ylab,las = 2)
	mtext("All Paths\nSorted By Component Membership",line = 5,side = 2)

	mtext(g,side = 1,at = xpos,cex = 0.6,las =2,line = 0.5,
			col = c(1:length(pp))%% 5 + 1)

	layout(1)
	par(mar = mpar)

}


#' Plots the structure of all path clusters
#'
#' Plots the structure of all path clusters
#'
#' @param ybinpaths The training paths computed by \code{\link{pathsToBinary}}.
#' @param clusters The pathway cluster model trained by \code{\link{pathCluster}} or \code{\link{pathClassifier}}.
#' @param col Colors for each path cluster.
#' @param grid A logical, whether to add a \code{\link[graphics]{grid}} to the plot
#' @param ... Extra paramaters passed to \code{plotClusterMatrix}
#'
#' @return
#' \code{plotClusterMatrix} 	plots an image of all paths the training dataset. Rows are the paths and columns
#' are the genes (features) included within each path. Paths are colored according to cluster membership.
#'
#' @author Ahmed Mohamed
#' @family Path clustering & classification methods
#' @family Plotting methods
#' @rdname plotClusters
#' @export
plotClusterMatrix <- function(ybinpaths, clusters, col=rainbow(clusters$params$M), grid=TRUE){
	if(missing(col))
		col=rainbow(clusters$params$M)

	mat <- as.matrix(ybinpaths$paths)
	if(is.null(clusters$y))
		mat.labels <- matrix(rep(clusters$labels, ncol(mat)), nrow(mat), ncol(mat))
	else
		mat.labels <- matrix(rep(clusters$component, ncol(mat)), nrow(mat), ncol(mat))

	mat[mat!=0] <- mat.labels[mat!=0]
	if(is.null(clusters$y))
		mat <- mat[order(clusters$labels),]
	else{
		mat <- mat[order(clusters$y,clusters$component),]
		mat <- cbind(mat, 'y=1'=clusters$y*(clusters$params$M+1))
		col <- c(col,"black")
	}
	mat[mat==0] <- NA

	mat <- t(mat)

	image(mat, col=col, axes=FALSE, ylab="Paths")
	# X axis
	axis(side=1, at=seq(from=0, to=1, length.out=nrow(mat)), labels= rownames(mat), las=2, cex.axis=0.6)
	# Y axis (ticks)
	axis(side=2, at=seq(from=0, to=1, length.out=ncol(mat)), tck=-0.01,labels=FALSE,las=2)

	if(!is.null(ybinpaths$y)){
		y.coor <- seq(from=0, to=1, length.out=ncol(mat)+1)
		axis(side=2, at=y.coor[which(!duplicated(ybinpaths$pidx[,2]))], labels=FALSE,, tck=1,las=2)
		axis(side=2, at=y.coor[which(!duplicated(ybinpaths$pidx[,2]))], labels=levels(ybinpaths$y) , tck=-0.03,las=2)
	}

	if(grid)
		grid(nx=nrow(mat), ny=ncol(mat))
}

#' @return
#' \code{plotClusterProbs} 	The training set posterior probabilities for each path belonging to a 3M component. \cr
#'
#' @rdname plotClusters
#' @export
#'
plotClusterProbs <- function(clusters, col=rainbow(clusters$params$M)){
	if(missing(col))
		col=rainbow(clusters$params$M)


	if(is.null(clusters$y))
		matplot(x=clusters$h[order(clusters$labels),], y=1:length(clusters$labels),
			type="l", lty=1, col=col, xlab="P(M|X)", ylab="Paths", yaxs="i")
	else
		matplot(x=clusters$path.probabilities[order(clusters$y,clusters$component),], y=1:length(clusters$component),
				type="l", lty=1, col=col, xlab="P(M|X)", ylab="Paths", yaxs="i")


}

#' @return
#' \code{plotClusters}:	combines the two plots produced by \code{plotClusterProbs} and \code{plotClusterMatrix}.
#'
#' @rdname plotClusters
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
#' 					K=20, minPathSize=8)
#'
#' 	## Convert paths to binary matrix.
#' 	ybinpaths <- pathsToBinary(ranked.p)
#' 	p.cluster <- pathCluster(ybinpaths, M=2)
#' 	plotClusters(ybinpaths, p.cluster, col=c("red", "blue") )
#'
plotClusters <- function(ybinpaths, clusters, col,...){
	if(missing(col))
		col <- rainbow(clusters$params$M)

	rect <-strwidth("M1", units="inches")
	layout(matrix(c(1,1,2),1,3,byrow=TRUE))
	par(omi=c(0,0,0, rect+0.5))

	plotClusterMatrix(ybinpaths, clusters, col=col, ...)
	plotClusterProbs(clusters, col=col)

	rect2<-legend(x=par("usr")[2]*1.05, y=nrow(clusters$h), legend=paste("M",1:clusters$params$M, sep=""),
			lty=1, lwd=2, col=col, cex=0.8, seg.len=1, xpd=NA)


	par(omi=c(0,0,0,0))
	layout(1)
}
