#' Converts the result from pathRanker into something suitable for pathClassifier or pathCluster.
#' 
#' Converts the result from pathRanker into something suitable for pathClassifier or pathCluster.    
#' 
#' Converts a set of pathways from \code{\link{pathRanker}} 
#' into a list of binary pathway matrices. If the pathways are grouped by a response label then the 
#' \emph{pathsToBinary} returns a list labeled by response class where each element is the binary 
#' pathway matrix for each class. If the pathways are from \code{\link{pathRanker}} then a list wiht
#' a single element containing the binary pathway matrix is returned. To look up the structure of a 
#' specific binary path in the corresponding \code{ypaths} object simply use matrix index by calling
#' \code{ypaths[[ybinpaths\$pidx[i,]]]}, where \code{i} is the row in the binary paths object you 
#' wish to reference.
#' 
#' @param ypaths The result of \code{\link{pathRanker}}.
#' 
#' @return A list with the following elements.
#' \item{paths}{All paths within ypaths converted to a binary string and concatenated into the one matrix.}
#' \item{y}{The response variable.}
#' \item{pidx}{An matrix where each row specifies the location of that path within the \code{ypaths} object.}
#' 
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
#' 		weight.method = "cor", use.attr="miriam.uniprot", 
#' 		y=factor(colnames(ex_microarray)), bootstrap = FALSE)
#' 
#' 	## Get ranked paths using probabilistic shortest paths.
#'  ranked.p <- pathRanker(rgraph, method="prob.shortest.path", 
#' 					K=20, minPathSize=6)
#' 	
#' 	## Convert paths to binary matrix. 
#' 	ybinpaths <- pathsToBinary(ranked.p)
#' 	p.cluster <- pathCluster(ybinpaths, M=3)
#' 	plotClusters(ybinpaths, p.cluster, col=c("red", "green", "blue") )
#'  
pathsToBinary <- function(ypaths) {
  makeBin <- function(pathGenes,allGenes) return(as.numeric(allGenes %in% pathGenes$genes))

  if(length(ypaths$path)==0){
    stop("ypaths is a an empty list. Please rerun pathRanker with different parameters.")
  }
  # if there are response labels
  if (!is.null(names(ypaths$paths))) { 
    all.genes <- c()
    path.data <- NULL

    all.genes <- unique(unlist(lapply(unlist(ypaths$paths,FALSE),"[[","genes")))
    
    resp <- c()
    for (p in 1:length(ypaths$paths)) {
      if(length(ypaths$paths[[p]])==0) next;
      
      binpaths <- data.frame(t(sapply(ypaths$paths[[p]],makeBin,allGenes = all.genes)))
      names(binpaths) <- all.genes

      if (is.null(path.data)) path.data <- binpaths
      else path.data <- rbind(path.data,binpaths)

      resp <- c(resp,rep(names(ypaths$paths)[p],length(ypaths$paths[[p]])))
    }
    pl <- sapply(ypaths$paths,length)
    m.idx <- as.matrix(cbind(rep(1,sum(pl)),
                         rep(1:length(ypaths$paths),pl),
                         as.numeric(sequence(pl))))

    return(list(paths = path.data,y = as.factor(resp), pidx = m.idx))
  } else {
    all.genes <- unique(unlist(lapply(ypaths$paths,"[[","genes")))
    
    binpaths <- data.frame(t(sapply(ypaths$paths,makeBin,allGenes = all.genes)))

    pl <- sapply(ypaths$paths,length)
    m.idx <- as.matrix(cbind(rep(1,sum(pl)),
                         as.numeric(sequence(pl))))

    names(binpaths) <- all.genes
    return(list(paths = binpaths,pidx = m.idx))
  }
}

#' HME3M Markov pathway classifier.
#' 
#' HME3M Markov pathway classifier.    
#' 
#' Take care with selection of lambda and alpha - make sure you check that the likelihood 
#' is always increasing. 
#' 
#' @param paths The training paths computed by \code{\link{pathsToBinary}}
#' @param target.class he label of the targe class to be classified.  This label must be present 
#' as a label within the \code{paths\$y} object
#' @param M Number of components within the paths to be extracted.
#' @param alpha The PLR learning rate. (between 0 and 1).
#' @param lambda The PLR regularization parameter. (between 0 and 2)
#' @param hme3miter Maximum number of HME3M iterations.  It will stop when likelihood change is < 0.001.
#' @param plriter Maximum number of PLR iteractions. It will stop when likelihood change is < 0.001.
#' @param init Specify whether to initialize the HME3M responsibilities with the 3M model - random is recommended.
#' 
#' @return A list with the following elements.
#' A list with the following values
#' \item{h}{A dataframe with the EM responsibilities.}
#' \item{theta}{A dataframe with the Markov parameters for each component.}
#' \item{beta}{A dataframe with the PLR coefficients for each component.}
#' \item{proportions}{The probability of each HME3M component.}
#' \item{posterior.probs}{The HME3M posterior probability.}
#' \item{likelihood}{The likelihood convergence history.}
#' \item{plrplr}{The posterior predictions from each components PLR model.}
#' \item{path.probabilities}{The 3M probabilities for each path belonging to each component.}
#' \item{params}{The parameters used to build the model.}
#' \item{y}{The binary response variable used by HME3M. A 1 indicates the location of the target.class labels in \code{paths\$y}}
#' \item{perf}{The training set ROC curve AUC.}
#' \item{label}{The HME3M predicted label for each path.}
#' \item{component}{The HME3M component assignment for each path.}
#' 
#' @author Timothy Hancock and Ichigaku Takigawa
#' @references Hancock, Timothy, and Mamitsuka, Hiroshi: A Markov Classification Model for Metabolic Pathways, Workshop on Algorithms in Bioinformatics (WABI) , 2009
#' @references Hancock, Timothy, and Mamitsuka, Hiroshi: A Markov Classification Model for Metabolic Pathways, Algorithms for Molecular Biology 2010
#' 
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
#' 		weight.method = "cor", use.attr="miriam.uniprot", 
#' 		y=factor(colnames(ex_microarray)), bootstrap = FALSE)
#' 
#' 	## Get ranked paths using probabilistic shortest paths.
#'  ranked.p <- pathRanker(rgraph, method="prob.shortest.path", 
#' 					K=20, minPathSize=6)
#' 	
#' 	## Convert paths to binary matrix. 
#' 	ybinpaths <- pathsToBinary(ranked.p)
#' 	p.class <- pathClassifier(ybinpaths, target.class = "BCR/ABL", M = 3)
#' 
#' 	## Contingency table of classification performance
#' 	table(ybinpaths$y,p.class$label)
#' 
#' 	## Plotting the classifier results.
#' 	plotClassifierROC(p.class)
#' 	plotClusters(ybinpaths, p.class)
#' 
pathClassifier <- function(paths,target.class,M,alpha=1,lambda=2,hme3miter = 100,plriter = 1,init = "random") {
    if ((target.class %in% levels(paths$y)) == FALSE) stop(paste("Cannot find",target.class,"in paths$y object"))
    y <- ifelse(paths$y == target.class,1,0) 
    x <- paths$paths

    # remove constant columns to train HME3M
    varying.cols <- which(sapply(x, sd) != 0)
    tr.x <- x[varying.cols]

    if (init == "3M") {
        message("Running initial 3M model")
        # initialize with a 3M model
        pclust <- pathCluster(list(paths = tr.x),M)
        pk <- pclust$proportions
        theta <- as.matrix(pclust$theta)
        beta <- matrix(0,nrow(theta),ncol(theta))
        pmx <- as.matrix(pclust$h) 
        fits <- matrix(0.5,nrow(tr.x),ncol = M)
    
        pkm <- matrix(pk,nrow=nrow(tr.x),ncol = M,byrow = TRUE)
        hij <- pkm*pmx*fits/rowSums(pkm*pmx*fits)
    } else {
        # random initialization
        clusters <- sample(1:M,nrow(x),replace = TRUE)
        pk <- rep(1/M,M)

        theta <- as.matrix(aggregate(tr.x,by = list(clusters),mean)[-1])
        pmx <- matrix(0,nrow(tr.x),M) 
        for (k in 1:nrow(theta)) {
            res <- as.matrix(tr.x) %*% diag(as.double(theta[k,])) 
            res[tr.x == 0] <- NA
            pmx[,k] <- apply(res,1,prod,na.rm = TRUE)
        }

        beta <- matrix(0,nrow(theta),ncol(theta))
        fits <- matrix(0.5,nrow(tr.x),ncol = M)

        pkm <- matrix(pk,nrow=nrow(tr.x),ncol = M,byrow = TRUE)
        hij <- pkm*pmx*fits/rowSums(pkm*pmx*fits)
    }

	fit <- .C("hme3m_R",
		y = as.double(y),
		x = as.double(as.matrix(tr.x)),
		m = as.integer(M),
		lambda = as.double(lambda),
		alpha = as.double(alpha),
		nrow = as.integer(nrow(tr.x)),
		ncol = as.integer(ncol(tr.x)),
		hme3miter = as.integer(hme3miter),
		plriter = as.integer(plriter),
		H = as.double(hij),
		PATHPROBS = as.double(pmx),
		PLRPRE = as.double(fits),
		THETA = as.double(theta),
		BETA = as.double(beta),
		PROPORTIONS = as.double(pk),
		HMEPRE = double(nrow(tr.x)),
		LIKELIHOOD = double(hme3miter))

    theta <- matrix(NA,nrow = M,ncol = ncol(x))
    theta[,c(1,ncol(theta))] <- 1
	t.complete <- matrix(as.double(fit$THETA),nrow = M,ncol = ncol(tr.x) ,byrow = TRUE)
    theta[,varying.cols] <- t.complete # add back the constant columns
    theta <- data.frame(theta)

    beta <- matrix(NA,nrow = M,ncol = ncol(x))
    b.complete <- matrix(as.double(fit$BETA),nrow = M,ncol = ncol(tr.x) ,byrow = TRUE)
    beta[,varying.cols] <- b.complete # add back the constant columns
    beta <- data.frame(beta)
	names(beta) <- names(theta) <- names(x)

	hij <- matrix(fit$H,nrow(x),M) 
	fits <- matrix(fit$PLRPRE,nrow(x)) 
	#perf <- compROC(y,fit$HMEPRE)$auc

    # cluster labels
    zm <- apply(hij,1,max)
    cl <- apply(hij == zm,1,which)
    if (is.list(cl)) {
        message("Multiple cluster labels can be assigned to some paths: some clusters might not be valid")
        clusters <-  paste("M",sapply(cl,"[[",1),sep = "")
    } else clusters <- cl

    pmx <- matrix(fit$PATHPROBS,nrow(x),M) 
    pmx <- pmx/rowSums(pmx)
	output <- list(h = hij,
		theta = theta,
		beta = beta,
		proportions = fit$PROPORTIONS,
		posterior.probs = fit$HMEPRE,
		likelihood = fit$LIKELIHOOD[1:fit$hme3miter],
		plr.probabilities = fits,
        path.probabilities = pmx,
		params = list(alpha = alpha,lambda = lambda,M = M),
		y = y,
		#perf = perf,
        labels = ifelse(fit$HMEPRE > 0.5,1,0),
        component = clusters)

	return(output)
}

#' Predicts new paths given a pathClassifier model.
#' 
#' Predicts new paths given a pathClassifier model.    
#' 
#' @param mix The result from \code{\link{pathClassifier}}.
#' @param newdata A data.frame containing the new paths to be classified.
#' 
#' @return A list with the following elements.
#' \item{h}{The posterior probabilities for each HME3M component.}
#' \item{posterior.probs}{The posterior probabilities for HME3M model to classify the response.}
#' \item{label}{A vector indicating the HME3M cluster membership.}
#' \item{component}{The HME3M component membership for each pathway.}
#' \item{path.probabilities}{The 3M path probabilities.}
#' \item{plr.probabilities}{The PLR predictions for each component.} 
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
#' 		weight.method = "cor", use.attr="miriam.uniprot", 
#' 		y=factor(colnames(ex_microarray)), bootstrap = FALSE)
#' 
#' 	## Get ranked paths using probabilistic shortest paths.
#'  ranked.p <- pathRanker(rgraph, method="prob.shortest.path", 
#' 					K=20, minPathSize=6)
#' 	
#' 	## Convert paths to binary matrix. 
#' 	ybinpaths <- pathsToBinary(ranked.p)
#' 	p.class <- pathClassifier(ybinpaths, target.class = "BCR/ABL", M = 3)
#' 
#' 	## Just an example of how to predict cluster membership
#'  pclass.pred <- predictPathCluster(p.class, ybinpaths$paths)
#' 
predictPathClassifier <- function(mix,newdata) {
    pexp <- matrix(0,nrow=nrow(newdata),ncol = nrow(mix$theta))
    pmx <- matrix(0,nrow=nrow(newdata),ncol = nrow(mix$theta))
    pk <- mix$proportions

    tt <- mix$theta
    tt[is.na(mix$theta)] <- 1
    tb <- mix$beta
    tb[is.na(mix$beta)] <- 0
    for (k in 1:nrow(mix$theta)) {
        pred <- as.matrix(newdata) %*% t(tb[k,])
        pexp[,k] <- 1/(1+exp(-pred))
        res <- sweep(as.matrix(newdata),2,t(tt[k,]),FUN = "*")
		res[newdata == 0] <- NA
		pmx[,k] <- apply(res,1,prod,na.rm = TRUE)
    }
    h <- sweep(pmx*pexp,2,t(pk),FUN = "*")
    h <- h/rowSums(h)

    # cluster labels
    zm <- apply(h,1,max)
    cl <- apply(h == zm,1,which)
    if (is.list(cl)) {
        message("Multiple cluster labels for some paths: some clusters might not be valid")
        clusters <-  paste("M",sapply(cl,"[[",1),sep = "")
    } else clusters <- cl

    post <- rowSums(pmx*pexp/rowSums(pmx))

    return(list(h = h,
                posterior.probs = post,
                label = ifelse(post > 0.5,1,0),
                component = clusters,
                path.probabilities = pmx/rowSums(pmx),
                plr.probabilities = pexp))
}

#' Diagnostic plots for pathClassifier.
#' 
#' Diagnostic plots for \code{\link{pathClassifier}}.    
#' 
#' @param mix The result from \code{\link{pathClassifier}}.
#' 
#' @return Diagnostic plots of the result from pathClassifier.
#' item{Top}{ROC curves for the posterior probabilities (\code{mix\$posterior.probs}) 
#' and for each HME3M component (\code{mix\$h}).  This gives information about what response 
#' label each relates to. A ROC curve with an \code{AUC < 0.5} relates to \code{y = 0}. 
#' Conversely ROC curves with \code{AUC > 0.5} relate to \code{y = 1}. }
#' item{Bottom}{The likelihood convergence history for the HME3M model.  If the parameters 
#' \code{alpha} or \code{lambda} are set too large then the likelihood may decrease.}
#' 
#' @author Timothy Hancock and Ichigaku Takigawa
#' @family Path clustering & classification methods
#' @family Plotting methods
#' @export
#' 
plotClassifierROC <- function(mix) {   
    palette("default") 

    layout(c(1,2), widths=c(1,1), heights=c(0.7,0.3))
    plotPathROC(mix)

    plot(na.omit(mix$likelihood),type="l",col=2,
		xlab="EM Iteration",ylab="Conditional\nLog-Likelihood",main="Likelihood Convergence",cex = 0.7)
}

#' Plots the structure of specified path found by pathClassifier.
#' 
#' Plots the structure of specified path found by pathClassifier.
#' 
#' @param ybinpaths The training paths computed by \code{\link{pathsToBinary}}
#' @param obj The pathClassifier \code{\link{pathClassifier}}.
#' @param m The path component to view.
#' @param tol A tolerance for 3M parameter \code{theta} which is the probability for each edge within each cluster.
#' If the tolerance is set all edges with a \code{theta} below that tolerance will be removed from the plot.
#' 
#' @return Produces a plot of the paths with the path probabilities and prediction probabilities and ROC curve overlayed.
#' \item{Center Plot}{An image of all paths the training dataset.  Rows are the paths and columns are the genes (vertices) 
#' included within each pathway.  A colour within image indicates if a particular gene (vertex) is included within a specific path.  
#' Colours flag whether a path belongs to the current HME3M component (P > 0.5).}
#' \item{Center Right}{The training set posterior probabilities for each path belonging to the current 3M component.}
#' \item{Center Top}{The ROC curve for this HME3M component.}
#' \item{Top Bar Plots}{\code{Theta}: The 3M component probabilities - indicates the importance of each edge is to a path.
#' \code{Beta}: The PLR coefficient - the magnitude indicates the importance of the edge to the classify the response.}
#' 
#' @author Timothy Hancock and Ichigaku Takigawa
#' 
#' @family Path clustering & classification methods
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
#' 	## Convert paths to binary matrix. 
#' 	ybinpaths <- pathsToBinary(ranked.p)
#' 	p.class <- pathClassifier(ybinpaths, target.class = "BCR/ABL", M = 3)
#' 
#' 	## Plotting the classifier results.
#' 	plotClassifierROC(p.class)
#' 	plotClusters(ybinpaths, p.class)
#' 
plotPathClassifier <- function(ybinpaths,obj,m,tol = NULL) {
    pp <- obj$theta[m,]
    fidx <- 1:length(pp)
    if (!is.null(tol)) fidx <- which((pp >= tol) | is.na(pp))
    pp <- pp[fidx]

    g <- names(pp)
	x <- ybinpaths$paths[fidx]
		
#    gc <- strsplit(g,":")[-c(1,length(g))]
#    frt <- c("",paste(sapply(gc,"[[",2),sapply(gc,"[[",4),sapply(gc,"[[",3),sep = "-"),"")
#    gn <- c("s",sapply(gc,"[[",1),"t")
#    pname <- sapply(gc,"[[",5)
#    pcol <- c(0,as.numeric(as.factor(pname)),0)
#    
    y <- ybinpaths$y
    h <- obj$h[,m]    
    bp <- obj$beta[m,][fidx]

    palette("default") 
    mpar <- par()$mar

    sy <- sort(h,index.return = TRUE)
 	ypred <- as.numeric(h > 0.5)

    n <- layout(matrix(c(1,2,4,5,5,3),3,2),
		heights = c(.15,0.15,.7,.33,.7),
		widths =  c(.7,.3,.7,.7,.3), TRUE)

    par(mar = c(0,8,2,0),pty = "m")
    xpos <- barplot(height = as.numeric(bp),space = 0,xlim = c(0,ncol(x)),xaxs = "i",
                ylab = "beta", col = c(1:length(pp))%% 5 + 1)
    abline(h = 0,lwd = 2)

	par(mar = c(0,8,1,0),pty = "m")
	xpos <- barplot(height = as.numeric(pp),space = 0,xlim = c(0,ncol(x)),xaxs = "i",ylim = c(0,1),
	ylab = "theta", col = c(1:length(pp))%% 5 + 1)
    abline(h = 0,lwd = 2)

	par(mar = c(8,0,0,4))	
	plot(NA,xlim = c(0,1),ylim = c(1,nrow(x)),
		xlab = paste("Probability P(M=",m,"|x,y)",sep =""),ylab = "",yaxs = "i",type = "l",axes = FALSE)
    points(x = h[sy$ix],y = 1:nrow(x),type = "l",col = 1)
	abline(v = 0.5,col = 2,lwd = 2)
	abline(v = 0,col = 1,lwd = 2)
	axis(side = 1,at = seq(0,1,1/10),las = 1)
	abline(h = nrow(x)+1)

	ystats <- as.numeric(summary(as.factor(ypred)))
	ylab <- c("","","","Not a\nMember","Member")
	ytick <- c(1,nrow(x),ystats[1],ystats[1]/2,ystats[2]/2 + ystats[1])
	
	par(mar = c(8,8,0,0))
	px <- x[sy$ix,]
    px[px == 0] <- NA
    ry <- matrix(rep(sy$x < 0.5,ncol(x)),nrow(x),ncol(x))
	px[ry] <- px[ry] + 1
	image(x = xpos,y = 1:nrow(x),as.matrix(t(px)),axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i",col = cm.colors(2))
	axis(side=2, at=ytick, labels=ylab, las=2)
    mtext(g, side=1, at=xpos, cex=0.6, las=2, line=0.5,
        col = c(1:length(pp))%% 5 + 1)
#    mtext(frt,side = 1,at = xpos,cex = 0.6,las = 2,line = 5,
#        col = c(1:length(pp))%% 5 + 1)
#    pd <- cbind(xpos[-c(1,length(xpos))],pcol[-c(1,length(pcol))])
#    pht <- aggregate(pd,list(pname),median)
#    mtext(pht[,1],side = 1,at = pht[,2],cex = 0.6,las = 2,line = 15,col = as.numeric(as.factor(pht[,3])) %% 5 + 1)
    mtext("All Paths\nSorted By Component Membership",line = 5,side = 2)

	par(mar = c(1,1,6,6))
    rocs <- compROC(obj$y,obj$h[,m]) 
    plot(NA,xlim = c(0,1),ylim = c(0,1),axes = FALSE)
    abline(0,1)
    lines(x = rocs$fnr,y = rocs$tpr,type = "l",col = 2,lwd = 2)
    axis(side = 3,at = seq(0,1,1/5),labels = seq(0,1,1/5),cex= 0.5)
    mtext("FNR",3,line = 2,cex = 0.7)
    axis(side = 4,at = seq(0,1,1/5),labels = seq(0,1,1/5),cex = 0.5)
    mtext("TPR",4,line= 2,cex = 0.7)
    mtext(paste("ROC for Path",m,"\nAUC = ",round(rocs$auc,3)),line = 3,cex = 0.7)

    mtext(paste("Structure of Path",m),side = 3,outer = TRUE,line = -2,cex = 2)

    layout(1)
    par(mar = mpar)
}

compROC <- function(y,yprob) {
	tpr <- c(0)
	fnr <- c(0)
	for (i in seq(0,1,1/1000)) {
		tpr <- c(tpr, sum(y == 1 & yprob >= i)  / sum(y == 1) ) 
		fnr <- c(fnr, sum(y == 0 & yprob >= i) /  sum(y == 0) )
	}
	trap.rule <- function(x, y) {
		idx <- 2:length(x)
		dx <- x[idx] - x[idx-1] 
		dy <- y[idx] + y[idx-1] 
		auc <- dx%*%dy / 2
		return( auc )
	}
	# The probabilites can get sparse and actual 0's and 1's occur!  
	# This causes problems with the AUC.
	if (sum(tpr == 0) > 0) fnr[tpr == 0] <- max(fnr[tpr == 0])
	if (sum(fnr == 0) > 0)tpr[fnr == 0] <- max(tpr[fnr == 0])
	
	fnr <- sort(fnr,index.return = TRUE)
	tpr <- tpr[fnr$ix]
	auc <- trap.rule(x = fnr$x,y = tpr)
	
	return(list(tpr = tpr,fnr = fnr$x,auc = auc))
}

plotPathROC <- function(pfit) {
	palette("default")
	
	plot(NA,xlim = c(0,2),ylim = c(0,1),axes = FALSE,
			ylab = "True Positive Rate", xlab = "False Positive Rate",
			main ="ROC Curve For Each HME3M Component")
	
	y <- pfit$y
	
	auc <- c()
	mpre <- pfit$h
	yprob <- cbind(mpre,pfit$posterior.probs)
	auc <- c()
	for (j in 1:ncol(data.frame(yprob))) {  
		zroc <- compROC(y,yprob[,j])
		lines(x = zroc$fnr, y = zroc$tpr,col = j,lwd = 2)
		auc <- c(auc,zroc$auc)
	}
	axis(1,seq(0,1,0.1),seq(0,1,0.1),pos = 0)
	axis(2,seq(0,1,0.1),seq(0,1,0.1),pos = 0)
	lines(x = c(0,1),y = c(1,1))
	lines(x = c(1,1),y = c(0,1))
	lines(x = c(0,1),y = c(0,1),lwd = 1,lty = "dashed")
	mtext("False Positive Rate",side = 1,line = 2,adj = .2)
	lnames <- c(paste("M = ",1:ncol(pfit$h)),"Complete")
	lnames <- paste(lnames,"(AUC =",round(auc,3),")",sep = "")
	legend(x = 1.05,y = .95,lnames,col = 1:ncol(yprob),lty = rep(1,ncol(yprob)),lwd = 2,bg = "white",cex = 0.8)
}

