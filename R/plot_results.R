## functions for plotting

#' extend model with node representing logic gate
#' @export
#' @param child define the child
#' @param logic define the logical gate
#' @param model normal model
#' @examples
#' model <- CreateRandomGraph(c("Ikk1", "Ikk2", "RelA"))
#' model2 <- AddLogicGates("RelA", "OR", model)
#' @return model list with additional logic gate
AddLogicGates <- function(child, logic, model) {
    ## Extends model with a node representing a logic gate for plotting
    parents <- names(which(model[,child]==1))
    redirectEdges <- function(model, parents, child) {
        ## Redirect edges via most recently added logic gate
        model[parents, child] <- 0
        model[nrow(model), child]   <- 1
        model[parents, ncol(model)] <- 1
        return(model)
    }
    model <- cbind(model, rep(0, dim(model)[1]))
    model <- rbind(model, rep(0, dim(model)[2]))
    rownames(model)[dim(model)[1]] <- logic
    colnames(model)[dim(model)[1]] <- logic
    model <- redirectEdges(model, parents, child)
    return(model)
}

#' Plots the winning pathway structure
#' @param x object of class epiNEM
#' @param ... other arguments
#' @export
#' @method plot epiNEM
#' @importFrom igraph graph.adjacency plot.igraph layout.fruchterman.reingold
#' @examples
#' data <- matrix(sample(c(0,1), 100*4, replace = TRUE), 100, 4)
#' colnames(data) <- c("A", "A.B", "B", "C")
#' rownames(data) <- paste("E", 1:100, sep = "_")
#' res <- epiNEM(data, method = "exhaustive")
#' plot(res)
#' @return plot of the logical network
plot.epiNEM <- function(x, ...) {
    results <- x
    result <- results$origModel
    logics <- results$logics
    ## logics <- paste("NOT ", gsub("masks.*", "", logics), sep = "")
    Egeneset <- results$EGeneset
    diag(result) <- 0
    nGenes    <- ncol(result)
    geneColor <- rep("lightpink", nGenes)
    geneShape <- rep("circle", nGenes)
    EGenes     <- rownames(Egeneset)
    nEGenes    <- nrow(Egeneset)
    EGeneShape <- rep("circle", nEGenes)
    EGeneColor <- rep("lightgrey", nEGenes)
    child <- colnames(result)[results$column]
    for (i in 1:length(logics)){
        result <- AddLogicGates(child[i], logics[i], result)
    }
    for (i in 1:nEGenes){
        result <- cbind(result, rep(0, dim(result)[1]))
        result <- rbind(result, rep(0, dim(result)[2]))
        rownames(result)[dim(result)[1]] <- Egeneset[i]
        colnames(result)[dim(result)[1]] <- Egeneset[i]
        ## if for unconnected node summarizing non-attached E-Genes
        if(i < nEGenes) {
            result[rownames(Egeneset)[i], ncol(result)]   <- 1
        }
    }
    labels     <- colnames(result)
    nLogics    <- length(logics)
    logicColor <- rep("lightgreen", nLogics)
    logicShape <- rep("rectangle", nLogics)
    geneSize   <- rep(30, nGenes)
    logicSize  <- unlist(lapply(logics, function(l) (nchar(l) + 5) * 3))
    EGeneSize  <- rep(30, nEGenes)
    ## lay3 as layout
    g <- graph.adjacency(result)
    plot.igraph(g, vertex.color=c(geneColor, logicColor, EGeneColor),
                vertex.size=c(geneSize, logicSize, EGeneSize),
                vertex.size2=c(rep(30,nGenes), rep(15,nLogics),
                               rep(30, nEGenes)),
                vertex.shape=c(geneShape, logicShape, EGeneShape),
                edge.arrow.size=c(rep(0.3, nGenes+nLogics+nEGenes)),
                vertex.label=labels,
                layout = layout.fruchterman.reingold)
}

#' Plots the sresults of a systematic knock-out screen
#' @param x object of class epiScreen
#' @param global plot global distribution or for each pair (FALSE)
#' @param ind index of pairs to plot
#' @param ... other arguments
#' @export
#' @method plot epiScreen
#' @examples
#' data <- matrix(sample(c(0,1), 100*9, replace = TRUE), 100, 9)
#' colnames(data) <- c("A.B", "A.C", "B.C", "A", "B", "C", "D", "E", "G")
#' rownames(data) <- paste("E", 1:100, sep = "_")
#' res <- epiScreen(data)
#' plot(res)
#' plot(res, global = FALSE, ind = 1:3)
#' @return plot(s) of an epiNEM screen analysis
plot.epiScreen <- function(x, global = TRUE, ind = NULL, ...) {

    doubles <- x$doubles

    if (is.null(ind)) { ind <- 1:length(doubles) }

    if (global) {

        distmat <- x$logic

        distmat[which(distmat == 0)] <- 7

        distmat[which(distmat %in% "AND")] <- 1
        distmat[which(distmat %in% "OR")] <- 2
        distmat[which(distmat %in% "XOR")] <- 3
        distmat[which(distmat %in% "NOEPI")] <- 6
        distmat[which(distmat %in% c("UNCON", "NOINFO", "NOINF"))] <- 7

        for (i in 1:ncol(distmat)) {

            genes <- unlist(strsplit(colnames(distmat)[i], "\\."))
            
            distmat[which(distmat[, i] %in%
                          paste(genes[1], " masks the effect of ", genes[2],
                                sep = "")), i] <- 4

            
            distmat[which(distmat[, i] %in%
                          paste(genes[2], " masks the effect of ", genes[1],
                                sep = "")), i] <- 5

        }

        distmat <- apply(distmat, c(1,2), as.numeric)

        for (i in 1:ncol(distmat)) {
            distmat[, i] <- rev(sort(distmat[, i]))
        }

        y <- distmat

        distmat <-
            distmat[, order(apply(distmat, 2, function(x) {
                return(sum(x == 1))
            }))]

        y[which(y == 5)] <- 4
        if (nrow(y) > 20) {
            rownames(distmat) <- NULL
        }
        
        labeltext <- c("", "no information\n\n\n", "no epistasis\n\n\n",
                       "masking (NOT B)\n\n\n",
                       "masking (NOT A)\n\n\n", "XOR\n\n\n",
                       "OR\n\n\n", "AND\n\n\n")

        HeatmapOP(distmat, Colv = FALSE, Rowv = FALSE,
                  main = "logic gate distribution",
                  sub = "", col = "Paired", breaks =
                                                seq(0.5,7.5, length.out = 8),
                  cexRow = 0, cexCol = 0.4, aspect = "fill",
                  colorkey = list(space = "right",
                                  labels = rev(labeltext), width = 1,
                                  at = seq(1.5,7.5, length.out = 8)),
                  xlab = "double knock-outs",
                  ylab = "modulators\n(different order for each pair)",
                  xrot = 45, bordercol = "transparent")

    } else {

        logicmat0 <- x$logic

        llmat0 <- x$ll

        for (i in 1:length(x$doubles)) {

            if (!(i %in% ind)) { next() }
            
            logicvec <- logicmat0[, i]

            llvec <- llmat0[, i]

            logicvec <- logicvec[order(llvec, decreasing = TRUE)]

            llvec <- llvec[order(llvec, decreasing = TRUE)]

            parents <- unlist(strsplit(doubles[i], "\\."))

            pchvec <- numeric(length(llvec))

            pchvec[which(logicvec %in% "AND")] <- 1
            pchvec[which(logicvec %in% "OR")] <- 2
            pchvec[which(logicvec %in% "XOR")] <- 3
            pchvec[grep(paste("^", parents[1], sep = ""), logicvec)] <- 4
            pchvec[grep(paste("^", parents[2], sep = ""), logicvec)] <- 5
            pchvec[which(logicvec %in% "NOEPI")] <- 6
            pchvec[which(logicvec %in% c("NOINFO", "NOINF"))] <- 7
            
            logicvec <- logicvec[-which(logicvec %in% "0")]
            pchvec <- pchvec[-which(pchvec == 0)]
            llvec <- llvec[-which(llvec == 0)]

            colvec <- pchvec
            
            thetop <- sum(!(logicvec %in% c("UNCON", "NOINFO", "NOINF")))
            
            if (all(is.infinite(llvec) == TRUE)) {

                llvec[1:length(llvec)] <- -1000

                margin <- 100

                donames <- 30

            } else {
                
                range <- max(llvec[1:thetop]) - min(llvec[1:thetop])

                offset <- range*0.05

                margin <- range*0.25

            }

            mark <- ""
            legendx <- length(llvec[1:thetop])
            p2max <- max(llvec[1:thetop])
            if (p2max == min(llvec[1:thetop])) {
                p2max <- p2max+margin*0.2
            }
            legendtext <- c("AND", "OR", "XOR",
                            paste(parents[1]," masks ", parents[2], sep = ""),
                            paste(parents[2], " masks ", parents[1], sep = ""),
                            "no epistasis")
            if (thetop == 0) { next() }
            plot = plot(llvec[1:thetop], pch = pchvec[1:thetop],
                        col = colvec[1:thetop],
                        ylab = "likelihood", xlab = "ranked single knockouts",
                        ylim = c(min(llvec[1:thetop]),
                                 max(llvec[1:thetop])+margin*0.2),
                        xlim = c(1, thetop+(thetop/100)),
                        main = paste(unlist(strsplit(doubles[i], "\\.")),
                                     collapse = " and "))
            text = text((1:thetop)+(thetop/100),
                        llvec[1:thetop]+offset,
                        labels = names(llvec)[1:thetop], cex = 0.6,
                        srt = 45, pos = 3,
                        offset = 0)
            mtext = mtext(mark, side = 3, line = 1, outer = FALSE,
                          cex = 4, adj = 0)
            legend = legend(legendx, p2max,
                            legend = legendtext,
                            col = 1:6, pch = 1:6, xjust = 1, yjust = 1,
                            cex = 0.7)
            
        }

    }

}

#' Plots the simulation results
#' @param x object of class epiSim
#' @param ... other arguments
#' @export
#' @method plot epiSim
#' @examples
#' res <- SimEpiNEM(runs = 1)
#' plot(res)
#' @return plot(s) of an epiNEM simulation analysis
plot.epiSim <- function(x, ...) {
    
    sens <- sim$sens
    spec <- sim$spec
    sens2 <- sim$sens2
    spec2 <- sim$spec2
    logics <- sim$logics
    time <- sim$time
    noiselvls <- sim$noiselvls

    acc <- (sens + spec)/2

    acc2 <- (sens2 + spec2)/2

    m <- rbind(c(1,1), c(2,2), c(3,4))

    layout(m)

    do <- sim$do

    do2 <- gsub("^b", "B-NEM",
                gsub("^e", "epiNEM",
                     gsub("^n", "NEM", gsub("^p", "PC algorithm",
                                            gsub("^a", "Aracne", sim$do)))))
    
    colvec <- accframe2 <- timeframe <- accframe <- logicsframe <- NULL

    colvecorg <- c("orange", "blue", "darkgreen", "brown", "darkgrey")

    for (i in 1:dim(spec2)[1]) {

        colvec <- c(colvec, rep(colvecorg[i], dim(spec2)[3]))

        timeframe[[do2[i]]] <- data.frame(time[i,,])

        colnames(timeframe[[do2[i]]]) <- noiselvls
        
        accframe2[[do2[i]]] <- data.frame(acc2[i,,])

        colnames(accframe2[[do2[i]]]) <- noiselvls

    }

    timeframe <- as.data.frame(timeframe)
    
    accframe2 <- as.data.frame(accframe2)

    colnames(timeframe) <- colnames(accframe2) <- rep(noiselvls, length(do2))
    
    boxplot(as.data.frame(timeframe), col = colvec, main = "running time",
            ylab = "seconds")

    abline(v=(1:(length(do)-1)*length(noiselvls) + 0.5), col = "black",
           lty = 6)

    axis(1, 4+(0:(length(do2)-1))*8,
         do2,
         tick = FALSE, pos = -25)

    boxplot(accframe2, col = colvec,
            main = "accuracy of the inferred edges", ylim = c(0,1))

    abline(v=(1:(length(do)-1)*length(noiselvls) + 0.5), col = "black",
           lty = 6)

    axis(1, 4+(0:(length(do2)-1))*8,
         do2,
         tick = FALSE, pos = -0.2)

    if (any(do %in% c("b", "e"))) {

        colvec2 <- unique(colvec)[which(do %in% c("b", "e"))]

        colvec2 <- rep(colvec2, each = length(noiselvls))
        
        for (i in 1:sum(do %in% c("b", "e"))) {
            
            accframe[[do2[i]]] <- data.frame(acc[i,,])
            
            colnames(accframe[[do2[i]]]) <- noiselvls
            
            logicsframe[[do2[i]]] <- data.frame(logics[i,,])
            
            colnames(logicsframe[[do2[i]]]) <- noiselvls
            
        }

        accframe <- as.data.frame(accframe)
        
        logicsframe <- as.data.frame(logicsframe)
        
        colnames(logicsframe) <- colnames(accframe) <- rep(noiselvls, 2)

        boxplot(logicsframe, col = colvec2,
                main = "accuracy of the inferred logic gate",
                ylim = c(0,1))
        
        abline(v=length(noiselvls)+0.5, col = "black", lty = 6)
        
        axis(1, 4+(0:(sum(do %in% c("b", "e"))-1))*8,
             do2[which(do %in% c("b", "e"))],
             tick = FALSE, pos = -0.2)

        boxplot(accframe, col = colvec2,
                main = "accuracy of the expected data",
                ylim = c(0,1))
        
        abline(v=length(noiselvls)+0.5, col = "black", lty = 6)
        
        axis(1, 4+(0:(sum(do %in% c("b", "e"))-1))*8,
             do2[which(do %in% c("b", "e"))],
             tick = FALSE, pos = -0.2)

    }

}


