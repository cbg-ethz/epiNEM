#' AUC permutation test
#'
#' computes the area under the rank enrichment score curve and does 
#' a permutation test to compute the p-value
#' @param x numeric vector of ranks
#' @param y numeric vector of the superset of x
#' @param alternative character for test type: 'less','greater','two.sided'
#' @param iter integer number of iterations
#' @return p-value
#' @author Martin Pirkl
#' @export
#' @examples
#' x <- 1:10
#' y <- 1:100
#' perm.rank.test(x,y,alternative='less')
#' perm.rank.test(x,y,alternative='greater')
perm.rank.test <- function(x,y=NULL,
                 alternative=c("two.sided", "less", "greater"),iter=1000,...) {
    xscore0 <- c(1,sort(x),length(y))
    scorel <- sort(x)/length(y)-
        (1:length(x))/length(x)
    score0 <- c(0,scorel,0)
    aucr <-
        sum(((score0[-length(score0)]+
                  score0[-1])/2)*(xscore0[-1]-
                                      xscore0[-length(xscore0)]))
    scorer <- (1:length(x))/length(x)-
        sort(x)/length(y)
    score0 <- c(0,scorer,0)
    aucl <- sum(((score0[-length(score0)]+
                      score0[-1])/2)*(xscore0[-1]-
                                          xscore0[-length(xscore0)]))
    count <- c(0,0)
    for (i in seq_len(iter)) {
        x.p <- sample(y,length(x))
        xscore0 <- c(1,sort(x.p),length(y))
        scorel <- sort(x.p)/length(y)-
            (1:length(x.p))/length(x.p)
        score0 <- c(0,scorel,0)
        aucr.p <-
            sum(((score0[-length(score0)]+
                      score0[-1])/2)*(xscore0[-1]-
                                          xscore0[-length(xscore0)]))
        scorer <- (1:length(x.p))/length(x.p)-
            sort(x.p)/length(y)
        score0 <- c(0,scorer,0)
        aucl.p <- sum(((score0[-length(score0)]+
                            score0[-1])/2)*(xscore0[-1]-
                                                xscore0[-length(xscore0)]))
        if (aucl<=aucl.p) {
            count[1] <- count[1]+1
        }
        if (aucr<=aucr.p) {
            count[2] <- count[2]+1
        }
    }
    p.values <- count/iter
    if (alternative[1]=='less') {
        return(list(p.value=p.values[1]))
    } else if (alternative[1]=='greater') {
        return(list(p.value=p.values[2]))
    } else if (alternative[1]=='two.sided') {
        p.tmp <- 2*min(p.values)
        return(list(p.value=p.tmp))
    }
}
#' Rank enrichment
#'
#' Infers a signalling pathway from peerturbation experiments.
#' @param data m times l matrix with m observed genes and l
#' variables with numeric values to rank the genes
#' @param list list of of vectors of genes
#' @param list2 optional list with same length as list
#' @param n length of the gradient (maximum: m)
#' @param main character string for main header; if NULL uses
#' the column names of data by default
#' @param col1 color of the gradient
#' @param col2 color of the first list
#' @param col3 color of the second list2
#' @param blim numeric vector of length two with the lower
#' and upper bounds for the gradient
#' @param p numeric adjustment (length four) of the left side of the
#' gradient (low means more to the left, high more to the right)
#' the right side of the enrichment lines and the top positions
#' of the additional matrices in case of vis='matrices'
#' @param lwd line width of the enrichment lines
#' @param test test function for the enrichment p-value; must
#' have input argument and output values same as perm.rank.test;
#' e.g., wilcox.test or ks.test (here 'less' and 'greater' are switched!)
#' @param vis method for visualisation: 'matrix' uses one
#' matrix heatmap for; 'matrices' uses several matrices
#' (experimental), 'colside'
#' uses the colSideColors argument for the ticks of genes in
#' list/list2 (can use a lot of memory; experimental)
#' @param verbose if TRUE gives prints additional output
#' @param ... additional arguments for epiNEM::HeatmapOP
#' @return transitively closed matrix or graphNEL
#' @author Martin Pirkl
#' @export
#' @import lattice
#' @importFrom latex2exp TeX
#' @examples
#' data <- matrix(rnorm(100*2),100,2)
#' rownames(data) <- 1:100
#' colnames(data) <- LETTERS[1:2]
#' list <- list(first = as.character(sample(1:100, 10)), second = as.character(sample(1:100, 20)))
#' rank.enrichment(data,list)
rank.enrichment <- function(data,list,list2=NULL,n=1000,main=NULL,col1="RdBu",
                            col2=rgb(1,0,0,0.75),col3=rgb(0,0,1,0.75),
                            blim=NULL,p=NULL,lwd=3,
                            test=wilcox.test,vis="matrix",
                            verbose=FALSE,...) {
    if (is.null(p)) {
        p <- c(0,1,0.466,0.333)
    }
    pvals <- matrix(NA,ncol(data),length(list))
    ylim <- c(-0,1)
    if (is.null(list2)) {
        pvals <- matrix(NA,length(list)*ncol(data),3)
        colnames(pvals) <- c("less","greater","two.sided")
    } else {
        pvals <- matrix(NA,length(list)*ncol(data),2*3)
        colnames(pvals) <- rep(c("less","greater","two.sided"),2)
    }
    rownames(pvals) <- rep(colnames(data),each=length(list))
    for (i in 1:ncol(data)) {
        data.tmp <- sort(data[,i],decreasing=FALSE)
        if (is.null(blim)) {
            blim2 <- c(-max(abs(data.tmp)),max(abs(data.tmp)))
        } else {
            blim2 <- blim
        }
        for (j in 1:length(list)) {
            tmp <- match(list[[j]],names(data.tmp))
            tmp <- tmp[which(is.na(tmp) == FALSE)]
            pvals[j+length(list)*(i-1),1] <-
                test(tmp,rank(data.tmp),alternative="less")$p.value
            pvals[j+length(list)*(i-1),2] <-
                test(tmp,rank(data.tmp),alternative="greater")$p.value
            pvals[j+length(list)*(i-1),3] <-
                2*min(pvals[j+length(list)*(i-1),1:2])
            mat.tmp <- rbind(data.tmp,NA)
            mat.tmp[2,tmp] <- blim2[1]
            if (vis == 'colside') {
                csc1 <- mat.tmp[2,]
                csc1[tmp] <- col2
                csc1[is.na(csc1)] <- "transparent"
            }
            xscore0 <- c(1,sort(tmp),length(data.tmp))
            ## left:
            scorel <- sort(tmp)/length(data.tmp)-
                (1:length(tmp))/length(tmp)
            score0 <- c(0,scorel,0)
            aucl <-
                sum(((score0[-length(score0)]+
                          score0[-1])/2)*(xscore0[-1]-
                                              xscore0[-length(xscore0)]))
            ## right:
            scorer <- (1:length(tmp))/length(tmp)-
                sort(tmp)/length(data.tmp)
            score0 <- c(0,scorer,0)
            aucr <- sum(((score0[-length(score0)]+
                              score0[-1])/2)*(xscore0[-1]-
                                                  xscore0[-length(xscore0)]))
            if (aucl >= aucr) {
                score <- scorel
            } else {
                score <- scorer
            }
            ## two-sided?:
            score <- abs(score)
            ## all from here:
            score0 <- c(0,score,0)
            xscore0 <- c(1,sort(tmp),length(data.tmp))
            ## auc:
            auc <- sum(((score0[-length(score0)]+
                             score0[-1])/2)*(xscore0[-1]-
                                                 xscore0[-length(xscore0)]))
            if (!is.null(list2)) {
                tmp2 <- tmp
                tmp <- match(list2[[j]],names(data.tmp))
                tmp <- tmp[which(is.na(tmp) == FALSE)]
                pvals[j+length(list)*(i-1),4] <-
                    test(tmp,rank(data.tmp),alternative="less")$p.value
                pvals[j+length(list)*(i-1),5] <-
                    test(tmp,rank(data.tmp),alternative="greater")$p.value
                pvals[j+length(list)*(i-1),6] <-
                    test(tmp,rank(data.tmp))$p.value
                mat.tmp <- rbind(mat.tmp,NA)
                mat.tmp[3,tmp] <- blim2[2]
                if (vis == 'colside') {
                    csc2 <- mat.tmp[3,]
                    csc2[tmp] <- col3
                    csc2[is.na(csc2)] <- "transparent"
                }
                score <- sort(tmp)/length(data.tmp)-
                    (1:length(tmp))/length(tmp)
                score <- c(0,abs(score),0)
                xscore <- c(1,sort(tmp),length(data.tmp))
                tmp <- c(tmp,tmp2)
            } else {
                score <- score0
                xscore <- xscore0
                col3 <- "transparent"
            }
            select <- sort(c(tmp,(1:n)*length(data.tmp)/n))
            mat.tmp <- mat.tmp[,select]
            csc0 <- NULL
            if (vis == 'colside') {
                csc1 <- csc1[select]
                csc2 <- csc2[select]
                csc0 <- rep("transparent",length(csc1))
                csc0[which(csc1 != "transparent")] <- col2
                csc0[which(csc2 != "transparent")] <- col3
                csc0[which(csc1 != "transparent" & csc2 != "transparent")] <-
                    "black"
            }
            colnames(mat.tmp) <- NULL
            rownames(mat.tmp) <- NULL
            if (is.null(main)) {
                main2 <- colnames(data)[i]
            } else {
                main2 <- main
            }
            if (is.null(list2)) {
                p.text <- c(paste0(gsub("_"," ",names(list)[j]), " (p-value: ",
                                   format(pvals[j+length(list)*(i-1),2],digits=3),
                                   " (greater), ",
                                   format(pvals[j+length(list)*(i-1),1],
                                          digits=3),
                                   " (less))"))
            } else {
                p.text <- c(paste0(gsub("_"," ",names(list)[j]), " (p-value: ",
                                   format(pvals[j+length(list)*(i-1),2],
                                          digits=3),
                                   " (greater), ",
                                   format(pvals[j+length(list)*(i-1),1],
                                          digits=3),
                                   " (less))"),
                            paste0(gsub("_"," ",names(list2)[j]),
                                   " (p-value: ",
                                   format(pvals[j+length(list)*(i-1),5],
                                          digits=3), " (greater), ",
                                   format(pvals[j+length(list)*(i-1),4],
                                          digits=3)," (less))"))
            }
    a <- lattice::xyplot(score0~xscore0,
                         strip=FALSE,type="l",ylim=ylim,
                         ylab="enrichment score",xlab="",col=col2,
                         key=list(corner=c(0.5,1),
                                  lines=list(col=c(col2,col3),
                                             lty=c(1,1),lwd=10),
    text=list(p.text)),
    lwd=lwd,par.settings=list(axis.line=list(lwd=0)),
    scales=list(x=list(draw=FALSE)),
    panel = function(...) {
        panel.xyplot(...)
        panel.abline(h=seq(-0.8,0.8,length.out=9),lty=3,col=rgb(0,0,0,0.75))
        panel.xyplot(x=xscore,y=score,type="l",lwd=lwd,col=col3)
    },
    main =  main2
            )
            if (!is.null(list2)) {
                sub <- ""
            } else {
                sub <- NULL
            }
    more <- FALSE
            if (vis == 'colside') {
                mat.tmp <- mat.tmp[1,]
            } else if (vis == 'matrices') {
                mat.tmp1 <- mat.tmp[2:3,]
                mat.tmp1[which(is.na(mat.tmp1)==TRUE)] <- blim2[2]
                mat.tmp1[2,] <- NA
                if (!is.null(list2)) {
                    mat.tmp2 <- mat.tmp[3,,drop=FALSE]
                    mat.tmp2[which(is.na(mat.tmp2)==TRUE)] <- blim2[1]
                }
                more <- TRUE
                mat.tmp[2:3,] <- NA
            }
            b <- HeatmapOP(mat.tmp,Colv=0,Rowv=0,
                                   bordercol="transparent",col=col1,
                                   borderwidth=0,colNA="transparent",
                                   breaks=seq(blim2[1],blim2[2],
                                              length.out=100),
                                   sub = sub, colSideColors = csc0, ...)
            print(a,position=c(0,0.5,p[2],1),more=TRUE)
            print(b,position=c(p[1],0,1,0.6),more=more)
            if (vis == 'matrices') {
                c1 <- HeatmapOP(mat.tmp1,Colv=0,Rowv=0,
                                bordercol="transparent",
                                col=c(col2,"white"),
                                borderwidth=0,colNA="transparent",
                                breaks=seq(blim2[1],blim2[2],
                                           length.out=100),
                                sub = sub, ...)
                more2 <- FALSE
                if (!is.null(list2)) {
                    c2 <- HeatmapOP(mat.tmp2,Colv=0,Rowv=0,
                                    bordercol="transparent",
                                    col=c("white",col3),
                                    borderwidth=0,colNA="transparent",
                                    breaks=seq(blim2[1],blim2[2],
                                               length.out=100),
                                    sub = sub, ...)
                    more2 <- TRUE
                }
                print(c1,position=c(p[1],0,1,p[3]),more=more2)
                if (!is.null(list2)) {
                    print(c2,position=c(p[1],0,1,p[4]))
                }
            }
            if (verbose) {
                print("Proliferative:")
                print(pvals[j+length(list)*(i-1),1:3])
                if (!is.null(list2)) {
                    print("Invasive:")
                    print(pvals[j+length(list)*(i-1),4:6])
                }
            }
        }
    }
    fdrs <- pvals
    fdrs[1:length(fdrs)] <- p.adjust(as.vector(pvals),method="BH")
    return(list(pvals=pvals,fdrs=fdrs))
}
#' Add logic.
#'
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

#' Plot pathway.
#'
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

#' Plot screen.
#'
#' Plots the sresults of a systematic knock-out screen
#' @param x object of class epiScreen
#' @param global plot global distribution or for each pair (FALSE)
#' @param ind index of pairs to plot
#' @param colorkey if TRUE prints colorkey
#' @param cexGene size of modulator annotation
#' @param off relative distance from the gene names to the respective
#' likelihoods
#' @param cexLegend font size of the legend
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
plot.epiScreen <- function(x, global = TRUE, ind = NULL, colorkey = TRUE,
                           cexGene = 1, off = 0.05, cexLegend = 1, ...) {

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
        ## if (nrow(y) > 20) {
            rownames(distmat) <- NULL
        ## }

        labeltext <- c("", "no information\n\n\n", "no epistasis\n\n\n",
                       "masking (NOT B)\n\n\n",
                       "masking (NOT A)\n\n\n", "XOR\n\n\n",
                       "OR\n\n\n", "AND\n\n\n")

        if (colorkey) {
            colorkey <- list(space = "right",
                                  labels = rev(labeltext), width = 1,
                             at = seq(1.5,7.5, length.out = 8))
        } else {
            colorkey <- NULL
        }

        HeatmapOP(distmat, Colv = FALSE, Rowv = FALSE,
                  main = "logic gate distribution",
                  sub = "", col = "Paired", breaks =
                                                seq(0.5,7.5, length.out = 8),
                  colorkey = colorkey,
                  xlab = "double knock-outs",
                  ylab = "modulators\n(different order for each pair)", ...)

    } else {

        palette(c("#4444cc", "#77aa77", "#009933",
                  "#ff0000", "#dd8811", "#aa44bb", "#999900"))

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

                if (range == 0) {

                    range <- 10

                    margin <- range*0.25

                    offset <- range*off

                    ylim <- c(llvec[1], llvec[1]+10)

                } else {

                    margin <- range*0.25

                    offset <- range*off

                    ylim <- c(min(llvec[1:thetop]),
                              max(llvec[1:thetop])+margin*0.2+offset+margin*3/5)

                }

                offset <- range*off

            }

            mark <- ""
            legendx <- length(llvec[1:thetop])
            p2max <- max(llvec[1:thetop])
            if (p2max == min(llvec[1:thetop])) {
                p2max <- p2max+margin*0.2
            }
            legendtext <- c("AND", "OR", "XOR",
                            paste(tolower(parents[1])," masks ",
                                  tolower(parents[2]), sep = ""),
                            paste(tolower(parents[2]), " masks ",
                                  tolower(parents[1]), sep = ""),
                            "no epistasis")
            if (thetop == 0) { next() }
            plot = plot(llvec[1:thetop], pch = pchvec[1:thetop],
                        col = colvec[1:thetop],
                        ylab = "likelihood", xlab = "ranked single knockouts",
                        ylim = ylim,
                        xlim = c(1, thetop+(thetop/100)),
                        main = tolower(paste(unlist(strsplit(doubles[i],
                                                             "\\.")),
                                             collapse = " and ")))
            text = text((1:thetop)+(thetop/100),
                        llvec[1:thetop]+offset,
                        labels = tolower(names(llvec)[1:thetop]),
                        srt = 45, pos = 3,
                        offset = 0, cex = cexGene, ...)
            mtext = mtext(mark, side = 3, line = 1, outer = FALSE,
                          cex = 4, adj = 0, ...)
            legend = legend("topright",
                ## legendx, p2max,
                legend = legendtext,
                col = 1:6, pch = 1:6, xjust = 1, yjust = 1,
                cex = cexLegend)

        }

    }

}

#' Plot simulations.
#'
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


