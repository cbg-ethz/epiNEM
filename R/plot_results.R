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
#' res <- epiNEM()
#' plot(res)
#' @return plot of the logical network
plot.epiNEM <- function(x, ...) {
    results <- x
    result <- results$origModel
    logics <- results$logics
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
