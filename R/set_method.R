## optional functions

#' @noRd
GreedyHillClimber <- function(n, experiments, data2,
                              data0, mutants, init = NULL) {
    cat("\n")
    cat(n)
    maxLlh <- -Inf
    startModel <- CreateRandomGraph(experiments)
    nonSelfEdges <- which(diag(nrow(startModel)) < 1)

    ##--- Iteratively update model and find highest scoring logics ---#
    extendedModels <- includeLogic(startModel, experiments, mutants)
    while (length(extendedModels)==0) {
        startModel <- CreateRandomGraph(experiments)
        extendedModels <- includeLogic(startModel, experiments, mutants)
    }
    extendedModels <- unlist(extendedModels, recursive=FALSE)
    result <- sapply(extendedModels, Mll, data2, data0)
    mLLscores <- unlist(result["mLL",])
    score <- max(mLLscores)
    if (is.null(init)) {
        model <- startModel
    } else {
        model <- init
    }

    while (score > maxLlh) {
        cat('.')
        maxLlh <- score
        prevModel <- model
        nextGen <- FindNeighbours(model, nonSelfEdges)
        extendedModels <- lapply(nextGen, includeLogic, experiments, mutants)
        extendedModels2 <- unlist(unlist(extendedModels, recursive=FALSE),
                                  recursive=FALSE)
        result <- sapply(extendedModels2, Mll, data2, data0)
        mll <- unlist(result["mLL",])
        index <- which.max(mll)

        score <- mll[index]
        model <- extendedModels2[[index]]$origModel
    }
    return(list(maxLlh=maxLlh, model=prevModel))
}

#' @noRd
FindNeighbours <- function(model, edges) {
    ## Return all models that differ only one edge from the current model by
    ## removing or adding a non-self edge.
    ChangeEdge <- function(edge, model) {
        if (model[edge] == 1) {
            model[edge] = 0
        } else if (model[edge] == 0) {
            model[edge] = 1
        }
        return(model)
    }
    models <- lapply(edges, ChangeEdge, model)
    models <- unique(models)
    return(models)
}

#' @noRd
EnumerateModels <- function(size, nodes=NULL) {
    ## Enumerates all possible adjacency matrices for a given size.
    ## Returns a list of
    ## all unique transitively closed adjacency matrices.
    ## Function adapted from NEM package
    if (!size %in% 2:5) {
        stop("Enumeration only feasible for networks up to 5 nodes.\n")
    }
    if (length(nodes) == 0) {
        nodes <- LETTERS[1:size]
    }

    createModels <- function(bincom, size, nodes) {
        model <- diag(size)
        model[which(model == 0)] <- bincom
        dimnames(model) <- list(nodes, nodes)
        return(list(model))
    }

    asMatrix <- function(model, size, nodes) {
        model <- matrix(model, size)
        dimnames(model) <- list(nodes, nodes)
        return(list(model))
    }

    bincom <- bincombinations(size * (size - 1))
    models <- apply(bincom, 1, createModels, size, nodes)
    models <- unique(matrix(unlist(models), ncol=size * size, byrow=TRUE))
    models <- unlist(apply(models, 1, asMatrix, size, nodes), recursive=FALSE)

    ## cat("\nGenerated", length(models),
    ## "unique models out of", nrow(bincom), "\n")
    return(models)
}
