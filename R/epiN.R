#' A binary effect matrix for a given set of knockouts.
#'
#' A publicly availiable example data set for demonstrating the use of the epiNEM package.
#' For data preprocessing steps please check vignette.
#'
#' @format A binary matrix with 2274 rows and 5 variables:
#' \describe{
#'   \item{row}{knockouts, E-Genes}
#'   \item{column}{S-Genes}
#'   ...
#' }
#' @source \url{http://www.holstegelab.nl/publications/sv/signaling_redundancy/downloads/DataS1.txt}
                                        # "epiNEMdata"


###--- MAIN SCRIPT ---###
#' Epistatic NEMs - main function
#' @param filename A binary, tab-delimited matrix. Columns: single and double knockdowns. Rows: genes showing effect or not?
#' Default: random; artificial data is generated according to 'random' specifications
#' @param method greedy or exhaustive search. Default: greedy
#' @param nIterations number of iterations. Default: 10
#' @param nModels number of Models. Default: 0
#' @param random list specifying how the data should be generated: no. of single mutants, no. of double mutants, no. of reporterGenes,
#' FP-rate, FN-rate, no. of replicates
#' @author Madeline Diekmann
#' @seealso nem
#' @export
epiNEM <- function(filename="random", method="greedy", nIterations=10, nModels=0,
                   random=list(single=4, double=1, reporters=100, FPrate=0.1,
                               FNrate=0.1, replicates=1), plotsy=TRUE) {

                                        #--- Sanity checks ---#
    methods <- c("greedy", "exhaustive")
    if (!method %in% methods)
        stop("Please enter a valid method: 'greedy' or 'exhaustive'.\n")

    if ((method == "greedy") && (nIterations < 1))
        stop("Please enter a positive number of iterations\n")

                                        #--- Process input ---#
    if (is.character(filename)){
        if (filename == "random") {
                                        # Generate artificial data
            topology    <- CreateTopology(random$single, random$double)
            topology    <- unlist(topology, recursive=FALSE)
            extTopology <- ExtendTopology(topology$model, random$reporters)
            sortedData  <- GenerateData(topology$model, extTopology, random$FPrate,
                                        random$FNrate, random$replicates)

            mutants   <- rownames(topology$model)
            singleKOs <- sort(colnames(topology$model))
            experiments <- singleKOs
            doubleKOs <- mutants[!mutants %in% singleKOs]
        } else {
                                        # Use data from file
            data        <- read.table(filename)
            sortedData  <- data[, order(unlist(sortedKOs))]
        }
    } else sortedData = filename

    mutants <- colnames(sortedData)
    singleKOs <- sort(unique(unlist(strsplit(mutants, '[.]'))))
    experiments <- singleKOs
    doubleKOs <- mutants[!mutants %in% singleKOs]

                                        #     mutants     = unique(sapply(strsplit(colnames(data), "_"), function(x){x[[1]]}))
                                        #     experiments = mutants[!grepl("\\.", mutants)]
                                        #     reporters   <- unique(rownames(data))
                                        #     dataNoWT    <- data[, (experiments[!experiments %in% grep('WT', experiments,
                                        #                                                               value=TRUE)])]
                                        #     singleKOs   <- sort(unique(unlist(strsplit(mutants, '[.]'))))
                                        #     parents = mutants[which(!mutants %in% experiments)]
                                        #     parents = unlist(strsplit(parents, ".", fixed=TRUE))
                                        #     #--- Sort mutant names (important for later use) ---#
                                        #     mutants[grep('[.]', mutants)] <- strsplit(mutants[grep('[.]', mutants)], '[.]')
                                        #     sortedKOs   <- lapply(mutants, function(d) do.call(paste, c(as.list(sort(d)),
                                        #                                                                 sep="")))
                                        #     data <- data.frame(data)
                                        #     #sortedData  <- data %>% select(order(unlist(sortedKOs)))
                                        #     sortedData  <- data[, order(unlist(sortedKOs))]
                                        #     mutants     <- sort(unlist(sortedKOs))
                                        #     doubleKOs   <- mutants[!mutants %in% singleKOs]
                                        #     #colnames(sortedData) <- mutants

                                        #data_inv <- 1-data

                                        #--- Initialize data matrices (this takes care of replicate experiments) ---#
    D1 <- sapply(mutants, function(s)
        rowSums(sortedData[, colnames(sortedData) == s, drop=FALSE]))
    D0 <- sapply(mutants, function(s) sum(colnames(sortedData) == s)) - D1

                                        #--- Optimize models ---#
    if (method == "exhaustive" & length(singleKOs) > 5) {
        cat("\nSorry, too many nodes for exhaustive search,\n")
        cat("doing greedy hill climbing instead.\n")
        method = "greedy"
    }

    else if (method == "exhaustive" & length(singleKOs) == 5) {
        cat("\nThis is going to take a while...\n")
    }

    if (method == "greedy") {
        print(paste("Progress (of ", nIterations, " iterations):", sep = ''))
                                        #iterations <- sapply(nIterations, GreedyHillClimber, experiments, dataset, data_inv)
        iterations <- sapply(1:nIterations, GreedyHillClimber, experiments, D1, D0, mutants)

        index <- which.max(iterations["maxLlh",])
        reconstruct <- includeLogic(iterations[,index]$model, experiments, mutants)
        reconstruct <- unlist(reconstruct, recursive=FALSE)
        score <- sapply(reconstruct, MLL, D1, D0)
        mll <- unlist(score["mLL",])
        index <- which.max(mll)
        posterior <- score[,index]$posterior

        results <- unlist(reconstruct[index], recursive=FALSE)
        results$score <- mll[index]
        results$EGeneset <- AttachEGenes(posterior, experiments)

        cat("\n")
    }

    else if (method == "exhaustive") {
        basicModels    <- EnumerateModels(length(singleKOs), singleKOs)
        extendedModels <- lapply(basicModels, includeLogic, experiments, unique(mutants))
        extendedModels <- unlist(unlist(extendedModels, recursive=FALSE), recursive = FALSE)
        uniqueModels   <- unique(lapply(extendedModels, function(e)
            identity((e["model"]))))
        uniqueModels <- uniqueModels[1:length(uniqueModels)-1]

        if (length(doubleKOs) == 0)
            cat("No double perturbations available --> computing NEM")
                                        #else
                                        #cat("Extended to", length(extendedModels), "models, ")
                                        #cat("of which", length(uniqueModels), "are unique")

        score <- sapply(uniqueModels, MLL, D1, D0)
        mll  <- unlist(score["mLL",])

        bestModel <- uniqueModels[[which.max(mll)]]$model
        allModels <- lapply(1:length(extendedModels), function(i)
            identity(extendedModels[[i]]$model))
        isBest    <- lapply(allModels, IsBestModel, bestModel)
        results   <- extendedModels[isBest == TRUE]
        results <- lapply(results, modifyList, list(score=mll[which.max(mll)]))
        results <- results[[1]]
        posterior <- AttachEGenes(score[,which.max(mll)]$posterior, experiments)

        results$EGeneset <- posterior
                                        #results <- lapply(results, modifyList, list(EGeneset = posterior))

                                        #     index <- which.max(mll)
                                        #     posterior <- score[,index]$posterior
                                        #
                                        #     results <- unlist(uniqueModels[index], recursive=FALSE)
                                        #
                                        #     results2 <- unlist(extendedModels[index], recursive=FALSE)
                                        #     results$score <- mll[index]
                                        #     results$EGeneset <- AttachEGenes(posterior, experiments)
                                        #
                                        #if (plotsy==TRUE){
                                        #       PlotResults(results, EGeneset)
                                        #     }
                                        #best_model <- index
                                        #best_score <- mll[index]
                                        #logic <- results$logic

    }
    return(results)#results)
}

###--- HELPER FUNCTIONS ---###

#' @noRd
#' @export
CreateTopology <- function(single, double) {
                                        # Returns the adjacency matrix of a randomly generated pathway topology
    extendedModels <- list()
    singleKOs <- LETTERS[1:single]
    experiments <- singleKOs
    doubleKOs <- lapply(1:double, function(d) GenerateDoubleKO(singleKOs))
    doubleKOs <- unlist(unique(doubleKOs))
    mutants   <- sort(c(singleKOs, doubleKOs))

    while (length(extendedModels)==0){
        startModel   <- CreateRandomGraph(singleKOs)
        extendedModels <- includeLogic(startModel, experiments, mutants)
                                        #extendedModels <- unlist(extendedModels, recursive=FALSE)
    }
    selectedModel <- sample(1:length(extendedModels), 1)
    topology      <- extendedModels[[selectedModel]]
    return(topology)
}

#' @noRd
#' @export
ExtendTopology <- function(topology, nReporters) {
                                        # Returns an extended topology in which reporters are linked to pathway genes.
                                        # The reporter genes are uniformly distributed over the pathway genes.

    reporters     <- unlist(lapply(1:nReporters, function(n) paste("reporter", n,
                                                                   sep="-")))
    linkedEffects <- sample(1:ncol(topology), nReporters, replace=TRUE)
    extTopology   <- sapply(linkedEffects, function(e) lapply(1:ncol(topology),
                                                              function(c) ifelse(e == c, 1, 0)))
    extTopology   <- matrix(unlist(extTopology), ncol=ncol(topology), byrow=TRUE)
    dimnames(extTopology) <- list(reporters, colnames(topology))
    return(extTopology)
}

#' @noRd
#' @export
GenerateData <- function(model, extTopology, FPrate, FNrate, replicates) {
                                        # Returns an artificial noisy data matrix

    MakeSomeNoise <- function(d) {
        if (d == 0)
            d <- sample(0:1, 1, prob=c(1 - FPrate, FPrate))
        else if (d == 1)
            d <- sample(0:1, 1, prob=c(FNrate, 1 - FNrate))
    }

    perfectData <- extTopology %*% t(model)
    perfectData <- perfectData[, rep(seq_len(ncol(perfectData)), replicates)]
    noisyData   <- t(apply(perfectData, 1, function(data) sapply(data, function(d)
        MakeSomeNoise(d))))

    return(noisyData)
}

#' @noRd
#' @export
IsBestModel <- function(thisModel, bestModel) {
                                        # Returns whether or not thisModel equals bestModel
    if (any(dim(thisModel) != dim(bestModel))) return(FALSE)
    else if (any(thisModel != bestModel)) return(FALSE)
    return(TRUE)
}

#' @noRd
#' @export
epi2bg <- function(t) {
    adj2dnf <- function(A) {

        dnf <- NULL

        for (i in 1:ncol(A)) {
            for (j in 1:nrow(A)) {
                if (i %in% j) { next() }
                if (A[i, j] == 1) {
                    dnf <- c(dnf, paste(colnames(A)[i], rownames(A)[j], sep = "="))
                }
                if (A[i, j] == -1) {
                    dnf <- c(dnf, paste("!", colnames(A)[i], "=", rownames(A)[j], sep = ""))
                }
            }
        }

        dnf <- unique(dnf)

        return(dnf)

    }

    if (is.matrix(t)) {

        colnames(t) <- paste("S_vs_S_", gsub("\\.", "_", colnames(t)), sep = "")

        return(t)

    } else {

        tmp <- apply(t$origModel, 2, sum)

        stim <- rownames(t$origModel)[which(tmp == min(tmp))]

        graph <- NULL

        for (i in 1:length(t$column)) {

            parents <- sort(rownames(t$origModel)[which(t$origModel[, t$column[i]] == 1)])

            child <- colnames(t$origModel)[t$column[i]]

            if (t$logics[i] %in% "OR") {

                graph <- unique(c(graph, convertGraph(adj2dnf(t$origModel))))

            }

            if (t$logics[i] %in% "AND") {

                graph <- unique(c(graph, transRed(convertGraph(adj2dnf(t$origModel)))))


                graph <- c(graph, convertGraph(graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]))

                graph <- graph[-grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]

            }

            if (t$logics[i] %in% paste(parents[2], " masks the effect of ", parents[1], sep = "")) {

                graph <- c(graph, unique(convertGraph(adj2dnf(t$origModel))))

                graph <- c(graph, gsub(parents[2], paste("!", parents[2], sep = ""), gsub("\\+\\+|^\\+", "", gsub(parents[1], "", graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]))), gsub("\\+=", "=", gsub("\\+\\+|^\\+", "", gsub(parents[2], "", graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]))))

                graph <- graph[-grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]

            }

            if (t$logics[i] %in% paste(parents[1], " masks the effect of ", parents[2], sep = "")) {

                graph <- c(graph, unique(convertGraph(adj2dnf(t$origModel))))

                graph <- c(graph, gsub(parents[1], paste("!", parents[1], sep = ""), gsub("\\+\\+|^\\+", "", gsub(parents[2], "", graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]))), gsub("\\+=", "=", gsub("\\+\\+|^\\+", "", gsub(parents[1], "", graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]))))

                graph <- graph[-grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]

            }

            if (t$logics[i] %in% "XOR") {

                graph <- unique(c(graph, convertGraph(adj2dnf(t$origModel))))

                edge <-  graph[grep(paste(paste(parents, collapse = ".*\\+.*"), child, sep = "="), graph)]

                for (j in parents) {

                    edge <- gsub(j, paste("!", j, sep = ""), edge)

                }

                graph <- unique(c(graph, edge))

            }

        }

        all <- rownames(t$origModel)

        children2 <- unique(gsub(".*=", "", graph))

        if (sum(!(all %in% children2)) > 0) {
            graph <- c(graph, paste("S", all[which(!(all %in% children2))], sep = "="))
        }

        return(transRed(unique(graph)))

    }

}

###--- END OF HELPER FUNCTIONS ---###
