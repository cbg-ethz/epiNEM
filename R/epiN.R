#' Example data: simulation results
#' 
#' Contains simulation results
#' 
#' @docType data
#' @examples 
#' data(sim)
#' @name sim
NA

#' Example data: epiNEM results for the Sameith et al., 2015 knock-out screen
#' 
#' @docType data
#' @examples 
#' data(sameith)
#' @name sameith
NA

#' Example data: epiNEM results for the Wageningen et al., 2010 knock-out screen
#' 
#' @docType data
#' @examples 
#' data(wageningen)
#' @name wageningen
NA

#' sig. of string interaction scores for Sameith et al., 2015 data
#' 
#' @docType data
#' @examples 
#' data(sameith_string)
#' @name sameith_string
NA

#' sig. of string interaction scores for van Wageningen et al., 2010 data
#' 
#' @docType data
#' @examples 
#' data(wageningen_string)
#' @name wageningen_string
NA

###--- MAIN SCRIPT ---###
#' Epistatic NEMs - main function
#' @param filename A binary, tab-delimited matrix.
#' Columns: single and double knockdowns. Rows: genes showing effect or not?
#' Default: random; artificial data is generated to 'random' specifications
#' @param method greedy or exhaustive search. Default: greedy
#' @param nIterations number of iterations. Default: 10
#' @param nModels number of Models. Default: 0
#' @param random list specifying how the data should be generated:
#' no. of single mutants, no. of double mutants, no. of reporterGenes,
#' FP-rate, FN-rate, no. of replicates
#' @param plotsy atm not used
#' @param ltype likelihood either "marginal" or "maximum"
#' @param para false positive and false negative rates
#' @author Madeline Diekmann
#' @seealso nem
#' @export
#' @import
#' igraph
#' e1071
#' nem
#' utils
#' @return optimized logical network
#' @examples
#' data <- matrix(sample(c(0,1), 100*4, replace = TRUE), 100, 4)
#' colnames(data) <- c("A", "A.B", "B", "C")
#' rownames(data) <- paste("E", 1:100, sep = "_")
#' res <- epiNEM(data, method = "exhaustive")
#' plot(res)
epiNEM <- function(filename="random", method="greedy",
                   nIterations=10, nModels=0,
                   random=list(single=4, double=1, reporters=100,
                               FPrate=0.1, FNrate=0.1, replicates=1),
                   plotsy=TRUE, ltype = "marginal", para = c(0.13, 0.05)) {

    ##--- Sanity checks ---#
    methods <- c("greedy", "exhaustive")
    if (!method %in% methods)
        stop("Please enter a valid method: 'greedy' or 'exhaustive'.\n")

    if ((method == "greedy") && (nIterations < 1))
        stop("Please enter a positive number of iterations\n")

    ##--- Process input ---#
    if (is.character(filename)){
        if (filename == "random") {
            ## Generate artificial data
            topology    <- CreateTopology(random$single, random$double)
            topology    <- unlist(topology, recursive=FALSE)
            extTopology <- ExtendTopology(topology$model, random$reporters)
            sortedData <-
                 GenerateData(topology$model, extTopology, random$FPrate,
                                        random$FNrate, random$replicates)

            mutants   <- rownames(topology$model)
            singleKOs <- sort(colnames(topology$model))
            experiments <- singleKOs
            doubleKOs <- mutants[!mutants %in% singleKOs]
        } else {
            ## Use data from file
            ## data        <- read.table(filename)
            ## sortedData  <- data[, order(unlist(sortedKOs))]
        }
    } else sortedData = filename

    mutants <- colnames(sortedData)
    singleKOs <- sort(unique(unlist(strsplit(mutants, '[.]'))))
    experiments <- singleKOs
    doubleKOs <- mutants[!mutants %in% singleKOs]
    ##-- Initialize data matrices (this takes care of replicate experiments) --
    D1 <- sapply(mutants, function(s)
        rowSums(sortedData[, colnames(sortedData) == s, drop=FALSE]))
    D0 <- sapply(mutants, function(s) sum(colnames(sortedData) == s)) - D1

    ##--- Optimize models ---#
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
        iterations <- sapply(1:nIterations, GreedyHillClimber,
                             experiments, D1, D0, mutants)

        index <- which.max(iterations["maxLlh",])
        reconstruct <- includeLogic(iterations[,index]$model,
                                    experiments, mutants)
        reconstruct <- unlist(reconstruct, recursive=FALSE)
        score <- sapply(reconstruct, Mll, D1, D0, ltype, para)
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
        extendedModels <- lapply(basicModels, includeLogic,
                                 experiments, unique(mutants))
        extendedModels <- unlist(unlist(extendedModels, recursive=FALSE),
                                 recursive = FALSE)
        uniqueModels   <- unique(lapply(extendedModels, function(e)
            identity((e["model"]))))
        uniqueModels <- uniqueModels#[1:length(uniqueModels)-1]

        if (length(doubleKOs) == 0) {
            print("No double perturbations available --> computing NEM")
            if (method == "exhaustive") { inference <- "search" }
            if (method == "greedy") { inference <- "nem.greedy" }
            options <-
                set.default.parameters(setdiff(unique(colnames(sortedData)),
                                               "time"))
            options$para <- para
            return(nem::nem(sortedData, inference = inference, control=options))
        }

        score <- sapply(uniqueModels, Mll, D1, D0, ltype, para)
        mll  <- unlist(score["mLL",])

        bestModel <- uniqueModels[[which.max(mll)]]$model
        allModels <- lapply(1:length(extendedModels), function(i)
            identity(extendedModels[[i]]$model))
        isBest    <- lapply(allModels, IsBestModel, bestModel)
        results   <- extendedModels[isBest == TRUE]
        results <- lapply(results, utils::modifyList,
                          list(score=mll[which.max(mll)]))
        results <- results[[1]]
        posterior <- AttachEGenes(score[,which.max(mll)]$posterior, experiments)

        results$EGeneset <- posterior

    }

    class(results) <- "epiNEM"
    
    return(results)
}

###--- HELPER FUNCTIONS ---###

#' Extending topology of normal "nem"
#' @param topology model of a topology from CreateTopology
#' @param nReporters number of effects reporters
#' @author Madeline Diekmann
#' @seealso CreateTopology
#' @export
#' @examples
#' topology <- CreateTopology(3, 1, force = TRUE)
#' topology <- unlist(unique(topology), recursive = FALSE)
#' extTopology <- ExtendTopology(topology$model, 100)
#' @return extended topology in which reporters are linked to pathway genes
ExtendTopology <- function(topology, nReporters) {
    reporters     <- unlist(lapply(1:nReporters,
                                   function(n) paste("reporter", n, sep="-")))
    linkedEffects <- sample(1:ncol(topology), nReporters, replace=TRUE)
    extTopology   <- sapply(linkedEffects,
                            function(e) lapply(1:ncol(topology),
                                               function(c) ifelse(e == c, 1, 0)))
    extTopology   <- matrix(unlist(extTopology), ncol=ncol(topology), byrow=TRUE)
    dimnames(extTopology) <- list(reporters, colnames(topology))
    return(extTopology)
}

#' Generate data from extended model
#' @param model model of a topology from CreateTopology
#' @param extTopology extended topology
#' @param FPrate false positive rate
#' @param FNrate false negative rate
#' @param replicates number of replicates
#' @author Madeline Diekmann
#' @seealso CreateTopology
#' @export
#' @examples
#' topology <- CreateTopology(3, 1, force = TRUE)
#' topology <- unlist(unique(topology), recursive = FALSE)
#' extTopology <- ExtendTopology(topology$model, 100)
#' sortedData <- GenerateData(topology$model, extTopology, 0.05, 0.13, 3)
#' @return data matrix
GenerateData <- function(model, extTopology, FPrate, FNrate, replicates) {
    perfectData <- extTopology %*% t(model)
    perfectData <- perfectData[, rep(seq_len(ncol(perfectData)), replicates)]
    fps <- sample(which(perfectData == 0), floor(sum(perfectData == 0)*FPrate))
    fns <- sample(which(perfectData == 1), floor(sum(perfectData == 1)*FNrate))
    noisyData   <- perfectData
    noisyData[c(fps, fns)] <- 1 - noisyData[c(fps, fns)]

    return(noisyData)
}

#' @noRd
IsBestModel <- function(thisModel, bestModel) {
    if (any(dim(thisModel) != dim(bestModel))) return(FALSE)
    else if (any(thisModel != bestModel)) return(FALSE)
    return(TRUE)
}

###--- END OF HELPER FUNCTIONS ---###
