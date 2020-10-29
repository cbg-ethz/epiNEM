#' Example data: simulation results
#'
#' Contains simulation results. How they were
#' aquired is explained in the vignette.
#' The data conists of a list of data matrices holding
#' sensitivity and specificity
#' (spec, sens) of network edges for the variious methods compared to
#' the ground truth, sensitivity and specificity (sens2, spec2)
#' of the expected data for epiNEM and Boolean NEMs and accuracy
#' of the inferred logics for both. The different methods are in the
#' rows and the columns denote the different independent simulation runs.
#' @docType data
#' @examples
#' data(sim)
#' @name sim
NA

#' Example data: epiNEM results for
#' the Sameith et al., 2015 knock-out screen
#'
#' The result of the epiNEM analysis of the data from
#' "http://www.holstegelab.nl/publications/
#' sv/signaling_redundancy/downloads/DataS1.txt".
#' The data consists of a list of matrices with the likelihoods (ll)
#' for each analysed triple of signalling genes and the inferred logic
#' (logic) for each triple. The signalling genes or modulators C are the
#' rows and the signalling genes from the double knock-downs are in the columns.
#' For details see the vignette.
#' @docType data
#' @examples
#' data(samscreen)
#' @name samscreen
NA

#' Example data: epiNEM results for
#' the Wageningen et al., 2010 knock-out screen
#' "http://www.holstegelab.nl/publications/GSTF_geneticinteractions/
#' downloads/del_mutants_limma.txt"
#'
#' The data consists of a list of matrices with the likelihoods (ll)
#' for each analysed triple of signalling genes and the inferred logic
#' (logic) for each triple. The signalling genes or modulators C are the
#' rows and the signalling genes from the double knock-downs are in the columns.
#' For details see the vignette.
#' @docType data
#' @examples
#' data(wagscreen)
#' @name wagscreen
NA

#' sig. of string interaction scores
#' for Sameith et al., 2015 data
#'
#' The data consists of a list including a vectors of pairs (for interactions)
#' and a corresponding list of interaction scores derived form the
#' string database.
#' For details see the vignette.
#' @docType data
#' @examples
#' data(sameith_string)
#' @name sameith_string
NA

#' sig. of string interaction scores
#' for van Wageningen et al., 2010 data
#'
#' The data consists of a list including a vectors of pairs (for interactions)
#' and a corresponding list of interaction scores derived form the
#' string database.
#' For details see the vignette.
#' @docType data
#' @examples
#' data(wageningen_string)
#' @name wageningen_string
NA

#' graph-based GO similarity scores, string GO annotations
#' for Sameith et al., 2015 data
#'
#' The data consists of lists including epiNEM identified and
#' general similarity scores and GO annotations for each triple.
#' For details see the vignette.
#' @docType data
#' @examples
#' data(sameith_GO)
#' @name sameith_GO
NA

#' graph-based GO similarity scores, string GO annotations
#' for van Wageningen et al., 2015 data
#'
#' The data consists of lists including epiNEM identified and
#' general similarity scores and GO annotations for each triple.
#' For details see the vignette.
#' @docType data
#' @examples
#' data(wageningen_GO)
#' @name wageningen_GO
NA

###--- MAIN SCRIPT ---###
#' Epistatic NEMs - main function.
#'
#' This function contains the inference
#' algorithm to learn logical networks from knock-down data including
#' double knock-downs.
#' @param filename A binary, tab-delimited matrix.
#' Columns: single and double knockdowns.
#' Rows: genes showing effect or not?
#' Default: random; artificial data is
#' generated to 'random' specifications
#' @param method greedy or exhaustive search.
#' Default: greedy
#' @param nIterations number of iterations.
#' Default: 10
#' @param nModels number of Models. Default: 0
#' @param random list specifying how the
#' data should be generated:
#' no. of single mutants, no. of
#' double mutants, no. of reporterGenes,
#' FP-rate, FN-rate, no. of replicates
#' @param ltype likelihood either
#' "marginal" or "maximum"
#' @param para false positive and
#' false negative rates
#' @param init adjacency matrix to initialise the greedy search
#' @author Madeline Diekmann
#' @seealso nem
#' @export
#' @import
#' stats
#' e1071
#' utils
#' @importFrom mnem mnem
#' @return List object with an adjacency matrix denoting the network,
#' the model of the silencing scheme (rows are knock-downs, columns
#' are signalling genes), a string with the inferred logial gates,
#' a column indices denoting position of logical gates, the log transformed
#' likelihood and the effect reporter distribution (rows are the signalling
#' genes including the null node).
#' @examples
#' data <- matrix(sample(c(0,1), 100*4, replace = TRUE), 100, 4)
#' colnames(data) <- c("A", "A.B", "B", "C")
#' rownames(data) <- paste("E", 1:100, sep = "_")
#' res <- epiNEM(data, method = "exhaustive")
#' plot(res)
epiNEM <- function(filename="random",
                   method="greedy",
                   nIterations=10,
                   nModels=0,
                   random=list(single=4,double=1,reporters=100,
                               FPrate=0.1,FNrate=0.1,replicates=1),
                   ltype = "marginal",
                   para = c(0.13, 0.05),
                   init = NULL) {

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
    } else if (method == "exhaustive" & length(singleKOs) == 5) {
        cat("\nThis is going to take a while...\n")
    }
    if (length(doubleKOs) == 0) {
        print("No double perturbations available --> computing NEM")
        return(mnem::mnem(sortedData, search=inference, fpfn=para, k=1))
    } else {
        if (method == "greedy") {
            print(paste("Progress (of ", nIterations, " iterations):",
                        sep = ''))
            iterations <- sapply(1:nIterations, GreedyHillClimber,
                                 experiments, D1, D0, mutants, init)

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
            if (ltype %in% "maximum") {
                names(results$score) <- "max"
            } else {
                names(results$score) <- "mLL"
            }
            results$EGeneset <- AttachEGenes(posterior, experiments)
            results$PostScore <- score[,which.max(mll)]$posterior
            cat("\n")
        } else if (method == "exhaustive") {
            basicModels    <- EnumerateModels(length(singleKOs), singleKOs)
            extendedModels <- lapply(basicModels, includeLogic,
                                     experiments, unique(mutants))
            extendedModels <- unlist(unlist(extendedModels, recursive=FALSE),
                                     recursive = FALSE)
            uniqueModels   <- unique(lapply(extendedModels, function(e)
                identity((e["model"]))))
            uniqueModels <- uniqueModels#[1:length(uniqueModels)-1]

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
            if (ltype %in% "maximum") {
                names(results$score) <- "max"
            } else {
                names(results$score) <- "mLL"
            }
            savescore <- score[,which.max(mll)]$posterior
            posterior <- AttachEGenes(score[,which.max(mll)]$posterior,
                                      experiments)
            results$EGeneset <- posterior
            results$PostScore <- savescore

        }

        class(results) <- "epiNEM"

        return(results)
    }
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
#' @return extended topology in which reporters
#' are linked to pathway genes
ExtendTopology <- function(topology, nReporters) {
    reporters     <- unlist(lapply(1:nReporters,
                                   function(n) paste("reporter", n, sep="-")))
    linkedEffects <- sample(1:ncol(topology), nReporters, replace=TRUE)
    extTopology   <-
        sapply(linkedEffects,
               function(e) lapply(1:ncol(topology),
                                  function(c) ifelse(e == c, 1, 0)))
    extTopology <-
        matrix(unlist(extTopology), ncol=ncol(topology), byrow=TRUE)
    dimnames(extTopology) <- list(reporters, colnames(topology))
    return(extTopology)
}

#' Generate data from extended model.
#'
#' Given a model created from
#' CreateTopology and ExtendTopology, this function creeates acorresponding
#' artificial data matrix, which is used as a ground truth for simulation
#' studies.
#' @param model model of a topology
#' from CreateTopology
#' @param extTopology extended topology
#' @param FPrate false positive rate
#' @param FNrate false negative rate
#' @param replicates number of replicates
#' @author Madeline Diekmann
#' @seealso CreateTopology
#' @export
#' @examples
#' topology <-
#' CreateTopology(3, 1, force = TRUE)
#' topology <-
#' unlist(unique(topology), recursive = FALSE)
#' extTopology <-
#' ExtendTopology(topology$model, 100)
#' sortedData <-
#' GenerateData(topology$model, extTopology, 0.05, 0.13, 3)
#' @return data matrix with effect reporters as rows and knock-downs
#' (including double kock-downs) as columns.
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

#' Compare algorithms.
#'
#' Compares different network reconstruction algorithm on
#' simulated data.
#' @param runs number simulation runs
#' @param do string vector of algorithms to
#' compare: e (epiNEM), n (Nested Effects Models),
#' b (B-NEM), p (PC algorithm), a (Aracne), e.g. c("e", "n", "p")
#' @param random list of false poitive rate FPrate, false negative rates
#' FNrate, number of single knock-downs single, number of double
#' knock-downs double, number of effect reporters reporters and number
#' of replicates replicates
#' @param maxTime TRUE if the algorithms are bound to a maximum running
#' time in respect to epiNEM
#' @param forcelogic if TRUE the randomly sampled ground truth network
#' includes a complex logic with probability 1
#' @param epinemsearch greedy or exhaustive search for epiNEM
#' @param bnemsearch genetic or greedy search for B-NEM
#' @param ... additional parameters
#' @author Martin Pirkl
#' @return returns list of specificity and sensitivity of inferred edges
#' (spec, sens) and inferred expected data (spec2, sens2) and accuracy
#' of logics (logics) and running time (time)
#' @export
#' @import
#' pcalg
#' minet
#' @importFrom mnem mnem transitive.reduction
#' @importFrom graph adj
#' @examples
#' res <- SimEpiNEM(runs = 1)
SimEpiNEM <- function(runs = 10, do = c("n", "e"),
                      random = list(FPrate = 0.1,
                                    FNrate = c(0.1, 0.5),
                                    single = 3, double = 1,
                                    reporters = 10,
                                    replicates = 2),
                      maxTime = FALSE,
                      forcelogic = TRUE,
                      epinemsearch = "greedy",
                      bnemsearch = "genetic", ...) {

    noiselvls <- random$FNrate

    spec <- sens  <- logics <- array(0, dim = c(2, runs, length(noiselvls)))

    sens2 <- spec2 <- time <- array(0, dim = c(5, runs, length(noiselvls)))

    logicgate <- matrix("", runs, length(noiselvls))

    edgenr <- matrix(0, runs, length(noiselvls))

    for (i in 1:runs) {

        print(paste("run ", i, sep = ""))

        for (j in 1:length(noiselvls)) {

            print(paste("noiselvl ", j, sep = ""))

            topology <- CreateTopology(random$single, random$double,
                                       force = forcelogic)

            topology <- unlist(unique(topology), recursive = FALSE)

            extTopology <- ExtendTopology(topology$model, random$reporters)

            sortedData <- GenerateData(topology$model, extTopology,
                                       random$FPrate, random$FNrate[j],
                                       random$replicates)

            logicgate[i, j] <- paste(topology$logics, collapse = "_")

            edgenr[i, j] <- sum(topology$origModel == 1)

            if ("e" %in% do) {
                print("epiNEM")

                start <- Sys.time()
                TriplModel <- epiNEM(filename = sortedData,
                                     method = epinemsearch,
                                     ...)
                time[1, i, j] <- difftime(Sys.time(), start, units = "secs")
                print(time[1, i, j])

                tp <- sum(topology$model == 1 & TriplModel$model == 1)
                tn <- sum(topology$model == 0 & TriplModel$model == 0)
                fp <- sum(topology$model == 0 & TriplModel$model == 1)
                fn <- sum(topology$model == 1 & TriplModel$model == 0)
                sens[1, i, j] <- tp/(tp+fn)
                spec[1, i, j] <- tn/(tn+fp)
                tp <- sum(topology$origModel == 1 & TriplModel$origModel == 1)
                tn <- sum(topology$origModel == 0 & TriplModel$origModel == 0)
                fp <- sum(topology$origModel == 0 & TriplModel$origModel == 1)
                fn <- sum(topology$origModel == 1 & TriplModel$origModel == 0)
                sens2[1, i, j] <- tp/(tp+fn)
                spec2[1, i, j] <- tn/(tn+fp)
                tp <- 0
                for (k in 1:length(topology$column)) {
                    for (l in 1:length(TriplModel$column)) {
                        if (topology$column[k] == TriplModel$column[l]) {
                            if (topology$logics[k] %in% TriplModel$logics[l]) {
                                tp <- tp + 1
                            }
                        }
                    }
                }
                logics[1, i, j] <- tp/(length(topology$logics) +
                                       length(TriplModel$logics) - tp)
                print(sens[1, i, j])
                print(spec[1, i, j])
                print(sens2[1, i, j])
                print(spec2[1, i, j])
                print(logics[1, i, j])

            }

            if ("b" %in% do) {
                print("B-NEM")

                gtn <- epiNEM2Bg(topology)

                fc <- cbind(Ctrl_vs_S = -1, epiNEM2Bg(sortedData))*(-1)

                bnemnoise <-
                    sample(1:nrow(fc), floor(nrow(fc)*random$FNrate[j]))

                fc[bnemnoise, 1] <- 0

                ers <- t(topology$model)*(-1)
                colnames(ers) <-
                    paste("S_vs_S_", gsub("\\.", "_", colnames(ers)), sep = "")
                ers <- cbind(Ctrl_vs_S = 1, ers)
                ers <- ers[, order(colnames(ers))]

                CNOlist <- dummyCNOlist(stimuli = "S",
                                        inhibitors = LETTERS[1:random$single],
                                        maxStim = 1, maxInhibit = 2,
                                        signals = LETTERS[1:random$single])

                parents <-
                    unique(unlist(
                        strsplit(colnames(sortedData)[
                            grep("\\.", colnames(sortedData))], "\\.")))

                nodes <-
                    unique(colnames(sortedData)[
                        -grep("\\.", colnames(sortedData))])

                child <- nodes[-which(nodes %in% parents)]

                sifMatrix <- NULL
                for (k in LETTERS[1:random$single]) {
                    sifMatrix <- rbind(sifMatrix, c("S", "1", k))
                    ##, c("S", "-1", k))
                    for (l in LETTERS[1:random$single]) {
                        if (k %in% l) { next() }
                        if (k %in% parents) {
                            sifMatrix <-
                                rbind(sifMatrix, c(k, "1", l), c(k, "-1", l))
                        } else {
                            sifMatrix <- rbind(sifMatrix, c(k, "1", l))
                        }
                    }
                }
                randfile <- paste("pkn_", as.numeric(Sys.time()), sep = "")
                write.table(sifMatrix, file = randfile, sep = "\t",
                            row.names = FALSE, col.names = FALSE, quote = FALSE)
                PKN <- readSIF(randfile)
                unlink(randfile)

                model <- preprocessing(CNOlist, PKN)

                initBstring <- absorption(rep(1, length(model$reacID)), model)

                if (maxTime) {
                    maxTime2 <- time[1, i, j]
                } else { maxTime2 <- Inf }

                start <- Sys.time()
                bga <- bnem(search = bnemsearch,
                            fc=fc,
                            CNOlist=CNOlist,
                            model=model,
                            initBstring=initBstring,
                            draw = FALSE,
                            verbose = FALSE,
                            popSize = popSize,
                            maxTime = maxTime2,
                            parallel = parallel
                            )
                time[2, i, j] <- difftime(Sys.time(), start, units = "secs")
                print(time[2, i, j])

                ers2 <-
                    computeFc(CNOlist,
                              t(simulateStatesRecursive(CNOlist,
                                                        model, bga$bString)))
                ers2 <- ers2[, unique(colnames(fc))]
                ers2 <- ers2[, order(colnames(ers2))]

                tp <- sum(ers == -1 & ers2 == -1)
                tn <- sum(ers == 0 & ers2 == 0)
                fn <- sum(ers == -1 & ers2 == 0)
                fp <- sum(ers == 0 & ers2 == -1)
                sens[2, i, j] <- tp/(tp+fn)
                spec[2, i, j] <- tn/(tn+fp)
                gtn2 <- abs(dnf2adj(gtn))
                if (length(grep("S", rownames(gtn2))) > 0) {
                    gtn2 <-
                        gtn2[-grep("S", rownames(gtn2)),
                             -grep("S", colnames(gtn2))]
                }
                gtn2 <- gtn2[order(rownames(gtn2)), order(colnames(gtn2))]
                res <- abs(dnf2adj(bga$graph))
                if (length(grep("S", rownames(res))) > 0) {
                    res <- as.matrix(res[-grep("S", rownames(res)),
                                         -grep("S", colnames(res))])
                }
                if (dim(res)[1] == 1) {
                    colnames(res) <- rownames(res) <- gsub(".*=", "", bga$graph)
                } else {
                    res <- res[order(rownames(res)), order(colnames(res))]
                }
                if (nrow(res) < nrow(gtn2)) {
                    res2 <-
                        rbind(cbind(res,
                                    matrix(0, nrow(res),
                                           nrow(gtn2) - nrow(res))),
                              matrix(0, nrow(gtn2) - nrow(res), ncol(gtn2)))
                    colnames(res2)[(ncol(res)+1):ncol(res2)] <-
                        colnames(gtn2)[which(!(colnames(gtn2)
                            %in% colnames(res)))]
                    rownames(res2)[(nrow(res)+1):nrow(res2)] <-
                        rownames(gtn2)[which(!(rownames(gtn2)
                            %in% rownames(res)))]
                    res2 <- res2[order(rownames(res2)), order(colnames(res2))]
                    res <- res2
                }
                diag(gtn2) <- diag(res) <- 0
                tp <- sum(gtn2 == 1 & res == 1)
                tn <- sum(gtn2 == 0 & res == 0)
                fn <- sum(gtn2 == 1 & res == 0)
                fp <- sum(gtn2 == 0 & res == 1)
                sens2[2, i, j] <- tp/(tp+fn)
                spec2[2, i, j] <- tn/(tn+fp)
                tp <- sum(bga$graph %in% gtn)
                logics[2, i, j] <- tp/(length(gtn) + length(bga$graph) - tp)
                print(sens[2, i, j])
                print(spec[2, i, j])
                print(sens2[2, i, j])
                print(spec2[2, i, j])
                print(logics[2, i, j])

                print(bga$graph)
                print(gtn)

            }

            if (any(c("n", "p", "a") %in% do)) {

                reddata <- sortedData[, -grep("\\.", colnames(sortedData))]
                gtnadj <- topology$origModel
                gtnadj <-
                    gtnadj[order(apply(gtnadj, 1, sum), decreasing = TRUE),
                           order(apply(gtnadj, 2, sum), decreasing = FALSE)]
                gtnadj[lower.tri(gtnadj)] <- gtnadj[upper.tri(gtnadj)]
                gtnadj <- gtnadj[order(rownames(gtnadj)),
                                 order(colnames(gtnadj))]
                eadj <- topology$origModel
                eadj <- eadj[order(rownames(eadj)), order(colnames(eadj))]
                reddata2 <- matrix(0, nrow(reddata)*random$replicates,
                                   length(unique(colnames(reddata))))
                for (k in 1:length(unique(colnames(reddata)))) {
                    reddata2[, k] <-
                        as.vector(reddata[,
                                          which(colnames(reddata) %in%
                                                unique(colnames(reddata))[k])])
                }
                colnames(reddata2) <- unique(colnames(reddata))

            }

            if ("n" %in% do) {
                print("NEM")

                start <- Sys.time()
                if (epinemsearch %in% "greedy") {
                    nemres <- mnem::mnem(reddata, k = 1)
                } else {
                    nemres <- mnem::mnem(reddata, k = 1, search = "exhaustive")
                }
                nadj <- mnem::transitive.reduction(nemres$comp[[1]]$phi)
                time[3, i, j] <- difftime(Sys.time(), start, units = "secs")
                print(time[3, i, j])

                tp <- sum(eadj  == 1 & nadj == 1)
                tn <- sum(eadj == 0 & nadj == 0)
                fp <- sum(eadj == 0 & nadj == 1)
                fn <- sum(eadj == 1 & nadj == 0)
                sens2[3, i, j] <- tp/(tp+fn)
                spec2[3, i, j] <- tn/(tn+fp)
                print(sens2[3, i, j])
                print(spec2[3, i, j])

            }

            if ("p" %in% do) {
                print("PCalg")

                start <- Sys.time()
                pc.fit <- pc(suffStat = list(C = cor(reddata2),
                                             n = nrow(reddata2)),
                             indepTest = gaussCItest,
                             ## indep.test: partial correlations
                             alpha=0.05, labels = colnames(reddata2),
                             verbose = FALSE)
                graph2adj <- function(gR) {
                    adj.matrix <- matrix(0,
                                         length(graph::nodes(gR)),
                                         length(graph::nodes(gR))
                                         )
                    rownames(adj.matrix) <- graph::nodes(gR)
                    colnames(adj.matrix) <- graph::nodes(gR)
                    for (i in 1:length(nodes(gR))) {
                        adj.matrix[graph::nodes(gR)[i],
                                   adj(gR,graph::nodes(gR)[i])[[1]]] <- 1
                    }

                    return(adj.matrix)
                }
                pcadj <- graph2adj(pc.fit@graph)
                time[4, i, j] <- difftime(Sys.time(), start, units = "secs")
                print(time[4, i, j])

                tp <- sum(gtnadj == 1 & pcadj == 1)
                tn <- sum(gtnadj  == 0 & pcadj == 0)
                fp <- sum(gtnadj == 0 & pcadj == 1)
                fn <- sum(gtnadj == 1 & pcadj == 0)
                sens2[4, i, j] <- tp/(tp+fn)
                spec2[4, i, j] <- tn/(tn+fp)
                print(sens2[4, i, j])
                print(spec2[4, i, j])

            }

            if ("a" %in% do) {
                print("Aracne")

                start <- Sys.time()
                ares <- build.mim(reddata2)
                ares <- aracne(ares)
                ares[which(ares > 0)] <- 1
                ares[which(ares < 0)] <- -1
                ares <- ares[order(rownames(ares)), order(colnames(ares))]
                nas <- which(is.na(ares) == TRUE)
                ares[nas] <- 0
                diag(ares) <- 0
                time[5, i, j] <- difftime(Sys.time(), start, units = "secs")
                print(time[5, i, j])

                tp <- sum(gtnadj == 1 & ares == 1)
                tn <- sum(gtnadj == 0 & ares == 0)
                fp <- sum(gtnadj == 0 & ares == 1)
                fn <- sum(gtnadj == 1 & ares == 0)
                sens2[5, i, j] <- tp/(tp+fn)
                spec2[5, i, j] <- tn/(tn+fp)
                print(sens2[5, i, j])
                print(spec2[5, i, j])

            }

        }

    }

    result <- list(sens = sens, spec = spec, sens2 = sens2,
                   spec2 = spec2, logics = logics, time = time)

    class(result) <- "epiSim"

    return(result)

}

#' Heatmap.
#'
#' Heatmap function based on the lattice package
#' more information: ?xyplot
#' @param x Matrix.
#' @param col Color. See brewer.pal.info for all available
#' color schemes.
#' @param colNA color for NAs; defaul is grey
#' @param coln Number of colors.
#' @param bordercol Border color.
#' @param borderwidth Border width.
#' @param breaks Defines the breaks in the color range. "sym"
#' makes the breaks symmetric around 0.
#' @param main Main title.
#' @param sub Subtitle.
#' @param dendrogram Draw dendrogram with "both", "col" or
#' "row", or do not draw with "none".
#' @param colorkey Draw colorkey list(space="left") or
#' list(space="right").
#' @param Colv Cluster columns (TRUE) or not (FALSE).
#' @param Rowv Cluster rows (TRUE) or not (FALSE).
#' @param xrot Rotate the column names by degree.
#' @param yrot Rotate the row names by degree.
#' @param shrink c(x,y) defines a range of size for the data
#' boxes from low to high.
#' @param cexCol Font size of column names.
#' @param cexRow Font size of row names.
#' @param cexMain Font size of main title.
#' @param cexSub Font size of subtitle.
#' @param colSideColors Defines a numeric vector to annotate
#' columns with different colors.
#' @param aspect "iso" for quadratic boxes or "fill" for
#' streched boxes.
#' @param contour TRUE adds a contour plot.
#' @param useRaster TRUE to add raster visuals
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param colSideColorsPos Place colSideColors at the "top" or "bottom".
#' @param clust p, s, or k for correlation clustering
#' @param clusterx Optional data matrix y with the same dimensions
#' as x. x is columns or rows are sorted by the cluster information of y.
#' @param \dots Optional arguments.
#' @author Martin Pirkl & Oscar Perpinan
#' at http://oscarperpinan.github.io/rastervis/
#' @return lattice object/matrix
#' @export
#' @import
#' lattice
#' latticeExtra
#' RColorBrewer
#' grDevices
#' @examples
#' x <- matrix(rnorm(50), 10, 5)
#' HeatmapOP(x, dendrogram = "both", aspect = "iso", xrot = 45)
HeatmapOP <-
    function(x, col = "RdYlGn", colNA = "grey", coln = 11, bordercol = "grey",
             borderwidth = 0.1, breaks = "sym",
             main = "",
             sub = "",
             dendrogram = "none", colorkey = list(space = "right"), Colv = TRUE,
             Rowv = TRUE, xrot = 90, yrot = 0, shrink = c(1,1), cexCol = 1,
             cexRow = 1, cexMain = 1, cexSub = 1,
             colSideColors = NULL, aspect = "fill",
             contour = FALSE, useRaster = FALSE, xlab = NULL, ylab = NULL,
             colSideColorsPos = "top", clust = NULL, clusterx = NULL, ...) {
        if (nrow(x)==1) {
            Rowv <- FALSE
        }
        if (ncol(x)==1) {
            Colv <- FALSE
        }
        if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
            x <- matrix(rnorm(100), 10, 10)
            main <- "max value equals min value"
            sub <- "random matrix plotted"
        }
        if (is.null(breaks)) {
            breaks <- seq(min(x, na.rm = TRUE),max(x, na.rm = TRUE),
            (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/45)
        }
        if ("sym" %in% breaks) {
            breaks <- max(abs(x), na.rm = TRUE)
            breaks2 <- breaks/45
            breaks <- seq(-breaks,breaks,breaks2)
        }
        if (length(breaks) == 1) {
            at <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE),
            (max(x, na.rm = TRUE)-min(x, na.rm = TRUE))/breaks)
        } else {
            at <- breaks
            x[x < breaks[1]] <- breaks[1]
            x[x > breaks[length(breaks)]] <- breaks[length(breaks)]
        }
        na.idx <- NULL
        if ((Colv | Rowv) & any(is.na(x) == TRUE)) {
            na.idx <- which(is.na(x) == TRUE)
            x[na.idx] <- mean(x, na.rm = TRUE)

        }
        dd.col <- NULL
        if (dendrogram == "row" | dendrogram == "both" & !is.null(colorkey)) {
            colorkey = list(space = "left")
        }
        if (Colv) {
            if (is.null(clust)) {
                if (is.null(clusterx)) {
                    dd.col <- as.dendrogram(hclust(dist(t(x))))
                } else {
                    dd.col <- as.dendrogram(hclust(dist(t(clusterx))))
                }
            } else {
                if (is.null(clusterx)) {
                    cor.tmp <- cor(x, method = clust)
                } else {
                    cor.tmp <- cor(clusterx, method = clust)
                }
                cor.tmp[which(is.na(cor.tmp) == TRUE)] <- 0
                dd.col <- as.dendrogram(hclust(as.dist(1-abs(cor.tmp))))
            }
            col.ord <- order.dendrogram(dd.col)
        } else {
            if (dendrogram %in% "both") {
                dendrogram <- "row"
            }
            if (dendrogram %in% "col") {
                dendrogram <- "none"
            }
            col.ord <- 1:ncol(x)
            legend = NULL
        }
        if (Rowv) {
            if (is.null(clust)) {
                if (is.null(clusterx)) {
                    dd.row <- as.dendrogram(hclust(dist(x)))
                } else {
                    dd.row <- as.dendrogram(hclust(dist(clusterx)))
                }
            } else {
                if (is.null(clusterx)) {
                    cor.tmp <- cor(t(x), method = clust)
                } else {
                    cor.tmp <- cor(t(clusterx), method = clust)
                }
                cor.tmp[which(is.na(cor.tmp) == TRUE)] <- 0
                dd.row <- as.dendrogram(hclust(as.dist(1-abs(cor.tmp))))
            }
            row.ord <- rev(order.dendrogram(dd.row))
        } else {
            if (dendrogram %in% "both") {
                dendrogram <- "col"
            }
            if (dendrogram %in% "row") {
                dendrogram <- "none"
            }
            row.ord <- 1:nrow(x)
            legend = NULL
        }
        if (!is.null(na.idx)) x[na.idx] <- NA
        add <- list(rect = list(col = "transparent",
                                fill = colSideColors[sort(col.ord)]))
        myTheme <- custom.theme(region=RColorBrewer::brewer.pal(n=coln, col))
        if (dendrogram != "none") {
            if (dendrogram == "both") {
                if (colSideColorsPos %in% "bottom") {
                    if (!is.null(colSideColors)) {
                        size <- 2
                    } else {
                        size <- 10
                    }
                    legend <- list(bottom =
                                       list(fun = dendrogramGrob,
                                            args =
                                                list(x = dd.col, ord = col.ord,
                                                     side = "top",
                                                     size = size,
                                                     add = add,
                                                     type = "rectangle")),
                                   right =
                                       list(fun = dendrogramGrob,
                                            args =
                                                list(x = dd.row,
                                                     side = "right",
                                                     size = 10)))
                }
                if (colSideColorsPos %in% "top") {
                    legend <- list(top =
                                       list(fun = dendrogramGrob,
                                            args =
                                                list(x = dd.col, ord = col.ord,
                                                     side = "top",
                                                     size = 10,
                                                     add = add,
                                                     type = "rectangle")),
                                   right =
                                       list(fun = dendrogramGrob,
                                            args =
                                                list(x = dd.row,
                                                     side = "right",
                                                     size = 10)))
                }
            }
            if (dendrogram == "row") {
                if (!is.null(colSideColors)) {
                    if (is.null(dd.col)) {
                        col.ord <- 1:length(colSideColors)
                    }
                    if (is.null(clusterx)) {
                        dd.col <- as.dendrogram(hclust(dist(t(x))*0))
                    } else {
                        dd.col <- as.dendrogram(hclust(dist(t(clusterx))*0))
                    }
                    if (colSideColorsPos %in% "bottom") {
                        legend <-
                            list(right =
                                     list(fun = dendrogramGrob,
                                          args =
                                              list(x = dd.row,
                                                   side = "right",
                                                   size = 10,
                                                   size.add= 0.5)),
                                 bottom =
                                     list(fun = dendrogramGrob,
                                          args = list(x = dd.col, ord = col.ord,
                                                      side = "top", size = 1,
                                                      size.add= 1,
                                                      add = add,
                                                      type = "rectangle")))
                    }
                    if (colSideColorsPos %in% "top") {
                        legend <-
                            list(right =
                                     list(fun = dendrogramGrob,
                                          args =
                                              list(x = dd.row,
                                                   side = "right",
                                                   size = 10, size.add= 0.5)),
                                 top =
                                     list(fun = dendrogramGrob,
                                          args = list(x = dd.col, ord = col.ord,
                                                      side = "top", size = 1,
                                                      size.add= 1,
                                                      add = add,
                                                      type = "rectangle")))
                    }
                } else {
                    legend <-
                        list(right =
                                 list(fun = dendrogramGrob,
                                      args =
                                          list(x = dd.row,
                                               side = "right",
                                               size = 10, size.add= 0.5)))
                }
            }
            if (dendrogram == "col") {
                if (colSideColorsPos %in% "bottom") {
                    if (!is.null(colSideColors)) {
                        size <- 2
                    } else {
                        size <- 10
                    }
                    legend <-
                        list(bottom =
                                 list(fun = dendrogramGrob,
                                      args = list(x = dd.col, ord = col.ord,
                                                  side = "top", size = size,
                                                  size.add= 0.5,
                                                  add = add,
                                                  type = "rectangle")))
                }
                if (colSideColorsPos %in% "top") {
                    legend <-
                        list(top =
                                 list(fun = dendrogramGrob,
                                      args = list(x = dd.col, ord = col.ord,
                                                  side = "top", size = 10,
                                                  size.add= 0.5,
                                                  add = add,
                                                  type = "rectangle")))
                }
            }
        } else {
            if (!is.null(colSideColors)) {
                if (is.null(dd.col)) {
                    col.ord <- 1:length(colSideColors)
                }
                if (is.null(clusterx)) {
                    dtx <- dist(t(x))
                    dtx[is.na(dtx)] <- 0
                    dd.col <- as.dendrogram(hclust(dtx*0))
                } else {
                    dist(t(clusterx))
                    dtcx[is.na(dtcx)] <- 0
                    dd.col <- as.dendrogram(hclust(dtcx*0))
                }
                if (colSideColorsPos %in% "bottom") {
                    legend <-
                        list(bottom =
                                 list(fun = dendrogramGrob,
                                      args = list(x = dd.col, ord = col.ord,
                                                  side = "top", size = 1,
                                                  size.add= 1,
                                                  add = add,
                                                  type = "rectangle"))
                             )
                }
                if (colSideColorsPos %in% "top") {
                    legend <-
                        list(top =
                                 list(fun = dendrogramGrob,
                                      args = list(x = dd.col, ord = col.ord,
                                                  side = "top", size = 1,
                                                  size.add= 1,
                                                  add = add,
                                                  type = "rectangle"))
                             )
                }
            } else {
                legend <- NULL
            }
        }
        d <- t(x[row.ord, col.ord,drop=FALSE])
        d <- d[, ncol(d):1, drop=FALSE]
        if (contour) {
            region <- TRUE
            col.regions <- terrain.colors(100)
        }
        ##  print(p, position=c(0,ypct-0.05,1,1), more=TRUE)
        ##  print(p2, position=c(0,0,1,ypct+0.05))
        if (is.null(rownames(x))) {
            ytck <- list(cex = 0, rot = 0, at = NULL)
        } else {
            ytck <- list(cex = cexRow, rot = yrot)
        }
        if (is.null(colnames(x))) {
            xtck <- list(cex = 0, rot = 0, at = NULL)
        } else {
            xtck <- list(cex = cexCol, rot = xrot)
        }
        levelplot(d,
                  main = list(label = main, cex = cexMain),
                  sub = list(label = sub, cex = cexSub),
                  aspect = aspect,
                  xlab=xlab,
                  ylab=ylab,
                  scales = list(x = xtck, y = ytck, tck = c(1,0)),
                  par.settings=myTheme,
                  border=bordercol,
                  border.lwd=borderwidth,
                  shrink=shrink,
                  legend = legend,
                  at = at,
                  colorkey = colorkey,
                  contour = contour,
                  panel = if (useRaster) {
                              function(...) {
                                  panel.fill(col = colNA)
                                  panel.levelplot.raster(...)
                              }
                          } else {
                              function(...) {
                                  panel.fill(col = colNA)
                                  panel.levelplot(...)
                              }
                          }
                  )
    }

#' Analyse large double knock-out screen.
#'
#' This function is used to analyse knock-out screens with multiple
#' double and single knock-outs combined in one data set.
#' @param data data matrix containing multiple single and double kock-downs
#' in columns and effect reporters in the rows
#' @param ... additional parameters, e.g. for the main epiNEM function
#' @export
#' @author Martin Pirkl
#' @return list object with vectors of double knock-downs, single knock-downs
#' and two matrices with doubles in the columns and singles in the rows. The
#' first matrix denotes the respective logical gate for the triple and the
#' second matrix the log-likelihood
#' @examples
#' data <- matrix(sample(c(0,1), 100*9, replace = TRUE), 100, 9)
#' colnames(data) <- c("A.B", "A.C", "B.C", "A", "B", "C", "D", "E", "G")
#' rownames(data) <- paste("E", 1:100, sep = "_")
#' res <- epiScreen(data)
epiScreen <- function(data, ...) {

    dataBin <- data

    doubles <- colnames(dataBin)[grep("\\.", colnames(dataBin))]

    if (length(grep("vs", doubles)) > 0) {
        doubles <- sort(doubles[-grep("vs", doubles)])
    } else { doubles <- sort(doubles) }

    doubles.genes <- unique(unlist(strsplit(doubles, "\\.")))

    if (length(grep("\\.", colnames(dataBin))) > 0) {
        singles <- colnames(dataBin)[-grep("\\.", colnames(dataBin))]
    } else { singles <- sort(singles) }

    singles <- unique(sort(singles))

    llmat <- logicmat <- matrix(0, length(singles), length(doubles))

    rownames(llmat) <- rownames(logicmat) <- singles

    colnames(llmat) <- colnames(logicmat) <- doubles

    globalgenes <- which(apply(dataBin, 1, max) == 1)

    targets <- list()

    for (i in doubles) {
        targets[[i]] <- list()
        print(i)
        doubles.singles <- unlist(strsplit(i, "\\."))
        egenes <-
            which(apply(dataBin[, which(colnames(dataBin) %in%
                                        c(i, doubles.singles))], 1, max) == 1)
        for (j in singles) {
            print(j)
            if (j %in% doubles.singles) { next() }

            dataTmp <- dataBin[, grep(paste(
                paste("^", c(i, j, doubles.singles), "$", sep = ""),
                collapse = "|"),
                colnames(dataBin))]

            dataTmp <- dataTmp[egenes, ]

            i1 <- which(singles %in% j)
            i2 <- which(doubles %in% i)

            if (!(is.null(dim(dataTmp)))) {

                if (any(dataTmp[, j] != 0)) {

                    epires <- epiNEM(dataTmp, method = "exhaustive")

                    targets[[i]][[j]] <-
                        rownames(dataTmp)[which(apply(epires$PostScore, 1,
                                                      which.max) ==
                                                which(colnames(epires$PostScore)
                                                      %in% j))]

                    tmp <- epires$logics
                    if ("OR" %in% tmp) {
                        if (sum(epires$origModel[, j]) != 2) {
                            tmp <- "NOEPI"
                        } else {
                            if (all(tmp %in% "OR")) {
                                tmp <- "OR"
                            } else {
                                tmp <- tmp[which(!(tmp %in% "OR"))]
                            }
                        }
                    }

                    logicmat[i1, i2] <- tmp
                    llmat[i1, i2] <- epires$score

                } else {

                    logicmat[i1, i2] <- "UNCON"
                    llmat[i1, i2] <- -Inf

                }

            } else {

                logicmat[i1, i2] <- "UNCON"
                llmat[i1, i2] <- -Inf

            }

        }

    }

    results <- list(doubles = doubles, singles = singles,
                    logic = logicmat, ll = llmat, targets = targets)

    class(results) <- "epiScreen"

    return(results)

}

#' Gate visualisation.
#'
#' Plots logical gate data annotation. The 8 heatmaps visualize what perfect
#' data would look like in respective to each logical gate. Perfect data is
#' equivalent to Boolean truth tables.
#' @references \url{https://en.wikipedia.org/wiki/Boolean_algebra}
#' @export
#' @author Martin Pirkl
#' @return plot of heatmaps showing the silencing scheme (=expected data,
#' truth tables)
#' @examples
#' epiAnno()
epiAnno <- function() {
    oldw <- getOption("warn")
    options(warn=-1)
    a1 <- HeatmapOP(matrix(c(1,-1,1,-1,1,1, 1, 1, 1), 3, 3,
                           dimnames = list(c("A", "B", "A.B"), LETTERS[1:3])),
                    Colv = FALSE, Rowv = FALSE,
                    main = "OR", col = "Greys", sub = "", colorkey = NULL)
    a2 <- HeatmapOP(matrix(c(1,-1,1,-1,1,1, -1, -1, 1), 3, 3,
                           dimnames = list(c("A", "B", "A.B"), LETTERS[1:3])),
                    Colv = FALSE, Rowv = FALSE,
                    main = "AND", col = "Greys", sub = "", colorkey = NULL)
    a3 <- HeatmapOP(matrix(c(1,-1,1,-1,1,1, 1, -1, -1), 3, 3,
                           dimnames = list(c("A", "B", "A.B"), LETTERS[1:3])),
                    Colv = FALSE, Rowv = FALSE,
                    main = "B masks effect of A", col = "Greys", sub = "",
                    colorkey = NULL)
    a4 <- HeatmapOP(matrix(c(1,-1,1,-1,1,1, -1, 1, -1), 3, 3,
                           dimnames = list(c("A", "B", "A.B"), LETTERS[1:3])),
                    Colv = FALSE, Rowv = FALSE,
                    main = "A masks effect of B", col = "Greys", sub = "",
                    colorkey = NULL)
    a5 <- HeatmapOP(matrix(c(1,-1,1,-1,1,1, 1, 1, -1), 3, 3,
                           dimnames = list(c("A", "B", "A.B"), LETTERS[1:3])),
                    Colv = FALSE, Rowv = FALSE,
                    main = "XOR", col = "Greys", sub = "", colorkey = NULL)
    a6 <- HeatmapOP(matrix(c(1,-1,1,-1,1,1, -1, 1, 1), 3, 3,
                           dimnames = list(c("A", "B", "A.B"), LETTERS[1:3])),
                    Colv = FALSE, Rowv = FALSE,
                    main = "No epistasis", col = "Greys", sub = "",
                    colorkey = NULL)
    a7 <- HeatmapOP(matrix(c(1,-1,1,-1,1,1, 1, -1, 1), 3, 3,
                           dimnames = list(c("A", "B", "A.B"), LETTERS[1:3])),
                    Colv = FALSE, Rowv = FALSE,
                    main = "No epistasis", col = "Greys", sub = "",
                    colorkey = NULL)
    a8 <- HeatmapOP(matrix(c(1,-1,1,-1,1,1, -1, -1, -1), 3, 3,
                           dimnames = list(c("A", "B", "A.B"), LETTERS[1:3])),
                    Colv = FALSE, Rowv = FALSE,
                    main = "No epistasis (unconnected)", col = "Greys",
                    sub = "", colorkey = NULL)
    print(a5, position = c(0,0, .25, .5), more = TRUE)
    print(a6, position = c(.25,0, .5, .5), more = TRUE)
    print(a7, position = c(.5,0, .75, .5), more = TRUE)
    print(a8, position = c(.75,0, 1, .5), more = TRUE)
    print(a1, position = c(0,.5, .25, 1), more = TRUE)
    print(a2, position = c(.25,.5, .5, 1), more = TRUE)
    print(a3, position = c(.5,.5, .75, 1), more = TRUE)
    print(a4, position = c(.75,.5, 1, 1))
    options(warn=oldw)
}

###--- END OF HELPER FUNCTIONS ---###


col = "RdYlGn"; colNA = "grey"; coln = 11; bordercol = "grey";
             borderwidth = 0.1; breaks = "sym";
             main = "";
             sub = "";
             dendrogram = "none"; colorkey = list(space = "right"); Colv = TRUE;
             Rowv = TRUE; xrot = 90; yrot = 0; shrink = c(1,1); cexCol = 1;
             cexRow = 1; cexMain = 1; cexSub = 1;
             colSideColors = NULL; aspect = "fill";
             contour = FALSE; useRaster = FALSE; xlab = NULL; ylab = NULL;
             colSideColorsPos = "top"; clust = NULL; clusterx = NULL;
