## functions needed in order to run LogicNEM

#' Create a random graph
#' @description Returns a model graph with randomly sampled edges.
#' Every possible edge has a probability to exist in the graph.
#' @param pathwayGenes vector of genes in the pathway
#' @param edgeProb probability of random edge
#' @export
#' @examples
#' graph <- CreateRandomGraph(c("Ikk1", "Ikk2", "RelA"))
#' @return adjacency matrix
CreateRandomGraph <- function(pathwayGenes, edgeProb=0.5) {
    size     <- length(pathwayGenes)
    model    <- diag(size)
    nEdges   <- size * (size - 1)
    model[which(model == 0)] <-
        sample(0:1, nEdges, TRUE, c(1 - edgeProb, edgeProb))
    dimnames(model) <- list(pathwayGenes, pathwayGenes)
    return(model)
}

get_col <- function(i, m) ((i-1) %/% nrow(m)) + 1
get_row <- function(i, m) (i-1) %% nrow(m) + 1

#' Evaluation of graphs
#' @param Phi model to be evaluated
#' @param D1 observed data matrix
#' @param D0 complementary D1
#' @param ltype likelihood type either "marginal" or "maximum"
#' @param para false positive and false negative rates
#' @description Computes marginal log-likelihood
#' for model Phi given observed data matrix D1
#' @export
#' @examples
#' Phi <- matrix(sample(c(0,1), 9, replace = TRUE), 3, 3)
#' data <- matrix(sample(c(0,1), 3*10, replace = TRUE), 10, 3)
#' rownames(Phi) <- colnames(Phi) <- colnames(data) <- c("Ikk1", "Ikk2", "RelA")
#' score <- Mll(Phi, D1 <- data, D0 <- 1 - data)
#' @return list with likelihood poster probability, egene positions
Mll <- function(Phi, D1, D0, ltype = "marginal", para = c(0.13, 0.05)) {
    ## Computes marginal log-likelihood for model Phi,
    ## observed data matrix D1, and
    ## complementary data matrix D0
    ## Function adapted from NEM package
    if (is.matrix(Phi)) {
        Phi2 <- Phi
    } else {
        Phi2 <- as.matrix(Phi$model)
    }
    Phi2 <- Phi2[colnames(D1), ]
    L <- para[1]^(D1 %*% (1 - Phi2)) *
                     (1 - para[1])^(D0 %*% (1 - Phi2)) *
                                       (1 - para[2])^(D1 %*% Phi2) *
                                                         para[2]^(D0 %*% Phi2)
    posterior <- L / (rowSums(L))
    LLperGene <- log(rowSums(L))
    mLL       <- sum(LLperGene)
    theta     <- apply(posterior, 1, function(e) e == max(e))
    mappos    <- apply(theta, 1, which)
    ## map instead of mll
    if (ltype %in% "maximum") {
        mLL <- sum(apply(L, 1, max))
    }
    return(list(posterior=posterior,
                LLperGene=LLperGene,
                mLL=mLL,
                mappos=mappos,
                para=para))
}

#' @noRd
getStarters <- function(mutants, experiments){
    ## get positions of doublemutants
    mutantslist <- strsplit(mutants, ".", fixed=TRUE)
    doublepos <- c()
    for (i in 1:length(mutantslist)){
        if (length(unlist(mutantslist[i])) == 2)
            doublepos <- c(doublepos, i)
    }
    starters <- c()
    dbvec <- c()
    ## create startStates for each knockout
    for (i in 1:length(experiments)){
        if (i %in% doublepos){
            for (j in experiments){
                if (j %in% unlist(mutantslist[doublepos]))
                    dbvec <- c(dbvec, 1)
                else dbvec <- c(dbvec, 0)
            }
            starters <- c(starters, dbvec)
        }
        z1 <- i-1
        z1 <- rep(0,z1)
        z2 <- length(experiments)-i
        z2 <- rep(0, z2)
        vector <- c(z1, 1, z2)
        starters <- c(starters, vector)
    }
    if (length(dbvec)==0){
        for (j in experiments){
            if (j %in% unlist(mutantslist[length(mutantslist)]))
                dbvec <- c(dbvec, 1)
            else dbvec <- c(dbvec, 0)
        }
        starters <- c(starters, dbvec)
    }
    return(starters)
}

#' Create an extended adjacency matrix
#' @description extend adjacency matrices taking cycles and logics into account.
#' For every given start state, the final state
#' is computed yu using BoolNet.
#' @param network network created by BoolNet from file
#' @param mutants vector of single knockouts
#' @param experiments vector of all knockouts
#' @importFrom BoolNet getPathToAttractor
#' @importFrom stats runif
#' @export
#' @examples
#' library(BoolNet)
#' data(cellcycle)
#' extModel <- CreateExtendedAdjacency(cellcycle,
#' c(cellcycle$genes, "CycD.Rb"), cellcycle$genes)
#' @return extended adjacency matrix
CreateExtendedAdjacency <- function(network, mutants, experiments){
    starters <- matrix(getStarters(mutants, experiments),
                       length(mutants), byrow=TRUE)
    rownames(starters) <- mutants
    colnames(starters) <- experiments

    ## get steady state of each startState
    getMyAttractor <- function(e, network, includeAttractorStates) {
        network$fixed[which(e == 1)] <- 1
        return(getPathToAttractor(network, e, includeAttractorStates))
    }
    data <- apply(starters, 1, getMyAttractor, network,
                  includeAttractorStates="all")

    ## get observed data by taking final state from all the random startStates.
    extadj <- data[[1]][nrow(data[[1]]),]
    for (i in 2:length(data)){
        d <- data[[i]][nrow(data[[i]]),]
        extadj <- rbind(extadj, d);
    }
    rownames(extadj) <- mutants
    return(extadj)
}

#' @noRd
includeLogic <- function(adj, experiments, mutants){
    if (is.list(adj)) {
        adj=matrix(unlist(adj),length(experiments), byrow=FALSE)
    }
    rownames(adj) <- experiments
    colnames(adj) <- rownames(adj)
    diag(adj)=0
    adj <- adj[order(rownames(adj)), order(colnames(adj))]
    mutantslist <- strsplit(mutants, ".", fixed=TRUE)
    doublepos <- c()
    for (i in 1:length(mutantslist)) {
        if (length(unlist(mutantslist[i])) == 2) {
            doublepos <- c(doublepos, i)
        }
    }
    ## look for suitable triples and create logic vector depending
    ## on the relationship of the parents
    notriples <- 0
    logic <- c()
    liste <- list()
    column <- c()
    for (c in 1:length(experiments)) {
        singles <- c()
        parents <- names(which(adj[,c]==1))
        if ((length(parents) == 2) &&
            ##check whether double mutant is available in the data
            (parents[1] %in% unlist(mutantslist[doublepos])) &&
            parents[2] %in% unlist(mutantslist[doublepos])) {
            for (i in 1:length(parents)) {
                singles <- cbind(singles, parents[i])
            }
            notriples <- notriples+1
            relation <- adj[singles, singles]
            diag(relation) <- 0
            NOT2 <- paste(parents[2], "masks the effect of", parents[1])
            NOT1 <- paste(parents[1], "masks the effect of", parents[2])
            if (sum(relation) == 0) {
                logic <- c("OR", "XOR", "AND", NOT1, NOT2)
            }
            if (relation[2] == 1) {
                logic <- c(NOT2)
            }
            if (relation[3] == 1) {
                logic <- c(NOT1)
            }
            liste[[notriples]] <- logic
            column <- cbind(column, c)
        }
        if (length(parents) == 2 &
            !(all(parents %in% unlist(mutantslist[doublepos])))) {
            notriples <- notriples + 1
            liste[[notriples]] <- "OR"
            column <- cbind(column, c)
        }
        if (length(parents) > 2) { 
            notriples <- notriples + 1
            liste[[notriples]] <- "OR"
            column <- cbind(column, c)
        }
        if (length(parents) == 1) {
            notriples <- notriples + 1
            liste[[notriples]] <- "OR"
            column <- cbind(column, c)
        }
    }
    if (length(liste) > 0) {
        logicmatrix <- as.matrix(expand.grid(liste))
        ## create logics file from adjacency matrix
        ## using logics provided by logic vector
        ## ready for using BoolNet
        randomnames <- sort(runif(nrow(logicmatrix)))
        for (modelno in 1:nrow(logicmatrix)) {
            lo <- 0
            if (!dir.exists("temp")) {
                dir.create("temp")
            }
            ## change that in the future:
            path <- paste("temp/outfile_", randomnames[modelno], ".txt", sep="")
            network <- character()
            countline <- 1
            network[countline] <- "targets, factors"
            for (c in experiments) {
                tmp <- NULL
                count <- 1
                if (sum(adj[, which(colnames(adj) %in% c)]) == 0) {
                    tmp <- paste(tmp, paste(c, ", ", c, sep=""), sep = "")
                } else {
                    tmp <- paste(tmp, paste(c, ", ", sep=""), sep = "")
                }
                count2 <- 0
                count3 <- 0
                for (r in experiments) {
                    if (adj[r,c]==1) {
                        if ((count==1) &&
                            (which(rownames(adj)==c) %in% column)) { 
                            help <- r
                            count <- count+1
                            if (sum(adj[, c]) == 1) {
                                tmp <- paste(tmp, paste(r, sep=""), sep = "")
                                lo <- lo + 1
                            }
                        }
                        else if ((count == 2) &&
                                 (which(rownames(adj)==c) %in% column)) {
                            NOT2 <- paste(r, "masks the effect of", help)
                            NOT1 <- paste(help, "masks the effect of", r)
                            if (count3 < 1) {
                                lo <- lo+1
                            }
                            if (logicmatrix[modelno, lo]=="OR" & count3 == 0) {
                                tmp <-
                                    paste(tmp, paste(help, " | ", r, sep=""),
                                          sep = "")
                                count3 <- count3 + 1
                            } else if
                            (logicmatrix[modelno, lo]=="OR" & count3 > 0) {
                                tmp <-
                                    paste(tmp,
                                          paste(" | ", help, " | ", r, sep=""),
                                          sep = "")
                            } else if
                            (logicmatrix[modelno, lo]=="AND") {
                                tmp <-
                                    paste(tmp,
                                          paste("(", help, " & ", r, ")",
                                                sep=""),
                                          sep = "")
                            } else if
                            (logicmatrix[modelno, lo]=="XOR") {
                                tmp <-
                                    paste(tmp,
                                          paste("( ", help, " & ! ", r
                                               ,") | (", r, " & ! ", help, ")",
                                                sep = ""),
                                          sep = "")
                                ## help refers to the first element
                            } else if
                            (logicmatrix[modelno, lo]==NOT2) {
                                tmp <-
                                    paste(tmp,
                                          paste("(", help, " & ! ", r, ")",
                                                sep=""),
                                          sep = "")
                            } else if
                            (logicmatrix[modelno, lo]==NOT1) {
                                tmp <-
                                    paste(tmp,
                                          paste("(", r, " & ! ", help, ")",
                                                sep=""),
                                          sep = "")
                            } else {
                                count2 <- count2 + 1
                                if
                                (count2 <
                                 sum(adj[, which(colnames(adj) %in% c)])) {
                                    tmp <- paste(tmp,
                                                 paste(r, " | ",
                                                       sep=""),
                                                 sep = "")
                                } else {
                                    tmp <- paste(tmp, paste(r,
                                                            sep=""),
                                                 sep = "")
                                    lo <- lo + 1
                                }
                            }
                        }
                    }
                }
                countline <- countline + 1
                network[countline] <- tmp
            }
            write(network, file = path)
        }
        test <- lapply(1:nrow(logicmatrix),
                       function(x) getExtendedAdjacency(x,
                                                        logicmatrix,
                                                        column,
                                                        adj,
                                                        mutants,
                                                        experiments,
                                                        sort(randomnames)))
        return(test)
    }
}

## to do: sehr unschoen!!!
#' create with logics extended adjacency matrix
#' @importFrom BoolNet loadNetwork
#' @noRd
getExtendedAdjacency <-function(modelno, logicmatrix,
                                column, adj, mutants,
                                experiments, randomnames) {
    randomnames <- sort(randomnames)
    path <- paste("temp/outfile_", randomnames[modelno], ".txt", sep="")
    network <- loadNetwork(path)
    extadj2 <- CreateExtendedAdjacency(network, unique(mutants), experiments)
    unlink(path)
    return(list(list(origModel=adj, model=extadj2,
                     logics=logicmatrix[modelno,], column=column)))
}

#' @noRd
AttachEGenes <- function(posterior, experiments){
    maxpost <-
        lapply(1:nrow(posterior),
               function(x) length(which(posterior[x,]==max(posterior[x,]))))
    attachedEs <- matrix()
    names <- c()
    attachedEsADD <- 0
    namesADD <- c()
    Egeneset <- matrix()
    Egeneset <- NULL
    EgeneADD <- 0
    EgeneADD <- NULL
    j <- 1
    k <- 1
    Epos <- c()
    for (i in 1:length(maxpost)){
        if (maxpost[i]==1){
            Epos <- cbind(Epos, i)
            attachedEs[j] <- names(which.max(posterior[i,]))
            names[j] <- rownames(posterior)[i]
            j <- j+1
        } else {
            ## for (e in 1:unlist(maxpost[i])){
            attachedEsADD <- attachedEsADD+1
            ##   namesADD[k] <- rownames(posterior)[i]
            ##   k <- k+1
            ## }
        }
    }
    attachedEs=as.matrix(attachedEs)
    rownames(attachedEs) <- names
    attachedEsADD=as.matrix(attachedEsADD)
    rownames(attachedEsADD) <- namesADD

    for (i in experiments){
        if (is.na(table(attachedEs)[i])){
            tablt=0
        } else tablt=table(attachedEs)[i]
        ## if (is.na(table(attachedEsADD)[i])){
        Egeneset <- rbind(Egeneset, tablt)
        ## } else Egeneset <- rbind(Egeneset,
        ## paste(tablt, "+", table(attachedEsADD)[i], sep=""))
    }
    Egeneset <- rbind(Egeneset, attachedEsADD)

    rownames(Egeneset)= c(experiments, "null")
    colnames(Egeneset)="noE"
    return(Egeneset)
}

#' create topology for a randomly generated pathway topology
#' @param single number of single knockouts
#' @param double number of double knockouts
#' @param force if true the random model will have a sophisticated logical gate
#' @export
#' @examples
#' model <- CreateTopology(3, 1)
#' @return adjacency matrix
CreateTopology <- function(single, double, force = TRUE) {
    extendedModels <- list()
    singleKOs <- LETTERS[1:single]
    experiments <- singleKOs
    doubleKOs <- lapply(1:double, GenerateDoubleKO, singleKOs)
    doubleKOs <- unlist(unique(doubleKOs))
    mutants   <- sort(c(singleKOs, doubleKOs))
    donotextend <- FALSE

    while (length(extendedModels)==0){
        startModel   <- CreateRandomGraph(singleKOs)
        startModel <- startModel[order(apply(startModel, 1, sum),
                                       decreasing = TRUE),
                                 order(apply(startModel, 1, sum),
                                       decreasing = TRUE)]
        startModel[lower.tri(startModel)] <- 0
        startModel <- startModel[order(rownames(startModel)),
                                 order(colnames(startModel))]
        diag(startModel) <- 0
        if (force) {
            if (sum(apply(startModel, 2, sum) >= 2) > 0) {
                if (sum(startModel[,which(
                    apply(startModel, 2, sum) >= 2)[1]] == 1) > 0) {
                    if (!(paste(
                             rownames(
                                 startModel)[which(
                                         startModel[, which(
                                             apply(startModel, 2, sum) >= 2)[1]]
                                         == 1)[1:2]], collapse = ".") %in%
                          mutants)) {
                        mutants <-
                            sort(c(singleKOs,
                                   paste(rownames(
                                       startModel)[which(startModel[,which(
                                                     apply(startModel,
                                                           2,
                                                           sum) >= 2)[1]]
                                                     == 1)[1:2]],
                                       collapse = ".")))
                    }
                    donotextend <- FALSE
                } else {
                    donotextend <- TRUE
                }
            } else {
                donotextend <- TRUE
            }
        }
        if (!donotextend) {
            extendedModels <- includeLogic(startModel, experiments, mutants)
            ## extendedModels <- unlist(extendedModels, recursive=FALSE)
        }
    }
    if (length(extendedModels) == 5) {
        prob <- c(0.018, 0.49, 0.49, 0.001, 0.001)
    } else {
        prob <- rep(1/length(extendedModels), length(extendedModels))
    }

    selectedModel <- sample(1:length(extendedModels), 1, prob = prob)
    topology      <- extendedModels[[selectedModel]]
    
    return(topology)
}

#' generate a random double mutant for artificial data
#' @param singleKOs: vector of single mutants
#' @importFrom gtools combinations
#' @noRd
GenerateDoubleKO <- function(d, singleKOs) {
    allDoubles    <- combinations(length(singleKOs), 2, singleKOs)
    doubleKO      <- allDoubles[d, ]
    doubleKO      <- do.call(paste, as.list(c(doubleKO, sep=".")))
    return(doubleKO)
}
