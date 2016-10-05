## functions needed in order to run LogicNEM

#' Create a random graph
#' @description Returns a model graph with randomly sampled edges. Every possible edge has a probability to exist in the graph.
#' @param pathwayGenes: vector of genes in the pathway
#' @param edgeProb: probability of random edge
#' @export
CreateRandomGraph <- function(pathwayGenes, edgeProb=0.5) {
    size     <- length(pathwayGenes)
    model    <- diag(size)
    nEdges   <- size * (size - 1)
    model[which(model == 0)] <- sample(0:1, nEdges, TRUE, c(1 - edgeProb, edgeProb))
    dimnames(model) <- list(pathwayGenes, pathwayGenes)
    return(model)
}

get_col <- function(i, m) ((i-1) %/% nrow(m)) + 1
get_row <- function(i, m) (i-1) %% nrow(m) + 1

#' Evaluation of graphs
#' @param Phi: model to be evaluated
#' @param D1: observed data matrix
#' @param D0: complementary D1
#' @description Computes marginal log-likelihood for model Phi given observed data matrix D1
#' @export
MLL <- function(Phi, D1, D0) {
                                        # Computes marginal log-likelihood for model Phi, observed data matrix D1, and
                                        # complementary data matrix D0
                                        # Function adapted from NEM package
    Phi2 <- as.matrix(Phi$model)
    Phi2 <- Phi2[colnames(D1), ]
    para <- c(0.13,0.05)
    L    <- para[1]^(D1 %*% (1 - Phi2)) * (1 - para[1])^(D0 %*% (1 - Phi2)) *
                                                            (1 - para[2])^(D1 %*% Phi2) * para[2]^(D0 %*% Phi2)
    posterior <- L / (rowSums(L))
    LLperGene <- log(rowSums(L))
    mLL       <- sum(LLperGene)
    theta     <- apply(posterior, 1, function(e) e == max(e))
    mappos    <- apply(theta, 1, which)
    return(list(posterior=posterior, LLperGene=LLperGene, mLL=mLL, mappos=mappos,
                para=para))
}

#' @noRd
colours <- function(logic, parents2){
                                        #defines colorcode for results
                                        #logic=as.character(logic[2])
    if (logic==NOT2){
        col <- "green3"
        pch <- 21
        log <- NOT2
    }else if (logic=="AND"){
        col <- "red"
        pch <- 22
        log <- "and"
    }else if (logic=="OR"){
        col <- "royalblue1"
        pch <- 23
        log <- "or"
    }else if (logic=="XOR"){
        col <- "orange"
        pch <- 24
        log <- "or"
    }else if (logic==NOT1){
        col <- "purple"
        pch <- 25
        log <- NOT1
    }
    return(list(col=col, pch=pch, log=log))
}

#' get starting vector
#' @description get mutant vectors for simulating boolean network
#' @param mutants: vector of single knockouts
#' @param experiments: vector of all knockouts
#' @export
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
#' @description extend adjacency matrices taking cycles and logics into account. For every given start state, the final state
#' is computed yu using BoolNet.
#' @param network: network created by BoolNet from file
#' @param mutants: vector of single knockouts
#' @param experiments: vector of all knockouts
#' @importFrom BoolNet getPathToAttractor
#' @export
createExtendedAdjacency <- function(network, mutants, experiments){
    starters <- matrix(getStarters(mutants, experiments), length(mutants), byrow=TRUE)
    rownames(starters) <- mutants
    colnames(starters) <- experiments

    ## get steady state of each startState
    getMyAttractor <- function(e, network, includeAttractorStates) {
        network$fixed[which(e == 1)] <- 1
        return(getPathToAttractor(network, e, includeAttractorStates))
    }
    data <- apply(starters, 1, getMyAttractor, network, includeAttractorStates="all")

    ## get observed data by taking final state from all the random startStates.
    extadj <- data[[1]][nrow(data[[1]]),]
    for (i in 2:length(data)){
        d <- data[[i]][nrow(data[[i]]),]
        extadj <- rbind(extadj, d);
    }
    rownames(extadj) <- mutants
    return(extadj)
}

#' include logics into network. Returns an extended adjacency matrix.
#' @param adj: adjacency matrix of graph
#' @param experiments: vector of all knockouts
#' @param mutants: vector of single knockouts
#' @export
includeLogic <- function(adj, experiments, mutants){
    adj=matrix(unlist(adj),length(experiments), byrow=FALSE)
    rownames(adj) <- experiments
    colnames(adj) <- rownames(adj)
    diag(adj)=0
    # adj <- adj[order(apply(adj, 1, sum), decreasing = T), order(apply(adj, 1, sum), decreasing = T)]
    # adj[lower.tri(adj)] <- 0
    adj <- adj[order(rownames(adj)), order(colnames(adj))]
    mutantslist <- strsplit(mutants, ".", fixed=TRUE)
    doublepos <- c()
    for (i in 1:length(mutantslist)){
        if (length(unlist(mutantslist[i])) == 2)
            doublepos <- c(doublepos, i)
    }
    ## look for suitable triples and create logic vector depending
    ## on the relationship of the parents
    notriples <- 0
    logic <- c()
    liste <- list()
    column <- c()
    parents <- experiments[which(experiments %in% unlist(mutantslist[doublepos], recursive=FALSE))]
    NOT2 <- paste(parents[2], "masks the effect of", parents[1])
    NOT1 <- paste(parents[1], "masks the effect of", parents[2])
    for (c in 1:length(experiments)){
        singles <- c()
        parents <- names(which(adj[,c]==1))
        if ((length(parents) == 2) &&
            ##check whether double mutant is abvailable in the data
            (parents[1] %in% unlist(mutantslist[doublepos])) && parents[2] %in% unlist(mutantslist[doublepos])){
            for (i in 1:length(parents)) {
                singles <- cbind(singles, parents[i])
            }
            notriples <- notriples+1
            relation <- adj[singles, singles]
            diag(relation) <- 0
            if (sum(relation) == 0){
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
        ## the next two ifs are from yours truly and yours truly should double check that with florian
        if (length(parents) >= 2 & !(all(parents %in% unlist(mutantslist[doublepos])))) {
            notriples <- notriples + 1
            liste[[notriples]] <- "OR" # because this is normal NEM; seems like a hack that will crash the universe
            column <- cbind(column, c)
        }
        if (length(parents) == 1) {
            notriples <- notriples + 1
            liste[[notriples]] <- "OR" # because this is normal NEM; seems like a hack that will crash the universe
            column <- cbind(column, c)
        }
    }
    if (length(liste) > 0){
        logicmatrix <- as.matrix(expand.grid(liste))
        ## create logics file from adjacency matrix using logics provided by logic vector
        ## ready for using BoolNet
        randomnames <- runif(nrow(logicmatrix)) # why is a vector of 5 not enough and i end up getting a outfile_NA? Because, you dummy, it is a combination of all possible logics for all possible gates.
        for (modelno in 1:nrow(logicmatrix)){
            lo <- 0
            path <- paste("outfile_", randomnames[modelno], ".txt", sep="") # change that !!! how? i don't know, think of something! yea, later. boy, i hope the rest is correct...
            sink(path)
            cat("targets, factors")
            cat("\n")
            for (c in experiments){
                count <- 1
                if (sum(adj[, which(colnames(adj) %in% c)]) == 0) {
                    cat(paste(c, ", ", c, sep=""))
                } else {
                    cat(paste(c, ", ", sep=""))
                }
                count2 <- 0
                count3 <- 0
                for (r in experiments){
                    if (adj[r,c]==1) {
                        if ((count==1) && (which(rownames(adj)==c) %in% column)){
                            help <- r
                            count <- count+1
                            if (sum(adj[, c]) == 1) {
                                cat(paste(r, sep=""))
                            }
                        }
                        else if ((count == 2) && (which(rownames(adj)==c) %in% column)){
                            lo <- lo+1
                            if (logicmatrix[modelno, lo]=="OR" & count3 == 0) { cat(paste(help, " | ", r, sep="")); count3 <- count3 + 1 }
                            else if (logicmatrix[modelno, lo]=="OR" & count3 > 0) cat(paste(" | ", help, " | ", r, sep=""))
                            else if (logicmatrix[modelno, lo]=="AND") cat(paste("(", help, " & ", r, ")", sep=""))
                            else if (logicmatrix[modelno, lo]=="XOR") cat(paste("( ", help, " & ! ", r ,") | (", r, " & ! ", help, ")"))
                            ## help refers to the first element
                            else if (logicmatrix[modelno, lo]==NOT2) cat(paste("(", help, " & ! ", r, ")", sep=""))
                            else if (logicmatrix[modelno, lo]==NOT1) cat(paste("(", r, " & ! ", help, ")", sep=""))
                        }
                        else {
                            count2 <- count2 + 1
                            if (count2 < sum(adj[, which(colnames(adj) %in% c)])) {
                                cat(paste(r, " | ", sep=""))
                            } else {
                                cat(paste(r, sep=""))
                            }
                        }
                    }
                }
                cat("\n")
            }
            sink()
        }
        test <- lapply(1:nrow(logicmatrix), function(x) getExtendedAdjacency(x, logicmatrix, column, adj, mutants, experiments, randomnames))
        return(test)
    }
}

## to do: sehr unschoen!!!
#' create with logics extended adjacency matrix
#' @importFrom BoolNet loadNetwork
#' @export
getExtendedAdjacency <- function(modelno, logicmatrix, column, adj, mutants, experiments, randomnames){
                                        #creates file read by boolNet
    path <- paste("outfile_", randomnames[modelno], ".txt", sep="")
    network <- loadNetwork(path)
    extadj2 <- createExtendedAdjacency(network, unique(mutants), experiments)
    unlink(path)
    return(list(list(origModel=adj, model=extadj2, logics=logicmatrix[modelno,], column=column)))
}

#' @noRd
#' @export
AttachEGenes <- function(posterior, experiments){
    maxpost <- lapply(1:nrow(posterior), function(x) length(which(posterior[x,]==max(posterior[x,]))))
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
                                        # for (e in 1:unlist(maxpost[i])){
            attachedEsADD <- attachedEsADD+1
                                        #   namesADD[k] <- rownames(posterior)[i]
                                        #   k <- k+1
                                        #}
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
                                        # if (is.na(table(attachedEsADD)[i])){
        Egeneset <- rbind(Egeneset, tablt)
                                        #}else Egeneset <- rbind(Egeneset, paste(tablt, "+", table(attachedEsADD)[i], sep=""))
    }
    Egeneset <- rbind(Egeneset, attachedEsADD)

    rownames(Egeneset)= c(experiments, "null")
    colnames(Egeneset)="noE"
    return(Egeneset)
}

                                        # AttachEGenes <- function(posterior, experiments){
                                        #   maxpost <- lapply(1:nrow(posterior), function(x) length(which(posterior[x,]==max(posterior[x,]))))
                                        #   attachedEs <- matrix()
                                        #   names <- c()
                                        #   attachedEsADD <- matrix()
                                        #   namesADD <- c()
                                        #   Egeneset <- matrix()
                                        #   Egeneset <- NULL
                                        #   EgeneADD <- matrix()
                                        #   EgeneADD <- NULL
                                        #   j <- 1
                                        #   k <- 1
                                        #   Epos <- c()
                                        #   for (i in 1:length(maxpost)){
                                        #     if (maxpost[i]==1){
                                        #       Epos <- cbind(Epos, i)
                                        #       attachedEs[j] <- names(which.max(posterior[i,]))
                                        #       names[j] <- rownames(posterior)[i]
                                        #       j <- j+1
                                        #     } else {
                                        #       for (e in 1:unlist(maxpost[i])){
                                        #         attachedEsADD[k] <- names(which(posterior[i,]==max(posterior[i,]))[e])
                                        #         namesADD[k] <- rownames(posterior)[i]
                                        #         k <- k+1
                                        #       }
                                        #     }
                                        #   }
                                        #   attachedEs=as.matrix(attachedEs)
                                        #   rownames(attachedEs) <- names
                                        #   attachedEsADD=as.matrix(attachedEsADD)
                                        #   rownames(attachedEsADD) <- namesADD
                                        #
                                        #   for (i in experiments){
                                        #     if (is.na(table(attachedEs)[i])){
                                        #       tablt=0
                                        #     } else tablt=table(attachedEs)[i]
                                        #     if (is.na(table(attachedEsADD)[i])){
                                        #       Egeneset <- rbind(Egeneset, tablt)
                                        #     }else Egeneset <- rbind(Egeneset, paste(tablt, "+", table(attachedEsADD)[i], sep=""))
                                        #   }
                                        #   rownames(Egeneset)=experiments
                                        #   colnames(Egeneset)="noE"
                                        #   return(Egeneset)
                                        # }

#' create topology for a randomly generated pathway topology
#' @param single: number of single knockouts
#' @param double: number of double knockouts
#' @export
CreateTopology <- function(single, double, force = T) {
    extendedModels <- list()
    singleKOs <- LETTERS[1:single]
    experiments <- singleKOs
    doubleKOs <- lapply(1:double, function(d) GenerateDoubleKO(singleKOs))
    doubleKOs <- unlist(unique(doubleKOs))
    mutants   <- sort(c(singleKOs, doubleKOs))
    donotextend <- FALSE

    while (length(extendedModels)==0){
        startModel   <- CreateRandomGraph(singleKOs)
        startModel <- startModel[order(apply(startModel, 1, sum), decreasing = T), order(apply(startModel, 1, sum), decreasing = T)]
        startModel[lower.tri(startModel)] <- 0
        startModel <- startModel[order(rownames(startModel)), order(colnames(startModel))]
        diag(startModel) <- 0
        if (force) {
            if (sum(apply(startModel, 2, sum) >= 2) > 0) {
                if (sum(startModel[, which(apply(startModel, 2, sum) >= 2)[1]] == 1) > 0) {
                    mutants <- sort(c(singleKOs, paste(rownames(startModel)[which(startModel[, which(apply(startModel, 2, sum) >= 2)[1]] == 1)[1:2]], collapse = ".")))
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
        prob <- c(0.018, 0.49, 0.49, 0.001, 0.001) # or = [1] is normal nem so has a higher prob for networks > 3, and the nots = [4:5] have sure probs if parents are related
    } else {
        prob <- rep(1/length(extendedModels), length(extendedModels)) # this does not make sense rbecause length should be equal to 1, right?
    }

    selectedModel <- sample(1:length(extendedModels), 1, prob = prob)
    topology      <- extendedModels[[selectedModel]]
    topology
    
    return(topology)
}

#' generate a random double mutant for artificial data
#' @param singleKOs: vector of single mutants
#' @importFrom gtools combinations
#' @export
GenerateDoubleKO <- function(singleKOs) {
    allDoubles    <- combinations(length(singleKOs), 2, singleKOs)
    randomDoubles <- sample(1:choose(length(singleKOs), 2), 1)
    doubleKO      <- allDoubles[randomDoubles, ]
    doubleKO      <- do.call(paste, as.list(c(doubleKO, sep=".")))
    return(doubleKO)
}

#' @noRd
#' @export
epiNEM.Simulations <- function(random, nSim){
                                        #perform simulation experiments for LogicNEM, Correlation and Mutual Information
    count <- rep(0, length(random$FNrate))
    log <- rep(0, length(random$FNrate))
    errorLOG <- list(0)
    errorCOR1 <- list(0)
    errorCOR2 <- list(0)
    errorCOR3 <- list(0)
    errorMUT1 <- list(0)
    errorMUT2 <- list(0)
    errorNEM <- list(0)
    good <- c()
    avgErrorLOG <- list()
    avgErrorMUT1 <- list()
    avgErrorMUT2 <- list()
    avgErrorCOR1 <- list()
    avgErrorCOR2 <- list()
    avgErrorCOR3 <- list()
    avgErrorNEM <- list()
    for (j in 1:length(random$FNrate)){
        for(i in 1:nSim){
            topology <- CreateTopology(random$single, random$double)
            topology <- unlist(unique(topology), recursive=FALSE)
            extTopology <- ExtendTopology(topology$model, random$reporters)
            sortedData  <- GenerateData(topology$model, extTopology, random$FPrate,
                                        random$FNrate[j], random$replicates)
            ##LogicNEM
                                        #       scoreTrip1 <- sapply(tripl12, MLL, sortedData, 1-sortedData)
                                        #       mll=unlist(scoreTrip1["mLL",])
                                        #       index <- which.max(mll)
                                        #       TriplModel <- tripl12[index]
            TriplModel <- LogicNEM(sortedData, method = "exhaustive")
                                        #       if (all.equal(TriplModel[[1]], topology)==TRUE){
                                        #         count[j]=count[j]+1
                                        #         good[count]=TriplModel[[1]]$logics
                                        #         errorLOG[[1]][i]=0
                                        #         log[j]=log[j]+1
                                        #      } else {
            if (TriplModel$logics==topology$logics){
                log[j]=log[j]+1
            }
            errorLOG[[1]][i]=sum(1-(as.numeric(topology$origModel==TriplModel$origModel)))
                                        #}
            ##Correlation
            cori <- cor(sortedData[,1:3])
            diag(cori) <- 0
            cori1 <- cori
            cori2 <- cori
            cori3 <- cori
            cori1[cori1>abs(0.2)] <- 1
            cori2[cori2>abs(0.5)] <- 1
            cori3[cori3>abs(0.8)] <- 1
            errorCOR1[[1]][i] <- sum(1-(as.numeric(topology$origModel==cori1)))/2
            errorCOR2[[1]][i] <- sum(1-(as.numeric(topology$origModel==cori2)))/2
            errorCOR3[[1]][i] <- sum(1-(as.numeric(topology$origModel==cori3)))/2
            ## Mutual Information
            mut <- minet(sortedData[,1:3], method="mrnet", estimator="spearman", disc="none", nbins=sqrt(NROW(dataset)))
            mut1 <- mut
            mut2 <- mut
            mut1[mut1>0.5] <- 1
            mut2[mut2>0.8] <- 1
            errorMUT1[[1]][i] <- sum(1-(as.numeric(topology$origModel==mut1)))/2
            errorMUT2[[1]][i]=sum(1-(as.numeric(topology$origModel==mut2)))/2

                                        #classic NEM
            nemi <- nem(sortedData, inference = "search")
            nemi <- igraph.from.graphNEL(nemi$graph)
            nemi <- as.matrix(get.adjacency(nemi))
            nemi <- nemi[which(rownames(nemi) %in% c("A", "B", "C")),which(colnames(nemi) %in% c("A", "B", "C"))]
            errorNEM[[1]][i] <- sum(1-(as.numeric(topology$origModel==nemi)))

        }
        avgErrorLOG <- c(avgErrorLOG, errorLOG)
        avgErrorCOR1 <- c(avgErrorCOR1, errorCOR1)
        avgErrorCOR2 <- c(avgErrorCOR2, errorCOR2)
        avgErrorCOR3 <- c(avgErrorCOR3, errorCOR3)
        avgErrorMUT1 <- c(avgErrorMUT1, errorMUT1)
        avgErrorMUT2 <- c(avgErrorMUT2, errorMUT2)
        avgErrorNEM <- c(avgErrorNEM, errorNEM)
    }
    return(list(count = count, errorCOR02 = avgErrorCOR1, errorCOR05 = avgErrorCOR2, errorCOR08 = avgErrorCOR3,
                errorMUT05 = avgErrorMUT1, errorMUT08 = avgErrorMUT2, errorNEM = avgErrorNEM, errorLOG=avgErrorLOG, logs=log))
                                        #, LogicNEM=avgErrorLOG, MutualInfo=avgErrorMUT, Correlation=avgErrorCOR, alpha=alpha))#, log=log))#, modTip=modTip, modTop=modTop))
}

#' @noRd
#' @export
getGeneName <- function(symbol){
    name <- as.character(unlist(xx[which(xx==symbol)], recursive=FALSE))
    if (length(name)==0){
        return(symbol)
    }else return(name)
}
