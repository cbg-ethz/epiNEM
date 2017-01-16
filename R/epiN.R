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
epiNEM <- function(filename="random", method="greedy", nIterations=10, nModels=0, random=list(single=4, double=1, reporters=100, FPrate=0.1, FNrate=0.1, replicates=1), plotsy=TRUE, ltype = "marginal", para = c(0.13, 0.05)) {

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
            sortedData  <- GenerateData(topology$model, extTopology, random$FPrate,
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
        iterations <- sapply(1:nIterations, GreedyHillClimber, experiments, D1, D0, mutants)

        index <- which.max(iterations["maxLlh",])
        reconstruct <- includeLogic(iterations[,index]$model, experiments, mutants)
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
        extendedModels <- lapply(basicModels, includeLogic, experiments, unique(mutants))
        extendedModels <- unlist(unlist(extendedModels, recursive=FALSE) , recursive = FALSE)
        uniqueModels   <- unique(lapply(extendedModels, function(e)
            identity((e["model"]))))
        uniqueModels <- uniqueModels#[1:length(uniqueModels)-1] # what is this?

        if (length(doubleKOs) == 0) {
            print("No double perturbations available --> computing NEM")
            if (method == "exhaustive") { inference <- "search" }
            if (method == "greedy") { inference <- "nem.greedy" }
            options <- set.default.parameters(setdiff(unique(colnames(D)),"time"))
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
        results <- lapply(results, utils::modifyList, list(score=mll[which.max(mll)]))
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
    reporters     <- unlist(lapply(1:nReporters, function(n) paste("reporter", n, sep="-")))
    linkedEffects <- sample(1:ncol(topology), nReporters, replace=TRUE)
    extTopology   <- sapply(linkedEffects, function(e) lapply(1:ncol(topology), function(c) ifelse(e == c, 1, 0)))
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

#' Convert epiNEM model into general Boolean graph
#' @param t full epiNEM model
#' @author Madeline Diekmann
#' @seealso CreateTopology
#' @export
#' @examples
#' topology <- CreateTopology(3, 1, force = TRUE)
#' topology <- unlist(unique(topology), recursive = FALSE)
#' extTopology <- ExtendTopology(topology$model, 100)
#' b <- EpiNEM2BooleanGraph(extTopology)
#' @return boolean hyper-graph
EpiNEM2BooleanGraph <- function(t) {
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
            if (length(parents) == 2) {
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
                    graph <- gsub("\\+=", "=", graph)
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
            if (length(parents) > 2 | length(parents) == 1) {
                graph <- c(graph, paste(sort(parents), child, sep = "="))
            }
        }
        all <- rownames(t$origModel)
        children2 <- unique(gsub(".*=", "", graph))
        if (sum(!(all %in% children2)) > 0) {
            graph <- c(graph, paste("S", all[which(!(all %in% children2))], sep = "="))
        }
        return(unique(graph))
    }
}

#' @noRd
convertGraph <- function(g) {
  g <- sort(g)
  targets <- gsub(".*=", "", g)
  g.new <- NULL
  for (i in unique(targets)) {
    dnf <- list()
    count <- 1
    for (j in g[grep(paste("=", i, sep = ""), g)]) {
      dnf[[count]] <- sort(unique(unlist(strsplit(gsub("=.*", "", j), "\\+"))))
      count <- count + 1
    }
    cnf <- expand.grid(dnf)
    dnf <- NULL
    for (j in 1:dim(cnf)[1]) {
      dnf <- c(dnf, paste(sort(unique(unlist(cnf[j, ]))), collapse = "+"))
    }
    dnf <- paste(sort(dnf), "=", i, sep = "")
    g.new <- c(g.new, dnf)
  }
  vertices <- sort(unique(unlist(strsplit(unlist(strsplit(g.new, "=")), "\\+"))))
  for (i in vertices) {
    if (length(grep(paste(i, ".*", i, ".*=", sep = ""), g.new)) > 0) {
      g.new <- g.new[-grep(paste(i, ".*", i, ".*=", sep = ""), g.new)]
    }
  }
  return(g.new)
}

#' @noRd
dnf2adj <- function(dnf, closed = FALSE) {
  if (length(dnf) == 0) {
    return(NULL)
  } else {
    nodes <- character()
    for (i in dnf) {
      tmp <- unlist(strsplit(i, "="))
      nodes <- c(nodes, tmp[2], unlist(strsplit(tmp[1], "\\+")))
    }
    nodes <- unique(gsub("!", "", nodes))
    adjmat <- matrix(0, length(nodes), length(nodes))
    colnames(adjmat) <- nodes
    rownames(adjmat) <- nodes
    for (i in dnf) {
      tmp <- unlist(strsplit(i, "="))
      child <- tmp[2]
      parents <- unlist(strsplit(tmp[1], "\\+"))
      for (j in parents) {
        if (gsub("!", "", j) %in% j) {
          adjmat[which(rownames(adjmat) %in% j), which(colnames(adjmat) %in% child)] <- 1
        } else {
          adjmat[which(rownames(adjmat) %in% gsub("!", "", j)), which(colnames(adjmat) %in% child)] <- -1
        }
      }
    }
    diag(adjmat) <- 1
    if (closed) {
      stop <- FALSE
      cons <- c(TRUE, rep(FALSE, (length(adjmat) - 1)))
      while(!stop) {
        adjmat <- adjmat%*%adjmat
        if (all(cons == (adjmat != 0))) {
          stop <- TRUE
        } else {
          cons <- (adjmat != 0)
        }
      }
      adjmat[adjmat > 1] <- 1
      adjmat[adjmat < -1] <- -1
    }
    return(adjmat)
  }
}

#' @noRd
getHierarchy <- function(graph) {
  adj <- dnf2adj(graph)
  dnf <- adj2dnf(adj)
  g <- plotDnf(dnf, draw = FALSE)
  Ypos <- g@renderInfo@nodes$labelY
  Ynames <- names(g@renderInfo@nodes$labelY)
  ## Ypos <- Ypos[-grep("and", Ynames)]
  ## Ynames <- Ynames[-grep("and", Ynames)]
  hierarchy <- list()
  count <- 0
  for (i in sort(unique(Ypos), decreasing = TRUE)) {
    count <- count + 1
    hierarchy[[count]] <- Ynames[which(Ypos == i)]
  }
  return(hierarchy)
}

#' @noRd
transClose <- function(g, max.iter = NULL, verbose = FALSE) {
  v <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(g, "=")), "\\+"))))
  if (is.null(max.iter)) {
    h <- getHierarchy(g)
    max.iter <- length(h) - 2 # should be sufficient
  }
  if (verbose) {
    print(paste("maximum iterations: ", max.iter, sep = ""))
  }
  g.out <- unique(gsub(".*=", "", g))
  g.closed <- g
  for (iter in 1:max.iter) {
    g.old <- g.closed
    
    if (verbose) {
      cat('\r', paste("iteration: ", iter, sep = ""))
      flush.console()
    }
    for (i in g.closed) {
      input <- gsub("!", "", unlist(strsplit(unlist(strsplit(i, "="))[1], "\\+")))
      input <- intersect(input, g.out)
      output <- gsub(".*=", "", i)
      if (length(input) == 0) { next() }
      for (j in unique(input)) {
        if (j %in% unlist(strsplit(unlist(strsplit(i, "="))[1], "\\+"))) {
          for (k in gsub("=.*", "", g[grep(paste("=", j, sep = ""), g)])) {
            g.closed <- c(g.closed, gsub(j, k, i))
          }
        } else {
          literals <- list()
          count <- 0
          for (k in unique(gsub("=.*", "", g[grep(paste("=", j, sep = ""), g)]))) {
            count <- count + 1
            literals[[count]] <- gsub("!!", "", paste("!", unlist(strsplit(k, "\\+")), sep = ""))
          }
          combis <- expand.grid(literals)
          combis <- apply(combis, c(1,2), as.character)
          for (k in 1:nrow(combis)) {
            g.closed <- c(g.closed, gsub(paste("!", j, sep = ""), paste(combis[k, ], collapse = "+"), i))
          }
        }
      }
    }
    g.closed <- unique(g.closed)
    if (all(g.closed %in% g.old)) {
      if (verbose) {
        cat("\n")
        print(paste("successfull convergence", sep = ""))
      }
      break()
    }
  }
  if (verbose) {
    cat("\n")
  }
  g.closed <- unique(g.closed)
  for (j in 1:length(g.closed)) {
    i <- g.closed[j]
    input <- unlist(strsplit(i, "="))
    output <- input[2]
    input <- unlist(strsplit(input[1], "\\+"))
    input <- unique(input)
    g.closed[j] <- paste(paste(input, collapse = "+"), output, sep = "=")
  }
  if (length(grep(paste(paste(v, ".*", v, ".*=", sep = ""), collapse = "|"), g.closed)) > 0) {
    g.closed <- g.closed[-grep(paste(paste(v, ".*", v, ".*=", sep = ""), collapse = "|"), g.closed)]
  }
  return(g.closed)
}

#' @noRd
transRed <- function(g, max.iter = NULL, verbose = FALSE) {
  v <- unique(gsub("!", "", unlist(strsplit(unlist(strsplit(g, "=")), "\\+"))))
  if (is.null(max.iter)) {
    max.iter <- length(v) - 2
  }
  a <- dnf2adj(g)
  g2 <- g
  h <- getHierarchy(g2)
  if (length(h) > 2) {
    for (i in 1:(length(h)-2)) {
      for (j in h[[i]]) {
        for (k in (i+2):length(h)) {
          for (l in h[[k]]) {
            if (length(grep(paste(".*", j, ".*=", l, sep = ""), g2)) != 0) {
              if (length(grep(paste(".*=", l, sep = ""), g2)) > 1) {
                g2 <- g2[-grep(paste(".*", j, ".*=", l, sep = ""), g2)]
              }
            }
          }
        }
      }
    }
  }
  g3 <- transClose(g2, max.iter)
  if (sum(g %in% g3) > 0) {
    g4 <- g[-which(g %in% g3)]
  }
  g5 <- unique(c(g2, g4))
  return(g5)
}

#' @noRd
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

#' @noRd
#' @import Rgraphviz
plotDnf <- function(dnf = NULL, freq = NULL, stimuli = c(), signals = c(), inhibitors = c(), connected = TRUE,  CNOlist = NULL, cex = NULL, fontsize = NULL, labelsize = NULL, type = 2, lwd = 2, edgelwd = 2, legend = 0, x = 0, y = 0, xjust = 0, yjust = 0, width = 1.5, height = 1, rankdir = "TB", rank = "same", layout = "dot", main = "", sub = "", cex.main = 1.5, cex.sub = 1, col.sub = "grey", fontcolor = NULL, nodestates = NULL, simulate = NULL, andcolor = "transparent", edgecol = NULL, labels = NULL, labelcol = "blue", nodelabel = NULL, nodecol = NULL, bordercol = NULL, nodeshape = NULL, verbose = FALSE, edgestyle = NULL, nodeheight = NULL, nodewidth = NULL, edgewidth = NULL, lty = NULL, hierarchy = NULL, showall = FALSE, nodefontsize = NULL, edgehead = NULL, edgelabel = NULL, edgetail = NULL, bool = TRUE, draw = TRUE, ...) {
    ## see graphvizCapabilities()$layoutTypes for supported layouts

    if (!bool & length(grep("\\+", dnf)) > 0) {
        dnf <- dnf[-grep("\\+", dnf)]
    }

    graph <- dnf

    if (!is.null(hierarchy)) {
        if (!showall) {
            nodes <- unique(unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+")))
            for (i in 1:length(hierarchy)) {
                hierarchy[[i]] <- intersect(hierarchy[[i]], nodes)
            }
        }
        graph2 <- NULL
        for (i in graph) {
            input <- unlist(strsplit(i, "="))
            output <- input[2]
            input <- gsub("!", "", unlist(strsplit(input[1], "\\+")))
            for (j in input) {
                graph2 <- c(graph2, paste(j, output, sep = "="))
            }
            graph2 <- unique(graph2)
        }
        hgraph <- NULL
        hgraph2 <- NULL
        for (i in 1:(length(hierarchy)-1)) {
            for (j in hierarchy[[i]]) {
                for (k in hierarchy[[i+1]]) {
                    hgraph <- c(hgraph, paste(j, k, sep = "="))
                    hgraph2 <- c(hgraph2, paste(k, j, sep = "="))
                }
            }
        }
        if (sum(hgraph %in% graph2) > 0) {
            hgraph <- hgraph[-which(hgraph %in% graph2)]
        }
        if (sum(hgraph2 %in% graph2) > 0) {
            hgraph <- hgraph[-which(hgraph2 %in% graph2)]
        }
        dnf <- c(graph, hgraph)
        ## update all the parameters e.g. edgestyle...
        if (is.null(edgecol)) {
            edgecol <- c(rep("black", length(grep("\\+", graph))+length(unlist(strsplit(graph, "\\+")))), rep("transparent", length(hgraph)))
            dnf2 <- dnf
            if (length(grep("\\+", dnf2)) > 0) {
                dnf2[-grep("\\+", dnf2)] <- gsub("=", "", dnf2[-grep("\\+", dnf2)])
            } else {
                dnf2 <- gsub("=", "", dnf2)
            }
            edgecol[grep("!", unlist(strsplit(unlist(strsplit(dnf2, "\\+")), "=")))] <- "red"
        } else {
            if (length(edgecol) == 1) {
                edgecol <- c(rep(edgecol, length(grep("\\+", graph))+length(unlist(strsplit(graph, "\\+")))), rep("transparent", length(hgraph)))
            } else {
                edgecol <- c(edgecol, rep("transparent", length(hgraph)))
            }
        }
    } else {
        if (is.null(lty) & !is.null(dnf)) {
            lty <- c(rep("solid", length(grep("\\+", graph))+length(unlist(strsplit(graph, "\\+")))))
        }
    }

    graph <- dnf

    dolegend <- FALSE
    if (is.null(dnf)) {
        dnf <- c("A=B")
        dolegend <- TRUE
    }
    
    if (!is.null(simulate)) {
        nodestates <- simulateDnf(graph, stimuli = simulate$stimuli, inhibitors = simulate$inhibitors)
    }
    
    if (is.null(freq)) {
        use.freq = FALSE
    } else {
        use.freq = TRUE
        if (is.null(labels)) {
            labels <- as.character(round(freq, 2)*100)
        }
    }
    if (is.null(labels)) {
        labels <- rep("", length(dnf))
    }
    
    if (is.null(fontsize)) {
        fontsize <- ""
    }
    if (is.null(labelsize)) {
        labelsize <- fontsize
    }

    if (!is.null(CNOlist)) {
        if (length(stimuli) == 0) {
            stimuli <- colnames(CNOlist@stimuli)
        }
        if (length(signals) == 0) {
            signals <- colnames(CNOlist@signals[[1]])
        }
        if(length(inhibitors) == 0) {
            inhibitors <- colnames(CNOlist@inhibitors)
        }
    }

    if (connected) {
        Vneg <- unique(c(c(unique(unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+"))))))
        if (length(grep("\\+", dnf)) > 0) {
            Vneg <- c(Vneg, paste("and", 1:length(grep("\\+", dnf)), sep = ""))
        }
        V <- unique(gsub("!", "", Vneg))
        stimuli <- intersect(stimuli, V)
        signals <- intersect(signals, V)
        inhibitors <- intersect(inhibitors, V)
    } else {
        Vneg <- unique(c(c(unique(unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+")))), stimuli, signals, inhibitors))
        if (length(grep("\\+", dnf)) > 0) {
            Vneg <- c(Vneg, paste("and", 1:length(grep("\\+", dnf)), sep = ""))
        }
        V <- unique(gsub("!", "", Vneg))
    }

    V <- sort(V)
    
    Vneg <- c(V, Vneg[grep("!", Vneg)])

    if (!is.null(nodecol)) {
        if (length(nodecol) == 1 & !is.list(nodecol)) {
            nodecol.tmp <- nodecol
            nodecol <- list()
            for (i in V) {
                if (length(grep("and", i)) == 0) {
                    nodecol[[i]] <- nodecol.tmp
                }
            }
        }
    }

    if (!is.null(bordercol)) {
        if (length(bordercol) == 1 & !is.list(bordercol)) {
            bordercol.tmp <- bordercol
            bordercol <- list()
            for (i in V) {
                if (length(grep("and", i)) == 0) {
                    bordercol[[i]] <- bordercol.tmp
                }
            }
        }
    }

    E <- list()

    for (i in V) {
        E[[i]] <- list()
    }

    Eneg <- list()

    for (i in Vneg) {
        Eneg[[i]] <- list()
    }

    count <- 0
    
    for (i in dnf) {
        tmp <- unlist(strsplit(i, "="))
        if (length(tmp)==1) {
            Eneg[[tmp]][["edges"]] <- c(Eneg[[tmp]][["edges"]], NULL)
            tmp <- gsub("!", "", tmp)
            E[[tmp]][["edges"]] <- c(E[[tmp]][["edges"]], NULL)
        } else {
            tmp2 <- unlist(strsplit(tmp[1], "\\+"))
            if (length(tmp2) > 1) {
                count <- count + 1
                Eneg[[paste("and", count, sep = "")]][["edges"]] <- c(Eneg[[paste("and", count, sep = "")]][["edges"]], which(Vneg %in% tmp[2]))
                E[[paste("and", count, sep = "")]][["edges"]] <- c(E[[paste("and", count, sep = "")]][["edges"]], which(V %in% tmp[2]))
                for (j in tmp2) {
                    Eneg[[j]][["edges"]] <- c(Eneg[[j]][["edges"]], which(Vneg %in% paste("and", count, sep = "")))
                    j <- gsub("!", "", j)
                    E[[j]][["edges"]] <- c(E[[j]][["edges"]], which(V %in% paste("and", count, sep = "")))
                }
            } else {
                Eneg[[tmp2]][["edges"]] <- c(Eneg[[tmp2]][["edges"]], which(Vneg %in% tmp[2]))
                tmp2 <- gsub("!", "", tmp2)
                E[[tmp2]][["edges"]] <- c(E[[tmp2]][["edges"]], which(V %in% tmp[2]))
            }
        }
    }

    g <- new("graphNEL",nodes=V,edgeL=E,edgemode="directed")

    gneg <- new("graphNEL",nodes=Vneg,edgeL=Eneg,edgemode="directed")

    nodes <- buildNodeList(g)

    edges <- buildEdgeList(g)
    
    nodesneg <- buildNodeList(gneg)

    edgesneg <- buildEdgeList(gneg)

    edgesnew <- list()
    
    for (i in sort(names(edges))) {
        edgesnew <- c(edgesnew, edges[[i]])
    }

    names(edgesnew) <- sort(names(edges))

    edges <- edgesnew

    edgesnegnew <- list()
    
    for (i in names(edgesneg)) {
        edgesnegnew <- c(edgesnegnew, edgesneg[[i]])
    }

    names(edgesnegnew) <- names(edgesneg)

    edgesneg <- edgesnegnew

    if (verbose) {
        print(paste("order of nodes: ", paste(names(nodes), collapse = ", "), sep = ""))
        print(paste("order of edges: ", paste(names(edges), collapse = ", "), sep = ""))
    }

    edges <- edgesneg
    names(edges) <- gsub("!", "", names(edges))

    for (i in 1:length(edges)) {
        edges[[i]]@from <- gsub("!", "", edges[[i]]@from)
    }
    
    nodeshape2 <- nodeshape
    nodeshape <- list()
    if (length(nodeshape2) == 1 & !(is.list(nodeshape2))) {
        for (i in 1:length(nodes)) {
            nodeshape[[nodes[[i]]@name]] <- nodeshape2
        }
    } else {
        nodeshape <- nodeshape2
    }
    
    for (i in 1:length(nodes)) {
        nodes[[i]]@attrs$height <- height
        nodes[[i]]@attrs$width <- width
        if (!is.null(nodelabel[[nodes[[i]]@name]])) {
            nodes[[i]]@attrs$name <- nodelabel[[nodes[[i]]@name]]
        }
        if (!is.null(nodeheight[[nodes[[i]]@name]])) {
            nodes[[i]]@attrs$height <- nodeheight[[nodes[[i]]@name]]
        }
        if (!is.null(nodewidth[[nodes[[i]]@name]])) {
            nodes[[i]]@attrs$width <- nodewidth[[nodes[[i]]@name]]
        }
        if (length(grep("and", nodes[[i]]@name)) > 0) {
            if (is.null(nodelabel)) {
                nodelabel <- list()
            }
            nodelabel[[nodes[[i]]@name]] <- "AND"
            nodes[[i]]@attrs$label <- ""
            nodes[[i]]@attrs$fontcolor <- andcolor
            if (is.null(bordercol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$color <- "grey"
            } else {
                nodes[[i]]@attrs$color <- bordercol[[nodes[[i]]@name]]
            }
            if (is.null(nodeshape[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$shape <- "box"
            } else {
                nodes[[i]]@attrs$shape <- nodeshape[[nodes[[i]]@name]]
            }
            if (is.null(nodecol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$fillcolor <- "grey"
            } else {
                nodes[[i]]@attrs$fillcolor <- nodecol[[nodes[[i]]@name]]
            }
            nodes[[i]]@attrs$width <- "0.5"
            nodes[[i]]@attrs$height <- "0.5"
            if (type == 2) {
                nodes[[i]]@attrs$fontsize <- "0"
            } else {
                nodes[[i]]@attrs$fontsize <- "0"
            }
        } else {
            nodes[[i]]@attrs$fontsize <- as.character(fontsize)
            if (is.null(nodecol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$fillcolor <- "white"
            } else {
                nodes[[i]]@attrs$fillcolor <- nodecol[[nodes[[i]]@name]]
            }
            if (is.null(nodeshape[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$shape <- "ellipse"
            } else {
                nodes[[i]]@attrs$shape <- nodeshape[[nodes[[i]]@name]]
            }
            if (is.null(bordercol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$color <- "black"
            } else {
                nodes[[i]]@attrs$color <- bordercol[[nodes[[i]]@name]]
            }
            if (names(nodes)[i] %in% stimuli & is.null(nodeshape[[nodes[[i]]@name]])) {
                if (type == 2) {
                    nodes[[i]]@attrs$shape <- "diamond"
                } else {
                    nodes[[i]]@attrs$shape <- "box"
                }
            }
            if (names(nodes)[i] %in% signals & is.null(nodecol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$fillcolor <- "lightblue"
            }
            if (names(nodes)[i] %in% inhibitors & is.null(bordercol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$color <- "red"
            }
        }
        if (!is.null(nodestates)) {
            if (sum(names(nodestates) %in% nodes[[i]]@name) == 1) {
                if (nodestates[which(names(nodestates) %in% nodes[[i]]@name)] == 0) {
                    if (is.null(nodecol[[nodes[[i]]@name]])) {
                        nodes[[i]]@attrs$fillcolor <- "white"
                    }
                    if (is.null(bordercol[[nodes[[i]]@name]])) {
                        nodes[[i]]@attrs$color <- "black"
                    }
                }
                if (nodestates[which(names(nodestates) %in% nodes[[i]]@name)] == 1) {
                    if (is.null(nodecol[[nodes[[i]]@name]])) {
                        nodes[[i]]@attrs$fillcolor <- "green"
                    }
                    if (is.null(bordercol[[nodes[[i]]@name]])) {
                        nodes[[i]]@attrs$color <- "black"
                    }
                }
            }
        }
    }
    if (length(edges) > 0) {
        for (i in 1:length(edges)) {
            edges[[i]]@attrs$fontsize <- as.character(labelsize)
            if (length(grep("and", names(edges)[i])) > 0) {
                tmp <- unlist(strsplit(names(edges)[i], "~"))
                k <- grep("and", tmp)
                inputN <- length(grep(tmp[k], edges))
                k <- as.numeric(gsub("and", "", tmp[k]))
                ## try to get the index of the correct edgecol:
                k2 <- grep("\\+", dnf)[k]
                if (grep("and", tmp) == 2) {
                    inputN2 <- which(gsub("!", "", unlist(strsplit(gsub("=.*", "", dnf[k2]), "\\+"))) %in% tmp[1])
                } else {
                    inputN2 <- length(unlist(strsplit(gsub("=.*", "", dnf[k2]), "\\+"))) + 1
                }
                if (k2 == 1) {
                    edgecolindex <- inputN2
                } else {
                    if (length(grep("\\+", graph[1:(k2-1)])) == 0) {
                        edgecolindex <- length(graph[1:(k2-1)]) + inputN2
                    } else {
                        edgecolindex <- length(unlist(strsplit(dnf[1:(k2-1)], "\\+"))) + length(grep("\\+", dnf[1:(k2-1)])) + inputN2
                    }
                }
                ## end
                inputN2 <- grep(tmp[1], unlist(strsplit(dnf[grep("\\+", dnf)[k]], "\\+")))-1
                edges[[i]]@attrs$style <- lty[grep("\\+", dnf)[k]]
                edges[[i]]@attrs$label <- labels[grep("\\+", dnf)[k]]
                if (use.freq) {
                    edges[[i]]@attrs$weight <- freq[grep("\\+", dnf)[k]]
                    edges[[i]]@attrs$fontcolor <- "blue"
                }
                if (!is.null(edgewidth)) {
                    edges[[i]]@attrs$weight <- edgewidth[grep("\\+", dnf)[k]]
                }
                if (!is.null(edgestyle)) {
                    if (!is.na(edgestyle[grep("\\+", dnf)[k]])) {
                        edges[[i]]@attrs$style <- edgestyle[edgecolindex]
                    }
                }
                if (!is.null(edgelabel)) {
                    if (!is.na(edgelabel[grep("\\+", dnf)[k]])) {
                        edges[[i]]@attrs$label <- edgelabel[edgecolindex]
                    }
                }
                if (length(grep("!", names(edgesneg)[i])) > 0) {
                    edges[[i]]@attrs$arrowhead <- "tee"
                    edges[[i]]@attrs$color <- "red"
                    if (!is.null(edgecol)) {
                        if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
                            edges[[i]]@attrs$color <- edgecol[edgecolindex]
                        }
                    }
                    if (!is.null(edgehead)) {
                        if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
                            edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex]
                        }
                    }
                    if (!is.null(edgetail)) {
                        if (!is.na(edgetail[grep("\\+", dnf)[k]])) {
                            edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex]
                        }
                    }
                } else {
                    edges[[i]]@attrs$arrowhead <- "open"
                    edges[[i]]@attrs$color <- "black"
                    if (gsub("and.*", "and", tmp[1]) %in% "and") {
                        if (!is.null(edgecol)) {
                            if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$color <- edgecol[edgecolindex]
                            }
                        }
                        if (!is.null(edgehead)) {
                            if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex]
                            }
                        }
                        if (!is.null(edgetail)) {
                            if (!is.na(edgetail[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex]
                            }
                        }
                    } else {
                        if (!is.null(edgecol)) {
                            if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$color <- edgecol[edgecolindex]
                            }
                        }
                        if (!is.null(edgehead)) {
                            if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex]
                            }
                        }
                        if (!is.null(edgetail)) {
                            if (!is.na(edgetail[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex]
                            }
                        }
                    }
                }
            } else {
                tmp <- unlist(strsplit(names(edges)[i], "~"))
                ## try to get the index of the correct edgecol:
                if (length(grep("!", names(edgesneg)[i])) == 0) {
                    k2 <- grep(paste("^", tmp[1], "=", tmp[2], sep = ""), dnf)
                } else {
                    k2 <- grep(paste("^!", tmp[1], "=", tmp[2], sep = ""), dnf)
                }
                if (k2 == 1) {
                    edgecolindex <- k2
                } else {
                    if (length(grep("\\+", dnf[1:(k2-1)])) == 0) {
                        edgecolindex <- k2
                    } else {
                        edgecolindex <- length(unlist(strsplit(dnf[1:(k2-1)], "\\+"))) + length(grep("\\+", dnf[1:(k2-1)])) + 1
                    }
                }
                ## end
                if (length(grep("!", names(edgesneg)[i])) > 0) {
                    edges[[i]]@attrs$style <- lty[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    edges[[i]]@attrs$label <- labels[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    if (use.freq) {
                        edges[[i]]@attrs$weight <- freq[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                        edges[[i]]@attrs$fontcolor <- "blue"
                    }
                    if (!is.null(edgewidth)) {
                        edges[[i]]@attrs$weight <- edgewidth[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    }
                    if (!is.null(edgestyle)) {
                        if (!is.na(edgestyle[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$style <- edgestyle[edgecolindex]
                        }
                    }
                    if (!is.null(edgelabel)) {
                        if (!is.na(edgelabel[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$label <- edgelabel[edgecolindex]
                        }
                    }
                    edges[[i]]@attrs$arrowhead <- "tee"
                    edges[[i]]@attrs$color <- "red"
                    if (!is.null(edgecol)) {
                        if (!is.na(edgecol[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$color <- edgecol[edgecolindex]
                        }
                    }
                    if (!is.null(edgehead)) {
                        if (!is.na(edgehead[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex]
                        }
                    }
                    if (!is.null(edgetail)) {
                        if (!is.na(edgetail[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex]
                        }
                    }
                } else {
                    edges[[i]]@attrs$style <- lty[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    edges[[i]]@attrs$label <- labels[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    if (use.freq) {
                        edges[[i]]@attrs$weight <- freq[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                        edges[[i]]@attrs$fontcolor <- "blue"
                    }
                    if (!is.null(edgewidth)) {
                        edges[[i]]@attrs$weight <- edgewidth[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    }
                    if (!is.null(edgestyle)) {
                        if (!is.na(edgestyle[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$style <- edgestyle[edgecolindex]
                        }
                    }
                    if (!is.null(edgelabel)) {
                        if (!is.na(edgelabel[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$label <- edgelabel[edgecolindex]
                        }
                    }
                    edges[[i]]@attrs$arrowhead <- "open"
                    edges[[i]]@attrs$color <- "black"
                    if (!is.null(edgecol)) {
                        if (!is.na(edgecol[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$color <- edgecol[edgecolindex]
                        }
                    }
                    if (!is.null(edgehead)) {
                        if (!is.na(edgehead[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex]
                        }
                    }
                    if (!is.null(edgetail)) {
                        if (!is.na(edgetail[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex]
                        }
                    }
                }
            }
        }
    }
    if (type == 1) {
        
        g2 <- agopen(name="boolean", nodes=nodes, recipEdges = "distinct", edges=edges, edgeMode="undirected", attrs=list(edge = list(), graph = list(lwd = lwd, rankdir = rankdir), node=list(lwd = lwd, fixedsize=FALSE)))

        plot(g2, "dot", lwd = lwd, ...)

    } else {

        arrowheads <- character()
        arrowtails <- character()
        arrowcolors <- character()
        arrowlabels <- character()
        arrowlwd <- character()
        arrowlty <- character()
        arrowfontsize <- character()
        if (length(edges) > 0) {
            for (i in 1:length(edges)) {
                if (length(edges[[i]]@attrs$style) == 0) { edges[[i]]@attrs$style <- "solid" }
                arrowlty <- c(arrowlty, edges[[i]]@attrs$style)
                arrowheads <- c(arrowheads, edges[[i]]@attrs$arrowhead)
                if (!is.null(edgetail)) {
                    arrowtails <- c(arrowtails, edges[[i]]@attrs$arrowtail)
                } else {
                    arrowtails <- c(arrowtails, "none")
                }
                arrowcolors <- c(arrowcolors, edges[[i]]@attrs$color)
                arrowfontsize <- c(arrowfontsize, edges[[i]]@attrs$fontsize)
                arrowlwd <- c(arrowlwd, edges[[i]]@attrs$weight)
                arrowlabels <- c(arrowlabels, edges[[i]]@attrs$label)
            }
        }
        arrowlwd <- as.numeric(arrowlwd)
        
        graph.trans <- NULL
        and.count <- 0
        for (i in dnf) {
            if (length(grep("\\+", i)) > 0) {
                and.count <- and.count + 1
                output <- unlist(strsplit(i, "="))
                input <- unlist(strsplit(output[1], "\\+"))
                output <- output[2]
                for (i in input) {
                    graph.trans <- c(graph.trans, paste(gsub("!", "", i), "~", paste("and", and.count, sep = ""), sep = ""))
                }
                graph.trans <- c(graph.trans, paste(paste("and", and.count, sep = ""), "~", output, sep = ""))
            } else {
                put <- unlist(strsplit(i, "="))
                graph.trans <- c(graph.trans, paste(gsub("!", "", put[1]), "~", put[2], sep = ""))
            }
        }

        if (length(edgecol) == length(arrowcolors)) {
            edgecol <- edgecol[order(match(graph.trans, names(edges)))]
            arrowcolors <- edgecol
        }
        
        nodeshapes <- character()
        nodecolors <- character()
        nodeheight <- character()
        nodewidth <- character()
        nodecolor <- character()
        
        for (i in 1:length(nodes)) {
            nodeshapes <- c(nodeshapes, nodes[[i]]@attrs$shape)
            nodecolors <- c(nodecolors, nodes[[i]]@attrs$fillcolor)
            nodeheight <- c(nodeheight, nodes[[i]]@attrs$height)
            nodewidth <- c(nodewidth, nodes[[i]]@attrs$width)
            nodecolor <- c(nodecolor, nodes[[i]]@attrs$color)
        }

        nodeheight[which(nodeheight == "0.4")] <- "0.2"

        if (is.null(lty) & is.null(edgestyle)) {
            arrowlty <- rep("solid", length(edges))
        }
        
        if (use.freq) {
            if (is.null(lty)) {
                if (edgestyle) {
                    arrowlty[which(as.numeric(arrowlabels) < 66)] <- "dashed"
                    arrowlty[which(as.numeric(arrowlabels) < 33)] <- "dotted"
                }
            }
            arrowlwd <- arrowlwd - min(arrowlwd)
            arrowlwd <- as.character((arrowlwd/max(arrowlwd)+0.1)*2*edgelwd)
        } else {
            if (is.null(edgewidth)) {
                arrowlwd <- rep(edgelwd, length(edges))
            }
        }

        if (is.null(edgewidth) & is.null(edgelabel)) {
            arrowlabels <- rep("", length(edges))
        }
        
        if (length(arrowlty) == 0) {
            arrowlty <- rep("solid", length(edges))
        }
        if (length(arrowlwd) == 0) {
            arrowlwd <- rep(lwd, length(edges))
        }
        
        names(arrowfontsize) <- names(arrowheads) <- names(arrowtails) <- names(arrowcolors) <- names(arrowlwd) <- names(arrowlty) <- names(arrowlabels) <- names(edges)

        names(nodecolor) <- names(nodewidth) <- names(nodeheight) <- names(nodeshapes) <- names(nodecolors) <- names(nodes)

        if (length(unique(names(edges))) < length(names(edges))) {
            for (i in names(edges)[-which(duplicated(names(edges)) == TRUE)]) {
                getpos <- grep(paste("^", i, "$", sep = ""), names(edges))
                if (length(getpos) > 1) {
                    if (use.freq) {
                        if (arrowheads[getpos[1]] %in% "tee") {
                            arrowlabels[getpos[1]] <- paste(paste(c("-", "+"), arrowlabels[getpos], sep = ""), collapse = "\n")
                        } else {
                            arrowlabels[getpos[1]] <- paste(paste(c("+", "-"), arrowlabels[getpos], sep = ""), collapse = "\n")
                        }
                    } else {
                        if (is.null(edgelabel)) {
                            arrowlabels[getpos[1]] <- ""
                        }
                    }
                    arrowheads[getpos[1]] <- "odiamond"
                    if (is.null(edgecol)) {
                        arrowcolors[getpos[1]] <- "black"
                    } else {
                        if (is.na(edgecol[getpos[1]])) {
                            arrowcolors[getpos[1]] <- "black"
                        }
                    }
                    arrowlwd[getpos[1]] <- as.character(mean(as.numeric(arrowlwd[getpos])))
                }
            }
        }
        if (length(labels) == length(arrowlabels) & is.null(edgelabel)) {
            arrowlabels[!is.na(labels)] <- labels[!is.na(labels)]
        }
        if (length(edgecol) == 1) {
            arrowcolors <- rep(edgecol, length(arrowcolors))
            names(arrowcolors) <- names(arrowlabels)
        }
        
        if (legend == 1 | legend == 3) {
            if (dolegend) {
                start <- 1
                g@nodes <- c("LEGEND:", "STIMULUS", "INHIBITOR", "SIGNAL", "NOTHING", "active", "inactive")
                g@edgeL <- list()
                g@edgeData@data <- list()
            } else {
                start <- length(g@nodes) + 1
                g@nodes <- c(g@nodes, "LEGEND:", "STIMULUS", "INHIBITOR", "SIGNAL", "NOTHING", "active", "inactive")
            }
            g@edgeL[["LEGEND:"]] <- list()
            g@edgeL[["STIMULUS"]] <- list()
            g@edgeL[["INHIBITOR"]] <- list()
            g@edgeL[["SIGNAL"]] <- list()
            g@edgeL[["NOTHING"]] <- list()
            g@edgeL[["active"]] <- list()
            g@edgeL[["inactive"]] <- list()
            g@edgeL[["LEGEND:"]][["edges"]] <- as.integer(start+1)
            g@edgeL[["STIMULUS"]][["edges"]] <- as.integer(start+2)
            g@edgeL[["INHIBITOR"]][["edges"]] <- as.integer(start+3)
            g@edgeL[["SIGNAL"]][["edges"]] <- as.integer(start+4)
            g@edgeL[["NOTHING"]][["edges"]] <- as.integer(start+5)
            g@edgeL[["active"]][["edges"]] <- c(as.integer(start+6), as.integer(start+4))
            g@edgeL[["inactive"]][["edges"]] <- as.integer(start+5)
            g@edgeData@data[["LEGEND:|STIMULUS"]] <- list()
            g@edgeData@data[["STIMULUS|INHIBITOR"]] <- list()
            g@edgeData@data[["INHIBITOR|SIGNAL"]] <- list()
            g@edgeData@data[["SIGNAL|NOTHING"]] <- list()
            g@edgeData@data[["NOTHING|active"]] <- list()
            g@edgeData@data[["active|inactive"]] <- list()
            g@edgeData@data[["active|NOTHING"]] <- list()
            g@edgeData@data[["inactive|active"]] <- list()
            g@edgeData@data[["LEGEND:|STIMULUS"]][["weight"]] <- 1
            g@edgeData@data[["STIMULUS|INHIBITOR"]][["weight"]] <- 1
            g@edgeData@data[["INHIBITOR|SIGNAL"]][["weight"]] <- 1
            g@edgeData@data[["SIGNAL|NOTHING"]][["weight"]] <- 1
            g@edgeData@data[["NOTHING|active"]][["weight"]] <- 1
            g@edgeData@data[["active|inactive"]][["weight"]] <- 1
            g@edgeData@data[["active|NOTHING"]][["weight"]] <- 1
            g@edgeData@data[["inactive|active"]][["weight"]] <- 1
            arrowheads <- c(arrowheads, "LEGEND:~STIMULUS" = "none", "STIMULUS~INHIBITOR" = "open", "INHIBITOR~SIGNAL" = "tee", "SIGNAL~NOTHING" = "odiamond", "NOTHING~active" = "open", "active~inactive" = "tee", "active~NOTHING" = "tee", "inactive~active" = "open")
            arrowcolors <- c(arrowcolors, "LEGEND:~STIMULUS" = "transparent", "STIMULUS~INHIBITOR" = "black", "INHIBITOR~SIGNAL" = "red", "SIGNAL~NOTHING" = "blue", "NOTHING~active" = "black", "active~inactive" = "red", "active~NOTHING" = "red", "inactive~active" = "black")
            arrowlabels <- c(arrowlabels, "LEGEND:~STIMULUS" = "", "STIMULUS~INHIBITOR" = "    positive", "INHIBITOR~SIGNAL" = "    negative", "SIGNAL~NOTHING" = "    ambiguous\npositive\nnegative", "NOTHING~active" = "    bidirectional\ndifferent", "active~inactive" = "    bidirectional\ndifferent", "active~NOTHING" = "", "inactive~active" = "")
            nodecolors <- c(nodecolors, "LEGEND:" = "white", "STIMULUS" = "white", "INHIBITOR" = "white", "SIGNAL" = "lightblue", "NOTHING" = "white", "active" = "green", "inactive" = "white")
            nodeheight <- c(nodeheight, "LEGEND:" = 0, "STIMULUS" = as.character(max(nodeheight)), "INHIBITOR" = as.character(max(nodeheight)), "SIGNAL" = as.character(max(nodeheight)), "NOTHING" = as.character(max(nodeheight)), "active" = as.character(max(nodeheight)), "inactive" = as.character(max(nodeheight)))
            nodewidth <- c(nodewidth, "LEGEND:" = as.character(max(nodewidth)), "STIMULUS" = as.character(max(nodewidth)), "INHIBITOR" = as.character(max(nodewidth)), "SIGNAL" = as.character(max(nodewidth)), "NOTHING" = as.character(max(nodewidth)), "active" = as.character(max(nodewidth)), "inactive" = as.character(max(nodewidth)))
            if (type == 2) {
                nodeshapes <- c(nodeshapes, "LEGEND:" = "box", "STIMULUS" = "diamond", "INHIBITOR" = "ellipse", "SIGNAL" = "ellipse", "NOTHING" = "ellipse", "active" = "ellipse", "inactive" = "ellipse")
            } else {
                nodeshapes <- c(nodeshapes, "LEGEND:" = "box", "STIMULUS" = "box", "INHIBITOR" = "ellipse", "SIGNAL" = "ellipse", "NOTHING" = "ellipse", "active" = "ellipse", "inactive" = "ellipse")
            }
            nodecolor <- c(nodecolor, "LEGEND:" = "white", "STIMULUS" = "black", "INHIBITOR" = "red", "SIGNAL" = "black", "NOTHING" = "black", "active" = "black", "inactive" = "black")
            dnf <- c(dnf, "NOTHING=active", "!active=NOTHING", "!active=inactive", "inactive=active")
        }
        nodelabels <- names(nodecolor)
        names(nodelabels) <- nodelabels
        for (i in 1:length(nodelabel)) {
            nodelabels[which(names(nodelabels) %in% names(nodelabel)[i])] <- nodelabel[i]
        }
        nodefontsizes <- NULL
        if (!is.null(nodefontsize)) {
            nodefontsizes <- rep(14, length(nodelabels))
            names(nodefontsizes) <- names(nodelabels)
            for (i in 1:length(nodefontsize)) {
                nodefontsizes[which(names(nodefontsizes) %in% names(nodefontsize)[i])] <- nodefontsize[[i]]
            }
        }
        g <- layoutGraph(g, edgeAttrs = list(arrowhead = arrowheads, color = arrowcolors, label = arrowlabels, arrowtail = arrowtails), nodeAttrs = list(labels = nodelabels, color = nodecolor, height = nodeheight, width = nodewidth, shape = nodeshapes, fillcolor = nodecolors), layoutType=layout)
        graph.par(list(graph=list(main = main, sub = sub, cex.main = cex.main, cex.sub = cex.sub, col.sub = col.sub), edges=list(textCol = labelcol, lwd = edgelwd, fontsize = labelsize), nodes=list(lwd = lwd, fontsize = fontsize, cex = cex)))
        edgeRenderInfo(g) <- list(lty = arrowlty, lwd = arrowlwd, label = arrowlabels)
        if (length(edges) > 0) {
            for (i in names(g@renderInfo@edges$direction)) {
                input <- unlist(strsplit(i, "~"))
                output <- input[2]
                input <- input[1]
                ambig <- FALSE
                if (paste("!", input, "=", output, sep = "") %in% dnf & paste("", input, "=", output, sep = "") %in% dnf) {
                    ambig <- TRUE
                }
                if ((length(grep("and", i)) == 0 & g@renderInfo@edges$direction[[i]] == "both") | ambig) {
                    pos <- which(names(g@renderInfo@edges$arrowhead) %in% i)
                    if (is.null(edgehead)) {
                        if (paste("!", input, "=", output, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowhead[pos] <- "tee"
                        }
                        if (paste(input, "=", output, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowhead[pos] <- "open"
                        }
                        if (paste("!", output, "=", input, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowtail[pos] <- "tee"
                        }
                        if (paste(output, "=", input, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowtail[pos] <- "open"
                        }
                        if (paste("!", output, "=", input, sep = "") %in% dnf & paste("", output, "=", input, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowtail[pos] <- "odiamond"
                        }
                        if (paste("!", input, "=", output, sep = "") %in% dnf & paste("", input, "=", output, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowhead[pos] <- "odiamond"
                        }
                    }
                    if (is.null(edgecol)) {
                        if (g@renderInfo@edges$arrowtail[pos] == "open" & g@renderInfo@edges$arrowhead[pos] == "open") {
                            g@renderInfo@edges$col[pos] <- "black"
                        }
                        if (g@renderInfo@edges$arrowtail[pos] == "tee" & g@renderInfo@edges$arrowhead[pos] == "tee") {
                            g@renderInfo@edges$col[pos] <- "red"
                        }
                        if (g@renderInfo@edges$arrowtail[pos] != g@renderInfo@edges$arrowhead[pos]) {
                            g@renderInfo@edges$col[pos] <- "brown"
                        }
                        if (g@renderInfo@edges$arrowtail[pos] == "odiamond" | g@renderInfo@edges$arrowhead[pos] == "odiamond") {
                            g@renderInfo@edges$col[pos] <- "blue"
                        }
                    } else {
                        if (is.null(edgecol)) { # is.na(edgecol[pos])
                            if (g@renderInfo@edges$arrowtail[pos] == "open" & g@renderInfo@edges$arrowhead[pos] == "open") {
                                g@renderInfo@edges$col[pos] <- "black"
                            }
                            if (g@renderInfo@edges$arrowtail[pos] == "tee" & g@renderInfo@edges$arrowhead[pos] == "tee") {
                                g@renderInfo@edges$col[pos] <- "red"
                            }
                            if (g@renderInfo@edges$arrowtail[pos] != g@renderInfo@edges$arrowhead[pos]) {
                                g@renderInfo@edges$col[pos] <- "brown"
                            }
                            if (g@renderInfo@edges$arrowtail[pos] == "odiamond" | g@renderInfo@edges$arrowhead[pos] == "odiamond") {
                                g@renderInfo@edges$col[pos] <- "blue"
                            }
                        }
                    }
                }
            }
        }
        if (!is.null(simulate$draw)) {
            for (i in simulate$inhibitors) {
                ## add the inhibiting node
                g@nodes <- c(g@nodes, paste(i, "_inhibited", sep = ""))
                g@renderInfo@nodes$nodeX <- c(g@renderInfo@nodes$nodeX, g@renderInfo@nodes$nodeX[which(names(g@renderInfo@nodes$nodeX) %in% i)])
                g@renderInfo@nodes$nodeY <- c(g@renderInfo@nodes$nodeY, g@renderInfo@nodes$nodeY[which(names(g@renderInfo@nodes$nodeY) %in% i)])
                names(g@renderInfo@nodes$nodeX)[length(g@renderInfo@nodes$nodeX)] <- paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$nodeY)[length(g@renderInfo@nodes$nodeY)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelX <- c(g@renderInfo@nodes$labelX, g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + 200)
                g@renderInfo@nodes$labelY <- c(g@renderInfo@nodes$labelY, g@renderInfo@nodes$labelY[which(names(g@renderInfo@nodes$labelY) %in% i)])
                names(g@renderInfo@nodes$labelX)[length(g@renderInfo@nodes$labelX)] <- paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$labelY)[length(g@renderInfo@nodes$labelY)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelJust <- c(g@renderInfo@nodes$labelJust, "n")
                names(g@renderInfo@nodes$labelJust)[length(g@renderInfo@nodes$labelJust)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$col <- c(g@renderInfo@nodes$col, "transparent")
                names(g@renderInfo@nodes$col)[length(g@renderInfo@nodes$col)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$fill <- c(g@renderInfo@nodes$fill, "transparent")
                names(g@renderInfo@nodes$fill)[length(g@renderInfo@nodes$fill)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$shape <- c(g@renderInfo@nodes$shape, "box")
                names(g@renderInfo@nodes$shape)[length(g@renderInfo@nodes$shape)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$style <- c(g@renderInfo@nodes$style, "")
                names(g@renderInfo@nodes$style)[length(g@renderInfo@nodes$style)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$height <- c(g@renderInfo@nodes$height, 2)
                names(g@renderInfo@nodes$height)[length(g@renderInfo@nodes$height)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$rWidth <- c(g@renderInfo@nodes$rWidth, 1)
                names(g@renderInfo@nodes$rWidth)[length(g@renderInfo@nodes$rWidth)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$lWidth <- c(g@renderInfo@nodes$lWidth, 1)
                names(g@renderInfo@nodes$lWidth)[length(g@renderInfo@nodes$lWidth)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$label[[paste(i, "_inhibited", sep = "")]] <- ""#"X            X\n X       X\n  X  X\n   X\n  X  X\n X       X\nX            X"
                g@renderInfo@nodes$labelWidth <- c(g@renderInfo@nodes$labelWidth, g@renderInfo@nodes$labelWidth[which(names(g@renderInfo@nodes$labelWidth) %in% i)])
                names(g@renderInfo@nodes$labelWidth)[length(g@renderInfo@nodes$labelWidth)] <- paste(i, "_inhibited", sep = "")

                ## add the inhibiting edge
                tmp.name <- paste(i, "_inhibited", sep = "")
                g@edgeL[[tmp.name]] <- list()
                g@edgeL[[tmp.name]][["edges"]] <- which(g@nodes %in% i)
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]] <- list()
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]]$weight <- 1
                g@renderInfo@edges$enamesFrom <- c(g@renderInfo@edges$enamesFrom, tmp.name)
                names(g@renderInfo@edges$enamesFrom)[length(g@renderInfo@edges$enamesFrom)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$enamesTo <- c(g@renderInfo@edges$enamesTo, i)
                names(g@renderInfo@edges$enamesTo)[length(g@renderInfo@edges$enamesTo)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelJust <- c(g@renderInfo@edges$labelJust, NA)
                names(g@renderInfo@edges$labelJust)[length(g@renderInfo@edges$labelJust)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelX <- c(g@renderInfo@edges$labelX, NA)
                names(g@renderInfo@edges$labelX)[length(g@renderInfo@edges$labelX)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelY <- c(g@renderInfo@edges$labelY, NA)
                names(g@renderInfo@edges$labelY)[length(g@renderInfo@edges$labelY)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelWidth <- c(g@renderInfo@edges$labelWidth, NA)
                names(g@renderInfo@edges$labelWidth)[length(g@renderInfo@edges$labelWidth)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$label <- c(g@renderInfo@edges$label, "")
                names(g@renderInfo@edges$label)[length(g@renderInfo@edges$label)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowhead <- c(g@renderInfo@edges$arrowhead, "tee")
                names(g@renderInfo@edges$arrowhead)[length(g@renderInfo@edges$arrowhead)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowtail <- c(g@renderInfo@edges$arrowtail, "odot")
                names(g@renderInfo@edges$arrowtail)[length(g@renderInfo@edges$arrowtail)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$col <- c(g@renderInfo@edges$col, "firebrick")
                names(g@renderInfo@edges$col)[length(g@renderInfo@edges$col)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lwd <- c(g@renderInfo@edges$lwd, lwd[1]*1)
                names(g@renderInfo@edges$lwd)[length(g@renderInfo@edges$lwd)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lty <- c(g@renderInfo@edges$lty, "solid")
                names(g@renderInfo@edges$lty)[length(g@renderInfo@edges$lty)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$direction <- c(g@renderInfo@edges$direction, "forward")
                names(g@renderInfo@edges$direction)[length(g@renderInfo@edges$direction)] <- paste(tmp.name, "~", i, sep = "")
                ## calculate splines
                tmp.splines <- rep(g@renderInfo@nodes$labelY[which(names(g@renderInfo@nodes$labelY) %in% i)], 8)
                tmp.splines[c(7,5,3,1)] <- round(seq(g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + g@renderInfo@nodes$rWidth[[i]] + 10,
                                                     g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + 200 - g@renderInfo@nodes$rWidth[[i]] - 10,
                                                     length.out = 4))
                tmp.splines[2] <- tmp.splines[2] + 50
                tmp.splines <- as.integer(tmp.splines)
                g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]] <- g@renderInfo@edges$splines[[1]]
                for (j in 1:4) {
                    g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]][[1]]@cPoints[[j]]@x <- tmp.splines[c(1,3,5,7)][j]
                    g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]][[1]]@cPoints[[j]]@y <- tmp.splines[c(2,4,6,8)][j]
                }
            }
            for (i in simulate$stimuli) {
                ## add the stimulating node
                g@nodes <- c(g@nodes, paste(i, "_inhibited", sep = ""))
                g@renderInfo@nodes$nodeX <- c(g@renderInfo@nodes$nodeX, g@renderInfo@nodes$nodeX[which(names(g@renderInfo@nodes$nodeX) %in% i)])
                g@renderInfo@nodes$nodeY <- c(g@renderInfo@nodes$nodeY, g@renderInfo@nodes$nodeY[which(names(g@renderInfo@nodes$nodeY) %in% i)])
                names(g@renderInfo@nodes$nodeX)[length(g@renderInfo@nodes$nodeX)] <- paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$nodeY)[length(g@renderInfo@nodes$nodeY)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelX <- c(g@renderInfo@nodes$labelX, g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + 200)
                g@renderInfo@nodes$labelY <- c(g@renderInfo@nodes$labelY, g@renderInfo@nodes$labelY[which(names(g@renderInfo@nodes$labelY) %in% i)])
                names(g@renderInfo@nodes$labelX)[length(g@renderInfo@nodes$labelX)] <- paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$labelY)[length(g@renderInfo@nodes$labelY)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelJust <- c(g@renderInfo@nodes$labelJust, "n")
                names(g@renderInfo@nodes$labelJust)[length(g@renderInfo@nodes$labelJust)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$col <- c(g@renderInfo@nodes$col, "transparent")
                names(g@renderInfo@nodes$col)[length(g@renderInfo@nodes$col)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$fill <- c(g@renderInfo@nodes$fill, "transparent")
                names(g@renderInfo@nodes$fill)[length(g@renderInfo@nodes$fill)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$shape <- c(g@renderInfo@nodes$shape, "box")
                names(g@renderInfo@nodes$shape)[length(g@renderInfo@nodes$shape)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$style <- c(g@renderInfo@nodes$style, "")
                names(g@renderInfo@nodes$style)[length(g@renderInfo@nodes$style)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$height <- c(g@renderInfo@nodes$height, 2)
                names(g@renderInfo@nodes$height)[length(g@renderInfo@nodes$height)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$rWidth <- c(g@renderInfo@nodes$rWidth, 1)
                names(g@renderInfo@nodes$rWidth)[length(g@renderInfo@nodes$rWidth)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$lWidth <- c(g@renderInfo@nodes$lWidth, 1)
                names(g@renderInfo@nodes$lWidth)[length(g@renderInfo@nodes$lWidth)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$label[[paste(i, "_inhibited", sep = "")]] <- ""#"X            X\n X       X\n  X  X\n   X\n  X  X\n X       X\nX            X"
                g@renderInfo@nodes$labelWidth <- c(g@renderInfo@nodes$labelWidth, g@renderInfo@nodes$labelWidth[which(names(g@renderInfo@nodes$labelWidth) %in% i)])
                names(g@renderInfo@nodes$labelWidth)[length(g@renderInfo@nodes$labelWidth)] <- paste(i, "_inhibited", sep = "")

                ## add the stimulating edge
                tmp.name <- paste(i, "_inhibited", sep = "")
                g@edgeL[[tmp.name]] <- list()
                g@edgeL[[tmp.name]][["edges"]] <- which(g@nodes %in% i)
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]] <- list()
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]]$weight <- 1
                g@renderInfo@edges$enamesFrom <- c(g@renderInfo@edges$enamesFrom, tmp.name)
                names(g@renderInfo@edges$enamesFrom)[length(g@renderInfo@edges$enamesFrom)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$enamesTo <- c(g@renderInfo@edges$enamesTo, i)
                names(g@renderInfo@edges$enamesTo)[length(g@renderInfo@edges$enamesTo)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelJust <- c(g@renderInfo@edges$labelJust, NA)
                names(g@renderInfo@edges$labelJust)[length(g@renderInfo@edges$labelJust)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelX <- c(g@renderInfo@edges$labelX, NA)
                names(g@renderInfo@edges$labelX)[length(g@renderInfo@edges$labelX)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelY <- c(g@renderInfo@edges$labelY, NA)
                names(g@renderInfo@edges$labelY)[length(g@renderInfo@edges$labelY)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelWidth <- c(g@renderInfo@edges$labelWidth, NA)
                names(g@renderInfo@edges$labelWidth)[length(g@renderInfo@edges$labelWidth)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$label <- c(g@renderInfo@edges$label, "")
                names(g@renderInfo@edges$label)[length(g@renderInfo@edges$label)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowhead <- c(g@renderInfo@edges$arrowhead, "open")
                names(g@renderInfo@edges$arrowhead)[length(g@renderInfo@edges$arrowhead)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowtail <- c(g@renderInfo@edges$arrowtail, "odot")
                names(g@renderInfo@edges$arrowtail)[length(g@renderInfo@edges$arrowtail)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$col <- c(g@renderInfo@edges$col, "limegreen")
                names(g@renderInfo@edges$col)[length(g@renderInfo@edges$col)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lwd <- c(g@renderInfo@edges$lwd, lwd[1]*1)
                names(g@renderInfo@edges$lwd)[length(g@renderInfo@edges$lwd)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lty <- c(g@renderInfo@edges$lty, "solid")
                names(g@renderInfo@edges$lty)[length(g@renderInfo@edges$lty)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$direction <- c(g@renderInfo@edges$direction, "forward")
                names(g@renderInfo@edges$direction)[length(g@renderInfo@edges$direction)] <- paste(tmp.name, "~", i, sep = "")
                ## calculate splines
                tmp.splines <- rep(g@renderInfo@nodes$labelY[which(names(g@renderInfo@nodes$labelY) %in% i)], 8)
                tmp.splines[c(7,5,3,1)] <- round(seq(g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + g@renderInfo@nodes$rWidth[[i]] + 10,
                                                     g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + 200 - g@renderInfo@nodes$rWidth[[i]] - 10,
                                                     length.out = 4))
                tmp.splines[2] <- tmp.splines[2] + 50
                tmp.splines <- as.integer(tmp.splines)
                g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]] <- g@renderInfo@edges$splines[[1]]
                for (j in 1:4) {
                    g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]][[1]]@cPoints[[j]]@x <- tmp.splines[c(1,3,5,7)][j]
                    g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]][[1]]@cPoints[[j]]@y <- tmp.splines[c(2,4,6,8)][j]
                }
            }
            g@renderInfo@graph$bbox[2,1] <- g@renderInfo@graph$bbox[2,1] + 150
            g@renderInfo@graph$bbox[2,2] <- g@renderInfo@graph$bbox[2,2] + 25
        }
        g <- g
        if (draw) {
            renderGraph(g, lwd = lwd, recipEdges = "distinct", ...)
        }
    }
    if (legend == 2 | legend == 3) {
        legend(x = x, y = y, legend = c("signals are blue", "stimuli are diamonds/boxes", "inhibitors have a red border", "positive regulation is green ->", "negative regulation is red -|", "ambiguous regulation is black -o"), fill = c("lightblue", "white", "red", "green", "red", "black"), col = c("lightblue", "white", "red", "green", "red", "black"), yjust = yjust, xjust = xjust)
    }
    return(g)
}


###--- END OF HELPER FUNCTIONS ---###
