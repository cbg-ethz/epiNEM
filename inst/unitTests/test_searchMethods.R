data <- matrix(sample(c(0,1), 100*4, replace = TRUE), 100, 4)
colnames(data) <- c("A", "A.B", "B", "C")
rownames(data) <- paste("E", 1:100, sep = "_")
res1 <- epiNEM(data, method = "exhaustive")
res2 <- epiNEM(data, method = "greedy", nIterations = 1,
               init = res1$origModel)
checkEquals(res1, res2, checkNames = TRUE)
