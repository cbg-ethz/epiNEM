NEMlist <- load('NEMlist.RData')

parallel <- NULL
stimuli <- "Stim"
inhibitors <- c("PPH3", "DUN1", "SSK2")
signals <- c("PPH3", "DUN1", "SSK2")

CNOlist <- dummyCNOlist(stimuli = stimuli, inhibitors = inhibitors, maxStim = 1, maxInhibit = 3)

Sgenes <- c(stimuli, inhibitors)

sifMatrix <- numeric()
for (i in Sgenes) {
  for (j in Sgenes) {
    if (i %in% j) { next() }
    sifMatrix <- rbind(sifMatrix, c(i, 1, j))
    sifMatrix <- rbind(sifMatrix, c(i, -1, j)) # falls du auch negative regulation zulassen möchtest
  }
}
write.table(sifMatrix, file = "temp.sif", sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE)
PKN <- readSIF("temp.sif")
unlink("temp.sif")

model <- preprocessing(CNOlist, PKN, maxInputsPerGate=100)

initBstring <- rep(0, length(model$reacID))

locRun <- localSearch(
  CNOlist=CNOlist,
  NEMlist=NEMlist,
  model=model,
  parallel=parallel,
  initSeed=initBstring,
  draw = TRUE # FALSE does not draw the network evolution and can be faster
)