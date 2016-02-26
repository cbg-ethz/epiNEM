## read in genes from user and prepare for further analysis
readInGenes <- function(genes){
  return(unlist(lapply(genes, function(x) {paste(x, '.del.vs..wt.1', sep='')})))
}

## binarise data: p-value < 0.05 significant, else not
## set user specific p-value?!
binarise <- function(x) {
  ifelse((x < 0.05 & x > 0), 1, 0)
}

## numericalize factors
asNumericFactor <- function(x) {as.numeric(levels(x))[x]}

binarizeData <- function(data) {
  binarized_data <- data %>%
    #select_(.dots = single) %>%
    select(-c(systematicName, geneSymbol)) %>%
    select(matches('.DEL.VS..WT.1')) %>%
    filter(-1) %>%
    mutate_each(funs(droplevels)) %>%
    mutate_each(funs(asNumericFactor)) %>%
    mutate_each(funs(binarise)) 
  return(binarized_data)
}

capitalize <- function(x){
  setnames(x, names(x), toupper(names(x)) %>% 
             gsub('.DEL.VS..WT.1', '', .) %>% 
             #gsub('.DEL.VS..WT', '', .) %>%
             gsub('.DEL.', '.', .)); 
  return(x)
}

findZeroRows <- function(binary_matrix){
  row_sums <- lapply(1:dim(binary_matrix)[1], function(x){rowSums(binary_matrix[x,])})
  zero_row <- which(row_sums==0)
  notzero_row <- which(!row_sums==0)
  result_list <- list(zerorow=zero_row, notzero=notzero_row)
  return(result_list)
}

removeZeroRows <- function(matrix, zero_row_index){
  matrix %<>%
    filter(-zero_row_index) 
  return(matrix)
}

addGeneName <-function(matrix, caption_matrix, not_zero_entry){
  matrix %>%
  mutate(geneName = caption_matrix$systematicName[not_zero_entry])
}