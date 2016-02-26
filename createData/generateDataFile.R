rm(list = ls())

library(magrittr)
library(data.table)
library(dplyr)

source("createData/auxFunction.R")

address <- "http://www.holstegelab.nl/publications/sv/signaling_redundancy/downloads/DataS1.txt"
data <- data.table(read.csv(url(address), sep="\t", header=TRUE))

pathway_data <- binarizeData(data)
pathway_data_capitalized <- capitalize(pathway_data) 

zero_row_non_zero <- findZeroRows(pathway_data_capitalized)

binary_matrix_without_rows <-removeZeroRows(pathway_data_capitalized, zero_row_non_zero$zerorow)

data_without_caption <- data %>%  filter(-1) 
full_matrix <- addGeneName(binary_matrix_without_rows, data_without_caption , zero_row_non_zero$notzero)
