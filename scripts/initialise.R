# can be run from either within the scripts or the examples directory
library(scater)
library(scran)
library(ggplot2)
library(patchwork)
library(batchelor)
library(matrixStats)

library(StabMap)

source("../scripts/adaptiveKNN.R")
source("../scripts/batch_wrap_functions.R")
source("../scripts/convenience_functions.R")
source("../scripts/comparison_functions.R")
source("../scripts/subsetUncertainty.R")
source("../scripts/tfidf.R")

# give each method it's own colour
library(nord)
method_colours = c(nord(palette = "lumina")[c(4,2)],
                   nord(palette = "rocky_mountain")[c(1)],
                   nord(palette = "lake_superior")[c(3)],
                   nord(palette = "lumina")[c(1)])
names(method_colours) = c("StabMap", "PCA", "MultiMAP", "UINMF", "PC_separate")
