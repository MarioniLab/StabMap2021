# can be run from either within the scripts or the examples directory
library(scater)
library(scran)
library(ggplot2)
library(patchwork)
library(batchelor)
library(matrixStats)
source("../scripts/adaptiveKNN.R")
source("../scripts/StabMap_functions.R")
source("../scripts/StabMap_generalised_experimental.R")
source("../scripts/subsetUncertainty.R")

# give each method it's own colour
library(nord)
method_colours = c(nord(palette = "lumina")[c(4,2)],
                   nord(palette = "rocky_mountain")[c(1)],
                   nord(palette = "lake_superior")[c(3)],
                   nord(palette = "lumina")[c(1)],
                   nord(palette = "lumina")[c(5)],
                   nord(palette = "lumina")[c(4)])
names(method_colours) = c("StabMap", "PCA", "MultiMAP", "UINMF", "PC_separate", "StabMap_labels",
                          "StabMap_unprojected")
