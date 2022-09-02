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
method_colours = c(nord(palette = "lumina")[c(4,4,4,4,4,2,2,2,2,2,2)],
                   nord(palette = "rocky_mountain")[c(1,1)],
                   nord(palette = "lake_superior")[c(3,3)],
                   nord(palette = "lumina")[c(1)],
                   nord(palette = "rocky_mountain")[c(2)])
names(method_colours) = c("StabMap", "StabMap_Harmony","StabMap_fastMNN","StabMap_Seurat","StabMap_embedding_corrected_MNN",
                          "fastMNN", "PCA", "Naive_Harmony", "Naive_Seurat", "Naive_fastMNN","PC_embedding_corrected_MNN",
                          "MultiMAP","MultiMAP_embedding",
                          "UINMF","UINMF_embedding",
                          "PC_separate", "Seurat")
