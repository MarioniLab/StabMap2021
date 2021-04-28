# wt is wild type chimera data
library(MouseGastrulationData)
library(SpatialUtils)
library(scater)
source("/Users/ghazan01/Dropbox/Backup/StabMap/StabMap2021/scripts/adaptiveKNN.R")
source("/Users/ghazan01/Dropbox/Backup/StabMap/StabMap2021/scripts/StabMap_functions.R")

##### data

wt = WTChimeraData(type = "processed", samples = c("6", "8", "10"))
setwd("/Users/ghazan01/Dropbox/Backup/SpatialEmbryos/scripts")
atlas_sub = readRDS("../analysis_output/E8.5/atlas_sub.Rds")
atlas = readRDS("../analysis_output/E8.5/atlas_for_celltype_mapping.Rds")
sce_filt = readRDS("../analysis_output/E8.5/E8.5_sce_filt_unlabelled.Rds")
sce_filt
sce_filt <- reNormalise(sce_filt)
mapping_dt = updateMesenchymeLabel(readRDS("../analysis_output/E8.5/celltype_annotation_refined.Rds"))
colData(sce_filt) <- cbind(colData(sce_filt), mapping_dt[colnames(sce_filt), setdiff(colnames(mapping_dt),colnames(colData(sce_filt)))])
source("celltype_colours.R")

ct_all = c(as.character(sce_filt$celltype_mapped_refined), as.character(atlas_sub$celltype_parsed))
names(ct_all) <- c(colnames(sce_filt), colnames(atlas_sub))



#### build params

nGenes_all = rep(c(100, 250, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000),
                 each = 3)

res = NULL

for (nGenes in nGenes_all) {
  
  print(nGenes)
  
  doUMAP = FALSE
  # referenceSCE = atlas_sub
  referenceSCE = atlas[,atlas$sample == 17]
  querySCE = atlas[,atlas$sample == 29]
  # querySCE = sce_filt
  # grouping = "celltype_parsed"
  grouping = "celltype"
  # groupingQuery = "celltype_mapped_refined"
  # genes = rownames(sce_filt)[1:50]
  assayNameReference = "logcounts"
  assayNameQuery = "logcounts"
  
  # genes = intersect(sample(rownames(referenceSCE), nGenes),
  #                   rownames(referenceSCE)[rowVars(assay(referenceSCE, assayNameReference)) > 0])
  genes = sample(rownames(referenceSCE)[rowVars(assay(referenceSCE, assayNameReference)) > 0])[seq_len(nGenes)]
  
  type = factor(ifelse(c(colnames(referenceSCE), colnames(querySCE)) %in% colnames(referenceSCE),
                       "reference", "query"))
  names(type) <- c(colnames(referenceSCE), colnames(querySCE))
  
  
  # rm(LD_embedding)
  # rm(PC_embedding)
  # rm(SC_embedding)
  
  LD_embedding = stabMapLabelled(referenceSCE = referenceSCE,
                                   querySCE = querySCE,
                                   genes = genes,
                                   grouping = grouping,
                                   assayNameReference = assayNameReference,
                                   assayNameQuery = assayNameQuery,
                                   prop_explained = 1)
  dim(LD_embedding)
  
  PC_embedding = mapPCA(referenceSCE = referenceSCE,
                                   querySCE = querySCE,
                                   genes = genes,
                                   assayNameReference = assayNameReference,
                                   assayNameQuery = assayNameQuery)
  dim(PC_embedding)
  
  #stable continuous
  SC_embedding = stabMapContinuous(referenceSCE = referenceSCE,
                                   querySCE = querySCE,
                                   genes = genes,
                                   assayNameReference = assayNameReference,
                                   assayNameQuery = assayNameQuery,
                                   sparse = FALSE)
  dim(SC_embedding)
  
  
  if (doUMAP) {
    LD_embedding_UMAP = calculateUMAP(t(LD_embedding))
    rownames(LD_embedding_UMAP) <- rownames(LD_embedding)
    PC_embedding_UMAP = calculateUMAP(t(PC_embedding))
    rownames(PC_embedding_UMAP) <- rownames(PC_embedding)
    SC_embedding_UMAP = calculateUMAP(t(SC_embedding))
    rownames(SC_embedding_UMAP) <- rownames(SC_embedding)
    
    
    # graph uncorrected
    
    par(mfrow = c(4,3))
    plot(rev(PC_embedding_UMAP[,1]), rev(PC_embedding_UMAP[,2]),
         col = type[rev(rownames(PC_embedding_UMAP))],
         main = "PCA uncorrected")
    plot(rev(LD_embedding_UMAP[,1]), rev(LD_embedding_UMAP[,2]),
         col = type[rev(rownames(LD_embedding_UMAP))],
         main = "LDA uncorrected")
    plot(rev(SC_embedding_UMAP[,1]), rev(SC_embedding_UMAP[,2]),
         col = type[rev(rownames(SC_embedding_UMAP))],
         main = "SC uncorrected")
    
    ct <- ct_all[rownames(PC_embedding_UMAP)]
    
    plot(rev(PC_embedding_UMAP[,1]), rev(PC_embedding_UMAP[,2]),
         col = rev(celltype_colours[ct]), main = "PCA uncorrected", pch = 16)
    plot(rev(LD_embedding_UMAP[,1]), rev(LD_embedding_UMAP[,2]),
         col = rev(celltype_colours[ct]), main = "LDA uncorrected", pch = 16)
    plot(rev(SC_embedding_UMAP[,1]), rev(SC_embedding_UMAP[,2]),
         col = rev(celltype_colours[ct]), main = "SC uncorrected", pch = 16)
    
  }
  
  # batch correct using reducedMNN
  
  batchFactor_unordered = rep(c("Reference", "Query"),
                              times = c(ncol(referenceSCE), ncol(querySCE)))
  names(batchFactor_unordered) <- c(colnames(referenceSCE), colnames(querySCE))
  batchFactor = batchFactor_unordered[rownames(LD_embedding)]
  # batchFactor = factor(substring(rownames(LD_embedding),1,1))
  # names(batchFactor) <- rownames(LD_embedding)
  
  LD_embedding_corrected = reducedMNN_batchFactor(LD_embedding,
                                                batchFactor)
  PC_embedding_corrected = reducedMNN_batchFactor(PC_embedding,
                                                batchFactor)
  SC_embedding_corrected = reducedMNN_batchFactor(as.matrix(SC_embedding),
                                                batchFactor)
  
  if (doUMAP) {
    LD_embedding_corrected_UMAP = calculateUMAP(t(LD_embedding_corrected))
    rownames(LD_embedding_corrected_UMAP) <- rownames(LD_embedding_corrected)
    PC_embedding_corrected_UMAP = calculateUMAP(t(PC_embedding_corrected))
    rownames(PC_embedding_corrected_UMAP) <- rownames(PC_embedding_corrected)
    SC_embedding_corrected_UMAP = calculateUMAP(t(SC_embedding_corrected))
    rownames(SC_embedding_corrected_UMAP) <- rownames(SC_embedding_corrected)
    
    # graph
    
    # par(mfrow = c(2,2))
    plot(rev(PC_embedding_corrected_UMAP[,1]), rev(PC_embedding_corrected_UMAP[,2]),
         col = type[rev(rownames(PC_embedding_corrected_UMAP))],
         main = "PCA corrected")
    plot(rev(LD_embedding_corrected_UMAP[,1]), rev(LD_embedding_corrected_UMAP[,2]),
         col = type[rev(rownames(LD_embedding_corrected_UMAP))],
         main = "LDA corrected")
    plot(rev(SC_embedding_corrected_UMAP[,1]), rev(SC_embedding_corrected_UMAP[,2]),
         col = type[rev(rownames(SC_embedding_corrected_UMAP))],
         main = "SC corrected")
    
    ct <- ct_all[rownames(PC_embedding_corrected_UMAP)]
    
    plot(rev(PC_embedding_corrected_UMAP[,1]), rev(PC_embedding_corrected_UMAP[,2]),
         col = rev(celltype_colours[ct]), main = "PCA corrected", pch = 16)
    plot(rev(LD_embedding_corrected_UMAP[,1]), rev(LD_embedding_corrected_UMAP[,2]),
         col = rev(celltype_colours[ct]), main = "LDA corrected", pch = 16)
    plot(rev(SC_embedding_corrected_UMAP[,1]), rev(SC_embedding_corrected_UMAP[,2]),
         col = rev(celltype_colours[ct]), main = "SC corrected", pch = 16)
    
    
  }
  
  # predict cell type labels using knn with k = 5
  referenceLabels = colData(referenceSCE)[,grouping]
  names(referenceLabels) = colnames(referenceSCE)
  
  knn_LD = embeddingKNN(LD_embedding_corrected,
                        referenceLabels,
                        type = "uniform_fixed",
                        k_values = 5)
  
  knn_PC = embeddingKNN(PC_embedding_corrected,
                        referenceLabels,
                        type = "uniform_fixed",
                        k_values = 5)
  
  knn_SC = embeddingKNN(SC_embedding_corrected,
                        referenceLabels,
                        type = "uniform_fixed",
                        k_values = 5)
  
  
  queryLabels = colData(querySCE)[,grouping]
  names(queryLabels) = colnames(querySCE)
  
  mean(isEqual(knn_PC[names(queryLabels),"predicted_labels"], queryLabels))
  mean(isEqual(knn_SC[names(queryLabels),"predicted_labels"], queryLabels))
  mean(isEqual(knn_LD[names(queryLabels),"predicted_labels"], queryLabels))
  
  res = rbind(res,
              data.frame(genes = length(genes),
                         type = c("PC", "SC", "LD"),
                         Accuracy = c(mean(isEqual(knn_PC[names(queryLabels),"predicted_labels"], queryLabels)),
                                      mean(isEqual(knn_SC[names(queryLabels),"predicted_labels"], queryLabels)),
                                      mean(isEqual(knn_LD[names(queryLabels),"predicted_labels"], queryLabels)))
              ))
  
  saveRDS(res, file = "res.Rds")
  
  
  g = ggplot(res, aes(x = genes, y = Accuracy)) + 
    theme_classic() +
    geom_point(aes(colour = type)) + 
    geom_smooth(aes(group = type, colour = type), fill = NA, method = "loess") + 
    NULL
  print(g)
  
}







if (FALSE) {
# predict spatial labels using SVM
LD_svm_pred = getSVMPrediction(atlas_sub, sce_filt, LD_embedding_corrected)

sort(table(LD_svm_pred))


### wrapper function

LDA_classify = function(atlas_sub, sce_filt, batchFactor) {
  
  atlas_sub <- add_quickSubClustering(atlas_sub, min.ncells = 100)
  
  print("add_quickSubClustering done")
  
  LD_embedding = stabMapLabelled(atlas_sub, sce_filt)
  
  print("stabMapLabelled done")
  
  LD_embedding_corrected = reducedMNN_batchFactor(LD_embedding, batchFactor)
  
  print("reducedMNN_batchFactor done")
  
  LD_svm_pred = getSVMPrediction(atlas_sub, sce_filt, LD_embedding_corrected)
  
  print("getSVMPrediction done")
  
  return(LD_svm_pred)
}


LD_svm_pred = LDA_classify(atlas_sub, sce_filt, batchFactor)
}


# wt chimera



if (FALSE) {
  # wt chimera comparison
  #genes = intersect(intersect(sample(rownames(wt))[1:500], names(which(rowSums(counts(wt))>0))), rownames(sce_filt))
  genes = intersect(sample(rownames(wt))[1:500], names(which(rowSums(counts(wt))>0)))
  length(genes)
  #genes = rownames(wt_sub)
  
  LD_embedding = stabMapLabelled(wt, NULL,
                                   genes = genes,
                                   grouping = "celltype.mapped",
                                   assayNameReference = "logcounts",
                                   assayNameQuery = "logcounts",
                                   prop_explained = 1)
  PC_embedding = mapPCA(wt, NULL,
                                   genes = genes,
                                   grouping = "celltype.mapped",
                                   assayNameReference = "logcounts",
                                   assayNameQuery = "logcounts",
                                   prop_explained = 0.99)
  
  LD_embedding_UMAP = calculateUMAP(t(LD_embedding))
  rownames(LD_embedding_UMAP) <- rownames(LD_embedding)
  PC_embedding_UMAP = calculateUMAP(t(PC_embedding))
  rownames(PC_embedding_UMAP) <- rownames(PC_embedding)
  
  par(mfrow = c(1,2))
  plot(PC_embedding_UMAP[colnames(wt),1], PC_embedding_UMAP[colnames(wt),2], col = celltype_colours[as.character(wt$celltype.mapped)], pch = 16, main = "PCA")
  plot(LD_embedding_UMAP[colnames(wt),1], LD_embedding_UMAP[colnames(wt),2], col = celltype_colours[as.character(wt$celltype.mapped)], pch = 16, main = "LDA")
}



# examine stabMapContinuous

emb = stabMapContinuous(referenceSCE,
                        querySCE,
                        genes = genes,
                        assayNameReference = "logcounts",
                        assayNameQuery = "logcounts",
                        ncomponentsFull = 50,
                        ncomponentsSubset = 50)

emb_umap = calculateUMAP(t(emb))
plot(emb_umap, col = celltype_colours[ct[rownames(emb)]])


# build some motivating data for stabMapComparative

nGenes_X = 10000
nGenes_Y = 5000
nGenes_xy = 100
genes_X = sample(rownames(atlas))[seq_len(nGenes_X)]
genes_Y = c(genes_X[seq_len(nGenes_xy)], sample(setdiff(rownames(atlas), genes_X))[seq_len(nGenes_Y - nGenes_xy)])

length(genes_X); length(genes_Y); length(intersect(genes_X, genes_Y))

SCE_X = atlas[genes_X,atlas$sample == 17]
SCE_Y = atlas[genes_Y,atlas$sample == 29]
assayNameX = "logcounts"
assayNameY = "logcounts"

emb = stabMapComparative(SCE_X,
                         SCE_Y,
                         genes = intersect(rownames(SCE_X), rownames(SCE_Y)),
                         assayNameX = assayNameX,
                         assayNameY = assayNameY,
                         ncomponentsFull = 50,
                         ncomponentsSubset = 50,
                         sparse = FALSE)
emb_UMAP = calculateUMAP(t(emb))
rownames(emb_UMAP) <- rownames(emb)

emb_PC = mapPCA(referenceSCE = SCE_X,
                           querySCE = SCE_Y,
                           genes = intersect(rownames(SCE_X), rownames(SCE_Y)),
                           assayNameReference = assayNameX,
                           assayNameQuery = assayNameY)
emb_PC_UMAP = calculateUMAP(t(emb_PC))
rownames(emb_PC_UMAP) <- rownames(emb_PC)


# visual comparison
par(mfrow = c(1,2))
plot(emb_UMAP, col = celltype_colours[ct_all[rownames(emb_UMAP)]], main = "stabMap")
plot(emb_PC_UMAP, col = celltype_colours[ct_all[rownames(emb_PC_UMAP)]], main = "PCA")


# predict cell type labels using knn with k = 5
LabelsX = colData(SCE_X)[,grouping]
names(LabelsX) = colnames(SCE_X)

LabelsY = colData(SCE_Y)[,grouping]
names(LabelsY) = colnames(SCE_Y)

knn_SC_X = embeddingKNN(emb,
                        LabelsX,
                        type = "uniform_fixed",
                        k_values = 5)

knn_SC_Y = embeddingKNN(emb,
                        LabelsY,
                        type = "uniform_fixed",
                        k_values = 5)

knn_PC_X = embeddingKNN(emb_PC,
                        LabelsX,
                      type = "uniform_fixed",
                      k_values = 5)

knn_PC_Y = embeddingKNN(emb_PC,
                        LabelsY,
                      type = "uniform_fixed",
                      k_values = 5)


# calculate accuracy for cell type prediction as a quantitative readout
mean(c(isEqual(knn_SC_X[names(LabelsY),"predicted_labels"], LabelsY),
       isEqual(knn_SC_Y[names(LabelsX),"predicted_labels"], LabelsX)))
mean(c(isEqual(knn_PC_X[names(LabelsY),"predicted_labels"], LabelsY),
       isEqual(knn_PC_Y[names(LabelsX),"predicted_labels"], LabelsX)))
