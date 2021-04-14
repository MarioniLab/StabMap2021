# wt is wild type chimera data
library(MouseGastrulationData)
library(SpatialUtils)
wt = WTChimeraData(type = "processed", samples = c("6", "8", "10"))
# wd "/Users/ghazan01/Dropbox/Backup/SpatialEmbryos/scripts"
atlas_sub = readRDS("../analysis_output/E8.5/atlas_sub.Rds")
sce_filt = readRDS("../analysis_output/E8.5/E8.5_sce_filt_unlabelled.Rds")
sce_filt
sce_filt <- reNormalise(sce_filt)
mapping_dt = updateMesenchymeLabel(readRDS("../analysis_output/E8.5/celltype_annotation_refined.Rds"))
colData(sce_filt) <- cbind(colData(sce_filt), mapping_dt[colnames(sce_filt), setdiff(colnames(mapping_dt),colnames(colData(sce_filt)))])


# get LDA embedding
get_LDA_embedding = function(sce, spatial,
                             genes = rownames(spatial),
                             grouping = "celltype_parsed_sub",
                             assayNameAtlas = "cosineNorm",
                             assayNameSpatial = "lm_cosineNorm",
                             prop_explained = 1) {
  require(MASS)
  
  fit = lda(t(assay(sce, assayNameAtlas)[genes,]),
            grouping = colData(sce)[,grouping])
  
  # proportion of between group variation explained by the linear
  # discriminants
  prop_rel = fit$svd^2/sum(fit$svd^2)
  n_prop_rel = which(cumsum(prop_rel) >= prop_explained)[1]
  
  if (is.null(spatial)) {
    all_assay = cbind(assay(sce, assayNameAtlas)[genes,])
  } else {
    all_assay = cbind(assay(sce, assayNameAtlas)[genes,],
                      assay(spatial, assayNameSpatial)[genes,])
  }
  
  resub = predict(object = fit, newdata = t(all_assay))
  
  if (prop_explained == 1) {
    resub_out = resub$x
  } else {
    resub_out = resub$x[,seq_len(n_prop_rel)]  
  }
  return(resub_out)
}

# get PCA embedding
get_PCA_embedding = function(sce, spatial,
                             genes = rownames(spatial),
                             grouping = "celltype_parsed_sub",
                             assayNameAtlas = "cosineNorm",
                             assayNameSpatial = "lm_cosineNorm",
                             prop_explained = 0.99) {
  require(irlba)
  
  if (is.null(spatial)) {
    all_assay = cbind(assay(sce, assayNameAtlas)[genes,])
  } else {
    all_assay = cbind(assay(sce, assayNameAtlas)[genes,],
                      assay(spatial, assayNameSpatial)[genes,])
  }
  
  pca_all <- irlba::prcomp_irlba(t(all_assay), n = 50, scale. = TRUE)$x
  rownames(pca_all) <- colnames(all_assay)
  
  return(pca_all)
}

genes = intersect(intersect(sample(rownames(wt))[1:500], names(which(rowSums(counts(wt))>0))), rownames(sce_filt))
length(genes)
genes = rownames(wt_sub)

LD_embedding = get_LDA_embedding(wt, NULL,
                                 genes = genes,
                                 grouping = "celltype.mapped",
                                 assayNameAtlas = "logcounts",
                                 assayNameSpatial = "logcounts",
                                 prop_explained = 1)
PC_embedding = get_PCA_embedding(wt, NULL,
                                 genes = genes,
                                 grouping = "celltype.mapped",
                                 assayNameAtlas = "logcounts",
                                 assayNameSpatial = "logcounts",
                                 prop_explained = 0.99)

LD_embedding_UMAP = calculateUMAP(t(LD_embedding))
rownames(LD_embedding_UMAP) <- rownames(LD_embedding)
PC_embedding_UMAP = calculateUMAP(t(PC_embedding))
rownames(PC_embedding_UMAP) <- rownames(PC_embedding)

par(mfrow = c(1,2))
plot(PC_embedding_UMAP[colnames(wt),1], PC_embedding_UMAP[colnames(wt),2], col = celltype_colours[as.character(wt$celltype.mapped)], pch = 16, main = "PCA")
plot(LD_embedding_UMAP[colnames(wt),1], LD_embedding_UMAP[colnames(wt),2], col = celltype_colours[as.character(wt$celltype.mapped)], pch = 16, main = "LDA")



LD_embedding = get_LDA_embedding(atlas_sub, sce_filt,
                                 genes = rownames(sce_filt),
                                 grouping = "celltype_parsed",
                                 assayNameAtlas = "logcounts",
                                 assayNameSpatial = "logcounts",
                                 prop_explained = 1)
dim(LD_embedding)

PC_embedding = get_PCA_embedding(atlas_sub, sce_filt,
                                 genes = rownames(sce_filt),
                                 grouping = "celltype_parsed",
                                 assayNameAtlas = "logcounts",
                                 assayNameSpatial = "logcounts",
                                 prop_explained = 1)
dim(PC_embedding)

LD_embedding_UMAP = calculateUMAP(t(LD_embedding))
rownames(LD_embedding_UMAP) <- rownames(LD_embedding)
PC_embedding_UMAP = calculateUMAP(t(PC_embedding))
rownames(PC_embedding_UMAP) <- rownames(PC_embedding)

par(mfrow = c(1,2))
plot(rev(PC_embedding_UMAP[,1]), rev(PC_embedding_UMAP[,2]),
     col = rev(factor(substring(rownames(PC_embedding_UMAP),1,1))), main = "PCA")
plot(rev(LD_embedding_UMAP[,1]), rev(LD_embedding_UMAP[,2]),
     col = rev(factor(substring(rownames(LD_embedding_UMAP),1,1))), main = "LDA")

ct = c(as.character(sce_filt$celltype_mapped_refined), as.character(atlas_sub$celltype_parsed))
names(ct) <- c(colnames(sce_filt), colnames(atlas_sub))
ct <- ct[rownames(PC_embedding_UMAP)]
par(mfrow = c(1,2))
plot(rev(PC_embedding_UMAP[,1]), rev(PC_embedding_UMAP[,2]),
     col = rev(celltype_colours[ct]), main = "PCA", pch = 16)
plot(rev(LD_embedding_UMAP[,1]), rev(LD_embedding_UMAP[,2]),
     col = rev(celltype_colours[ct]), main = "LDA", pch = 16)


# batch correct within this embedding
get_LD_MNN_corrected = function(LD_embedding,
                                batchFactor) {
  # batchFactor is a named vector that is matched
  
  require(batchelor)
  
  batchFactor_used = batchFactor[rownames(LD_embedding)]
  
  out = reducedMNN(LD_embedding, batch = batchFactor_used)
  resub_corrected = out$corrected
  
  return(resub_corrected)
}

batchFactor = factor(substring(rownames(LD_embedding),1,1))
names(batchFactor) <- rownames(LD_embedding)

LD_embedding_corrected = get_LD_MNN_corrected(LD_embedding,
                                              batchFactor)
PC_embedding_corrected = get_LD_MNN_corrected(PC_embedding,
                                              batchFactor)

LD_embedding_corrected_UMAP = calculateUMAP(t(LD_embedding_corrected))
rownames(LD_embedding_corrected_UMAP) <- rownames(LD_embedding_corrected)
PC_embedding_corrected_UMAP = calculateUMAP(t(PC_embedding_corrected))
rownames(PC_embedding_corrected_UMAP) <- rownames(PC_embedding_corrected)

par(mfrow = c(1,2))
plot(rev(PC_embedding_corrected_UMAP[,1]), rev(PC_embedding_corrected_UMAP[,2]),
     col = rev(factor(substring(rownames(PC_embedding_corrected_UMAP),1,1))), main = "PCA")
plot(rev(LD_embedding_corrected_UMAP[,1]), rev(LD_embedding_corrected_UMAP[,2]),
     col = rev(factor(substring(rownames(LD_embedding_corrected_UMAP),1,1))), main = "LDA")

ct = c(as.character(sce_filt$celltype_mapped_refined), as.character(atlas_sub$celltype_parsed))
names(ct) <- c(colnames(sce_filt), colnames(atlas_sub))
ct <- ct[rownames(PC_embedding_corrected_UMAP)]
par(mfrow = c(1,2))
plot(rev(PC_embedding_corrected_UMAP[,1]), rev(PC_embedding_corrected_UMAP[,2]),
     col = rev(celltype_colours[ct]), main = "PCA", pch = 16)
plot(rev(LD_embedding_corrected_UMAP[,1]), rev(LD_embedding_corrected_UMAP[,2]),
     col = rev(celltype_colours[ct]), main = "LDA", pch = 16)


# predict spatial labels using SVM
LD_svm_pred = getSVMPrediction(atlas_sub, sce_filt, LD_embedding_corrected)

sort(table(LD_svm_pred))




### wrapper function

LDA_classify = function(atlas_sub, sce_filt, batchFactor) {
  
  atlas_sub <- add_quickSubClustering(atlas_sub, min.ncells = 100)
  
  print("add_quickSubClustering done")
  
  LD_embedding = get_LDA_embedding(atlas_sub, sce_filt)
  
  print("get_LDA_embedding done")
  
  LD_embedding_corrected = get_LD_MNN_corrected(LD_embedding, batchFactor)
  
  print("get_LD_MNN_corrected done")
  
  LD_svm_pred = getSVMPrediction(atlas_sub, sce_filt, LD_embedding_corrected)
  
  print("getSVMPrediction done")
  
  return(LD_svm_pred)
}


LD_svm_pred = LDA_classify(atlas_sub, sce_filt, batchFactor)