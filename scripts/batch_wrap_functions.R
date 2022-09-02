reducedMNN_batchFactor = function(embedding,
                                  batchFactor,
                                  ...) {
  # batch correct within this embedding, wrapper around reducedMNN
  # batchFactor is a named vector that is matched
  # ... passed to reducedMNN
  
  require(batchelor)
  
  # for some reason reducedMNN only takes dense matrices
  embedding <- as.matrix(embedding)
  
  batchFactor_used = batchFactor[rownames(embedding)]
  
  out = reducedMNN(embedding, batch = batchFactor_used, ...)
  resub_corrected = out$corrected
  
  return(resub_corrected)
}

ComBat_batchFactor = function(embedding,
                              batchFactor,
                              ...) {
  # batch correct within this embedding, wrapper around sva's ComBat
  # batchFactor is a named vector that is matched
  # ... passed to reducedMNN
  
  require(sva)
  
  batchFactor_used = batchFactor[rownames(embedding)]
  
  out = ComBat(t(embedding), batch = batchFactor_used, ...)
  resub_corrected = t(out)
  
  return(resub_corrected)
}

Harmony_batchFactor = function(embedding,
                               batchFactor, ...) {
  # batch correct within this embedding, wrapper around sva's ComBat
  # batchFactor is a named vector that is matched
  # ... passed to reducedMNN
  
  require(harmony)
  
  batchFactor_used = batchFactor[rownames(embedding)]
  
  out = HarmonyMatrix(as.matrix(t(embedding)), batchFactor_used, do_pca = FALSE, ...)
  resub_corrected = t(out)
  
  return(resub_corrected)
}

batchCorrect_batchFactor = function(embedding,
                                    batchFactor,
                                    PARAM) {
  # batch correct within this embedding, wrapper around sva's ComBat
  # batchFactor is a named vector that is matched
  # PARAM passed to batchCorrect
  
  require(batchelor)
  
  batchFactor_used = batchFactor[rownames(embedding)]
  
  out = batchCorrect(t(embedding), batch = batchFactor_used, PARAM = PARAM)
  resub_corrected = t(assays(out)[[1]])
  
  return(resub_corrected)
}

Seurat_batchFactor = function(embedding,
                              batchFactor) {
  
  # Use the Seurat CCA approach to perform batch integration
  require(Seurat)
  
  seurat.list <- sapply(unique(batchFactor), function(x){
    CreateSeuratObject(t(embedding[names(batchFactor[batchFactor == x]),]), project = as.character(x))
  }, simplify = FALSE)
  
  # reduc.list = sapply(unique(batchFactor), function(x){
  #   CreateDimReducObject(embeddings = embedding[names(batchFactor[batchFactor == x]),])
  # }, simplify = FALSE)
  
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = Inf)
  })
  
  features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = Inf)
  
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features)
  # seurat.combined <- tryCatch(
  #   IntegrateData(anchorset = seurat.anchors),
  #   error = function(cond) IntegrateData(anchorset = seurat.anchors, k.weight = 10)
  # )
  # tryCatch(seurat.combined <- IntegrateData(anchorset = seurat.anchors, k.weight = 10)
  seurat.combined <- IntegrateData(anchorset = seurat.anchors, k.weight = 10)
  
  DefaultAssay(seurat.combined) <- "integrated"
  
  int = t(seurat.combined@assays$integrated@data)
  
  return(int)
}
