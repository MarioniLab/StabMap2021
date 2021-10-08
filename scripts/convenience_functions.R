calculateUMAP_rnames = function(embedding, ...) {
  # pass to calculateUMAP from scater
  # but keep rownames
  require(scater)
  message("calculating UMAP...")
  embedding_UMAP = calculateUMAP(t(embedding), ...)
  rownames(embedding_UMAP) <- rownames(embedding)
  return(embedding_UMAP)
}