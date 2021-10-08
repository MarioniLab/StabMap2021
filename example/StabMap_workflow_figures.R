# build a fake graph

library(igraph)

graph = graph.edgelist(rbind(c("D1", "D2"),
                             c("D1", "D3")),
                       directed = FALSE)

# graph = featureNetwork(assay_list)
V(graph)$size = 50
V(graph)$color = c("lightblue", "lightpink", "lightblue")
V(graph)$label.color = "black"

E(graph)$label <- c("F1", "F2")
E(graph)$weight = c(50, 20)

# pdf(paste0("../../Figures/raw/PBMC_nonoverlapping_multiome_", probval,
#            "_network.pdf"), 
#     width = 6, height = 6)
plot(graph, layout = layout.drl)
# dev.off()
# g
