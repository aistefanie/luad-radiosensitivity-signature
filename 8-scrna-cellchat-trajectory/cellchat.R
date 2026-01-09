# Figure11A、B: 细胞类型间的相互作用网络图 ####
seu_filtered <- subset(seu, subset = !is.na(Cell_type)) # 去掉 NA 分组细胞
cellchat_obj <- createCellChat(seu_filtered, group.by = "Cell_type.y")
cellchat_obj <- setIdent(cellchat_obj, ident.use = "Cell_type.y")
CellChatDB <- CellChatDB.human
cellchat_obj@DB <- CellChatDB
cellchat_obj <- subsetData(cellchat_obj)
cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
cellchat_obj <- computeCommunProb(cellchat_obj)
cellchat_obj <- aggregateNet(cellchat_obj)

# 网络图
groupSize <- as.numeric(table(cellchat_obj@idents))
pdf("CellChat_Net.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), xpd = TRUE)
p3 <- netVisual_circle(cellchat_obj@net$count, vertex.weight = groupSize,
                       weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
p4 <- netVisual_circle(cellchat_obj@net$weight, vertex.weight = groupSize,
                       weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights")
dev.off()
par(mfrow = c(1, 1), xpd = FALSE)
