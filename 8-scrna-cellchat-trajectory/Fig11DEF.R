library(CellChat)
library(ggplot2)
library(dplyr)

pathways.target <- c("MIF", "SPP1", "CXCL", "VEGF", "FN1")
pathways.show <- intersect(pathways.target, cellchat_obj@netP$pathways)

pdf("Fig11D_Bubble_Fibro_Immune.pdf", width=7.5, height=9)
p <- netVisual_bubble(
  cellchat_obj,
  sources.use = "Fibroblasts",
  targets.use = c("T lymphocytes", "Myeloid cells", "NK cells"),
  signaling = pathways.show,
  remove.isolate = TRUE
)
p + theme_bw(base_size = 15) +
  theme(plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x      = element_text(size = 12, color = "black", angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 11, color = "black"),
    legend.text      = element_text(size = 14), legend.title     = element_text(size = 16, face = "bold"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border     = element_rect(color="black", fill=NA, linewidth=0.7),
    plot.margin      = margin(10, 10, 10, 10)
  )
dev.off()

##############  Fig11E  ###########
df_all <- data.frame()
for (pw in pw.list) {
  res <- netAnalysis_contribution(cellchat_obj, signaling = pw)
  df <- res$data
  df$pair <- df$name
  df$Pathway <- pw
  
  df_all <- rbind(df_all, df)
}
pdf("Fig11E_PathwayContribution.pdf", width = 7, height = 8)
ggplot(df_all, aes(x = reorder(pair, contribution), y = contribution)) +
  geom_col(fill = "#4DBBD5", width = 0.7, color = "black") +
  coord_flip() +
  facet_wrap(~ Pathway, ncol = 1, scales = "free_y") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "Ligand–receptor pair", y = "Contribution",
       title = "Fibroblast-Related Pathway Contribution"
  ) +
  theme_classic(base_size = 14) +
  theme(strip.text   = element_text(face = "bold", size = 18),
        axis.text.y  = element_text(size = 14, color = "black"),
        axis.text.x  = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        plot.title   = element_text(size = 20, face = "bold", hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA)
  )
dev.off()

##############  Fig11F  ##############

cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)
cellchat_obj <- computeCommunProbPathway(cellchat_obj)
cellchat_obj <- aggregateNet(cellchat_obj)
cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")
cellchat_obj@netP$pathways

netP_prob <- cellchat_obj@netP$prob
sending_strength    <- apply(netP_prob, 1, sum)
receiving_strength  <- apply(netP_prob, 2, sum)

df_role <- data.frame(
  CellGroup = levels(cellchat_obj@idents),
  Sending = sending_strength,
  Receiving = receiving_strength
)

pdf("Fig11F_SignalingRole.pdf", width=6, height=4.5)
ggplot(df_role, aes(x=CellGroup)) +
  geom_col(aes(y=Sending, fill="Sending")) +
  geom_col(aes(y=-Receiving, fill="Receiving")) +
  coord_flip() +
  scale_fill_manual(values=c("Sending"="#D55E00","Receiving"="#0072B2")) +
  labs(title="Signaling Role (sender vs receiver)", 
       y="Communication Strength", x="", fill="Role") +
  theme_classic(base_size=14) +
  theme(plot.title = element_text(hjust=0.5, face="bold"), 
        axis.text.y = element_text(size = 14, color = "black"))
dev.off()

