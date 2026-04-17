if (!require(readxl)) install.packages("readxl")
if (!require(pheatmap)) install.packages("pheatmap")
if (!require(dplyr)) install.packages("dplyr")

library(readxl)
library(pheatmap)
library(dplyr)

data <- read_excel("Klon1vsKontrolle_deg.xlsx", sheet = 1)


data <- data %>% 
  filter(
    rowMeans(select(., K1_1, K1_2, K1_3), na.rm = TRUE) > 0 &
      rowMeans(select(., Ctr_1, Ctr_2, Ctr_3), na.rm = TRUE) > 0
  )

required_cols <- c("K1_1","K1_2","K1_3",
                   "Ctr_1","Ctr_2","Ctr_3",
                   "log2FoldChange","pvalue","padj","gene_name")

missing <- required_cols[!required_cols %in% colnames(data)]
if(length(missing) > 0){
  stop(paste("Fehlende Spalten:", paste(missing, collapse=", ")))
}

data <- data %>% 
  mutate(
    log2FoldChange = as.numeric(log2FoldChange),
    pvalue = as.numeric(pvalue),
    padj = as.numeric(padj),
    K1_1 = as.numeric(K1_1),
    K1_2 = as.numeric(K1_2),
    K1_3 = as.numeric(K1_3),
    Ctr_1 = as.numeric(Ctr_1),
    Ctr_2 = as.numeric(Ctr_2),
    Ctr_3 = as.numeric(Ctr_3)
  )

data <- data %>% 
  mutate(
    negLog10Padj = -log10(padj),
    volcano_score = abs(log2FoldChange) * negLog10Padj
  )


volcano_filtered <- data %>% 
  filter(
    !is.na(padj),
    abs(log2FoldChange) >= 1,   
    padj <= 0.05                
  )

top20_up <- volcano_filtered %>% 
  filter(log2FoldChange > 0) %>% 
  arrange(desc(volcano_score)) %>% 
  slice_head(n = 20)

top20_down <- volcano_filtered %>% 
  filter(log2FoldChange < 0) %>% 
  arrange(desc(volcano_score)) %>% 
  slice_head(n = 20)


top_genes <- bind_rows(top20_up, top20_down)

top_genes <- top_genes %>% 
  arrange(desc(log2FoldChange))

count_mat <- as.matrix(
  top_genes[, c("K1_1","K1_2","K1_3","Ctr_1","Ctr_2","Ctr_3")]
)
rownames(count_mat) <- top_genes$gene_name

count_mat_log <- log2(count_mat + 1)
count_mat_scaled <- t(scale(t(count_mat_log)))


annotation_col <- data.frame(
  Condition = c("Treatment","Treatment","Treatment",
                "Control","Control","Control")
)
rownames(annotation_col) <- colnames(count_mat_scaled)

ann_colors <- list(
  Condition = c(
    Treatment = "#B2182B",
    Control   = "#2166AC"
  )
)

pheatmap(
  count_mat_scaled,
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  cluster_rows = FALSE,   
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize_row = 9,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 10,
  border_color = NA,
  main = "Top 20 up- and downregulated genes (Volcano-based selection)"
)

png("Heatmap.png", width = 2000, height = 3500, res = 300)
pheatmap(
  count_mat_scaled,
  color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize_row = 9,
  fontsize_col = 10,
  cellwidth = 15,
  cellheight = 10,
  border_color = NA
)
dev.off()

