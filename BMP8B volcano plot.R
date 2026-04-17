if (!require(ggplot2)) install.packages("ggplot2")
if (!require(readxl)) install.packages("readxl")
if (!require(ggrepel)) install.packages("ggrepel")

library(ggplot2)
library(readxl)
library(ggrepel)

data <- read_excel("Klon1vsKontrolle_deg.xlsx", sheet = 1)

data$negLog10Padj <- -log10(data$padj)

data$Regulation <- "Not significant"
data$Regulation[data$padj < 0.05 & data$log2FoldChange > 1]  <- "Upregulated"
data$Regulation[data$padj < 0.05 & data$log2FoldChange < -1] <- "Downregulated"

data$Regulation <- factor(
  data$Regulation,
  levels = c("Upregulated", "Downregulated", "Not significant")
)
cols <- c(
  "Upregulated"     = "#D73027", 
  "Downregulated"   = "#4575B4", 
  "Not significant" = "black"
)


genes_to_label <- c(
)

label_data <- data[data$gene_name %in% genes_to_label, ]


volcano <- ggplot(
  data,
  aes(x = log2FoldChange, y = negLog10Padj)
) +
  geom_point(
    aes(color = Regulation),
    size = 1.4,
    alpha = 0.8
  ) +
  scale_color_manual(values = cols) +
  

  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4,
    color = "black"
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.4,
    color = "black"
  ) +
  
  geom_text_repel(
    data = label_data,
    aes(label = gene_name),
    size = 3.2,
    color = "black",
    box.padding = 0.5,
    point.padding = 0.4,
    segment.size = 0.3,
    max.overlaps = Inf
  ) +
  

  labs(
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~adjusted~italic(p))
  ) +
  

  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold")
  )

print(volcano)


ggsave(
  filename = "Volcano.png",
  plot = volcano,
  width = 6,
  height = 5,
  dpi = 600,
  bg = "white"
)
ggsave("Volcano.pdf", volcano, width = 6, height = 5)

