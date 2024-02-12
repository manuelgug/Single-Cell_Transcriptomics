
library(edgeR)
library(ggplot2)
library(ggrepel)


### IMPORT DATA (downloaded from https://geo.metadataplus.biothings.io/geo/query/acc.cgi?acc=GSE151940;)

#counts
raw_counts <- read.csv("GSM4594094_M11_dcm.csv", header = F)
raw_counts <- as.data.frame(t(raw_counts))
dim(raw_counts)

sum(rowSums(raw_counts) == 0) # no empty rows. GOOD

#genes
genes <- read.csv("GSM4594094_M11_genes.csv", header=F)
dim(genes)

#conditions (in this case, barcodes were used to represent sample identity (ie: OD))
annot <- read.csv("GSM4594094_M11_barcodes.csv", header =F)
annot$V1 <- sub("^[^M]*_", "", annot$V1)
dim(annot)
unique(annot)


# use gene names as rownames for raw_counts
rownames(raw_counts) <- genes$V1


### ANALYZE DATA

# Create a DGEList object
dge <- DGEList(counts = raw_counts, group = factor(as.matrix(annot)), remove.zeros = T)

# Filter shitty samples
keep <- filterByExpr(y = dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Carculate normalization factors and normalize through TMM
dge <- calcNormFactors(object = dge)

# Estimate dispersions
dge <- estimateDisp(y = dge)

# Calculate average expression levels
ave_expr <- rowMeans(cpm(dge))

# Perform an exact test and adjust p-values thorugh FDR
et <- exactTest(object = dge)
top_degs = topTags(et, n = Inf)$table
top_degs

write.csv(top_degs, "top_degs.csv")


### VISUALIZE RESULTS

# Define significance cutoffs
logFC_cutoff <- 0.1    # Laughably lax, but whatever...
FDR_cutoff <- 0.05

# Define colors for upregulated and downregulated genes
top_degs$color <- ifelse(top_degs$logFC > logFC_cutoff & top_degs$FDR < FDR_cutoff, "red",
                         ifelse(top_degs$logFC < -logFC_cutoff & top_degs$FDR < FDR_cutoff, "cyan3", "gray"))

# Vvolcano plot BOOM
volcano_plot <- ggplot(top_degs, aes(x = logFC, y = -log10(PValue), color = color, label = rownames(top_degs))) +
  geom_point(size = 4) +
  geom_hline(yintercept = -log10(FDR_cutoff), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "grey") +
  geom_text_repel(size = 3, force = 10) +
  scale_color_identity() +
  labs(x = "log2 Fold Change", y = "-log10(PValue)", title = "Volcano Plot") +
  theme_minimal()

print(volcano_plot)

# MA plot
ma_plot <- ggplot(top_degs, aes(x = ave_expr, y = logFC, color = color, label = rownames(top_degs))) +
  geom_point() +
  geom_text_repel(size = 3, force = 10) +  # Add dot labels
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_identity() +
  labs(x = "Average Expression", y = "log2 Fold Change", title = "MA Plot") +
  theme_minimal()

print(ma_plot)


ggsave("volcano_plot.png", plot = volcano_plot, width = 8, height = 6, units = "in", bg = "white")
ggsave("ma_plot.png", plot = ma_plot, width = 8, height = 6, units = "in", bg = "white")
