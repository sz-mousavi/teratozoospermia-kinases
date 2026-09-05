# Dysregulated kinases in teratozoospermia
# Microarray meta-analysis: GSE6967, GSE6968, GSE6872
# Paper: JRI 2025;26(4):224
# DOI: https://doi.org/10.18502/jri.v26i4.21088

# Put GEO series matrices and GPL annotation files in data/
# GSE6967_series_matrix.txt + GPL2507.txt
# GSE6968_series_matrix.txt + GPL2700.txt
# GSE6872_series_matrix.txt + GPL570.txt

library(data.table)
library(limma)
library(ggplot2)
library(sva)
library(pheatmap)
library(dplyr)

data_dir <- "data"
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

data1 <- read.delim(file.path(data_dir, "GSE6967_series_matrix.txt"))
data2 <- read.delim(file.path(data_dir, "GSE6968_series_matrix.txt"))
data3 <- read.delim(file.path(data_dir, "GSE6872_series_matrix.txt"))

data1 <- na.omit(data1)
data2 <- na.omit(data2)

data1[, -1] <- log2(data1[, -1] - min(data1[, -1]) + 1)
data2[, -1] <- log2(data2[, -1] - min(data2[, -1]) + 1)
data3[, -1] <- log2(1 + data3[, -1])

annot1 <- fread(file.path(data_dir, "GPL2507.txt"), data.table = FALSE)
annot2 <- fread(file.path(data_dir, "GPL2700.txt"), data.table = FALSE)
annot3 <- fread(file.path(data_dir, "GPL570.txt"), data.table = FALSE)

annot1$GB_ACC <- gsub("\\..*", "", annot1$GB_ACC)
annot2$GB_ACC <- gsub("\\..*", "", annot2$GB_ACC)

annot1 <- annot1[, c("ID", "GB_ACC")]
annot2 <- annot2[, c("ID", "GB_ACC")]
annot3 <- annot3[, c("ID", "GB_ACC")]

rownames(annot1) <- annot1$ID
rownames(annot2) <- annot2$ID
rownames(annot3) <- annot3$ID

data1$GB_ACC <- annot1[as.character(data1$ID_REF), "GB_ACC"]
data2$GB_ACC <- annot2[as.character(data2$ID_REF), "GB_ACC"]
data3$GB_ACC <- annot3[as.character(data3$ID_REF), "GB_ACC"]

data1 <- data1[, -1]
data2 <- data2[, -1]
data3 <- data3[, -1]

data1 <- aggregate(. ~ GB_ACC, data1, mean)
data2 <- aggregate(. ~ GB_ACC, data2, mean)
data3 <- aggregate(. ~ GB_ACC, data3, mean)

all1 <- merge(data1, data2, by = "GB_ACC")
all <- merge(all1, data3, by = "GB_ACC")
all <- subset(all, GB_ACC != "")

pdf(file.path(out_dir, "boxplot_merged.pdf"))
boxplot(all[, -1])
dev.off()

rownames(all) <- all$GB_ACC
all <- subset(all, select = -GB_ACC)
all <- data.matrix(all)

batch <- factor(c(
  rep(1, ncol(data1) - 1),
  rep(2, ncol(data2) - 1),
  rep(3, ncol(data3) - 1)
))

allm <- all - rowMeans(all)
pc <- prcomp(allm)
pcr <- data.frame(pc$rotation[, 1:3], batch)

pdf(file.path(out_dir, "PCAbatch.pdf"))
ggplot(pcr, aes(PC1, PC2, color = batch)) +
  geom_point() +
  theme_bw()
dev.off()

allc <- ComBat(allm, batch)
allc <- normalizeQuantiles(allc)

png(file.path(out_dir, "boxplot_afterComBat.png"), width = 720, height = 950)
boxplot(allc)
dev.off()

png(file.path(out_dir, "boxplotAfternormalization.png"), width = 720, height = 950)
boxplot(allc)
dev.off()

allm <- allc - rowMeans(allc)
pc <- prcomp(allm)
pcr <- data.frame(pc$rotation[, 1:3], batch)

pdf(file.path(out_dir, "PCAbatch_after_Normalization.pdf"))
ggplot(pcr, aes(PC1, PC2, color = batch)) +
  geom_point() +
  theme_bw()
dev.off()

allm <- as.data.frame(allm)
write.table(allm, file.path(out_dir, "allm.txt"))

# Sample order must match columns of the merged matrix
gr <- c(
  rep("terato", 5), rep("normo", 1),
  rep("terato", 3), rep("normo", 4),
  rep("terato", 6), rep("normo", 4),
  rep("normo", 13), rep("terato", 8)
)

pdf(file.path(out_dir, "boxplotallm.pdf"))
boxplot(allm)
dev.off()

exs <- t(scale(t(allm), center = TRUE, scale = FALSE))
pc <- prcomp(exs)
pcr <- data.frame(pc$rotation[, 1:2], Group = gr)

pdf(file.path(out_dir, "PCAallm.pdf"))
ggplot(pcr, aes(PC1, PC2, color = Group, label = colnames(exs))) +
  geom_point(size = 7) +
  geom_text(size = 0.8, color = "white") +
  theme_bw()
dev.off()

png(file.path(out_dir, "Heatmap_allm.png"), width = 720, height = 950)
pheatmap(
  cor(allm),
  labels_row = gr,
  labels_col = gr,
  border_color = NA,
  fontsize_row = 15,
  fontsize_col = 19
)
dev.off()

gr <- factor(gr)
design <- model.matrix(~ gr + 0)
colnames(design) <- levels(gr)
fit <- lmFit(allm, design)
cont.matrix <- makeContrasts(terato - normo, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit3 <- eBayes(fit2, 0.01)
tT <- topTable(fit3, adjust = "fdr", sort.by = "B", number = Inf)
tT <- subset(tT, select = c("GB_ACC", "adj.P.Val", "logFC"))

write.csv(tT, file.path(out_dir, "terato-normo.csv"))

updif <- subset(tT, logFC > 1 & adj.P.Val < 0.01)
downdif <- subset(tT, logFC < (-1) & adj.P.Val < 0.01)
alldif <- rbind(updif, downdif)

write.csv(updif, file.path(out_dir, "Up-DEGs1.csv"))
write.csv(downdif, file.path(out_dir, "Down-DEGs1.csv"))
write.csv(alldif, file.path(out_dir, "all-DEGs1.csv"))

allm_t <- t(allm)
write.csv(allm_t, file.path(out_dir, "transposed_allm.csv"))

allm_t <- as.data.frame(allm_t)
allm_t$id <- rownames(allm_t)
alldif$id <- rownames(alldif)
plotd <- subset(allm_t, id %in% alldif$id)
write.csv(plotd, file.path(out_dir, "DEGs-Matrix.csv"))
