library(edgeR)
library(VennDiagram)
library(pheatmap)
library(GenomicFeatures)
library(corrplot)
library(gridExtra)

#functions####
###heatmap function
heatmap <- function(dataframe1, dataframe2, dataframe3, title) {
  all_DEGs <- Reduce(union, list(
    rownames(subset(dataframe1, FDR < 0.05)),
    rownames(subset(dataframe2, FDR < 0.05)),
    rownames(subset(dataframe3, FDR < 0.05))
  ))
  pheatmap( logcpm_dge[all_DEGs,] , scale = "row" , cluster_cols = F, color = colorRampPalette(c("navy","white","firebrick"))(90), show_rownames=F, main = title, border_color = "black")
}

#volcano plot
get_upregulated_and_downregulated <- function(FDR, logFC) { 
  if (FDR < 0.05 && logFC > 0) { 
    return("UP")
  } else if (FDR < 0.05 && logFC < 0) {
    return("DOWN")
  } else {
    return("NO")  
  }
}

#PARENTAL####
#Feature count from reads mapped against RPE-1 haplotype1
count_data_hap1 <- read.table("rpe1v1.1_hap1_feature_counts.txt", header=TRUE, row.names=1, sep="\t")
count_data_hap1_f <- count_data_hap1[, c(43,42,41)]
fpkms_hap1 = rpkm(count_data_hap1_f, gene.length=count_data_hap1$Length)

#Feature count from reads mapped against RPE-1 haplotype2
count_data_hap2 <- read.table("rpe1v1.1_hap2_feature_counts.txt", header=TRUE, row.names=1, sep="\t")
count_data_hap2_f <- count_data_hap2[, c(43,42,41)]
fpkms_hap2 = rpkm(count_data_hap2_f, gene.length=count_data_hap2$Length)

#Feature count from reads mapped against hg38
count_data_hg38 <- read.table("hg38_feature_counts.txt", header=TRUE, row.names=1, sep="\t")
count_data_hg38_f <- count_data_hg38[, c(43,42,41)]
fpkms_hg38 = rpkm(count_data_hg38_f, gene.length=count_data_hg38$Length)

fpkms_finale <- as.data.frame(cbind(fpkms_hap1, fpkms_hap2, fpkms_hg38))
keep_fpkms <- rowSums(fpkms_finale > 1) >= 4
fpkms_finale <- fpkms_finale[keep_fpkms,]

fpkm.table <- data.frame(A1=rowMeans(fpkms_finale[,1:3]), 
                         A2=rowMeans(fpkms_finale[,4:6]),
                         A3=rowMeans(fpkms_finale[,7:9]))

rownames(fpkm.table) = rownames(fpkms_finale)

count_data_hap1_f <- count_data_hap1_f[rownames(fpkms_finale),]
colnames(count_data_hap1_f) <- c( "parent_1_hap1", "parent_2_hap1", "parent_3_hap1")

count_data_hap2_f <- count_data_hap2_f[rownames(fpkms_finale),]
colnames(count_data_hap2_f) <- c( "parent_1_hap2", "parent_2_hap2", "parent_3_hap2")

count_data_hg38_f <- count_data_hg38_f[rownames(fpkms_finale),]
colnames(count_data_hg38_f) <- c( "parent_1_hg38", "parent_2_hg38", "parent_3_hg38")

counts_final <- cbind(count_data_hap1_f, count_data_hap2_f, count_data_hg38_f)

condition_final <- factor(c("B1","B1","B1",
                            "B2","B2","B2",
                            "B3","B3","B3"))

dge_final = DGEList(counts=counts_final,group=condition_final)
cpms_dge = cpm(dge_final)
logcpm_dge <- log2(cpms_dge + 1)
#cpms = cpm(counts_final + 1)
#logcpm_f <- log2(cpms)
norm_counts <- cpm(dge_final, log = TRUE)
#heatmap to show all the logCPM of all genes
pheatmap(logcpm_dge, cluster_cols = F, color = colorRampPalette(c("navy","white","firebrick"))(90), show_rownames=F, main = "logCPM All Genes")

dge_final <- calcNormFactors(dge_final)

plotMDS(dge_final) ####standard method (based on the fold change)
plotMDS(dge_final, method="bcv") ####biological coefficient of variation method (based on the bcv <- more able to separate samples)

M = cor(logcpm_dge) ###Create a correlation matrix which calculate the pearson correlation coefficient for each pair of samples based on the expression
corrplot(M, method = "number") ####create a correlation plot 

####Based on this correlation matrix you can calculate a distance matrix
####I subtract this correlation values to 1, lower the value, higher the similarity
####I convert it in to a distance object using as.dist function
distCor <- as.dist(1-M) 

####Hierarchical clustering of samples based on pearson correlation and the distance between samples 
hc <- hclust(distCor)
plot(hc)

design <- model.matrix(~0+condition_final, data=dge_final$samples)
rownames(design) <- colnames(counts_final) ####assign the column names of the count
design

dge_final_final <- estimateDisp(dge_final, design)
fit <- glmFit(dge_final_final, design)

my.contrasts <- makeContrasts(B1vsB2=condition_finalB1-condition_finalB2,
                              B1vsB3=condition_finalB1-condition_finalB3,
                              B2vsB3=condition_finalB2-condition_finalB3,
                              levels=design)

lrt.B1vsB2 <- glmLRT(fit, contrast=my.contrasts[,"B1vsB2"]) #HAP1-HAP2 left-right
lrt.B1vsB3 <- glmLRT(fit, contrast=my.contrasts[,"B1vsB3"]) #HAP1-HG38
lrt.B2vsB3 <- glmLRT(fit, contrast=my.contrasts[,"B2vsB3"]) #HAP2-HG38

tp_lrt.B1vsB2 = topTags(lrt.B1vsB2, n=nrow(lrt.B1vsB2$table))
tp_lrt.B1vsB3 = topTags(lrt.B1vsB3, n=nrow(lrt.B1vsB3$table))
tp_lrt.B2vsB3 = topTags(lrt.B2vsB3, n=nrow(lrt.B2vsB3$table))

tp_lrt.B1vsB2_frame = as.data.frame(tp_lrt.B1vsB2)
tp_lrt.B1vsB3_frame = as.data.frame(tp_lrt.B1vsB3)
tp_lrt.B2vsB3_frame = as.data.frame(tp_lrt.B2vsB3)

tp_lrt.B1vsB2_frame$diffexpressed <- mapply(get_upregulated_and_downregulated, tp_lrt.B1vsB2_frame$FDR, tp_lrt.B1vsB2_frame$logFC)
tp_lrt.B1vsB2_frame$genes <- rownames(tp_lrt.B1vsB2_frame)
tp_lrt.B1vsB3_frame$diffexpressed <- mapply(get_upregulated_and_downregulated, tp_lrt.B1vsB3_frame$FDR, tp_lrt.B1vsB3_frame$logFC)
tp_lrt.B1vsB3_frame$genes <- rownames(tp_lrt.B1vsB3_frame)
tp_lrt.B2vsB3_frame$diffexpressed <- mapply(get_upregulated_and_downregulated, tp_lrt.B2vsB3_frame$FDR, tp_lrt.B2vsB3_frame$logFC)
tp_lrt.B2vsB3_frame$genes <- rownames(tp_lrt.B2vsB3_frame)

#Volcano
#you do p1 p2 and p3 by changing the df and then grid.arrange to merge the three volcano plots
p1 <- ggplot(data = tp_lrt.B1vsB2_frame, aes(x = logFC, y = -log10(PValue), col = diffexpressed)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 200), xlim = c(-15, 15)) +
  ggtitle("Parental Hap1 vs. Hap2") +
  theme_bw() +
  theme(
    plot.margin = margin(t = 15, r = 2, b = 15, l = 2), legend.position = "none") +
  scale_x_continuous(breaks = seq(-15, 15, 5)) 

p2 <- ggplot(data = tp_lrt.B1vsB3_frame, aes(x = logFC, y = -log10(PValue), col = diffexpressed)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 200), xlim = c(-15, 15)) +
  ggtitle("Parental Hap1 vs. HG38") +
  theme_bw() +
  theme(
    plot.margin = margin(t = 15, r = 2, b = 15, l = 2), legend.position = "none") +
  scale_x_continuous(breaks = seq(-15, 15, 5)) 

p3 <- ggplot(data = tp_lrt.B2vsB3_frame, aes(x = logFC, y = -log10(PValue), col = diffexpressed)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 200), xlim = c(-15, 15)) +
  ggtitle("Parental Hap2 vs. HG38") +
  theme_bw() +
  theme(
    plot.margin = margin(t = 15, r = 2, b = 15, l = 2), legend.position = "none") +
  scale_x_continuous(breaks = seq(-15, 15, 5)) 
grid.arrange(p1, p2, p3, nrow = 1)

#Heatmap function to show the locpm of only the DGEs in the three comparisons
heatmap(tp_lrt.B1vsB2_frame, tp_lrt.B1vsB3_frame, tp_lrt.B2vsB3_frame, "HAP1vsHAP2, HAP1vsHG38, HAP2vsHG38")

#Single heatmap of each DEGs in each comparison and creating files for venn diagram
all_DEGs_B1vsB2 <- rownames(subset(tp_lrt.B1vsB2_frame, FDR<.05 ))
pheatmap(logcpm_dge[all_DEGs_B1vsB2,], scale = "row" , color = colorRampPalette(c("navy","white","firebrick"))(90), show_rownames=F, main = "Hap1 vs Hap2 DEGs")
all_DEGs_B1vsB3 <- rownames(subset(tp_lrt.B1vsB3_frame, FDR<.05 ))
pheatmap(logcpm_dge[all_DEGs_B1vsB3,], scale = "row" , color = colorRampPalette(c("navy","white","firebrick"))(90), show_rownames=F, main = "Hap1 vs HG38 DEGs")
all_DEGs_B2vsB3 <- rownames(subset(tp_lrt.B2vsB3_frame, FDR<.05 ))
pheatmap(logcpm_dge[all_DEGs_B2vsB3,], scale = "row" , color = colorRampPalette(c("navy","white","firebrick"))(90), show_rownames=F, main = "Hap2 vs HG38 DEGs")

#Creating files for FDR plot
all_DEGs_hap1_hap2_parental_unique <- setdiff(all_DEGs_B1vsB2, union(all_DEGs_B2vsB3, all_DEGs_B1vsB3))
all_DEGs_hap1_hg38_parental_unique <- setdiff(all_DEGs_B1vsB3, union(all_DEGs_B2vsB3, all_DEGs_B1vsB2))
all_DEGs_hap2_hg38_parental_unique <- setdiff(all_DEGs_B2vsB3, union(all_DEGs_B1vsB3, all_DEGs_B1vsB2))

all_unique <- c(all_DEGs_hap1_hap2_parental_unique,all_DEGs_hap1_hg38_parental_unique,all_DEGs_hap2_hg38_parental_unique)

B1vsB2 <- tp_lrt.B1vsB2_frame[all_unique,]
B1vsB3 <- tp_lrt.B1vsB3_frame[all_unique,]
B2vsB3 <- tp_lrt.B2vsB3_frame[all_unique,]

B1vsB2$genome <- "Parental_HAP1vsHAP2"
B1vsB2$gene <- rownames(B1vsB2)
B1vsB3$genome <- "Parental_HAP1vsHG38"
B1vsB3$gene <- rownames(B1vsB3)
B2vsB3$genome <- "Parental_HAP2vsHG38"
B2vsB3$gene <- rownames(B2vsB3)

FDR_unique <- rbind(B1vsB2, B1vsB3, B2vsB3)
FDR_unique$gene <- gsub("gene:", "", FDR_unique$gene)