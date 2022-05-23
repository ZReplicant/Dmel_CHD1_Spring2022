# This is the example of pair-wise DEA

library(stringi)
library(GenomicRanges)
library(DESeq2)
library(data.table)
library(dplyr)
library(stringr)
library(gplots)
library(tximport)
library(openxlsx)
library(tidyverse)

setwd("/home/guest/Документы/rnaseq/main_experiment/dif/dif_male-mut_male-wt")

# Genes were taken from FlyBase, 6.41 release. It is possible to make the process faster with .gtf
# E.g. ftp://ftp.flybase.net/releases/FB2019_05/dmel_r6.30/gff/dmel-all-r6.30.gff.gz
genes <- fread("dmel-all-r6.41_alt.gff", nThread = 10, showProgress = TRUE, verbose = TRUE) %>% 
  dplyr::select(c(1, 4, 5, 7, 9)) %>% 
  setNames(c("chr", "start", "end", "strand", "attr")) %>% 
  mutate(id = sub('^ID=(FBgn[0-9]+);.*', '\\1', attr), chr = paste0("chr", chr), gene_name = sub('.*Name=([^;]+);.*', '\\1', attr), tss = ifelse(strand == "+", start, end)) %>% 
  dplyr::select(-attr)

save(genes, file = "./genes.RData")

genes.gr <- GRanges(
  seqnames = Rle(genes$chr),
  ranges = IRanges(
    start = genes$start,
    end = genes$end,
    names = genes$gene_name
  )
)

# Import Salmon output
files <- paste0("salmonout/",
                c("transcript_V300093791_L01_89", "transcript_V300093791_L01_90", "transcript_V300093791_L01_94"))
files <- file.path(files, "quant.genes.sf")
names(files) <- c("mut1", "mut2", "wt2")

txi.mut.wt <- tximport(files, "salmon", txOut = T, importer = read.delim)

head(txi.mut.wt$counts)

# DESEQ2
samples <- data.frame(row.names = names(files), conditions = c("mut", "mut", "wt"), orig = "evrogen")

ddsTxi <- DESeqDataSetFromTximport(txi.mut.wt, colData = samples, design =~ conditions)

keep <- rowSums(counts(ddsTxi)) >= 10

ddsTxi <- ddsTxi[keep,]
ddsTxi$conditions <- relevel(ddsTxi$conditions, ref = "wt")

dds <- DESeq(ddsTxi)
res <- results(dds)
# plotDispEsts(dds)

# These commands prevent my local RStudio abortion error
rm(genes.gr)
gc()


resultsNames(dds)

resLFC <- lfcShrink(dds, coef = "conditions_mut_vs_wt", type = "apeglm",
                    lfcThreshold = 1)
pdf("./plots/MA.plot.hp1.bothkd.pdf")
DESeq2::plotMA(resLFC, alpha = 0.05, ylim=c(-6, 6))
dev.off()

summary(res, alpha = 0.05)

res.df <- as.data.frame(res) %>% mutate(id = rownames(res))
res.df <- merge(genes, res.df, by = "id", all.y = T)
write.xlsx(res.df, file = "./tables/DESeq2.mut.wt.results.xlsx", dec = ",")
#res.df %>% filter(is.na(gene_name)) %>% View

# Filter DEGs to get min and max log2foldchanges for our samples
ranking <- filter(res.df, padj < 0.05, abs(log2FoldChange) >= 1)
ordered_ranking <- ranking[order(ranking[,log2FoldChange], decreasing = TRUE)]
write.xlsx(ordered_ranking, file = "./tables/ordered_ranking.xlsx", dec = ",")


res.df.tss.gr <- makeGRangesFromDataFrame(res.df %>% 
                                            mutate(start = tss,
                                                   end = tss) %>% 
                                            filter(-tss | padj < 0.05),#here instead of | may be &, but this variable isn't used anywhere anyway
                                          keep.extra.columns = T)

# get TPM values and filter out those that differ greatly between replicas (also no zeroes)
mut.wt.tpm <- as.data.frame(txi.mut.wt$abundance) %>% mutate(id = rownames(txi.mut.wt$abundance))
#                                                      ,wt.var = abs(WT1-WT2)/(WT1 + WT2),
#                                                      kd.var = abs(HP1KD1 - HP1KD2)/(HP1KD1 + HP1KD2)) %>% 
# filter(wt.var < 0.3, kd.var < 0.3)
# make df with average values for kd and wt
mut.wt.tpm.avg <- mut.wt.tpm %>% mutate(wt = wt2,
                                              mut = (mut1 + mut2)/2) %>% 
  dplyr::select(id, wt, mut)

boxplot(mut.wt.tpm.avg$mut, mut.wt.tpm.avg$wt, outline = F)

cor.mut.wt <- cor(mut.wt.tpm.avg$wt, mut.wt.tpm.avg$mut)
write.xlsx(cor.mut.wt, file = "correlation.xls", dec = ",")
# hp1.tpm.avg <- merge(genes, hp1.tpm.avg, by = "id")

# Pseudocounts for averaged TPM values
pco <- min(mut.wt.tpm.avg$mut[mut.wt.tpm.avg$mut > 0])
# hp1.tpm.avg <- merge(genes, hp1.tpm.avg, by = "id") %>%
#   mutate(WT = ifelse(WT == 0, pco, WT),
#                      KD = ifelse(KD == 0, pco, KD))

mut.wt.tpm.avg <- merge(genes, mut.wt.tpm.avg, by = "id") %>%
  mutate(wt = ifelse(wt == 0, pco, wt),
         mut = ifelse(mut == 0, pco, mut))

mut.wt.tpm.tss.gr <- makeGRangesFromDataFrame(mut.wt.tpm.avg %>% 
                                                   mutate(start = tss,
                                                          end = tss) %>% 
                                                   dplyr::select(-tss), keep.extra.columns = T)

mut.wt.tpm.gb.gr <- makeGRangesFromDataFrame(mut.wt.tpm.avg, keep.extra.columns = T)

# which(is.na(hp1.tpm.avg$tss.loc))

mut.wt.tpm.avg.reg <- mut.wt.tpm.avg %>% mutate(
  reg = NA
)
mut.wt.tpm.avg.reg$reg[mut.wt.tpm.avg.reg$id %in% (res.df %>%
                                                           filter(padj < 0.05, log2FoldChange > 0))$id &
                            log2(mut.wt.tpm.avg.reg$mut/mut.wt.tpm.avg.reg$wt) > 1] <- "up"

mut.wt.tpm.avg.reg$reg[mut.wt.tpm.avg.reg$id %in% (res.df %>%
                                                           filter(padj < 0.05, log2FoldChange < 0))$id &
                            log2(mut.wt.tpm.avg.reg$mut/mut.wt.tpm.avg.reg$wt) < -1] <- "down"
mut.wt.tpm.avg.reg$reg[is.na(mut.wt.tpm.avg.reg$reg)] <- "aaaaaaaa"
mut.wt.tpm.avg.reg$reg <- factor(mut.wt.tpm.avg.reg$reg)

table(mut.wt.tpm.avg.reg$reg)

pdf("./plots/TPM.based.mut.wt.difex.pdf")
ggplot(mut.wt.tpm.avg.reg %>% filter(wt > 0.1, mut > 0.1) %>% 
         arrange(reg), aes(x = log2(wt), y = log2(mut)))+
  geom_point(aes(col = reg))+
  stat_function(fun = function(x) x, col = "blue")+
  theme_bw()+
  theme(line = element_blank())+
  scale_color_manual(breaks = c("down", "up"),
                     values=c("#BBBBBB", "#0010DD", "#DD1010"),
                     name = "regulated")
dev.off()

save(cor.mut.wt, mut.wt.tpm.avg.reg, file = "mut.wt.tpm.avg.reg.RData")
