rm(list = ls())
library(ggplot2)
library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pheatmap)
library(stringr)
####38um

deseq38 = read.csv("D:\\gzy/38um_remove0_annotation_delmin_DESeq2.csv",header = T)
deseq50 = read.csv("D:\\gzy/50um_remove0_annotation_delmin_DESeq2.csv",header = T)
deseq75 = read.csv("D:\\gzy/75um_remove0_annotation_delmin_DESeq2.csv",header = T)

all38 = unique(c(deseq38$X,deseq50$X,deseq75$X))
entrezid_38all = bitr(all38, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrezid_38all = as.character(entrezid_38all[,2])

deseq38 = read.csv("D:\\gzy/38um_remove0_annotation_delmin_DESeq2.csv",header = T)


deseq38$Group <- "not-significant"
deseq38$Group[which(deseq38$padj <= 0.05 & deseq38$log2FoldChange >= 0.58)] = "up-regulate"
deseq38$Group[which(deseq38$padj <= 0.05 & deseq38$log2FoldChange <= (-0.58))] = "down-regulate"

up38 = deseq38$X[deseq38$padj <= 0.05 & deseq38$log2FoldChange >= 0.58]
down38 = deseq38$X[deseq38$padj <=0.05 & deseq38$log2FoldChange <= -0.58]
upvalue = length(up38)
downvalue = length(down38)

ma38 = ggplot(deseq38,aes(x = log2(baseMean),y = log2FoldChange ,colour = Group))+
  geom_point(alpha = 0.5, size=1)+
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  annotate(geom="text", x=15, y=6, label=upvalue,  fontface="bold",colour="#ff4757", size=5)+
  annotate(geom="text", x=15, y=-5, label=downvalue,  fontface="bold",colour="#546de5", size=5)+
  theme()+
  theme_set(theme_bw())+
  theme(panel.grid.major=element_line(colour=NA))
ggsave("D:\\gzy/ma38.tiff",ma38,width = 8)

#all38 = read.table("D:\\gzy/control_75.txt",header = T)
#all38 = all38[rowSums(all38[c(1,2,3),])>=3 | rowSums(all38[c(4,5,6),])>=3,]
#all38 = rownames(all38)

entrezid_38up = bitr(up38, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrezid_38up = as.character(entrezid_38up[,2])
go38up <- enrichGO(gene   = entrezid_38up,
                   OrgDb  = org.Hs.eg.db,
                   universe = entrezid_38all,
                   ont   = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable = T,
                   qvalueCutoff  = 0.05)
go38up = data.frame(go38up)
if (nrow(go38up)){
  go38up$termRatioStr =paste0(go38up$Count,'/', unlist(strsplit(go38up$BgRatio, '/'))[seq(1,nrow(go38up)*2,2)])
  go38up$termRatio = round(go38up$Count/as.numeric(unlist(strsplit(go38up$BgRatio, '/'))[seq(1,nrow(go38up)*2,2)]),2)
}
go38up = go38up[order(go38up$pvalue, decreasing=F),]
write.csv(go38up,"D:\\gzy/go38up.csv")
write.table(go38up,"D:\\gzy/goup.txt",quote = F)
if (nrow(go38up)>10){
  go38up = go38up[1:10,]
}
go38up$Description = factor(go38up$Description,levels = as.character(go38up$Description))
rownames(go38up) = go38up[,3]
re = rev(as.character(go38up$Description))
go38up = go38up[re,]

p_go38up <- ggplot(go38up, aes(x=-log10( pvalue ), y=Description)) +
  geom_point(aes(size= termRatio , colour = Count)  ) + 
  scale_y_discrete(limits=go38up$Description,labels = function(x) str_wrap(x, width = 40))+
  ggtitle("GO enrichment")  +  scale_color_gradient(low = 'blue', high = 'red') + #xlim(range(log10(out$pvalue))) +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )
  theme(axis.text.x=element_text(angle=0,size=10, vjust=0.7), axis.text.y=element_text(angle=0,size=12, vjust=0.7),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16),
        panel.background = element_rect(fill="white", colour='gray'),
        panel.grid.major = element_line(size = 0.05, colour = "gray"),
        panel.grid.minor.y = element_line(size=0.05, colour="gray"),
        panel.grid.minor.x = element_line(size=0.05, colour="gray"))
ggsave("D:\\gzy/go38up.tiff",p_go38up,width = 6.4)

kegg38up = enrichKEGG(entrezid_38up, universe = entrezid_38all)
kegg38up = data.frame(kegg38up)
if (nrow(kegg38up)){
  kegg38up$termRatioStr =paste0(kegg38up$Count,'/', unlist(strsplit(kegg38up$BgRatio, '/'))[seq(1,nrow(kegg38up)*2,2)])
  kegg38up$termRatio =round(kegg38up$Count/as.numeric(unlist(strsplit(kegg38up$BgRatio, '/'))[seq(1,nrow(kegg38up)*2,2)]),2)
}
kegg38up = kegg38up[order(kegg38up$pvalue, decreasing = F),]
write.csv(kegg38up, "D:\\gzy/38um_UP_KEGGanalysis.csv", row.names = F)
if (nrow(kegg38up)>10){
  kegg38up = kegg38up[1:10,]
}
kegg38up$Description = factor(kegg38up$Description,levels = as.character(kegg38up$Description))
rownames(kegg38up) = kegg38up[,2]
re = rev(as.character(kegg38up$Description))
kegg38up = kegg38up[re,]
p_kegg38up <- ggplot(kegg38up, aes(x=-log10( pvalue ), y=Description)) +
  geom_point(aes( size= termRatio, colour = Count)  ) + scale_y_discrete(limits=kegg38up$Description,labels = function(x) str_wrap(x, width = 40))+
  ggtitle("KEGG enrichment")  +  scale_color_gradient(low = 'blue', high = 'red') + #xlim(range(log10(out$pvalue))) +
  theme(axis.text.x=element_text(angle=0,size=10, vjust=0.7), axis.text.y=element_text(angle=0,size=12, vjust=0.7),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16),
        panel.background = element_rect(fill="white", colour='gray'),
        panel.grid.major = element_line(size = 0.05, colour = "gray"),
        panel.grid.minor.y = element_line(size=0.05, colour="gray"),
        panel.grid.minor.x = element_line(size=0.05, colour="gray"))
ggsave("D:\\gzy/kegg38up.tiff",p_kegg38up)
##38 down
entrezid_38down = bitr(down38, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrezid_38down = as.character(entrezid_38down[,2])
go38down <- enrichGO( gene   = entrezid_38down,
                      OrgDb  = org.Hs.eg.db,
                      universe = entrezid_38all,
                      ont   = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      readable = T,
                      qvalueCutoff  = 0.05)
go38down = data.frame(go38down)
if (nrow(go38down)){
  go38down$termRatioStr =paste0(go38down$Count,'/', unlist(strsplit(go38down$BgRatio, '/'))[seq(1,nrow(go38down)*2,2)])
  go38down$termRatio = round(go38down$Count/as.numeric(unlist(strsplit(go38down$BgRatio, '/'))[seq(1,nrow(go38down)*2,2)]),2)
}
go38down = go38down[order(go38down$pvalue, decreasing=F),]
write.csv(go38down,"D:\\gzy/go38down.csv",quote = F)
if (nrow(go38down)>10){
  go38down = go38down[1:10,]
}
go38down$Description = factor(go38down$Description,levels = as.character(go38down$Description))
rownames(go38down) = go38down[,3]
re = rev(as.character(go38down$Description))
go38down = go38down[re,]

p_go38down <- ggplot(go38down, aes(x=-log10( pvalue ), y=Description)) +
  geom_point(aes(size= termRatio , colour = Count)  ) + scale_y_discrete(limits=go38down$Description,labels = function(x) str_wrap(x, width = 40))+
  ggtitle("GO enrichment")  +  scale_color_gradient(low = 'blue', high = 'red') + #xlim(range(log10(out$pvalue))) +
  #scale_y_discrete(labels = function(x) str_wrap(x, width = 40) )
  theme(axis.text.x=element_text(angle=0,size=10, vjust=0.7), axis.text.y=element_text(angle=0,size=12, vjust=0.7),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16),
        panel.background = element_rect(fill="white", colour='gray'),
        panel.grid.major = element_line(size = 0.05, colour = "gray"),
        panel.grid.minor.y = element_line(size=0.05, colour="gray"),
        panel.grid.minor.x = element_line(size=0.05, colour="gray"))
ggsave("D:\\gzy/go38down.tiff",p_go38down,width = 6.4)

kegg38down = enrichKEGG(entrezid_38down, universe = entrezid_38all)
kegg38down = data.frame(kegg38down)
if (nrow(kegg38down)){
  kegg38down$termRatioStr =paste0(kegg38down$Count,'/', unlist(strsplit(kegg38down$BgRatio, '/'))[seq(1,nrow(kegg38down)*2,2)])
  kegg38down$termRatio =round(kegg38down$Count/as.numeric(unlist(strsplit(kegg38down$BgRatio, '/'))[seq(1,nrow(kegg38down)*2,2)]),2)
}
kegg38down = kegg38down[order(kegg38up$pvalue, decreasing = F),]
write.csv(kegg38down,"D:\\gzy/kegg38down.csv",quote = F)
#write.csv(kegg, "G:\\gzy/38um_UP_KEGGanalysis.csv", row.names = F)
if (nrow(kegg38down)>10){
  kegg38down = kegg38down[1:10,]
}
kegg38down$Description = factor(kegg38down$Description,levels = as.character(kegg38down$Description))
rownames(kegg38down) = kegg38down[,2]
re = rev(as.character(kegg38down$Description))
kegg38down = kegg38down[re,]
p_kegg38down <- ggplot(kegg38down, aes(x=-log10( pvalue ), y=Description)) +
  geom_point(aes( size= termRatio, colour = Count)  ) + scale_y_discrete(limits=kegg38down$Description,labels = function(x) str_wrap(x, width = 40))+
  ggtitle("KEGG enrichment")  +  scale_color_gradient(low = 'blue', high = 'red') + #xlim(range(log10(out$pvalue))) +
  theme(axis.text.x=element_text(angle=0,size=10, vjust=0.7), axis.text.y=element_text(angle=0,size=12, vjust=0.7),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16),
        panel.background = element_rect(fill="white", colour='gray'),
        panel.grid.major = element_line(size = 0.05, colour = "gray"),
        panel.grid.minor.y = element_line(size=0.05, colour="gray"),
        panel.grid.minor.x = element_line(size=0.05, colour="gray"))
ggsave("D:\\gzy/kegg38down.tiff",p_kegg38down)
####response to endoplasmic reticulum stress

deseq38 = read.csv("D:\\gzy/38um_remove0_annotation_delmin_DESeq2.csv",header = T)
deseq50 = read.csv("D:\\gzy/50um_remove0_annotation_delmin_DESeq2.csv",header = T)
deseq75 = read.csv("D:\\gzy/75um_remove0_annotation_delmin_DESeq2.csv",header = T)

up38 = deseq38$X[deseq38$padj <= 0.05 & deseq38$log2FoldChange >= 0.58]
up50 = deseq50$X[deseq50$padj <= 0.05 & deseq50$log2FoldChange >= 0.58]
up75 = deseq75$X[deseq75$padj <= 0.05 & deseq75$log2FoldChange >= 0.58]
#up38 = c(up38,"ENSG00000116717")

entrezid_38up = bitr(up38, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrezid_38up = as.character(entrezid_38up[,2])

entrezid_50up = bitr(up50, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrezid_50up = as.character(entrezid_50up[,2])

entrezid_75up = bitr(up75, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
entrezid_75up = as.character(entrezid_75up[,2])

go38up <- enrichGO(gene   = entrezid_38up,
                   OrgDb  = org.Hs.eg.db,
                   universe = NULL,
                   ont   = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable = T,
                   qvalueCutoff  = 0.05)
go38up = data.frame(go38up)
row.names(go38up) = go38up[,3]
go38up = go38up["response to endoplasmic reticulum stress",]
go38up = go38up$geneID

go50up <- enrichGO(gene   = entrezid_50up,
                   OrgDb  = org.Hs.eg.db,
                   universe = NULL,
                   ont   = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable = T,
                   qvalueCutoff  = 0.05)
go50up = data.frame(go50up)
row.names(go50up) = go50up[,3]
go50up = go50up["response to endoplasmic reticulum stress",]
go50up = go50up$geneID

go75up <- enrichGO(gene   = entrezid_75up,
                   OrgDb  = org.Hs.eg.db,
                   universe = NULL,
                   ont   = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable = T,
                   qvalueCutoff  = 0.05)
go75up = data.frame(go75up)
row.names(go75up) = go75up[,3]
go75up = go75up["response to endoplasmic reticulum stress",]
go75up = go75up$geneID

go38up = unlist(strsplit(go38up,"/"))
go50up = unlist(strsplit(go50up,"/"))
go75up = unlist(strsplit(go75up,"/"))
ergene =c()
ergene = union(go38up,go50up)
ergene = union(ergene,go75up)

#gadd = bitr("ENSG00000116717", fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#gadd = gadd$ENTREZID
ergene = union(ergene,"GADD45A")

deseq38 = deseq38[deseq38$Symbol %in% ergene,]
deseq38 = deseq38[,c("Symbol","log2FoldChange")]
rownames(deseq38) = NULL
colnames(deseq38) = c("Symbol","38")

deseq50 = deseq50[deseq50$Symbol %in% ergene,]
deseq50 = deseq50[,c("Symbol","log2FoldChange")]
rownames(deseq50) = NULL
colnames(deseq50) = c("Symbol","50")

deseq75 = deseq75[deseq75$Symbol %in% ergene,]
deseq75 = deseq75[,c("Symbol","log2FoldChange")]
rownames(deseq75) = NULL
colnames(deseq75) = c("Symbol","75")

deseq = merge(deseq38,deseq50,by = "Symbol")
deseq = merge(deseq,deseq75,by = "Symbol")
rownames(deseq) = deseq[,1]
deseq = deseq[,-1]
p = pheatmap(t(deseq),color = colorRampPalette(c( "white", "red"))(50),
             cluster_row = FALSE,
             cluster_rows = FALSE, 
             cellwidth =9, cellheight =9)
ggsave("D:\\gzy/heatmap.tiff",p,width = 10)

gn=rownames(t(deseq))[p$tree_row[["order"]]]
sn=colnames(t(deseq))[p$tree_col[["order"]]]
new_deseq=deseq[gn,sn]


wwcgene = "AKT2,BIRC5,CCND1,FASN,GADD45A,GADD45B,GADD45G,JMY,JUN,MDM4,PIK3R3,PMAIP1,RPRM,SFN,TNFRSF10B,TP53INP1,TP63,TP73"
wwcgene = unlist(strsplit(wwcgene,","))

deseq38 = read.csv("D:\\gzy/38um_remove0_annotation_delmin_DESeq2.csv",header = T)
#deseq50 = read.csv("D:\\gzy/50um_remove0_annotation_delmin_DESeq2.csv",header = T)
#deseq75 = read.csv("D:\\gzy/75um_remove0_annotation_delmin_DESeq2.csv",header = T)
wwcgene = deseq38$X[deseq38$Symbol %in% wwcgene]

all38 = read.table("D:\\gzy/control_38.txt",header = T)
all50 = read.table("D:\\gzy/control_50.txt",header = T)
all75 = read.table("D:\\gzy/control_75.txt",header = T)

all38 = all38[wwcgene,]
all50 = all50[wwcgene,]
all75 = all75[wwcgene,]

all50 = all50[,c(4,5,6)]
all75 = all75[,c(4,5,6)]

allgene = merge(all38,all50,by = "row.names")
rownames(allgene) = allgene[,1]
allgene = allgene[,-1]

allgene = merge(allgene, all75,by = "row.names")
write.csv(allgene,"D:\\gzy/gene18.csv")
##########
mito = read.csv("D:\\gzy/go75down.csv",header = T)
mito_gene = mito$geneID[grepl("mito",mito$Description)]

mito = c()
for (i in mito_gene){
  a = unlist(strsplit(i,"/"))
  mito = c(mito,a)
}
mito = unique(mito)

deseq38 = read.csv("D:\\gzy/38um_remove0_annotation_delmin_DESeq2.csv",header = T)
deseq50 = read.csv("D:\\gzy/50um_remove0_annotation_delmin_DESeq2.csv",header = T)
deseq75 = read.csv("D:\\gzy/75um_remove0_annotation_delmin_DESeq2.csv",header = T)

deseq38 = deseq38[deseq38$Symbol %in% mito,]
deseq38 = deseq38[,c("Symbol","log2FoldChange")]
rownames(deseq38) = NULL
colnames(deseq38) = c("Symbol","38")

deseq50 = deseq50[deseq50$Symbol %in% mito,]
deseq50 = deseq50[,c("Symbol","log2FoldChange")]
rownames(deseq50) = NULL
colnames(deseq50) = c("Symbol","50")

deseq75 = deseq75[deseq75$Symbol %in% mito,]
deseq75 = deseq75[,c("Symbol","log2FoldChange")]
rownames(deseq75) = NULL
colnames(deseq75) = c("Symbol","75")

deseq = merge(deseq38,deseq50,by = "Symbol")
deseq = merge(deseq,deseq75,by = "Symbol")
rownames(deseq) = deseq[,1]
deseq = deseq[,-1]
p = pheatmap(t(deseq),color = colorRampPalette(c(  "blue","white"))(50),
             cluster_row = FALSE,
             cluster_rows = FALSE, 
             cellwidth =9, cellheight =9)

deseq38 = read.csv("D:\\gzy/38um_remove0_annotation_delmin_DESeq2.csv",header = T)
deseq50 = read.csv("D:\\gzy/50um_remove0_annotation_delmin_DESeq2.csv",header = T)
deseq75 = read.csv("D:\\gzy/75um_remove0_annotation_delmin_DESeq2.csv",header = T)

deseq38up = deseq38[deseq38$padj <= 0.05 & deseq38$log2FoldChange >= 0.58,]
deseq38down = deseq38[deseq38$padj <= 0.05 & deseq38$log2FoldChange <= -0.58,]

deseq50up = deseq50[deseq50$padj <= 0.05 & deseq50$log2FoldChange >= 0.58,]
deseq50down = deseq50[deseq50$padj <= 0.05 & deseq50$log2FoldChange <= -0.58,]

deseq75up = deseq75[deseq75$padj <= 0.05 & deseq75$log2FoldChange >= 0.58,]
deseq75down = deseq75[deseq75$padj <= 0.05 & deseq75$log2FoldChange <= -0.58,]

write.csv(deseq38up,"D:\\gzy/deseq38up.csv")
write.csv(deseq38down,"D:\\gzy/deseq38down.csv")
write.csv(deseq50up,"D:\\gzy/deseq50up.csv")
write.csv(deseq50down,"D:\\gzy/deseq50down.csv")
write.csv(deseq75up,"D:\\gzy/deseq75up.csv")
write.csv(deseq75down,"D:\\gzy/deseq75dwn.csv")
