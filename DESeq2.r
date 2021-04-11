rm(list = ls())
library(DESeq2)
mydata <- read.table("E:\\gzy/control_50um_remove0.txt", header = TRUE, quote = '\t')

condition <- factor(c("control","control","control","treat","treat","treat"))
colData <- data.frame(row.names=colnames(mydata), condition)

dds <- DESeqDataSetFromMatrix(round(mydata), DataFrame(condition), design= ~ condition )
dds2 <- DESeq(dds)

res <- results(dds2)
write.csv(res,"E:\\gzy/50um_remove0_DESeq2.csv")
#table(res$padj<0.05)

library("biomaRt")
library("curl")
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
my_ensembl_gene_id<-row.names(res)
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = ensembl)
ensembl_gene_id<-rownames(res)
res<-cbind(ensembl_gene_id,res)
colnames(res)[1]<-c("ensembl_gene_id")
diff_name<-merge(res,mms_symbols,by="ensembl_gene_id")
rownames(diff_name)<-diff_name[,1]
diff_name<-diff_name[,-1]
write.csv(diff_name,"E:\\gzy/50um_remove0_annotation_DESeq2.csv")

library(ggpubr)
jpeg("E:\\gzy/50um_remove0.jpeg",width=800,height=600,res=120)
ggmaplot(diff_name, 
               fdr = 0.05, fc = 1.5, size = 0.7,
               palette = c("#B31B21", "#1465AC", "darkgray"),
               genenames = as.vector(diff_name$external_gene_name),
               xlab = "BaseMean",ylab = "log2(FC)",
               legend = "top", top = 20,
               font.label = c("bold", 11),
               font.legend = "bold",
               font.main = "bold",
               ggtheme = ggplot2::theme_minimal())
dev.off()
