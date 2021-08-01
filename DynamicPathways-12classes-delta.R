library(ggplot2)
library(reshape2)
library(mdscore)
library(org.Hs.eg.db)
library(stats)
library(clusterProfiler)
library(MAGeCKFlute)
library(ReactomePA)
library(ggpubr)
library(ggsci)
library(matrixStats)
source('D:\\wwc/noEnrichPlot.R')

source('D:\\wwc/GOKEGGPlot.R')
setwd('D:\\wwc2/')
sabDB = loadDb('org.H')

H2 = read.csv('D:\\gzy/38um_remove0_annotation_delmin_DESeq2.csv')
H4 = read.csv('D:\\gzy/50um_remove0_annotation_delmin_DESeq2.csv')
H6 = read.csv('D:\\gzy/75um_remove0_annotation_delmin_DESeq2.csv')

outdir = 'D:\\wwc2/'
#dir.create(outdir)
setwd(outdir)

data1 = data.frame('Gene' = H2$X,'c38.LFC' = H2$log2FoldChange, 'c38.padj' = H2$padj)
data2 = data.frame('Gene' = H4$X,'c50.LFC' = H4$log2FoldChange, 'c50.padj' = H4$padj)
data3 = data.frame('Gene' = H6$X,'c75.LFC' = H6$log2FoldChange, 'c75.padj' = H6$padj)
##reserve genes expressed(TPM>1) in at least one time point
data = merge(data1,data2, by.x = 'Gene', by.y = 'Gene', all.x = T, all.y = T)
data = merge(data, data3, by.x = 'Gene',by.y = 'Gene', all.x=T, all.y=T)
#remove all 'na' values by give them 0 (for LFC) or 1 (for padj)
data$c38.LFC[is.na(data["c38.LFC"])] = 0
data$c50.LFC[is.na(data["c50.LFC"]) ] = 0
data$c75.LFC[is.na(data['c75.LFC'])] = 0
data$c38.padj[is.na(data$c38.padj)] = 1
data$c50.padj[is.na(data$c50.padj)] = 1
data$c75.padj[is.na(data$c75.padj)] = 1

#samples = str_split(as.character(colnames(data)[grep('LFC', colnames(data))]), '.LFC', simplify = T)[,1]
samples = as.character(colnames(data)[grep('LFC', colnames(data))])
padjs = as.character(colnames(data)[grep('padj', colnames(data))])

#expdata = data[!is.na(rowMeans(data[,as.character(samples)])),]
##data = data[(rowMins(as.matrix(data[,as.character(padjs)]))<0.05) &(!is.na(rowMaxs(as.matrix(data[,as.character(padjs)]))<0.05)),]
thred = 0.58
data = data[(rowMaxs(as.matrix(abs(data[,as.character(samples)]))) > thred)&(rowMins(as.matrix(data[,as.character(padjs)]))<0.05) &(!is.na(rowMaxs(as.matrix(data[,as.character(padjs)]))<0.05)),]


data$c38.LFC.slope = data$c38.LFC
for (i in 1:(length(samples)-1)){
  pairsample = samples[i:(i+1)]
  timeseq = c(i,i+1)
  #title = paste0(sample, 'vs', pre)
  title = samples[i+1]
  relation = apply(data[as.character(pairsample)], 1, function(x) glm(x ~ timeseq))
  waldstat = apply(data[as.character(pairsample)], 1, function(x) wald.test(glm(x ~ timeseq), terms =2))
  data[paste0(title, '.slope')] = sapply(relation, function(x) x[1]$coefficients[2])
  data[paste0(title, '.intercept')] = sapply(relation, function(x) x[1]$coefficients[1])
  #data[paste0(title, '.wald')] = sapply(waldstat, function(x) x$W)
  #data[paste0(title, '.waldpvalue')] = sapply(waldstat, function(x) x$pvalue)
}
# timeseq = c(1,2,3)
# #title = paste0(sample, 'vs', pre)
# title = 'All'
# relation = apply(data[as.character(samples)], 1, function(x) glm(x ~ timeseq))
# data[paste0(title, '.slope')] = sapply(relation, function(x) x[1]$coefficients[2])
# data[paste0(title, '.intercept')] = sapply(relation, function(x) x[1]$coefficients[1])

ggplot(data) + geom_histogram(aes(c38.LFC.slope), binwidth = 0.5)
ggsave('c38.slope.hist.pdf')
ggplot(data) + geom_histogram(aes(c50.LFC.slope), binwidth = 0.5)
ggplot(data) + geom_histogram(aes(c75.LFC.slope), binwidth = 0.5)
#ggplot(data) + geom_histogram(aes(All.slope), binwidth = 0.5)
nrow(data[data$c50.LFC.slope > 0 &(data$c50 < 1),])

subdata0 = data[, c('Gene', 'c38.LFC.slope')]
colnames(subdata0) = c('Gene', 'Slope')
subdata0$Group = '0-38m'
subdata1 = data[, c('Gene', 'c50.LFC.slope')]
colnames(subdata1) = c('Gene', 'Slope')
subdata1$Group = '38-50m'
subdata2 = data[, c('Gene', 'c75.LFC.slope')]
colnames(subdata2) = c('Gene', 'Slope')
subdata2$Group = '50-75m'
#subdata3 = data[, c('Gene', 'All.slope')]
#colnames(subdata3) = c('Gene', 'Slope')
#subdata3$Group = '3h-24h'
slopes = rbind(subdata0,subdata1, subdata2)

slopes$Group = factor(slopes$Group,levels= unique(slopes$Group))
ggplot(slopes) + geom_violin(aes(Group, Slope)) + geom_abline(slope = 0,intercept = 1, linetype = 'dotted') + 
  geom_abline(slope = 0,intercept = -1, linetype = 'dotted') +  
  theme(panel.background = element_rect(fill="white", colour='gray'))
ggsave('violin_timepoint_slope.pdf', units = 'in', dpi = 600, width = 3, height = 3)
slopes$Group = factor(slopes$Group,levels= unique(slopes$Group))
ggplot(slopes,aes(Group, Slope)) + geom_violin() + 
  geom_jitter(aes(color=Group)) + #geom_jitter()+
  geom_abline(slope = 0,intercept = 0.58, linetype = 'dotted') + 
  geom_abline(slope = 0,intercept = -0.58, linetype = 'dotted') +  
  theme(panel.background = element_rect(fill="white", colour='gray'))
ggsave('violin_thmepoint_slope_cutoff0.58.pdf', units = 'in', dpi = 600, width = 3, height = 3)
t.test(data$c50.LFC.slope,conf.level = 0.90)
t.test(data$c75.LFC.slope,conf.level = 0.99)
#t.test(data$All.slope)
t.test(data$c38.LFC.slope, conf.level = 0.99)

data$group = 'No'
data$group[data$All.slope>1] = 'Up'
data$group[data$All.slope< -1] = 'Down'
subdata = data[,c('Gene', 'c38.LFC', 'c50.LFC', 'c75.LFC', 'group')]
subdata4 = melt(subdata, id.vars = c('Gene','group'))


ggplot(subdata4, aes(variable, value,color=group, group = Gene)) + geom_point() +
  geom_line(alpha= 0.5) + labs(y = 'Log2FoldChange') 
ggsave('pointplot_dynamic.pdf', dpi=600, width = 5, height = 6)

cutoff = 0
data$classone = 0
data$classone[data$c38.LFC.slope > thred] = 2
data$classone[data$c38.LFC.slope < -thred] = 1
data$classtwo = 0
data$classtwo[data$c50.LFC.slope > cutoff] = 6
data$classtwo[data$c50.LFC.slope < -cutoff] = 3
data$classthree = 0
data$classthree[data$c75.LFC.slope > cutoff] = 18
data$classthree[data$c75.LFC.slope < -cutoff] = 9

data$class = data$classone + data$classtwo + data$classthree

subdata = data[,c('Gene', 'c38.LFC', 'c50.LFC', 'c75.LFC', 'class')]
subdata5 = melt(subdata, id.vars = c('Gene','class'))
ggplot(subdata5, aes(variable, value,color=class, group = Gene)) + geom_point() +
  geom_line(alpha= 0.5) + labs(y = 'Log2FoldChange') 
ggsave('pointplot_dynamic_27classes.pdf', dpi=600, width = 5, height = 6)

# outdir = 'G:/01.projects/COVID-19/DESeqResults/DESeqPipeline/EnrichAnalysis-27classes-cutoff0.58'
# dir.create(outdir)
# setwd(outdir)

symbols = bitr(as.character(data$Gene), fromType = 'ENSEMBL', toType = c('ENSEMBL','SYMBOL', 'ENTREZID'), OrgDb = "org.Hs.eg.db", drop = FALSE)
#symbols = symbols[!duplicated(symbols$'ENSEMBL'),]
print(nrow(symbols))

symbols$'Gene' = symbols[, which(colnames(symbols) == 'ENSEMBL')]

data = merge(data, symbols, by.x = 'Gene', by.y = 'Gene', all.x = T)
#res$Gene = res$hgnc_symbol
colnames(data)[colnames(data)=='ENTREZID'] = 'EntrezID'
colnames(data)[colnames(data)=='SYMBOL'] = 'Symbol'

write.table(data, 'delta_slopeSummary_pfiltered.csv', row.names = F, quote = F,sep = ',')
#expdata = merge(expdata, symbols, by.x='Gene', by.y='Gene', all.x = T)

for (i in unique(data$class)){
  subdata = data[data$class == i,]
  n = nrow(subdata)
  print(n)
  ggplot(subdata5[subdata5$class == i,], aes(variable, value, color= class, group = Gene)) + geom_point() +
    geom_line(alpha= 0.5) + labs(y = 'Log2FoldChange', title = paste(n, 'Genes')) 
  ggsave(paste0('pointplot_dynamic_delta_class-', i,'.pdf'), dpi=600, width = 4, height = 5)
  
  bg = na.omit(data$ENTREZID) ###change background genes list from expdata to data, more strigent
  outdir = paste0('class-',i)
  dir.create(outdir)
  genelist = na.omit(subdata$EntrezID)
  GOKEGGPlot(genelist,outdir=outdir, main = outdir, kegg_organism = 'hsa',universe=bg ,GO_orgdb = "org.Hs.eg.db")
}

#################### plot several classes in one fig ##############
mycolors = c("#8dd3c7", 
             "#fb8072", "#80b1d3", "#fdb462", 
             "#bc80bd", "#b3de69", "#bebada")
for(i in 1:round(max(data$class)/4)){
  subdata = data[data$class %in% sort(unique(data$class))[-1][((i-1)*4+1):(i*4)],]
  nclass = length(unique(subdata$class))
  n = nrow(subdata)
  names(mycolors) = unique(subdata$class)
  ggplot(subdata5[subdata5$class %in% subdata$class,], aes(variable, value, color= as.character(class), group = Gene)) + geom_point() +
    scale_color_manual('Categories',values = mycolors[1:nclass])  +
    geom_line(alpha = 0.6) + labs(x= 'Varibles', y = 'Log2FoldChange', title = paste(n, 'Genes')) + scale_fill_npg() + theme_bw() +
    theme(axis.text.x=element_text(angle = 20,vjust =0.7,size=10) )
  ggsave(paste0('pointplot_dynamic_delta_subclasses-', i,'.pdf'), dpi=600, width = 3, height = 4)
  
}



######### check one specific dynamic group (class-15)
# goterms = read.csv('class-15/class-15group_GO_enrich.csv',header = T)
# goterms = goterms[order(goterms$p.adjust),]
# genes = goterms$geneID[1]
# genes = unlist(strsplit(as.character(genes), '/'))
# ensembles = data$Gene[data$Symbol %in% genes]
# class15 = subdata5[subdata5$class == 15,]
# n=nrow(class15)
# class15$highlight = 'None'
# class15$highlight[class15$Gene %in% ensembles] = 'inflammatory'
# ggplot(class15, aes(variable, value, color= highlight, group = Gene)) + geom_point() +
#   geom_line(alpha= 0.5) + labs(y = 'Log2FoldChange', title = paste(n, 'Genes')) 
# ggsave('class15_highlighted_pointplot.pdf')
# 
# 
# keggterms = read.csv('class-15/class-15group_KEGG_enrich.csv',header = T)
# keggterms = keggterms[order(keggterms$p.adjust),]
# genes = keggterms$geneID[1]
# genes = unlist(strsplit(as.character(genes), '/'))
# ensembles = data$Gene[data$EntrezID %in% genes]
# class15 = subdata5[subdata5$class == 15,]
# n=nrow(class15)
# class15$highlight = 'None'
# class15$highlight[class15$Gene %in% ensembles] = 'KEGG'
# ggplot(class15, aes(variable, value, color= highlight, group = Gene)) + geom_point() +
#   geom_line(alpha= 0.5) + labs(y = 'Log2FoldChange', title = paste(n, 'Genes')) 
# ggsave('class15_highlighted_pointplot_kegg.pdf')
# 
# 
