# Loading package: TCGAbiolinks,maftools
library(TCGAbiolinks)
library(maftools)
# Mutation data
PAAD.mutect.maf <- GDCquery_Maf('PAAD',pipelines = 'mutect2')
laml <- read.maf(maf = PAAD.mutect.maf, isTCGA = TRUE)
snpdata <- laml@data
snpdata <- snpdata[,c("Hugo_Symbol","Variant_Type","Variant_Classification","Tumor_Sample_Barcode")]
TP53.snp <- subset(snpdata,Hugo_Symbol=='TP53')
silent <- laml@maf.silent
silent <- silent[,c("Hugo_Symbol","Variant_Type","Variant_Classification","Tumor_Sample_Barcode")]
TP53.silent <- subset(silent,Hugo_Symbol=='TP53')
tidydata <- PAAD.mutect.maf[,c("Hugo_Symbol","Variant_Type",
                               "Variant_Classification","Tumor_Sample_Barcode")]
library(reshape2)
rsdata <- dcast(tidydata, Hugo_Symbol~Tumor_Sample_Barcode)
rownames(rsdata) <- rsdata$Hugo_Symbol
rsdata <- rsdata[,-1]
Gene <- rsdata['TP53',]
names(Gene) <-  sapply(strsplit(names(Gene),'-'),function(x) paste0(x[1:3],collapse="-"))
MAF_barcode <-  sapply(strsplit(names(Gene),'-'),function(x) paste0(x[1:3],collapse="-"))

# Clinical data
library(dplyr)
clinical <- GDCquery_clinic(project = "TCGA-PAAD", type = "clinical")
clinical <- clinical[,c("submitter_id","vital_status","days_to_death","days_to_last_follow_up")]
which(is.na(clinical$vital_status),arr.ind = T)
clinical[clinical$vital_status=='Dead',]$vital_status <- 1
clinical[clinical$vital_status=='Alive',]$vital_status <- 0
clinical$vital_status <- as.numeric(clinical$vital_status)
clinical$time <- clinical$days_to_death
clinical$time[which(is.na(clinical$time))] <- clinical$days_to_last_follow_up[which(is.na(clinical$time))]
clinical <- filter(clinical,time != "--")
clinical$time <- as.numeric(clinical$time)

# Data fusion: mutation data and clinical data
interbarcode <- intersect(MAF_barcode,clinical$submitter_id)
Element <- c('TP53')
geneMut <- Gene[Element,interbarcode]
geneMut <- t(geneMut) %>% as.data.frame() %>% data.matrix() %>% as.data.frame()
geneMut$submitter_id <- rownames(geneMut)
gene <- c('ATF3','TP53')
mergdata <- merge(clinical,geneMut,by = "submitter_id")
mergdata$Group <- ifelse(mergdata[,'TP53']==0,"Wild Type","Mutant")
# Survival analysis
library(survival)
library(survminer)
fit <- survfit(Surv(time, vital_status) ~ Group, data = mergdata)
table(mergdata$Group)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "darkblue", "bold"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           legend.labs = c("Mutant(n=112)", "Wild Type(n=66)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw()
)
