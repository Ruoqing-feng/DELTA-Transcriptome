# Loading package: TCGAbiolinks,maftools
library(TCGAbiolinks)
library(maftools)
# Gene expression data of BRCA
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")
query <- GDCquery(project = "TCGA-BRCA", 
                  experimental.strategy = "RNA-Seq",
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - FPKM")
GDCdownload(query)
BRCARnaseq <- GDCprepare(query)
BRCA_FPKM <- assay(BRCARnaseq)
write.csv(BRCA_FPKM, file = 'TCGA-BRCA-Counts.csv')
BRCA_FPKM <- read.csv('TCGA-BRCA-Counts.csv',header = TRUE,sep = ',',check.names = FALSE)
rownames(BRCA_FPKM) <- BRCA_FPKM[,1]
View(BRCA_FPKM)
BRCA_FPKM <- BRCA_FPKM[,-1]
library("SummarizedExperiment")
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(BRCA_FPKM),typesample = c("TP"))
BRCA_FPKM <- BRCA_FPKM[c('ENSG00000162772','ENSG00000141510'),samplesTP]
row.names(BRCA_FPKM) <- c("ATF3","TP53")
names(BRCA_FPKM) <-  sapply(strsplit(names(BRCA_FPKM),'-'),function(x) paste0(x[1:3],collapse="-"))
barcode <- sapply(strsplit(names(BRCA_FPKM),'-'),function(x) paste0(x[1:3],collapse="-"))
# Clinical data
clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
library(dplyr)
clinical <- clinical[,c("submitter_id","vital_status","days_to_death","days_to_last_follow_up")]
which(is.na(clinical$vital_status),arr.ind = T)
clinical <- clinical[-1098,]
clinical[clinical$vital_status=='Dead',]$vital_status <- 1
clinical[clinical$vital_status=='Alive',]$vital_status <- 0
clinical$vital_status <- as.numeric(clinical$vital_status)
clinical$time <- clinical$days_to_death
clinical$time[which(is.na(clinical$time))] <- clinical$days_to_last_follow_up[which(is.na(clinical$time))]
clinical <- filter(clinical,time != "--")
clinical$time <- as.numeric(clinical$time)

# Mutation data
BRCA.mutect.maf <- GDCquery_Maf('BRCA',pipelines = 'mutect2')
laml <- read.maf(maf = BRCA.mutect.maf, isTCGA = TRUE)
snpdata <- laml@data
snpdata <- snpdata[,c("Hugo_Symbol","Variant_Type",
                      "Variant_Classification","Tumor_Sample_Barcode")]
TP53.snp <- subset(snpdata,Hugo_Symbol=='TP53')
silent <- laml@maf.silent
silent <- silent[,c("Hugo_Symbol","Variant_Type",
                    "Variant_Classification","Tumor_Sample_Barcode")]
TP53.silent <- subset(silent,Hugo_Symbol=='TP53')
interbarcode <- intersect(TP53.snp$Tumor_Sample_Barcode,TP53.silent$Tumor_Sample_Barcode)
tidydata <- BRCA.mutect.maf[,c("Hugo_Symbol","Variant_Type",
                               "Variant_Classification","Tumor_Sample_Barcode")]
library(reshape2)
rsdata <- dcast(tidydata, Hugo_Symbol~Tumor_Sample_Barcode)
rownames(rsdata) <- rsdata$Hugo_Symbol
rsdata <- rsdata[,-1]
Gene <- rsdata['TP53',]
names(Gene) <-  sapply(strsplit(names(Gene),'-'),function(x) paste0(x[1:3],collapse="-"))
MAF_barcode <-  sapply(strsplit(names(Gene),'-'),function(x) paste0(x[1:3],collapse="-"))
Gene[,'TCGA-BH-A0BW'] <- 0
Gene[,'TCGA-E2-A1B4'] <- 0
Gene[,'TCGA-A8-A09Z'] <- 0
Gene[,'TCGA-E9-A1RF'] <- 0

# Survival Analysis: Mutation Status
# data fusion: mutation data and clinical data
library(survival)
library(survminer)
Gene <- t(Gene) %>% as.data.frame()
Gene$submitter_id <- rownames(Gene)
interbarcode <- intersect(Gene$submitter_id,clinical$submitter_id)
Element <- c('TP53')
Gene <- t(Gene) %>% as.data.frame()
Gene <- Gene[Element,interbarcode]
Gene <- t(Gene) %>% as.data.frame()
Gene$submitter_id <- rownames(Gene)
mergdata <- merge(clinical,Gene,by = "submitter_id")
mergdata$Group <- ifelse(mergdata[,'TP53']==0,"Wild Type","Mutant")
table(mergdata$Group)
fit <- survfit(Surv(time, vital_status) ~ Group, data=mergdata)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "darkblue", "bold"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           ylim = c(0,1),
           legend.labs = c("Mutant(n=338)", "Wild Type(n=648)"),
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)

# 将ATF3、TP53表达矩阵与TP53突变数据、临床数据融合
Gene_wt <- subset(Gene,TP53=='0')
Gene_wt$Variant_Classification <- 'WT'
Gene_wt$Tumor_Sample_Barcode <- rownames(Gene_wt)
Gene_mut <- subset(Gene,TP53!='0')
Mut <- subset(snpdata,Hugo_Symbol=='TP53')
Mut <- Mut[,c("Variant_Classification","Tumor_Sample_Barcode")]
Gene_mut$Tumor_Sample_Barcode <- rownames(Gene_mut)
mergdata <- merge(Mut,Gene_mut,by = "Tumor_Sample_Barcode")
Mut <- rbind(Gene_wt,mergdata)
Mut <- Mut[,-1]
interbarcode <- intersect(barcode,Mut$Tumor_Sample_Barcode)
interbarcode <- intersect(interbarcode,clinical$submitter_id)
gene <- c('ATF3','TP53')
geneExp <- BRCA_FPKM[gene,interbarcode]
geneExp <- t(geneExp) %>% as.data.frame() %>% data.matrix() %>% as.data.frame()
geneExp$submitter_id <- rownames(geneExp)
Mut <- t(Mut) %>% as.data.frame()
View(Mut)
colnames(Mut) <- Mut[1,]
Element <- c('Variant_Classification')
Mut <- Mut[Element,interbarcode]
Mut <- t(Mut) %>% as.data.frame()
Mut$submitter_id <- rownames(Mut)
mergdata <- merge(geneExp,Mut,by = "submitter_id")
mergdata <- merge(clinical,mergdata,by = "submitter_id")
table(mergdata$Variant_Classification)
p <- ggplot(data = mergdata,aes(x = Variant_Classification,y = TP53,fill = Variant_Classification))+  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  coord_cartesian(ylim = c(0,80))
p
mutant <- mergdata[-which(mergdata[,'Variant_Classification'] %in% c('WT','Splice_Site')),]
mutant$Group <- ifelse(mutant$Variant_Classification %in% c("Missense_Mutation","In_Frame_Del"),"Non_truncation","Truncation")
table(mutant$Group)
truncation <- mutant[which(mutant[,'Group'] %in% c('Truncation')),]
non_truncation <- mutant[which(mutant[,'Group'] %in% c('Non_truncation')),]
### Expression: ATF3 in BRCA-TP53-truncation                       
sur.cut <- surv_cutpoint(truncation,time= 'time',
                         event = 'vital_status' ,
                         variables = 'ATF3' )
truncation1 <- surv_categorize(sur.cut)
table(truncation1$ATF3)
fit <- survfit(Surv(time, vital_status) ~ ATF3, data = truncation1)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "bold", "darkblue"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'ATF3',
           ylim = c(0,1),
           legend.labs = c("High(n=29)", "Low(n=76)"),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### Expression: TP53 in BRCA-TP53-truncation
sur.cut <- surv_cutpoint(truncation,time= 'time',
                         event = 'vital_status' ,
                         variables = 'TP53' )
truncation2 <- surv_categorize(sur.cut)
table(truncation2$TP53)
fit <- survfit(Surv(time, vital_status) ~ TP53, data = truncation2)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "bold", "darkblue"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           legend.labs = c("High(n=91)", "Low(n=14)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### Expression: ATF3 in BRCA-TP53-non_truncation
sur.cut <- surv_cutpoint(non_truncation,time= 'time',
                         event = 'vital_status' ,
                         variables = 'ATF3' )
non_truncation1 <- surv_categorize(sur.cut)
table(non_truncation1$ATF3)
fit <- survfit(Surv(time, vital_status) ~ ATF3, data = non_truncation1)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "bold", "darkblue"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'ATF3',
           ylim = c(0,1),
           legend.labs = c("High(n=24)", "Low(n=181)"),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### Expression: TP53 in BRCA-TP53-non_truncation
sur.cut <- surv_cutpoint(non_truncation,time= 'time',
                         event = 'vital_status' ,
                         variables = 'TP53' )
non_truncation2 <- surv_categorize(sur.cut)
table(non_truncation2$TP53)
fit <- survfit(Surv(time, vital_status) ~ TP53, data = non_truncation2)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "bold", "darkblue"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           legend.labs = c("High(n=62)", "Low(n=143)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)                       
data1 <- mergdata[which(mergdata[,'Variant_Classification'] %in% c('WT','Missense_Mutation','In_Frame_Del')),]
data2 <- mergdata[-which(mergdata[,'Variant_Classification'] %in% c('Missense_Mutation','In_Frame_Del','Splice_Site')),]
data1$Group <- ifelse(data1$Variant_Classification %in% c("WT"),"Wild Type","Non-truncation")
table(data1$Group)
data2$Group <- ifelse(data2$Variant_Classification %in% c("WT"),"Wild Type","Truncation")
table(data2$Group)
### Expression: TP53-WT vs TP53-non_truncation                       
fit <- survfit(Surv(time, vital_status) ~ Group, data=data1)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "darkblue", "bold"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           legend.labs = c("Non-truncation(n=205)", "Wild Type(n=644)"),
           ylim = c(0,1),
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### Expression: TP53-WT vs TP53-truncation
fit <- survfit(Surv(time, vital_status) ~ Group, data=data2)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "darkblue", "bold"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           legend.labs = c("Truncation(n=105)", "Wild Type(n=644)"),
           ylim = c(0,1),
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### Expression: TP53-truncation vs TP53-non_truncation
fit <- survfit(Surv(time, vital_status) ~ Group, data=mutant)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "darkblue", "bold"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'Status',
           legend.labs = c("Non-truncation(n=205)", "Truncation(n=105)"),
           ylim = c(0,1),
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)

## cutpoint
### Expression: ATF3 in TP53 wt
View(mergdata)
wt <- mergdata[which(mergdata[,'Variant_Classification'] %in% c('WT')),]
View(wt)
mut <- mergdata[-which(mergdata[,'Variant_Classification'] %in% c('WT')),]
View(mut)
sur.cut <- surv_cutpoint(wt,time= 'time',
                         event = 'vital_status' ,
                         variables = 'ATF3' )
wt1 <- surv_categorize(sur.cut)
table(wt1$ATF3)
fit <- survfit(Surv(time, vital_status) ~ ATF3, data = wt1)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "darkblue", "bold"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'ATF3',
           legend.labs = c("High(n=325)", "Low(n=319)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### Expession; TP53 in TP53 wt
sur.cut <- surv_cutpoint(wt,time= 'time',
                         event = 'vital_status' ,
                         variables = 'TP53' )
wt2 <- surv_categorize(sur.cut)
table(wt2$TP53)
fit <- survfit(Surv(time, vital_status) ~ TP53, data = wt2)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "bold", "darkblue"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           legend.labs = c("High(n=521)", "Low(n=123)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### ATF3 in TP53 mutant
sur.cut <- surv_cutpoint(mut,time= 'time',
                         event = 'vital_status' ,
                         variables = 'ATF3' )
mutant1 <- surv_categorize(sur.cut)
table(mutant1$ATF3)
fit <- survfit(Surv(time, vital_status) ~ ATF3, data = mutant1)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "bold", "darkblue"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'ATF3',
           legend.labs = c("High(n=34)", "Low(n=302)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### TP53 in TP53 mutant
sur.cut <- surv_cutpoint(mut,time= 'time',
                         event = 'vital_status' ,
                         variables = 'TP53' )
mutant2 <- surv_categorize(sur.cut)
table(mutant2$TP53)
fit <- survfit(Surv(time, vital_status) ~ TP53, data = mutant2)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "bold", "darkblue"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           legend.labs = c("High(n=187)", "Low(n=149)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)

# 将ATF3、TP53表达矩阵与临床数据融合
interbarcode <- intersect(barcode,clinical$submitter_id)
geneExp <- BRCA_FPKM[gene,interbarcode]
geneExp <- t(geneExp) %>% as.data.frame() %>% data.matrix() %>% as.data.frame()
geneExp$submitter_id <- rownames(geneExp)
mergdata <- merge(geneExp,clinical,by = "submitter_id")
# 生存分析
## cutpoint
### Expression: ATF3 in BRCA
sur.cut <- surv_cutpoint(mergdata,time= 'time',
                         event = 'vital_status' ,
                         variables = 'ATF3' )
mergdata1 <- surv_categorize(sur.cut)
table(mergdata1$ATF3)
fit <- survfit(Surv(time, vital_status) ~ ATF3, data = mergdata1)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "bold", "darkblue"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'ATF3',
           ylim = c(0,1),
           legend.labs = c("High(n=488)", "Low(n=602)"),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### Expression: TP53 in BRCA
sur.cut <- surv_cutpoint(mergdata,time= 'time',
                         event = 'vital_status' ,
                         variables = 'TP53' )
mergdata2 <- surv_categorize(sur.cut)
table(mergdata2$TP53)
fit <- survfit(Surv(time, vital_status) ~ TP53, data = mergdata2)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "bold", "darkblue"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           legend.labs = c("High(n=587)", "Low(n=503)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)

##### Drug

# 下载BRCA临床药物信息
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Clinical",
                  file.type = "xml")
query_download <- GDCdownload(query)
clinical.drug <- GDCprepare_clinic(query,"drug")

# 将TP53突变数据与临床用药数据和ATF3、TP53表达谱数据融合
interbarcode <- intersect(barcode,MAF_barcode)
interbarcode <- intersect(interbarcode,clinical.drug$bcr_patient_barcode)
geneExp <- BRCA_FPKM[gene,interbarcode]
geneExp <- t(geneExp) %>% as.data.frame() %>% data.matrix() %>% as.data.frame()
geneExp$bcr_patient_barcode <- rownames(geneExp)
Element <- c('TP53')
Gene <- t(Gene) %>% as.data.frame()
geneMut <- Gene[Element,interbarcode]
geneMut <- t(geneMut) %>% as.data.frame()
geneMut$bcr_patient_barcode <- rownames(geneMut)
mergdata <- merge(geneExp,geneMut,by = "bcr_patient_barcode")
mergdata <- merge(clinical.drug,mergdata,by = "bcr_patient_barcode")
names(mergdata)[26] <- 'TP53'
names(mergdata)[27] <- 'status'
mergdata$Group <- ifelse(mergdata[,'status']==0,"Wild Type","Mutant")
response <- mergdata[,c('bcr_patient_barcode','measure_of_response','ATF3','TP53','Group')]
response <- response[which(response[,'measure_of_response'] %in% c('Complete Response','Partial Response','Clinical Progressive Disease','Stable Disease')),]
response <- response[!duplicated(response), ]

# responder:Complete Response&Partial Response; 
# non-responder:Clinical Progressive Disease&Stable Disease
responder <- response[which(response[,'measure_of_response'] %in% c('Complete Response','Partial Response')),]
responder <- responder[!duplicated(responder$bcr_patient_barcode), ]
non_responder <- response[which(response[,'measure_of_response'] %in% c('Clinical Progressive Disease','Stable Disease')),]
non_responder <- non_responder[!duplicated(non_responder$bcr_patient_barcode), ]
interbarcode <- intersect(responder$bcr_patient_barcode,non_responder$bcr_patient_barcode)
responder <- responder[-which(responder[,'bcr_patient_barcode'] %in% interbarcode),]
non_responder <- non_responder[-which(non_responder[,'bcr_patient_barcode'] %in% interbarcode),]
remerg <- merge(responder,non_responder,all = TRUE)
remerg$response <- ifelse(remerg[,'measure_of_response'] %in% c('Complete Response','Partial Response'),"Responder","Non-responder")
wt <- remerg[which(remerg[,'Group'] %in% 'Wild Type'),]
mutant <- remerg[which(remerg[,'Group'] %in% 'Mutant'),]
## Wild Type
### ATF3
with(wt, shapiro.test(ATF3[response == "Responder"]))
with(wt, shapiro.test(ATF3[response == "Non-responder"]))
var.test(ATF3 ~ response, data = wt)
wilcox.test(ATF3 ~ response, data = wt, var.equal = TRUE)
library(ggplot2)
library(ggpubr)
p1 <- ggplot(data = wt,aes(x = response,y = ATF3,fill = response))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  coord_cartesian(ylim = c(0,30))
p1+stat_compare_means(method = 'wilcox.test',label.x = 1, label.y = 30)
### TP53
with(wt, shapiro.test(TP53[response == "Responder"]))
with(wt, shapiro.test(TP53[response == "Non-responder"]))
var.test(TP53 ~ response, data = wt)
wilcox.test(TP53 ~ response, data = wt, var.equal = TRUE)
p1 <- ggplot(data = wt,aes(x = response,y = TP53,fill = response))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  coord_cartesian(ylim = c(0,45))
p1+stat_compare_means(method = 'wilcox.test',label.x = 1, label.y = 45)

## Mutant
### ATF3
with(mutant, shapiro.test(ATF3[response == "Responder"]))
with(mutant, shapiro.test(ATF3[response == "Non-responder"]))
var.test(ATF3 ~ response, data = mutant)
wilcox.test(ATF3 ~ response, data = mutant, var.equal = FALSE)
p1 <- ggplot(data = mutant,aes(x = response,y = ATF3,fill = response))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  coord_cartesian(ylim = c(0,25))
p1+stat_compare_means(method = 'wilcox.test',label.x = 1, label.y = 25)
### TP53
with(mutant, shapiro.test(TP53[response == "Responder"]))
with(mutant, shapiro.test(TP53[response == "Non-responder"]))
var.test(TP53 ~ response, data = mutant)
wilcox.test(TP53 ~ response, data = mutant, var.equal = TRUE)
p1 <- ggplot(data = mutant,aes(x = response,y = TP53,fill = response))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  coord_cartesian(ylim = c(0,75))
p1+stat_compare_means(method = 'wilcox.test',label.x = 1, label.y = 75)

# Non-truncation mutation
View(Mut)
interbarcode <- intersect(Mut$submitter_id,mutant$bcr_patient_barcode)
View(interbarcode)
Mut <- t(Mut) %>% as.data.frame()
Element <- c('Variant_Classification')
Mut <- Mut[Element,interbarcode]
Mut <- t(Mut) %>% as.data.frame()
Mut$bcr_patient_barcode <- rownames(Mut)
View(Mut)
mergdata <- merge(mutant,Mut,by = "bcr_patient_barcode")
nontruncation <- mergdata[which(mergdata[,'Variant_Classification'] %in% c('Missense_Mutation')),]
### TP53
with(nontruncation, shapiro.test(TP53[response == "Responder"]))
with(nontruncation, shapiro.test(TP53[response == "Non-responder"]))
var.test(TP53 ~ response, data = nontruncation)
wilcox.test(TP53 ~ response, data = nontruncation, var.equal = TRUE)
p1 <- ggplot(data = nontruncation,aes(x = response,y = TP53,fill = response))+
  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  coord_cartesian(ylim = c(0,75))
p1+stat_compare_means(method = 'wilcox.test',label.x = 1, label.y = 75)

