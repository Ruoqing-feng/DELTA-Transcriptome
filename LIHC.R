# Loading the package: TCGAbiolinks,maftools
library(TCGAbiolinks)
library(maftools)
# Mutation data
LIHC.mutect.maf <- GDCquery_Maf('LIHC',pipelines = 'mutect2')
laml <- read.maf(maf = LIHC.mutect.maf, isTCGA = TRUE)
snpdata <- laml@data
snpdata <- snpdata[,c("Hugo_Symbol","Variant_Type",
                      "Variant_Classification","Tumor_Sample_Barcode")]
TP53.snp <- subset(snpdata,Hugo_Symbol=='TP53')
silent <- laml@maf.silent
silent <- silent[,c("Hugo_Symbol","Variant_Type",
                    "Variant_Classification","Tumor_Sample_Barcode")]
TP53.silent <- subset(silent,Hugo_Symbol=='TP53')
interbarcode <- intersect(TP53.snp$Tumor_Sample_Barcode,TP53.silent$Tumor_Sample_Barcode)
tidydata <- LIHC.mutect.maf[,c("Hugo_Symbol","Variant_Type",
                               "Variant_Classification","Tumor_Sample_Barcode")]
TP53 <- subset(tidydata,Hugo_Symbol=='TP53')
library(reshape2)
rsdata <- dcast(tidydata, Hugo_Symbol~Tumor_Sample_Barcode)
rownames(rsdata) <- rsdata$Hugo_Symbol
rsdata <- rsdata[,-1]
Gene <- rsdata['TP53',]
names(Gene) <-  sapply(strsplit(names(Gene),'-'),function(x) paste0(x[1:3],collapse="-"))
MAF_barcode <-  sapply(strsplit(names(Gene),'-'),function(x) paste0(x[1:3],collapse="-"))
Gene[,'TCGA-DD-A4NQ'] <- 0
Gene[,'TCGA-G3-AAV4'] <- 0
Gene[,'TCGA-CC-5259'] <- 0

# Gene expression data
TCGAbiolinks:::getProjectSummary("TCGA-LIHC")
query <- GDCquery(project = "TCGA-LIHC", 
                  experimental.strategy = "RNA-Seq",
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - FPKM")
GDCdownload(query)
LIHCRnaseq <- GDCprepare(query)
library("SummarizedExperiment")
LIHC_FPKM <- assay(LIHCRnaseq)
write.csv(LIHC_FPKM, file = 'TCGA-LIHC-Counts.csv')
LIHC_FPKM <- read.csv('TCGA-LIHC-Counts.csv',header = TRUE,sep = ',',check.names = FALSE)
rownames(LIHC_FPKM) <- LIHC_FPKM[,1]
LIHC_FPKM <- LIHC_FPKM[,-1]
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(LIHC_FPKM),typesample = c("TP"))
LIHC_FPKM <- LIHC_FPKM[c('ENSG00000162772','ENSG00000141510'),samplesTP]
row.names(LIHC_FPKM) <- c("ATF3","TP53")
names(LIHC_FPKM) <-  sapply(strsplit(names(LIHC_FPKM),'-'),function(x) paste0(x[1:3],collapse="-"))
barcode <- sapply(strsplit(names(LIHC_FPKM),'-'),function(x) paste0(x[1:3],collapse="-"))

# Clinical data
clinical <- GDCquery_clinic(project = "TCGA-LIHC", type = "clinical")
library(dplyr)
clinical <-clinical[,c("submitter_id","vital_status","days_to_death","days_to_last_follow_up")]
View(clinical)
clinical <- clinical[-13,]
clinical[clinical$vital_status=='Dead',]$vital_status <- 1
clinical[clinical$vital_status=='Alive',]$vital_status <- 0
clinical$vital_status <- as.numeric(clinical$vital_status)
clinical$time <- clinical$days_to_death
clinical$time[which(is.na(clinical$time))] <- clinical$days_to_last_follow_up[which(is.na(clinical$time))]
clinical <- filter(clinical,time != "--")
clinical$time <- as.numeric(clinical$time)

# data fusion: gene expression data, mutation data and clinical data
interbarcode <- intersect(barcode,MAF_barcode)
interbarcode <- intersect(interbarcode,clinical$submitter_id)
Element <- c('TP53')
geneMut <- Gene[Element,interbarcode]
geneMut <- t(geneMut) %>% as.data.frame() %>% data.matrix() %>% as.data.frame()
geneMut$submitter_id <- rownames(geneMut)
library(dplyr)
gene <- c('ATF3','TP53')
geneExp <- LIHC_FPKM[gene,interbarcode]
geneExp <- t(geneExp) %>% as.data.frame() %>% data.matrix() %>% as.data.frame()
geneExp$submitter_id <- rownames(geneExp)
mergdata <- merge(geneExp,geneMut,by = "submitter_id")
mergdata <- merge(clinical,mergdata,by = "submitter_id")
names(mergdata)[7] <- 'TP53'
names(mergdata)[8] <- 'condition'               
# Divide the data into TP53-WT and TP53-mutation
mergdata$Group <- ifelse(mergdata[,'condition']==0,"Wild Type","Mutant")
wt <- mergdata[which(mergdata[,'Group'] %in% 'Wild Type'),]
mutant <- mergdata[which(mergdata[,'Group'] %in% 'Mutant'),]

# Survival analysis
library(survival)
library(survminer)
## cutpoint
### Expression: ATF3 in TP53 wt
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
           legend.labs = c("High(n=116)", "Low(n=140)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### Expression: TP53 in TP53 wt
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
           legend.labs = c("High(n=227)", "Low(n=29)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### Expression: ATF3 in TP53 mutant
sur.cut <- surv_cutpoint(mutant,time= 'time',
                         event = 'vital_status' ,
                         variables = 'ATF3' )
mutant1 <- surv_categorize(sur.cut)
table(mutant1$ATF3)
fit <- survfit(Surv(time, vital_status) ~ ATF3, data = mutant1)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "bold", "darkblue"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'ATF3',
           legend.labs = c("High(n=57)", "Low(n=44)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### Expression: TP53 in TP53 mutant
sur.cut <- surv_cutpoint(mutant,time= 'time',
                         event = 'vital_status' ,
                         variables = 'TP53' )
mutant2 <- surv_categorize(sur.cut)
table(mutant2$TP53)
fit <- survfit(Surv(time, vital_status) ~ TP53, data = mutant2)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "bold", "darkblue"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           legend.labs = c("High(n=44)", "Low(n=57)"),
           title='TP53 Expression in TP53 Mutant',
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)

# Data fusion: gene expression data and clinical data
interbarcode <- intersect(barcode,clinical$submitter_id)
geneExp <- LIHC_FPKM[gene,interbarcode]
geneExp <- t(geneExp) %>% as.data.frame() %>% data.matrix() %>% as.data.frame()
geneExp$submitter_id <- rownames(geneExp)
mergdata <- merge(geneExp,clinical,by = "submitter_id")
# Survival Analysis
## cutpoint
### Expression: ATF3 in LIHC
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
           legend.labs = c("High(n=230)", "Low(n=139)"),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)

### Expression: TP53 in LIHC
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
           legend.labs = c("High(n=293)", "Low(n=76)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)

# Data fusion: mutation data and clinical data
interbarcode <- intersect(MAF_barcode,clinical$submitter_id)
GENE <- 'TP53'
geneMut <- Gene[GENE,interbarcode]
geneMut <- t(geneMut) %>% as.data.frame() %>% data.matrix() %>% as.data.frame()
geneMut$submitter_id <- rownames(geneMut)
mergdata <- merge(clinical,geneMut,by = "submitter_id")
mergdata$Group <- ifelse(mergdata[,'TP53']==0,"Wild Type","Mutant")
table(mergdata$Group)
# Survival analysis: TP53-WT vs TP53-mutation                  
fit <- survfit(Surv(time, vital_status) ~ Group, data=mergdata)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "darkblue", "bold"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           ylim = c(0,1),
           legend.labs = c("Mutant(n=102)", "Wild Type(n=260)"),
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)

# data fusion: gene expression and mutation data
Gene <- t(Gene) %>% as.data.frame()
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
geneExp <- LIHC_FPKM[gene,interbarcode]
geneExp <- t(geneExp) %>% as.data.frame() %>% data.matrix() %>% as.data.frame()
geneExp$submitter_id <- rownames(geneExp)
Mut$submitter_id <- Mut[,'Tumor_Sample_Barcode']
mergdata <- merge(geneExp,Mut,by = "submitter_id")
# TP53 expression distribution
p <- ggplot(data = mergdata,aes(x = Variant_Classification,y = TP53,fill = Variant_Classification))+  geom_boxplot(outlier.shape = NA)+
  theme_bw()+
  coord_cartesian(ylim = c(0,50))
p

# Data fusion: mutation data and clinical data
interbarcode <- intersect(Mut$Tumor_Sample_Barcode,clinical$submitter_id)
Mut <- t(Mut) %>% as.data.frame()
colnames(Mut) <- Mut[3,]
Element <- c('Variant_Classification')
Mut <- Mut[Element,interbarcode]
Mut <- t(Mut) %>% as.data.frame()
Mut$submitter_id <- rownames(Mut)
mergdata <- merge(clinical,Mut,by = "submitter_id")      
data1 <- mergdata[which(mergdata[,'Variant_Classification'] %in% c('WT','Missense_Mutation','In_Frame_Del')),]
data2 <- mergdata[-which(mergdata[,'Variant_Classification'] %in% c('Missense_Mutation','In_Frame_Del','In_Frame_Ins','Splice_Site')),]
data1$Group <- ifelse(data1$Variant_Classification %in% c("WT"),"Wild Type","Non-truncation")
table(data1$Group)
data2$Group <- ifelse(data2$Variant_Classification %in% c("WT"),"Wild Type","Truncation")
table(data2$Group)
fit <- survfit(Surv(time, vital_status) ~ Group, data=data1)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "darkblue", "bold"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           legend.labs = c("Non-truncation(n=62)", "Wild Type(n=260)"),
           ylim = c(0,1),
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
fit <- survfit(Surv(time, vital_status) ~ Group, data=data2)
ggsurvplot(fit, main = "Survival curve",
           font.main = c(10, "darkblue", "bold"),
           font.tickslab = c(12, "plain", "darkgreen"),
           legend.title = 'TP53',
           legend.labs = c("Truncation(n=31)", "Wild Type(n=260)"),
           ylim = c(0,1),
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)

# data fusion: gene expression data, mutation data and clinical data
Mut$Tumor_Sample_Barcode <- rownames(Mut)
interbarcode <- intersect(barcode,Mut$Tumor_Sample_Barcode)
interbarcode <- intersect(interbarcode,clinical$submitter_id)
geneExp <- LIHC_FPKM[gene,interbarcode]
geneExp <- t(geneExp) %>% as.data.frame() %>% data.matrix() %>% as.data.frame()
geneExp$submitter_id <- rownames(geneExp)
Mut <- t(Mut) %>% as.data.frame()
colnames(Mut) <- Mut[3,]
Element <- c('Variant_Classification')
Mut <- Mut[Element,interbarcode]
Mut <- t(Mut) %>% as.data.frame()
Mut$submitter_id <- rownames(Mut)
mergdata <- merge(geneExp,Mut,by = "submitter_id")
mergdata <- merge(clinical,mergdata,by = "submitter_id")
mutant <- mergdata[-which(mergdata[,'Variant_Classification'] %in% c('WT','Splice_Site','In_Frame_Ins')),]
mutant$Group <- ifelse(mutant$Variant_Classification %in% c("Missense_Mutation","In_Frame_Del"),"Non_truncation","Truncation")
truncation <- mutant[which(mutant[,'Group'] %in% c('Truncation')),]
non_truncation <- mutant[which(mutant[,'Group'] %in% c('Non_truncation')),]
### Expression: ATF3 in LIHC-TP53-truncation
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
           legend.labs = c("High(n=14)", "Low(n=16)"),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)

### Expression: TP53 in LIHC-TP53-truncation
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
           legend.labs = c("High(n=12)", "Low(n=18)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)
### Expression: ATF3 in LIHC-TP53-non_truncation
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
           legend.labs = c("High(n=10)", "Low(n=52)"),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)

### Expression: TP53 in LIHC-TP53-non_truncation
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
           legend.labs = c("High(n=40)", "Low(n=22)"),
           ylim = c(0,1),
           xlab = 'Time(day)',
           pval = TRUE,
           censor = FALSE,
           palette = c("#FF0000", "#2200FF"),
           ggtheme = theme_bw() 
)    

