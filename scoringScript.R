library("superheat")
library(dplyr)
library(xlsx)


################################################
#scoring for DE understudied kinases  
################################################
#load the DE genes from each cancer 

IDG_kinase <- read.csv("/projects/ccs/schurerlab/Rimpi/kinases_IDG.txt", sep="\t")
updated_IDG_kinases <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/Understudiedkinases_2019June.csv")

CanType <- c("BLCA",  "KIRP","KIRC", "BRCA", "COAD", "ESCA",  "HNSC" , "LIHC", "LUAD", "LUSC", "STAD", "THCA", "CESC", "GBM" , "PRAD", "CHOL", "KICH", "UCEC" ,  "READ" ,"PCPG")   #"CHOL", "KICH", "UCEC" ,  "READ" ,"PCPG" #


kk <- data.frame()
for(bb in 1:length(CanType) ){
  print(CanType[bb])
#load the DE genes with p-value from each cancer 
load(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[bb], "/AllGeneswithFC.RData", sep =""))  #dataDEGs


# tt <- dataDEGs[grep('\\c', rownames(dataDEGs), value = TRUE),]
# tt1 <- gsub(".*\\((.*)\\).*", "\\1", tt$genes)
# 
# 
# tt3 <- gsub(" ", "|", tt2, fixed=TRUE)
# 
# # tt1 <-  unlist(strsplit(tt$genes,"\\c|\\("))
# tt1 <- sapply(strsplit(tt$genes, ","), "[", 1)
# # 
# tt2 <-  unlist(strsplit(tt1,"\\c|\\("))
# 
# tt3 <- gsub('list|[[:punct:]]', "", tt2)
# 
# tt4 <- grep(pattern = "[a-zA-Z0-9]", tt3, value = TRUE)
# 
# tt$genes1 <- tt2
# rownames(tt) <- tt$genes1
# 
# tt4 <- tt[,c(1:5)]
# 
# dataDEGs1 <- dataDEGs[!grepl('\\c', rownames(dataDEGs)),]
# dataDEGs1 <- rbind(dataDEGs1, tt4)
# 
# dim(dataDEGs1) == dim(dataDEGs)




Underst <- subset(dataDEGs, rownames(dataDEGs)  %in% IDG_kinase$HGNC.Sym)

gg <- dataDEGs[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(dataDEGs)),]
gg
gg <- gg[-3,]
gg$genes <- c("PAK5", "EEF2K", "COQ8B", "UCK2", "NRBP2")
rownames(gg) <- gg$genes
gg
gg <- gg[ ,c(1:5)]


Underst <- rbind(Underst, gg)


library(dplyr)
Underst_adj.P <- select(Underst, FDR)
colnames(Underst_adj.P) <- paste("FDR_",CanType[19], sep="")
Underst_adj.P$Row.names <- rownames(Underst_adj.P)

kk <- merge(as.data.frame(kk), as.data.frame(Underst_adj.P), by='Row.names', all=TRUE)
head(kk)
dim(kk)
#kk <- cbind(kk, Underst_adj.P)

#rownames(kk) == rownames(Underst_adj.P)
}

write.csv(kk, file= "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/DE_Pvalue.csv", row.names = FALSE)







################################################
#scoring for stages in each cancer 
################################################

#read the files

IDG_kinase <- read.csv("/projects/ccs/schurerlab/Rimpi/kinases_IDG.txt", sep="\t")
Dir <- "/projects/ccs/schurerlab/DerekJ/TCGA/ClinicalMerged/"

#no stages "CESC", "GBM" ,"PCPG", "PRAD","UCEC"

CanType <- c("BRCA",  "CHOL", "COAD", "ESCA",  "HNSC" , "LIHC", "LUAD", "LUSC",  "READ", "STAD", "THCA")     #"BLCA", "KICH", "KIRP"


 



#Save files in 

save <- "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/stagewisePvalueForeachCancer/"


################################################
#for loop to call all the files and and calculate the P-value for each stage using 
################################################


Anova_stage <- data.frame()
Anova_grade <- data.frame()
for(i in 1:length(CanType)){
 print(CanType[i] )

  a <- read.delim(paste(Dir, CanType[1], ".Normalize.expression_IDG_clinicalmerged.csv", sep="")) 
  table(a$ajcc_pathologic_tumor_stage)
  
  b <- a[, colnames(a) %in% IDG_kinase$HGNC.Sym]
  c <- select(a, TumorType, bcr_patient_barcode, ajcc_pathologic_tumor_stage, histological_grade)
  c <- cbind(c, b)
  f <- as.character(c$ajcc_pathologic_tumor_stage)
  f[f == "Stage IA" ] <- "Stage I" 
  f[f == "Stage IB" ] <- "Stage I"
  f[f == "Stage IIA" ] <- "Stage II" 
  f[f == "Stage IIB" ] <- "Stage II" 
  f[f == "Stage IIIB" ] <- "Stage III" 
  f[f == "Stage IIIA" ] <- "Stage III" 
  f[f == "Stage IIIC" ] <- "Stage III"
  c$ajcc_pathologic_tumor_stage1 <- f
  
  d <- subset(c, subset = ajcc_pathologic_tumor_stage1 == "Stage I" | ajcc_pathologic_tumor_stage1 == "Stage II" | ajcc_pathologic_tumor_stage1 == "Stage III"| ajcc_pathologic_tumor_stage1 == "Stage IV" )
 e <- subset(c, subset = histological_grade == "G1" | histological_grade == "G2"| histological_grade == "G3" | histological_grade == "G4"  )
genes <- as.character(colnames(b))

for(j in 5:ncol(d)){
  
  print(colnames(d[j]))
  
  
  formula <- as.formula(paste(colnames(d)[j], " ~ ajcc_pathologic_tumor_stage1", sep=""))
  
  res.aov <- aov(formula, data = d)
  #res.aov <- aov(colnames(c[j]) ~ ajcc_pathologic_tumor_stage, data = c)
  
  
  
  #res.aov_grade <- aov(genes[,j] ~ histological_grade, data = c)
  
  
  
  print(paste("Anova for the gene == ", colnames(c[j])))
    print(TukeyHSD(res.aov))
    
    res.aov_result <-  TukeyHSD(res.aov)[1]
    res.aov_result <- data.frame(res.aov_result$ajcc_pathologic_tumor_stage1)
    
    res.aov_result$gene <- colnames(d[j])
    Anova_s <- res.aov_result
    
    Anova_stage <- rbind(Anova_stage, Anova_s)
    
    
    for(k in 5:ncol(e)){
    #for the grade 
    formula1 <- as.formula(paste(colnames(e)[k], " ~ histological_grade", sep=""))
    res.aov_grade <- aov(formula1, data = e)
  
    res.aov_grade_result <-  TukeyHSD(res.aov_grade)[1]
    res.aov_grade_result <-  data.frame(res.aov_grade_result$histological_grade)
    
    res.aov_grade_result$gene <- colnames(e[k]) 
    Anova_g <- res.aov_grade_result
    Anova_grade <- rbind(Anova_grade, Anova_g)
    
    }
                                        
    
}

write.csv(Anova_stage, file = paste(save, CanType[11], "Stage.csv"))
length(unique(Anova_stage$gene)) == length(genes)
        
write.csv(Anova_grade, file = paste(save, CanType[10], "Grade.csv"))
length(unique(Anova_grade$gene)) == length(genes)
  
}

  
################################################
#survival analysis of all the understudied kinases in each cancer


colnames(BRCA.surv_rnaseq.cat) <- c("times", "patient.vital_status", "expr", "cohart")
BRCA.surv_rnaseq.cat$gene <- "WNK2"

################################################
library(survival)
library(RTCGA.rnaseq)
library(RTCGA.clinical)
library(dplyr)
library(tidyr)

CanType <- CanType <- c("BLCA",  "KIRP", "BRCA", "COAD", "ESCA",  "HNSC" , "LIHC", "LUAD", "LUSC", "STAD", "THCA", "CESC", "GBM" , "PRAD", "CHOL", "KICH", "UCEC" ,  "READ" ,"PCPG")   #"CHOL", "KICH", "UCEC" ,  "READ" ,"PCPG" #
understudied <- as.character(IDG_kinase$HGNC.Sym)
understudied <- understudied[-c(75,79)]  #PAK7='PAK7', PAK6='PAK6',



#BLCA
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-BLCA/TCGA-BLCA_clinical.rds")


#BLCA
BLCA.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/BLCA_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#BLCA
survivalTCGA(BLCA.clin) -> BLCA.clinical.surv
dd <- BLCA.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(BLCA.Normalize.expression)),]  #COQ8B changed to ADCK4
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
BLCA.Normalize.expression <- BLCA.Normalize.expression[-grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(BLCA.Normalize.expression)),]
BLCA.Normalize.expression <- rbind(BLCA.Normalize.expression, dd)

BLCA.Normalized <- data.frame(t(BLCA.Normalize.expression))
BLCA.Normalized <- cbind(bcr_patient_barcode = rownames(BLCA.Normalized), BLCA.Normalized)
rownames(BLCA.Normalized) <- NULL
BLCA.rnaseq <- BLCA.Normalized



expressionsTCGA(
 BLCA.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> BLCA.rnaseq11  



BLCA.clinical.surv %>%
  left_join(BLCA.rnaseq11,
            by = "bcr_patient_barcode") ->
  BLCA.surv_rnaseq


table(BLCA.surv_rnaseq$cohort, useNA = "always")


BLCA.surv_rnaseq <- BLCA.surv_rnaseq %>%
  filter(!is.na(cohort))



BLCA.surv_rnaseq<- BLCA.surv_rnaseq[!duplicated(BLCA.surv_rnaseq$bcr_patient_barcode),]
dim(BLCA.surv_rnaseq)
######################
#GENES for survival analysis
#######################
CDKL4
library(survminer)

gen <- as.character(colnames(BLCA.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
# geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#            "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#            "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#            "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

# dataSurvTable <- data.frame()
# dataSurvplot <- data.frame()
# pvalue <- data.frame()
# fittable <- data.frame()
# surmedian <- data.frame()
BLCA.surv_rnaseq <- BLCA.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
 BLCA.surv_rnaseq.cut <- surv_cutpoint(
    BLCA.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(BLCA.surv_rnaseq.cut)
  
  #distri <- plot(LIHC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  BLCA.surv_rnaseq.cat <- surv_categorize(BLCA.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = BLCA.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  surv_median(fit)
 
  p <- rbind(p, pval) 
  
  save(pval,fit,BLCA.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-BLCA/survival/", geneid[i], "BLCA_Pvalue.RData", sep ="") )
  
  png(file = paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-BLCA/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= BLCA.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
    
  )  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/BLCA_Pvalue.csv")
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/BLCA_Pvalue.csv") 
  
 

########
#KIRP
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRP/TCGA-KIRP_clinical.rds")


#LIHC
KIRP.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/KIRP_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#KIRP
survivalTCGA(KIRP.clin) -> KIRP.clinical.surv
dd <- KIRP.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(KIRP.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
KIRP.Normalize.expression <- rbind(KIRP.Normalize.expression, dd)

KIRP.Normalized <- data.frame(t(KIRP.Normalize.expression))
KIRP.Normalized <- cbind(bcr_patient_barcode = rownames(KIRP.Normalized), KIRP.Normalized)
rownames(KIRP.Normalized) <- NULL
KIRP.rnaseq <- KIRP.Normalized



expressionsTCGA(
  KIRP.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K ='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> KIRP.rnaseq1  



KIRP.clinical.surv %>%
  left_join(KIRP.rnaseq1,
            by = "bcr_patient_barcode") ->
  KIRP.surv_rnaseq


table(KIRP.surv_rnaseq$cohort, useNA = "always")


KIRP.surv_rnaseq <- KIRP.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(KIRP.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(KIRP.surv_rnaseq))
head(gen)
gen <- gen[-c(1:4)]

geneid <- gen
# #geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#            "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#            "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#            "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

KIRP.surv_rnaseq<- KIRP.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  KIRP.surv_rnaseq.cut <- surv_cutpoint(
    KIRP.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(KIRP.surv_rnaseq.cut)
  
  #distri <- plot(LIHC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  KIRP.surv_rnaseq.cat <- surv_categorize(KIRP.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = KIRP.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,KIRP.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRP/survival/", geneid[i], "KIRP_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRP/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= KIRP.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/KIRP_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/KIRP_Pvalue.RData") 

############

#BRCA
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-BRCA/TCGA-BRCA_clinical.rds")


#LIHC
BRCA.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/BRCA_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#BRCA
survivalTCGA(BRCA.clin) -> BRCA.clinical.surv
dd <- BRCA.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(BRCA.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
BRCA.Normalize.expression <- rbind(BRCA.Normalize.expression, dd)

BRCA.Normalized <- data.frame(t(BRCA.Normalize.expression))
BRCA.Normalized <- cbind(bcr_patient_barcode = rownames(BRCA.Normalized), BRCA.Normalized)
rownames(BRCA.Normalized) <- NULL
BRCA.rnaseq <- BRCA.Normalized



expressionsTCGA(
  BRCA.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> BRCA.rnaseq1  



BRCA.clinical.surv %>%
  left_join(BRCA.rnaseq1,
            by = "bcr_patient_barcode") ->
  BRCA.surv_rnaseq


table(BRCA.surv_rnaseq$cohort, useNA = "always")


BRCA.surv_rnaseq <- BRCA.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(BRCA.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(BRCA.surv_rnaseq))
colnames(gen)
gen <- gen[-c(1:4)]

geneid <- gen
# #geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#            "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#            "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#            "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

BRCA.surv_rnaseq<- BRCA.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  BRCA.surv_rnaseq.cut <- surv_cutpoint(
    BRCA.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(BRCA.surv_rnaseq.cut)
  
  #distri <- plot(LIHC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  BRCA.surv_rnaseq.cat <- surv_categorize(BRCA.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = BRCA.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  
  save(pval,fit,BRCA.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-BRCA/survival/", geneid[i], "C", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-BRCA/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= BRCA.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
    
    
  )
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/BRCA_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/BRCA_Pvalue.RData") 


##############
#


#COAD
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-COAD/TCGA-COAD_clinical.rds")


#LIHC
COAD.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/COAD_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#COAD
survivalTCGA(COAD.clin) -> COAD.clinical.surv
dd <- COAD.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(COAD.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
COAD.Normalize.expression <- rbind(COAD.Normalize.expression, dd)

COAD.Normalized <- data.frame(t(COAD.Normalize.expression))
COAD.Normalized <- cbind(bcr_patient_barcode = rownames(COAD.Normalized), COAD.Normalized)
rownames(COAD.Normalized) <- NULL
COAD.rnaseq <- COAD.Normalized



expressionsTCGA(
  COAD.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> COAD.rnaseq1  



COAD.clinical.surv %>%
  left_join(COAD.rnaseq1,
            by = "bcr_patient_barcode") ->
  COAD.surv_rnaseq


table(COAD.surv_rnaseq$cohort, useNA = "always")


COAD.surv_rnaseq <- COAD.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(COAD.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(COAD.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
# #geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
# "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
# "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
# "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

COAD.surv_rnaseq<- COAD.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  COAD.surv_rnaseq.cut <- surv_cutpoint(
    COAD.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(COAD.surv_rnaseq.cut)
  
  #distri <- plot(LIHC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  COAD.surv_rnaseq.cat <- surv_categorize(COAD.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = COAD.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  
  save(pval,fit,COAD.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-COAD/survival/", geneid[i], "COAD_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-COAD/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= COAD.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/COAD_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/COAD_Pvalue.RData") 



#######

#ESCA
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-ESCA/TCGA-ESCA_clinical.rds")


#LIHC
ESCA.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/ESCA_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#ESCA
survivalTCGA(ESCA.clin) -> ESCA.clinical.surv
dd <- ESCA.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(ESCA.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
ESCA.Normalize.expression <- rbind(ESCA.Normalize.expression, dd)

ESCA.Normalized <- data.frame(t(ESCA.Normalize.expression))
ESCA.Normalized <- cbind(bcr_patient_barcode = rownames(ESCA.Normalized), ESCA.Normalized)
rownames(ESCA.Normalized) <- NULL
ESCA.rnaseq <- ESCA.Normalized



expressionsTCGA(
  ESCA.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> ESCA.rnaseq1  



ESCA.clinical.surv %>%
  left_join(ESCA.rnaseq1,
            by = "bcr_patient_barcode") ->
  ESCA.surv_rnaseq


table(ESCA.surv_rnaseq$cohort, useNA = "always")


ESCA.surv_rnaseq <- ESCA.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(ESCA.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(ESCA.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
# "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
# "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
# "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

ESCA.surv_rnaseq<- ESCA.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  ESCA.surv_rnaseq.cut <- surv_cutpoint(
    ESCA.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(ESCA.surv_rnaseq.cut)
  
  #distri <- plot(LIHC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  ESCA.surv_rnaseq.cat <- surv_categorize(ESCA.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = ESCA.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  
  save(pval,fit,ESCA.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-ESCA/survival/", geneid[i], "ESCA_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-ESCA/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= ESCA.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/ESCA_Pvalue.csv" , row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/ESCA_Pvalue.RData") 



############
#HNSC
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-HNSC/TCGA-HNSC_clinical.rds")


#LIHC
HNSC.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/HNSC_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#HNSC
survivalTCGA(HNSC.clin) -> HNSC.clinical.surv
dd <- HNSC.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(HNSC.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
HNSC.Normalize.expression <- rbind(HNSC.Normalize.expression, dd)

HNSC.Normalized <- data.frame(t(HNSC.Normalize.expression))
HNSC.Normalized <- cbind(bcr_patient_barcode = rownames(HNSC.Normalized), HNSC.Normalized)
rownames(HNSC.Normalized) <- NULL
HNSC.rnaseq <- HNSC.Normalized



expressionsTCGA(
  HNSC.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> HNSC.rnaseq1  



HNSC.clinical.surv %>%
  left_join(HNSC.rnaseq1,
            by = "bcr_patient_barcode") ->
  HNSC.surv_rnaseq


table(HNSC.surv_rnaseq$cohort, useNA = "always")


HNSC.surv_rnaseq <- HNSC.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(HNSC.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(HNSC.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
# "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
# "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
# "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

HNSC.surv_rnaseq<- HNSC.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  HNSC.surv_rnaseq.cut <- surv_cutpoint(
    HNSC.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(HNSC.surv_rnaseq.cut)
  
  #distri <- plot(LIHC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  HNSC.surv_rnaseq.cat <- surv_categorize(HNSC.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = HNSC.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,HNSC.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-HNSC/survival/", geneid[i], "HNSC_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-HNSC/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= HNSC.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/HNSC_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/HNSC_Pvalue.RData") 



#########
#LIHC
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/TCGA-LIHC_clinical.rds")


#LIHC
LIHC.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LIHC_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#LIHC
survivalTCGA(LIHC.clin) -> LIHC.clinical.surv
dd <- LIHC.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(LIHC.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
LIHC.Normalize.expression <- rbind(LIHC.Normalize.expression, dd)

LIHC.Normalized <- data.frame(t(LIHC.Normalize.expression))
LIHC.Normalized <- cbind(bcr_patient_barcode = rownames(LIHC.Normalized), LIHC.Normalized)
rownames(LIHC.Normalized) <- NULL
LIHC.rnaseq <- LIHC.Normalized



expressionsTCGA(
  LIHC.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> LIHC.rnaseq1  



LIHC.clinical.surv %>%
  left_join(LIHC.rnaseq1,
            by = "bcr_patient_barcode") ->
  LIHC.surv_rnaseq


table(LIHC.surv_rnaseq$cohort, useNA = "always")


LIHC.surv_rnaseq <- LIHC.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(LIHC.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(LIHC.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
# "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
# "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
# "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

LIHC.surv_rnaseq<- LIHC.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  LIHC.surv_rnaseq.cut <- surv_cutpoint(
    LIHC.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(LIHC.surv_rnaseq.cut)
  
  #distri <- plot(LIHC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  LIHC.surv_rnaseq.cat <- surv_categorize(LIHC.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = LIHC.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,LIHC.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], "LIHC_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= LIHC.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/LIHC_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/LIHC_Pvalue.RData") 



##############

#LUAD
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUAD/TCGA-LUAD_clinical.rds")


#LUAD
LUAD.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LUAD_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#LUAD
survivalTCGA(LUAD.clin) -> LUAD.clinical.surv
dd <- LUAD.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(LUAD.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
LUAD.Normalize.expression <- rbind(LUAD.Normalize.expression, dd)

LUAD.Normalized <- data.frame(t(LUAD.Normalize.expression))
LUAD.Normalized <- cbind(bcr_patient_barcode = rownames(LUAD.Normalized), LUAD.Normalized)
rownames(LUAD.Normalized) <- NULL
LUAD.rnaseq <- LUAD.Normalized



expressionsTCGA(
  LUAD.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> LUAD.rnaseq1  



LUAD.clinical.surv %>%
  left_join(LUAD.rnaseq1,
            by = "bcr_patient_barcode") ->
  LUAD.surv_rnaseq


table(LUAD.surv_rnaseq$cohort, useNA = "always")


LUAD.surv_rnaseq <- LUAD.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(LUAD.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(LUAD.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
# #geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
# "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
# "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
# "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

LUAD.surv_rnaseq<- LUAD.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  LUAD.surv_rnaseq.cut <- surv_cutpoint(
    LUAD.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(LUAD.surv_rnaseq.cut)
  
  #distri <- plot(LUAD.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUAD/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  LUAD.surv_rnaseq.cat <- surv_categorize(LUAD.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = LUAD.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,LUAD.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUAD/survival/", geneid[i], "LUAD_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUAD/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= LUAD.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/LUAD_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/LUAD_Pvalue.RData") 


#######
#LUSC
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUSC/TCGA-LUSC_clinical.rds")


#LUSC
LUSC.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LUSC_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#LUSC
survivalTCGA(LUSC.clin) -> LUSC.clinical.surv
dd <- LUSC.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(LUSC.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
LUSC.Normalize.expression <- rbind(LUSC.Normalize.expression, dd)

LUSC.Normalized <- data.frame(t(LUSC.Normalize.expression))
LUSC.Normalized <- cbind(bcr_patient_barcode = rownames(LUSC.Normalized), LUSC.Normalized)
rownames(LUSC.Normalized) <- NULL
LUSC.rnaseq <- LUSC.Normalized



expressionsTCGA(
  LUSC.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> LUSC.rnaseq1  



LUSC.clinical.surv %>%
  left_join(LUSC.rnaseq1,
            by = "bcr_patient_barcode") ->
  LUSC.surv_rnaseq


table(LUSC.surv_rnaseq$cohort, useNA = "always")


LUSC.surv_rnaseq <- LUSC.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(LUSC.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(LUSC.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
# #geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
# "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
# "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
# "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

LUSC.surv_rnaseq<- LUSC.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  LUSC.surv_rnaseq.cut <- surv_cutpoint(
    LUSC.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(LUSC.surv_rnaseq.cut)
  
  #distri <- plot(LUSC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUSC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  LUSC.surv_rnaseq.cat <- surv_categorize(LUSC.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = LUSC.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,LUSC.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUSC/survival/", geneid[i], "LUSC_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUSC/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= LUSC.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/LUSC_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/LUSC_Pvalue.RData") 



###############
#STAD
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-STAD/TCGA-STAD_clinical.rds")


#STAD
STAD.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/STAD_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#STAD
survivalTCGA(STAD.clin) -> STAD.clinical.surv
dd <- STAD.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(STAD.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
STAD.Normalize.expression <- rbind(STAD.Normalize.expression, dd)

STAD.Normalized <- data.frame(t(STAD.Normalize.expression))
STAD.Normalized <- cbind(bcr_patient_barcode = rownames(STAD.Normalized), STAD.Normalized)
rownames(STAD.Normalized) <- NULL
STAD.rnaseq <- STAD.Normalized



expressionsTCGA(
  STAD.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> STAD.rnaseq1  



STAD.clinical.surv %>%
  left_join(STAD.rnaseq1,
            by = "bcr_patient_barcode") ->
  STAD.surv_rnaseq


table(STAD.surv_rnaseq$cohort, useNA = "always")


STAD.surv_rnaseq <- STAD.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(STAD.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(STAD.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#           "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#          "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#          "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

STAD.surv_rnaseq<- STAD.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  STAD.surv_rnaseq.cut <- surv_cutpoint(
    STAD.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(STAD.surv_rnaseq.cut)
  
  #distri <- plot(STAD.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-STAD/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  STAD.surv_rnaseq.cat <- surv_categorize(STAD.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = STAD.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,STAD.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-STAD/survival/", geneid[i], "STAD_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-STAD/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= STAD.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/STAD_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/STAD_Pvalue.RData") 


##############
#THCA
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-THCA/TCGA-THCA_clinical.rds")


#THCA
THCA.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/THCA_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#THCA
survivalTCGA(THCA.clin) -> THCA.clinical.surv
dd <- THCA.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(THCA.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
THCA.Normalize.expression <- rbind(THCA.Normalize.expression, dd)

THCA.Normalized <- data.frame(t(THCA.Normalize.expression))
THCA.Normalized <- cbind(bcr_patient_barcode = rownames(THCA.Normalized), THCA.Normalized)
rownames(THCA.Normalized) <- NULL
THCA.rnaseq <- THCA.Normalized



expressionsTCGA(
  THCA.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> THCA.rnaseq1  



THCA.clinical.surv %>%
  left_join(THCA.rnaseq1,
            by = "bcr_patient_barcode") ->
  THCA.surv_rnaseq


table(THCA.surv_rnaseq$cohort, useNA = "always")


THCA.surv_rnaseq <- THCA.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(THCA.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(THCA.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#           "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#          "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#          "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

THCA.surv_rnaseq<- THCA.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  THCA.surv_rnaseq.cut <- surv_cutpoint(
    THCA.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(THCA.surv_rnaseq.cut)
  
  #distri <- plot(THCA.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-THCA/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  THCA.surv_rnaseq.cat <- surv_categorize(THCA.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = THCA.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,THCA.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-THCA/survival/", geneid[i], "THCA_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-THCA/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= THCA.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/THCA_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/THCA_Pvalue.RData") 


########################
########################

#CESC
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CESC/TCGA-CESC_clinical.rds")


#CESC
CESC.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/CESC_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#CESC
survivalTCGA(CESC.clin) -> CESC.clinical.surv
dd <- CESC.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(CESC.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
CESC.Normalize.expression <- rbind(CESC.Normalize.expression, dd)

CESC.Normalized <- data.frame(t(CESC.Normalize.expression))
CESC.Normalized <- cbind(bcr_patient_barcode = rownames(CESC.Normalized), CESC.Normalized)
rownames(CESC.Normalized) <- NULL
CESC.rnaseq <- CESC.Normalized



expressionsTCGA(
  CESC.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> CESC.rnaseq1  



CESC.clinical.surv %>%
  left_join(CESC.rnaseq1,
            by = "bcr_patient_barcode") ->
  CESC.surv_rnaseq


table(CESC.surv_rnaseq$cohort, useNA = "always")


CESC.surv_rnaseq <- CESC.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(CESC.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(CESC.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#           "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#          "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#          "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

CESC.surv_rnaseq<- CESC.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  CESC.surv_rnaseq.cut <- surv_cutpoint(
    CESC.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(CESC.surv_rnaseq.cut)
  
  #distri <- plot(CESC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CESC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  CESC.surv_rnaseq.cat <- surv_categorize(CESC.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = CESC.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,CESC.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CESC/survival/", geneid[i], "CESC_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CESC/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= CESC.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/CESC_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/CESC_Pvalue.RData") 




########################
########################

#GBM
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-GBM/TCGA-GBM_clinical.rds")


#GBM
GBM.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                       patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                       patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/GBM_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#GBM
survivalTCGA(GBM.clin) -> GBM.clinical.surv
dd <- GBM.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(GBM.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
GBM.Normalize.expression <- rbind(GBM.Normalize.expression, dd)

GBM.Normalized <- data.frame(t(GBM.Normalize.expression))
GBM.Normalized <- cbind(bcr_patient_barcode = rownames(GBM.Normalized), GBM.Normalized)
rownames(GBM.Normalized) <- NULL
GBM.rnaseq <- GBM.Normalized



expressionsTCGA(
  GBM.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> GBM.rnaseq1  



GBM.clinical.surv %>%
  left_join(GBM.rnaseq1,
            by = "bcr_patient_barcode") ->
  GBM.surv_rnaseq


table(GBM.surv_rnaseq$cohort, useNA = "always")


GBM.surv_rnaseq <- GBM.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(GBM.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(GBM.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#           "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#          "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#          "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

GBM.surv_rnaseq<- GBM.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  GBM.surv_rnaseq.cut <- surv_cutpoint(
    GBM.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(GBM.surv_rnaseq.cut)
  
  #distri <- plot(GBM.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-GBM/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  GBM.surv_rnaseq.cat <- surv_categorize(GBM.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = GBM.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,GBM.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-GBM/survival/", geneid[i], "GBM_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-GBM/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= GBM.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/GBM_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/GBM_Pvalue.RData") 



########################
########################

#PRAD
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PRAD/TCGA-PRAD_clinical.rds")


#PRAD
PRAD.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/PRAD_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#PRAD
survivalTCGA(PRAD.clin) -> PRAD.clinical.surv
dd <- PRAD.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(PRAD.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
PRAD.Normalize.expression <- rbind(PRAD.Normalize.expression, dd)

PRAD.Normalized <- data.frame(t(PRAD.Normalize.expression))
PRAD.Normalized <- cbind(bcr_patient_barcode = rownames(PRAD.Normalized), PRAD.Normalized)
rownames(PRAD.Normalized) <- NULL
PRAD.rnaseq <- PRAD.Normalized



expressionsTCGA(
  PRAD.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> PRAD.rnaseq1  



PRAD.clinical.surv %>%
  left_join(PRAD.rnaseq1,
            by = "bcr_patient_barcode") ->
  PRAD.surv_rnaseq


table(PRAD.surv_rnaseq$cohort, useNA = "always")


PRAD.surv_rnaseq <- PRAD.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(PRAD.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(PRAD.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#           "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#          "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#          "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

PRAD.surv_rnaseq<- PRAD.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  PRAD.surv_rnaseq.cut <- surv_cutpoint(
    PRAD.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(PRAD.surv_rnaseq.cut)
  
  #distri <- plot(PRAD.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PRAD/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  PRAD.surv_rnaseq.cat <- surv_categorize(PRAD.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = PRAD.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval)
  save(pval,fit,PRAD.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PRAD/survival/", geneid[i], "PRAD_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PRAD/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= PRAD.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/PRAD_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/PRAD_Pvalue.RData") 



########################
########################

#CHOL
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CHOL/TCGA-CHOL_clinical.rds")


#CHOL
CHOL.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/CHOL_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#CHOL
survivalTCGA(CHOL.clin) -> CHOL.clinical.surv
dd <- CHOL.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(CHOL.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
CHOL.Normalize.expression <- rbind(CHOL.Normalize.expression, dd)

CHOL.Normalized <- data.frame(t(CHOL.Normalize.expression))
CHOL.Normalized <- cbind(bcr_patient_barcode = rownames(CHOL.Normalized), CHOL.Normalized)
rownames(CHOL.Normalized) <- NULL
CHOL.rnaseq <- CHOL.Normalized



expressionsTCGA(
  CHOL.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> CHOL.rnaseq1  



CHOL.clinical.surv %>%
  left_join(CHOL.rnaseq1,
            by = "bcr_patient_barcode") ->
  CHOL.surv_rnaseq


table(CHOL.surv_rnaseq$cohort, useNA = "always")


CHOL.surv_rnaseq <- CHOL.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(CHOL.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(CHOL.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#           "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#          "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#          "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

CHOL.surv_rnaseq<- CHOL.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  CHOL.surv_rnaseq.cut <- surv_cutpoint(
    CHOL.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(CHOL.surv_rnaseq.cut)
  
  #distri <- plot(CHOL.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CHOL/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  CHOL.surv_rnaseq.cat <- surv_categorize(CHOL.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = CHOL.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval)
  save(pval,fit,CHOL.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CHOL/survival/", geneid[i], "CHOL_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CHOL/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= CHOL.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/CHOL_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/CHOL_Pvalue.RData") 



########################
########################

#KICH
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KICH/TCGA-KICH_clinical.rds")


#KICH
KICH.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/KICH_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#KICH
survivalTCGA(KICH.clin) -> KICH.clinical.surv
dd <- KICH.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(KICH.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
KICH.Normalize.expression <- rbind(KICH.Normalize.expression, dd)

KICH.Normalized <- data.frame(t(KICH.Normalize.expression))
KICH.Normalized <- cbind(bcr_patient_barcode = rownames(KICH.Normalized), KICH.Normalized)
rownames(KICH.Normalized) <- NULL
KICH.rnaseq <- KICH.Normalized



expressionsTCGA(
  KICH.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> KICH.rnaseq1  



KICH.clinical.surv %>%
  left_join(KICH.rnaseq1,
            by = "bcr_patient_barcode") ->
  KICH.surv_rnaseq


table(KICH.surv_rnaseq$cohort, useNA = "always")


KICH.surv_rnaseq <- KICH.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(KICH.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(KICH.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#           "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#          "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#          "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

KICH.surv_rnaseq<- KICH.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  KICH.surv_rnaseq.cut <- surv_cutpoint(
    KICH.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(KICH.surv_rnaseq.cut)
  
  #distri <- plot(KICH.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KICH/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  KICH.surv_rnaseq.cat <- surv_categorize(KICH.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = KICH.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,KICH.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KICH/survival/", geneid[i], "KICH_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KICH/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= KICH.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/KICH_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/KICH_Pvalue.RData") 



########################
########################

#UCEC
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-UCEC/TCGA-UCEC_clinical.rds")


#UCEC
UCEC.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/UCEC_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#UCEC
survivalTCGA(UCEC.clin) -> UCEC.clinical.surv
dd <- UCEC.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(UCEC.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
UCEC.Normalize.expression <- rbind(UCEC.Normalize.expression, dd)

UCEC.Normalized <- data.frame(t(UCEC.Normalize.expression))
UCEC.Normalized <- cbind(bcr_patient_barcode = rownames(UCEC.Normalized), UCEC.Normalized)
rownames(UCEC.Normalized) <- NULL
UCEC.rnaseq <- UCEC.Normalized



expressionsTCGA(
  UCEC.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> UCEC.rnaseq1  



UCEC.clinical.surv %>%
  left_join(UCEC.rnaseq1,
            by = "bcr_patient_barcode") ->
  UCEC.surv_rnaseq


table(UCEC.surv_rnaseq$cohort, useNA = "always")


UCEC.surv_rnaseq <- UCEC.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(UCEC.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(UCEC.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#           "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#          "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#          "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

UCEC.surv_rnaseq<- UCEC.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  UCEC.surv_rnaseq.cut <- surv_cutpoint(
    UCEC.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(UCEC.surv_rnaseq.cut)
  
  #distri <- plot(UCEC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-UCEC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  UCEC.surv_rnaseq.cat <- surv_categorize(UCEC.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = UCEC.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,UCEC.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-UCEC/survival/", geneid[i], "UCEC_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-UCEC/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= UCEC.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/UCEC_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/UCEC_Pvalue.RData") 


########################
########################

#READ
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-READ/TCGA-READ_clinical.rds")


#READ
READ.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/READ_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#READ
survivalTCGA(READ.clin) -> READ.clinical.surv
dd <- READ.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(READ.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
READ.Normalize.expression <- rbind(READ.Normalize.expression, dd)

READ.Normalized <- data.frame(t(READ.Normalize.expression))
READ.Normalized <- cbind(bcr_patient_barcode = rownames(READ.Normalized), READ.Normalized)
rownames(READ.Normalized) <- NULL
READ.rnaseq <- READ.Normalized



expressionsTCGA(
  READ.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> READ.rnaseq1  



READ.clinical.surv %>%
  left_join(READ.rnaseq1,
            by = "bcr_patient_barcode") ->
  READ.surv_rnaseq


table(READ.surv_rnaseq$cohort, useNA = "always")


READ.surv_rnaseq <- READ.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(READ.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(READ.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#           "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#          "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#          "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

READ.surv_rnaseq<- READ.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  READ.surv_rnaseq.cut <- surv_cutpoint(
    READ.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(READ.surv_rnaseq.cut)
  
  #distri <- plot(READ.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-READ/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  READ.surv_rnaseq.cat <- surv_categorize(READ.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = READ.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,READ.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-READ/survival/", geneid[i], "READ_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-READ/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= READ.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/READ_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/READ_Pvalue.RData") 




########################
########################

#PCPG
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PCPG/TCGA-PCPG_clinical.rds")


#PCPG
PCPG.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/PCPG_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#PCPG
survivalTCGA(PCPG.clin) -> PCPG.clinical.surv
dd <- PCPG.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(PCPG.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
PCPG.Normalize.expression <- rbind(PCPG.Normalize.expression, dd)

PCPG.Normalized <- data.frame(t(PCPG.Normalize.expression))
PCPG.Normalized <- cbind(bcr_patient_barcode = rownames(PCPG.Normalized), PCPG.Normalized)
rownames(PCPG.Normalized) <- NULL
PCPG.rnaseq <- PCPG.Normalized



expressionsTCGA(
  PCPG.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> PCPG.rnaseq1  



PCPG.clinical.surv %>%
  left_join(PCPG.rnaseq1,
            by = "bcr_patient_barcode") ->
  PCPG.surv_rnaseq


table(PCPG.surv_rnaseq$cohort, useNA = "always")


PCPG.surv_rnaseq <- PCPG.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(PCPG.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(PCPG.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#           "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#          "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#          "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

PCPG.surv_rnaseq<- PCPG.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  PCPG.surv_rnaseq.cut <- surv_cutpoint(
    PCPG.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(PCPG.surv_rnaseq.cut)
  
  #distri <- plot(PCPG.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PCPG/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  PCPG.surv_rnaseq.cat <- surv_categorize(PCPG.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = PCPG.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,PCPG.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PCPG/survival/", geneid[i], "PCPG_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PCPG/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= PCPG.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/PCPG_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/PCPG_Pvalue.RData") 



#############
#KIRC
clinical <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/TCGA-KIRC_clinical.rds")


#KIRC
KIRC.clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)


load("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/KIRC_Normalize_filter.RData")



##################
#Survival 
############

#gdc_cases.samples.portions.analytes.aliquots.submitter_id

#KIRC
survivalTCGA(KIRC.clin) -> KIRC.clinical.surv
dd <- KIRC.Normalize.expression[grep("NRBP2|UCK2|EEF2K|COQ8B|PAK5", rownames(KIRC.Normalize.expression)),]
rownames(dd)
dd <- dd[-3,]
rownames(dd) <- c("PAK5", "EEF2K", "ADCK4", "UCK2", "NRBP2")
KIRC.Normalize.expression <- rbind(KIRC.Normalize.expression, dd)

KIRC.Normalized <- data.frame(t(KIRC.Normalize.expression))
KIRC.Normalized <- cbind(bcr_patient_barcode = rownames(KIRC.Normalized), KIRC.Normalized)
rownames(KIRC.Normalized) <- NULL
KIRC.rnaseq <- KIRC.Normalized



expressionsTCGA(
  KIRC.rnaseq,
  extract.cols = understudied) %>%
  rename(cohort = dataset,
         CAMKK1 = 'CAMKK1', CAMK1G = 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223', PRKAG1 = 'PRKAG1', PRKAB1='PRKAB1', ADCK2='ADCK2',  ADCK4= 'ADCK4', ALPK2= 'ALPK2', BCKDK='BCKDK', CDKL3='CDKL3', 
         CDKL4 = 'CDKL4', CDKL1=  'CDKL1', CDKL2 = 'CDKL2',CDK15 = 'CDK15',CDK17 = 'CDK17',    CDK11B = 'CDK11B',   CDK14 ='CDK14', CLK4= 'CLK4', 
         CDK10='CDK10', CLK3='CLK3', CSNK2A2='CSNK2A2', DGKH='DGKH', DYRK3='DYRK3', DYRK1B='DYRK1B', DYRK4='DYRK4', EEF2K='EEF2K', FASTK='FASTK', 
         GK2='GK2', HIPK4='HIPK4', CSNK1G3='CSNK1G3', CSNK1G2 ='CSNK1G2', PRKACG='PRKACG', CSNK1G1='CSNK1G1', ITPK1='ITPK1', PSKH2='PSKH2',
         CAMK1D ='CAMK1D', TK2='TK2', PRKACB='PRKACB', PSKH1='PSKH1', PHKA1='PHKA1', RPS6KC1='RPS6KC1', PRKCQ='PRKCQ', MARK4='MARK4', MAST3='MAST3',
         LRRK1='LRRK1' ,MAP3K14= 'MAP3K14', LTK='LTK', MAST4='MAST4', MAST2='MAST2', MAPK4='MAPK4', MKNK2='MKNK2', CDC42BPB='CDC42BPB', NEK11= 'NEK11', 
         NEK5='NEK5', NEK10='NEK10', NIM1K='NIM1K', NEK4='NEK4', NEK6='NEK6', NRK='NRK' ,NEK7 = 'NEK7', SCYL1='SCYL1', NRBP2 ='NRBP2', PIK3C2G ='PIK3C2G',
         PANK3= 'PANK3', PIK3C2B='PIK3C2B',  SCYL3='SCYL3', PDIK1L='PDIK1L', PIP4K2C='PIP4K2C', PIP5K1A='PIP5K1A', PHKG1='PHKG1',
         PHKG2='PHKG2' ,PI4KA = 'PI4KA', PIP5K1B='PIP5K1B', PKN3='PKN3',PRPF4B ='PRPF4B', TP53RK='TP53RK', PRKRA='PRKRA', PRKY='PRKY', PXK='PXK', 
         RIOK3='RIOK3',SBK2 ='SBK2', SCYL2='SCYL2', SBK3='SBK3', POMK='POMK', STK40='STK40' ,STK32A = 'STK32A', STK38L='STK38L', STKLD1='STKLD1', 
         TESK2='TESK2',TESK1 ='TESK1', TBCK='TBCK', TLK1='TLK1', TPK1='TPK1', TLK2='TLK2', TSSK1B='TSSK1B', TTBK2='TTBK2', TSSK4='TSSK4', TTBK1='TTBK1', 
         UCK1='UCK1' ,ULK4 = 'ULK4', UCK2 ='UCK2', VRK2 ='VRK2', WEE2='WEE2') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> KIRC.rnaseq1  



KIRC.clinical.surv %>%
  left_join(KIRC.rnaseq1,
            by = "bcr_patient_barcode") ->
  KIRC.surv_rnaseq


table(KIRC.surv_rnaseq$cohort, useNA = "always")


KIRC.surv_rnaseq <- KIRC.surv_rnaseq %>%
  filter(!is.na(cohort))


dim(KIRC.surv_rnaseq)
######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(KIRC.surv_rnaseq))
gen <- gen[-c(1:4)]

geneid <- gen
#geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
#           "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
#          "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
#          "PLK5","RPS6KL1", "SGK223")
geneid<- geneid[-c(31, 39, 97, 122)]

KIRC.surv_rnaseq<- KIRC.surv_rnaseq[,-c(35, 43, 101, 126)]
######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
for(i in 1: length(geneid)){
  
  KIRC.surv_rnaseq.cut <- surv_cutpoint(
    KIRC.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(KIRC.surv_rnaseq.cut)
  
  #distri <- plot(KIRC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  KIRC.surv_rnaseq.cat <- surv_categorize(KIRC.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = KIRC.surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  save(pval,fit,KIRC.surv_rnaseq.cat, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/survival/", geneid[i], "KIRC_Pvalue.RData", sep ="") )
  
  png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/survival/", geneid[i], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= KIRC.surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
  )  
  
  print(gg)
  dev.off()
}

write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/KIRC_Pvalue.csv", row.names = FALSE)
save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/KIRC_Pvalue.RData") 




###############
#concatenate survival P-value for each cancer
########

p <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/KIRC_Pvalue.csv")

p <- p[,c(1,4)]
p1 <- separate(data = p, col = variable, into = c("left", "right"), sep = "\\+")
p1 <- separate(data = p1, col = pval.txt, into = c("left1", "right1"), sep = "\\=")
which(p$pval.txt == "p < 0.0001")
p2 <- p1[,c(1,4)]

p2$right1[is.na(p2$right1)] <- "< 0.0001"


colnames(p2) <- c("Genes", "Survival_P_KIRC")

#p3 <- p2

Survival <- cbind(Survival, Survival_P_KIRC = p2$Survival_P_KIRC)

head(Survival)

write.csv(Survival, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/FinalTableSurvivalwithPvalue.csv", row.names = FALSE)



###########
#build a dataframe 
#######
tt <- data.frame()
for(bb in 1:length(geneid)){
  for(cc in 1:length(CanType)){
load(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[cc], "/survival/",  geneid[bb], CanType[cc],"_Pvalue.RData", sep=""))
gg <- get(paste(CanType[cc],".surv_rnaseq.cat", sep=""))

colnames(gg) <- c("times", "patient.vital_status", "expr", "cohart")
gg$gene <- paste(geneid[bb])
gg$cohart <- paste(CanType[cc])

tt <- rbind(tt, gg)
     
     }
}

write.csv(tt, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/data/Survival.csv")


###########################################################
# mapping kinases data with DE data
#get the DE 
########################################
#################

#CanID <- c("TCGA-LGG","TCGA-ACC","TCGA-SARC","TCGA-THYM","TCGA-OV","TCGA-UVM","TCGA-SKCM","TCGA-PAAD","TCGA-TGCT","TCGA-UCS")

Kinases <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/listofkinases.csv")
Kinases <- Kinases$sym

kk <- data.frame()
for(i in 2:length(CanType)){
  print(paste("starting", CanType[i]))
load(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[i], "/AllGeneswithFC.RData", sep ="")) #dataDEGs
#Kinases_filter <- rownames(dataDEGs)[Kinases,]


Kinases_filter <- dataDEGs[rownames(dataDEGs) %in% Kinases,]

Filter <- dataDEGs[grep("AMHR2|CDK9|CDKL4|CIT|DAPK3|DDR1|EEF2K|GRK7|GK2|GUCY2F|HASPIN|PRKACG|PSKH2|IRAK1|IRAK3|CSNK1E|MAP3K21|MAP3K19|MAP3K20|MOS|NME2P1|NRBP2|NMRK2|PDXK|PAK6|PDPK2P|PINK1|PGK2|PLK5|PRAG1|GRK1|SBK2|SBK3|SPEG|TAOK1|UCK2|UCKL1|WEE2", rownames(dataDEGs)),]

#part_kinases <- Kinases_filter[grep("^c", rownames(Kinases_filter)),]

#Kinases_filter1 <- Kinases_filter[-which(rownames(Kinases_filter) %in% rownames(part_kinases)), ]

part_kinases <- Filter[grep("^c", rownames(Filter)),]

split1 <- unlist(strsplit(as.character(rownames(part_kinases)), "[()]"))
n <- length(split1) 
split2 <- split1[seq(2, n, 2)] 

#remove everything after , 

split3 <- gsub(",.*","",split2)
dd <- gsub('[[:punct:]]',"",split3)
rownames(part_kinases) <- dd


Kinases_filter1 <- rbind(Kinases_filter, part_kinases)


sc <- dim(Kinases_filter)[1] + dim(part_kinases)[1] 
sc == dim(Kinases_filter1)[1]

Underst_adj.P <- select(Kinases_filter1, logFC, FDR)

print(setdiff(Kinases, rownames(Underst_adj.P)))

length(setdiff(Kinases, rownames(Underst_adj.P))) + dim(Underst_adj.P)[1] == length(Kinases)

colnames(Underst_adj.P) <- c(paste("logFC_",CanType[i], sep=""),  paste("FDR_",CanType[i], sep=""))
Underst_adj.P$Row.names <- rownames(Underst_adj.P)

#kk <- Underst_adj.P
kk <- merge(as.data.frame(kk), as.data.frame(Underst_adj.P), by='Row.names', all=TRUE)
head(kk)
dim(kk)
#kk <- cbind(kk, Underst_adj.P)

print(paste("the foldchange and P_values are extracted for::::::", CanType[i]))
#rownames(kk) == rownames(Underst_adj.P)
}

write.csv(kk, file= "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/Kinases_all_Pvalue_FC.csv", row.names = FALSE)


###########################################################
# Survival analysis of kinases data
#plus which do not have control sample ("LGG","ACC","SARC","THYM","OV","UVM","SKCM","PAAD","TGCT","UCS")
########################################
library(survival)
library(RTCGA.rnaseq)
library(RTCGA.clinical)
library(dplyr)
library(tidyr)
library(survminer)

#
CanType <-  c("BLCA",  "KIRP","KIRC", "BRCA", "COAD", "ESCA",  "HNSC" , "LIHC", "LUAD", "LUSC", "CESC", "GBM" , "PRAD", "CHOL", "KICH", "UCEC" ,  "READ" ,"PCPG", "STAD", "THCA", "LGG","ACC","SARC","THYM","OV","UVM","SKCM","PAAD","TGCT","UCS")   #"CHOL", "KICH", "UCEC" ,  "READ" ,"PCPG" #

Kinases <- as.character(Kinases)

for(i in 1:length(CanType)){
  
print(paste("Processing the Survival Analysis for", CanType[i]))
clinical <- readRDS(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[i], "/TCGA-", CanType[i],"_clinical.rds", sep=""))




#extract attributes needed for survival analysis
clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                        patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = clinical$xml_days_to_death)

load(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/",CanType[i], "_Normalize_filter.RData", sep=""))

#clinical.surv <- survivalTCGA(get(paste(CanType[i],".clin", sep="")))
clinical.surv <- survivalTCGA(clin)

norm <- get(paste(CanType[i],".Normalize.expression", sep =""))

norm_filter <- norm[grep("^c", rownames(norm)),]

norm_filter1 <- norm[-which(rownames(norm) %in% rownames(norm_filter)), ]

#split the string with 
split1 <- unlist(strsplit(as.character(rownames(norm_filter)), "[()]"))
n <- length(split1) 
split2 <- split1[seq(2, n, 2)] 
split3 <- gsub(",.*","",split2)
split3 <- gsub('[[:punct:]]',"",split3) 
#split3 <- gsub("[^[:alnum:][:blank:]+]", "", split2)


#split4 <- gsub("[ ]", "/", split3)

# split2 <- noquote(split2)
# split3 <- gsub(",.*","",split2)
#split3 <- gsub("\\,", "/", split2)
#remove everything after , 
# split3 <- gsub('[[:punct:]]',"",split3) 

rownames(norm_filter) <- split3
norm_filter2 <- rbind(norm_filter1, norm_filter)
dim(norm_filter2) == dim(norm)


rnaseq <- data.frame(t(norm_filter2))
rnaseq$bcr_patient_barcode <- rownames(rnaseq)
rownames(rnaseq) <- NULL

#"NME1-NME2" = "NME1-NME2" , "HASPIN" = "HASPIN" ,"MAP3K21" = "MAP3K21" ,"PDPK2P" = "PDPK2P" ,"PRAG1" = "PRAG1" ,"NME2P1" = "NME2P1" ,"MAP3K20" = "MAP3K20" ,

expressionsTCGA(
  rnaseq,
  extract.cols = Kinases) %>%
  rename(cohort = dataset,
         "ACVR1" = "ACVR1" , "ADCK2" = "ADCK2" , "ADCK5" = "ADCK5" , "ADK" = "ADK" , "ABL2" = "ABL2" , "AAK1" = "AAK1" , 
         "TNK2" = "TNK2" , "ACVR1C" = "ACVR1C" , "ALK" = "ALK" , "ACVR1B" = "ACVR1B" , "ACVRL1" = "ACVRL1" , "PRKAB1" = "PRKAB1" ,
         "PRKAA1" = "PRKAA1" , "ABL1" = "ABL1" , "AMHR2" = "AMHR2" , "ADCK1" = "ADCK1" , "AKT1" = "AKT1" , "PRKAG1" = "PRKAG1" ,
         "PRKAA2" = "PRKAA2" , "ARAF" = "ARAF" , "ALPK3" = "ALPK3" , "NPR1" = "NPR1" , "ADPGK" = "ADPGK" , "ANKK1" = "ANKK1" , 
         "NPR2" = "NPR2" , "ALPK1" = "ALPK1" , "ALPK2" = "ALPK2" , "GRK2" = "GRK2" , "GRK3" = "GRK3" , "AKT2" = "AKT2" , 
         "AKT3" = "AKT3" , "ATM" = "ATM" , "ATR" = "ATR" , "AURKC" = "AURKC" , "ACVR2A" = "ACVR2A" , "ACVR2B" = "ACVR2B" ,
         "BLK" = "BLK" , "BCKDK" = "BCKDK" , "AURKB" = "AURKB" , "AURKA" = "AURKA" , "BRAF" = "BRAF" , "BMPR1B" = "BMPR1B" , 
         "BMX" = "BMX" , "BTK" = "BTK" , "CAMKV" = "CAMKV" , "BMPR2" = "BMPR2" , "BMPR1A" = "BMPR1A" , "BRSK1" = "BRSK1" , 
         "BRSK2" = "BRSK2" , "BMP2K" = "BMP2K" , "BUB1B" = "BUB1B" , "BUB1" = "BUB1" , "CDC7" = "CDC7" , "CDK17" = "CDK17" , 
         "CDK18" = "CDK18" , "CDK9" = "CDK9" , "CHKA" = "CHKA" , "CLK2" = "CLK2" , "CLK4" = "CLK4" , "CMPK2" = "CMPK2" , 
         "CHEK1" = "CHEK1" , "CHKB" = "CHKB" , "CHEK2" = "CHEK2" , "CDK2" = "CDK2" , "CDK7" = "CDK7" , "CDK8" = "CDK8" , 
         "CDK10" = "CDK10" , "CDK11A" = "CDK11A" , "CDK19" = "CDK19" , "CDK1" = "CDK1" , "CDK5" = "CDK5" , "CDKL4" = "CDKL4" , 
         "CDK11B" = "CDK11B" , "CDK12" = "CDK12" , "CDK13" = "CDK13" , "CDK14" = "CDK14" , "CDK15" = "CDK15" , "CDK16" = "CDK16" , 
         "CDK20" = "CDK20" , "CDK3" = "CDK3" , "CDK4" = "CDK4" , "CDK6" = "CDK6" , "CDKL1" = "CDKL1" , "CDKL2" = "CDKL2" , 
         "CDKL3" = "CDKL3" , "CDKL5" = "CDKL5" , "CERK" = "CERK" , "CSNK2A2" = "CSNK2A2" , "COQ8A" = "COQ8A" , "COQ8B" = "COQ8B" ,
         "CLK3" = "CLK3" , "CLK1" = "CLK1" , "CSNK2A1" = "CSNK2A1" , "CSNK2A3" = "CSNK2A3" , "CSNK2B" = "CSNK2B" , "CIT" = "CIT" , "CSK" = "CSK" , "CSF1R" = "CSF1R" , "CASK" = "CASK" , "DCLK3" = "DCLK3" , "DAPK2" = "DAPK2" , "DCLK1" = "DCLK1" , "DDR2" = "DDR2" , "DAPK1" = "DAPK1" , "DAPK3" = "DAPK3" , "DCK" = "DCK" , "DCLK2" = "DCLK2" , "DGKQ" = "DGKQ" ,
         "DGKH" = "DGKH" , "DGUOK" = "DGUOK" , "DMPK" = "DMPK" , "DGKA" = "DGKA" , "DDR1" = "DDR1" , "DYRK1A" = "DYRK1A" , "EIF2AK1" = "EIF2AK1" , "EIF2AK4" = "EIF2AK4" , "DYRK1B" = "DYRK1B" , "DYRK3" = "DYRK3" , "EIF2AK2" = "EIF2AK2" , "EGFR" = "EGFR" , "ETNK1" = "ETNK1" , "ETNK2" = "ETNK2" , "EPHB2" = "EPHB2" , "ERN2" = "ERN2" , "DSTYK" = "DSTYK" ,
         "DYRK2" = "DYRK2" , "DYRK4" = "DYRK4" , "EIF2AK3" = "EIF2AK3" , "ERBB3" = "ERBB3" , "EPHA7" = "EPHA7" , "EPHB1" = "EPHB1" , "ERN1" = "ERN1" , "EEF2K" = "EEF2K" , "EFNA1" = "EFNA1" , "EPHA1" = "EPHA1" , "EPHA8" = "EPHA8" , "EPHA10" = "EPHA10" , "PFKFB2" = "PFKFB2" , "PTK2" = "PTK2" , "EPHA6" = "EPHA6" , "ERBB2" = "ERBB2" , "ERBB4" = "ERBB4" ,
         "EPHA2" = "EPHA2" , "EPHA3" = "EPHA3" , "EPHA4" = "EPHA4" , "EPHA5" = "EPHA5" , "EPHB3" = "EPHB3" , "EPHB4" = "EPHB4" , "EPHB6" = "EPHB6" , "FGFR2" = "FGFR2" , "FGFR3" = "FGFR3" , "FGFR4" = "FGFR4" , "PTK2B" = "PTK2B" , "FGFR1" = "FGFR1" , "FGGY" = "FGGY" , "PFKFB1" = "PFKFB1" , "PFKFB3" = "PFKFB3" , "FASTK" = "FASTK" , "FER" = "FER" ,
         "FN3K" = "FN3K" , "FES" = "FES" , "FUK" = "FUK" , "FYN" = "FYN" , "GAK" = "GAK" , "FLT3" = "FLT3" , "PIKFYVE" = "PIKFYVE" , "FRK" = "FRK" , "FGR" = "FGR" , "GALK1" = "GALK1" , "GALK2" = "GALK2" , "GK" = "GK" , "GLYCTK" = "GLYCTK" , "GRK7" = "GRK7" , "IDNK" = "IDNK" , "GUCY2C" = "GUCY2C" , 
         "GK5" = "GK5" , "GSK3A" = "GSK3A" , "GSK3B" = "GSK3B" , "GK3P" = "GK3P" , "GK2" = "GK2" , "GRK6" = "GRK6" , "GUCY2F" = "GUCY2F" ,  "GRK4" = "GRK4" , "GRK5" = "GRK5" , "HIPK4" = "HIPK4" , "GUCY2D" = "GUCY2D" , "MASTL" = "MASTL" , "HCK" = "HCK" , "HIPK1" = "HIPK1" , "HIPK2" = "HIPK2" , "HIPK3" = "HIPK3" , "HKDC1" = "HKDC1" ,
         "HUNK" = "HUNK" , "HK2" = "HK2" , "HK3" = "HK3" , "ICK" = "ICK" , "IGF1R" = "IGF1R" , "ILK" = "ILK" , "ITPKA" = "ITPKA" , "ITPKC" = "ITPKC" , "IP6K1" = "IP6K1" , "IP6K2" = "IP6K2" , "HK1" = "HK1" , "CHUK" = "CHUK" , "IKBKB" = "IKBKB" , "INSR" = "INSR" , "ITK" = "ITK" , "PRKACG" = "PRKACG" , "CSNK1A1L" = "CSNK1A1L" , "CSNK1A1" = "CSNK1A1" ,
         "CSNK1D" = "CSNK1D" , "CAMK2D" = "CAMK2D" , "CKB" = "CKB" , "IKBKE" = "IKBKE" , "IP6K3" = "IP6K3" , "JAK2" = "JAK2" , "JAK3" = "JAK3" , "AK1" = "AK1" , "AK2" = "AK2" , "KALRN" = "KALRN" , "JAK1" = "JAK1" , "AK3" = "AK3" , "AK5" = "AK5" , "AK6" = "AK6" , "AK8" = "AK8" , "AK9" = "AK9" , "PRKACA" = "PRKACA" , "CKM" = "CKM" , "CKMT1A" = "CKMT1A" , "TK1" = "TK1" , "AK4" = "AK4" , "CAMK1D" = "CAMK1D" , "CAMK1G" = "CAMK1G" , "CMPK1" = "CMPK1" , "GUK1" = "GUK1" , "PRKD3" = "PRKD3" , "PSKH2" = "PSKH2" , "INSRR" = "INSRR" , "IRAK1" = "IRAK1" , "IRAK4" = "IRAK4" , "PRKACB" = "PRKACB" , "CSNK1G1" = "CSNK1G1" , "CSNK1G3" = "CSNK1G3" , "CAMK2B" = "CAMK2B" , "CAMK4" = "CAMK4" , "GCK" = "GCK" , "IPMK" = "IPMK" , "IRAK3" = "IRAK3" , "ITPKB" = "ITPKB" , "IRAK2" = "IRAK2" , "ITPK1" = "ITPK1" , "CSNK1E" = "CSNK1E" , "CSNK1G2" = "CSNK1G2" , "CAMK1" = "CAMK1" , 
         "PNCK" = "PNCK" , "CAMK2A" = "CAMK2A" , "CAMK2G" = "CAMK2G" , "PRKG1" = "PRKG1" , "PRKG2" = "PRKG2" , "PSKH1" = "PSKH1" , "RPS6KA1" = "RPS6KA1" , "RPS6KB1" = "RPS6KB1" , "LMTK3" = "LMTK3" , "AATK" = "AATK" , "MVK" = "MVK" , "TK2" = "TK2" , "KIT" = "KIT" , "CAMKK2" = "CAMKK2" , "LCK" = "LCK" , "RPS6KA2" = "RPS6KA2" , "RPS6KA3" = "RPS6KA3" , "RPS6KA4" = "RPS6KA4" , "RPS6KA6" = "RPS6KA6" , "DTYMK" = "DTYMK" , "LMTK2" = "LMTK2" , "KHK" = "KHK" , "CAMKK1" = "CAMKK1" , "PHKA2" = "PHKA2" , "PRKCD" = "PRKCD" , "PRKCI" = "PRKCI" , "PRKCH" = "PRKCH" , "PRKCQ" = "PRKCQ" , "LATS2" = "LATS2" , "PRKCA" = "PRKCA" , "PRKCB" = "PRKCB" , "PRKCE" = "PRKCE" , "PRKCG" = "PRKCG" , "PRKCZ" = "PRKCZ" , "PKM" = "PKM" , "PKLR" = "PKLR" , "RPS6KA5" = "RPS6KA5" , "RPS6KC1" = "RPS6KC1" , "KSR2" = "KSR2" , "SYK" = "SYK" , "FN3KRP" = "FN3KRP" , "LATS1" = "LATS1" , "LIMK1" = "LIMK1" , "LIMK2" = "LIMK2" , "CKMT2" = "CKMT2" , "PHKA1" = "PHKA1" , "PRKD1" = "PRKD1" , "PRKD2" = "PRKD2" , "RPS6KB2" = "RPS6KB2" , "KSR1" = "KSR1" , "LRRK2" = "LRRK2" , "LYN" = "LYN" , "MAP3K15" = "MAP3K15" , "MAP3K9" = "MAP3K9" , "MATK" = "MATK" , "MET" = "MET" , "MAP3K10" = "MAP3K10" , "MAP3K11" = "MAP3K11" , "MAP3K13" = "MAP3K13" , "MAP3K14" = "MAP3K14" ,  "MAP3K4" = "MAP3K4" , "MAP3K5" = "MAP3K5" , "MAP3K7" = "MAP3K7" , "MAP4K1" = "MAP4K1" , "MAP4K5" = "MAP4K5" , "MAPKAPK3" = "MAPKAPK3" , "MAST1" = "MAST1" , "MAST2" = "MAST2" , "MAST4" = "MAST4" , "MAP3K12" = "MAP3K12" , "MAP3K2" = "MAP3K2" , "MAP4K3" = "MAP4K3" , "MAK" = "MAK" , "MAP3K19" = "MAP3K19" , "MAP3K1" = "MAP3K1" ,  "MAP3K3" = "MAP3K3" , "MAP3K6" = "MAP3K6" , "MAP4K2" = "MAP4K2" , "MARK1" = "MARK1" , "MARK3" = "MARK3" , "MAPKAPK2" = "MAPKAPK2" , "MAPKAPK5" = "MAPKAPK5" , "MAST3" = "MAST3" , "LRRK1" = "LRRK1" , "MAP4K4" = "MAP4K4" , "LTK" = "LTK" , "MAP3K8" = "MAP3K8" , "MARK2" = "MARK2" , "MARK4" = "MARK4" , "MAPK14" = "MAPK14" , "MLKL" = "MLKL" , "MOS" = "MOS" , "MOK" = "MOK" , "MAP2K2" = "MAP2K2" , "MAP2K6" = "MAP2K6" , "MAPK1" = "MAPK1" , "MAPK3" = "MAPK3" , "MAPK4" = "MAPK4" , "MAPK7" = "MAPK7" , "MAPK9" = "MAPK9" , "MAP2K1" = "MAP2K1" , "MAP2K3" = "MAP2K3" , "MAP2K4" = "MAP2K4" , "MAP2K5" = "MAP2K5" , "MAP2K7" = "MAP2K7" , "MINK1" = "MINK1" , "MAPK15" = "MAPK15" , "MKNK2" = "MKNK2" , "CDC42BPG" = "CDC42BPG" , "MELK" = "MELK" , "MERTK" = "MERTK" , "MAPK6" = "MAPK6" , "MAPK8" = "MAPK8" , "MAPK10" = "MAPK10" , "MAPK11" = "MAPK11" , "MAPK12" = "MAPK12" , "MAPK13" = "MAPK13" , "MKNK1" = "MKNK1" , "NME6" = "NME6" , "MYLK4" = "MYLK4" , "MYLK" = "MYLK" ,  "NADK" = "NADK" , "MYLK2" = "MYLK2" , "MYLK3" = "MYLK3" , "NADK2" = "NADK2" , "NME3" = "NME3" , "NEK3" = "NEK3" , "NEK5" = "NEK5" , "NEK7" = "NEK7" , "MUSK" = "MUSK" , "MYO3A" = "MYO3A" , "MYO3B" = "MYO3B" , "CDC42BPA" = "CDC42BPA" ,
         "CDC42BPB" = "CDC42BPB" , "MTOR" = "MTOR" , "NME4" = "NME4" , "NEK8" = "NEK8" , "NEK9" = "NEK9" , "IKBKG" = "IKBKG" , "NAGK" = "NAGK" , "NME5" = "NME5" , "NME1" = "NME1" ,  "NEK11" = "NEK11" , "NIM1K" = "NIM1K" , "NEK4" = "NEK4" , "NRBP2" = "NRBP2" , "NTRK2" = "NTRK2" , "NTRK3" = "NTRK3" , "NEK10" = "NEK10" , "NEK1" = "NEK1" , "NEK2" = "NEK2" , "NEK6" = "NEK6" , "NLK" = "NLK" , "NRBP1" = "NRBP1" , "NTRK1" = "NTRK1" , "NRK" = "NRK" , "NUAK2" = "NUAK2" , "PI4K2A" = "PI4K2A" , "PANK2" = "PANK2" , "PANK4" = "PANK4" , "NUAK1" = "NUAK1" , "PAK3" = "PAK3" , "PAK5" = "PAK5" , "PANK1" = "PANK1" , "NMRK1" = "NMRK1" , "NMRK2" = "NMRK2" , "OBSCN" = "OBSCN" , "SCYL3" = "SCYL3" , "OXSR1" = "OXSR1" , "PIK3C2A" = "PIK3C2A" , "PIK3C2B" = "PIK3C2B" , "PIK3C2G" = "PIK3C2G" , "PIK3R1" = "PIK3R1" , "PAK1" = "PAK1" , "PAK4" = "PAK4" , "PDIK1L" = "PDIK1L" , "PDK1" = "PDK1" , "PDK2" = "PDK2" , "PDK3" = "PDK3" , "PDK4" = "PDK4" , "PDXK" = "PDXK" , "PFKM" = "PFKM" , "PDGFRB" = "PDGFRB" , "PCK1" = "PCK1" , "PCK2" = "PCK2" , "PAK2" = "PAK2" , "PAK6" = "PAK6" , "PAN3" = "PAN3" , "PANK3" = "PANK3" ,  
         "PEAK1" = "PEAK1" , "PHKG1" = "PHKG1" , "PHKG2" = "PHKG2" , "PINK1" = "PINK1" , "PKDCC" = "PKDCC" , "PDPK1" = "PDPK1" , "PFKL" = "PFKL" , "PIM3" = "PIM3" , "PASK" = "PASK" , "PDGFRA" = "PDGFRA" , "PGK1" = "PGK1" , "PGK2" = "PGK2", 
         "PRPS1" = "PRPS1" , "PIM2" = "PIM2" , "PIK3C3" = "PIK3C3" , "PIK3CB" = "PIK3CB" , "PIK3CD" = "PIK3CD" , "PKN1" = "PKN1" , "PLK3" = "PLK3" , "PI4KB" = "PI4KB" , "PIP5K1C" = "PIP5K1C" , "PLK2" = "PLK2" , "PKMYT1" = "PKMYT1" , 
         "PIK3R4" = "PIK3R4" , "PIP5K1A" = "PIP5K1A" , "PIM1" = "PIM1" , "PKN2" = "PKN2" , "PLK1" = "PLK1" , "PLK5" = "PLK5" , "PIP4K2B" = "PIP4K2B" , "PIP4K2C" = "PIP4K2C" , "PIK3CA" = "PIK3CA" , "PFKP" = "PFKP" , "PIP4K2A" = "PIP4K2A" , "PI4KA" = "PI4KA" , "PIP5K1B" = "PIP5K1B" , "PIK3CG" = "PIK3CG" , "PKN3" = "PKN3" , "PLK4" = "PLK4" , "PRKX" = "PRKX" , 
          "PRKRA" = "PRKRA" , "PRPS2" = "PRPS2" , "PTK6" = "PTK6" , "PRKDC" = "PRKDC" , "PRKY" = "PRKY" , "PRPF4B" = "PRPF4B" , "PRPS1L1" = "PRPS1L1" , "PTK7" = "PTK7" , "PSTK" = "PSTK" , "TP53RK" = "TP53RK" , "PXK" = "PXK" , "RAF1" = "RAF1" , "RIOK1" = "RIOK1" , "RIOK2" = "RIOK2" , "RBKS" = "RBKS" , "RIOK3" = "RIOK3" , "RIPK2" = "RIPK2" , "RIPK4" = "RIPK4" , "RET" = "RET" , "RIPK1" = "RIPK1" , "RIPK3" = "RIPK3" , "GRK1" = "GRK1" , "ROR1" = "ROR1" , "RNASEL" = "RNASEL" , "MST1R" = "MST1R" , "ROS1" = "ROS1" , "ROCK1" = "ROCK1" , "ROCK2" = "ROCK2" , "ROR2" = "ROR2" , "RPS6KL1" = "RPS6KL1" , 
         "SCYL1" = "SCYL1" , "SCYL2" = "SCYL2" , "SIK1" = "SIK1" , "RYK" = "RYK" , "SHPK" = "SHPK" , "SBK1" = "SBK1" , "SBK2" = "SBK2" , "SBK3" = "SBK3" , "SNRK" = "SNRK" , "POMK" = "POMK" , "SGK494" = "SGK494" , "SGK3" = "SGK3" , "SPHK2" = "SPHK2" , "SLK" = "SLK" , "SRC" = "SRC" , "SPHK1" = "SPHK1" , "SRPK3" = "SRPK3" , "STK32C" = "STK32C" , "STK26" = "STK26" , "STK35" = "STK35" , "STKLD1" = "STKLD1" , "SRMS" = "SRMS" , "SRPK1" = "SRPK1" , "SRPK2" = "SRPK2" , "SGK1" = "SGK1" , "SGK2" = "SGK2" , "SPEG" = "SPEG" , "STK38L" = "STK38L" , "SIK2" = "SIK2" , "SIK3" = "SIK3" , "SMG1" = "SMG1" , "STK11" = "STK11" , "STK31" = "STK31" , "STRADA" = "STRADA" , "STK17A" = "STK17A" , "STK38" = "STK38" , "STK39" = "STK39" , "STRADB" = "STRADB" , "STYK1" = "STYK1" , "STK24" = "STK24" , "STK25" = "STK25" , "STK33" = "STK33" , "STK40" = "STK40" , "STK4" = "STK4" , "TAOK2" = "TAOK2" , "TAOK3" = "TAOK3" , "STK17B" = "STK17B" , "STK32A" = "STK32A" , "STK32B" = "STK32B" , "STK10" = "STK10" , "STK16" = "STK16" , "STK36" = "STK36" , "STK3" = "STK3" , "TAB1" = "TAB1" , "TEX14" = "TEX14" , "TAOK1" = "TAOK1" , "TEC" = "TEC" , "TTN" = "TTN" , "TLK1" = "TLK1" , "TAF1" = "TAF1" , "TBCK" = "TBCK" , "TGFBR1" = "TGFBR1" , "TESK1" = "TESK1" , "TESK2" = "TESK2" , "TIE1" = "TIE1" , "TEK" = "TEK" , "TBK1" = "TBK1" , "TNK1" = "TNK1" , "TSSK3" = "TSSK3" , "TSSK4" = "TSSK4" , "TYRO3" = "TYRO3" , "TGFBR2" = "TGFBR2" , "TRIB2" = "TRIB2" , "TNIK" = "TNIK" , "TRIB1" = "TRIB1" , "TRIO" = "TRIO" , "TSSK1B" = "TSSK1B" , "TSSK6" = "TSSK6" , "TTBK1" = "TTBK1" , "TTBK2" = "TTBK2" , "TLK2" = "TLK2" , "TPK1" = "TPK1" , "PBK" = "PBK" , "TRIB3" = "TRIB3" , "TRPM6" = "TRPM6" , "TNNI3K" = "TNNI3K" , "ULK2" = "ULK2" , "ULK3" = "ULK3" , "ULK4" = "ULK4" , "KDR" = "KDR" , "UCK2" = "UCK2" , "UCK1" = "UCK1" , "TYK2" = "TYK2" , "AXL" = "AXL" , "UCKL1" = "UCKL1" , "TXK" = "TXK" , "UHMK1" = "UHMK1" , "TSSK2" = "TSSK2" , "TTK" = "TTK" , "VRK1" = "VRK1" , "VRK3" = "VRK3" , "WEE2" = "WEE2" , "WNK4" = "WNK4" , "VRK2" = "VRK2" , "WEE1" = "WEE1" , "WNK3" = "WNK3" ,
         "ULK1" = "ULK1" , "FLT1" = "FLT1" , "FLT4" = "FLT4" , "WNK1" = "WNK1" , 
         "WNK2" = "WNK2" , "XYLB" = "XYLB" , "YES1" = "YES1" , "ZAP70" = "ZAP70") %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> rnaseq1 




#`CDK9`, `CIT`, `DAPK3`, `DDR1`, `EEF2K`, `HASPIN`, `IRAK1`, `IRAK3`, `CSNK1E`, `MAP3K21`, `MAP3K20`, `NME2P1`, `NME1-NME2`, `NRBP2`, `PDXK`, `PAK6`, `PDPK2P`, `PINK1`, `PRAG1`, `SPEG`, `TAOK1`, `UCK2`, `UCKL1`

clinical.surv %>%
  left_join(rnaseq1,
            by = "bcr_patient_barcode") ->
  surv_rnaseq


table(surv_rnaseq$cohort, useNA = "always")


surv_rnaseq <- surv_rnaseq %>%
  filter(!is.na(cohort))


dim(surv_rnaseq)


######################
#GENES for survival analysis
#######################

library(survminer)

gen <- as.character(colnames(surv_rnaseq))
gen <- gen[-c(1:4)]

gen <- gen[-grep("^GK2$|^PSKH2$|^MOS$", gen)]

setwd(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[i], "/survival/",  sep=""  ))
#dir.create(file.path(paste0(getwd(), "/Kinase/")))

######################
#Plot SURVIVAL CURVES FOR GENES
#######################
p <- data.frame()
#BLCA- #GK2 -180, PSKH2-237, 346-MOS
#KIRP - GK2 -180, PSKH2-237, 346-MOS

for(j in 1:length(gen)){
  
  print(paste("Survival Analysis for", gen[j]))
  
  surv_rnaseq.cut <- surv_cutpoint(
   surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(gen[j]), "cohort")
  )
  summary(surv_rnaseq.cut)
  
  #distri <- plot(LIHC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
 surv_rnaseq.cat <- surv_categorize(surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , gen[j], "+ cohort")),
                 data = surv_rnaseq.cat)
  
  pval <- surv_pvalue(fit)
  
  p <- rbind(p, pval) 
  
  
  save(pval,fit,surv_rnaseq.cat, rnaseq1,clinical.surv, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[i], "/survival/Kinase/", gen[j], "_", CanType[i], "_Pvalue.RData", sep ="") )
  
  png(file = paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[i], "/survival/Kinase/", gen[j], ".survival.jpg", sep=""))
  gg <-ggsurvplot(
    fit, 
    data= surv_rnaseq.cat,
    break.time.by = 500, 
    ggtheme = theme_bw(), 
    risk.table = TRUE,
    xlim = c(0,3000),
    xlab = "Time in days",
    pval = TRUE,
    font.x = 16, 
    font.y=16,
    risk.table.height = 0.25,
    fontsize=3, 
    risk.table.col= "strata",
    size=1, palette = c("#2E9FDF", "red"),
    risk.table.fontsize = 5, 
    risk.table.y.text = FALSE
    
  )  
  print(gg)
  dev.off()
  
  print(paste("done for", gen[j]))
}

print(paste("done for Cancer", CanType[i]))
write.csv(p, file=paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases/", CanType[i], "_Pvalue.csv", sep = ""))
save(p, file=paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases/", CanType[i], "_Pvalue.RData", sep = "")) 
}
#write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/BLCA_Pvalue.csv")
#save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/BLCA_Pvalue.RData)

###########################################################
# get all the Pvalues for kinases together in one csv file
########################################

file <- "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases"
Survival <- data_frame()
for(i in 2:length(CanType)){
p <- read.csv(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases/", CanType[i], "_Pvalue.csv", sep = ""))

p <- p %>% select(variable, pval.txt)
p1 <- separate(data = p, col = variable, into = c("left", "right"), sep = "\\+")
p1 <- separate(data = p1, col = pval.txt, into = c("left1", "right1"), sep = "\\=")
which(p$pval.txt == "p < 0.0001")
p2 <- p1[,c(1,4)]

p2$right1[is.na(p2$right1)] <- "< 0.0001"


colnames(p2) <- c("Genes", paste0("Survival_P_", CanType[i]))

#p3 <- p2

#Survival <- p2

#p2$Row.names <- p2$Genes
#Survival <- p2
#kk <- Underst_adj.P
Survival <- merge(as.data.frame(Survival), as.data.frame(p2), by='Genes', all=TRUE)

#Survival <- cbind(Survival, p2[,2,drop = FALSE])

print(head(Survival))

}

write.csv(Survival, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases/FinalTableSurvivalwithPvalue.csv", row.names = FALSE)




###########################################################
# mapping Pseudo- kinases data with DE data
#get the DE [56 psedokinases 54 are present the list of kinases]
########################################
#################

pseudoKinases <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/pseudoKinase_by_stephan.csv")
pseudoKinases <- pseudoKinases$gene

pseudoKinases <- gsub("_.*","",pseudoKinases)
pseudoKinases <- gsub("-.*","",pseudoKinases) #56

kk <- data.frame()
for(i in 2:length(CanType)){
  print(paste("starting", CanType[i]))
  load(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[i], "/AllGeneswithFC.RData", sep ="")) #dataDEGs
  #Kinases_filter <- rownames(dataDEGs)[Kinases,]
  

  pseudoKinases_filter <- dataDEGs[rownames(dataDEGs) %in% pseudoKinases,]
  
  
}






###########################################################
# Mutation Data anaysis   
#
########################################
#################

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")

library(maftools)
library(TCGAmutations)


TCGAmutations::tcga_load(study = "KIRC")


#Typing tcga_kirc prints summary of the object
tcga_kirc_mc3

tcga_mc3 <- read.maf(maf = tcga_kirc_mc3@data, clinicalData = tcga_kirc_mc3@clinical.data)
#Shows sample summary
getSampleSummary(tcga_kirc_mc3)

#Shows gene summary
getGeneSummary(tcga_kirc_mc3)

#Clinical data; printing only first ten columns for display convenience
getClinicalData(tcga_kirc_mc3)[1:10, 1:10]


#Plotting MAF summary
png(file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/Mutation_New/mafSummary_kirc.png", width = 1000, height = 1000, units = "px")

plotmafSummary(maf = tcga_mc3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = T)
dev.off()

tcga_load(study = "BRCA")
#We will draw oncoplots for top ten mutated genes.
#oncoplot(maf = tcga_kirc_mc3, top = 10, fontSize = 12)


#draw any number of genes of ur choice understudied 
png(file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/Mutation/understudied_kirc.jpg", width = 1000, height = 1000, units = "px")
oncostrip(maf = tcga_kirc_mc3, genes = c('PRKAG1','PRKAB1','ADCK2','ADCK5','ALPK3','ADCK4','ALPK2','BCKDK','CAMKV','CDKL3','CDKL4','CDKL1','CDKL2','CDK15','CDK17' , 'CDK11B','CDK14','CDK18','CLK4','CDK10','CLK3','CSNK2A2' , 'DGKH','DYRK3','DYRK1B','DYRK4','EEF2K','DYRK2','ERN2','FASTK' , 'GK2','HIPK4','ITPKA','CSNK1G3','CSNK1G2','PRKACG','CSNK1G1','ITPK1','PSKH2','CAMK1D','TK2','PRKACB','PNCK','CAMK1G','CAMKK1' , 'PSKH1','PHKA1','RPS6KC1' , 'PRKCQ','LMTK3','MARK4','MAST3','LRRK1','MAP3K14','LTK','MAP3K10' , 'MAST4','MAST2','MAPK15','MAPK4' , 'MKNK2','CDC42BPG','CDC42BPB','NEK11','NEK5','NEK10','NIM1K','NEK4','NEK6','NRK','NEK7','SCYL1','NRBP2','PIK3C2G' , 'PAK7' , 'PANK3','PIK3C2B' , 'PAK3','PAK6','SCYL3','PDIK1L','PIP4K2C','PIP5K1A','PLK5','PHKG1','PHKG2','PI4KA','PIP5K1B','PKMYT1','PKN3' , 'PRPF4B','TP53RK','PRKRA','PRKY','PXK','RIOK1','RIOK3','RPS6KL1' , 'SBK2','SCYL2','SBK3','SGK494','POMK','SGK223','SRPK3' , 'STK32C','STK33','STK32B','STK40','STK17A','STK36','STK32A','STK38L','STKLD1','STK3','STK31','TESK2','TESK1','TBCK','TLK1' , 'TPK1','TLK2','TSSK6','TSSK1B','TTBK2','TSSK4','TTBK1','TSSK3','UCK1','ULK4','UCK2','WNK2','VRK2','WEE2' ),
          fontSize = 8)
dev.off()



#subset the data  




gg <- subsetMaf(maf = tcga_kirc_mc3, genes = c('ACVR1' , 'ADCK2' , 'ADCK5' , 'ADK' , 'ABL2' , 'AAK1' , 'TNK2' , 'ACVR1C' , 'ALK' , 'ACVR1B' , 'ACVRL1' , 'PRKAB1' , 'PRKAA1' , 'ABL1' , 'AMHR2' , 'ADCK1' , 'AKT1' , 
                                               'PRKAG1' , 'PRKAA2' , 'ARAF' , 'ALPK3' , 'NPR1' , 'ADPGK' , 'ANKK1' , 'NPR2' , 'ALPK1' , 'ALPK2' , 'GRK2' , 'GRK3' , 'AKT2' , 'AKT3' , 'ATM' , 'ATR' , 'AURKC' , 'ACVR2A' ,
                                               'ACVR2B' , 'BLK' , 'BCKDK' , 'AURKB' , 'AURKA' , 'BRAF' , 'BMPR1B' , 'BMX' , 'BTK' , 'CAMKV' , 'BMPR2' , 'BMPR1A' , 'BRSK1' , 'BRSK2' , 'BMP2K' , 'BUB1B' , 'BUB1' , 'CDC7' ,
                                               'CDK17' , 'CDK18' , 'CDK9' , 'CHKA' , 'CLK2' , 'CLK4' , 'CMPK2' , 'CHEK1' , 'CHKB' , 'CHEK2' , 'CDK2' , 'CDK7' , 'CDK8' , 'CDK10' , 'CDK11A' , 'CDK19' , 'CDK1' , 'CDK5' , 
                                               'CDKL4' , 'CDK11B' , 'CDK12' , 'CDK13' , 'CDK14' , 'CDK15' , 'CDK16' , 'CDK20' , 'CDK3' , 'CDK4' , 'CDK6' , 'CDKL1' , 'CDKL2' , 'CDKL3' , 'CDKL5' , 'CERK' , 'CSNK2A2' , 
                                               'COQ8A' , 'COQ8B' , 'CLK3' , 'CLK1' , 'CSNK2A1' , 'CSNK2A3' , 'CSNK2B' , 'CIT' , 'CSK' , 'CSF1R' , 'CASK' , 'DCLK3' , 'DAPK2' , 'DCLK1' , 'DDR2' , 'DAPK1' , 'DAPK3' , 
                                               'DCK' , 'DCLK2' , 'DGKQ' , 'DGKH' , 'DGUOK' , 'DMPK' , 'DGKA' , 'DDR1' , 'DYRK1A' , 'EIF2AK1' , 'EIF2AK4' , 'DYRK1B' , 'DYRK3' , 'EIF2AK2' , 'EGFR' , 'ETNK1' , 'ETNK2' , 
                                               'EPHB2' , 'ERN2' , 'DSTYK' , 'DYRK2' , 'DYRK4' , 'EIF2AK3' , 'ERBB3' , 'EPHA7' , 'EPHB1' , 'ERN1' , 'EEF2K' , 'EFNA1' , 'EPHA1' , 'EPHA8' , 'EPHA10' , 'PFKFB2' , 'PTK2' , 
                                               'EPHA6' , 'ERBB2' , 'ERBB4' , 'EPHA2' , 'EPHA3' , 'EPHA4' , 'EPHA5' , 'EPHB3' , 'EPHB4' , 'EPHB6' , 'FGFR2' , 'FGFR3' , 'FGFR4' , 'PTK2B' , 'FGFR1' , 'FGGY' , 'PFKFB1' ,
                                               'PFKFB3' , 'FASTK' , 'FER' , 'FN3K' , 'FES' , 'FUK' , 'FYN' , 'GAK' , 'FLT3' , 'PIKFYVE' , 'FRK' , 'FGR' , 'GALK1' , 'GALK2' , 'GK' , 'GLYCTK' , 'GRK7' , 'IDNK' , 'GUCY2C' , 'GK5' , 'GSK3A' , 'GSK3B' , 'GK3P' , 'GK2' , 'GRK6' , 'GUCY2F' , 'HASPIN' , 'GRK4' , 'GRK5' , 'HIPK4' , 'GUCY2D' , 'MASTL' , 'HCK' , 'HIPK1' , 'HIPK2' , 'HIPK3' , 'HKDC1' , 'HUNK' , 'HK2' , 'HK3' , 'ICK' , 'IGF1R' , 'ILK' , 'ITPKA' , 'ITPKC' , 'IP6K1' , 'IP6K2' , 'HK1' , 'CHUK' , 'IKBKB' , 'INSR' , 'ITK' , 'PRKACG' , 'CSNK1A1L' , 'CSNK1A1' , 'CSNK1D' , 'CAMK2D' , 'CKB' , 'IKBKE' , 'IP6K3' , 'JAK2' , 'JAK3' , 'AK1' , 'AK2' , 'KALRN' , 'JAK1' , 'AK3' , 'AK5' , 'AK6' , 'AK8' , 'AK9' , 'PRKACA' , 'CKM' , 'CKMT1A' , 'TK1' , 'AK4' , 'CAMK1D' , 'CAMK1G' , 'CMPK1' , 'GUK1' , 'PRKD3' , 'PSKH2' , 'INSRR' , 'IRAK1' , 'IRAK4' , 'PRKACB' , 'CSNK1G1' , 'CSNK1G3' , 'CAMK2B' , 'CAMK4' , 'GCK' , 'IPMK' , 'IRAK3' , 'ITPKB' , 'IRAK2' , 'ITPK1' , 'CSNK1E' , 'CSNK1G2' , 'CAMK1' , 'PNCK' , 'CAMK2A' , 'CAMK2G' , 'PRKG1' , 'PRKG2' , 'PSKH1' , 'RPS6KA1' , 'RPS6KB1' , 'LMTK3' , 'AATK' , 'MVK' , 'TK2' , 'KIT' , 'CAMKK2' , 'LCK' , 'RPS6KA2' , 'RPS6KA3' , 'RPS6KA4' , 'RPS6KA6' , 'DTYMK' , 'LMTK2' , 'KHK' , 'CAMKK1' , 'PHKA2' , 'PRKCD' , 'PRKCI' , 'PRKCH' , 'PRKCQ' , 'LATS2' , 'PRKCA' , 'PRKCB' , 'PRKCE' , 'PRKCG' , 'PRKCZ' , 'PKM' , 'PKLR' , 'RPS6KA5' , 'RPS6KC1' , 'KSR2' , 'SYK' , 'FN3KRP' , 'LATS1' , 'LIMK1' , 'LIMK2' , 'CKMT2' , 'PHKA1' , 'PRKD1' , 'PRKD2' , 'RPS6KB2' , 'KSR1' , 'LRRK2' , 'LYN' , 'MAP3K15' , 'MAP3K9' , 'MATK' , 'MET' , 'MAP3K10' , 'MAP3K11' , 'MAP3K13' , 'MAP3K14' , 'MAP3K21' , 'MAP3K4' , 'MAP3K5' , 'MAP3K7' , 'MAP4K1' , 'MAP4K5' , 'MAPKAPK3' , 'MAST1' , 'MAST2' , 'MAST4' , 'MAP3K12' , 'MAP3K2' , 'MAP4K3' , 'MAK' , 'MAP3K19' , 'MAP3K1' , 'MAP3K20' , 'MAP3K3' , 'MAP3K6' , 'MAP4K2' , 'MARK1' , 'MARK3' , 'MAPKAPK2' , 'MAPKAPK5' , 'MAST3' , 'LRRK1' , 'MAP4K4' , 'LTK' , 'MAP3K8' , 'MARK2' , 'MARK4' , 'MAPK14' , 'MLKL' , 'MOS' , 'MOK' , 'MAP2K2' , 'MAP2K6' , 'MAPK1' , 'MAPK3' , 'MAPK4' , 'MAPK7' , 'MAPK9' , 'MAP2K1' , 'MAP2K3' , 'MAP2K4' , 'MAP2K5' , 'MAP2K7' , 'MINK1' , 'MAPK15' , 'MKNK2' , 'CDC42BPG' , 'MELK' , 'MERTK' , 'MAPK6' , 'MAPK8' , 'MAPK10' , 'MAPK11' , 'MAPK12' , 'MAPK13' , 'MKNK1' , 'NME6' , 'MYLK4' , 'MYLK' , 'NME2P1' , 'NADK' , 'MYLK2' , 'MYLK3' , 'NADK2' , 'NME3' , 'NEK3' , 'NEK5' , 'NEK7' , 'MUSK' , 'MYO3A' , 'MYO3B' , 'CDC42BPA' , 'CDC42BPB' , 'MTOR' , 'NME4' , 'NEK8' , 'NEK9' , 'IKBKG' , 'NAGK' , 'NME5' , 'NME1' , 'NME1-NME2' , 'NEK11' , 'NIM1K' , 'NEK4' , 'NRBP2' , 'NTRK2' , 'NTRK3' , 'NEK10' , 'NEK1' , 'NEK2' , 'NEK6' , 'NLK' , 'NRBP1' , 'NTRK1' , 'NRK' , 'NUAK2' , 'PI4K2A' , 'PANK2' , 'PANK4' , 'NUAK1' , 'PAK3' , 'PAK5' , 'PANK1' , 'NMRK1' , 'NMRK2' , 'OBSCN' , 'SCYL3' , 'OXSR1' , 'PIK3C2A' , 'PIK3C2B' , 'PIK3C2G' , 'PIK3R1' , 'PAK1' , 'PAK4' , 'PDIK1L' , 'PDK1' , 'PDK2' , 'PDK3' , 'PDK4' , 'PDXK' , 'PFKM' , 'PDGFRB' , 'PCK1' , 'PCK2' , 'PAK2' , 'PAK6' , 'PAN3' , 'PANK3' , 'PDPK2P' , 'PEAK1' , 'PHKG1' , 'PHKG2' , 'PINK1' , 'PKDCC' , 'PDPK1' , 'PFKL' , 'PIM3' , 'PASK' , 'PDGFRA' , 'PGK1' , 'PGK2' , 'PRPS1' , 'PIM2' , 'PIK3C3' , 'PIK3CB' , 'PIK3CD' , 'PKN1' , 'PLK3' , 'PI4KB' , 'PIP5K1C' , 'PLK2' , 'PKMYT1' , 'PIK3R4' , 'PIP5K1A' , 'PIM1' , 'PKN2' , 'PLK1' , 'PLK5' , 'PIP4K2B' , 'PIP4K2C' , 'PIK3CA' , 
                                               'PFKP' , 'PIP4K2A' , 'PI4KA' , 'PIP5K1B' , 'PIK3CG' , 'PKN3' , 'PLK4' , 'PRKX' , 'PRAG1' , 'PRKRA' , 'PRPS2' , 'PTK6' , 'PRKDC' , 'PRKY' , 'PRPF4B' , 'PRPS1L1' , 
                                               'PTK7' , 'PSTK' , 'TP53RK' , 'PXK' , 'RAF1' , 'RIOK1' , 'RIOK2' , 'RBKS' , 'RIOK3' , 'RIPK2' , 'RIPK4' , 'RET' , 'RIPK1' , 'RIPK3' , 'GRK1' , 'ROR1' , 'RNASEL' , 
                                               'MST1R' , 'ROS1' , 'ROCK1' , 'ROCK2' , 'ROR2' , 'RPS6KL1' , 'SCYL1' , 'SCYL2' , 'SIK1' , 'RYK' , 'SHPK' , 'SBK1' , 'SBK2' , 'SBK3' , 'SNRK' , 'POMK' , 'SGK494' , 
                                               'SGK3' , 'SPHK2' , 'SLK' , 'SRC' , 'SPHK1' , 'SRPK3' , 'STK32C' , 'STK26' , 'STK35' , 'STKLD1' , 'SRMS' , 'SRPK1' , 'SRPK2' , 'SGK1' , 'SGK2' , 'SPEG' , 'STK38L' ,
                                               'SIK2' , 'SIK3' , 'SMG1' , 'STK11' , 'STK31' , 'STRADA' , 'STK17A' , 'STK38' , 'STK39' , 'STRADB' , 'STYK1' , 'STK24' , 'STK25' , 'STK33' , 'STK40' , 'STK4' , 'TAOK2' , 
                                               'TAOK3' , 'STK17B' , 'STK32A' , 'STK32B' , 'STK10' , 'STK16' , 'STK36' , 'STK3' , 'TAB1' , 'TEX14' , 'TAOK1' , 'TEC' , 'TTN' , 'TLK1' , 'TAF1' , 'TBCK' , 'TGFBR1' , 
                                               'TESK1' , 'TESK2' , 'TIE1' , 'TEK' , 'TBK1' , 'TNK1' , 'TSSK3' , 'TSSK4' , 'TYRO3' , 'TGFBR2' , 'TRIB2' , 'TNIK' , 'TRIB1' , 'TRIO' , 'TSSK1B' , 'TSSK6' , 'TTBK1' , 
                                               'TTBK2' , 'TLK2' , 'TPK1' , 'PBK' , 'TRIB3' , 'TRPM6' , 'TNNI3K' , 'ULK2' , 'ULK3' , 'ULK4' , 'KDR' , 'UCK2' , 'UCK1' , 'TYK2' , 'AXL' , 'UCKL1' , 'TXK' , 'UHMK1' ,
                                               'TSSK2' , 'TTK' , 'VRK1' , 'VRK3' , 'WEE2' , 'WNK4' , 'VRK2' , 'WEE1' , 'WNK3' , 'ULK1' , 'FLT1' , 'FLT4' , 'WNK1' , 'WNK2' , 'XYLB' , 'YES1' , 'ZAP70' ))



sort(table(gg$Hugo_Symbol))

gg1 <- data.frame(table(gg$Hugo_Symbol))

gg2 <- gg1[which(gg1$Freq > 5),]

#draw any number of genes of ur choice kinases
png(file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/Mutation/kinases_kirc.jpg", width = 1000, height = 1000, units = "px")
oncostrip(maf = tcga_kirc_mc3, genes = as.character(gg2$Var1), fontSize = 8)
dev.off()





#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
#lollipopPlot(maf = tcga_kirc_mc3, gene = 'DNMT3A', AACol = 'GMAF', showMutationRate = TRUE)  #Protein_change for other cancers



#rainfallPlot(maf = tcga_kirc_mc3, detectChangePoints = TRUE, fontSize = 12, pointSize = 0.6)


#Detecting cancer driver genes based on positional clustering.
sig = oncodrive(maf = tcga_mc3, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')

head(sig)

sig1 <- select(sig, Hugo_Symbol, clusterScores, protLen, zscore,  pval, fdr, fract_muts_in_clusters)



sig1 <- sig1[which(sig1$Hugo_Symbol  %in% gg1$Var1), ]

sig1$cancer <- "KIRC"

jj <- select(sig1, Hugo_Symbol, pval, fdr)
colnames(jj) <- c("Hugo_Symbol", "pval_KIRC", "fdr_KIRC")

mergeJJ <- jj


png(file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/Mutation_New/plotOncodrive_kirc.jpg", width = 1000, height = 1000, units = "px")
plotOncodrive(res = sig, fdrCutOff = 0.1, useFraction = TRUE)
dev.off()

save(tcga_mc3, sig, gg1, jj, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/Mutation_New/AllDatawithsig.RData")

#######completed for KIRC 
#######loop to start with other cancers



###################
# Mutation analysis for other cancers 
#################

CanType <- c("BLCA",  "KIRP","BRCA", "COAD", "ESCA",  "HNSC" , "LIHC", "LUAD", "LUSC", "STAD", "THCA", "CESC", "GBM" ,  "KICH", "UCEC" ,  "READ" ,"PCPG", "PRAD", "CHOL")   #"CHOL", "KICH", "UCEC" ,  "READ" ,"PCPG" #
#KICH has only one signature gene "TP53" and only one gene above 5 mutation "TTN" so could not plot  plotOncodrive_KICH.png, kinases_KICH.png
#CHOL has only one gene "TTN above 5 mutations so  plotOncodrive_CHOL.png is with all the genes

for(qq in 1:length(CanType) ) {
  print(paste("Started with", CanType[qq], sep ="-------"))
TCGAmutations::tcga_load(study = paste(CanType[qq]))

tcga_mc3 <- get(paste("tcga_",tolower(CanType[qq]),"_mc3", sep =""))

#Typing tcga_kirc prints summary of the object
#tcga_kirc_mc3

#Shows sample summary
getSampleSummary(tcga_mc3)

#Shows gene summary
getGeneSummary(tcga_mc3)

#Clinical data; printing only first ten columns for display convenience
getClinicalData(tcga_mc3)[1:10, 1:10]

#
dir.create(file.path(filePath, paste("TCGA-", CanType[qq], sep = ""), "Mutation"))
#setwd(file.path(filePath, paste("TCGA-", CanType[3], sep =""), "Mutation/" ))

#Plotting MAF summary
png(file =paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[qq], "/Mutation/MAFSummary.png", sep =""),  width = 1000, height = 1000, units = "px")
plotmafSummary(maf = tcga_mc3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()


#We will draw oncoplots for top ten mutated genes.
#oncoplot(maf = tcga_mc3, top = 10, fontSize = 12)


#draw any number of genes of ur choice understudied 
png(file =paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[qq], "/Mutation/understudied_", CanType[qq], ".png", sep =""), width = 1000, height = 1000, units = "px")
oncostrip(maf = tcga_mc3, genes = c('PRKAG1','PRKAB1','ADCK2','ADCK5','ALPK3','ADCK4','ALPK2','BCKDK','CAMKV','CDKL3','CDKL4','CDKL1','CDKL2','CDK15','CDK17' , 'CDK11B','CDK14','CDK18','CLK4','CDK10','CLK3','CSNK2A2' , 'DGKH','DYRK3','DYRK1B','DYRK4','EEF2K','DYRK2','ERN2','FASTK' , 'GK2','HIPK4','ITPKA','CSNK1G3','CSNK1G2','PRKACG','CSNK1G1','ITPK1','PSKH2','CAMK1D','TK2','PRKACB','PNCK','CAMK1G','CAMKK1' , 'PSKH1','PHKA1','RPS6KC1' , 'PRKCQ','LMTK3','MARK4','MAST3','LRRK1','MAP3K14','LTK','MAP3K10' , 'MAST4','MAST2','MAPK15','MAPK4' , 'MKNK2','CDC42BPG','CDC42BPB','NEK11','NEK5','NEK10','NIM1K','NEK4','NEK6','NRK','NEK7','SCYL1','NRBP2','PIK3C2G' , 'PAK7' , 'PANK3','PIK3C2B' , 'PAK3','PAK6','SCYL3','PDIK1L','PIP4K2C','PIP5K1A','PLK5','PHKG1','PHKG2','PI4KA','PIP5K1B','PKMYT1','PKN3' , 'PRPF4B','TP53RK','PRKRA','PRKY','PXK','RIOK1','RIOK3','RPS6KL1' , 'SBK2','SCYL2','SBK3','SGK494','POMK','SGK223','SRPK3' , 'STK32C','STK33','STK32B','STK40','STK17A','STK36','STK32A','STK38L','STKLD1','STK3','STK31','TESK2','TESK1','TBCK','TLK1' , 'TPK1','TLK2','TSSK6','TSSK1B','TTBK2','TSSK4','TTBK1','TSSK3','UCK1','ULK4','UCK2','WNK2','VRK2','WEE2' ),
          fontSize = 8)
dev.off()



#subset the data  




gg <- subsetMaf(maf = tcga_mc3, genes = c('ACVR1' , 'ADCK2' , 'ADCK5' , 'ADK' , 'ABL2' , 'AAK1' , 'TNK2' , 'ACVR1C' , 'ALK' , 'ACVR1B' , 'ACVRL1' , 'PRKAB1' , 'PRKAA1' , 'ABL1' , 'AMHR2' , 'ADCK1' , 'AKT1' , 
                                               'PRKAG1' , 'PRKAA2' , 'ARAF' , 'ALPK3' , 'NPR1' , 'ADPGK' , 'ANKK1' , 'NPR2' , 'ALPK1' , 'ALPK2' , 'GRK2' , 'GRK3' , 'AKT2' , 'AKT3' , 'ATM' , 'ATR' , 'AURKC' , 'ACVR2A' ,
                                               'ACVR2B' , 'BLK' , 'BCKDK' , 'AURKB' , 'AURKA' , 'BRAF' , 'BMPR1B' , 'BMX' , 'BTK' , 'CAMKV' , 'BMPR2' , 'BMPR1A' , 'BRSK1' , 'BRSK2' , 'BMP2K' , 'BUB1B' , 'BUB1' , 'CDC7' ,
                                               'CDK17' , 'CDK18' , 'CDK9' , 'CHKA' , 'CLK2' , 'CLK4' , 'CMPK2' , 'CHEK1' , 'CHKB' , 'CHEK2' , 'CDK2' , 'CDK7' , 'CDK8' , 'CDK10' , 'CDK11A' , 'CDK19' , 'CDK1' , 'CDK5' , 
                                               'CDKL4' , 'CDK11B' , 'CDK12' , 'CDK13' , 'CDK14' , 'CDK15' , 'CDK16' , 'CDK20' , 'CDK3' , 'CDK4' , 'CDK6' , 'CDKL1' , 'CDKL2' , 'CDKL3' , 'CDKL5' , 'CERK' , 'CSNK2A2' , 
                                               'COQ8A' , 'COQ8B' , 'CLK3' , 'CLK1' , 'CSNK2A1' , 'CSNK2A3' , 'CSNK2B' , 'CIT' , 'CSK' , 'CSF1R' , 'CASK' , 'DCLK3' , 'DAPK2' , 'DCLK1' , 'DDR2' , 'DAPK1' , 'DAPK3' , 
                                               'DCK' , 'DCLK2' , 'DGKQ' , 'DGKH' , 'DGUOK' , 'DMPK' , 'DGKA' , 'DDR1' , 'DYRK1A' , 'EIF2AK1' , 'EIF2AK4' , 'DYRK1B' , 'DYRK3' , 'EIF2AK2' , 'EGFR' , 'ETNK1' , 'ETNK2' , 
                                               'EPHB2' , 'ERN2' , 'DSTYK' , 'DYRK2' , 'DYRK4' , 'EIF2AK3' , 'ERBB3' , 'EPHA7' , 'EPHB1' , 'ERN1' , 'EEF2K' , 'EFNA1' , 'EPHA1' , 'EPHA8' , 'EPHA10' , 'PFKFB2' , 'PTK2' , 
                                               'EPHA6' , 'ERBB2' , 'ERBB4' , 'EPHA2' , 'EPHA3' , 'EPHA4' , 'EPHA5' , 'EPHB3' , 'EPHB4' , 'EPHB6' , 'FGFR2' , 'FGFR3' , 'FGFR4' , 'PTK2B' , 'FGFR1' , 'FGGY' , 'PFKFB1' ,
                                               'PFKFB3' , 'FASTK' , 'FER' , 'FN3K' , 'FES' , 'FUK' , 'FYN' , 'GAK' , 'FLT3' , 'PIKFYVE' , 'FRK' , 'FGR' , 'GALK1' , 'GALK2' , 'GK' , 'GLYCTK' , 'GRK7' , 'IDNK' , 'GUCY2C' , 'GK5' , 'GSK3A' , 'GSK3B' , 'GK3P' , 'GK2' , 'GRK6' , 'GUCY2F' , 'HASPIN' , 'GRK4' , 'GRK5' , 'HIPK4' , 'GUCY2D' , 'MASTL' , 'HCK' , 'HIPK1' , 'HIPK2' , 'HIPK3' , 'HKDC1' , 'HUNK' , 'HK2' , 'HK3' , 'ICK' , 'IGF1R' , 'ILK' , 'ITPKA' , 'ITPKC' , 'IP6K1' , 'IP6K2' , 'HK1' , 'CHUK' , 'IKBKB' , 'INSR' , 'ITK' , 'PRKACG' , 'CSNK1A1L' , 'CSNK1A1' , 'CSNK1D' , 'CAMK2D' , 'CKB' , 'IKBKE' , 'IP6K3' , 'JAK2' , 'JAK3' , 'AK1' , 'AK2' , 'KALRN' , 'JAK1' , 'AK3' , 'AK5' , 'AK6' , 'AK8' , 'AK9' , 'PRKACA' , 'CKM' , 'CKMT1A' , 'TK1' , 'AK4' , 'CAMK1D' , 'CAMK1G' , 'CMPK1' , 'GUK1' , 'PRKD3' , 'PSKH2' , 'INSRR' , 'IRAK1' , 'IRAK4' , 'PRKACB' , 'CSNK1G1' , 'CSNK1G3' , 'CAMK2B' , 'CAMK4' , 'GCK' , 'IPMK' , 'IRAK3' , 'ITPKB' , 'IRAK2' , 'ITPK1' , 'CSNK1E' , 'CSNK1G2' , 'CAMK1' , 'PNCK' , 'CAMK2A' , 'CAMK2G' , 'PRKG1' , 'PRKG2' , 'PSKH1' , 'RPS6KA1' , 'RPS6KB1' , 'LMTK3' , 'AATK' , 'MVK' , 'TK2' , 'KIT' , 'CAMKK2' , 'LCK' , 'RPS6KA2' , 'RPS6KA3' , 'RPS6KA4' , 'RPS6KA6' , 'DTYMK' , 'LMTK2' , 'KHK' , 'CAMKK1' , 'PHKA2' , 'PRKCD' , 'PRKCI' , 'PRKCH' , 'PRKCQ' , 'LATS2' , 'PRKCA' , 'PRKCB' , 'PRKCE' , 'PRKCG' , 'PRKCZ' , 'PKM' , 'PKLR' , 'RPS6KA5' , 'RPS6KC1' , 'KSR2' , 'SYK' , 'FN3KRP' , 'LATS1' , 'LIMK1' , 'LIMK2' , 'CKMT2' , 'PHKA1' , 'PRKD1' , 'PRKD2' , 'RPS6KB2' , 'KSR1' , 'LRRK2' , 'LYN' , 'MAP3K15' , 'MAP3K9' , 'MATK' , 'MET' , 'MAP3K10' , 'MAP3K11' , 'MAP3K13' , 'MAP3K14' , 'MAP3K21' , 'MAP3K4' , 'MAP3K5' , 'MAP3K7' , 'MAP4K1' , 'MAP4K5' , 'MAPKAPK3' , 'MAST1' , 'MAST2' , 'MAST4' , 'MAP3K12' , 'MAP3K2' , 'MAP4K3' , 'MAK' , 'MAP3K19' , 'MAP3K1' , 'MAP3K20' , 'MAP3K3' , 'MAP3K6' , 'MAP4K2' , 'MARK1' , 'MARK3' , 'MAPKAPK2' , 'MAPKAPK5' , 'MAST3' , 'LRRK1' , 'MAP4K4' , 'LTK' , 'MAP3K8' , 'MARK2' , 'MARK4' , 'MAPK14' , 'MLKL' , 'MOS' , 'MOK' , 'MAP2K2' , 'MAP2K6' , 'MAPK1' , 'MAPK3' , 'MAPK4' , 'MAPK7' , 'MAPK9' , 'MAP2K1' , 'MAP2K3' , 'MAP2K4' , 'MAP2K5' , 'MAP2K7' , 'MINK1' , 'MAPK15' , 'MKNK2' , 'CDC42BPG' , 'MELK' , 'MERTK' , 'MAPK6' , 'MAPK8' , 'MAPK10' , 'MAPK11' , 'MAPK12' , 'MAPK13' , 'MKNK1' , 'NME6' , 'MYLK4' , 'MYLK' , 'NME2P1' , 'NADK' , 'MYLK2' , 'MYLK3' , 'NADK2' , 'NME3' , 'NEK3' , 'NEK5' , 'NEK7' , 'MUSK' , 'MYO3A' , 'MYO3B' , 'CDC42BPA' , 'CDC42BPB' , 'MTOR' , 'NME4' , 'NEK8' , 'NEK9' , 'IKBKG' , 'NAGK' , 'NME5' , 'NME1' , 'NME1-NME2' , 'NEK11' , 'NIM1K' , 'NEK4' , 'NRBP2' , 'NTRK2' , 'NTRK3' , 'NEK10' , 'NEK1' , 'NEK2' , 'NEK6' , 'NLK' , 'NRBP1' , 'NTRK1' , 'NRK' , 'NUAK2' , 'PI4K2A' , 'PANK2' , 'PANK4' , 'NUAK1' , 'PAK3' , 'PAK5' , 'PANK1' , 'NMRK1' , 'NMRK2' , 'OBSCN' , 'SCYL3' , 'OXSR1' , 'PIK3C2A' , 'PIK3C2B' , 'PIK3C2G' , 'PIK3R1' , 'PAK1' , 'PAK4' , 'PDIK1L' , 'PDK1' , 'PDK2' , 'PDK3' , 'PDK4' , 'PDXK' , 'PFKM' , 'PDGFRB' , 'PCK1' , 'PCK2' , 'PAK2' , 'PAK6' , 'PAN3' , 'PANK3' , 'PDPK2P' , 'PEAK1' , 'PHKG1' , 'PHKG2' , 'PINK1' , 'PKDCC' , 'PDPK1' , 'PFKL' , 'PIM3' , 'PASK' , 'PDGFRA' , 'PGK1' , 'PGK2' , 'PRPS1' , 'PIM2' , 'PIK3C3' , 'PIK3CB' , 'PIK3CD' , 'PKN1' , 'PLK3' , 'PI4KB' , 'PIP5K1C' , 'PLK2' , 'PKMYT1' , 'PIK3R4' , 'PIP5K1A' , 'PIM1' , 'PKN2' , 'PLK1' , 'PLK5' , 'PIP4K2B' , 'PIP4K2C' , 'PIK3CA' , 
                                               'PFKP' , 'PIP4K2A' , 'PI4KA' , 'PIP5K1B' , 'PIK3CG' , 'PKN3' , 'PLK4' , 'PRKX' , 'PRAG1' , 'PRKRA' , 'PRPS2' , 'PTK6' , 'PRKDC' , 'PRKY' , 'PRPF4B' , 'PRPS1L1' , 
                                               'PTK7' , 'PSTK' , 'TP53RK' , 'PXK' , 'RAF1' , 'RIOK1' , 'RIOK2' , 'RBKS' , 'RIOK3' , 'RIPK2' , 'RIPK4' , 'RET' , 'RIPK1' , 'RIPK3' , 'GRK1' , 'ROR1' , 'RNASEL' , 
                                               'MST1R' , 'ROS1' , 'ROCK1' , 'ROCK2' , 'ROR2' , 'RPS6KL1' , 'SCYL1' , 'SCYL2' , 'SIK1' , 'RYK' , 'SHPK' , 'SBK1' , 'SBK2' , 'SBK3' , 'SNRK' , 'POMK' , 'SGK494' , 
                                               'SGK3' , 'SPHK2' , 'SLK' , 'SRC' , 'SPHK1' , 'SRPK3' , 'STK32C' , 'STK26' , 'STK35' , 'STKLD1' , 'SRMS' , 'SRPK1' , 'SRPK2' , 'SGK1' , 'SGK2' , 'SPEG' , 'STK38L' ,
                                               'SIK2' , 'SIK3' , 'SMG1' , 'STK11' , 'STK31' , 'STRADA' , 'STK17A' , 'STK38' , 'STK39' , 'STRADB' , 'STYK1' , 'STK24' , 'STK25' , 'STK33' , 'STK40' , 'STK4' , 'TAOK2' , 
                                               'TAOK3' , 'STK17B' , 'STK32A' , 'STK32B' , 'STK10' , 'STK16' , 'STK36' , 'STK3' , 'TAB1' , 'TEX14' , 'TAOK1' , 'TEC' , 'TTN' , 'TLK1' , 'TAF1' , 'TBCK' , 'TGFBR1' , 
                                               'TESK1' , 'TESK2' , 'TIE1' , 'TEK' , 'TBK1' , 'TNK1' , 'TSSK3' , 'TSSK4' , 'TYRO3' , 'TGFBR2' , 'TRIB2' , 'TNIK' , 'TRIB1' , 'TRIO' , 'TSSK1B' , 'TSSK6' , 'TTBK1' , 
                                               'TTBK2' , 'TLK2' , 'TPK1' , 'PBK' , 'TRIB3' , 'TRPM6' , 'TNNI3K' , 'ULK2' , 'ULK3' , 'ULK4' , 'KDR' , 'UCK2' , 'UCK1' , 'TYK2' , 'AXL' , 'UCKL1' , 'TXK' , 'UHMK1' ,
                                               'TSSK2' , 'TTK' , 'VRK1' , 'VRK3' , 'WEE2' , 'WNK4' , 'VRK2' , 'WEE1' , 'WNK3' , 'ULK1' , 'FLT1' , 'FLT4' , 'WNK1' , 'WNK2' , 'XYLB' , 'YES1' , 'ZAP70' ))



sort(table(gg$Hugo_Symbol))

gg1 <- data.frame(table(gg$Hugo_Symbol))

gg2 <- gg1[which(gg1$Freq > 5),]

print(dim(gg2))

#draw any number of genes of ur choice kinases
png(file =paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[qq], "/Mutation/kinases_", CanType[qq],".png", sep =""), width = 1000, height = 1000, units = "px")
oncostrip(maf = tcga_mc3, genes = as.character(gg2$Var1), fontSize = 8)
dev.off()





#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
#lollipopPlot(maf = tcga_kirc_mc3, gene = 'DNMT3A', AACol = 'GMAF', showMutationRate = TRUE)  #Protein_change for other cancers



#rainfallPlot(maf = tcga_kirc_mc3, detectChangePoints = TRUE, fontSize = 12, pointSize = 0.6)


#Detecting cancer driver genes based on positional clustering
tcga_mc3@maf.silent

if("HGVSp_Short" %in% names(tcga_mc3@maf.silent)){
sig = oncodrive(maf = tcga_mc3, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')

head(sig)
} else {
  
sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')

}


sig1 <- select(sig, Hugo_Symbol, clusterScores, protLen, zscore,  pval, fdr, fract_muts_in_clusters)



sig1 <- sig1[which(sig1$Hugo_Symbol  %in% gg1$Var1), ]

sig1$cancer <- CanType[qq]

jj <- select(sig1, Hugo_Symbol, pval, fdr)
colnames(jj) <- c("Hugo_Symbol", paste("pval_",CanType[qq], sep = ""), paste("fdr_", CanType[qq], sep =""))

mergeJJ <- merge(as.data.frame(mergeJJ), as.data.frame(jj), by='Hugo_Symbol', all=TRUE)





png(file = paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[qq], "/Mutation/plotOncodrive_", CanType[qq], ".png", sep =""), width = 1000, height = 1000, units = "px")
plotOncodrive(res = sig, fdrCutOff = 0.1, useFraction = TRUE)
dev.off()

save(tcga_mc3, sig, gg1, jj, file = paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[qq], "/Mutation/AllDatawithsig.RData", sep =""))


print(paste("######################## Ended with ", CanType[qq], "######################", sep = "    "))
}

write.csv(mergeJJ, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/Mutation_Pvalue.csv")



################################################
#scoring for stages in each cancer 
################################################

#read the files

#IDG_kinase <- read.csv("/projects/ccs/schurerlab/Rimpi/kinases_IDG.txt", sep="\t")
Dir <- "/projects/ccs/schurerlab/DerekJ/TCGA/Cell_ClinicalMerged/"
head(Kinases)

#no stages "CESC", "GBM" ,"PCPG", "PRAD","UCEC"

CanType <- c("BLCA",  "KIRP","KIRC", "BRCA",  "LUAD", "LUSC",  "THCA", "CESC",  "CHOL", "KICH",   "READ" , "LGG","ACC","SKCM","PAAD","TGCT","UCS", "UVM")     #"BLCA", "KICH", "KIRP"


#"COAD", "ESCA",  "HNSC" , "LIHC","STAD","GBM" ,"PRAD", "UCEC" ,"PCPG","SARC","THYM",

#does noit have pathlogical stage data
#GBM, PRAD,PCPG, SARC, THYM, OV", "UVM",

#clinical stages
#UCEC - clinical stage (does not have pathlogical stage)
#THYM
#OV
#UVM
#CESC



#TGCT different type of stage 1S


#Save files in 

save <- "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/stagewisePvalueForeachCancer/Kinases/"


################################################
#for loop to call all the files and and calculate the P-value for each stage using 
################################################


Anova_stage <- data.frame()
Anova_grade <- data.frame()
for(i in 1:length(CanType)){
  print(CanType[i] )
  
  a <- read.csv(paste(Dir, "TCGA_", CanType[i], "_cellclinicalmerge1.csv", sep="")) 
  table(a$ajcc_pathologic_tumor_stage)
  
  b <- a[, colnames(a) %in% Kinases]
  bcr_patient_barcode <- data.frame(bcr_patient_barcode = paste(a$TCGA_ID, "-", a$Sample, sep =""))
  
  c <- select(a, type, ajcc_pathologic_tumor_stage , histological_grade)
  c$bcr_patient_barcode <- bcr_patient_barcode

  f <- as.character(c$ajcc_pathologic_tumor_stage)
  f[f == "Stage IA" ] <- "Stage I" 
  f[f == "Stage IB" ] <- "Stage I"
  f[f == "Stage IIA" ] <- "Stage II" 
  f[f == "Stage IIB" ] <- "Stage II" 
  f[f == "Stage IIIB" ] <- "Stage III" 
  f[f == "Stage IIIA" ] <- "Stage III" 
  f[f == "Stage IIIC" ] <- "Stage III"
  f[f == "Stage IVA" ] <- "Stage IV"
  f[f == "Stage IVB" ] <- "Stage IV"
  c$ajcc_pathologic_tumor_stage1 <- f
  c <- cbind(c, b)
  d <- subset(c, subset = ajcc_pathologic_tumor_stage1 == "Stage I" | ajcc_pathologic_tumor_stage1 == "Stage II" | ajcc_pathologic_tumor_stage1 == "Stage III"| ajcc_pathologic_tumor_stage1 == "Stage IV" )
  e <- subset(c, subset = histological_grade == "G1" | histological_grade == "G2"| histological_grade == "G3" | histological_grade == "G4"  )
  genes <- as.character(colnames(b))
  
  
  Anova_stage <- data.frame()
  Anova_grade <- data.frame()
  for(j in 6:ncol(d)){
    
    print(colnames(d[j]))
    
    
    formula <- as.formula(paste(colnames(d)[j], " ~ ajcc_pathologic_tumor_stage1", sep=""))
    
    res.aov <- aov(formula, data = d)
    #res.aov <- aov(colnames(c[j]) ~ ajcc_pathologic_tumor_stage, data = c)
    
    
    
    #res.aov_grade <- aov(genes[,j] ~ histological_grade, data = c)
    
    
    
    print(paste("Anova for the gene == ", colnames(c[j])))
    print(TukeyHSD(res.aov))
    
    res.aov_result <-  TukeyHSD(res.aov)[1]
    res.aov_result <- data.frame(res.aov_result$ajcc_pathologic_tumor_stage1)
    
    res.aov_result$gene <- colnames(d[j])
    Anova_s <- res.aov_result
    
    Anova_stage <- rbind(Anova_stage, Anova_s)
    
    print(paste("*************gene is done ", colnames(d[j]), "*****************"))
    
    if(length(e$histological_grade != 0)){
    #if(e$histological_grade == 0){
      
    for(k in 6:ncol(e)){
      #for the grade 
      formula1 <- as.formula(paste(colnames(e)[k], " ~ histological_grade", sep=""))
      res.aov_grade <- aov(formula1, data = e)
      
      res.aov_grade_result <-  TukeyHSD(res.aov_grade)[1]
      res.aov_grade_result <-  data.frame(res.aov_grade_result$histological_grade)
      
      res.aov_grade_result$gene <- colnames(e[k]) 
      Anova_g <- res.aov_grade_result
      Anova_grade <- rbind(Anova_grade, Anova_g)
    }
    
    }else{
      
      print(paste("*************histological_grade is not Avaiable:::::", CanType[i], "***************" , sep =""))
    } 
    print(paste("*************the canceris done::::::", CanType[i], "****************", sep ="" ))
  }
  
  write.csv(Anova_stage, file = paste(save, CanType[i], "Stage.csv"))
  length(unique(Anova_stage$gene)) == length(genes)
  
  write.csv(Anova_grade, file = paste(save, CanType[i], "Grade.csv"))
  length(unique(Anova_grade$gene)) == length(genes)
  
}






################################################
#Survival analysis for those cancers which do not have control sample 
################################################

library(survival)
library(RTCGA.rnaseq)
library(RTCGA.clinical)
library(dplyr)
library(tidyr)
library(survminer)

#
CanType <-  c("LGG","ACC","SARC","THYM","OV","UVM","SKCM","PAAD","TGCT","UCS")   #"CHOL", "KICH", "UCEC" ,  "READ" ,"PCPG" #

Kinases <- as.character(Kinases)

for(i in 7:length(CanType)){
  
  print(paste("Processing the Survival Analysis for", CanType[i]))
  clinical <- readRDS(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[i], "/TCGA-", CanType[i],"_clinical.rds", sep=""))
  
  
  
  
  #extract attributes needed for survival analysis
  clin <- data.frame(patient.bcr_patient_barcode = clinical$gdc_cases.submitter_id,
                     patient.vital_status = clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = clinical$cgc_case_days_to_last_follow_up, 
                     patient.days_to_death = clinical$xml_days_to_death)
  
  load(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/",CanType[i], "_Normalize_filter.RData", sep=""))
  
  #clinical.surv <- survivalTCGA(get(paste(CanType[i],".clin", sep="")))
  clinical.surv <- survivalTCGA(clin)
  
  norm <- get(paste(CanType[i],".Normalize.expression", sep =""))
  
  norm_filter <- norm[grep("^c", rownames(norm)),]
  
  norm_filter1 <- norm[-which(rownames(norm) %in% rownames(norm_filter)), ]
  
  #split the string with 
  split1 <- unlist(strsplit(as.character(rownames(norm_filter)), "[()]"))
  n <- length(split1) 
  split2 <- split1[seq(2, n, 2)] 
  split3 <- gsub(",.*","",split2)
  split3 <- gsub('[[:punct:]]',"",split3) 
  #split3 <- gsub("[^[:alnum:][:blank:]+]", "", split2)
  
  
  #split4 <- gsub("[ ]", "/", split3)
  
  # split2 <- noquote(split2)
  # split3 <- gsub(",.*","",split2)
  #split3 <- gsub("\\,", "/", split2)
  #remove everything after , 
  # split3 <- gsub('[[:punct:]]',"",split3) 
  
  rownames(norm_filter) <- split3
  norm_filter2 <- rbind(norm_filter1, norm_filter)
  dim(norm_filter2) == dim(norm)
  
  
  rnaseq <- data.frame(t(norm_filter2))
  rnaseq$bcr_patient_barcode <- rownames(rnaseq)
  rownames(rnaseq) <- NULL
  
  #"NME1-NME2" = "NME1-NME2" , "HASPIN" = "HASPIN" ,"MAP3K21" = "MAP3K21" ,"PDPK2P" = "PDPK2P" ,"PRAG1" = "PRAG1" ,"NME2P1" = "NME2P1" ,"MAP3K20" = "MAP3K20" ,
  
  expressionsTCGA(
    rnaseq,
    extract.cols = Kinases) %>%
    rename(cohort = dataset,
           "ACVR1" = "ACVR1" , "ADCK2" = "ADCK2" , "ADCK5" = "ADCK5" , "ADK" = "ADK" , "ABL2" = "ABL2" , "AAK1" = "AAK1" , 
           "TNK2" = "TNK2" , "ACVR1C" = "ACVR1C" , "ALK" = "ALK" , "ACVR1B" = "ACVR1B" , "ACVRL1" = "ACVRL1" , "PRKAB1" = "PRKAB1" ,
           "PRKAA1" = "PRKAA1" , "ABL1" = "ABL1" , "AMHR2" = "AMHR2" , "ADCK1" = "ADCK1" , "AKT1" = "AKT1" , "PRKAG1" = "PRKAG1" ,
           "PRKAA2" = "PRKAA2" , "ARAF" = "ARAF" , "ALPK3" = "ALPK3" , "NPR1" = "NPR1" , "ADPGK" = "ADPGK" , "ANKK1" = "ANKK1" , 
           "NPR2" = "NPR2" , "ALPK1" = "ALPK1" , "ALPK2" = "ALPK2" , "GRK2" = "GRK2" , "GRK3" = "GRK3" , "AKT2" = "AKT2" , 
           "AKT3" = "AKT3" , "ATM" = "ATM" , "ATR" = "ATR" , "AURKC" = "AURKC" , "ACVR2A" = "ACVR2A" , "ACVR2B" = "ACVR2B" ,
           "BLK" = "BLK" , "BCKDK" = "BCKDK" , "AURKB" = "AURKB" , "AURKA" = "AURKA" , "BRAF" = "BRAF" , "BMPR1B" = "BMPR1B" , 
           "BMX" = "BMX" , "BTK" = "BTK" , "CAMKV" = "CAMKV" , "BMPR2" = "BMPR2" , "BMPR1A" = "BMPR1A" , "BRSK1" = "BRSK1" , 
           "BRSK2" = "BRSK2" , "BMP2K" = "BMP2K" , "BUB1B" = "BUB1B" , "BUB1" = "BUB1" , "CDC7" = "CDC7" , "CDK17" = "CDK17" , 
           "CDK18" = "CDK18" , "CDK9" = "CDK9" , "CHKA" = "CHKA" , "CLK2" = "CLK2" , "CLK4" = "CLK4" , "CMPK2" = "CMPK2" , 
           "CHEK1" = "CHEK1" , "CHKB" = "CHKB" , "CHEK2" = "CHEK2" , "CDK2" = "CDK2" , "CDK7" = "CDK7" , "CDK8" = "CDK8" , 
           "CDK10" = "CDK10" , "CDK11A" = "CDK11A" , "CDK19" = "CDK19" , "CDK1" = "CDK1" , "CDK5" = "CDK5" , "CDKL4" = "CDKL4" , 
           "CDK11B" = "CDK11B" , "CDK12" = "CDK12" , "CDK13" = "CDK13" , "CDK14" = "CDK14" , "CDK15" = "CDK15" , "CDK16" = "CDK16" , 
           "CDK20" = "CDK20" , "CDK3" = "CDK3" , "CDK4" = "CDK4" , "CDK6" = "CDK6" , "CDKL1" = "CDKL1" , "CDKL2" = "CDKL2" , 
           "CDKL3" = "CDKL3" , "CDKL5" = "CDKL5" , "CERK" = "CERK" , "CSNK2A2" = "CSNK2A2" , "COQ8A" = "COQ8A" , "COQ8B" = "COQ8B" ,
           "CLK3" = "CLK3" , "CLK1" = "CLK1" , "CSNK2A1" = "CSNK2A1" , "CSNK2A3" = "CSNK2A3" , "CSNK2B" = "CSNK2B" , "CIT" = "CIT" , "CSK" = "CSK" , "CSF1R" = "CSF1R" , "CASK" = "CASK" , "DCLK3" = "DCLK3" , "DAPK2" = "DAPK2" , "DCLK1" = "DCLK1" , "DDR2" = "DDR2" , "DAPK1" = "DAPK1" , "DAPK3" = "DAPK3" , "DCK" = "DCK" , "DCLK2" = "DCLK2" , "DGKQ" = "DGKQ" ,
           "DGKH" = "DGKH" , "DGUOK" = "DGUOK" , "DMPK" = "DMPK" , "DGKA" = "DGKA" , "DDR1" = "DDR1" , "DYRK1A" = "DYRK1A" , "EIF2AK1" = "EIF2AK1" , "EIF2AK4" = "EIF2AK4" , "DYRK1B" = "DYRK1B" , "DYRK3" = "DYRK3" , "EIF2AK2" = "EIF2AK2" , "EGFR" = "EGFR" , "ETNK1" = "ETNK1" , "ETNK2" = "ETNK2" , "EPHB2" = "EPHB2" , "ERN2" = "ERN2" , "DSTYK" = "DSTYK" ,
           "DYRK2" = "DYRK2" , "DYRK4" = "DYRK4" , "EIF2AK3" = "EIF2AK3" , "ERBB3" = "ERBB3" , "EPHA7" = "EPHA7" , "EPHB1" = "EPHB1" , "ERN1" = "ERN1" , "EEF2K" = "EEF2K" , "EFNA1" = "EFNA1" , "EPHA1" = "EPHA1" , "EPHA8" = "EPHA8" , "EPHA10" = "EPHA10" , "PFKFB2" = "PFKFB2" , "PTK2" = "PTK2" , "EPHA6" = "EPHA6" , "ERBB2" = "ERBB2" , "ERBB4" = "ERBB4" ,
           "EPHA2" = "EPHA2" , "EPHA3" = "EPHA3" , "EPHA4" = "EPHA4" , "EPHA5" = "EPHA5" , "EPHB3" = "EPHB3" , "EPHB4" = "EPHB4" , "EPHB6" = "EPHB6" , "FGFR2" = "FGFR2" , "FGFR3" = "FGFR3" , "FGFR4" = "FGFR4" , "PTK2B" = "PTK2B" , "FGFR1" = "FGFR1" , "FGGY" = "FGGY" , "PFKFB1" = "PFKFB1" , "PFKFB3" = "PFKFB3" , "FASTK" = "FASTK" , "FER" = "FER" ,
           "FN3K" = "FN3K" , "FES" = "FES" , "FUK" = "FUK" , "FYN" = "FYN" , "GAK" = "GAK" , "FLT3" = "FLT3" , "PIKFYVE" = "PIKFYVE" , "FRK" = "FRK" , "FGR" = "FGR" , "GALK1" = "GALK1" , "GALK2" = "GALK2" , "GK" = "GK" , "GLYCTK" = "GLYCTK" , "GRK7" = "GRK7" , "IDNK" = "IDNK" , "GUCY2C" = "GUCY2C" , 
           "GK5" = "GK5" , "GSK3A" = "GSK3A" , "GSK3B" = "GSK3B" , "GK3P" = "GK3P" , "GK2" = "GK2" , "GRK6" = "GRK6" , "GUCY2F" = "GUCY2F" ,  "GRK4" = "GRK4" , "GRK5" = "GRK5" , "HIPK4" = "HIPK4" , "GUCY2D" = "GUCY2D" , "MASTL" = "MASTL" , "HCK" = "HCK" , "HIPK1" = "HIPK1" , "HIPK2" = "HIPK2" , "HIPK3" = "HIPK3" , "HKDC1" = "HKDC1" ,
           "HUNK" = "HUNK" , "HK2" = "HK2" , "HK3" = "HK3" , "ICK" = "ICK" , "IGF1R" = "IGF1R" , "ILK" = "ILK" , "ITPKA" = "ITPKA" , "ITPKC" = "ITPKC" , "IP6K1" = "IP6K1" , "IP6K2" = "IP6K2" , "HK1" = "HK1" , "CHUK" = "CHUK" , "IKBKB" = "IKBKB" , "INSR" = "INSR" , "ITK" = "ITK" , "PRKACG" = "PRKACG" , "CSNK1A1L" = "CSNK1A1L" , "CSNK1A1" = "CSNK1A1" ,
           "CSNK1D" = "CSNK1D" , "CAMK2D" = "CAMK2D" , "CKB" = "CKB" , "IKBKE" = "IKBKE" , "IP6K3" = "IP6K3" , "JAK2" = "JAK2" , "JAK3" = "JAK3" , "AK1" = "AK1" , "AK2" = "AK2" , "KALRN" = "KALRN" , "JAK1" = "JAK1" , "AK3" = "AK3" , "AK5" = "AK5" , "AK6" = "AK6" , "AK8" = "AK8" , "AK9" = "AK9" , "PRKACA" = "PRKACA" , "CKM" = "CKM" , "CKMT1A" = "CKMT1A" , "TK1" = "TK1" , "AK4" = "AK4" , "CAMK1D" = "CAMK1D" , "CAMK1G" = "CAMK1G" , "CMPK1" = "CMPK1" , "GUK1" = "GUK1" , "PRKD3" = "PRKD3" , "PSKH2" = "PSKH2" , "INSRR" = "INSRR" , "IRAK1" = "IRAK1" , "IRAK4" = "IRAK4" , "PRKACB" = "PRKACB" , "CSNK1G1" = "CSNK1G1" , "CSNK1G3" = "CSNK1G3" , "CAMK2B" = "CAMK2B" , "CAMK4" = "CAMK4" , "GCK" = "GCK" , "IPMK" = "IPMK" , "IRAK3" = "IRAK3" , "ITPKB" = "ITPKB" , "IRAK2" = "IRAK2" , "ITPK1" = "ITPK1" , "CSNK1E" = "CSNK1E" , "CSNK1G2" = "CSNK1G2" , "CAMK1" = "CAMK1" , 
           "PNCK" = "PNCK" , "CAMK2A" = "CAMK2A" , "CAMK2G" = "CAMK2G" , "PRKG1" = "PRKG1" , "PRKG2" = "PRKG2" , "PSKH1" = "PSKH1" , "RPS6KA1" = "RPS6KA1" , "RPS6KB1" = "RPS6KB1" , "LMTK3" = "LMTK3" , "AATK" = "AATK" , "MVK" = "MVK" , "TK2" = "TK2" , "KIT" = "KIT" , "CAMKK2" = "CAMKK2" , "LCK" = "LCK" , "RPS6KA2" = "RPS6KA2" , "RPS6KA3" = "RPS6KA3" , "RPS6KA4" = "RPS6KA4" , "RPS6KA6" = "RPS6KA6" , "DTYMK" = "DTYMK" , "LMTK2" = "LMTK2" , "KHK" = "KHK" , "CAMKK1" = "CAMKK1" , "PHKA2" = "PHKA2" , "PRKCD" = "PRKCD" , "PRKCI" = "PRKCI" , "PRKCH" = "PRKCH" , "PRKCQ" = "PRKCQ" , "LATS2" = "LATS2" , "PRKCA" = "PRKCA" , "PRKCB" = "PRKCB" , "PRKCE" = "PRKCE" , "PRKCG" = "PRKCG" , "PRKCZ" = "PRKCZ" , "PKM" = "PKM" , "PKLR" = "PKLR" , "RPS6KA5" = "RPS6KA5" , "RPS6KC1" = "RPS6KC1" , "KSR2" = "KSR2" , "SYK" = "SYK" , "FN3KRP" = "FN3KRP" , "LATS1" = "LATS1" , "LIMK1" = "LIMK1" , "LIMK2" = "LIMK2" , "CKMT2" = "CKMT2" , "PHKA1" = "PHKA1" , "PRKD1" = "PRKD1" , "PRKD2" = "PRKD2" , "RPS6KB2" = "RPS6KB2" , "KSR1" = "KSR1" , "LRRK2" = "LRRK2" , "LYN" = "LYN" , "MAP3K15" = "MAP3K15" , "MAP3K9" = "MAP3K9" , "MATK" = "MATK" , "MET" = "MET" , "MAP3K10" = "MAP3K10" , "MAP3K11" = "MAP3K11" , "MAP3K13" = "MAP3K13" , "MAP3K14" = "MAP3K14" ,  "MAP3K4" = "MAP3K4" , "MAP3K5" = "MAP3K5" , "MAP3K7" = "MAP3K7" , "MAP4K1" = "MAP4K1" , "MAP4K5" = "MAP4K5" , "MAPKAPK3" = "MAPKAPK3" , "MAST1" = "MAST1" , "MAST2" = "MAST2" , "MAST4" = "MAST4" , "MAP3K12" = "MAP3K12" , "MAP3K2" = "MAP3K2" , "MAP4K3" = "MAP4K3" , "MAK" = "MAK" , "MAP3K19" = "MAP3K19" , "MAP3K1" = "MAP3K1" ,  "MAP3K3" = "MAP3K3" , "MAP3K6" = "MAP3K6" , "MAP4K2" = "MAP4K2" , "MARK1" = "MARK1" , "MARK3" = "MARK3" , "MAPKAPK2" = "MAPKAPK2" , "MAPKAPK5" = "MAPKAPK5" , "MAST3" = "MAST3" , "LRRK1" = "LRRK1" , "MAP4K4" = "MAP4K4" , "LTK" = "LTK" , "MAP3K8" = "MAP3K8" , "MARK2" = "MARK2" , "MARK4" = "MARK4" , "MAPK14" = "MAPK14" , "MLKL" = "MLKL" , "MOS" = "MOS" , "MOK" = "MOK" , "MAP2K2" = "MAP2K2" , "MAP2K6" = "MAP2K6" , "MAPK1" = "MAPK1" , "MAPK3" = "MAPK3" , "MAPK4" = "MAPK4" , "MAPK7" = "MAPK7" , "MAPK9" = "MAPK9" , "MAP2K1" = "MAP2K1" , "MAP2K3" = "MAP2K3" , "MAP2K4" = "MAP2K4" , "MAP2K5" = "MAP2K5" , "MAP2K7" = "MAP2K7" , "MINK1" = "MINK1" , "MAPK15" = "MAPK15" , "MKNK2" = "MKNK2" , "CDC42BPG" = "CDC42BPG" , "MELK" = "MELK" , "MERTK" = "MERTK" , "MAPK6" = "MAPK6" , "MAPK8" = "MAPK8" , "MAPK10" = "MAPK10" , "MAPK11" = "MAPK11" , "MAPK12" = "MAPK12" , "MAPK13" = "MAPK13" , "MKNK1" = "MKNK1" , "NME6" = "NME6" , "MYLK4" = "MYLK4" , "MYLK" = "MYLK" ,  "NADK" = "NADK" , "MYLK2" = "MYLK2" , "MYLK3" = "MYLK3" , "NADK2" = "NADK2" , "NME3" = "NME3" , "NEK3" = "NEK3" , "NEK5" = "NEK5" , "NEK7" = "NEK7" , "MUSK" = "MUSK" , "MYO3A" = "MYO3A" , "MYO3B" = "MYO3B" , "CDC42BPA" = "CDC42BPA" ,
           "CDC42BPB" = "CDC42BPB" , "MTOR" = "MTOR" , "NME4" = "NME4" , "NEK8" = "NEK8" , "NEK9" = "NEK9" , "IKBKG" = "IKBKG" , "NAGK" = "NAGK" , "NME5" = "NME5" , "NME1" = "NME1" ,  "NEK11" = "NEK11" , "NIM1K" = "NIM1K" , "NEK4" = "NEK4" , "NRBP2" = "NRBP2" , "NTRK2" = "NTRK2" , "NTRK3" = "NTRK3" , "NEK10" = "NEK10" , "NEK1" = "NEK1" , "NEK2" = "NEK2" , "NEK6" = "NEK6" , "NLK" = "NLK" , "NRBP1" = "NRBP1" , "NTRK1" = "NTRK1" , "NRK" = "NRK" , "NUAK2" = "NUAK2" , "PI4K2A" = "PI4K2A" , "PANK2" = "PANK2" , "PANK4" = "PANK4" , "NUAK1" = "NUAK1" , "PAK3" = "PAK3" , "PAK5" = "PAK5" , "PANK1" = "PANK1" , "NMRK1" = "NMRK1" , "NMRK2" = "NMRK2" , "OBSCN" = "OBSCN" , "SCYL3" = "SCYL3" , "OXSR1" = "OXSR1" , "PIK3C2A" = "PIK3C2A" , "PIK3C2B" = "PIK3C2B" , "PIK3C2G" = "PIK3C2G" , "PIK3R1" = "PIK3R1" , "PAK1" = "PAK1" , "PAK4" = "PAK4" , "PDIK1L" = "PDIK1L" , "PDK1" = "PDK1" , "PDK2" = "PDK2" , "PDK3" = "PDK3" , "PDK4" = "PDK4" , "PDXK" = "PDXK" , "PFKM" = "PFKM" , "PDGFRB" = "PDGFRB" , "PCK1" = "PCK1" , "PCK2" = "PCK2" , "PAK2" = "PAK2" , "PAK6" = "PAK6" , "PAN3" = "PAN3" , "PANK3" = "PANK3" ,  
           "PEAK1" = "PEAK1" , "PHKG1" = "PHKG1" , "PHKG2" = "PHKG2" , "PINK1" = "PINK1" , "PKDCC" = "PKDCC" , "PDPK1" = "PDPK1" , "PFKL" = "PFKL" , "PIM3" = "PIM3" , "PASK" = "PASK" , "PDGFRA" = "PDGFRA" , "PGK1" = "PGK1" , "PGK2" = "PGK2", 
           "PRPS1" = "PRPS1" , "PIM2" = "PIM2" , "PIK3C3" = "PIK3C3" , "PIK3CB" = "PIK3CB" , "PIK3CD" = "PIK3CD" , "PKN1" = "PKN1" , "PLK3" = "PLK3" , "PI4KB" = "PI4KB" , "PIP5K1C" = "PIP5K1C" , "PLK2" = "PLK2" , "PKMYT1" = "PKMYT1" , 
           "PIK3R4" = "PIK3R4" , "PIP5K1A" = "PIP5K1A" , "PIM1" = "PIM1" , "PKN2" = "PKN2" , "PLK1" = "PLK1" , "PLK5" = "PLK5" , "PIP4K2B" = "PIP4K2B" , "PIP4K2C" = "PIP4K2C" , "PIK3CA" = "PIK3CA" , "PFKP" = "PFKP" , "PIP4K2A" = "PIP4K2A" , "PI4KA" = "PI4KA" , "PIP5K1B" = "PIP5K1B" , "PIK3CG" = "PIK3CG" , "PKN3" = "PKN3" , "PLK4" = "PLK4" , "PRKX" = "PRKX" , 
           "PRKRA" = "PRKRA" , "PRPS2" = "PRPS2" , "PTK6" = "PTK6" , "PRKDC" = "PRKDC" , "PRKY" = "PRKY" , "PRPF4B" = "PRPF4B" , "PRPS1L1" = "PRPS1L1" , "PTK7" = "PTK7" , "PSTK" = "PSTK" , "TP53RK" = "TP53RK" , "PXK" = "PXK" , "RAF1" = "RAF1" , "RIOK1" = "RIOK1" , "RIOK2" = "RIOK2" , "RBKS" = "RBKS" , "RIOK3" = "RIOK3" , "RIPK2" = "RIPK2" , "RIPK4" = "RIPK4" , "RET" = "RET" , "RIPK1" = "RIPK1" , "RIPK3" = "RIPK3" , "GRK1" = "GRK1" , "ROR1" = "ROR1" , "RNASEL" = "RNASEL" , "MST1R" = "MST1R" , "ROS1" = "ROS1" , "ROCK1" = "ROCK1" , "ROCK2" = "ROCK2" , "ROR2" = "ROR2" , "RPS6KL1" = "RPS6KL1" , 
           "SCYL1" = "SCYL1" , "SCYL2" = "SCYL2" , "SIK1" = "SIK1" , "RYK" = "RYK" , "SHPK" = "SHPK" , "SBK1" = "SBK1" , "SBK2" = "SBK2" , "SBK3" = "SBK3" , "SNRK" = "SNRK" , "POMK" = "POMK" , "SGK494" = "SGK494" , "SGK3" = "SGK3" , "SPHK2" = "SPHK2" , "SLK" = "SLK" , "SRC" = "SRC" , "SPHK1" = "SPHK1" , "SRPK3" = "SRPK3" , "STK32C" = "STK32C" , "STK26" = "STK26" , "STK35" = "STK35" , "STKLD1" = "STKLD1" , "SRMS" = "SRMS" , "SRPK1" = "SRPK1" , "SRPK2" = "SRPK2" , "SGK1" = "SGK1" , "SGK2" = "SGK2" , "SPEG" = "SPEG" , "STK38L" = "STK38L" , "SIK2" = "SIK2" , "SIK3" = "SIK3" , "SMG1" = "SMG1" , "STK11" = "STK11" , "STK31" = "STK31" , "STRADA" = "STRADA" , "STK17A" = "STK17A" , "STK38" = "STK38" , "STK39" = "STK39" , "STRADB" = "STRADB" , "STYK1" = "STYK1" , "STK24" = "STK24" , "STK25" = "STK25" , "STK33" = "STK33" , "STK40" = "STK40" , "STK4" = "STK4" , "TAOK2" = "TAOK2" , "TAOK3" = "TAOK3" , "STK17B" = "STK17B" , "STK32A" = "STK32A" , "STK32B" = "STK32B" , "STK10" = "STK10" , "STK16" = "STK16" , "STK36" = "STK36" , "STK3" = "STK3" , "TAB1" = "TAB1" , "TEX14" = "TEX14" , "TAOK1" = "TAOK1" , "TEC" = "TEC" , "TTN" = "TTN" , "TLK1" = "TLK1" , "TAF1" = "TAF1" , "TBCK" = "TBCK" , "TGFBR1" = "TGFBR1" , "TESK1" = "TESK1" , "TESK2" = "TESK2" , "TIE1" = "TIE1" , "TEK" = "TEK" , "TBK1" = "TBK1" , "TNK1" = "TNK1" , "TSSK3" = "TSSK3" , "TSSK4" = "TSSK4" , "TYRO3" = "TYRO3" , "TGFBR2" = "TGFBR2" , "TRIB2" = "TRIB2" , "TNIK" = "TNIK" , "TRIB1" = "TRIB1" , "TRIO" = "TRIO" , "TSSK1B" = "TSSK1B" , "TSSK6" = "TSSK6" , "TTBK1" = "TTBK1" , "TTBK2" = "TTBK2" , "TLK2" = "TLK2" , "TPK1" = "TPK1" , "PBK" = "PBK" , "TRIB3" = "TRIB3" , "TRPM6" = "TRPM6" , "TNNI3K" = "TNNI3K" , "ULK2" = "ULK2" , "ULK3" = "ULK3" , "ULK4" = "ULK4" , "KDR" = "KDR" , "UCK2" = "UCK2" , "UCK1" = "UCK1" , "TYK2" = "TYK2" , "AXL" = "AXL" , "UCKL1" = "UCKL1" , "TXK" = "TXK" , "UHMK1" = "UHMK1" , "TSSK2" = "TSSK2" , "TTK" = "TTK" , "VRK1" = "VRK1" , "VRK3" = "VRK3" , "WEE2" = "WEE2" , "WNK4" = "WNK4" , "VRK2" = "VRK2" , "WEE1" = "WEE1" , "WNK3" = "WNK3" ,
           "ULK1" = "ULK1" , "FLT1" = "FLT1" , "FLT4" = "FLT4" , "WNK1" = "WNK1" , 
           "WNK2" = "WNK2" , "XYLB" = "XYLB" , "YES1" = "YES1" , "ZAP70" = "ZAP70") %>%
    filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
    # only cancer samples
    mutate(bcr_patient_barcode = 
             substr(bcr_patient_barcode, 1, 12)) -> rnaseq1 
  
  
  
  
  #`CDK9`, `CIT`, `DAPK3`, `DDR1`, `EEF2K`, `HASPIN`, `IRAK1`, `IRAK3`, `CSNK1E`, `MAP3K21`, `MAP3K20`, `NME2P1`, `NME1-NME2`, `NRBP2`, `PDXK`, `PAK6`, `PDPK2P`, `PINK1`, `PRAG1`, `SPEG`, `TAOK1`, `UCK2`, `UCKL1`
  
  clinical.surv %>%
    left_join(rnaseq1,
              by = "bcr_patient_barcode") ->
    surv_rnaseq
  
  
  table(surv_rnaseq$cohort, useNA = "always")
  
  
  surv_rnaseq <- surv_rnaseq %>%
    filter(!is.na(cohort))
  
  
  dim(surv_rnaseq)
  
  
  ######################
  #GENES for survival analysis
  #######################
  
  library(survminer)
  
  gen <- as.character(colnames(surv_rnaseq))
  gen <- gen[-c(1:4)]
  
  gen <- gen[-grep("^GK2$|^PSKH2$|^MOS$", gen)]
  
  setwd(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[i], "/survival/",  sep=""  ))
  #dir.create(file.path(paste0(getwd(), "/Kinase/")))
  
  ######################
  #Plot SURVIVAL CURVES FOR GENES
  #######################
  p <- data.frame()
  #BLCA- #GK2 -180, PSKH2-237, 346-MOS
  #KIRP - GK2 -180, PSKH2-237, 346-MOS
  
  for(j in 1:length(gen)){
    
    print(paste("Survival Analysis for", gen[j]))
    
    surv_rnaseq.cut <- surv_cutpoint(
      surv_rnaseq,
      time = "times",
      event = "patient.vital_status",
      variables = c(paste(gen[j]), "cohort")
    )
    summary(surv_rnaseq.cut)
    
    #distri <- plot(LIHC.surv_rnaseq.cut, geneid[i], palette = "npg")
    
    #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
    #print(distri)
    
    
    #dev.off()
    
    surv_rnaseq.cat <- surv_categorize(surv_rnaseq.cut)
    
    
    
    fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , gen[j], "+ cohort")),
                   data = surv_rnaseq.cat)
    
    pval <- surv_pvalue(fit)
    
    p <- rbind(p, pval) 
    
    
    save(pval,fit,surv_rnaseq.cat, rnaseq1,clinical.surv, file= paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[i], "/survival/Kinase/", gen[j], "_", CanType[i], "_Pvalue.RData", sep ="") )
    
    png(file = paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[i], "/survival/Kinase/", gen[j], ".survival.jpg", sep=""))
    gg <-ggsurvplot(
      fit, 
      data= surv_rnaseq.cat,
      break.time.by = 500, 
      ggtheme = theme_bw(), 
      risk.table = TRUE,
      xlim = c(0,3000),
      xlab = "Time in days",
      pval = TRUE,
      font.x = 16, 
      font.y=16,
      risk.table.height = 0.25,
      fontsize=3, 
      risk.table.col= "strata",
      size=1, palette = c("#2E9FDF", "red"),
      risk.table.fontsize = 5, 
      risk.table.y.text = FALSE
      
    )  
    print(gg)
    dev.off()
    
    print(paste("done for", gen[j]))
  }
  
  print(paste("done for Cancer", CanType[i]))
  write.csv(p, file=paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases/", CanType[i], "_Pvalue.csv", sep = ""))
  save(p, file=paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases/", CanType[i], "_Pvalue.RData", sep = "")) 
}
#write.csv(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/BLCA_Pvalue.csv")
#save(p, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/BLCA_Pvalue.RData)

###########################################################
# get all the Pvalues for kinases together in one csv file
########################################

file <- "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases"
Survival <- data_frame()
for(i in 2:length(CanType)){
  p <- read.csv(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases/", CanType[i], "_Pvalue.csv", sep = ""))
  
  p <- p %>% select(variable, pval.txt)
  p1 <- separate(data = p, col = variable, into = c("left", "right"), sep = "\\+")
  p1 <- separate(data = p1, col = pval.txt, into = c("left1", "right1"), sep = "\\=")
  which(p$pval.txt == "p < 0.0001")
  p2 <- p1[,c(1,4)]
  
  p2$right1[is.na(p2$right1)] <- "< 0.0001"
  
  
  colnames(p2) <- c("Genes", paste0("Survival_P_", CanType[i]))
  
  #p3 <- p2
  
  #Survival <- p2
  
  #p2$Row.names <- p2$Genes
  #Survival <- p2
  #kk <- Underst_adj.P
  Survival <- merge(as.data.frame(Survival), as.data.frame(p2), by='Genes', all=TRUE)
  
  #Survival <- cbind(Survival, p2[,2,drop = FALSE])
  
  print(head(Survival))
  
}

write.csv(Survival, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases/SurvivalPvalue_cancersWithoutContl.csv", row.names = FALSE)








################################################
#scoring for stages in each cancer T_M_N for all the kinases
################################################

#read the files

#IDG_kinase <- read.csv("/projects/ccs/schurerlab/Rimpi/kinases_IDG.txt", sep="\t")
Dir <- "/projects/ccs/schurerlab/DerekJ/TCGA/Firehose_Clinical/"
head(Kinases)

#no stages "CESC", "GBM" ,"PCPG", "PRAD","UCEC"

CanType <- c("BLCA",  "KIRP","KIRC", "BRCA", "COAD", "ESCA",  "HNSC" , "LIHC", "LUAD", "LUSC", "STAD", "THCA", "CESC",  "PRAD", "CHOL", "KICH" ,  "READ" ,"ACC","SKCM","PAAD")     #"BLCA", "KICH", "KIRP"

#does noit have pathlogical stage data
#"GBM" , "PCPG", "LGG","SARC","OV","UVM",

#clinical stages
#UCEC - clinical stage (does not have pathlogical stage)
#THYM
#OV
#UVM
#CESC



#TGCT different type of stage 1S


#Save files in 

save <- "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/stagewisePvalueForeachCancer/Firehous_clinical/"


################################################
#for loop to call all the files and and calculate the P-value for each stage using 
################################################

for(i in 19:length(CanType)){
Anova_clinical_stage_t <- data.frame()
Anova_pathology_M_stage <- data.frame()
Anova_pathology_N_stage <- data.frame()
Anova_pathology_T_stage <- data.frame()



print(CanType[i] )
a <- read.csv(paste(Dir, "TCGA_", CanType[i], "_fhclinicalmerge1.csv", sep=""))
#print(table(a$clinical_stage_t))
print(table(a$pathology_M_stage))
print(table(a$pathology_N_stage))
print(table(a$pathology_T_stage))
b <- a[, colnames(a) %in% Kinases]
bcr_patient_barcode <- data.frame(bcr_patient_barcode = paste(a$TCGA_ID, "-", a$Sample, sep =""))
c <- select(a, pathology_M_stage, pathology_N_stage, pathology_T_stage)
c$bcr_patient_barcode <- bcr_patient_barcode

m <- as.character(c$pathology_M_stage)
m[m == "cm0 (i+)" ] <- "m0"
m[m == "m1a" ] <- "m1"
m[m == "m1b" ] <- "m1"
m[m == "m1c" ] <- "m1"
c$pathology_M_stage <- m



n <- as.character(a$pathology_N_stage)
n[n == "n0 (i-)" ] <- "n0"
n[n == "n0 (i+)" ] <- "n0"
n[n == "n0 (mol+)" ] <- "n0"
n[n == "n1a" ] <- "n1"
n[n == "n1b" ] <- "n1"
n[n == "n1c" ] <- "n1"
n[n == "n1mi" ] <- "n1"
n[n == "n2a" ] <- "n2"
n[n == "n2c" ] <- "n2"
n[n == "n2c" ] <- "n2"
n[n == "n3a" ] <- "n3"
n[n == "n3b" ] <- "n3"
n[n == "n3c" ] <- "n3"
c$pathology_N_stage <- n
f <- as.character(c$pathology_T_stage)
f[f == "tis" ] <- "t0"
f[f == "t1a" ] <- "t1"
f[f == "t1a1" ] <- "t1"
f[f == "t1b" ] <- "t1"
f[f == "t1b1" ] <- "t1"
f[f == "t1b2" ] <- "t1"
f[f == "t1c" ] <- "t1"
f[f == "t2a1" ] <- "t2"
f[f == "t2a2" ] <- "t2"
f[f == "t2b" ] <- "t2"
f[f == "t2c" ] <- "t2"
f[f == "t3a" ] <- "t3"
f[f == "t3b" ] <- "t3"
f[f == "t3c" ] <- "t3"
f[f == "t4a" ] <- "t4"
f[f == "t4b" ] <- "t4"
f[f == "t4c" ] <- "t4"
f[f == "t4d" ] <- "t4"
f[f == "t4e" ] <- "t4"
c$pathology_T_stage <- f


# g <- as.character(a$clinical_stage_t)
# g[g == "t2a" ] <- "t2"
# g[g == "t2b" ] <- "t2"
# g[g == "t3b" ] <- "t3"
# g[g == "t3b" ] <- "t3"
# g[g == "t4a" ] <- "t4"
# c$clinical_stage_t <- g
c <- cbind(c, b)
#stage_t <- subset(c, subset = clinical_stage_t == "t1" | clinical_stage_t == "t2" | clinical_stage_t == "t3"| clinical_stage_t == "t4" )
p_M_stage <- subset(c, subset = pathology_M_stage == "m0" | pathology_M_stage == "m1" )
p_N_stage <- subset(c, subset = pathology_N_stage == "n0" | pathology_N_stage == "n1" | pathology_N_stage == "n2" | pathology_N_stage == "n3")
p_T_stage <- subset(c, subset = pathology_T_stage == "t0" | pathology_T_stage == "t1" | pathology_T_stage == "t2" | pathology_T_stage == "t3" | pathology_T_stage == "t4")

genes <- as.character(colnames(b))
# if(length(stage_t$clinical_stage_t != 0)){
# for(j in 6:ncol(stage_t)){
#   print(colnames(stage_t[j]))
#   formula <- as.formula(paste(colnames(stage_t)[j], " ~ clinical_stage_t", sep=""))
#   res.aov <- aov(formula, data = stage_t)
#   #res.aov <- aov(colnames(c[j]) ~ ajcc_pathologic_tumor_stage, data = c)
#   #res.aov_grade <- aov(genes[,j] ~ histological_grade, data = c)
#   print(paste("Anova for the gene == ", colnames(c[j])))

#   print(TukeyHSD(res.aov))
#   res.aov_result <-  TukeyHSD(res.aov)[1]
#   res.aov_result <- data.frame(res.aov_result$clinical_stage_t)
#   res.aov_result$gene <- colnames(stage_t[j])

#   Anova_s <- res.aov_result
#   Anova_clinical_stage_t <- rbind(Anova_clinical_stage_t, Anova_s)
#   print(paste("*************gene is done ", colnames(stage_t[j]), "for", CanType[i], "*****************"))
# }
#   }else{

#     print(paste("*************stage_t is not Avaiable:::::", CanType[i], "***************" , sep =""))
if(length(p_M_stage$pathology_M_stage != 0)){
for(k in 5:ncol(p_M_stage)){
formula1 <- as.formula(paste(colnames(p_M_stage)[k], " ~ pathology_M_stage", sep=""))
res.aov_grade1 <- aov(formula1, data = p_M_stage)
res.aov_grade_result1 <-  TukeyHSD(res.aov_grade1)[1]
res.aov_grade_result1 <-  data.frame(res.aov_grade_result1$pathology_M_stage)
res.aov_grade_result1$gene <- colnames(p_M_stage[k])
Anova_g <- res.aov_grade_result1
Anova_pathology_M_stage  <- rbind(Anova_pathology_M_stage , Anova_g)
}

}else{

print(paste("*************histological_grade is not Avaiable:::::", CanType[i], "***************" , sep =""))
}
if(length(p_N_stage$pathology_N_stage != 0)){
for(l in 5:ncol(p_N_stage)){
print(colnames(p_N_stage[l]))
formula2 <- as.formula(paste(colnames(p_N_stage)[l], " ~ pathology_N_stage", sep=""))
res.aov2 <- aov(formula2, data = p_N_stage)
print(paste("Anova for the gene == ", colnames(p_N_stage[l])))
print(TukeyHSD(res.aov2))
res.aov_result2 <-  TukeyHSD(res.aov2)[1]
res.aov_result2 <- data.frame(res.aov_result2$pathology_N_stage)
res.aov_result2$gene <- colnames(p_N_stage[l])
Anova_N <- res.aov_result2
Anova_pathology_N_stage <- rbind(Anova_pathology_N_stage, Anova_N)

print(paste("*************gene is done ", colnames(p_N_stage[l]), "for", CanType[i], "*****************"))
}
}else{

print(paste("*************pathological_N_stage is not Avaiable:::::", CanType[i], "***************" , sep =""))
}
if(length(p_T_stage$pathology_M_stage != 0)){

#if(e$histological_grade == 0){
for(p in 5:ncol(p_M_stage)){
#for the grade
formula3 <- as.formula(paste(colnames(p_T_stage)[p], " ~ pathology_T_stage", sep=""))
res.aov_grade3 <- aov(formula3, data = p_T_stage)
res.aov_grade_result3 <-  TukeyHSD(res.aov_grade3)[1]
res.aov_grade_result3$gene <- colnames(p_T_stage[p])
Anova_T <- res.aov_grade_result3

Anova_pathology_T_stage  <- rbind(Anova_pathology_T_stage , Anova_T)
}
print(paste("*************gene is done ", colnames(p_T_stage[p]), "for", CanType[i], "*****************"))
}else{
print(paste("*************pathology_T_stage is not Avaiable:::::", CanType[i], "***************" , sep =""))
}
print(paste("*************the canceris done::::::", CanType[i], "****************", sep ="" ))
# write.csv(Anova_clinical_stage_t, file = paste(save, CanType[i], "Anova_clinical_stage_t.csv"))
# length(unique(Anova_stage$gene)) == length(genes)
write.csv(Anova_pathology_M_stage, file = paste(save, CanType[i], "Anova_pathology_M_stage.csv"))length(unique(Anova_pathology_M_stage$gene)) == length(genes)
write.csv(Anova_pathology_N_stage, file = paste(save, CanType[i], "Anova_pathology_N_stage.csv"))

length(unique(Anova_pathology_N_stage$gene)) == length(genes)
write.csv(Anova_pathology_T_stage, file = paste(save, CanType[i], "Anova_pathology_T_stage.csv"))
length(unique(Anova_pathology_T_stage$gene)) == length(genes)
}





###################
# merge tables for survival analysis (without controls and with controls) kinases for the scoring 
#################

kin_sur <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases/FinalTableSurvivalwithPvalue.csv")
kin_sur_no_control <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases/SurvivalPvalue_cancersWithoutContl.csv")
Merge_sur <- merge(kin_sur, kin_sur_no_control, by = "Genes")

#change colnames of the data frame 
col <- gsub("\\Survival_P_", "", colnames(Merge_sur))
colnames(Merge_sur) <- col

#Remove the understudied kinases from kinases survival 
gene_Com <- intersect(Merge_sur$Genes,understudied)
sur_Kinases <- Merge_sur[!Merge_sur$Genes %in% gene_Com,]


rownames(sur_Kinases) <- sur_Kinases$Genes

sur_Kinases <- sur_Kinases[,-1]

#Sort the colnames in alphabetic order

Survival_kinases <- sur_Kinases %>% 
  select(sort(names(.)))

write.csv(Survival_kinases, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/SurvivalPvalue/Kinases/Survival_Without_understudied_kinases.csv")
############
#data for the app 
############
geneid <- as.character(rownames(Survival_kinases))
CanType <- as.character(colnames(Survival_kinases))
tt <- data.frame()
for(bb in 1:length(geneid)){
  for(cc in 1:length(CanType)){
    
    FileName <- paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[cc], "/survival/Kinase/",  geneid[bb],"_", CanType[cc],"_Pvalue.RData", sep="")
    
    if (file.exists(FileName)) {
    load(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[cc], "/survival/Kinase/",  geneid[bb],"_", CanType[cc],"_Pvalue.RData", sep=""))
    gg <- F
    
    colnames(gg) <- c("times", "patient.vital_status", "expr", "cohart")
    gg$gene <- paste(geneid[bb])
    gg$cohart <- paste(CanType[cc])
    
    tt <- rbind(tt, gg)
    } else {
      cat("file does not exist")
    }
  }
}

write.csv(tt, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/data/Survival_Kinases.csv")

#################
#
###############

gene_final <- read.delim("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/624.txt")


#################
#merge each cancer scored file with IDG
###############

library("xlsx")


#read scored IDG file 
IDG <- read.delim("/projects/ccs/schurerlab/DerekJ/TCGA/Scored/listofkinases.csv", sep =",")

IDG <- IDG[,c(2,4)]
colnames(IDG) <- c("Gene", "tdl")

#read scored file for each cancer 
CancerID <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", "READ", "SKCM", "STAD", "THCA")
for(i in 1:length(CancerID)){
  paste(print(CancerID[i]))
fil <- read.xlsx(file = paste("/projects/ccs/schurerlab/DerekJ/TCGA/Scored/", CancerID[i] , "_ScoredFinal.xlsx", sep =""), sheetIndex =1, header=TRUE, colClasses=NA)


fil1 <- merge(fil, IDG, by = "Gene", all.x = TRUE)


write.csv(fil1, file = paste("/projects/ccs/schurerlab/DerekJ/TCGA/Scored/final_score/", CancerID[i] , "_ScoredFinal.csv", sep =""))

}

#################
#merge each cancer scored file with survival and DE files and Mutation files
###############

Sur  <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/SCORESFinalTableSurvivalwithPvalue.xlsx",  sheetIndex =2, header=TRUE, colClasses=NA)
DE <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/SCOREKinases_all_Pvalue_FC.xlsx",  sheetIndex =2, header=TRUE, colClasses=NA)
Mut <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/Mutation_New_PvalueSCORED.xlsx", sheetIndex =2, header=TRUE, colClasses=NA)


#####Change colnames of the files
#Survival
col <- gsub("\\Score_", "", colnames(Sur))
colnames(Sur) <- col

#DE
col <- gsub("\\Score_", "", colnames(DE))
colnames(DE) <- col

#Mut
col <- gsub("\\Score_", "", colnames(Mut))
colnames(Mut) <- col
##########
#MergeSurvival, mutation and DE 



Can <- as.character(intersect(intersect(colnames(Sur), colnames(DE)), colnames(Mut)))

for(i in 1:length(Can)){
  print(Can[i])
  
  sub <- subset(Mut, select= c("Gene", paste(Can[i])))
  colnames(sub) <- c("Gene", "Mut")
  sub1 <- subset(Sur, select= c("Genes", paste(Can[i])))
  colnames(sub1) <- c("Gene", "Survival")
  sub2 <- subset(DE, select= c("Row.names", paste(Can[i])))
  colnames(sub2) <- c("Gene", "DE")
  mer <- merge(sub, sub1, by = "Gene", all = TRUE)
  mer1 <- merge(mer, sub2, by = "Gene", all = TRUE)
  #colnames(mer) <- c("Gene", "Mutation", "Survival", "DE")
  write.csv(mer1, paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/F_De_Sur_Mut/merging/", Can[i], "_DE_sur_Mut.csv", sep =""), row.names = FALSE)
}


#LIHC has values in mut and Sur so merge those 

Can <- "LIHC"

sub1 <- subset(Sur, select= c("Genes", paste(Can[i])))
colnames(sub1) <- c("Gene", "Survival")
sub <- subset(Mut, select= c("Gene", paste(Can[i])))
colnames(sub) <- c("Gene", "Mut")
mer <- merge(sub, sub1, by = "Gene", all = TRUE)
write.csv(mer, paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/F_De_Sur_Mut/merging/", Can[i], "_DE_sur_Mut.csv", sep =""), row.names = FALSE)

#############
#merge clinical with Survival, mutaion and DE analysis
###########


Can <- as.character(intersect(intersect(colnames(Sur), colnames(DE)), colnames(Mut)))  #PCPG, UCEC, GBM, COAD

for(i in 1:length(Can)){
  print(Can[i])
fil  <- read.xlsx(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/", Can[i], "_ScoredFinal.xlsx", sep =""),  sheetIndex =1, header=TRUE, colClasses=NA)
fil2 <- read.delim(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/F_De_Sur_Mut/merging/", Can[i], "_DE_sur_Mut.csv", sep = ""), sep = ",")

All_Clin_smd <- merge(as.data.frame(fil), as.data.frame(fil2), by = "Gene", all=TRUE)  #merge all clinic, survival mutation and DE

write.csv(All_Clin_smd, paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/F_De_Sur_Mut/", Can[i], "_ScoredFinal_clin_DE_sur_Mut.csv", sep =""), row.names = FALSE)

}


#LIHC
Can <- "LIHC"
for(i in 1:length(Can)){
  print(Can[i])
  fil  <- read.delim(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/", Can[i], "_ScoredFinal.csv", sep =""),  sep = ",")
  fil2 <- read.delim(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/F_De_Sur_Mut/merging/", Can[i], "_DE_sur_Mut.csv", sep = ""), sep = ",")
  
  All_Clin_smd <- merge(as.data.frame(fil), as.data.frame(fil2), by = "Gene", all=TRUE)  #merge all clinic, survival mutation and DE
  
  write.csv(All_Clin_smd, paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/F_De_Sur_Mut/", Can[i], "_ScoredFinal_clin_DE_sur_Mut.csv", sep =""), row.names = FALSE)
  
}









######################
#DEAnalysis and Survival heatmaps 
######################

DE <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/SCOREKinases_all_Pvalue_FC.xlsx",  sheetIndex =2, header=TRUE, colClasses=NA)
col <- gsub("\\Score_", "", colnames(DE))
colnames(DE) <- col

DE_Ana <- DE
#Remove columns with all 0  Rows with all 0 

rownames(DE_Ana) <- DE_Ana[,1]
DE_Ana <- DE_Ana[,-c(1,22)]

colSums(DE_Ana) #624

#Final_Sur_kinase <- DE_Ana[,-which(colSums(DE_Ana) == 0) ]  #colSums== 0 


Final_DE_kinase <- DE_Ana[rowSums(DE_Ana) >1,]  #rowSums> 1  #425

Final_DE_kinase <- Final_DE_kinase[order(rowSums(Final_DE_kinase), decreasing = T),]

Final_DE_kinase <- Final_DE_kinase[,order(colSums(Final_DE_kinase))]
data_matrix <- data.matrix(Final_DE_kinase)   #394 genes 19 cancers


#heatmap(t(Final_Sur_kinase2), scale = "none", Rowv = NA, Colv = NA, col= c("grey", "black"))
#
library("superheat")

png("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/FiguresForPaper/DE_kinases.png",width = 1000, height = 800, units = "px")
superheat(X = t(data_matrix), # heatmap matrix
          title = "DE Kinases Analysis",
          # change the angle of the label text
          bottom.label.text.angle = 90,left.label.text.size = 5,
          left.label.text.alignment = "right",
          bottom.label.text.alignment = "right",
          # scale the matrix columns
          heat.pal = c("lightgrey", "red"),
          bottom.label.text.size = 2, force.bottom.label = T, bottom.label.size = 0.2, yr = rowSums(t(data_matrix)) ,
          yr.axis.name = "Total Genes in Cancer",yr.plot.type = "bar", yr.bar.col = "black", yr.obs.col = rep("red", nrow(t(data_matrix))), 
          yt = colSums(t(data_matrix)), yt.axis.name = "Genes in all Cancer", yt.obs.col = rep("red", ncol(t(data_matrix))))
dev.off()


######################
#Survival Analysis heapmaps for low vs high
######################
#Sur  <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/SCORESFinalTableSurvivalwithPvalue.xlsx",  sheetIndex =2, header=TRUE, colClasses=NA)
Sur <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/Survival_P_value_Low.xlsx", sheetIndex =1, header=TRUE, colClasses=NA)
Sur$variable <- str_remove_all(Sur$variable, "[+cohort]")
#col <- gsub("\\Score_", "", colnames(Sur))
#colnames(Sur) <- col


#Remove columns with all 0  Rows with all 0 
Sur_Ana <- Sur[,c(1,7,8)]
Sur_Ana$survival <- as.numeric(levels(Sur_Ana$survival))[Sur_Ana$survival]

Sur_Ana$survival[is.na(Sur_Ana$survival)] <- 0

Sur_Ana1 <-  dcast(Sur_Ana, variable ~ Cancer)
Sur_Ana1$survival[is.na(Sur_Ana1$survival)] <- 0


rownames(Sur_Ana1) <- Sur_Ana1[,1]
Sur_Ana1 <- Sur_Ana1[,-1]

colSums(Sur_Ana1)
rowSums(Sur_Ana1)

#Final_Sur_kinase <- DE_Ana[,-which(colSums(DE_Ana) == 0) ]  #colSums== 0 


Final_Sur_kinase <- Sur_Ana1[rowSums(Sur_Ana1) >1,]  #rowSums> 1

Final_Sur_kinase <- Final_Sur_kinase[order(rowSums(Final_Sur_kinase), decreasing = T),]

Final_Sur_kinase <- Final_Sur_kinase[,-which(colSums(Final_Sur_kinase) == 0) ]

Final_Sur_kinase <- Final_Sur_kinase[,order(colSums(Final_Sur_kinase))]
data_matrix <- data.matrix(Final_Sur_kinase)   #624 genes 19 cancers


#heatmap(t(Final_Sur_kinase2), scale = "none", Rowv = NA, Colv = NA, col= c("grey", "black"))
#
library(superheat)
#png("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Paper/Figures_paper/Survival_kinases.png",width = 1000, height = 800, units = "px")
png("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Paper/Figures_paper/Survival_kinases_highvslow.png",width = 1000, height = 800, units = "px")
superheat(X = t(data_matrix), # heatmap matrix
          title = "Survival Kinases Analysis",
          # change the angle of the label text
          bottom.label.text.angle = 90,left.label.text.size = 5,
          left.label.text.alignment = "right",
          bottom.label.text.alignment = "right",
          # scale the matrix columns
          heat.pal = c("lightgrey", "red"),
          bottom.label.text.size = 2, force.bottom.label = T, bottom.label.size = 0.2, yr = rowSums(t(data_matrix)) ,
          yr.axis.name = "Total Genes in Cancer",yr.plot.type = "bar", yr.bar.col = "black", yr.obs.col = rep("red", nrow(t(data_matrix))), 
          yt = colSums(t(data_matrix)), yt.axis.name = "Genes in all Cancer", yt.obs.col = rep("red", ncol(t(data_matrix))))
dev.off()



######################
#Mutation Analysis heapmaps 
######################
Mut <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/Mutation_New_PvalueSCORED.xlsx", sheetIndex =2, header=TRUE, colClasses=NA)

col <- gsub("\\Score_", "", colnames(Mut))
colnames(Mut) <- col

Mut_Ana <- Mut
#Remove columns with all 0  Rows with all 0 

rownames(Mut_Ana) <- Mut_Ana[,1]
Mut_Ana <- Mut_Ana[,-c(1,22)]

colSums(Mut_Ana)
rowSums(Mut_Ana)

Final_Mut_kinase <- Mut_Ana[,-which(colSums(Mut_Ana) == 0) ]  #colSums== 0 


Final_Mut_kinase <- Final_Mut_kinase[rowSums(Final_Mut_kinase) >1,]  #rowSums> 1

Final_Mut_kinase <- Final_Mut_kinase[order(rowSums(Final_Mut_kinase), decreasing = T),]

Final_Mut_kinase <- Final_Mut_kinase[,order(colSums(Final_Mut_kinase))]
data_matrix <- data.matrix(Final_Mut_kinase)   #394 genes 19 cancers


#heatmap(t(Final_Mut_kinase2), scale = "none", Rowv = NA, Colv = NA, col= c("grey", "black"))
#
library(Superheatmap)

png("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/FiguresForPaper/Mutation_kinases.png",width = 1000, height = 800, units = "px")
superheat(X = t(data_matrix), # heatmap matrix
          title = "Muttation Kinases Analysis",
          # change the angle of the label text
          bottom.label.text.angle = 90,left.label.text.size = 5,
          left.label.text.alignment = "right",
          bottom.label.text.alignment = "right",
          # scale the matrix columns
          heat.pal = c("lightgrey", "red"),
          bottom.label.text.size = 4, force.bottom.label = T, bottom.label.size = 0.3, yr = rowSums(t(data_matrix)) ,
          yr.axis.name = "Total Genes in Cancer",yr.plot.type = "bar", yr.bar.col = "black", yr.obs.col = rep("red", nrow(t(data_matrix))), 
          yt = colSums(t(data_matrix)), yt.axis.name = "Genes in all Cancer", yt.obs.col = rep("red", ncol(t(data_matrix))))
dev.off()


######################
#DEAnalysis updated IDG list heatmaps 
######################
updated_IDG_kinases <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/Understudiedkinases_2019June.csv")
DE <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/SCOREKinases_all_Pvalue_FC.xlsx",  sheetIndex =2, header=TRUE, colClasses=NA)
col <- gsub("\\Score_", "", colnames(DE))
colnames(DE) <- col

DE_Ana <- DE
#Remove columns with all 0  Rows with all 0 

rownames(DE_Ana) <- DE_Ana[,1]
DE_Ana <- DE_Ana[,-c(1,22)]

colSums(DE_Ana) #624

#Final_Sur_kinase <- DE_Ana[,-which(colSums(DE_Ana) == 0) ]  #colSums== 0 


Final_DE_kinase <- DE_Ana[rowSums(DE_Ana) >1,]  #rowSums> 1  #425

Final_DE_kinase <- Final_DE_kinase[order(rowSums(Final_DE_kinase), decreasing = T),]

Final_DE_kinase <- Final_DE_kinase[,order(colSums(Final_DE_kinase))]
data_matrix <- data.matrix(Final_DE_kinase)   #394 genes 19 cancers

saveRDS(Final_DE_kinase, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/data/Upregulated_Kinases_cancer.rds")  

#heatmap(t(Final_Sur_kinase2), scale = "none", Rowv = NA, Colv = NA, col= c("grey", "black"))
#
library("superheat")

png("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/FiguresForPaper/DE_kinases.png",width = 1000, height = 800, units = "px")
superheat(X = t(data_matrix), # heatmap matrix
          title = "DE Kinases Analysis",
          # change the angle of the label text
          bottom.label.text.angle = 90,left.label.text.size = 5,
          left.label.text.alignment = "right",
          bottom.label.text.alignment = "right",
          # scale the matrix columns
          heat.pal = c("lightgrey", "red"),
          bottom.label.text.size = 2, force.bottom.label = T, bottom.label.size = 0.2, yr = rowSums(t(data_matrix)) ,
          yr.axis.name = "Total Genes in Cancer",yr.plot.type = "bar", yr.bar.col = "black", yr.obs.col = rep("red", nrow(t(data_matrix))), 
          yt = colSums(t(data_matrix)), yt.axis.name = "Genes in all Cancer", yt.obs.col = rep("red", ncol(t(data_matrix))))
dev.off()






###########################
#upload new IDG kinase list and DE genes of all cancer and extract understudied kinases
##########################
# differentially expressed Kinases 
DE <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FINALSCORE/SCOREKinases_all_Pvalue_FC.xlsx",  sheetIndex =2, header=TRUE, colClasses=NA)
colnames(DE)[1] <- "Gene"
# updated IDG list 
updated_IDG_kinases <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/Understudiedkinases_2019June.csv")


Intersect_DE <- intersect(DE$Row.names, updated_IDG_kinases$Gene)


#extract the score for understudied kinases


unstudiedK <- join(DE, updated_IDG_kinases, by="Gene", type="inner")

#select only those understudied kinases present in atleast 2 cancers

unstudiedK1 <- unstudiedK[which(unstudiedK$Total_Score > 1),]  # 102
unstudiedK1 <- unstudiedK1[,1:21]

rownames(unstudiedK1) <- unstudiedK1[,1]

unstudiedK1 <- unstudiedK1[,-1]

#edit the columnames remove score

col <- gsub("\\Score_", "", colnames(unstudiedK1))

colnames(unstudiedK1) <- col

#generate a superheatmap for this 


unstud_DE <-unstudiedK1[order(rowSums(unstudiedK1), decreasing = T),]

unstud_DE <- unstud_DE[,order(colSums(unstud_DE))]

# save for running in the app
saveRDS(unstud_DE, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/data/Upregulated_understudiedKinases_cancer.rds")  
data_matrix <- data.matrix(unstud_DE)   #102 in 20 cancers


#heatmap(t(Final_Sur_kinase2), scale = "none", Rowv = NA, Colv = NA, col= c("grey", "black"))
#
library("superheat")

png("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/FiguresForPaper/understudied_DE_kinases.png",width = 1000, height = 800, units = "px")
superheat(X = t(data_matrix), # heatmap matrix
          title = "DE Kinases Analysis",
          # change the angle of the label text
          bottom.label.text.angle = 90,left.label.text.size = 5,
          left.label.text.alignment = "right",
          bottom.label.text.alignment = "right",
          # scale the matrix columns
          heat.pal = c("lightgrey", "red"),
          bottom.label.text.size = 2, force.bottom.label = T, bottom.label.size = 0.2, yr = rowSums(t(data_matrix)) ,
          yr.axis.name = "Total Genes in Cancer",yr.plot.type = "bar", yr.bar.col = "black", yr.obs.col = rep("red", nrow(t(data_matrix))), 
          yt = colSums(t(data_matrix)), yt.axis.name = "Genes in all Cancer", yt.obs.col = rep("red", ncol(t(data_matrix))))
dev.off()




###########################
#score for publication and app (rank the kinases by column "kinase score percent")
##########################

library(dplyr)
library(xlsx)



#files OV ACC are in csv so read them as csv
CanType <- c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "READ", "STAD", "THCA")

Score_percentage <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/BLCA_ScoredFinal_clin_DE_sur_Mut.xlsx", sheetIndex=1, header=TRUE, colClasses=NA)

Score_percentage <- select(Score_percentage, Gene, class)
colnames(Score_percentage) <- c("Genes", "class")
score <- data.frame()
for(i in 1:length(CanType)){
print(CanType[i])

scor <- paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/", CanType[i], "_ScoredFinal_clin_DE_sur_Mut.xlsx", sep="")
if(file.exists(scor)) {
  scor <- read.xlsx(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/", CanType[i], "_ScoredFinal_clin_DE_sur_Mut.xlsx", sep=""), sheetIndex=1, header=TRUE, colClasses=NA)
} else {
  scor <- read.csv(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/", CanType[i], "_ScoredFinal_clin_DE_sur_Mut.csv", sep=""))
}

scor <- scor %>% mutate(rank = dense_rank(desc(Kinase_score_percent)))
scor_fil <- select(scor,  Gene, Cancer, Kinase_score_percent, rank, class)

scor_filt <- select(scor,  Gene, Kinase_score_percent, class)
colnames(scor_filt) <- c("Genes", paste0("Kinase_score_percent_", CanType[i]), "class")

score <- rbind(score, scor_fil)

Score_percentage <- merge(Score_percentage, scor_filt , by = c("Genes", "class"), all=TRUE)
} 


write.csv(Score_percentage, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/Score_percentage.csv", row.names = FALSE)
write.csv(score, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/rank_for_App.csv", row.names = FALSE)


################
#merge the full kinase infirmation with rank 
##############
full_kinase <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/SupplementalInfo/Kinase_Full_Data_Final.csv")

col <- gsub("\\Score_", "", colnames(Mut))
col <- gsub("\\Kinase_score_percent_","",  colnames(full_kinase))

colnames(full_kinase) <- col
#remove kinase score average column

full_kinase <- full_kinase[,-20]
ful_mel <- melt(full_kinase)

names(ful_mel)[11] <- paste("Cancer")


###get the rank file and merge it with kinase full information file
rank <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/rank_for_App.csv")

data<- ful_mel %>% right_join(rank, by = c("Gene", "Cancer"))
data <- data[,-12]




write.csv(data, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/data/Kinase_full_information_for_App.csv", row.names = FALSE)
#merge moA data with data 
rank_data <- data.frame(read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/data/Kinase_full_information_for_App.csv"))
rank_data <- rank_data[,c(1,10,2,3,11,12, 5:8)]
new_dataset <- rank_data %>% right_join(TargetMOA, by=c("Gene", "Cancer"))

new_dataset <- new_dataset[,1:11]
colnames(new_dataset) <- c("Gene", "Cancer", "TDL",  "IDG_Status",  "Kinase_score_percent", "rank",     "Kinase_Type", "Kinase_Group", "Kinase_Family" ,"Kinase_protein","MOA")
levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_BLCA"] <- "yes"
levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_BRCA"] <- "yes"
levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_CHOL"] <- "yes"
> levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_HNSC"] <- "yes"
> levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_KICH"] <- "yes"
> levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_KIRC"] <- "yes"
> levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_KIRP"] <- "yes"
> levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_KIRP"] <- "yes"
> levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_LIHC"] <- "yes"
> levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_LUAD"] <- "yes"
> levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_LUSC"] <- "yes"
> levels(new_dataset$MOA)[levels(new_dataset$MOA) == "TMOA_THCA"] <- "yes"
levels(new_dataset$MOA)[levels(new_dataset$MOA) == "Tbio"] <- "no"
levels(new_dataset$MOA)[levels(new_dataset$MOA) == "Tchem"] <- "no"
levels(new_dataset$MOA)[levels(new_dataset$MOA) == "Tclin"] <- "no"
levels(new_dataset$MOA)[levels(new_dataset$MOA) == "Tdark"] <- "no"

write.csv(new_dataset, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/data/Kinase_full_information_for_App.csv", row.names = FALSE)
###################
#complete suplementary table /projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/SupplementalInfo/TCGA_clinicaldataoverview.xlsx
###################

CanType <- c("THCA", "CESC", "BLCA", "BRCA", "CHOL", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "READ", "STAD", "PRAD",  
             "PCPG", "UCEC", "PAAD","ACC","SKCM","GBM","OV", "SARC","LGG","TGCT","THYM","UCS","UVM","MESO" )

for(i in 1:length(CanType)){
  
  data1 <- readRDS(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[i], "/TCGA-", CanType[i], ".rds", sep = ""))
  
  print(paste("***********", CanType[i], "***********",   table(colData(data1)$Shortnametype)))
  
}


#####################
#cancatenate  all the FC and Pvalue data from DE analysis for the volcano plots 
####################


CanType <- c("BLCA",  "KIRP","KIRC", "BRCA", "COAD", "ESCA",  "HNSC" , "LIHC", "LUAD", "LUSC", "STAD", "THCA", "CESC", "GBM" , "PRAD", "CHOL", "KICH", "UCEC" ,  "READ" ,"PCPG")   #"CHOL", "KICH", "UCEC" ,  "READ" ,"PCPG" #


DEData <- data.frame()
for(bb in 1:length(CanType) ){
  print(CanType[bb])
  #load the DE genes with p-value from each cancer 
  load(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[bb], "/AllGeneswithFC.RData", sep =""))  #dataDEGs
  
  DE_d <- dataDEGs[,c(1,2,5)]
  DE_d$Genes <- rownames(DE_d)
  DE_d$Cancer <- CanType[bb]
  DE_d <- DE_d[, c(4,1:2,3,5)]
  rownames(DE_d) <- NULL
  
  DEData <- rbind(DEData, DE_d)
}

  write.csv(DEData, file= "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/All_cancers_genes_Pvalue_FC.csv", row.names = FALSE)
  write.csv(DEData, file= "/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/data/All_cancers_genes_Pvalue_FC.csv", row.names = FALSE)
  
  
  #####################
  #cancatenate  all the FC and Pvalue data from DE analysis for the volcano plots 
  ####################

############################################
PKMYT1 = 2.015556, FDR= 5.65091e-11

ITPKA= 2.009993, FDR = 0.0004100969

load(paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-", CanType[1], "/AllGeneswithFC.RData", sep =""))



ggplot(dd) + geom_point(aes(x=dd$logFC, y=-log10(dd$FDR), colour=threshold)) + 
  geom_text_repel(aes(x= dd$logFC, y = -log10(dd$FDR)))
  ggtitle("overexpression") + 
  theme(legend.position = "none", plot.title = element_text(size = rel(1.5), hjust = 0.5), axis.title = element_text(size= rel(1.25)))
  
  ###################
  # prepare staging files for app boxplots
  ##################
  
  #Read the files
  CanType <- c("BLCA",  "KIRP","KIRC", "BRCA", "COAD", "ESCA",  "HNSC" , "LIHC", "LUAD", "LUSC", "STAD", "THCA", "CESC",  "PRAD", "CHOL", "KICH" ,  "READ" ,"ACC","SKCM","PAAD")     #"BLCA", "KICH", "KIRP"
  
  #KIRCfh <- data.frame(read.csv("/projects/ccs/schurerlab/DerekJ/TCGA/Firehose_Clinical/TCGA_KIRC_fhclinicalmerge1.csv"))
  
  #ff <- KIRCfh %>% select(HIPK1, pathology_T_stage)
  
  #boxplot(HIPK1~pathology_T_stage, data = ff)
  #read the files
  
  #IDG_kinase <- read.csv("/projects/ccs/schurerlab/Rimpi/kinases_IDG.txt", sep="\t")
  Dir <- "/projects/ccs/schurerlab/DerekJ/TCGA/Firehose_Clinical/"

  
  #no stages "CESC", "GBM" ,"PCPG", "PRAD","UCEC"
  
  CanType <- c("BLCA",  "KIRP","KIRC", "BRCA", "COAD", "ESCA",  "HNSC" , "LIHC", "LUAD", "LUSC", "STAD", "THCA", "CESC",  "PRAD", "CHOL", "KICH" ,  "READ" ,"ACC","SKCM","PAAD")     #"BLCA", "KICH", "KIRP"
  
  #does noit have pathlogical stage data
  #"GBM" , "PCPG", "LGG","SARC","OV","UVM",
  
  #clinical stages
  #UCEC - clinical stage (does not have pathlogical stage)
  #THYM
  #OV
  #UVM
  #CESC
  
  
  
  #Save files in 
  
  save <- "/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/data/stage_expression/"
  
  
  ################################################
  #for loop to call all the files and and calculate the P-value for each stage using 
  ################################################
  
  for(i in 1:length(CanType)){
    
    a <- read.csv(paste(Dir, "TCGA_", CanType[i], "_fhclinicalmerge1.csv", sep=""))
    
    
    print(CanType[i]) 
          print(dim(a)) 
    
    
    #print(table(a$clinical_stage_t))
    print(table(a$pathology_M_stage))
    print(table(a$pathology_N_stage))
    print(table(a$pathology_T_stage))
    
    
    m <- as.character(a$pathology_M_stage)
    m[m == "cm0 (i+)" ] <- "m0"
    m[m == "m1a" ] <- "m1"
    m[m == "m1b" ] <- "m1"
    m[m == "m1c" ] <- "m1"
    a$pathology_M_stage <- m
    
    
    
    n <- as.character(a$pathology_N_stage)
    n[n == "n0 (i-)" ] <- "n0"
    n[n == "n0 (i+)" ] <- "n0"
    n[n == "n0 (mol+)" ] <- "n0"
    n[n == "n1a" ] <- "n1"
    n[n == "n1b" ] <- "n1"
    n[n == "n1c" ] <- "n1"
    n[n == "n1mi" ] <- "n1"
    n[n == "n2a" ] <- "n2"
    n[n == "n2c" ] <- "n2"
    n[n == "n2c" ] <- "n2"
    n[n == "n3a" ] <- "n3"
    n[n == "n3b" ] <- "n3"
    n[n == "n3c" ] <- "n3"
    a$pathology_N_stage <- n
    
    f <- as.character(a$pathology_T_stage)
    f[f == "tis" ] <- "t0"
    f[f == "t1a" ] <- "t1"
    f[f == "t1a1" ] <- "t1"
    f[f == "t1b" ] <- "t1"
    f[f == "t1b1" ] <- "t1"
    f[f == "t1b2" ] <- "t1"
    f[f == "t1c" ] <- "t1"
    f[f == "t2a" ] <- "t2"
    f[f == "t2a1" ] <- "t2"
    f[f == "t2a2" ] <- "t2"
    f[f == "t2b" ] <- "t2"
    f[f == "t2c" ] <- "t2"
    f[f == "t3a" ] <- "t3"
    f[f == "t3b" ] <- "t3"
    f[f == "t3c" ] <- "t3"
    f[f == "t4a" ] <- "t4"
    f[f == "t4b" ] <- "t4"
    f[f == "t4c" ] <- "t4"
    f[f == "t4d" ] <- "t4"
    f[f == "t4e" ] <- "t4"
    a$pathology_T_stage <- f
    
    

    
    
       
      
      print(paste("*************pathological_stage is not Avaiable:::::", CanType[i], "***************" , sep =""))
      
  write.csv(a, file = paste(save, CanType[i], "_stage.csv", sep =""))
  
  
  nam <- paste("stage", CanType[i], sep = ".")
  assign(nam, a)
  
  }
  
  save(stage.ACC,stage.BLCA, stage.BRCA, stage.CESC, stage.CHOL, stage.COAD, stage.ESCA, stage.HNSC, stage.KICH, stage.KIRC, 
       stage.KIRP, stage.LIHC, stage.LUAD, stage.LUSC, stage.PAAD, 
       file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/data/stage.RData")
    
  
  ################################################
  #figure5a ---/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/FiguresForPaper/Final_Figures_08_02_2019/Figure5/
  #for loop to call all the MSIGDB files and merge to make a corplot 
  ################################################
  CanType <- c("BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "HNSC" , "KICH" , "KIRC",  "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "READ",  "STAD", "THCA") 
  
  
  library(corrplot)
  top25_prognostic <- data.frame()
  for(i in 1:length(CanType)) {
    t <- print(i)
    print(CanType[i])
    #allprognostic <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/MSIGDB_pancancer_top25prognosticgenes.xlsx", sheetIndex=t, header=TRUE, colClasses=NA)
  allprognostic <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/MSIGDB_pancancer_allprognosticgenes.xlsx", sheetIndex=t, header=TRUE, colClasses=NA)
  allprognostic <- data.frame(allprognostic[7:106,1])
  allprognostic$cancer <- CanType[i]
  colnames(allprognostic) <- c("Gene_Set_Name", "Cancer")
  
  top25_prognostic <- rbind(top25_prognostic, allprognostic)
  
  } 
  
  #write.csv(c, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/MSIGDB_pancancer_allprognosticgenes.csv", row.names = FALSE)
  write.csv(top25_prognostic, "/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/MSIGDB_pancancer_top25prognosticgenes.csv", row.names = FALSE)
  
  library(colorRamp)
  library(grDevices)
  library(corrplot)
  
  prognostic <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/MSIGDB_pancancer_prognosticgenes.csv")
  don <- data.frame()
  for(i in 1:length(CanType)){
    for(j in 1:length(CanType)){
    test1 <- prognostic[which(prognostic$Cancer == CanType[i]), ]
    test2 <- prognostic[which(prognostic$Cancer == CanType[j]), ]
    len <- length(intersect(test1$Gene_Set_Name, test2$Gene_Set_Name))
    print(paste(CanType[i], "and", CanType[j] ))
    
    
    gen <- data.frame(cbind(CanType[i], CanType[j], len))
    colnames(gen) <- c("can1",  "can2",  "num")
    don <- data.frame(rbind(don,gen))
    
  }
  }  
  
  dcast_don <- dcast(don, formula = can1~can2)
  
  row.names(dcast_don) <- as.character(dcast_don$can1)
  
  dcast_don <- dcast_don[,-1]
  cols <- colnames(dcast_don)
  colnames(dcast_don) <- cols
  Matrix1 <- apply(dcast_don,2,as.numeric)
  
  
  row.names(Matrix1) <- row.names(dcast_don)
  colnames(Matrix1) <- colnames(dcast_don)

  col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                             "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                             "#4393C3", "#2166AC", "#053061"))
  
  
png("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/FiguresForPaper/Final_Figures_08_02_2019/Figure5/Overlap_gene_sets_cancers.jpg", width = 500, height = 600, units = "px")
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
#addCoef.col = rgb(0,0,0, alpha =0.6)
corrplot(Matrix1,  method = "color", type = "lower",  hclust.method="complete", is.corr = FALSE, order = "hclust",  tl.col = "black", tl.cex = 0.8,  col=col2(200), 
         addCoef.col = "black", mar=c(0,0,1,0), diag= TRUE, addgrid.col	="white", cl.lim = c(30, 100), number.cex=0.7, rect.lwd		=10)

#corrplot(Cor1, method = "color",type="lower",hclust.method="complete",order="hclust",col=col2(200),
         #diag = TRUE,addgrid.col	="white",rect.lwd		=10,tl.col = "black",addCoef.col = "black",number.cex=0.5)
  
dev.off()
  

################################################
#figure5b ---/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/FiguresForPaper/Final_Figures_08_02_2019/Figure5/
#spearman cor plot
################################################


library(xlsx)
Table <- read.xlsx("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/SupplementalInfo/Supplemental_Table_4_FullKinaseDataandScore.xlsx", sheetIndex =2, header=TRUE, colClasses=NA)
Table2 <- Table[,1:19]
Table2 <- Table2[,-2]


library(corrplot)
library(colorRamps)
library(grDevices)
row.names(Table2) <- as.character(Table2$Gene)
Table2 <- Table2[,-1]
cols <- colnames(Table2)
cols2 <- gsub("Kinase_score_","",cols)
colnames(Table2) <- cols2
Matrix1 <- apply(Table2,2,as.numeric)
row.names(Matrix1) <- row.names(Table2)
colnames(Matrix1) <- colnames(Table2)

Cor1 <- cor(Matrix1,method="spearman")
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
corrplot(Cor1, method = "color",type="lower",hclust.method="complete",order="hclust",col=col2(40))
corrplot(Cor1,type="lower",hclust.method="complete",order="hclust",col=col2(40),addgrid.col	="white")
png("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/FiguresForPaper/Final_Figures_08_02_2019/Figure5/figure5b.jpg", width = 550, height = 650, units = "px")

corrplot(Cor1, method = "color",type="lower",hclust.method="complete",order="hclust",col=col2(200),
         diag = TRUE,addgrid.col	="white",rect.lwd		=10,tl.col = "black",addCoef.col = "black",number.cex=0.8)

dev.off()


################################################
#TNM staging superheatmap  
################################################

TNM_stage <- read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/SupplementalInfo/Kinases_TNMClinicalScore_total.csv")
TNM_stage_sub <- TNM_stage[,c(1,2,3:19)]




col <- gsub("\\TNM_", "", colnames(TNM_stage_sub))
colnames(TNM_stage_sub) <- col

TNM_stage_sub_Ana <- TNM_stage_sub
#Remove columns with all 0  Rows with all 0 

rownames(TNM_stage_sub_Ana) <- TNM_stage_sub_Ana[,1]
TNM_stage_sub_Ana <- TNM_stage_sub_Ana[,c(3:19)]

colSums(TNM_stage_sub_Ana) #624

#Final_Sur_kinase <- DE_Ana[,-which(colSums(DE_Ana) == 0) ]  #colSums== 0 


#Final_DE_kinase <- DE_Ana[rowSums(DE_Ana) >1,]  #rowSums> 1  #425
TNM_stage_sub_Ana_kinase <- TNM_stage_sub_Ana[rowSums(TNM_stage_sub_Ana) >1,]

TNM_stage_sub_Ana_kinase <- TNM_stage_sub_Ana_kinase[order(rowSums(TNM_stage_sub_Ana_kinase), decreasing = T),]

TNM_stage_sub_Ana_kinase <- TNM_stage_sub_Ana_kinase[,order(colSums(TNM_stage_sub_Ana_kinase))]
data_matrix <- data.matrix(TNM_stage_sub_Ana_kinase)   #394 genes 19 cancers


#heatmap(t(Final_Sur_kinase2), scale = "none", Rowv = NA, Colv = NA, col= c("grey", "black"))
#
library("superheat")

png("/projects/ccs/schurerlab/Rimpi/TCGARecount2/Scoring/FinalScore_28-05-2019/ReadyforPublication/Analysis/FiguresForPaper/DE_kinases.png",width = 1000, height = 800, units = "px")
superheat(X = t(data_matrix), # heatmap matrix
          title = "DE Kinases Analysis",
          # change the angle of the label text
          bottom.label.text.angle = 90,left.label.text.size = 5,
          left.label.text.alignment = "right",
          bottom.label.text.alignment = "right",
          # scale the matrix columns
          heat.pal = c("lightgrey","darkgrey", "blue", "green"),
          bottom.label.text.size = 2, force.bottom.label = T, bottom.label.size = 0.2, yr = rowSums(t(data_matrix)) ,
          yr.axis.name = "Total score of Genes in Cancer",yr.plot.type = "bar", yr.bar.col = "black", yr.obs.col = rep("red", nrow(t(data_matrix))),
          yt = colSums(t(data_matrix)), yt.axis.name = "Genes in all Cancer", yt.obs.col = rep("red", ncol(t(data_matrix))))
          
dev.off()

#yt = colSums(t(data_matrix)), yt.axis.name = "Genes in all Cancer", yt.obs.col = rep("red", ncol(t(data_matrix))))

## get the counts
TNM_stage_stack <- TNM_stage[,c(1,22,42,60)]

TNM_stage_subset <- TNM_stage[TNM_stage$Gene  %in% rownames(data_matrix),]
TNM_stage_subset1 <- TNM_stage_subset[,c(1,4:20,24:40, 43:59)]
stagecan <- colnames(TNM_stage_subset1[,-1])
tda1 <- data.frame()


for(i in 1:length(stagecan)){

tb<-data.frame(table(TNM_stage_subset$stagecan[i]))
tb1 <- tb[apply(tb!=0, 1, all),]
tb2 <- sum(tb1$Freq)
t_dta <- data.frame(cbind(tb2, stagecan[i]))
tda1 <- rbind(t_dta1,t_dta)

}


############
# split DE analysis table by cancer and save in RDS file 
###########




DEAnlysis <- data.frame(read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/Appkin/data/All_cancers_genes_Pvalue_FC.csv"))
head(DEAnlysis)
Cantype <- unique(DEAnlysis$Cancer)

for(i in 1:length(Cantype)){
  print(Cantype[i])
  data_can <- DEAnlysis[which(DEAnlysis$Cancer == Cantype[i]),]
  
  saveRDS(data_can, paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/Appkin/data/DEAnalysis/",Cantype[i], "_DE.rds", sep = "" ))
  
}

############
# plotly for the score analysis  
###########

library(plotly)
packageVersion('plotly')

rank_data <- data.frame(read.csv("/projects/ccs/schurerlab/Rimpi/TCGARecount2/ShinyApp/Appkin/data/Kinase_full_information_for_App.csv"))


#subset the data
subset <- rank_data[which(rank_data$TDL == "Tbio"),]

p <- plot_ly(subset, x = ~Kinase_score_percent, y = ~Cancer, name = "Kinase_score_percent", type = 'scatter',
                          mode = "markers", marker = list(color = "lightgreen"), text = ~paste("Score: ", Kinase_score_percent, '<br>Gene:', Gene)) %>%
layout(
title = "Tbio",
xaxis = list(title = "kinase_score"),
 margin = list(l = 100)
)
p


g <- list(
  scope = 'usa',
  showland = T,
  landcolor = toRGB("gray90"),
  showcountries = F,
  subunitcolor = toRGB("white")
)


one_map <- function(dat) {
  plot_ly(dat) %>%
    add_markers(x = ~Kinase_score_percent, y = ~Cancer, color = I("blue"), alpha = 0.5) %>%
    add_text(x = 0, y = 47, text = ~unique(TDL), color = I("black")) %>%
    layout(geo = g)



p <- rank_data %>%
  group_by(TDL) %>%
  do(map = one_map(.)) %>%
  subplot(nrows = 9) %>%
  layout(
    showlegend = FALSE,
    title = 'New Walmart Stores per year 1962-2006<br> Source: <a href="http://www.econ.umn.edu/~holmes/data/WalMart/index.html">University of Minnesota</a>',
    width = 1000,
    height = 900,
    hovermode = FALSE
  )



vars <- unique(rank_data$TDL)
scaterplot <- lapply(vars, function(var) {
  plot_ly(x = state.x77[, var], y = state.name) %>%
    add_bars(orientation = "h", name = var) %>%
    layout(showlegend = FALSE, hovermode = "y",
           yaxis = list(showticklabels = FALSE))
  
  
  
  p <- ggplot(rank_data, aes(Kinase_score_percent, Cancer)) + 
    geom_point() +
    facet_wrap(~ TDL) +
    ggtitle("Diamonds dataset facetted by clarity")
  
  
  
  
  p = ggplot(rank_data, aes(x=Kinase_score_percent, y=Cancer, color = Kinase_Group, text = paste("gene:",Gene))) +
   geom_point(alpha = (1/3)) +
    facet_grid(. ~ TDL, scales = "free_x") +
   labs(x='') + 
    theme(axis.text.x = element_text(angle = 90, 
                                           hjust = 1,
                                           size=3))
  p <- ggplotly(p)
  p
  
  
