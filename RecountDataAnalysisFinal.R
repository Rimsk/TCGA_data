
source('http://bioconductor.org/biocLite.R')
biocLite('recount')
library('recount')
browseVignettes('recount')
library('dplyr')
library('recount')
library('magrittr')
library('limma')
library('edgeR')
library('ffpe')
library('RSkittleBrewer')
library('SummarizedExperiment')
library('devtools')
trop <- RSkittleBrewer::RSkittleBrewer('tropical')
library('TCGAbiolinks')
library('DESeq2')
library('tidyverse')

###base file path
filePath <- "/projects/ccs/schurerlab/Rimpi/"


####output Directory
outputDir  <- "TCGARecount2"

setwd(file.path(filePath))
###IDG Kinase
IDG_kinase <- read.csv("/projects/ccs/schurerlab/Rimpi/kinases_IDG.txt", sep="\t")

download_study( 'TCGA', type='rse-gene')
load(file.path('TCGA', 'rse_gene.Rdata'))

# extract breast cancer data 
#CanID <-  as.character(unique(rse_gene$cgc_file_investigation))
CanID <- c("TCGA-LIHC","TCGA-PRAD","TCGA-READ","TCGA-BLCA",  "TCGA-BRCA","TCGA-UCEC","TCGA-KIRC","TCGA-PCPG","TCGA-LUSC",   
           "TCGA-LUAD","TCGA-STAD","TCGA-GBM", "TCGA-THCA","TCGA-CESC","TCGA-COAD"
                       ,"TCGA-HNSC", "TCGA-KICH","TCGA-ESCA","TCGA-KIRP","TCGA-CHOL")

# ("TCGA-LGG" = 514 samples in TP
"TCGA-ACC" =TP-79 
"TCGA-SARC" = NT-2  TP-259
"TCGA-DLBC" = TP - 48 
"TCGA-MESO"=TP-87 
"TCGA-THYM"= NT-2  TP-120 
"TCGA-OV"=TP-422 
"TCGA-LAML"= 126 TB (Primary Blood Derived Cancer - Peripheral Blood) 
"TCGA-UVM"=TP-80 
"TCGA-SKCM"=NT-1  TP-103 
"TCGA-PAAD"=NT-4  TP -178
"TCGA-TGCT"=TP-150 
"TCGA-ACC"=TP-79 
"TCGA-UCS"=TP-57 )



AllUpregulated <- data.frame()
#6,11,12,15,19,24,27,31
for(i in 1:length(CanID)) {
print(paste(CanID[i]))
 
#brca.rec <- TCGAquery_recount2(project = "tcga", tissue = "adrenal_gland" )paste(tissue))
exp <- rse_gene[, rse_gene$cgc_file_investigation %in% c(paste(CanID[i]))]

#colData(brca.rec) <-colData(brca.rec)[ , ! apply( colData(brca.rec) , 2 , function(x) all(is.na(x)) ) ]

###remove duplicated  TCGA barcodes from  
clic.exp <- data.frame(colData(exp))
clic.exp <- clic.exp[ , ! apply( clic.exp , 2 , function(x) all(is.na(x)) ) ]
samples <- c("Primary Tumor", "Solid Tissue Normal")
clic.exp <- clic.exp[clic.exp$gdc_cases.samples.sample_type %in% samples,]
clic.exp$Shortnametype <- ifelse(clic.exp$gdc_cases.samples.sample_type == "Solid Tissue Normal", "NT", "TP")
# data_TP <- as.character(clic.exp$gdc_cases.samples.portions.analytes.aliquots.submitter_id[which(clic.exp$Shortnametype == "TP")])
# data_NT <- as.character(clic.exp$gdc_cases.samples.portions.analytes.aliquots.submitter_id[which(clic.exp$Shortnametype == "NT")])
# keepID <- setdiff(data_TP, data_NT)
# 
# #pick up the candidates in normal and tumor
# Primary_T <- clic.exp[which(clic.exp$Shortnametype == "TP"),]
# Normal <- clic.exp[which(clic.exp$Shortnametype == "NT"),]
# 
# #keep those the unique ones in tumor data 
# Primary_T1 <- subset(Primary_T, gdc_cases.samples.portions.analytes.aliquots.submitter_id %in% keepID)
# Primary_T1 <- Primary_T1[!duplicated(Primary_T1$gdc_cases.samples.portions.analytes.aliquots.submitter_id),]
# 
# #bind both the both normal and tumor data set
# clinical <- rbind(Primary_T1, Normal)
clinical <- clic.exp

save(clinical, file = "./TOMsimilarity.signed.RData")



clinical1 <- data.frame(barcode = clinical$gdc_cases.samples.portions.analytes.aliquots.submitter_id, caseId = rownames(clinical), bigfile = clinical$bigwig_file, tissuetype = clinical$cgc_file_investigation,
                        project = clinical$project, clinical$reads_downloaded,
                        experimentalStrategy = clinical$cgc_file_experimental_strategy,
                        pairedend = clinical$paired_end, mappedreads= clinical$mapped_read_count,
                        auc = clinical$auc, filesize= clinical$gdc_file_size ,
                        access =clinical$gdc_access, platform = clinical$gdc_platform, state = clinical$gdc_state, category = clinical$gdc_data_category,
                        center = clinical$gdc_center.short_name,centertype= clinical$gdc_center.center_type, gender = clinical$gdc_cases.demographic.gender,
                        birthyear = clinical$gdc_cases.demographic.year_of_birth, race = clinical$gdc_cases.demographic.race,
                        deathyear = clinical$gdc_cases.demographic.year_of_death, tissue = clinical$gdc_cases.project.primary_site,
                        stage = clinical$gdc_cases.diagnoses.tumor_stage, vitalstatus= clinical$gdc_cases.diagnoses.vital_status,
                        sampletype= clinical$gdc_cases.samples.sample_type, lastfollowupdays= clinical$cgc_case_days_to_last_follow_up,
                        followupID = clinical$cgc_follow_up_id,patientid= clinical$xml_bcr_patient_barcode,drugtherapy= clinical$cgc_drug_therapy_pharmaceutical_therapy_type,
                        site= clinical$xml_icd_o_3_site, followuptumorstatus= clinical$cgc_follow_up_tumor_status,
                        Daytstobirth= clinical$xml_days_to_birth, icdhistology= clinical$xml_icd_o_3_histology, tumorstatus= clinical$cgc_case_tumor_status,
                        ageofdiagnosis= clinical$cgc_case_age_at_diagnosis, tissuesourcecite= clinical$cgc_sample_tissue_source_site_code, samplecreationdatetime = clinical$gdc_cases.samples.portions.creation_datetime,
                        tissueoforigin= clinical$gdc_cases.diagnoses.tissue_or_organ_of_origin, diagnosisdaystobirth = clinical$gdc_cases.diagnoses.days_to_birth,
                        projectaccession = clinical$gdc_cases.project.program.dbgap_accession_number,
                        Shortnametype = clinical$Shortnametype, tissuesite= clinical$gdc_cases.tissue_source_site.project, daystodeath = clinical$xml_days_to_death)


clinical2 <- clinical1[!duplicated(clinical1$barcode),]








if(length(data_NT) > 3 & length(data_TP) > 3) {



###remove NAs from 
exp.genes <- data.frame(rowData(exp))
table(is.na(exp.genes$symbol))



exp.genes.sort <- exp.genes[which(apply(!is.na(exp.genes), 1, all)),]
table(is.na(exp.genes.sort$symbol))

dim(exp.genes.sort[duplicated(exp.genes.sort$symbol),])
exp.genes.sort.dup <- exp.genes.sort[!duplicated(exp.genes.sort$symbol),]






#breastCancer1[duplicated(breastCancer1$barcode),]
#breastCancer3 <- breastCancer2[!duplicated(breastCancer2$barcode),]
#breastCancer2 <- breastCancer[!duplicated(breastCancer$gdc_cases.samples.portions.analytes.aliquots.submitter_id),]
#breastCancer3 <- breastCancer[!duplicated(breastCancer$gdc_cases.submitter_id),]
table(clinical2$sampletype)


####get the expression counts
exp.counts <- assay(exp)

####Remove the expression counts id not present in genesN
exp.counts.sort <- exp.counts[which(rownames(exp.counts) %in% exp.genes.sort.dup$gene_id),]


####Remove the expression counts barcodes not present in breastCancer2
#which(colnames(expressioncounts1) == "137FDD43-DC3C-4027-ADED-DA9D69CF909C")
#which(colnames(expressioncounts1) == "1E347C7A-FBE8-408D-8014-EE33869116E0")

rm <- setdiff(colnames(exp.counts.sort),clinical2$caseId)
expr.count.fin <- exp.counts.sort[ , !(colnames(exp.counts.sort) %in% rm)]


#add rownames as gene symbol
rownames(expr.count.fin) <- exp.genes.sort.dup$symbol
colnames(expr.count.fin) <- clinical2$barcode

#
rownames(clinical2) <- clinical2$barcode

dir.create(file.path(filePath, outputDir, paste(CanID[i])))
setwd(file.path(filePath, outputDir, paste(CanID[i])))


saveRDS(clinical, file=paste(getwd(),'/', CanID[i],"_clinical.rds", sep = ""))
write.table(clinical, file=paste(getwd(),'/', CanID[i],"_clinical.csv", sep="" ), sep=",")
#saveRDS(breastCancer2, file = paste(getwd(),'/',"_clinic.Rdata", sep="" ))

saveRDS(expr.count.fin, file = paste(getwd(),'/', "_expression.Rdata", sep="" ))
write.table(expr.count.fin, paste(getwd(),'/', CanID[i], "_expression.csv", sep="" ), sep=",")

rownames(exp.genes.sort.dup) <- exp.genes.sort.dup$symbol

#rownames(clinical1) <- rownames(clinical)
# convert into SummarizedExperiment
data <- SummarizedExperiment(assays = list(counts=expr.count.fin),
                    rowData= exp.genes.sort.dup, colData=clinical2)
data1 <- as(data, "RangedSummarizedExperiment")
saveRDS(data1, file = paste(getwd(),'/',CanID[i],".rds", sep="" ))

readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/TCGA-LIHC.rds")
gdc_cases.samples.portions.submitter_id

#DE Analysis
dataTP <- as.character(clinical2$barcode[which(clinical2$Shortnametype == "TP")])
dataNT <- as.character(clinical2$barcode[which(clinical2$Shortnametype == "NT")])

dataPrep <- TCGAanalyze_Preprocenamesss(object = data1, cor.cut = 0.6)

#dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                      geneInfo = geneInfo,
#                                    method = "gcContent")   

rownames(exp.genes.sort.dup) <- exp.genes.sort.dup$gene_id
#normalize read counts
Normdd <- scale_counts(data1)
Normddval <- as.matrix(assays(Normdd)$counts)

if(paste(CanID[i]) =="TCGA-LIHC"){
  LIHC.clinical <- clinical
  LIHC.Normalize.expression <-Normddval
  print("completed LIHC clinical")
} else if (paste(CanID[i]) == "TCGA-PRAD"){
  PRAD.clinical <- clinical
  PRAD.Normalize.expression <-Normddval
  print("completed PRAD clinical")
  
} else if (paste(CanID[i]) == "TCGA-READ"){
  READ.clinical <- clinical
  READ.Normalize.expression <-Normddval
  print("completed READ clinical")
  
} else if (paste(CanID[i]) == "TCGA-BLCA"){
  BLCA.clinical <- clinical
  BLCA.Normalize.expression <-Normddval
  print("completed BLCA clinical")
  
}  else if (paste(CanID[i]) == "TCGA-BRCA"){
  BRCA.clinical <- clinical
  BRCA.Normalize.expression <-Normddval
  print("completed BRCA clinical")
  
} else if (paste(CanID[i]) == "TCGA-UCEC"){
  UCEC.clinical <- clinical
  UCEC.Normalize.expression <-Normddval
  print("completed UCEC clinical")
  
} else if (paste(CanID[i]) == "TCGA-KIRC"){
  KIRC.clinical <- clinical
  KIRC.Normalize.expression <-Normddval
  print("completed KIRC clinical")
  
} else if (paste(CanID[i]) == "TCGA-PCPG"){
  PCPG.clinical <- clinical
  PCPG.Normalize.expression <-Normddval
  print("completed PCPG clinical")
  
} else if (paste(CanID[i]) == "TCGA-LUSC"){
  LUSC.clinical <- clinical
  LUSC.Normalize.expression <-Normddval
  print("completed LUSC clinical")
  
} else if (paste(CanID[i]) == "TCGA-LUAD"){
  LUAD.clinical <- clinical
  LUAD.Normalize.expression <-Normddval
  print("completed LUAD clinical")
  
} else if (paste(CanID[i]) == "TCGA-STAD"){
  STAD.clinical <- clinical
  STAD.Normalize.expression <-Normddval
  print("completed STAD clinical")
  
} else if (paste(CanID[i]) == "TCGA-GBM"){
  GBM.clinical <- clinical
  GBM.Normalize.expression <-Normddval
  print("completed GBM clinical")
  
} else if (paste(CanID[i]) == "TCGA-THCA"){
  THCA.clinical <- clinical
  THCA.Normalize.expression <-Normddval
  print("completed THCA clinical")
  
} else if (paste(CanID[i]) == "TCGA-CESC"){
  CESC.clinical <- clinical
  CESC.Normalize.expression <-Normddval
  print("completed CESC clinical")
  
} else if (paste(CanID[i]) == "TCGA-COAD"){
  COAD.clinical <- clinical
  COAD.Normalize.expression <-Normddval
  print("completed COAD clinical")
  
} else if (paste(CanID[i]) == "TCGA-HNSC"){
  HNSC.clinical <- clinical
  HNSC.Normalize.expression <-Normddval
  print("completed HNSC clinical")
  
} else if (paste(CanID[i]) == "TCGA-KICH"){
  KICH.clinical <- clinical
  KICH.Normalize.expression <-Normddval
  print("completed KICH clinical")
  
} else if (paste(CanID[i]) == "TCGA-ESCA"){
  ESCA.clinical <- clinical
  ESCA.Normalize.expression <-Normddval
  print("completed ESCA clinical")
  
} else if (paste(CanID[i]) == "TCGA-KIRP"){
  KIRP.clinical <- clinical
  KIRP.Normalize.expression <-Normddval
  print("completed KIRP clinical")
  
} else if (paste(CanID[i]) == "TCGA-CHOL"){
  CHOL.clinical <- clinical
  CHOL.Normalize.expression <-Normddval
  print("completed CHOL clinical")
  
}else  {print("nothing is working") }






dataFilt <- TCGAanalyze_Filtering(tabDF = Normddval,
                                  method = "quantile", 
                                  qnt.cut =  0.25) 

#


#new_dataFilt <- t(dataFilt)

#write.csv(new_dataFilt, file= paste(getwd(),'/', CanID[i],"new_dataFilt.csv", sep="" ), sep=",")
#DE
#####################
#Normalized and filtered data for all cancers is saved in "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis" 
#####################

#save(LIHC.clinical, LIHC.Normalize.expression, LIHC.Normalize.expression_filtered , file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LIHC_Normalize_filter.RData")



#############
#differentially expressed data
#######

print(paste("Differentially Expressed", CanID[i]))
dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataNT],
                            mat2 = dataFilt[,dataTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT") 

print(dim(dataDEGs)) #4888

write.table(dataDEGs, file= paste(getwd(),'/', CanID[i],"_DE.csv", sep="" ), sep=",")

print(paste("Diffirrentially expressed genes in ",  CanID[i], " - ", sep=""))
print(paste(dim(dataDEGs)))
#nam<-paste("DE", CanID[6],sep=".")
#assign(CanID[i], dataDEGs)

print(paste("Number of common candidates in", CanID[i], "and IDG", sep="-"))
print(paste(table(rownames(dataDEGs) %in% IDG_kinase$HGNC.Sym)))

commomIDG <- dataDEGs[intersect(rownames(dataDEGs),IDG_kinase$HGNC.Sym),]
print(paste("common DE candidates in IDG and ", CanID[i], " = ", dim(commomIDG)[1], sep=""))

write.table(commomIDG, file= paste(getwd(),'/', CanID[i],"commomIDG.csv", sep="" ), sep=",")
print(paste("completed ", CanID[i], sep=""))

#up regulated genes 

Upregulated <- commomIDG[which(commomIDG$logFC > 1),]
Upregulated$genes <- rownames(Upregulated)
Upregulated$cancer <- CanID[i]
AllUpregulated <- rbind(AllUpregulated, Upregulated)

write.table(Upregulated, file= paste(getwd(),'/', CanID[i],"upregulatedIDG.csv", sep="" ), sep=",")



# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(commomIDG,"Tumor","Normal",
                                          dataFilt[,dataTP],dataFilt[,dataNT])

write.table(dataDEGsFiltLevel, file= paste(getwd(),'/', CanID[i],"dataDEGsFiltLevel.csv", sep="" ), sep=",")
} else { next }

}



####################
# Get the common upregulated genes in all the cancers 
################

keep <-  names(which(table(AllUpregulated$genes) < 2))
Overlap <- AllUpregulated[!AllUpregulated$genes %in% keep,]

write.table(Overlap, file= paste(getwd(),'/',"Overlap_understudied_upregulated_kinases.csv", sep="" ))
saveRDS(Overlap, file= paste(getwd(),'/',"Overlap_understudied_upregulated_kinases.RData", sep="" ))

cbPalette <- c("#F0A3FF", "#99d8c9", "#9ebcda","#fdae6b","#2BCE48","#FFCC99",
               "#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380",
               "#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF", "#8da0cb")

png("/projects/ccs/schurerlab/Rimpi/TCGARecount2/allGenes.jpg", width = 1000, height = 800)
p <- ggplot(Overlap, aes(x = fct_infreq(genes), fill=cancer, order= fct_infreq(genes), width=.25)) + geom_bar(lwd=1, color="white", width = 0.75)
p1 <- p+ scale_fill_manual(values=cbPalette)  + theme_classic(base_size = 5, base_family = "") #position = position_stack(reverse = TRUE))
p1 <- p1+theme( legend.title = element_text(size=20, face="bold",  ), 
                 legend.text = element_text(size=16),text = element_text(size=8), 
  axis.text.x = element_text(face="bold", color="black", 
                                   size=10, angle=90),  
        axis.text.y = element_text(face="bold", color="black", 
                                   size=10))
print(p1)
dev.off()

##################
#Survival of specific genes which are prsent in more than 4 times
############


Genes.cut <- Overlap[grep("ADCK5|ALPK2|ALPK3|CAMK1D|CAMK1G|CAMKK1|CAMKV|CDC42BPG|CDK18|DYRK2|ERN2|ITPKA|LMTK3|MAP3K10|MAPK15|NEK5|PAK3|PHKA1|PHKG2|PKMYT1|PLK5|PNCK|POMK|RIOK1|RPS6KL1|SGK223|SGK494|SRPK3|STK17A|STK3|STK31|STK32A|STK32C|STK36|TSSK3|TSSK6|TTBK1|WNK2", Overlap$genes),]

Genes.cut <- Overlap[grep("ADCK5|ALPK2|ALPK3|CAMK1D|CAMK1G|CAMKK1|CAMKV|CDC42BPG|CDK18|DYRK2|ERN2|ITPKA|LMTK3|MAP3K10|MAPK15|NEK5|PAK3|PHKA1|PHKG2|PKMYT1|PLK5|PNCK|POMK|RIOK1|RPS6KL1|SGK223|SGK494|SRPK3|STK17A|STK3|STK31|STK32A|STK32C|STK36|TSSK3|TSSK6|TTBK1|WNK2"  , Overlap$genes),]


#LIHC
LIHC.clin <- data.frame(patient.bcr_patient_barcode = LIHC.clinical$gdc_cases.submitter_id,
patient.vital_status = LIHC.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = LIHC.clinical$cgc_case_days_to_last_follow_up, 
patient.days_to_death = LIHC.clinical$xml_days_to_death)






##################
#Survival 
############
  
#gdc_cases.samples.portions.analytes.aliquots.submitter_id
  

survivalTCGA(LIHC.clin) -> LIHC.clinical.surv


LIHC.Normalized <- data.frame(t(LIHC.Normalize.expression))
LIHC.Normalized <- cbind(bcr_patient_barcode = rownames(LIHC.Normalized), LIHC.Normalized)
rownames(LIHC.Normalized) <- NULL
LIHC.rnaseq <- LIHC.Normalized

expressionsTCGA(
  LIHC.rnaseq,
  extract.cols = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
                   "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
                   "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
                   "PLK5","RPS6KL1", "SGK223")) %>%
  rename(cohort = dataset,
        CAMKK1 = 'CAMKK1', CAMK1G= 'CAMK1G', PAK3 = 'PAK3', STK3 = 'STK3', CDK18= 'CDK18',
         RIOK1 = 'RIOK1', DYRK2 = 'DYRK2', PKMYT1= 'PKMYT1', STK33= 'STK33' , MAP3K10 = 'MAP3K10', PNCK= 'PNCK',
         ERN2='ERN2', ALPK3= 'ALPK3', ITPKA='ITPKA', LMTK3= 'LMTK3', STK32B='STK32B', TSSK3='TSSK3', STK36='STK36',
         CAMKV= 'CAMKV', STK17A='STK17A', WNK2 = 'WNK2',  STK32C='STK32C' ,SGK494='SGK494',CDC42BPG='CDC42BPG',
         ADCK5='ADCK5', TSSK6='TSSK6', MAPK15='MAPK15', STK31='STK31',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',
         SGK223='SGK223') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("CAMKK1", "CAMK1G", "PAK3", "STK3", "CDK18", "RIOK1", "DYRK2", "PKMYT1","STK33",
  "MAP3K10", "PNCK", "ERN2", "ALPK3", "ITPKA", "LMTK3", "STK32B","TSSK3", "STK36", "CAMKV",
  "STK17A","WNK2", "STK32C", "SGK494", "CDC42BPG", "ADCK5", "TSSK6", "MAPK15","STK31", "SRPK3",
  "PLK5","RPS6KL1", "SGK223")

dataSurvTable <- data.frame()
dataSurvplot <- data.frame()
pvalue <- data.frame()
fittable <- data.frame()
surmedian <- data.frame()
######################
#Plot SURVIVAL CURVES FOR GENES
#######################

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
  
 
  
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LIHC/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
"n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                       "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  #dev.off()
  
  print(paste("completed", geneid[i]))
  
}

print("##################
#Survival ANALYSIS OF TCGA-PRAD GENES
##############################")

########################################

#CLINICAL DATA OF TCGA-PRAD

########################################
PRAD.clin <- data.frame(patient.bcr_patient_barcode = PRAD.clinical$gdc_cases.submitter_id,
                        patient.vital_status = PRAD.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = PRAD.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = PRAD.clinical$xml_days_to_death)









survivalTCGA(PRAD.clin) -> PRAD.clinical.surv
####################################
#get the normaliized gene expression
#####################################
PRAD.Normalized <- data.frame(t(PRAD.Normalize.expression))
PRAD.Normalized <- cbind(bcr_patient_barcode = rownames(PRAD.Normalized), PRAD.Normalized)
rownames(PRAD.Normalized) <- NULL
PRAD.rnaseq <- PRAD.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############
expressionsTCGA(
  PRAD.rnaseq,
  extract.cols = c("CAMKV", "NEK5")) %>%
  rename(cohort = dataset,
         CAMKV = `CAMKV`, NEK5='NEK5') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("CAMKV", "NEK5")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  PRAD.surv_rnaseq.cut <- surv_cutpoint(
    PRAD.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(PRAD.surv_rnaseq.cut)
  
  distri <- plot(PRAD.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PRAD/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  PRAD.surv_rnaseq.cat <- surv_categorize(PRAD.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = PRAD.surv_rnaseq.cat)
  
  
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PRAD/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  pval <- surv_pvalue(fit)
  
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}





print("##################
#Survival ANALYSIS OF TCGA-READ GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-READ

########################################
#READ
READ.clin <- data.frame(patient.bcr_patient_barcode = READ.clinical$gdc_cases.submitter_id,
                        patient.vital_status = READ.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = READ.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = READ.clinical$xml_days_to_death)


survivalTCGA(READ.clin) -> READ.clinical.surv

####################################
#get the normaliized gene expression
#####################################

READ.Normalized <- data.frame(t(READ.Normalize.expression))
READ.Normalized <- cbind(bcr_patient_barcode = rownames(READ.Normalized), READ.Normalized)
rownames(READ.Normalized) <- NULL
READ.rnaseq <- READ.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  READ.rnaseq,
  extract.cols = c("PHKA1","PKMYT1","MAP3K10","LMTK3","TTBK1","STK36",
                   "CAMKV","SGK494","ADCK5","TSSK6","MAPK15","STK31","NEK5","RPS6KL1")) %>%
  rename(cohort = dataset,
         PHKA1= 'PHKA1',PKMYT1='PKMYT1',MAP3K10='MAP3K10',LMTK3='LMTK3',TTBK1='TTBK1',STK36='STK36',
         CAMKV='CAMKV',SGK494='SGK494',ADCK5='ADCK5',TSSK6='TSSK6',MAPK15='MAPK15',STK31='STK31',NEK5='NEK5',
         RPS6KL1='RPS6KL1') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("PHKA1","PKMYT1","MAP3K10","LMTK3","TTBK1","STK36",
           "CAMKV","SGK494","ADCK5","TSSK6","MAPK15","STK31","NEK5","RPS6KL1")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  READ.surv_rnaseq.cut <- surv_cutpoint(
    READ.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(READ.surv_rnaseq.cut)
  
  distri <- plot(READ.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-READ/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  print(distri)
  
  
  #dev.off()
  
  READ.surv_rnaseq.cat <- surv_categorize(READ.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = READ.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-READ/survival/", geneid[i], ".survival.jpg", sep=""))
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
  pval <- surv_pvalue(fit)
  
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  #print(gg)
  #dev.off()
  
  print(paste("completed", geneid[i]))
  
}









#########################

print("##################
#Survival ANALYSIS OF TCGA-BLCA GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-BLCA

########################################

#BLCA
BLCA.clin <- data.frame(patient.bcr_patient_barcode = BLCA.clinical$gdc_cases.submitter_id,
                        patient.vital_status = BLCA.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = BLCA.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = BLCA.clinical$xml_days_to_death)


survivalTCGA(BLCA.clin) -> BLCA.clinical.surv
####################################
#get the normaliized gene expression
#####################################

BLCA.Normalized <- data.frame(t(BLCA.Normalize.expression))
BLCA.Normalized <- cbind(bcr_patient_barcode = rownames(BLCA.Normalized), BLCA.Normalized)
rownames(BLCA.Normalized) <- NULL
BLCA.rnaseq <- BLCA.Normalized

################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  BLCA.rnaseq,
  extract.cols = c("PKMYT1","ITPKA","TTBK1","CAMKV","STK32C",
                   "SGK494","ADCK5","TSSK6","MAPK15","STK31","SGK223")) %>%
  rename(cohort = dataset,
         PKMYT1='PKMYT1',ITPKA='ITPKA',TTBK1='TTBK1',CAMKV='CAMKV',STK32C='STK32C',
         SGK494='SGK494',ADCK5='ADCK5',TSSK6='TSSK6',MAPK15='MAPK15',STK31='STK31',SGK223='SGK223') %>%
  filter(substr(bcr_patient_barcode, 14, 15) == "01") %>% 
  # only cancer samples
  mutate(bcr_patient_barcode = 
           substr(bcr_patient_barcode, 1, 12)) -> BLCA.rnaseq1



BLCA.clinical.surv %>%
  left_join(BLCA.rnaseq1,
            by = "bcr_patient_barcode") ->
  BLCA.surv_rnaseq


table(BLCA.surv_rnaseq$cohort, useNA = "always")


BLCA.surv_rnaseq <- BLCA.surv_rnaseq %>%
  filter(!is.na(cohort))


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("PKMYT1","ITPKA","TTBK1","CAMKV","STK32C",
           "SGK494","ADCK5","TSSK6","MAPK15","STK31","SGK223")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  BLCA.surv_rnaseq.cut <- surv_cutpoint(
    BLCA.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(BLCA.surv_rnaseq.cut)
  
  distri <- plot(BLCA.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-BLCA/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  BLCA.surv_rnaseq.cat <- surv_categorize(BLCA.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = BLCA.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-BLCA/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  pval <- surv_pvalue(fit)
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}


##################################
print("##################
#Survival ANALYSIS OF TCGA-BRCA GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-BRCA

########################################
#BRCA
BRCA.clin <- data.frame(patient.bcr_patient_barcode = BRCA.clinical$gdc_cases.submitter_id,
                        patient.vital_status = BRCA.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = BRCA.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = BRCA.clinical$xml_days_to_death)


survivalTCGA(BRCA.clin) -> BRCA.clinical.surv

####################################
#get the normaliized gene expression
#####################################

BRCA.Normalized <- data.frame(t(BRCA.Normalize.expression))
BRCA.Normalized <- cbind(bcr_patient_barcode = rownames(BRCA.Normalized), BRCA.Normalized)
rownames(BRCA.Normalized) <- NULL
BRCA.rnaseq <- BRCA.Normalized

################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  BRCA.rnaseq,
  extract.cols = c("CAMK1G","PKMYT1","PNCK","ERN2","ITPKA","LMTK3","TTBK1","PHKG2","CAMKV","ADCK5","MAPK15","STK31","NEK5")) %>%
  rename(cohort = dataset,
         CAMK1G= 'CAMK1G', PKMYT1= 'PKMYT1', PNCK ='PNCK', ERN2='ERN2', ITPKA='ITPKA', LMTK3='LMTK3',
         TTBK1='TTBK1', PHKG2='PHKG2',CAMKV='CAMKV',ADCK5='ADCK5',MAPK15='MAPK15',STK31='STK31',
         NEK5='NEK5') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("CAMK1G","PKMYT1","PNCK","ERN2","ITPKA","LMTK3","TTBK1","PHKG2","CAMKV","ADCK5","MAPK15","STK31","NEK5")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  BRCA.surv_rnaseq.cut <- surv_cutpoint(
    BRCA.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(BRCA.surv_rnaseq.cut)
  
  distri <- plot(BRCA.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-BRCA/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  BRCA.surv_rnaseq.cat <- surv_categorize(BRCA.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = BRCA.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-BRCA/survival/", geneid[i], ".survival.jpg", sep=""))
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
 # print(gg)
 # dev.off()
  pval <- surv_pvalue(fit)
  
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}


###
#save the normalized filtered files in clinical file 
######


#save the data in file normalized and filtered and clinical data


LIHC.clinical 
LIHC.clinical$gdc_cases.annotations <- NULL
LIHC.clinical$gdc_cases.samples.portions.slides <- NULL
LIHC.clinical$gdc_annotations <- NULL
LIHC.clinical$gdc_cases.samples.portions.analytes.aliquots.annotations <- NULL
LIHC.Normalize.expression 
LIHC.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = LIHC.Normalize.expression,
                                                            method = "quantile", 
                                                            qnt.cut =  0.25) 
write.table(LIHC.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LIHC.clinical.csv", sep=",")
save(LIHC.clinical, LIHC.Normalize.expression, LIHC.Normalize.expression_filtered , file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LIHC_Normalize_filter.RData")


PRAD.clinical 
PRAD.clinical$gdc_annotations <- NULL
PRAD.clinical$gdc_cases.annotations <- NULL
PRAD.clinical$gdc_cases.samples.portions.slides <- NULL
PRAD.Normalize.expression
PRAD.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF =PRAD.Normalize.expression, method = "quantile", 
                                                            qnt.cut =  0.25) 
write.table(PRAD.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/PRAD.clinical.csv", sep=",")
save(PRAD.clinical, PRAD.Normalize.expression,PRAD.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/PRAD_Normalize_filter.RData")



READ.clinical 
t_READ.clinical <- t(READ.clinical)
READ.clinical$gdc_cases.samples.portions.slides <- NULL
READ.Normalize.expression 
READ.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF =READ.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(READ.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/READ.clinical.csv", sep=",")
save(READ.clinical, READ.Normalize.expression, READ.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/READ_Normalize_filter.RData" )

BLCA.clinical 
BLCA.clinical$gdc_annotations <- NULL
BLCA.clinical$gdc_cases.annotations <- NULL
BLCA.clinical$gdc_cases.samples.portions.slides <- NULL
BLCA.Normalize.expression 
BLCA.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF =BLCA.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(BLCA.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/BLCA.clinical.csv", sep=",")
save(BLCA.clinical, BLCA.Normalize.expression, BLCA.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/BLCA_Normalize_filter.RData")

BRCA.clinical
BRCA.clinical$gdc_annotations <- NULL
BRCA.clinical$gdc_cases.annotations <- NULL
BRCA.clinical$gdc_cases.samples.portions.slides <- NULL
BRCA.Normalize.expression 
BRCA.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = BRCA.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)

write.table(BRCA.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/BRCA.clinical.csv", sep=",")
save(BRCA.clinical,BRCA.Normalize.expression ,BRCA.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/BRCA_Normalize_filter.RData")


UCEC.clinical
#UCEC.clinical$gdc_annotations <- NULL
#UCEC.clinical$gdc_cases.annotations <- NULL
UCEC.clinical$gdc_cases.samples.portions.slides <- NULL
UCEC.Normalize.expression 
UCEC.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF =UCEC.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(UCEC.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/UCEC.clinical.csv", sep=",")
save(UCEC.clinical,UCEC.Normalize.expression, UCEC.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/UCEC_Normalize_filter.RData")


KIRC.clinical 
KIRC.clinical$gdc_cases.annotations <- NULL
KIRC.clinical$gdc_annotations <- NULL
KIRC.clinical$gdc_cases.samples.portions.slides <- NULL
KIRC.Normalize.expression 
KIRC.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = KIRC.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)

save(KIRC.clinical, KIRC.Normalize.expression, KIRC.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/KIRC_Normalize_filter.RData")

PCPG.clinical 
PCPG.clinical$gdc_annotations <- NULL
PCPG.clinical$gdc_cases.samples.portions.slides <- NULL 

PCPG.Normalize.expression 
PCPG.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = PCPG.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)

write.table(PCPG.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/PCPG.clinical.csv", sep=",")
save(PCPG.clinical,PCPG.Normalize.expression, PCPG.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/PCPG_Normalize_filter.RData")

LUSC.clinical
LUSC.clinical$gdc_annotations <- NULL
LUSC.clinical$gdc_cases.annotations <- NULL
LUSC.clinical$gdc_cases.samples.portions.slides <- NULL
LUSC.clinical$gdc_cases.samples.portions.analytes.annotations <- NULL
LUSC.clinical$gdc_cases.samples.portions.analytes.aliquots.annotations <- NULL
LUSC.Normalize.expression 
LUSC.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF =LUSC.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)

write.table(LUSC.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LUSC.clinical.csv", sep=",")
save(LUSC.clinical, LUSC.Normalize.expression ,LUSC.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LUSC_Normalize_filter.RData")   


LUAD.clinical 
t_LUAD.clinical <- t(LUAD.clinical)
LUAD.clinical$gdc_annotations <- NULL
LUAD.clinical$gdc_cases.annotations <- NULL
LUAD.clinical$gdc_cases.samples.portions.slides <- NULL
LUAD.Normalize.expression 
LUAD.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = LUAD.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(LUAD.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LUAD.clinical.csv", sep=",")
save(LUAD.clinical, LUAD.Normalize.expression, LUAD.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LUAD_Normalize_filter.RData") 

STAD.clinical 
t_STAD.clinical <- t(STAD.clinical)
STAD.clinical$gdc_cases.samples.portions.slides <- NULL

STAD.Normalize.expression 
STAD.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = STAD.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)

write.table(STAD.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/STAD.clinical.csv", sep=",")
save(STAD.clinical ,STAD.Normalize.expression, STAD.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/STAD_Normalize_filter.RData")

GBM.clinical 
t_GBM.clinical <- t(GBM.clinical)
GBM.clinical$gdc_annotations <- NULL
GBM.clinical$gdc_cases.annotations <- NULL
GBM.clinical$gdc_cases.samples.portions.slides <- NULL
GBM.clinical$gdc_cases.samples.portions.analytes.annotations <- NULL
GBM.Normalize.expression 
GBM.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = GBM.Normalize.expression,method = "quantile", 
                                                           qnt.cut =  0.25)
write.table(GBM.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/GBM.clinical.csv", sep=",")
save(GBM.clinical, GBM.Normalize.expression, GBM.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/GBM_Normalize_filter.RData")

THCA.clinical 
t_THCA.clinical <- t(THCA.clinical)
THCA.clinical$gdc_annotations <- NULL
THCA.clinical$gdc_cases.annotations <- NULL
THCA.clinical$gdc_cases.samples.portions.slides <- NULL
THCA.Normalize.expression 
THCA.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = THCA.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(THCA.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/THCA.clinical.csv", sep=",")
save(THCA.clinical,THCA.Normalize.expression, THCA.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/THCA_Normalize_filter.RData")

CESC.clinical 
t_CESC.clinical <- t(CESC.clinical)
#CESC.clinical$gdc_annotations <- NULL
#CESC.clinical$gdc_cases.annotations <- NULL
CESC.clinical$gdc_cases.samples.portions.slides <- NULL
CESC.Normalize.expression 
CESC.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = CESC.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(CESC.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/CESC.clinical.csv", sep=",")
save(CESC.clinical, CESC.Normalize.expression, CESC.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/CESC_Normalize_filter.RData")

COAD.clinical 
t_COAD.clinical <- t(COAD.clinical)
#COAD.clinical$gdc_annotations <- NULL
#COAD.clinical$gdc_cases.annotations <- NULL
COAD.clinical$gdc_cases.samples.portions.slides <- NULL
COAD.Normalize.expression 
COAD.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = COAD.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(COAD.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/COAD.clinical.csv", sep=",")
save(COAD.clinical, COAD.Normalize.expression, COAD.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/COAD_Normalize_filter.RData")

HNSC.clinical
t_HNSC.clinical <- t(HNSC.clinical)
#HNSC.clinical$gdc_annotations <- NULL
#HNSC.clinical$gdc_cases.annotations <- NULL
HNSC.clinical$gdc_cases.samples.portions.slides <- NULL
HNSC.Normalize.expression 
HNSC.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = HNSC.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(HNSC.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/HNSC.clinical.csv", sep=",")
save(HNSC.clinical, HNSC.Normalize.expression, HNSC.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/HNSC_Normalize_filter.RData")

KICH.clinical
t_KICH.clinical <- t(KICH.clinical)
KICH.clinical$gdc_annotations <- NULL
KICH.clinical$gdc_cases.annotations <- NULL
KICH.clinical$gdc_cases.samples.portions.slides <- NULL
KICH.Normalize.expression 
KICH.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = KICH.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(KICH.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/KICH.clinical.csv", sep=",")
save(KICH.clinical, KICH.Normalize.expression, KICH.Normalize.expression_filtered , file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/KICH_Normalize_filter.RData")   

ESCA.clinical 
t_ESCA.clinical <- t(ESCA.clinical)
ESCA.clinical$gdc_annotations <- NULL
ESCA.clinical$gdc_cases.annotations <- NULL
ESCA.clinical$gdc_cases.samples.portions.slides <- NULL
ESCA.Normalize.expression 
ESCA.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = ESCA.Normalize.expression ,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(ESCA.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/ESCA.clinical.csv", sep=",")
save(ESCA.clinical, ESCA.Normalize.expression, ESCA.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/ESCA_Normalize_filter.RData")                                    

KIRP.clinical
t_KIRP.clinical <- t(KIRP.clinical)
KIRP.clinical$gdc_annotations <- NULL
KIRP.clinical$gdc_cases.annotations <- NULL
KIRP.clinical$gdc_cases.samples.portions.slides <- NULL
KIRP.Normalize.expression
KIRP.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = KIRP.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(KIRP.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/KIRP.clinical.csv", sep=",")
save(KIRP.clinical, KIRP.Normalize.expression, KIRP.Normalize.expression_filtered , file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/KIRP_Normalize_filter.RData")                                     

CHOL.clinical 
t_CHOL.clinical <- t(CHOL.clinical)
#CHOL.clinical$gdc_annotations <- NULL
#CHOL.clinical$gdc_cases.annotations <- NULL
CHOL.clinical$gdc_cases.samples.portions.slides <- NULL
CHOL.Normalize.expression 
CHOL.Normalize.expression_filtered <- TCGAanalyze_Filtering(tabDF = CHOL.Normalize.expression,method = "quantile", 
                                                            qnt.cut =  0.25)
write.table(CHOL.clinical, file="/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/CHOL.clinical.csv", sep=",")
save(CHOL.clinical, CHOL.Normalize.expression, CHOL.Normalize.expression_filtered, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/CHOL_Normalize_filter.RData")

#######################

print("##################
#Survival ANALYSIS OF TCGA-UCEC GENES
      ##############################")

########################################

#CLINICAL DATA OF UCEC-PRAD

########################################
#UCEC
UCEC.clin <- data.frame(patient.bcr_patient_barcode = UCEC.clinical$gdc_cases.submitter_id,
                        patient.vital_status = UCEC.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = UCEC.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = UCEC.clinical$xml_days_to_death)


survivalTCGA(UCEC.clin) -> UCEC.clinical.surv
####################################
#get the normaliized gene expression
#####################################

UCEC.Normalized <- data.frame(t(UCEC.Normalize.expression))
UCEC.Normalized <- cbind(bcr_patient_barcode = rownames(UCEC.Normalized), UCEC.Normalized)
rownames(UCEC.Normalized) <- NULL
UCEC.rnaseq <- UCEC.Normalized

################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  UCEC.rnaseq,
  extract.cols = c("PKMYT1","ERN2","ALPK3","LMTK3","PHKG2","CAMKV","CDC42BPG","ADCK5","POMK","STK31","ALPK2")) %>%
  rename(cohort = dataset,
         PKMYT1='PKMYT1',ERN2='ERN2',ALPK3='ALPK3',LMTK3='LMTK3',PHKG2='PHKG2',CAMKV='CAMKV',
         CDC42BPG='CDC42BPG',ADCK5='ADCK5',POMK='POMK',STK31='STK31',ALPK2='ALPK2') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("PKMYT1","ERN2","ALPK3","LMTK3","PHKG2","CAMKV","CDC42BPG","ADCK5","POMK","STK31","ALPK2")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  UCEC.surv_rnaseq.cut <- surv_cutpoint(
    UCEC.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(UCEC.surv_rnaseq.cut)
  
  distri <- plot(UCEC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-UCEC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  UCEC.surv_rnaseq.cat <- surv_categorize(UCEC.surv_rnaseq.cut)
  
  
  
fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = UCEC.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-UCEC/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  pval <- surv_pvalue(fit)
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}




###############

print("##################
#Survival ANALYSIS OF TCGA-KIRC GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-KIRC

########################################
#KIRC
KIRC.clin <- data.frame(patient.bcr_patient_barcode = KIRC.clinical$gdc_cases.submitter_id,
                        patient.vital_status = KIRC.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = KIRC.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = KIRC.clinical$xml_days_to_death)


survivalTCGA(KIRC.clin) -> KIRC.clinical.surv

####################################
#get the normaliized gene expression
#####################################

KIRC.Normalized <- data.frame(t(KIRC.Normalize.expression))
KIRC.Normalized <- cbind(bcr_patient_barcode = rownames(KIRC.Normalized), KIRC.Normalized)
rownames(KIRC.Normalized) <- NULL
KIRC.rnaseq <- KIRC.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  KIRC.rnaseq,
  extract.cols = c("CAMKK1","CDK18","PKMYT1","PNCK","ITPKA","TSSK3","MAPK15","CAMK1D","SRPK3","PLK5","RPS6KL1","ALPK2")) %>%
  rename(cohort = dataset,
         CAMKK1= 'CAMKK1',CDK18='CDK18',PKMYT1='PKMYT1',PNCK='PNCK',ITPKA='ITPKA',TSSK3='TSSK3',
         MAPK15='MAPK15',CAMK1D='CAMK1D',SRPK3='SRPK3',PLK5='PLK5',RPS6KL1='RPS6KL1',ALPK2='ALPK2') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("CAMKK1","CDK18","PKMYT1","PNCK","ITPKA","TSSK3","MAPK15","CAMK1D","SRPK3","PLK5","RPS6KL1","ALPK2")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  KIRC.surv_rnaseq.cut <- surv_cutpoint(
    KIRC.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(KIRC.surv_rnaseq.cut)
  
  distri <- plot(KIRC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  KIRC.surv_rnaseq.cat <- surv_categorize(KIRC.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = KIRC.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  pval <- surv_pvalue(fit)
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  print(paste("completed", geneid[i]))
  
}

#################

print("##################
#Survival ANALYSIS OF TCGA-PCPG GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-PCPG

########################################

#PCPG
PCPG.clin <- data.frame(patient.bcr_patient_barcode = PCPG.clinical$gdc_cases.submitter_id,
                        patient.vital_status = PCPG.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = PCPG.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = PCPG.clinical$xml_days_to_death)


survivalTCGA(PCPG.clin) -> PCPG.clinical.surv
####################################
#get the normaliized gene expression
#####################################

PCPG.Normalized <- data.frame(t(PCPG.Normalize.expression))
PCPG.Normalized <- cbind(bcr_patient_barcode = rownames(PCPG.Normalized), PCPG.Normalized)
rownames(PCPG.Normalized) <- NULL
PCPG.rnaseq <- PCPG.Normalized

################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  PCPG.rnaseq,
  extract.cols = c("PAK3","LMTK3","TTBK1","CAMKV","RPS6KL1")) %>%
  rename(cohort = dataset,
         PAK3= 'PAK3',LMTK3='LMTK3',TTBK1='TTBK1',CAMKV='CAMKV',RPS6KL1='RPS6KL1') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid =c("PAK3","LMTK3","TTBK1","CAMKV","RPS6KL1")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  PCPG.surv_rnaseq.cut <- surv_cutpoint(
    PCPG.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(PCPG.surv_rnaseq.cut)
  
  distri <- plot(PCPG.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PCPG/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  PCPG.surv_rnaseq.cat <- surv_categorize(PCPG.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = PCPG.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-PCPG/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  pval <- surv_pvalue(fit)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}



#################

print("##################
#Survival ANALYSIS OF TCGA-LUSC GENES
##############################")

########################################

#CLINICAL DATA OF TCGA-LUSC

########################################
#LUSC
LUSC.clin <- data.frame(patient.bcr_patient_barcode = LUSC.clinical$gdc_cases.submitter_id,
                        patient.vital_status = LUSC.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = LUSC.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = LUSC.clinical$xml_days_to_death)


survivalTCGA(LUSC.clin) -> LUSC.clinical.surv
####################################
#get the normaliized gene expression
#####################################

LUSC.Normalized <- data.frame(t(LUSC.Normalize.expression))
LUSC.Normalized <- cbind(bcr_patient_barcode = rownames(LUSC.Normalized), LUSC.Normalized)
rownames(LUSC.Normalized) <- NULL
LUSC.rnaseq <- LUSC.Normalized

################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  LUSC.rnaseq,
  extract.cols = c("PHKA1","RIOK1","PKMYT1","MAP3K10","PNCK","ITPKA","TTBK1","CAMKV","WNK2","SGK494","ADCK5","TSSK6",
                   "SRPK3","POMK","STK31","ALPK2","STK38L")) %>%
  rename(cohort = dataset,
         PHKA1='PHKA1',RIOK1='RIOK1',PKMYT1='PKMYT1',MAP3K10='MAP3K10',PNCK='PNCK',ITPKA='ITPKA',
         TTBK1='TTBK1',CAMKV='CAMKV',WNK2='WNK2',SGK494='SGK494',ADCK5='ADCK5',TSSK6='TSSK6',
         SRPK3='SRPK3',POMK='POMK',STK31='STK31',ALPK2='ALPK2',STK38L='STK38L') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("PHKA1","RIOK1","PKMYT1","MAP3K10","PNCK","ITPKA","TTBK1","CAMKV","WNK2","SGK494","ADCK5","TSSK6",
           "SRPK3","POMK","STK31","ALPK2","STK38L")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  LUSC.surv_rnaseq.cut <- surv_cutpoint(
    LUSC.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(LUSC.surv_rnaseq.cut)
  
  distri <- plot(LUSC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUSC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  LUSC.surv_rnaseq.cat <- surv_categorize(LUSC.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = LUSC.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUSC/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  pval <- surv_pvalue(fit)
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}



###################
print("##################
#Survival ANALYSIS OF TCGA-LUAD GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-LUAD

########################################


#LUAD
LUAD.clin <- data.frame(patient.bcr_patient_barcode = LUAD.clinical$gdc_cases.submitter_id,
                        patient.vital_status = LUAD.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = LUAD.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = LUAD.clinical$xml_days_to_death)


survivalTCGA(LUAD.clin) -> LUAD.clinical.surv

####################################
#get the normaliized gene expression
#####################################

LUAD.Normalized <- data.frame(t(LUAD.Normalize.expression))
LUAD.Normalized <- cbind(bcr_patient_barcode = rownames(LUAD.Normalized), LUAD.Normalized)
rownames(LUAD.Normalized) <- NULL
LUAD.rnaseq <- LUAD.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  LUAD.rnaseq,
  extract.cols = c("PHKA1","PAK3","DYRK2","PKMYT1","PNCK","ERN2","ITPKA","TTBK1",
                   "STK32B","WNK2","SGK494","STK32A","ADCK5","SRPK3","POMK","STK31","RPS6KL1","ALPK2")) %>%
  rename(cohort = dataset,
         PHKA1='PHKA1',PAK3='PAK3',DYRK2='DYRK2',PKMYT1='PKMYT1',PNCK='PNCK',ERN2='ERN2',
         ITPKA='ITPKA',TTBK1='TTBK1',STK32B='STK32B',WNK2='WNK2',SGK494='SGK494',STK32A='STK32A',
         ADCK5='ADCK5',SRPK3='SRPK3',POMK='POMK',STK31='STK31',RPS6KL1='RPS6KL1',ALPK2='ALPK2') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid =c("PHKA1","PAK3","DYRK2","PKMYT1","PNCK","ERN2","ITPKA","TTBK1",
          "STK32B","WNK2","SGK494","STK32A","ADCK5","SRPK3","POMK","STK31","RPS6KL1","ALPK2")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  LUAD.surv_rnaseq.cut <- surv_cutpoint(
    LUAD.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(LUAD.surv_rnaseq.cut)
  
  distri <- plot(LUAD.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUAD/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  LUAD.surv_rnaseq.cat <- surv_categorize(LUAD.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = LUAD.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-LUAD/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  pval <- surv_pvalue(fit)
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}


################


###################

print("##################
#Survival ANALYSIS OF TCGA-STAD GENES
##############################")

########################################

#CLINICAL DATA OF TCGA-STAD

########################################
#STAD
STAD.clin <- data.frame(patient.bcr_patient_barcode = STAD.clinical$gdc_cases.submitter_id,
                        patient.vital_status = STAD.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = STAD.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = STAD.clinical$xml_days_to_death)


survivalTCGA(STAD.clin) -> STAD.clinical.surv
####################################
#get the normaliized gene expression
#####################################

STAD.Normalized <- data.frame(t(STAD.Normalize.expression))
STAD.Normalized <- cbind(bcr_patient_barcode = rownames(STAD.Normalized), STAD.Normalized)
rownames(STAD.Normalized) <- NULL
STAD.rnaseq <- STAD.Normalized

################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  STAD.rnaseq,
  extract.cols = c("STK3","PKMYT1","CAMKV","STK32C","SGK494","MAPK15","POMK","STK31",
                   "NEK5","RPS6KL1","SGK223")) %>%
  rename(cohort = dataset,
         STK3='STK3',PKMYT1='PKMYT1',CAMKV='CAMKV',STK32C='STK32C',SGK494='SGK494',MAPK15='MAPK15',
         POMK='POMK',STK31='STK31',NEK5='NEK5',RPS6KL1='RPS6KL1',SGK223='SGK223') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("STK3","PKMYT1","CAMKV","STK32C","SGK494","MAPK15","POMK","STK31",
           "NEK5","RPS6KL1","SGK223")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  STAD.surv_rnaseq.cut <- surv_cutpoint(
    STAD.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(STAD.surv_rnaseq.cut)
  
  distri <- plot(STAD.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-STAD/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  STAD.surv_rnaseq.cat <- surv_categorize(STAD.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = STAD.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-STAD/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  pval <- surv_pvalue(fit)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}




#####################
print("##################
#Survival ANALYSIS OF TCGA-GBM GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-GBM

########################################

#GBM
GBM.clin <- data.frame(patient.bcr_patient_barcode = GBM.clinical$gdc_cases.submitter_id,
                       patient.vital_status = GBM.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = GBM.clinical$cgc_case_days_to_last_follow_up, 
                       patient.days_to_death = GBM.clinical$xml_days_to_death)


survivalTCGA(GBM.clin) -> GBM.clinical.surv
####################################
#get the normaliized gene expression
#####################################

GBM.Normalized <- data.frame(t(GBM.Normalize.expression))
GBM.Normalized <- cbind(bcr_patient_barcode = rownames(GBM.Normalized), GBM.Normalized)
rownames(GBM.Normalized) <- NULL
GBM.rnaseq <- GBM.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  GBM.rnaseq,
  extract.cols = c("PHKA1","STK33","ERN2","STK36","STK17A","STK32A","STK38L")) %>%
  rename(cohort = dataset,
         PHKA1='PHKA1',STK33='STK33',ERN2='ERN2',STK36='STK36',STK17A='STK17A',STK32A='STK32A',STK38L='STK38L') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("PHKA1","STK33","ERN2","STK36","STK17A","STK32A","STK38L")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  GBM.surv_rnaseq.cut <- surv_cutpoint(
    GBM.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(GBM.surv_rnaseq.cut)
  
  distri <- plot(GBM.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-GBM/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  GBM.surv_rnaseq.cat <- surv_categorize(GBM.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = GBM.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-GBM/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  pval <- surv_pvalue(fit)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}


##############
print("##################
#Survival ANALYSIS OF TCGA-THCA GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-THCA

########################################

#THCA
THCA.clin <- data.frame(patient.bcr_patient_barcode = THCA.clinical$gdc_cases.submitter_id,
                        patient.vital_status = THCA.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = THCA.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = THCA.clinical$xml_days_to_death)


survivalTCGA(THCA.clin) -> THCA.clinical.surv

####################################
#get the normaliized gene expression
#####################################

THCA.Normalized <- data.frame(t(THCA.Normalize.expression))
THCA.Normalized <- cbind(bcr_patient_barcode = rownames(THCA.Normalized), THCA.Normalized)
rownames(THCA.Normalized) <- NULL
THCA.rnaseq <- THCA.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  THCA.rnaseq,
  extract.cols = c("CAMK1G","PKMYT1","ERN2","ITPKA","STK32A","CDC42BPG")) %>%
  rename(cohort = dataset,
         CAMK1G='CAMK1G',PKMYT1='PKMYT1',ERN2='ERN2',ITPKA='ITPKA',STK32A='STK32A',CDC42BPG='CDC42BPG') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("CAMK1G","PKMYT1","ERN2","ITPKA","STK32A","CDC42BPG")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  THCA.surv_rnaseq.cut <- surv_cutpoint(
    THCA.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(THCA.surv_rnaseq.cut)
  
  distri <- plot(THCA.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-THCA/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  THCA.surv_rnaseq.cat <- surv_categorize(THCA.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = THCA.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-THCA/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  pval <- surv_pvalue(fit)
  
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}




#######################

print("##################
#Survival ANALYSIS OF TCGA-CESC GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-CESC

########################################


#CESC
CESC.clin <- data.frame(patient.bcr_patient_barcode = CESC.clinical$gdc_cases.submitter_id,
                        patient.vital_status = CESC.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = CESC.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = CESC.clinical$xml_days_to_death)


survivalTCGA(CESC.clin) -> CESC.clinical.surv

####################################
#get the normaliized gene expression
#####################################

CESC.Normalized <- data.frame(t(CESC.Normalize.expression))
CESC.Normalized <- cbind(bcr_patient_barcode = rownames(CESC.Normalized), CESC.Normalized)
rownames(CESC.Normalized) <- NULL
CESC.rnaseq <- CESC.Normalized

################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  CESC.rnaseq,
  extract.cols = c("PKMYT1","LMTK3","PHKG2","CDC42BPG","ADCK5","SGK223")) %>%
  rename(cohort = dataset,
         PKMYT1='PKMYT1',LMTK3='LMTK3',PHKG2='PHKG2',CDC42BPG='CDC42BPG',ADCK5='ADCK5',SGK223='SGK223')%>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("PKMYT1","LMTK3","PHKG2","CDC42BPG","ADCK5","SGK223")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  CESC.surv_rnaseq.cut <- surv_cutpoint(
    CESC.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(CESC.surv_rnaseq.cut)
  
  distri <- plot(CESC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
 # png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CESC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  CESC.surv_rnaseq.cat <- surv_categorize(CESC.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = CESC.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CESC/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  pval <- surv_pvalue(fit)
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  print(paste("completed", geneid[i]))
  
}



##################
print("##################
#Survival ANALYSIS OF TCGA-COAD GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-COAD

########################################


#COAD
COAD.clin <- data.frame(patient.bcr_patient_barcode = COAD.clinical$gdc_cases.submitter_id,
                        patient.vital_status = COAD.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = COAD.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = COAD.clinical$xml_days_to_death)


survivalTCGA(COAD.clin) -> COAD.clinical.surv
####################################
#get the normaliized gene expression
#####################################

COAD.Normalized <- data.frame(t(COAD.Normalize.expression))
COAD.Normalized <- cbind(bcr_patient_barcode = rownames(COAD.Normalized), COAD.Normalized)
rownames(COAD.Normalized) <- NULL
COAD.rnaseq <- COAD.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  COAD.rnaseq,
  extract.cols = c("PHKA1","STK3","PKMYT1","ALPK3","LMTK3","CAMKV","STK32C","SGK494",
                   "MAPK15","POMK","PLK5","STK31","NEK5","RPS6KL1")) %>%
  rename(cohort = dataset,
         PHKA1='PHKA1',STK3='STK3',PKMYT1='PKMYT1',ALPK3='ALPK3',LMTK3='LMTK3',CAMKV='CAMKV',STK32C='STK32C',
         SGK494='SGK494',MAPK15='MAPK15',POMK='POMK',PLK5='PLK5',STK31='STK31',NEK5='NEK5',RPS6KL1='RPS6KL1') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("PHKA1","STK3","PKMYT1","ALPK3","LMTK3","CAMKV","STK32C","SGK494",
           "MAPK15","POMK","PLK5","STK31","NEK5","RPS6KL1")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  COAD.surv_rnaseq.cut <- surv_cutpoint(
    COAD.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(COAD.surv_rnaseq.cut)
  
  distri <- plot(COAD.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-COAD/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  COAD.surv_rnaseq.cat <- surv_categorize(COAD.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = COAD.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-COAD/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  
  
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  pval <- surv_pvalue(fit)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}

##################
print("##################
#Survival ANALYSIS OF TCGA-HNSC GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-HNSC

########################################

#HNSC
HNSC.clin <- data.frame(patient.bcr_patient_barcode = HNSC.clinical$gdc_cases.submitter_id,
                        patient.vital_status = HNSC.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = HNSC.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = HNSC.clinical$xml_days_to_death)


survivalTCGA(HNSC.clin) -> HNSC.clinical.surv

####################################
#get the normaliized gene expression
#####################################
HNSC.Normalized <- data.frame(t(HNSC.Normalize.expression))
HNSC.Normalized <- cbind(bcr_patient_barcode = rownames(HNSC.Normalized), HNSC.Normalized)
rownames(HNSC.Normalized) <- NULL
HNSC.rnaseq <- HNSC.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  HNSC.rnaseq,
  extract.cols = c("CAMK1G","STK3","CDK18","PKMYT1","PNCK","ITPKA","TTBK1","TSSK3","CAMKV",
                   "STK17A","SGK494","ADCK5","TSSK6","STK31")) %>%
  rename(cohort = dataset,
         CAMK1G='CAMK1G',STK3='STK3',CDK18='CDK18',PKMYT1='PKMYT1',PNCK='PNCK',ITPKA='ITPKA',TTBK1='TTBK1',
         TSSK3='TSSK3',CAMKV='CAMKV',STK17A='STK17A',SGK494='SGK494',ADCK5='ADCK5',TSSK6='TSSK6',STK31='STK31') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("CAMK1G","STK3","CDK18","PKMYT1","PNCK","ITPKA","TTBK1","TSSK3","CAMKV",
           "STK17A","SGK494","ADCK5","TSSK6","STK31")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  HNSC.surv_rnaseq.cut <- surv_cutpoint(
    HNSC.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(HNSC.surv_rnaseq.cut)
  
  distri <- plot(HNSC.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-HNSC/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  HNSC.surv_rnaseq.cat <- surv_categorize(HNSC.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = HNSC.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-HNSC/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  pval <- surv_pvalue(fit)
  
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}



##################
print("##################
#Survival ANALYSIS OF TCGA-KICH GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-KICH

########################################


#KICH
KICH.clin <- data.frame(patient.bcr_patient_barcode = KICH.clinical$gdc_cases.submitter_id,
                        patient.vital_status = KICH.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = KICH.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = KICH.clinical$xml_days_to_death)


survivalTCGA(KICH.clin) -> KICH.clinical.surv

####################################
#get the normaliized gene expression
#####################################

KICH.Normalized <- data.frame(t(KICH.Normalize.expression))
KICH.Normalized <- cbind(bcr_patient_barcode = rownames(KICH.Normalized), KICH.Normalized)
rownames(KICH.Normalized) <- NULL
KICH.rnaseq <- KICH.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  KICH.rnaseq,
  extract.cols = c("PHKA1","PKMYT1","PNCK","ALPK3","ITPKA","STK32A","CAMK1D","RPS6KL1")) %>%
  rename(cohort = dataset,
         PHKA1='PHKA1',PKMYT1='PKMYT1',PNCK='PNCK',ALPK3='ALPK3',ITPKA='ITPKA',STK32A='STK32A',
         CAMK1D='CAMK1D',RPS6KL1='RPS6KL1') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("PHKA1","PKMYT1","PNCK","ALPK3","ITPKA","STK32A","CAMK1D","RPS6KL1")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  KICH.surv_rnaseq.cut <- surv_cutpoint(
    KICH.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(KICH.surv_rnaseq.cut)
  
  distri <- plot(KICH.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KICH/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  KICH.surv_rnaseq.cat <- surv_categorize(KICH.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = KICH.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KICH/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  pval <- surv_pvalue(fit)
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}






##################
print("##################
#Survival ANALYSIS OF TCGA-ESCA GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-ESCA

########################################

#ESCA
ESCA.clin <- data.frame(patient.bcr_patient_barcode = ESCA.clinical$gdc_cases.submitter_id,
                        patient.vital_status = ESCA.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = ESCA.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = ESCA.clinical$xml_days_to_death)


survivalTCGA(ESCA.clin) -> ESCA.clinical.surv

####################################
#get the normaliized gene expression
#####################################

ESCA.Normalized <- data.frame(t(ESCA.Normalize.expression))
ESCA.Normalized <- cbind(bcr_patient_barcode = rownames(ESCA.Normalized), ESCA.Normalized)
rownames(ESCA.Normalized) <- NULL
ESCA.rnaseq <- ESCA.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  ESCA.rnaseq,
  extract.cols = c("CAMKK1","CAMK1G","STK3","RIOK1","PKMYT1","STK32C","ADCK5",
                   "MAPK15","STK31","RPS6KL1","SGK223")) %>%
  rename(cohort = dataset,
         CAMKK1= 'CAMKK1',CAMK1G='CAMK1G',STK3='STK3',RIOK1='RIOK1',PKMYT1='PKMYT1',STK32C='STK32C',
         ADCK5='ADCK5', MAPK15='MAPK15',STK31='STK31',RPS6KL1='RPS6KL1',SGK223='SGK223') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("CAMKK1","CAMK1G","STK3","RIOK1","PKMYT1","STK32C","ADCK5",
           "MAPK15","STK31","RPS6KL1","SGK223")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  ESCA.surv_rnaseq.cut <- surv_cutpoint(
    ESCA.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(ESCA.surv_rnaseq.cut)
  
  distri <- plot(ESCA.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-ESCA/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  ESCA.surv_rnaseq.cat <- surv_categorize(ESCA.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = ESCA.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-ESCA/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  pval <- surv_pvalue(fit)
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}





##################
print("##################
#Survival ANALYSIS OF TCGA-KIRP GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-KIRP

########################################

#KIRP
KIRP.clin <- data.frame(patient.bcr_patient_barcode = KIRP.clinical$gdc_cases.submitter_id,
                        patient.vital_status = KIRP.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = KIRP.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = KIRP.clinical$xml_days_to_death)


survivalTCGA(KIRP.clin) -> KIRP.clinical.surv

####################################
#get the normaliized gene expression
#####################################

KIRP.Normalized <- data.frame(t(KIRP.Normalize.expression))
KIRP.Normalized <- cbind(bcr_patient_barcode = rownames(KIRP.Normalized), KIRP.Normalized)
rownames(KIRP.Normalized) <- NULL
KIRP.rnaseq <- KIRP.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  KIRP.rnaseq,
  extract.cols = c("CAMKK1","DYRK2","PKMYT1","PNCK","ERN2","ITPKA","TSSK3",
                   "CAMKV","WNK2","MAPK15","CAMK1D","PLK5")) %>%
  rename(cohort = dataset,
         CAMKK1= 'CAMKK1',DYRK2='DYRK2',PKMYT1='PKMYT1',PNCK='PNCK',ERN2='ERN2',ITPKA='ITPKA',TSSK3='TSSK3',
         CAMKV='CAMKV',WNK2='WNK2',MAPK15='MAPK15',CAMK1D='CAMK1D',PLK5='PLK5') %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("CAMKK1","DYRK2","PKMYT1","PNCK","ERN2","ITPKA","TSSK3",
           "CAMKV","WNK2","MAPK15","CAMK1D","PLK5")


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  KIRP.surv_rnaseq.cut <- surv_cutpoint(
    KIRP.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(KIRP.surv_rnaseq.cut)
  
  distri <- plot(KIRP.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRP/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  KIRP.surv_rnaseq.cat <- surv_categorize(KIRP.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = KIRP.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRP/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  
  pval <- surv_pvalue(fit)
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}





##################
print("##################
#Survival ANALYSIS OF TCGA-CHOL GENES
      ##############################")

########################################

#CLINICAL DATA OF TCGA-CHOL

########################################

#CHOL
CHOL.clin <- data.frame(patient.bcr_patient_barcode = CHOL.clinical$gdc_cases.submitter_id,
                        patient.vital_status = CHOL.clinical$gdc_cases.diagnoses.vital_status, patient.days_to_last_followup = CHOL.clinical$cgc_case_days_to_last_follow_up, 
                        patient.days_to_death = CHOL.clinical$xml_days_to_death)


survivalTCGA(CHOL.clin) -> CHOL.clinical.surv

####################################
#get the normaliized gene expression
#####################################
CHOL.Normalized <- data.frame(t(CHOL.Normalize.expression))
CHOL.Normalized <- cbind(bcr_patient_barcode = rownames(CHOL.Normalized), CHOL.Normalized)
rownames(CHOL.Normalized) <- NULL
CHOL.rnaseq <- CHOL.Normalized


################
# EXTRACT EXPRESSION DATA OF GENES 
###############

expressionsTCGA(
  CHOL.rnaseq,
  extract.cols = c("CAMKK1","CAMK1G","PHKA1","PAK3","STK3","CDK18","RIOK1","DYRK2","PKMYT1" ,"STK33","MAP3K10","PNCK",
                   "ERN2","ALPK3","ITPKA","LMTK3","STK32B","PHKG2","TSSK3","STK36","CAMKV","STK17A","WNK2","STK32C",
                   "SGK494","STK32A","CDC42BPG","ADCK5","TSSK6","MAPK15","CAMK1D","SRPK3","POMK","PLK5","STK31","NEK5",    
                   "RPS6KL1" , "STK38L","SGK223" )) %>%
  rename(cohort = dataset,
         CAMKK1='CAMKK1',CAMK1G='CAMK1G',PHKA1='PHKA1', PAK3='PAK3',STK3='STK3',CDK18='CDK18',RIOK1='RIOK1',
         DYRK2='DYRK2',PKMYT1='PKMYT1' ,STK33='STK33',MAP3K10='MAP3K10',PNCK='PNCK',
         ERN2='ERN2',ALPK3='ALPK3',ITPKA='ITPKA',LMTK3='LMTK3',STK32B='STK32B',PHKG2='PHKG2',
         TSSK3='TSSK3',STK36='STK36',CAMKV='CAMKV',STK17A='STK17A',WNK2='WNK2',STK32C='STK32C',
         SGK494='SGK494',STK32A='STK32A',CDC42BPG='CDC42BPG',ADCK5='ADCK5',TSSK6='TSSK6',MAPK15='MAPK15',
         CAMK1D='CAMK1D',SRPK3='SRPK3',POMK='POMK',PLK5='PLK5',STK31='STK31',NEK5='NEK5',    
         RPS6KL1='RPS6KL1' , STK38L='STK38L',SGK223='SGK223' ) %>%
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


######################
#GENES for survival analysis
#######################

library(survminer)

geneid = c("CAMKK1","CAMK1G","PHKA1","PAK3","STK3","CDK18","RIOK1","DYRK2","PKMYT1" ,"STK33","MAP3K10","PNCK",
           "ERN2","ALPK3","ITPKA","LMTK3","STK32B","PHKG2","TSSK3","STK36","CAMKV","STK17A","WNK2","STK32C",
           "SGK494","STK32A","CDC42BPG","ADCK5","TSSK6","MAPK15","CAMK1D","SRPK3","POMK","PLK5","STK31","NEK5",    
           "RPS6KL1" , "STK38L","SGK223" )


######################
#Plot SURVIVAL CURVES FOR GENES
#######################

for(i in 1: length(geneid)){
  
  CHOL.surv_rnaseq.cut <- surv_cutpoint(
    CHOL.surv_rnaseq,
    time = "times",
    event = "patient.vital_status",
    variables = c(paste(geneid[i]), "cohort")
  )
  summary(CHOL.surv_rnaseq.cut)
  
  distri <- plot(CHOL.surv_rnaseq.cut, geneid[i], palette = "npg")
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CHOL/survival/", geneid[i], ".distri.png", sep=""), width = 1000, height = 800)
  #print(distri)
  
  
  #dev.off()
  
  CHOL.surv_rnaseq.cat <- surv_categorize(CHOL.surv_rnaseq.cut)
  
  
  
  fit <- survfit(as.formula(paste0("Surv(times, patient.vital_status) ~" , geneid[i], "+ cohort")),
                 data = CHOL.surv_rnaseq.cat)
  
  #png(paste0("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-CHOL/survival/", geneid[i], ".survival.jpg", sep=""))
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
  #print(gg)
  #dev.off()
  
  pval <- surv_pvalue(fit)
  surmed <- surv_median(fit)
  surmedian <- rbind(surmedian, surmed)
  dataSur <- gg$data.survtable
  colnames(dataSur) <- c("strata", "time", "n.risk", "pct.risk", "n.event", "cum.n.event",
                         "n.censor","cum.n.censor", "strata_size", "gene", "cohort")
  dataSurvTable <- rbind(dataSurvTable, dataSur)
  dataplot <- gg$data.survplot
  colnames(dataplot) <- c("time", "n.risk", "n.event n.censor", "surv","std.err","upper","lower",
                          "strata", "gene",  "cohort")
  dataSurvplot <- rbind(dataSurvplot, dataplot) 
  pvalue <- rbind(pvalue,  pval)
  high = fit$strata[1]
  low <- fit$strata[2]
  fitt1 <- cbind(high, low)
  fitt2 <- cbind(fitt1, pval)
  fittable <- rbind(fittable,fitt2)
  
  print(paste("completed", geneid[i]))
  
}

##############
#split the column into two

#########

aa <- separate(data = surmedian, col = strata, into = c("gene", "express", "cohart", "cancer" ))

aa%>% 
  group_by(gene,cancer) %>%
  arrange(desc(median)) %>%
  slice(1:2) -> aa1



#extract odd rows with higest median

aa1 %>% dplyr::filter(row_number() %% 2 == 1) -> high.med.aa1

#split the colum in fittable

fittable1 <- fittable
fittable1$strata <- rownames(fittable1)

fittable2 <- separate(data = fittable1, col = strata, into = c("gene", "express", "cohart", "cancer" ))

#merge the 2 data frames based on gene and cancer
data <- merge(high.med.aa1, fittable2, by=c("gene","cancer"))

data1 <- data[,c(1:3, 5,13)]

#create a boxplot 

#data2 <- melt(data1)
########remove missing data median has NA value
newdata <- na.omit(data1)
newdata1 <- newdata[!duplicated(newdata),]


png("/projects/ccs/schurerlab/Rimpi/TCGARecount2/boxplot.survival1.jpg", width = 8, height = 5, units = 'in', res=300)
 p <- ggplot(newdata1, aes(x= fct_infreq(gene), y=pval, color=express.x)) + coord_flip()
 p1 <-p + scale_color_manual(values = c("#2E9FDF", "red"))  +
    geom_jitter(position=position_jitter(0.2)) +
   theme(axis.text.x = element_text(face="bold", color="black", 
                                    size=5, angle=45),  
         axis.text.y = element_text(face="bold", color="black", 
                                    size=5)) + geom_hline(yintercept=0.05, colour="black",linetype="dashed")
 p1 
print(p1)

dev.off()


########################################
#Normalize the expression of cancer that does not have normal samples (LGG, SARC, DLBC, THYM, OV, SKCM, LAML, UVM, PAAD, TGCT, ACC, UCS) 
########################################
#TCGA-DLBC, TCGA-LAML
CanID <- c("TCGA-LGG","TCGA-ACC","TCGA-SARC","TCGA-THYM","TCGA-OV","TCGA-UVM","TCGA-SKCM","TCGA-PAAD","TCGA-TGCT","TCGA-UCS")

for(i in 8:length(CanID)){
print(CanID[i])
exp <- rse_gene[, rse_gene$cgc_file_investigation %in% c(paste(CanID[4]))]

###remove duplicated  TCGA barcodes from  
clic.exp <- data.frame(colData(exp))
clic.exp <- clic.exp[ , ! apply( clic.exp , 2 , function(x) all(is.na(x)) ) ]

clinical <- clic.exp


clinical1 <- data.frame(barcode = clinical$gdc_cases.samples.portions.analytes.aliquots.submitter_id, caseId = rownames(clinical), bigfile = clinical$bigwig_file, tissuetype = clinical$cgc_file_investigation,
                        project = clinical$project, clinical$reads_downloaded,
                        experimentalStrategy = clinical$cgc_file_experimental_strategy,
                        pairedend = clinical$paired_end, mappedreads= clinical$mapped_read_count,
                        auc = clinical$auc, filesize= clinical$gdc_file_size ,
                        access =clinical$gdc_access, platform = clinical$gdc_platform, state = clinical$gdc_state, category = clinical$gdc_data_category,
                        center = clinical$gdc_center.short_name,centertype= clinical$gdc_center.center_type, gender = clinical$gdc_cases.demographic.gender,
                        birthyear = clinical$gdc_cases.demographic.year_of_birth, race = clinical$gdc_cases.demographic.race,
                        deathyear = clinical$gdc_cases.demographic.year_of_death, tissue = clinical$gdc_cases.project.primary_site,
                        stage = clinical$gdc_cases.diagnoses.tumor_stage, vitalstatus= clinical$gdc_cases.diagnoses.vital_status,
                        sampletype= clinical$gdc_cases.samples.sample_type, lastfollowupdays= clinical$cgc_case_days_to_last_follow_up,
                        followupID = clinical$cgc_follow_up_id,patientid= clinical$xml_bcr_patient_barcode,drugtherapy= clinical$cgc_drug_therapy_pharmaceutical_therapy_type,
                        site= clinical$xml_icd_o_3_site, followuptumorstatus= clinical$cgc_follow_up_tumor_status,
                        Daytstobirth= clinical$xml_days_to_birth, icdhistology= clinical$xml_icd_o_3_histology, tumorstatus= clinical$cgc_case_tumor_status,
                        ageofdiagnosis= clinical$cgc_case_age_at_diagnosis, tissuesourcecite= clinical$cgc_sample_tissue_source_site_code, samplecreationdatetime = clinical$gdc_cases.samples.portions.creation_datetime,
                        tissueoforigin= clinical$gdc_cases.diagnoses.tissue_or_organ_of_origin, diagnosisdaystobirth = clinical$gdc_cases.diagnoses.days_to_birth,
                        projectaccession = clinical$gdc_cases.project.program.dbgap_accession_number,
                         tissuesite= clinical$gdc_cases.tissue_source_site.project, daystodeath = clinical$xml_days_to_death)
#Shortnametype = clinical$Shortnametype,

clinical2 <- clinical1[!duplicated(clinical1$barcode),]

###remove NAs from 
exp.genes <- data.frame(rowData(exp))
table(is.na(exp.genes$symbol))



exp.genes.sort <- exp.genes[which(apply(!is.na(exp.genes), 1, all)),]
table(is.na(exp.genes.sort$symbol))

dim(exp.genes.sort[duplicated(exp.genes.sort$symbol),])
exp.genes.sort.dup <- exp.genes.sort[!duplicated(exp.genes.sort$symbol),]






#breastCancer1[duplicated(breastCancer1$barcode),]
#breastCancer3 <- breastCancer2[!duplicated(breastCancer2$barcode),]
#breastCancer2 <- breastCancer[!duplicated(breastCancer$gdc_cases.samples.portions.analytes.aliquots.submitter_id),]
#breastCancer3 <- breastCancer[!duplicated(breastCancer$gdc_cases.submitter_id),]
table(clinical2$sampletype)


####get the expression counts
exp.counts <- assay(exp)

####Remove the expression counts id not present in genesN
exp.counts.sort <- exp.counts[which(rownames(exp.counts) %in% exp.genes.sort.dup$gene_id),]


####Remove the expression counts barcodes not present in breastCancer2
#which(colnames(expressioncounts1) == "137FDD43-DC3C-4027-ADED-DA9D69CF909C")
#which(colnames(expressioncounts1) == "1E347C7A-FBE8-408D-8014-EE33869116E0")

rm <- setdiff(colnames(exp.counts.sort),clinical2$caseId)
expr.count.fin <- exp.counts.sort[ , !(colnames(exp.counts.sort) %in% rm)]


#add rownames as gene symbol
rownames(expr.count.fin) <- exp.genes.sort.dup$symbol
colnames(expr.count.fin) <- clinical2$barcode

#
rownames(clinical2) <- clinical2$barcode

rownames(exp.genes.sort.dup) <- exp.genes.sort.dup$symbol

dir.create(file.path(filePath, outputDir, paste(CanID[i])))
setwd(file.path(filePath, outputDir, paste(CanID[i])))


saveRDS(clinical, file=paste(getwd(),'/', CanID[i],"_clinical.rds", sep = ""))
write.table(clinical, file=paste(getwd(),'/', CanID[i],"_clinical.csv", sep="" ), sep=",")
#saveRDS(breastCancer2, file = paste(getwd(),'/',"_clinic.Rdata", sep="" ))

saveRDS(expr.count.fin, file = paste(getwd(),'/', "_expression.Rdata", sep="" ))
write.table(expr.count.fin, paste(getwd(),'/', CanID[i], "_expression.csv", sep="" ), sep=",")


#rownames(clinical1) <- rownames(clinical)
# convert into SummarizedExperiment
data <- SummarizedExperiment(assays = list(counts=expr.count.fin),
                             rowData= exp.genes.sort.dup, colData=clinical2)
data1 <- as(data, "RangedSummarizedExperiment")
saveRDS(data1, file = paste(getwd(),'/',CanID[i],".rds", sep="" ))




#DE Analysis
# dataTP <- as.character(clinical2$barcode[which(clinical2$Shortnametype == "TP")])
# dataNT <- as.character(clinical2$barcode[which(clinical2$Shortnametype == "NT")])

#dataPrep <-TCGAanalyze_Preprocess(object = exp, cor.cut = 0.6)


#TCGAanalyze_Preprocenamesss
#dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                      geneInfo = geneInfo,
#                                    method = "gcContent")   


#normalize read counts
Normdd <- scale_counts(data1)
Normddval <- as.matrix(assays(Normdd)$counts)


dataFilt <- TCGAanalyze_Filtering(tabDF = Normddval,
                                  method = "quantile", 
                                  qnt.cut =  0.25) 






split <- unlist(strsplit(CanID[i], "[-]"))

assign(paste(split[2],".clinical", sep =""), clinical) 
assign(paste(split[2],".Normalize.expression", sep =""), Normddval)

assign(paste(split[2],".Normalize.expression_filter", sep =""), dataFilt)

write.csv(clinical, file= paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/", split[2], ".clinical.csv", sep=""))
write.csv(Normddval, file= paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/", split[2], ".Normalize.expression.csv", sep=""))
write.csv(dataFilt, file= paste("/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/", split[2], ".Normalize.expression_filtered.csv", sep=""))
}                  
save(LGG.clinical, LGG.Normalize.expression, LGG.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LGG_Normalize_filter.RData")
save(ACC.clinical, ACC.Normalize.expression, ACC.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/ACC_Normalize_filter.RData")
save(SARC.clinical, SARC.Normalize.expression, SARC.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/SARC_Normalize_filter.RData")
save(DLBC.clinical, DLBC.Normalize.expression, DLBC.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/DLBC_Normalize_filter.RData")

save(THYM.clinical, THYM.Normalize.expression, THYM.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/THYM_Normalize_filter.RData")
save(OV.clinical, OV.Normalize.expression, OV.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/OV_Normalize_filter.RData")
save(LAML.clinical, LAML.Normalize.expression, LAML.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/LAML_Normalize_filter.RData")
save(UVM.clinical, UVM.Normalize.expression, UVM.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/UVM_Normalize_filter.RData")
save(SKCM.clinical, SKCM.Normalize.expression, SKCM.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/SKCM_Normalize_filter.RData")
save(PAAD.clinical, PAAD.Normalize.expression, PAAD.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/PAAD_Normalize_filter.RData")
save(TGCT.clinical, TGCT.Normalize.expression, TGCT.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/TGCT_Normalize_filter.RData")
save(UCS.clinical, UCS.Normalize.expression, UCS.Normalize.expression_filter, file = "/projects/ccs/schurerlab/Rimpi/TCGARecount2/For_Derek_analysis/UCS_Normalize_filter.RData")




########################################
#get the DE  (LGG, SARC, DLBC, THYM, OV, SKCM, LAML, UVM, PAAD, TGCT, ACC, UCS) 
########################################


########################################

# 3rd analaysis WGCNA with DE genes in each cancer

#####################################



######################
#WGCNA
#####################

dim(dataDEGs)

#####################
# Remove non differentially-expressed genes
##################

# Filter out genes from the dataFilt file from TCGAbiolinks (the datFilt file is Genes with an average expression level <0.25 in all )
 dataFilt
Normddval
data1 <- readRDS("/projects/ccs/schurerlab/Rimpi/TCGARecount2/TCGA-KIRC/TCGA-KIRC.rds")
clic <- colData(data1)

#preprocess raw data
dataPrep <- TCGAanalyze_Preprocess(object = data1, cor.cut = 0.6)

#dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
#                                      geneInfo = geneInfo,
#                                    method = "gcContent")   

#normalize read counts
Normdd <- scale_counts(data1)
Normddval <- as.matrix(assays(Normdd)$counts)

#data filter according to the 
dataFilt <- TCGAanalyze_Filtering(tabDF = Normddval,
                                  method = "quantile", 
                                  qnt.cut =  0.25) 


#filter genes using  CV > 0.5 from genefilter package 
library(genefilter)
ffun <- filterfun(cv(0.5))
data_filt_cv_rem <- genefilter(dataFilt,ffun)


#Remove the auxiliary data and transpose the expression data for further analysis
datExpr0 = as.data.frame(t(dataFilt));
#datExpr0 <- dataFilt;

http://pklab.med.harvard.edu/scw2014/WGCNA.html

#grouping in high and NT and TP (trait)
trait <- as.character(clic$Shortnametype)
trait[trait == 'NT'] = 0
trait[trait == 'TP'] = 1

#calculate gene significance measure for lymphocyte score (lscore) - Welch's t-Test
GS_lscore = t(sapply(1:ncol(WGCNA_matrix2),function(x)c(t.test(WGCNA_matrix2[,x]~trait,var.equal=F)$p.value,
                          t.test(WGCNA_matrix2[,x]~trait,var.equal=F)$estimate[1],
                          t.test(WGCNA_matrix2[,x]~trait,var.equal=F)$estimate[2])))
GS_lscore = cbind(GS_lscore, abs(GS_lscore[,2] - GS_lscore[,3]))
colnames(GS_lscore) = c('p_value','mean_high_lscore','mean_low_lscore',
                        'effect_size(high-low score)'); rownames(GS_lscore) = colnames(WGCNA_matrix2)


#calculate module significance
MS.lscore = as.data.frame(cbind(GS_lscore,modules))
MS.lscore$log_p_value = -log10(as.numeric(MS.lscore$p_value))
MS.lscore = ddply(MS.lscore, .(modules), summarize, mean(log_p_value), sd(log_p_value))
colnames(MS.lscore) = c('modules','pval','sd')
MS.lscore.bar = as.numeric(MS.lscore[,2])
MS.lscore.bar[MS.lscore.bar<(-log10(0.05))] = 0
names(MS.lscore.bar) = GO.module.name

METree.GO = METree
label.order = match(METree$labels,paste0('ME',labels2colors(0:(length(unique(modules))-1))))
METree.GO$labels = GO.module.name[label.order]
plotTree.wBars(as.phylo(METree.GO), MS.lscore.bar, tip.labels = TRUE, scale = 0.2)







log_counts <- log2(data_filt_cv_rem + 1)

log_counts <- log_counts[apply(log_counts, 1, var) > 0,]

x = melt(as.matrix(log_counts))

colnames(x) = c('gene_id', 'sample', 'value')
ggplot(x, aes(x=value, color=sample)) + geom_density()
RNAseq_voom = voom(data_filt_cv_rem)$E






cordist <- function(dat) {
 cor_matrix  <- cor(t(dat))

dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
dist_matrix <- log1p(dist_matrix)
dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
}
sim_matrix <- cordist(log_counts)

# Construct adjacency matrix
adj_matrix <- adjacency.fromSimilarity(sim_matrix, power=12, type='signed')

# Delete similarity matrix to free up memory
rm(sim_matrix)
gc()


# Convert to matrix
gene_ids <- rownames(adj_matrix)

adj_matrix <- matrix(adj_matrix, nrow=nrow(adj_matrix))
rownames(adj_matrix) <- gene_ids
colnames(adj_matrix) <- gene_ids

heatmap.2(t(adj_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Adjacency matrix',
          density.info='none', revC=TRUE)


# Cluster gene expression profiles; the flashClust function from
# the authors of WGCNA is another options for larger datasets.
# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
gene_tree <- hclust(as.dist(1 - adj_matrix), method="average")

# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=15,
                                   deepSplit=TRUE)

# assign a color to each module for easier visualization and referencing
module_colors <- labels2colors(module_labels)


















##########################################
#geom_boxplot(color="black") +
######################
#genes below pval < 0.05

sig.genes <- subset(newdata1, newdata1$pval < 0.05) #39



################
p <- boxplot(pval ~ gene, data=data1)
p + scale_y_discrete(name = "pval", limits= c(0, 1)) 
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)


boxplot(pval ~ gene, data=data1, lwd = 2, ylab = 'pval')
stripchart(pval ~ gene, vertical = TRUE, data = data1, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

ggboxplot(data1, x = "gene", y = "pval",
          palette = "jco", color = "express.x",
          add = "jitter")


ggplot(data1, aes(gene, pval, fill=express.x)) + geom_boxplot()

ggplot(data1, aes(x=gene, y=pval, fill=express.x)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center')
##################









group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))

####################Survival rate
dataSurvIDG <- TCGAanalyze_SurvivalK(clinical_patient = breastCancer2,
                                   dataGE = dataFilt,
                                   Genelist = rownames(commomIDG),
                                   Survresult = FALSE,
                                   ThreshTop = 0.67,
                                   ThreshDown = 0.33,
                                   p.cut = 0.05, group1, group2)

tabSurvKMcomplete <- rbind(dataFilt,dataSurvIDG)


dataSurv <- TCGAanalyze_SurvivalK(clinical_patient = breastCancer2,
                                  dataGE = dataFilt,
                                  Genelist = rownames(dataDEGs),
                                  Survresult = FALSE,
                                  ThreshTop = 0.67,
                                  ThreshDown = 0.33,
                                  p.cut = 0.05, group1, group2)

################################














##################################################################

# get the CommonIDG data with each breast cancer
clinical <- breastCancer[!duplicated(breastCancer$gdc_cases.samples.portions.analytes.aliquots.submitter_id),]
clinical$IDs <- toupper(clinical$gdc_cases.samples.portions.analytes.aliquots.submitter_id)
sum(clinical$IDs %in% colnames(dataFilt))


# get the columns that contain data we can use: days to death, new tumor event, last day contact to....
ind_keep <- grep('follow_up',colnames(clinical))

# this is a bit tedious, since there are numerous follow ups, let's collapse them together and keep the first value (the higher one) if more than one is available
new_tum <- as.matrix(clinical[,ind_keep])
new_tum <- as.matrix(clinical[,ind_keep])
new_tum_collapsed <- c()
for (i in 1:dim(new_tum)[1]){
  if ( sum ( is.na(new_tum[i,])) < dim(new_tum)[2]){
    m <- min(new_tum[i,],na.rm=T)
    new_tum_collapsed <- c(new_tum_collapsed,m)
  } else {
    new_tum_collapsed <- c(new_tum_collapsed,'NA')
  }
}

# do the same to death
ind_keep <- grep('days_to_death',colnames(clinical))
death <- as.matrix(clinical[,ind_keep])
death_collapsed <- c()
for (i in 1:dim(death)[1]){
  if ( sum ( is.na(death[i,])) < dim(death)[2]){
    m <- max(death[i,],na.rm=T)
    death_collapsed <- c(death_collapsed,m)
  } else {
    death_collapsed <- c(death_collapsed,'NA')
  }
}


# and days last follow up here we take the most recent which is the max number
ind_keep <- grep('days_to_last_followup',colnames(clinical))
fl <- as.matrix(clinical[,ind_keep])
fl_collapsed <- c()
for (i in 1:dim(fl)[1]){
  if ( sum (is.na(fl[i,])) < dim(fl)[2]){
    m <- max(fl[i,],na.rm=T)
    fl_collapsed <- c(fl_collapsed,m)
  } else {
    fl_collapsed <- c(fl_collapsed,'NA')
  }
}

# and put everything together
all_clin <- data.frame(new_tum_collapsed,death_collapsed,fl_collapsed)
colnames(all_clin) <- c('new_tumor_days', 'death_days', 'followUp_days')

head(all_clin)




######################

table(duplicated(breastCancer$gdc_cases.samples.portions.analytes.aliquots.submitter_id))



gdc_cases.samples.portions.analytes.aliquots.submitter_id

rowData(brca.rec$tcga_breast)$symbol


gdc_cases.project.state  "legancy"


data <- readGeneExpressionQuantification(files = files, 
                                         cases = query$results[[1]]$cases, summarizedExperiment = summarizedExperiment, 
                                         experimental.strategy = unique(query$results[[1]]$experimental_strategy))

# breastCancer1 <- data.frame(barcode = breastCancer$gdc_cases.samples.portions.analytes.aliquots.submitter_id, caseId = rownames(breastCancer), bigfile = breastCancer$bigwig_file, tissuetype = breastCancer$cgc_file_investigation,
#                               project = breastCancer$project, breastCancer$reads_downloaded,
#                               experimentalStrategy = breastCancer$cgc_file_experimental_strategy,
#                               pairedend = breastCancer$paired_end, mappedreads= breastCancer$mapped_read_count,
#                               auc = breastCancer$auc, filesize= breastCancer$gdc_file_size ,
#                               access =breastCancer$gdc_access, platform = breastCancer$gdc_platform, state = breastCancer$gdc_state, category = breastCancer$gdc_data_category,
#                               center = breastCancer$gdc_center.short_name,centertype= breastCancer$gdc_center.center_type, gender = breastCancer$gdc_cases.demographic.gender,
#                               birthyear = breastCancer$gdc_cases.demographic.year_of_birth, race = breastCancer$gdc_cases.demographic.race,
#                               deathyear = breastCancer$gdc_cases.demographic.year_of_death, tissue = breastCancer$gdc_cases.project.primary_site,
#                               stage = breastCancer$gdc_cases.diagnoses.tumor_stage, vitalstatus= breastCancer$gdc_cases.diagnoses.vital_status,
#                               sampletype= breastCancer$gdc_cases.samples.sample_type, lastfollowupdays= breastCancer$cgc_case_days_to_last_follow_up,
#                             followupID = breastCancer$cgc_follow_up_id,patientid= breastCancer$xml_bcr_patient_barcode,drugtherapy= breastCancer$cgc_drug_therapy_pharmaceutical_therapy_type,
#                             site= breastCancer$xml_icd_o_3_site, followuptumorstatus= breastCancer$cgc_follow_up_tumor_status,
#                             Daytstobirth= breastCancer$xml_days_to_birth, icdhistology= breastCancer$xml_icd_o_3_histology, tumorstatus= breastCancer$cgc_case_tumor_status,
#                             ageofdiagnosis= breastCancer$cgc_case_age_at_diagnosis, tissuesourcecite= breastCancer$cgc_sample_tissue_source_site_code, samplecreationdatetime = breastCancer$gdc_cases.samples.portions.creation_datetime,
#                             tissueoforigin= breastCancer$gdc_cases.diagnoses.tissue_or_organ_of_origin, diagnosisdaystobirth = breastCancer$gdc_cases.diagnoses.days_to_birth,
#                             projectaccession = breastCancer$gdc_cases.project.program.dbgap_accession_number,
#                             Shortnametype = breastCancer$Shortnametype, tissuesite= breastCancer$gdc_cases.tissue_source_site.project, daystodeath = breastCancer$xml_days_to_death)

#ER.status= breastCancer$xml_breast_carcinoma_estrogen_receptor_status, 
#
#Tissues in recount2
#tissues <- c("adipose tissue", "adrenal", "gland", "bladder", 
"blood", "blood vessel", "bone marrow", "brain", "breast", 
"cervix uteri", "colon", "esophagus", "fallopian tube", 
"heart", "kidney", "liver", "lung", "muscle", "nerve", 
"ovary", "pancreas", "pituitary", "prostate", "salivary", 
"gland", "skin", "small intestine", "spleen", "stomach", 
"testis", "thyroid", "uterus", "vagina")

{
  isServeOK()
  if (missing(breastCancer)) 
    stop("Please set query parameter")
  if (any(duplicated(breastCancer$gdc_cases.samples.portions.analytes.aliquots.submitter_id) & breastCancer$data.type != 
      "Clinical data" & breastCancer$data.type != "Protein expression quantification") {
    dup <- breastCancer$cases[duplicated(breastCancer$gdc_cases.samples.portions.analytes.aliquots.submitter_id)]
    dup <- breastCancer[breastCancer$gdc_cases.samples.portions.analytes.aliquots.submitter_id %in% 
                                dup, c("tags", "gdc_cases.samples.portions.analytes.aliquots.submitter_id", "gdc_experimental_strategy")]
    dup <- dup[order(dup$gdc_cases.samples.portions.analytes.aliquots.submitter_id), ]
    print(knitr::kable(dup))
    stop("There are samples duplicated. We will not be able to preapre it")
  }