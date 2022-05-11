# ANALYSIS FOR SIMON: PROGNOSIS VALUES OF GENE SIGNATURES DERIEVED FROM IMMUNE EVADING CLONES
# 20201114
# Sunny Wu
# 
# qrsh -pe smp 2 -l mem_requested=10G
# source activate Renv
# Rdev
# 
# GENE SIGNATURES SUMMARY
## 6 files
### IE1_UP; up-reg in IE1 vs bulk
### IE1_DOWN; down-reg in IE2 vs bulk
### IE2_UP; up-reg in IE2 vs bulk
### IE2_DOWN; down-reg in IE1 vs bulk
### common_UP; common up-reg beween IE1 & IE2
### common_down; common down reg between IE1 & IE2
#
# 01: SETUP R-ENVIRONMENT -------------------------------------------------

library(survival)
library(survminer)
library(RColorBrewer)
library(org.Hs.eg.db)
library(dplyr)

ncol1 = brewer.pal(12,"Paired")
ncol2 = brewer.pal(8,"Dark2")

#heatmap
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#heatmap-as-raster-image

# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)


# 02: DIRECTORIES -------------------------------------------------------------

dir.create("/share/ScratchGeneral/sunwu/projects/immune_evasion_prognosis/run04_same_gene_list")
setwd("/share/ScratchGeneral/sunwu/projects/immune_evasion_prognosis/run04_same_gene_list")

# 03: LOAD COHORT EXPRESISON DATA ---------------------------------------------------------------

# METABRIC
temp_dir_new <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/12_survival_analysis/cohort_data_heloisa/"
## S1, read the expression table
load(paste0(temp_dir_new,"METABRIC_Discovery.RData")) # METABRIC dis
expr_METABRIC_dis <- Discovery
rm(Discovery)

load(paste0(temp_dir_new,"METABRIC_Validation.RData")) # METABRIC val
expr_METABRIC_val <- Validation
rm(Validation)

# TCGA
temp_dir <- "/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/12_survival_analysis/cohort_data/"
expr_TCGA <- readRDS(paste0(temp_dir,"TCGA.BRCA.RNAseq.RPKM.rds"))           # TCGA

# GSE21653 and GSE58812
expr_GSE21653 = as.matrix(read.delim(paste0(temp_dir,"GSE21653.self_subtract")))     # GEO BRCA dataset
expr_GSE58812 = as.matrix(read.delim(paste0(temp_dir,"GSE58812.self_subtract")))     # GEO BRCA dataset

rownames(expr_GSE21653) = as.vector(unlist(mget(rownames(expr_GSE21653), envir=org.Hs.egSYMBOL, ifnotfound=NA)))
rownames(expr_GSE58812) = as.vector(unlist(mget(rownames(expr_GSE58812), envir=org.Hs.egSYMBOL, ifnotfound=NA)))


# 04: LOAD COHORT CLINICAL DATA -----------------------------------------------

# METABRIC
clin_METABRIC = read.delim("/share/ScratchGeneral/sunwu/projects/MINI_ATLAS_PROJECT/Sept2019/15_CIBERSORTx_decon_prognosis//brca_metabric_clinical_data.tsv")
clin_METABRIC$Overall.Survival.Status <- gsub(":LIVING","",clin_METABRIC$Overall.Survival.Status)
clin_METABRIC$Overall.Survival.Status <- gsub(":DECEASED","",clin_METABRIC$Overall.Survival.Status)
clin_METABRIC$Overall.Survival.Status <- as.numeric(clin_METABRIC$Overall.Survival.Status)
clin_METABRIC$Overall.Survival..Months. <- as.numeric(clin_METABRIC$Overall.Survival..Months.)

clin_METABRIC_dis <- clin_METABRIC[clin_METABRIC$Patient.ID %in% as.vector(colnames(expr_METABRIC_dis)),,drop=F]
clin_METABRIC_val <- clin_METABRIC[clin_METABRIC$Patient.ID %in% as.vector(colnames(expr_METABRIC_val)),,drop=F]

# discovery cohort
expr_METABRIC_dis <- expr_METABRIC_dis[,colnames(expr_METABRIC_dis) %in% as.vector(clin_METABRIC_dis$Patient.ID),drop=F]
rownames(clin_METABRIC_dis) <- clin_METABRIC_dis$Patient.ID
clin_METABRIC_dis <- clin_METABRIC_dis[colnames(expr_METABRIC_dis),,drop=F]
all.equal(colnames(expr_METABRIC_dis),
          as.vector(clin_METABRIC_dis$Patient.ID))

# validation cohort
expr_METABRIC_val <- expr_METABRIC_val[,colnames(expr_METABRIC_val) %in% as.vector(clin_METABRIC_val$Patient.ID),drop=F]
rownames(clin_METABRIC_val) <- clin_METABRIC_val$Patient.ID
clin_METABRIC_val <- clin_METABRIC_val[colnames(expr_METABRIC_val),,drop=F]
all.equal(colnames(expr_METABRIC_val),
          as.vector(clin_METABRIC_val$Patient.ID))

# combined cohort
temp_intersect <- intersect(rownames(expr_METABRIC_val), rownames(expr_METABRIC_dis))
expr_METABRIC_val2 <- expr_METABRIC_val[rownames(expr_METABRIC_val) %in% temp_intersect,,drop=F]
expr_METABRIC_dis2 <- expr_METABRIC_dis[rownames(expr_METABRIC_dis) %in% temp_intersect,,drop=F]
all.equal(rownames(expr_METABRIC_val2),
          rownames(expr_METABRIC_dis2)
          )
expr_METABRIC <- cbind(expr_METABRIC_val2, expr_METABRIC_dis2)


clin_METABRIC <- clin_METABRIC[clin_METABRIC$Patient.ID %in% colnames(expr_METABRIC),,drop=F]
rownames(clin_METABRIC) <- clin_METABRIC$Patient.ID

# TCGA
load(paste0(temp_dir,"TCGAclin.RData"))
clin_TCGA = clin[which(clin$cancer=="BRCA"),];rownames(clin_TCGA) = clin_TCGA[,1]

# GSE21653 and GSE58812
clin_GSE21653 = read.delim(paste0(temp_dir,"GSE21653.clinical"))
clin_GSE58812 = read.delim(paste0(temp_dir,"GSE58812.OS_TN"))

determine_subtype = function(clin)
{
  clin$subtype = NA
  clin[which(clin$LumA==1),"subtype"] = "LumA"
  clin[which(clin$LumB==1),"subtype"] = "LumB"
  clin[which(clin$Her2==1),"subtype"] = "Her2"
  clin[which(clin$Basal==1),"subtype"] = "Basal"
  return(clin)
}
clin_GSE21653 = determine_subtype(clin_GSE21653)
clin_GSE58812$subtype = "Basal"


# 05: FILTER TCGA EXPRESSION DATA CORRECTLY ---------------------------------------------

# remove substring
head(colnames(expr_TCGA))
colnames(expr_TCGA) <- substr(colnames(expr_TCGA),1,12)
head(colnames(expr_TCGA))
length(colnames(expr_TCGA)[duplicated(colnames(expr_TCGA))])

# temp_duplicated <- colnames(expr_TCGA)[duplicated(colnames(expr_TCGA))]
# temp <- expr_TCGA[,colnames(expr_TCGA) %in% temp_duplicated,drop=F]
# temp <- temp[,sort(colnames(temp)),drop=F]

expr_TCGA <- expr_TCGA[, !duplicated(colnames(expr_TCGA))]
clin_TCGA <- clin_TCGA[clin_TCGA$patient %in% colnames(expr_TCGA),,drop=F]
clin_TCGA <- clin_TCGA[colnames(expr_TCGA),,drop=F]
print(all.equal(rownames(clin_TCGA),
                colnames(expr_TCGA)))

# 06: LOAD GENE SIGNATURES ----------------------------------------------------

# METABRIC 20058 genes
# TCGA 18863 genes
intersect_genes_TCGA_METABRIC <- intersect(rownames(expr_METABRIC),
                                           rownames(expr_TCGA))

# lengt intersect
# > length(intersect_genes_TCGA_METABRIC)
# [1] 16842


temp_files <- list.files("../gene_signatures")
# for(cutoff in c(50,100,200)){
  temp_genesig_df <- data.frame()
  for(file in temp_files){
    print(file)
    temp_csv <- read.csv(paste0("../gene_signatures/",file))
    # temp_csv <- temp_csv[1:cutoff,,drop=F]
    temp_csv <- temp_csv[,"H.Ensembl",drop=F]
    temp_name <- gsub(".csv","",file)
    colnames(temp_csv) <- "gene"
    temp_csv <- temp_csv[temp_csv$gene %in% intersect_genes_TCGA_METABRIC,,drop=F]
    rownames(temp_csv) <- c(1:nrow(temp_csv))
    
    
    
    
    temp_csv$cluster <- temp_name
    
    
    print(temp_name)
    
    
    if(nrow(temp_genesig_df)==0){
      temp_genesig_df <- temp_csv
    } else {
      temp_genesig_df <- rbind(temp_genesig_df, temp_csv)
    }
  }
  # n <- paste0("temp_genesig_df_",cutoff)
  # assign(n, temp_genesig_df)
  # rm(temp_genesig_df)
# }

# 07: GENERATE AVERAGE EXPRESSION VALUES --------------------------------------

# inmatrix <- expr_METABRIC_dis
# gene <- temp_genesig_df
# top <- 50

averaged_gene_signature <- function(inmatrix, gene, top = temp_number, cohort)
{
  outmatrix = NULL
  name = NULL
  
  temp_cluster_ids <- as.vector(unique(gene$cluster))
  
  # filter for shared genes between platforms
  gene <- gene[gene$gene %in% rownames(inmatrix),,drop=F]
  temp_df_combined <- data.frame(row.names = c(1:top),
                                 gene = c(1:top))
  for(i in temp_cluster_ids)
  {
    temp_length_intersect <- length(intersect(rownames(inmatrix), 
                                              as.character(na.omit(gene[gene$cluster==i,]$gene[1:top]))))
    # print(i)
    # print(temp_length_intersect)
    
    outmatrix = cbind(outmatrix, 
                      apply(inmatrix[intersect(rownames(inmatrix), 
                                               as.character(na.omit(gene[gene$cluster==i,]$gene[1:top]))),],
                            2,
                            mean))
    name = c(name, gsub("-","",gsub("\\+","",gsub("\\ ","",as.character(i)))))
    
    temp_df <- data.frame(cluster = intersect(rownames(inmatrix), 
                                              as.character(na.omit(gene[gene$cluster==i,]$gene[1:top]))))
    colnames(temp_df) <- i
    
    if(temp_length_intersect >= top){
      temp_df_combined <- cbind(temp_df_combined, temp_df) 
    } else {
      print(paste0(i," has only ",temp_length_intersect," genes from intersect with cohort matrix"))
      temp_df_combined <- rowr::cbind.fill(temp_df_combined, temp_df, fill=NA)
    }
    
    # print number of genes overlapped
    # print(i)
    # print(length(as.character(na.omit(gene[gene$cluster==i,]$gene[1:top]))))
    # print(length(intersect(rownames(inmatrix), 
    #                        as.character(na.omit(gene[gene$cluster==i,]$gene[1:top]))))/
    #         length(as.character(na.omit(gene[gene$cluster==i,]$gene[1:top])))
    # )
  }
  # write.csv(temp_df_combined, 
  #           paste0("genes_used_for_prognosis/temp_df_combined_",cohort,".csv"))
  
  colnames(outmatrix) = name
  return(outmatrix)
}

temp_cutoffs <- c(25)
for(cutoff in temp_cutoffs){
  print(cutoff)
  expr_METABRIC_mean <- averaged_gene_signature(expr_METABRIC, temp_genesig_df, top=cutoff)
  expr_METABRIC_mean = data.frame(clin_METABRIC[rownames(expr_METABRIC_mean),], expr_METABRIC_mean)
  
  colnames(expr_METABRIC_mean) <- gsub("Overall.Survival..Months.","OS",colnames(expr_METABRIC_mean))
  colnames(expr_METABRIC_mean) <- gsub("Overall.Survival.Status","EVENT",colnames(expr_METABRIC_mean))
  colnames(expr_METABRIC_mean) <- gsub("Age.at.Diagnosis","AGE",colnames(expr_METABRIC_mean))
  colnames(expr_METABRIC_mean) <- gsub("Pam50...Claudin.low.subtype","PAM50SUBTYPE",colnames(expr_METABRIC_mean))
  
  expr_TCGA_mean <-  averaged_gene_signature(expr_TCGA, temp_genesig_df, top=cutoff)
  expr_TCGA_mean = data.frame(clin_TCGA[substr(rownames(expr_TCGA_mean),1,12),], expr_TCGA_mean)
  colnames(expr_TCGA_mean) <- gsub("subtype","PAM50SUBTYPE",colnames(expr_TCGA_mean))
  colnames(expr_TCGA_mean) <- gsub("age","AGE",colnames(expr_TCGA_mean))
  
  expr_GSE21653_mean <-  averaged_gene_signature(expr_GSE21653, temp_genesig_df, top=cutoff)
  expr_GSE21653_mean = data.frame(clin_GSE21653[rownames(expr_GSE21653_mean),], expr_GSE21653_mean); colnames(expr_GSE21653_mean)[1] <- "OS"
  colnames(expr_GSE21653_mean) <- gsub("subtype","PAM50SUBTYPE",colnames(expr_GSE21653_mean))
  colnames(expr_GSE21653_mean) <- gsub("Age","AGE",colnames(expr_GSE21653_mean))
  
  expr_GSE58812_mean <-  averaged_gene_signature(expr_GSE58812, temp_genesig_df, top=cutoff)
  expr_GSE58812_mean = data.frame(clin_GSE58812[rownames(expr_GSE58812_mean),], expr_GSE58812_mean); colnames(expr_GSE58812_mean)[2] <- "EVENT"
  colnames(expr_GSE58812_mean) <- gsub("subtype","PAM50SUBTYPE",colnames(expr_GSE58812_mean))
  colnames(expr_GSE58812_mean) <- gsub("Age","AGE",colnames(expr_GSE58812_mean))
  
  # expr_METABRIC_dis_mean <- averaged_gene_signature(expr_METABRIC_dis, gene, top=cutoff)
  # expr_METABRIC_dis_mean = data.frame(clin_METABRIC_dis[rownames(expr_METABRIC_dis_mean),], expr_METABRIC_dis_mean)
  # 
  # expr_METABRIC_val_mean <- averaged_gene_signature(expr_METABRIC_val, gene, top=cutoff)
  # expr_METABRIC_val_mean = data.frame(clin_METABRIC_val[rownames(expr_METABRIC_val_mean),], expr_METABRIC_val_mean)
  
  n <- paste0("expr_METABRIC_mean_",cutoff)
  assign(n, expr_METABRIC_mean)
  
  n <- paste0("expr_TCGA_mean_",cutoff)
  assign(n, expr_TCGA_mean)
  
  n <- paste0("expr_GSE21653_mean_",cutoff)
  assign(n, expr_GSE21653_mean)
  
  n <- paste0("expr_GSE58812_mean_",cutoff)
  assign(n, expr_GSE58812_mean)
  
  # n <- paste0("expr_METABRIC_dis_mean_",cutoff)
  # assign(n, expr_METABRIC_dis_mean)
  # n <- paste0("expr_METABRIC_val_mean_",cutoff)
  # assign(n, expr_METABRIC_val_mean)
  
  rm(expr_METABRIC_mean, expr_TCGA_mean,expr_GSE21653_mean,expr_GSE58812_mean)
}

# 08: RUN SURVIVAL 30vs30 WITH YEAR CUTOFFS ---------------------------------------------------------

temp_celltypes_list <- colnames(expr_METABRIC_mean_25)[36:41]

for(cutoff in temp_cutoffs){
  dir.create(paste0("Clinical_Mean_with_year_cutoffs_",cutoff))
}

temp_percent_cutoff <- 0.3
temp_df_combined <- data.frame()
for(cutoff in temp_cutoffs){
  print(cutoff)
  for(cohort in c("TCGA", "METABRIC")){
    print(cohort)
    inmatrix_cohort <- get(paste0("expr_",cohort,"_mean_",cutoff))
    inmatrix_cohort$PAM50SUBTYPE <- as.character(inmatrix_cohort$PAM50SUBTYPE)
    
    if(length(inmatrix_cohort$PAM50SUBTYPE[is.na(inmatrix_cohort$PAM50SUBTYPE)])>0){
      inmatrix_cohort$PAM50SUBTYPE[is.na(inmatrix_cohort$PAM50SUBTYPE)] <- "no call"
    }
    
    # year cutoffs
    if(cohort == "METABRIC"){
      inmatrix_cohort <- inmatrix_cohort[inmatrix_cohort$OS < 241,,drop=F] # limit METABRIC to 20 years
    }
    if(cohort == "TCGA"){
      inmatrix_cohort <- inmatrix_cohort[inmatrix_cohort$OS < 5476,,drop=F] # limit TCGA to 15 years
    }
    
    inmatrix_cohort <- inmatrix_cohort[!grepl('NA', rownames(inmatrix_cohort)), ]
    
    if(!cohort == "GSE58812"){
      temp_subtypes <- c("Basal", "Her2", "LumA", "LumB", "All")
    } else {
      temp_subtypes <- c("Basal")
    }
    
    for(subtype in temp_subtypes){
      # print(subtype)
      if(!subtype == "All"){
        inmatrix <- inmatrix_cohort[inmatrix_cohort$PAM50SUBTYPE == subtype,,drop=F]
      } else {
        inmatrix <- inmatrix_cohort
      }
      
      # inmatrix <- inmatrix[!rownames(inmatrix) %in%,,drop=F]
      
      for(celltype in "common_up_GeneList_HO_collapsed"){
        # print(celltype)
        # M1, stratify by top 30% vs bot 30%
        
        inmatrix_subset = inmatrix[order(inmatrix[,celltype], decreasing=T),]
        
        temp_top = (round(nrow(inmatrix_subset)*temp_percent_cutoff))
        temp_bot = nrow(inmatrix_subset)-temp_top+1
        
        inmatrix_top = inmatrix_subset[1:temp_top,];inmatrix_top$CELLTYPE = "top30%"
        inmatrix_bottom = inmatrix_subset[temp_bot:nrow(inmatrix_subset),];inmatrix_bottom$CELLTYPE = "bot30%"
        
        inmatrix_subset = rbind(inmatrix_top,inmatrix_bottom)
        colnames(inmatrix_subset) <- toupper(colnames(inmatrix_subset))
        surv <- survfit(Surv(OS, EVENT) ~ CELLTYPE, data = inmatrix_subset)
        
        # make p-values plotted (which is from a log-rank test) = to p-value derived from coxprop model
        coxph.surv=Surv(inmatrix_subset[,'OS'],inmatrix_subset[,'EVENT'])
        temp_summary <- as.data.frame(summary(coxph(coxph.surv~inmatrix_subset[,'CELLTYPE']))$coefficients[1,c(2,5)])
        
        # temp_summary_save <- as.data.frame(summary(coxph(coxph.surv~inmatrix_subset[,'CELLTYPE']))$coefficients)
        temp_summary <- t(temp_summary)
        colnames(temp_summary) <- c("HR", "pval")
        rownames(temp_summary) <- NULL
        temp_summary <- as.data.frame(temp_summary)
        temp_summary$signature <- celltype
        temp_summary$cohort <- cohort
        temp_summary$subtype <- subtype
        temp_summary$cutoff <- cutoff
        
        # log rank p-value
        temp_summary_logrank <- surv_pvalue(surv)
        temp_summary_logrank <- temp_summary_logrank[,"pval"]
        temp_summary$logrankpval <- temp_summary_logrank
        
        if(nrow(temp_df_combined)==0){
          temp_df_combined <- temp_summary
        } else {
          temp_df_combined <- rbind(temp_df_combined, temp_summary)
        }
        
        # colnames(temp_summary_save) <- paste0("METABRIC_",cohort)
        # write.csv(temp_summary_save,
        #           paste0("Clinical_Mean_",cutoff,"/Coxprop_stats_01_",cohort,"_KMcomparison.csv"))
        
        temp_plot <- ggsurvplot(surv,
                                data = inmatrix_subset,
                                title = paste(celltype,":", cohort, " ",subtype),
                                font.title = c(8, "bold", "darkblue"),
                                # palette = ncol1,
                                pval=T,
                                # pval = round(temp_summary[2,],digits=8),
                                legend="top",
                                # legend.labs = c("bottom 20%","top 20%"),
                                risk.table = T)

        png(file = paste0("Clinical_Mean_with_year_cutoffs_",cutoff,"/",cohort,"_01_",celltype,"_",subtype,"_KMcomparison.png"),
            width = 6,
            height = 6,
            res=600,
            units="in")
        print(temp_plot)
        dev.off()

        pdf(file = paste0("Clinical_Mean_with_year_cutoffs_",cutoff,"/",cohort,"_01_",celltype,"_",subtype,"_KMcomparison.pdf"),
            width = 6,
            height =6,
            onefile = F)
        print(temp_plot)
        dev.off()
        
        write.csv(inmatrix_subset,
                  paste0("Clinical_Mean_with_year_cutoffs_",cutoff,"/",cohort,"_01_",celltype,"_",subtype,"_clinical_data.csv"))
        
      }
      
    }
  }
}

write.csv(temp_df_combined,
          "Clinical_Mean_with_year_cutoffs_25/pvalues_30_vs_30_yearcutoff.csv")






# 09: RUN SURVIVAL 30vs30 WITH MID & YEAR CUTOFFS ---------------------------------------------------------

temp_celltypes_list <- colnames(expr_METABRIC_mean_25)[36:41]

for(cutoff in temp_cutoffs){
  dir.create(paste0("Clinical_Mean_with_mid_year_cutoffs_",cutoff))
}

temp_percent_cutoff <- 0.3
temp_df_combined <- data.frame()
for(cutoff in temp_cutoffs){
  print(cutoff)
  for(cohort in c("TCGA", "METABRIC")){
    print(cohort)
    inmatrix_cohort <- get(paste0("expr_",cohort,"_mean_",cutoff))
    inmatrix_cohort$PAM50SUBTYPE <- as.character(inmatrix_cohort$PAM50SUBTYPE)

    if(length(inmatrix_cohort$PAM50SUBTYPE[is.na(inmatrix_cohort$PAM50SUBTYPE)])>0){
      inmatrix_cohort$PAM50SUBTYPE[is.na(inmatrix_cohort$PAM50SUBTYPE)] <- "no call"
    }

    inmatrix_cohort <- inmatrix_cohort[!is.na(inmatrix_cohort$OS),,drop=F]

    # year cutoffs
    if(cohort == "METABRIC"){
      inmatrix_cohort <- inmatrix_cohort[inmatrix_cohort$OS < 241,,drop=F] # limit METABRIC to 20 years
    }
    if(cohort == "TCGA"){
      inmatrix_cohort <- inmatrix_cohort[inmatrix_cohort$OS < 5476,,drop=F] # limit TCGA to 15 years
    }

    inmatrix_cohort <- inmatrix_cohort[!grepl('NA', rownames(inmatrix_cohort)), ]

    if(!cohort == "GSE58812"){
      temp_subtypes <- c("Basal", "Her2", "LumA", "LumB", "All")
    } else {
      temp_subtypes <- c("Basal")
    }

    for(subtype in temp_subtypes){
      # print(subtype)
      if(subtype == "All"){
        inmatrix <- inmatrix_cohort
        } else {
          inmatrix <- inmatrix_cohort[inmatrix_cohort$PAM50SUBTYPE == subtype,,drop=F]
      }

      # inmatrix <- inmatrix[!rownames(inmatrix) %in%,,drop=F]

      for(celltype in temp_celltypes_list){
        # print(celltype)
        # M1, stratify by top 30% vs bot 30%

        inmatrix_subset = inmatrix[order(inmatrix[,celltype], decreasing=T),]

        temp_top = (round(nrow(inmatrix_subset)*temp_percent_cutoff))
        temp_bot = nrow(inmatrix_subset)-temp_top+1

        inmatrix_top = inmatrix_subset[1:temp_top,];inmatrix_top$CELLTYPE = "top30%"
        inmatrix_bottom = inmatrix_subset[temp_bot:nrow(inmatrix_subset),];inmatrix_bottom$CELLTYPE = "bot30%"

        inmatrix_mid <- inmatrix_subset[!rownames(inmatrix_subset) %in% c(rownames(inmatrix_top),rownames(inmatrix_bottom)),,drop=F]
        inmatrix_mid$CELLTYPE = "mid40%"

        inmatrix_subset = rbind(inmatrix_top,inmatrix_mid,inmatrix_bottom)
        colnames(inmatrix_subset) <- toupper(colnames(inmatrix_subset))
        surv <- survfit(Surv(OS, EVENT) ~ CELLTYPE, data = inmatrix_subset)

        # make p-values plotted (which is from a log-rank test) = to p-value derived from coxprop model
        coxph.surv=Surv(inmatrix_subset[,'OS'],inmatrix_subset[,'EVENT'])
        temp_summary <- as.data.frame(summary(coxph(coxph.surv~inmatrix_subset[,'CELLTYPE']))$coefficients[1,c(2,5)])

        # temp_summary_save <- as.data.frame(summary(coxph(coxph.surv~inmatrix_subset[,'CELLTYPE']))$coefficients)
        temp_summary <- t(temp_summary)
        colnames(temp_summary) <- c("HR", "pval")
        rownames(temp_summary) <- NULL
        temp_summary <- as.data.frame(temp_summary)
        temp_summary$signature <- celltype
        temp_summary$cohort <- cohort
        temp_summary$subtype <- subtype
        temp_summary$cutoff <- cutoff

        # log rank p-value
        temp_summary_logrank <- surv_pvalue(surv)
        temp_summary_logrank <- temp_summary_logrank[,"pval"]
        temp_summary$logrankpval <- temp_summary_logrank

        if(nrow(temp_df_combined)==0){
          temp_df_combined <- temp_summary
        } else {
          temp_df_combined <- rbind(temp_df_combined, temp_summary)
        }

        # colnames(temp_summary_save) <- paste0("METABRIC_",cohort)
        # write.csv(temp_summary_save,
        #           paste0("Clinical_Mean_",cutoff,"/Coxprop_stats_01_",cohort,"_KMcomparison.csv"))

        temp_plot <- ggsurvplot(surv,
                                data = inmatrix_subset,
                                title = paste(celltype,":", cohort, " ",subtype),
                                font.title = c(8, "bold", "darkblue"),
                                # palette = ncol1,
                                pval=T,
                                # pval = round(temp_summary[2,],digits=8),
                                legend="top",
                                # legend.labs = c("bottom 20%","top 20%"),
                                risk.table = T)

        png(file = paste0("Clinical_Mean_with_mid_year_cutoffs_",cutoff,"/",cohort,"_01_",celltype,"_",subtype,"_KMcomparison.png"),
            width = 8,
            height = 6,
            res=600,
            units="in")
        print(temp_plot)
        dev.off()

        pdf(file = paste0("Clinical_Mean_with_mid_year_cutoffs_",cutoff,"/",cohort,"_01_",celltype,"_",subtype,"_KMcomparison.pdf"),
            width = 8,
            height =6,
            onefile = F)
        print(temp_plot)
        dev.off()

        n <- paste0("inmatrix_subset_",cutoff, "_", cohort, "_", celltype, "_",subtype)
        assign(n, inmatrix_subset)

        write.csv(inmatrix_subset,
                  paste0("Clinical_Mean_with_mid_year_cutoffs_",cutoff,"/",cohort,"_01_",celltype,"_",subtype,"_clinical_data.csv"))

      }

    }
  }
}

write.csv(temp_df_combined,
          "Clinical_Mean_with_mid_year_cutoffs_25/pvalues_30_vs_30_yearcutoff_with_mid_group.csv")

# 10: DISTRIBUTION SCORES OF EACH SIGNATURE -----------------------------------

dir.create("distribution_scores")

for(cutoff in temp_cutoffs){
  print(cutoff)
  
  for(cohort in c("METABRIC", "TCGA")){
    print(cohort)
    # inmatrix_cohort <- get(paste0("expr_",cohort,"_mean_",cutoff))
    # 
    # if(length(inmatrix_cohort$PAM50SUBTYPE[is.na(inmatrix_cohort$PAM50SUBTYPE)])>0){
    #   inmatrix_cohort$PAM50SUBTYPE[is.na(inmatrix_cohort$PAM50SUBTYPE)] <- "no call"
    # }
    
      for(signature in "common_up_GeneList_HO_collapsed"){
        inmatrix_subset <- get(paste0("inmatrix_subset_",cutoff, "_", cohort, "_", signature, "_","All"))
        inmatrix_subset$patient <- rownames(inmatrix_subset)
        inmatrix_subset <- inmatrix_subset[inmatrix_subset$PAM50SUBTYPE %in% c("Basal"),,drop=F]
        
        inmatrix_subset$patient <- factor(inmatrix_subset$patient,
                                                            levels=inmatrix_subset$patient)
        
        temp_min <- min(inmatrix_subset$COMMON_UP_GENELIST_HO_COLLAPSED)
        temp_max <- max(inmatrix_subset$COMMON_UP_GENELIST_HO_COLLAPSED)
        
        inmatrix_subset$CELLTYPE <- factor(inmatrix_subset$CELLTYPE,
                                           levels=c("top30%", "mid40%", "bot30%"))
        
        temp_ggplot <- ggplot(inmatrix_subset,
                              aes(x = patient,
                                  y = COMMON_UP_GENELIST_HO_COLLAPSED,
                                  fill = CELLTYPE)) + 
          geom_bar(stat = "identity") +
          xlab(" ") + ylab("signature score") + ggtitle(paste0(cohort, "(",cutoff,")", " ",signature)) +
          theme(axis.text.x = element_blank()) +
          facet_wrap(. ~ PAM50SUBTYPE,
                     # scales="free_x",
                     nrow=1) +
          coord_cartesian(ylim = c(temp_min, temp_max)) +
          scale_fill_manual(values=c("top30%" = "blue", "mid40%" = "green", "bot30%" = "red"))
        
        png(file = paste0("distribution_scores/cutoff_",cutoff, "_",cohort, "_",signature,".png"),
            width = 5,
            height = 3.5,
            res=300,
            units="in")
        print(temp_ggplot)
        dev.off()
    }
  }
}

# 11: CLUSTER GENES -------------------------------------------------------

dir.create("Cluster_genes")

for(cutoff in temp_cutoffs){
  print(cutoff)
  
  for(cohort in c("TCGA", "METABRIC")){
    print(cohort)
    expr <- get(paste0("expr_",cohort))
    
    # temp_df_combined <- data.frame(row.names = c(1:top),
    #                                gene = c(1:top))
    
    for(signature in "common_up_GeneList_HO_collapsed"){
      print(signature)
      gene <- temp_genesig_df[temp_genesig_df$cluster == signature,,drop=F]
      gene <- gene[gene$gene %in% rownames(expr),,drop=F]
      gene <- gene[1:cutoff,,drop=F]
      
      n <- paste0("gene_",cohort)
      assign(n, gene)
      
      inmatrix_subset <- get(paste0("inmatrix_subset_",cutoff, "_", cohort, "_", signature, "_","Basal"))
      inmatrix_subset$patient <- rownames(inmatrix_subset)
      # inmatrix_subset <- inmatrix_subset[inmatrix_subset$PAM50SUBTYPE %in% c("Basal"),,drop=F]
      inmatrix_subset$patient <- factor(inmatrix_subset$patient,
                                        levels=inmatrix_subset$patient)
      
      

      
      # filter top 25 genes for expression matrix
      expr_filtered <- expr[rownames(expr) %in% gene$gene,,drop=F]
      # filter for patients
      expr_filtered <- expr_filtered[,colnames(expr_filtered) %in% as.vector(inmatrix_subset$patient),drop=F]

      # compute z-scores
      # if(cohort == "TCGA"){
        expr_filtered_zscore <- expr_filtered
        expr_filtered_zscore <- scale(t(expr_filtered_zscore))
        # for(row in c(1:nrow(expr_filtered))){
          # expr_filtered_zscore[row,] <- as.numeric(expr_filtered_zscore[row,])
        #   expr_filtered_zscore[row,] <- as.vector(scale(expr_filtered_zscore[row,],center = T))
        # }
      # } else {
      #   expr_filtered_zscore <- expr_filtered
      # }
        
      expr_filtered_zscore <- expr_filtered_zscore[as.vector(inmatrix_subset$patient),,drop=F]  
      print(all.equal(rownames(expr_filtered_zscore),
                      rownames(inmatrix_subset)))
      
      expr_filtered_zscore_plot <- t(expr_filtered_zscore)
      expr_filtered_zscore_plot <- expr_filtered_zscore_plot[as.vector(gene$gene),,drop=F]  
      
      temp_heatmap_annotation <- HeatmapAnnotation(group = as.vector(inmatrix_subset$CELLTYPE),
                                                   col = list(group = c("top30%" = "blue", "mid40%" = "green", "bot30%" = "red"))
                                                   )
      for(cluster in c(T,F)){
      temp_heatmap <- Heatmap(expr_filtered_zscore_plot,
                              show_column_names = FALSE,
                              cluster_columns = F,
                              cluster_rows = cluster,
                              clustering_distance_rows = "pearson",
                              width = unit(10, "cm"), 
                              height = unit(14, "cm"),
                              top_annotation = temp_heatmap_annotation)
        png(file = paste0("Cluster_genes/01_heatmap_",cohort,"_cluster_",cluster, ".png"),
            width = 7,
            height = 9,
            res=300,
            units="in")
        
        # draw(temp_heatmap, temp_heatmap_annotation)
        
        print(temp_heatmap)
        dev.off()
        
        n <- paste0("inmatrix_subset_",cohort)
        assign(n, inmatrix_subset)
      }
      
    }
  }
}

# print(intersect(gene_METABRIC$gene,
#                 gene_TCGA$gene))
# 
# print(setdiff(gene_METABRIC$gene,
#                 gene_TCGA$gene))

gene <- temp_genesig_df[temp_genesig_df$cluster == signature,,drop=F]
# gene <- gene[gene$gene %in% rownames(expr),,drop=F]
gene <- gene[1:25,,drop=F]

write.csv(gene,
          "top25_gene_signature.csv")

# 12: CORRELATE WITH CTL (BASALS) --------------------------------------------------

temp_genesig_df_IMMUNE <- data.frame(gene = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1"))
temp_genesig_df_IMMUNE$cluster <- "CTL levels"

{
  cutoff <- 5
  print(cutoff)
  expr_METABRIC_mean <- averaged_gene_signature(expr_METABRIC, temp_genesig_df_IMMUNE, top=cutoff)
  expr_METABRIC_mean = data.frame(clin_METABRIC[rownames(expr_METABRIC_mean),], expr_METABRIC_mean)
  
  colnames(expr_METABRIC_mean) <- gsub("Overall.Survival..Months.","OS",colnames(expr_METABRIC_mean))
  colnames(expr_METABRIC_mean) <- gsub("Overall.Survival.Status","EVENT",colnames(expr_METABRIC_mean))
  colnames(expr_METABRIC_mean) <- gsub("Age.at.Diagnosis","AGE",colnames(expr_METABRIC_mean))
  colnames(expr_METABRIC_mean) <- gsub("Pam50...Claudin.low.subtype","PAM50SUBTYPE",colnames(expr_METABRIC_mean))
  
  expr_TCGA_mean <-  averaged_gene_signature(expr_TCGA, temp_genesig_df_IMMUNE, top=cutoff)
  expr_TCGA_mean = data.frame(clin_TCGA[substr(rownames(expr_TCGA_mean),1,12),], expr_TCGA_mean)
  colnames(expr_TCGA_mean) <- gsub("subtype","PAM50SUBTYPE",colnames(expr_TCGA_mean))
  colnames(expr_TCGA_mean) <- gsub("age","AGE",colnames(expr_TCGA_mean))
  
  expr_GSE21653_mean <-  averaged_gene_signature(expr_GSE21653, temp_genesig_df_IMMUNE, top=cutoff)
  expr_GSE21653_mean = data.frame(clin_GSE21653[rownames(expr_GSE21653_mean),], expr_GSE21653_mean); colnames(expr_GSE21653_mean)[1] <- "OS"
  colnames(expr_GSE21653_mean) <- gsub("subtype","PAM50SUBTYPE",colnames(expr_GSE21653_mean))
  colnames(expr_GSE21653_mean) <- gsub("Age","AGE",colnames(expr_GSE21653_mean))
  
  expr_GSE58812_mean <-  averaged_gene_signature(expr_GSE58812, temp_genesig_df_IMMUNE, top=cutoff)
  expr_GSE58812_mean = data.frame(clin_GSE58812[rownames(expr_GSE58812_mean),], expr_GSE58812_mean); colnames(expr_GSE58812_mean)[2] <- "EVENT"
  colnames(expr_GSE58812_mean) <- gsub("subtype","PAM50SUBTYPE",colnames(expr_GSE58812_mean))
  colnames(expr_GSE58812_mean) <- gsub("Age","AGE",colnames(expr_GSE58812_mean))
  
  # expr_METABRIC_dis_mean <- averaged_gene_signature(expr_METABRIC_dis, gene, top=cutoff)
  # expr_METABRIC_dis_mean = data.frame(clin_METABRIC_dis[rownames(expr_METABRIC_dis_mean),], expr_METABRIC_dis_mean)
  # 
  # expr_METABRIC_val_mean <- averaged_gene_signature(expr_METABRIC_val, gene, top=cutoff)
  # expr_METABRIC_val_mean = data.frame(clin_METABRIC_val[rownames(expr_METABRIC_val_mean),], expr_METABRIC_val_mean)
  
  n <- paste0("expr_METABRIC_mean_CTL_",cutoff)
  assign(n, expr_METABRIC_mean)
  
  n <- paste0("expr_TCGA_mean_CTL_",cutoff)
  assign(n, expr_TCGA_mean)
  
  n <- paste0("expr_GSE21653_mean_CTL_",cutoff)
  assign(n, expr_GSE21653_mean)
  
  n <- paste0("expr_GSE58812_mean_CTL_",cutoff)
  assign(n, expr_GSE58812_mean)
  
  # n <- paste0("expr_METABRIC_dis_mean_",cutoff)
  # assign(n, expr_METABRIC_dis_mean)
  # n <- paste0("expr_METABRIC_val_mean_",cutoff)
  # assign(n, expr_METABRIC_val_mean)
  
  rm(expr_METABRIC_mean, expr_TCGA_mean,expr_GSE21653_mean,expr_GSE58812_mean)
}

# append to IE dataframe
CTL_cutoff <- 5
IE_cutoff <- 25
for(cohort in c("TCGA", "METABRIC")){
  
  expr_mean_CTL <- get(paste0("expr_",cohort, "_mean_CTL_",CTL_cutoff))
  expr_mean <- get(paste0("inmatrix_subset_",cohort))
  
  expr_mean_CTL <- expr_mean_CTL[rownames(expr_mean_CTL) %in% rownames(expr_mean),,drop=F]
  expr_mean_CTL <- expr_mean_CTL[rownames(expr_mean),,drop=F]
  
  if(all.equal(rownames(expr_mean),
               rownames(expr_mean_CTL))){
    print("rows match")
    expr_mean$CTLlevels <- expr_mean_CTL$CTLlevels
  }
  
  n <- paste0("expr_",cohort, "_mean_CTL_",CTL_cutoff)
  assign(n, expr_mean)
  
}


# 13: PLOT CORRELATIONS (BASALS) ---------------------------------------------------
 
# Just Basals
dir.create("TIL_correlations")
for(cohort in c("TCGA", "METABRIC")){
  
  expr_mean <- get(paste0("expr_",cohort, "_mean_CTL_",CTL_cutoff))
  expr_mean <- expr_mean[,colnames(expr_mean) %in% c("PAM50SUBTYPE", "COMMON_UP_GENELIST_HO_COLLAPSED", "CTLlevels", "CELLTYPE"),drop=F]  
  # if(length(expr_mean$PAM50SUBTYPE[is.na(expr_mean$PAM50SUBTYPE)])>0){
  #   expr_mean$PAM50SUBTYPE <- as.vector(expr_mean$PAM50SUBTYPE)
  #   expr_mean$PAM50SUBTYPE[is.na(expr_mean$PAM50SUBTYPE)] <- "no call"
  # }
  # expr_mean <- expr_mean[expr_mean$PAM50SUBTYPE == "Basal",,drop=F]
  # expr_mean <- expr_mean[,colnames(expr_mean) %in% c("common_up_GeneList_HO_collapsed", "CTLlevels"),drop=F]  
  # expr_mean$patient <- rownames(expr_mean)
  # rownames(expr_mean) <- NULL
  # 
  # temp_linearMod <- lm(CTLlevels ~ COMMON_UP_GENELIST_HO_COLLAPSED,
  #                      data = expr_mean)
  # temp <- summary(temp_linearMod)
  # temp_Rsquared <- paste0("adj.r.squared = ", round(temp$adj.r.squared, digits=3))
  

  # temp_ggplot <- 
  #   ggplot(temp_new_df_combined, 
  #          aes(x = FT, y = CRYO, color = cluster, linetype=condition)) +
  #   stat_smooth(method="loess", 
  #               se=FALSE,
  #               level = 0.95, 
  #               fullrange = TRUE) +
  #   theme(
  #     plot.title = element_text(color = '#666666',face = 'bold',size = 25
  #     ))
  # 
  # temp_pdf_function(paste0("06_cluster_level_correlations/loess_",sampleID,"correlations_by_cluster.pdf"))
  # print(temp_ggplot)
  # dev.off()
  temp_ggplot <-
    ggscatter(expr_mean, x = "COMMON_UP_GENELIST_HO_COLLAPSED", y = "CTLlevels", 
              add = "reg.line", conf.int = TRUE, size=1,
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "COMMON_UP_GENELIST_HO_COLLAPSED", ylab = "CTL levels")
  
  pdf(file = paste0("TIL_correlations/CTLlevels_vs_common_up_GeneList_HO_collapsed_",cohort, ".pdf"),
      width = 4,
      height = 4,
      onefile = F)
  print(temp_ggplot)
  dev.off()
  
  
  
  
}

# 12: CORRELATE WITH CTL (ALL BRCA) --------------------------------------------------

temp_genesig_df_IMMUNE <- data.frame(gene = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1"))
temp_genesig_df_IMMUNE$cluster <- "CTL levels"

{
  cutoff <- 5
  print(cutoff)
  expr_METABRIC_mean <- averaged_gene_signature(expr_METABRIC, temp_genesig_df_IMMUNE, top=cutoff)
  expr_METABRIC_mean = data.frame(clin_METABRIC[rownames(expr_METABRIC_mean),], expr_METABRIC_mean)
  
  colnames(expr_METABRIC_mean) <- gsub("Overall.Survival..Months.","OS",colnames(expr_METABRIC_mean))
  colnames(expr_METABRIC_mean) <- gsub("Overall.Survival.Status","EVENT",colnames(expr_METABRIC_mean))
  colnames(expr_METABRIC_mean) <- gsub("Age.at.Diagnosis","AGE",colnames(expr_METABRIC_mean))
  colnames(expr_METABRIC_mean) <- gsub("Pam50...Claudin.low.subtype","PAM50SUBTYPE",colnames(expr_METABRIC_mean))
  
  expr_TCGA_mean <-  averaged_gene_signature(expr_TCGA, temp_genesig_df_IMMUNE, top=cutoff)
  expr_TCGA_mean = data.frame(clin_TCGA[substr(rownames(expr_TCGA_mean),1,12),], expr_TCGA_mean)
  colnames(expr_TCGA_mean) <- gsub("subtype","PAM50SUBTYPE",colnames(expr_TCGA_mean))
  colnames(expr_TCGA_mean) <- gsub("age","AGE",colnames(expr_TCGA_mean))
  
  expr_GSE21653_mean <-  averaged_gene_signature(expr_GSE21653, temp_genesig_df_IMMUNE, top=cutoff)
  expr_GSE21653_mean = data.frame(clin_GSE21653[rownames(expr_GSE21653_mean),], expr_GSE21653_mean); colnames(expr_GSE21653_mean)[1] <- "OS"
  colnames(expr_GSE21653_mean) <- gsub("subtype","PAM50SUBTYPE",colnames(expr_GSE21653_mean))
  colnames(expr_GSE21653_mean) <- gsub("Age","AGE",colnames(expr_GSE21653_mean))
  
  expr_GSE58812_mean <-  averaged_gene_signature(expr_GSE58812, temp_genesig_df_IMMUNE, top=cutoff)
  expr_GSE58812_mean = data.frame(clin_GSE58812[rownames(expr_GSE58812_mean),], expr_GSE58812_mean); colnames(expr_GSE58812_mean)[2] <- "EVENT"
  colnames(expr_GSE58812_mean) <- gsub("subtype","PAM50SUBTYPE",colnames(expr_GSE58812_mean))
  colnames(expr_GSE58812_mean) <- gsub("Age","AGE",colnames(expr_GSE58812_mean))
  
  # expr_METABRIC_dis_mean <- averaged_gene_signature(expr_METABRIC_dis, gene, top=cutoff)
  # expr_METABRIC_dis_mean = data.frame(clin_METABRIC_dis[rownames(expr_METABRIC_dis_mean),], expr_METABRIC_dis_mean)
  # 
  # expr_METABRIC_val_mean <- averaged_gene_signature(expr_METABRIC_val, gene, top=cutoff)
  # expr_METABRIC_val_mean = data.frame(clin_METABRIC_val[rownames(expr_METABRIC_val_mean),], expr_METABRIC_val_mean)
  
  n <- paste0("expr_METABRIC_mean_CTL_",cutoff)
  assign(n, expr_METABRIC_mean)
  
  n <- paste0("expr_TCGA_mean_CTL_",cutoff)
  assign(n, expr_TCGA_mean)
  
  n <- paste0("expr_GSE21653_mean_CTL_",cutoff)
  assign(n, expr_GSE21653_mean)
  
  n <- paste0("expr_GSE58812_mean_CTL_",cutoff)
  assign(n, expr_GSE58812_mean)
  
  # n <- paste0("expr_METABRIC_dis_mean_",cutoff)
  # assign(n, expr_METABRIC_dis_mean)
  # n <- paste0("expr_METABRIC_val_mean_",cutoff)
  # assign(n, expr_METABRIC_val_mean)
  
  rm(expr_METABRIC_mean, expr_TCGA_mean,expr_GSE21653_mean,expr_GSE58812_mean)
}

# All BrCas
for(cutoff in temp_cutoffs){
  print(cutoff)
  
  for(cohort in c("TCGA", "METABRIC")){
    print(cohort)
    expr <- get(paste0("expr_",cohort))
    
    # temp_df_combined <- data.frame(row.names = c(1:top),
    #                                gene = c(1:top))
    
    for(signature in "common_up_GeneList_HO_collapsed"){
      print(signature)
      gene <- temp_genesig_df[temp_genesig_df$cluster == signature,,drop=F]
      gene <- gene[gene$gene %in% rownames(expr),,drop=F]
      gene <- gene[1:cutoff,,drop=F]
    
      inmatrix_subset <- get(paste0("inmatrix_subset_",cutoff, "_", cohort, "_", signature, "_","All"))
      inmatrix_subset$patient <- rownames(inmatrix_subset)
      # inmatrix_subset <- inmatrix_subset[inmatrix_subset$PAM50SUBTYPE %in% c("Basal"),,drop=F]
      inmatrix_subset$patient <- factor(inmatrix_subset$patient,
                                        levels=inmatrix_subset$patient)
      
      
      n <- paste0("inmatrix_subset_allgenes_",cohort)
      assign(n, inmatrix_subset)
    }
    
  }
}

# append to IE dataframe
CTL_cutoff <- 5
IE_cutoff <- 25
for(cohort in c("TCGA", "METABRIC")){
  
  expr_mean_CTL <- get(paste0("expr_",cohort, "_mean_CTL_",CTL_cutoff))
  expr_mean <- get(paste0("inmatrix_subset_allgenes_",cohort))
  
  expr_mean_CTL <- expr_mean_CTL[rownames(expr_mean_CTL) %in% rownames(expr_mean),,drop=F]
  expr_mean_CTL <- expr_mean_CTL[rownames(expr_mean),,drop=F]
  
  if(all.equal(rownames(expr_mean),
               rownames(expr_mean_CTL))){
    print("rows match")
    expr_mean$CTLlevels <- expr_mean_CTL$CTLlevels
  }
  
  n <- paste0("expr_",cohort, "_mean_allgenes_CTL_",CTL_cutoff)
  assign(n, expr_mean)
  
}

# 13: PLOT CORRELATIONS (ALL BRCA) ---------------------------------------------------

for(cohort in c("TCGA", "METABRIC")){

  expr_mean <- get(paste0("expr_",cohort, "_mean_allgenes_CTL_",CTL_cutoff))
  expr_mean <- expr_mean[,colnames(expr_mean) %in% c("PAM50SUBTYPE", "COMMON_UP_GENELIST_HO_COLLAPSED", "CTLlevels", "CELLTYPE"),drop=F]
  expr_mean <- expr_mean[expr_mean$PAM50SUBTYPE %in% c("Basal", "LumA", "LumB", "Her2"),,drop=F]

  # if(length(expr_mean$PAM50SUBTYPE[is.na(expr_mean$PAM50SUBTYPE)])>0){
  #   expr_mean$PAM50SUBTYPE <- as.vector(expr_mean$PAM50SUBTYPE)
  #   expr_mean$PAM50SUBTYPE[is.na(expr_mean$PAM50SUBTYPE)] <- "no call"
  # }
  # expr_mean <- expr_mean[expr_mean$PAM50SUBTYPE == "Basal",,drop=F]
  # expr_mean <- expr_mean[,colnames(expr_mean) %in% c("common_up_GeneList_HO_collapsed", "CTLlevels"),drop=F]
  # expr_mean$patient <- rownames(expr_mean)
  # rownames(expr_mean) <- NULL
  # 
  # temp_linearMod <- lm(CTLlevels ~ COMMON_UP_GENELIST_HO_COLLAPSED,
  #                      data = expr_mean)
  # temp <- summary(temp_linearMod)
  # temp_Rsquared <- paste0("R = ", round(temp$adj.r.squared, digits=3))

  # cor.test(expr_mean$CTLlevels, expr_mean$COMMON_UP_GENELIST_HO_COLLAPSED)

  # temp_ggplot <-
  #   ggplot(temp_new_df_combined,
  #          aes(x = FT, y = CRYO, color = cluster, linetype=condition)) +
  #   stat_smooth(method="loess",
  #               se=FALSE,
  #               level = 0.95,
  #               fullrange = TRUE) +
  #   theme(
  #     plot.title = element_text(color = '#666666',face = 'bold',size = 25
  #     ))
  #
  # temp_pdf_function(paste0("06_cluster_level_correlations/loess_",sampleID,"correlations_by_cluster.pdf"))
  # print(temp_ggplot)
  # dev.off()
  # if(cohort == "METABRIC"){
  #   minx <- 5.95
  #   maxy <- 8.9
  # }
  # if(cohort == "TCGA"){
  #   minx <- 2
  #   maxy <- 55
  # }
  # temp_ggplot <-
  #   ggplot(expr_mean,
  #          aes(x = COMMON_UP_GENELIST_HO_COLLAPSED,
  #              y = CTLlevels)) +
  #   geom_point(aes(colour = PAM50SUBTYPE),
  #     stat = "identity") +
  #   stat_smooth(method="lm",
  #               formula = y ~ x) +
  #   theme(
  #     plot.title = element_text(color = '#666666',face = 'bold',size = 25
  #     )) + annotate(geom="text",
  #                   x=minx,
  #                   y=maxy,
  #                   label=temp_Rsquared,
  #                   color="red")

  temp_ggplot <-
       ggscatter(expr_mean, x = "COMMON_UP_GENELIST_HO_COLLAPSED", y = "CTLlevels", size = 0.5,
            add = "reg.line", conf.int = TRUE,facet.by="PAM50SUBTYPE", nrow=2,
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "COMMON_UP_GENELIST_HO_COLLAPSED", ylab = "CTL levels")
  
  pdf(file = paste0("TIL_correlations/CTLlevels_vs_common_up_GeneList_HO_collapsed_allbrcas_",cohort, ".pdf"),
      width = 6,
      height =6,
      onefile = F)
  print(temp_ggplot)
  dev.off()




}



