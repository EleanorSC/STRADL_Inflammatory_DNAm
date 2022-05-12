# Data cleaning 
install.packages("tidyverse")
install.packages("magrittr")

library(tidyverse)
library(magrittr)

# Load ID linkages
ID_link_GS_STRADL_DNAm <- readRDS("Methylation_Genotype_Imaging_linkageID.rds")

ID_link_attempt <- read.csv("STRADL_DNAm_target_REM_17April2020.csv")

BMI_covariates <- read.csv("Depression_Inflammation_DNAm_final.csv")

BMI_covariates %<>% select(ID, sex, age,
                          bmi, whr, body_fat,
                          diabetes_Y, osteo_arthritis_Y, rheum_arthritis_Y, depression_Y,
                          asthma_Y, COPD_Y, inflam_disease
                          )

# Load the DNAm CRP scores
CRP_DNAm <- readRDS("GS_DNAm_CRP_score.rds")

names(CRP_DNAm)[names(CRP_DNAm) == 'ID'] <- 'meth_ID'

DNAm_CRP_STRADL <- merge(ID_link_GS_STRADL_DNAm,CRP_DNAm,by = "meth_ID")

# Note we are missing 772 / 801 CRP values ? different assay?
skimr::skim(DNAm_CRP_STRADL$crp)

# Load CRP blood data
blood_data <- read.csv("Bloods.csv")

# Load DNAm EpiScores data

EpiScores <- read.csv("STRADL_scores_trained_LBC.csv")

names(EpiScores)[names(EpiScores) == 'ID'] <- 'meth_ID'

DNAm_STRADL <- merge(ID_link_attempt, 
                     EpiScores,
                     by = "meth_ID")

# Link with BMI attempt
names(BMI_covariates)[names(BMI_covariates) == 'ID'] <- 'GS_id'

#n.b drop sex from BMI covariates 

DNAm_STRADL_BMI <- merge(DNAm_STRADL, 
                         BMI_covariates,
                     by = "GS_id")


