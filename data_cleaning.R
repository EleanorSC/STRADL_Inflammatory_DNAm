# Data cleaning 
install.packages("tidyverse")
install.packages("magrittr")

library(tidyverse)
library(magrittr)

ID_link_GS_STRADL_DNAm <- readRDS("Methylation_Genotype_Imaging_linkageID.rds")

BMI_covariates <- read.csv("Depression_Inflammation_DNAm_final.csv")

BMI_covariates %<>% select(ID, sex, age,
                          bmi, whr, body_fat,
                          diabetes_Y, osteo_arthritis_Y, rheum_arthritis_Y, depression_Y,
                          asthma_Y, COPD_Y, inflam_disease
                          )

