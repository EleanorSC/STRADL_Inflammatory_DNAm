## ---------------------------
##
## Script Purpose: Analysing proteomic STRADL data trained in KORA with neuroimaging outcomes
##   
##   
##
## Notes: you will need both: 
##        (a) a .csv file containing: Master data (protein, neuroimaging, cognitive, demographics)
##        (b) a .csv file containing: 

##
## In the scirpt below, these files are listed as:
##   "protein_data_full.csv" = (protein, neuroimaging, cognitive, demographics data)
##   "STRADL_scores_trained_KORA.csv" = STRADL Episcores trained in KORA
##   ""Annotation_data.csv" = Linked UniProt IDs with SeqIds
##
## This script will
##        (a) Load the KORA trained EpiScores
##
##
## ---------------------------
protein_data_full <- read.csv("prot_file_150621.csv", check.names=FALSE)
names(protein_data_full)[names(protein_data_full) == 'SampleId'] <- 'stradl_ID'
# Remove columns we don't need
protein_data_full <- protein_data_full[ -c(1:6,8:32) ]

protein_data_only <- protein_data_full[c(1:4236)]
STRADL_main_data <- protein_data_full[c(1, 4237:4286)]




STRADL_KORA <- read.csv("STRADL_scores_trained_KORA.csv", check.names=FALSE)
#annotation_data <- read.csv("Annotation_data.csv", check.names=FALSE)

newnames <- c()
for(colname in names(STRADL_KORA)){
  if(colname %in% annotation_data$SeqId){
    egname <- annotation_data[annotation_data$SeqId == colname,]$`Entrez Gene Name`
    newnames <- append(newnames, egname)
  } else {
    newnames <- append(newnames, colname)
  }
}

print(newnames)

names(STRADL_KORA) <- newnames

# Now link these to STRADL ID
names(STRADL_KORA)[names(STRADL_KORA) == 'ID'] <- 'meth_ID'

# Load ID linkages
#ID_link_GS_STRADL_DNAm <- readRDS("Methylation_Genotype_Imaging_linkageID.rds")
ID_link_attempt <- read.csv("STRADL_DNAm_target_REM_17April2020.csv")
Meth_link <- ID_link_attempt %>% select("stradl_ID", "meth_ID")

STRADL_KORA <- merge(Meth_link, 
                     STRADL_KORA,
                     by = "meth_ID")

###########

DNAm_STRADL_KORA <- merge(STRADL_KORA, 
                          STRADL_main_data,
                                  by = "stradl_ID")

# Find out what proteins we have
KORA_DNAm_list <- ls(STRADL_KORA[c(3:86)])

KORA_DNAm_list <- list("ACY1",
                     "ADAMTS13",
                     "ADIPOQ",
                     "AFM",
                     "B2M",
                     "BCAM",
                     "BMP1",
                     "C4A|C4B",
                     "C5",
                     "C9",
                     "CCL17",
                     "CCL18",
                     "CCL21",
                     "CCL22",
                     "CCL25",
                     "CD163",
                     "CD209",
                     "CD48",
                     "CD5L",
                     "CHIT1",
                     "CLEC11A",
                     #"CLEC11A.1",
                     "CNTN4",
                     "CRP",
                     "CXCL10",
                     "CXCL11",
                     "EDA",
                     "ENPP7",
                     "ESM1",
                     "F7",
                     "FAP",
                     "FCER2",
                     "FCGR3B",
                     "GHR",
                     "GNLY",
                     "GP1BA",
                     "GZMA",
                     "HGFAC",
                     "ICAM5",
                     "IDUA",
                     "IGFBP1",
                     "IGFBP4",
                     "IL19",
                     "INSR",
                     "LGALS3BP", 
                     "LGALS4",
                     "LTA|LTB",
                     "LTF",
                     "LY9",
                     "LYZ",
                     "MIA",
                     "MMP1",
                     "MMP12",
                     "MMP2",
                     "MMP9",
                     "MPL",
                     "MPO",
                     "MRC2",
                     "MST1",
                     "NCAM1",
                     "NOTCH1",
                     "NTRK3",
                     "OMD",
                     "PAPPA",
                     "PIGR",
                     "PRSS2",
                     "RARRES2",
                     "RETN",
                     "S100A9",
                     "SELE",
                     "SELL",
                     "SEMA3E",
                     "SERPINA3",
                     "SERPIND1",
                     "SHBG",
                     "SLITRK5",
                     "SPOCK2",
                     "STC1",
                     "THBS2",
                     "TNFRSF17",
                     "TNFRSF1B",
                     "TPSB2",
                     "VCAM1",
                     "WFIKKN2")

DNAm_regressions <- list()

# First we will examine just GM volume

GM_null <- summary(lm(scale(Global_GM_Volume) ~ scale(st_age) + sex, 
                      data = DNAm_STRADL_KORA))

##

i <- 1
for (DNAm in KORA_DNAm_list) { 
  linear_mod_output <- summary(lm(scale(Global_GM_Volume) ~ scale(st_age) + sex + scale(DNAm_STRADL_KORA[DNAm]), 
                                  data = DNAm_STRADL_KORA))
  rownames(linear_mod_output$coefficients)[4] <- DNAm
  DNAm_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

#This function outputs my various regression summaries into a table
summary_table_new <- function(DNAm_regressions, lmNULL, KORA_DNAm_list, inflam_protein) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  r2d <- NULL
  
  for (regression in DNAm_regressions) {
    beta <- c(beta, regression$coefficients[4,1])
    se <- c(se, regression$coefficients[4,2])
    pvals <- c(pvals, regression$coefficients[4,4])
    r2 <- c(r2, regression$r.squared)
    r2d <- c(r2d, regression$r.squared - lmNULL$r.squared) 
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    r2d = r2d,
                    inflammation = unlist(KORA_DNAm_list),
                    protein = unlist(KORA_DNAm_list)
                    
  ))
}


#Here tell the function which regressions to use
GM <- summary_table_new(DNAm_regressions, GM_null, KORA_DNAm_list, KORA_DNAm_list)
GM <- mutate(GM, neuroimaging = "GM")
GM <- mutate(GM, Trait = "GM")
GM <- mutate(GM, proteomics = "DNAm")


#########
#Now for the proteins
# Find out what proteins we have
protein_list <- list("ACY1",
                   "ADAMTS13",
                   "ADIPOQ",
                   "AFM",
                  "B2M",
                   "BCAM",
                   "BMP1",
                  "C4A|C4B",
                  "C5",
                  "C9",
                  "CCL17",
                  "CCL18",
                  "CCL21",
                  "CCL22",
                  "CCL25",
                  "CD163",
                  "CD209",
                  "CD48",
                  "CD5L",
                  "CHIT1",
                  "CLEC11A",
                  #"CLEC11A.1",
                  "CNTN4",
                  "CRP",
                  "CXCL10",
                  "CXCL11",
                  "EDA",
                  "ENPP7",
                  "ESM1",
                  "F7",
                  "FAP",
                  "FCER2",
                  "FCGR3B",
                  "GHR",
                  "GNLY",
                  "GP1BA",
                  "GZMA",
                  "HGFAC",
                  "ICAM5",
                  "IDUA",
                  "IGFBP1",
                  "IGFBP4",
                  "IL19",
                  "INSR",
                  "LGALS3BP", 
                  "LGALS4",
                  "LTA|LTB",
                  "LTF",
                  "LY9",
                  "LYZ",
                  "MIA",
                  "MMP1",
                  "MMP12",
                  "MMP2",
                  "MMP9",
                  "MPL",
                  "MPO",
                  "MRC2",
                  "MST1",
                  "NCAM1",
                  "NOTCH1",
                  "NTRK3",
                  "OMD",
                  "PAPPA",
                  "PIGR",
                  "PRSS2",
                  "RARRES2",
                  "RETN",
                  "S100A9",
                  "SELE",
                  "SELL",
                  "SEMA3E",
                  "SERPINA3",
                  "SERPIND1",
                  "SHBG",
                  "SLITRK5",
                  "SPOCK2",
                  "STC1",
                  "THBS2",
                  "TNFRSF17",
                  "TNFRSF1B",
                  "TPSB2",
                  "VCAM1",
                  "WFIKKN2")

#protein_list <- KORA_DNAm_list
protein_regressions <- list()

# First we will examine just GM volume

GM_null <- summary(lm(scale(Global_GM_Volume) ~ scale(st_age) + sex, 
                      data = protein_data_full))

##
i <- 1
for (protein in protein_list) { 
  print(protein)
  linear_mod_output <- summary(lm(scale(Global_GM_Volume) ~ scale(st_age) + sex + scale(protein_data_full[protein]), 
                                  data = protein_data_full))
  rownames(linear_mod_output$coefficients)[4] <- protein
  protein_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

#This function outputs my various regression summaries into a table
summary_table_new2 <- function(protein_regressions, lmNULL, protein_list, inflam_protein) {
  
  beta <- NULL
  se <- NULL
  pvals <- NULL
  r2 <- NULL
  r2d <- NULL
  
  for (regression in protein_regressions) {
    beta <- c(beta, regression$coefficients[4,1])
    se <- c(se, regression$coefficients[4,2])
    pvals <- c(pvals, regression$coefficients[4,4])
    r2 <- c(r2, regression$r.squared)
    r2d <- c(r2d, regression$r.squared - lmNULL$r.squared) 
  }
  
  return(data.frame(beta = beta,
                    se = se,
                    pvals = pvals,
                    r2 = r2,
                    r2d = r2d,
                    inflammation = unlist(protein_list),
                    protein = unlist(protein_list)
                    
  ))
}


#Here tell the function which regressions to use
GM2 <- summary_table_new2(protein_regressions, GM_null, protein_list, protein_list)
GM2 <- mutate(GM2 , neuroimaging = "GM")
GM2 <- mutate(GM2 , Trait = "GM")
GM2 <- mutate(GM2, proteomics = "protein")

#### PLOT

plot <- rbind(GM, GM2)

plot <- plot  %>% 
  mutate(significant =
           case_when(pvals <= 0.05 ~ 'significant',
                     TRUE ~ 'null'))

plot <- plot %>% filter(significant == "significant")

#highlight the 25
#plot <- plot %>% mutate(highlight_flag = ifelse(protein %in% select_proteins, T, F)) 

##
beta_plot <- ggplot(data=plot,
                    aes(x=reorder(protein, -beta),
                        y=beta,
                         group = proteomics
                    )
) +
  
  geom_bar(aes(fill=reorder(protein, -beta),
                alpha = proteomics
  ),
  stat="identity", 
  colour="black",
  alpha = 0.7,
  position=position_dodge()) +
  
  geom_errorbar(aes(ymin=beta-se, ymax=beta+se), 
                position = position_dodge(0.9),
                width=0.3,
                colour="darkgrey", alpha=0.9, size=0.8) +
  coord_flip() +
  xlab("GM") + 
  ylab("Standardised effect size") 


mycolors <- wesanderson::wes_palette("Zissou1", 84, type = "continuous")

p <- beta_plot + 
  #facet_wrap(vars(proteomics),  scales = "free_y", ncol =3) +
  #facet_wrap(vars(positive_negative),  scales = "free_y", ncol =2) +
  scale_fill_manual(values = mycolors) +
  theme_bw() + 
  theme(legend.position= c(0.81, 0.93)) +
  guides(fill = FALSE) +
  labs(alpha = "inflammation")+
  theme(legend.position = "bottom") +
  #ylim(-0.16, 0.19) +
  theme(axis.text.x = element_text(size = 9), 
        axis.text.y = element_text(size = 9)) +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.title.y = element_text(face = "bold")) 

p






##### Do for the protein data
#newnames <- c()
#for(colname in names(protein_data_full)){
#  if(colname %in% annotation_data$SeqId){
#    egname <- annotation_data[annotation_data$SeqId == colname,]$`Entrez Gene Name`
#    newnames <- append(newnames, egname)
#  } else {
#    newnames <- append(newnames, colname)
#  }
#}
#
#print(newnames)
#
#names(protein_data_full) <- newnames
#