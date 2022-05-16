## ---------------------------
##
## Script Purpose: Analysing proteomic STRADL data trained in KORA with neuroimaging outcomes
##                 Plot for only specific protein / DNAm assocs 
##   
##   
##
# Load DNAm EpiScores data trained in STRADL
EpiScores <- read.csv("STRADL_scores_trained_LBC.csv")

names(EpiScores)[names(EpiScores) == 'ID'] <- 'meth_ID'

DNAm_STRADL <- merge(ID_link_attempt, 
                     EpiScores,
                     by = "meth_ID")

# Remove columns we don't need
DNAm_STRADL_only <- DNAm_STRADL[c(2,7:31)]

# Merge with other variables
DNAm_STRADL_LBC <- merge(DNAm_STRADL_only, 
                         STRADL_main_data,
                         by = "stradl_ID")

##
DNAm_LBC_list <- ls(DNAm_STRADL_only[-c(1)])

DNAm_regressions <- list()

# First we will examine just GM volume

GM_null <- summary(lm(scale(Global_GM_Volume) ~ scale(st_age) + sex, 
                      data = DNAm_STRADL_LBC))

##

i <- 1
for (DNAm in DNAm_LBC_list) { 
  linear_mod_output <- summary(lm(scale(Global_GM_Volume) ~ scale(st_age) + sex + scale(DNAm_STRADL_LBC[DNAm]), 
                                  data = DNAm_STRADL_LBC))
  rownames(linear_mod_output$coefficients)[4] <- DNAm
  DNAm_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

#This function outputs my various regression summaries into a table
summary_table_new <- function(DNAm_regressions, lmNULL, DNAm_LBC_list, inflam_protein) {
  
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
                    inflammation = unlist(DNAm_LBC_list),
                    protein = unlist(DNAm_LBC_list)
                    
  ))
}


#Here tell the function which regressions to use
GM_LBC_DNAm <- summary_table_new(DNAm_regressions, GM_null, DNAm_LBC_list, DNAm_LBC_list)
GM_LBC_DNAm <- mutate(GM_LBC_DNAm, neuroimaging = "GM")
GM_LBC_DNAm <- mutate(GM_LBC_DNAm, Trait = "GM")
GM_LBC_DNAm <- mutate(GM_LBC_DNAm, Trained = "LBC")
GM_LBC_DNAm <- mutate(GM_LBC_DNAm, proteomics = "DNAm")

#########
#Now for the proteins
# Find out what proteins we have
#protein_list <- DNAm_LBC_list

protein_list <- list(
  "CCL11",
  #"CD6",
  "CRTAM",
  "CXCL10",
  "CXCL11",
  "CXCL9",
  # "EN.RAGE",
  "EZR",
  #  "FcRL2",
  #  "FGF.21",
  #  "G.CSF",
  #  "GDF.8",
  "GZMA",
  "HGF",
  #  "MMP.1",
  #  "N.CDase",
  #  "NEP",
  "NMNAT1",
  "NTRK3",
  "OSM",
  #  "SIGLEC1",
  #  "SKR3",
  "SMPD1",
  #  "TGF.alpha",
  "VEGFA"
)


#protein_list <- KORA_DNAm_list
protein_regressions <- list()

# First we will examine just GM volume

GM_null <- summary(lm(scale(Global_GM_Volume) ~ scale(st_age) + sex, 
                      data = protein_data_full))

##
i <- 1
for (protein in protein_list) { 
  print(protein)
  # Build an if loop here where if protein not present it skips but sends a warning!
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
GM_LBC_protein <- summary_table_new2(protein_regressions, GM_null, protein_list, protein_list)
GM_LBC_protein <- mutate(GM_LBC_protein , neuroimaging = "GM")
GM_LBC_protein <- mutate(GM_LBC_protein , Trait = "GM")
GM_LBC_protein<- mutate(GM_LBC_protein, Trained = "LBC")
GM_LBC_protein <- mutate(GM_LBC_protein, proteomics = "protein")



######KORA trained
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
GM_KORA_DNAm <- summary_table_new2(DNAm_regressions, GM_null, KORA_DNAm_list, KORA_DNAm_list)
GM_KORA_DNAm <- mutate(GM_KORA_DNAm , neuroimaging = "GM")
GM_KORA_DNAm <- mutate(GM_KORA_DNAm , Trait = "GM")
GM_KORA_DNAm<- mutate(GM_KORA_DNAm, Trained = "KORA")
GM_KORA_DNAm <- mutate(GM_KORA_DNAm, proteomics = "DNAm")


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
GM_KORA_protein <- summary_table_new2(protein_regressions, GM_null, protein_list, protein_list)
GM_KORA_protein <- mutate(GM_KORA_protein , neuroimaging = "GM")
GM_KORA_protein <- mutate(GM_KORA_protein , Trait = "GM")
GM_KORA_protein<- mutate(GM_KORA_protein, Trained = "KORA")
GM_KORA_protein <- mutate(GM_KORA_protein, proteomics = "protein")




plot <- rbind(GM_LBC_protein, 
              GM_LBC_DNAm,
              GM_KORA_protein,
              GM_KORA_DNAm)



#plot <- plot %>% subset(-c("")
#                        
#drop <- c("GP1BA","GNLY", "SELL","SPOCK2",
#          "NOTCH1")
#
#

plot <- plot %>% filter(
  Trained == "KORA" &  protein == "STC1" |
  Trained == "KORA" &  protein == "RARRES2" |
  Trained == "KORA" &  protein == "CCL18" |
  Trained == "KORA" &  protein == "CRP" |
  Trained == "KORA" &  protein == "IGFBP4" |
  Trained == "KORA" &  protein == "HGFAC" |
  Trained == "KORA" &  protein == "S100A9" |
  Trained == "KORA" &  #  protein == "GP1BA" |
  Trained == "KORA" &  protein == "MMP1" |
  Trained == "KORA" &  # protein == "GNLY" |
  Trained == "KORA" &  protein == "CNTN4" |
  Trained == "KORA" &  # protein == "SELL" |
  Trained == "KORA" &  protein == "MRC2" |
  Trained == "KORA" &  protein == "FAP" |
  Trained == "KORA" &#  protein == "NTRK3" |
  Trained == "KORA" &  protein == "SEMA3E" |
  Trained == "KORA" &  protein == "ADIPOQ" |
  Trained == "KORA" &  protein == "CNTN4" |
  Trained == "KORA" &  protein == "MMP2" |
  Trained == "KORA" &  protein == "S100A9" |
  Trained == "KORA" &  protein == "THBS2" | 
    
    # LBC to keep
  Trained == "LBC" &   protein == "CCL11" |
  Trained == "LBC" &   protein == "CRTAM"|
# Trained == "LBC" & #  protein == "CXCL10"|
# Trained == "LBC" & #  protein == "CXCL11"|
  Trained == "LBC" &   protein == "CXCL9"|
  Trained == "LBC" &   protein == "EZR"|
  Trained == "LBC" &   protein == "GZMA"|
  Trained == "LBC" &   protein == "HGF"|
  Trained == "LBC" &   protein == "NMNAT1"|
  Trained == "LBC" &   protein == "NTRK3"|
  Trained == "LBC" &   protein == "OSM"|
  Trained == "LBC" &   protein == "SMPD1"|
  Trained == "LBC" &   protein == "VEGFA"  )
  
)





plot$concat <-  paste(plot$Trained, plot$proteomics)


plot <- plot %>% filter(Trained == "KORA")

##
beta_plot <- ggplot(data=plot,
                    aes(x=reorder(protein, -beta),
                        y=beta,
                        group = concat,
                        alpha = proteomics
                    )
) +
  
 geom_bar(aes(fill = Trained),
   #fill=reorder(protein, -beta),
              #alpha = Trained),
          
          stat="identity", 
          colour="black",
          #alpha = 0.7,
          position=position_dodge()) +
  
  geom_errorbar(aes(ymin=beta-se, ymax=beta+se), 
                position = position_dodge(0.9),
                width=0.3,
                colour="darkgrey", alpha=0.9, size=0.8) +
  coord_flip() +
  xlab("GM") + 
  ylab("Standardised effect size") 


#mycolors <- wesanderson::wes_palette("Zissou1", 32, type = "continuous")


p <- beta_plot + 
 # facet_wrap(vars(proteomics),  scales = "free_y", ncol =2) +
  #facet_wrap(vars(positive_negative),  scales = "free_y", ncol =2) +
#  scale_fill_manual(values = mycolors) +
  theme_bw() + 
  theme(legend.position= c(0.81, 0.93)) +
  guides(fill = FALSE) +
  #labs(alpha = "inflammation")+
  theme(legend.position = "bottom") +
  #ylim(-0.16, 0.19) +
  theme(axis.text.x = element_text(size = 10.5), 
        axis.text.y = element_text(size = 10.5)) +
  theme(axis.title.x = element_text(face = "bold",size = 11)) +
  theme(axis.title.y = element_text(face = "bold",size = 11)) 

p



####

# GENE2FUNC


# Proteins we're interested in

#STC1 
#RARRES2 
#CCL18 
#CRP 
#IGFBP4 
#HGFAC 
#S100A9 
#GP1BA 
#MMP1 
#GNLY 
#CNTN4 
#SELL 
#MRC2 
#FAP 
#NTRK3 
#SEMA3E 
#ADIPOQ 
#CNTN4 
#MMP2 
#S100A9 
#THBS2 




