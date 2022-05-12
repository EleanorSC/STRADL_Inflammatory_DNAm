## ---------------------------
##
## Script Purpose: Analysing proteomic STRADL data with neuroimaging outcomes
##   
##   
##
## Notes: you will need both: 
##        (a) a .csv file containing: Master data (protein, neuroimaging, cognitive, demographics)
##        (b) a .csv file containing: 

##
## In the scirpt below, these files are listed as:
##   "protein_data_full.csv" = (protein, neuroimaging, cognitive, demographics data)
##
## ---------------------------

setwd()

protein_data_full <- protein_data_full[ -c(1:6,8:32) ]

# try to subset this instead
protein_neuroimaging <- protein_data_full[complete.cases(protein_data_full$Global_GM_Volume),]

# Make a list of just the proteins 

protein_list <- list(annotation_data$`Entrez Gene Name`) 


protein_list <- list(
  "TPP1",
  "NEU1",
  "FABP1",
  "ASGR1",
  "TRAPPC3",
  "CANT1",
  "TREM1",
  "GDF15",
  "HEXB",
  "SVEP1",
  "PLXDC1",
  "ITIH4",
  "RBL2",
  "AFM",
  "ACY1",
  "CFHR1",
  "F9",
  "SERPIND1",
  "SIGLEC12",
  "LMAN2",
  "ADAMTSL2",
  "IGLON5",
  "NCAN",
  "GLIPR2",
  "RBL2",
  
  ### Ones Danni did not select:
  "POR",
  "PDK1",
  "OLFM2",
  "ECI2",
  "DDX58",
  "MX1",
  "MMAB",
  "STXBP6",
  "CST5",
  "IAPP",
  "MATN3"
)

#proteins_only <- protein_data_full[, c(2:157)]
#protein_list <- colnames(proteins_only)

protein_regressions <- list()

# First we will examine just GM volume

GM_null <- summary(lm(scale(Global_GM_Volume) ~ scale(st_age) + sex, 
                      data = protein_neuroimaging))

##

i <- 1

for (protein in protein_list) { 
  #print(protein)
  linear_mod_output <- summary(lm(scale(Global_GM_Volume) ~ scale(st_age) + sex + scale(protein_neuroimaging[protein]), 
                                  data = protein_neuroimaging))
  rownames(linear_mod_output$coefficients)[4] <- protein
  #rownames(linear_mod_output$coefficients)[2] <- age 
  protein_regressions[[i]] <- linear_mod_output
  i <- (i + 1)
}

#This function outputs my various regression summaries into a table
summary_table_new <- function(protein_regressions, lmNULL, protein_list, inflam_protein) {
  
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
                    protein = unlist(inflam_protein)
                    
  ))
}


#Here tell the function which regressions to use
GM <- summary_table_new(protein_regressions, GM_null, protein_list, protein_list)
GM <- mutate(GM, neuroimaging = "GM")
GM <- mutate(GM, Trait = "GM")
GM <- mutate(GM, proteomics = "protein")


## Plot data --------------------------------

plot <- GM


# If only selecting the ones Danni did:
#select_proteins <- list("TPP1", #17691-1
#                        "NEU1", #15426-5
#                        "FABP1", #18888-37
#                        "ASGR1", #5452-71
#                        "TRAPPC3", #14337-1
#                        "CANT1", #6480-1
#                        "TREM1", #
#                        "GDF15",
#                        "HEXB",
#                        "SVEP1",
#                        "PLXDC1",
#                        "ITIH4",
#                        "RBL2",
#                        "AFM",
#                        "ACY1",
#                        "CFHR1",
#                        "F9",
#                        "SERPIND1",
#                        "SIGLEC12",
#                        "LMAN2",
#                        "ADAMTSL2",
#                        "IGLON5",
#                        "NCAN",
#                        "GLIPR2",
#                        "RBL2")

# Drop some associations if they are too small

#plot <- plot  %>% 
#  mutate(what_to_drop =
#         case_when(beta <= -0.02 ~ 'Pass',
#                   beta >= 0.02 ~ 'Pass',
#                   TRUE ~ 'Fail'))
#
#plot <- plot %>% filter(what_to_drop == "Pass")

# positive or negative association ?

#plot <- plot  %>% 
#  mutate(positive_negative =
#           case_when(beta <= 0 ~ 'negative',
#                     TRUE ~ 'positive'))

# p-value < 0.05 

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
                        y=beta
                       # group = highlight_flag
                        )
) +
  
  geom_bar(aes(fill=reorder(protein, -beta)
              # alpha = highlight_flag
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


mycolors <- wesanderson::wes_palette("Zissou1", 30, type = "continuous")

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
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 8)) +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.title.y = element_text(face = "bold")) 

p


#skimr::skim(plot$highlight_flag)

