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
##   "Master_data.csv" = (protein, neuroimaging, cognitive, demographics data)
##
## ---------------------------

setwd()

# Make a list of just the proteins 

proteins_only <- Master_data[, c(2:157)]

protein_list <- colnames(proteins_only)

protein_regressions <- list()

# First we will examine just GM volume

GM_null <- summary(lm(scale(Global_GM_Volume) ~ scale(st_age) + sex, 
                      data = Master_data))

##

i <- 1

for (protein in protein_list) { 
  #print(protein)
  linear_mod_output <- summary(lm(scale(Global_GM_Volume) ~ scale(st_age) + sex + scale(Master_data[protein]), 
                                  data = Master_data))
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

##
beta_plot <- ggplot(data=plot,
                    aes(x=reorder(protein, -beta),
                        y=beta,
                        group = proteomics)
) +
  
  geom_bar(aes(fill=reorder(protein, -beta), 
               alpha = proteomics),
           stat="identity", 
           colour="black",
           position=position_dodge()) +
  
  geom_errorbar(aes(ymin=beta-se, ymax=beta+se), 
                position = position_dodge(0.9),
                width=0.3,
                colour="darkgrey", alpha=0.9, size=0.8) +
  coord_flip() +
  xlab("GM") + 
  ylab("Standardised effect size") 


mycolors <- wesanderson::wes_palette("Zissou1", 156, type = "continuous")

p <- beta_plot + 
  facet_wrap(vars(proteomics),  scales = "free_y", ncol =3) +
  scale_fill_manual(values = mycolors) +
  theme_bw() + 
  theme(legend.position= c(0.81, 0.93)) +
  guides(fill = FALSE) +
  labs(alpha = "inflammation")+
  theme(legend.position = "bottom") +
  ylim(-0.16, 0.19) +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.title.y = element_text(face = "bold")) 

p