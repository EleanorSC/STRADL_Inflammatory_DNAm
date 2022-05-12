## ---------------------------
##
## Script Purpose: aligning protein SeqID with Proteijn name
##   
##   
##
## Notes: you will need both: 
##        (a) a .csv file containing: Master data (protein, neuroimaging, cognitive, demographics)
##        (b) a .csv file containing the uniprot IDs

##
## 
##
## ---------------------------
protein_data_full <- read.csv("prot_file_150621.csv", check.names=FALSE)

annotation_data <- read.csv("Annotation_data.csv", check.names=FALSE)

newnames <- c()
for(colname in names(protein_data_full)){
  if(colname %in% annotation_data$SeqId){
    egname <- annotation_data[annotation_data$SeqId == colname,]$`Entrez Gene Name`
    newnames <- append(newnames, egname)
  } else {
    newnames <- append(newnames, colname)
  }
}

print(newnames)

names(protein_data_full) <- newnames