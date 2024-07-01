########################################################
# Libraries
########################################################
library(data.table)
library(dplyr)



########################################################
# STEP 1: Load data
########################################################

# Load eductaion-related phenotype data. besides the identifier (eid) column, the file has two columns:
# age_end_edu (field #845): Negative values set to missing; Average over instances;
# qualifications (field # 6138): Negative values set to missing; All values that are not 1 (= College or University degree) are set to 0; Maximum over instances retained 
pheno <- as.data.frame(fread("education_phenotypes.txt"))



########################################################
# STEP 2: Define composite EA phenotype
########################################################

# Modify EA
phenotypes$age_end_edu[which(phenotypes$age_end_edu < 0)] <- NA 
phenotypes$age_end_edu[which(phenotypes$age_end_edu < 14)] <- 14
phenotypes$age_end_edu[which(phenotypes$age_end_edu > 19)] <- 19 
phenotypes$age_end_edu[which(phenotypes$qualifications == 1)] <- 19
phenotypes <- select(phenotypes, -qualifications) 

# Save results
fwrite(phenotypes, "EA_WB_raw_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
