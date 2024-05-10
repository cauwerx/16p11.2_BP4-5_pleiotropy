########################################################
# Libraries
########################################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)



########################################################
# STEP 1: Load data & define parameters
########################################################

# Parameter definition: Number of controls per CNV carrier (25)
n <- 25

# Load unrelated white British individuals
# Load sample file (vector with identifiers, with header "IID")
samples <- as.data.frame(fread("samples_white_british_All.txt", col.names = "IID"))

# Load CNV file and subset samples (output from CNV calling pipeline, as described in Auwerx et al., 2022 (AJHG))
cnvs <- as.data.frame(fread("/Users/cauwerx/Documents/Work/PhD/Projects/CNV_GWA/UKBB_500K_CNV/GENERAL_DATA/ukb_cnv_global.gz"))
cnvs$Sample_Name <- sub("_.*", "", cnvs$Sample_Name)
cnvs <- cnvs[cnvs$Sample_Name %in% samples$eid, ]

# Load matching factor data (identifier + 1 column/factor, including age, sex, BMI, TDI, income class, EA; ethnicity is matched by only working on the "white British" data subset)
cov <- as.data.frame(fread("matching_factors.txt"))



########################################################
# STEP 2: Identify CNV carriers
########################################################

# Identify CNV carriers according to the following criteria; Distinguish deletion and duplication carriers;
# Start: 29.4-29.8 Mb; End: 30.05-30.4 Mb; |QS| > 0.5 
cnv_carriers <- cnvs[which(cnvs$Chromosome == 16 & cnvs$Start_Position_bp >= 29400000 & cnvs$Start_Position_bp <= 29800000 & cnvs$End_Position_bp >= 30050000 & cnvs$End_Position_bp <= 30400000), "Sample_Name"]

# High confidence CNV carriers
dup_carriers <- cnvs[which(cnvs$Copy_Number > 2), "Sample_Name"]
del_carriers <- cnvs[which(cnvs$Copy_Number < 2), "Sample_Name"]

# Split covariate file by CNV status
cov_ctrl <- cov[!cov$IID %in% cnv_carriers, ]
cov_dup <- cov[cov$IID %in% dup_carriers, ]
cov_del <- cov[cov$IID %in% del_carriers, ]

# Exclude CNV carriers with missing data & set a unique identifier for each carrier
cov_dup <- na.omit(cov_dup[,c(1:4,6:8)]) # 89 -> 61
cov_dup$ID <- seq(1:nrow(cov_dup))
cov_del <- na.omit(cov_del[,c(1:4,6:8)]) # 73 -> 58
cov_del$ID <- seq(1:nrow(cov_del))



########################################################
# STEP 3: Identify CNV matching controls
########################################################

### Duplication carriers ###############################

# Create an empty dataframe to store results
controls_dup <- data.frame()
controls_dup_count <- data.frame(IID = cov_dup$IID, ID = cov_dup$ID)

# Loop over duplication carriers
for (i in 1:nrow(cov_dup)) {

	# Start with all non-carriers
	df <- cov_ctrl

	# Select for age (+/- 2.5 years)
	age <- cov_dup[i, "age"]
	df <- df[which(df$age >= age - 2.5 & df$age <= age + 2.5), ]

	# Select for sex (identical)
	sex <- cov_dup[i, "sex"]
	df <- df[which(df$sex == sex), ]

	# Select for BMI (+/- 2.5 kg/m2)
	BMI <- cov_dup[i, "BMI"]
	df <- df[which(df$BMI >= BMI - 2.5 & df$BMI <= BMI + 2.5), ]

	# Select for TDI (+/- 2)
	TDI <- cov_dup[i, "TDI"]
	df <- df[which(df$TDI >= TDI - 2 & df$TDI <= TDI + 2), ]

	# Select for income (identical)
	income <- cov_dup[i, "income"]
	df <- df[which(df$income == income), ]

	# Select for EA (+/- 1)
	EA <- cov_dup[i, "age_end_edu"]
	df <- df[which(df$age_end_edu >= EA - 1 & df$age_end_edu <= EA + 1), ]

	# Add identified
	if(nrow(df) > 0) {df$ID <- i}
	
	# Assess the number of remaining controls & randomly select them
	print(paste0("Number of remaining controls for DUP #", i, ": ", nrow(df)))
	controls_dup_count[i, "count"] <- nrow(df)
	if (nrow(df) >= n) {
		# Sample
		df <- sample_n(df, n)
		# Merge dataframe
		controls_dup <- rbind(controls_dup, df)
	} else {print(paste0("Less than ", n, " samples remaining for DUP ", i, "! No selection performed."))}

}
rm(df, age, sex, BMI, TDI, income, EA)

print(paste0("Number of retained duplication carriers: ", length(unique(controls_dup$ID))))

# Save results
fwrite(controls_dup_count, paste0("/data/case_control/controls_dup_min_", n, "_count.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
fwrite(controls_dup, paste0("/data/case_control/controls_dup_min_", n, ".txt.gz"), col.names = T, row.names = F, quote = F, sep = "\t")


### Deletion carriers ##################################

# Create an empty dataframe to store results
controls_del <- data.frame()
controls_del_count <- data.frame(IID = cov_del$IID, ID = cov_del$ID)

# Loop over deletion carriers
for (i in 1:nrow(cov_del)) {

  # Start with all non-carriers
	df <- cov_ctrl

	# Select for age (+/- 2.5 years)
	age <- cov_del[i, "age"]
	df <- df[which(df$age >= age - 2.5 & df$age <= age + 2.5), ]

	# Select for sex (identical)
	sex <- cov_del[i, "sex"]
	df <- df[which(df$sex == sex), ]

	# Select for BMI (+/- 2.5 kg/m2)
	BMI <- cov_del[i, "BMI"]
	df <- df[which(df$BMI >= BMI - 2.5 & df$BMI <= BMI + 2.5), ]

	# Select for TDI (+/- 2)
	TDI <- cov_del[i, "TDI"]
	df <- df[which(df$TDI >= TDI - 2 & df$TDI <= TDI + 2), ]

	# Select for income (identical)
	income <- cov_del[i, "income"]
	df <- df[which(df$income == income), ]

	# Select for EA (EA +/- 1)
	EA <- cov_del[i, "age_end_edu"]
	df <- df[which(df$age_end_edu >= EA - 1 & df$age_end_edu <= EA + 1), ]

	# Add identified
	if(nrow(df) > 0) {df$ID <- i}
	
	# Assess the number of remaining controls & randomly select them
	print(paste0("Number of remaining controls for DEL #", i, ": ", nrow(df)))
	controls_del_count[i, "count"] <- nrow(df)
	if (nrow(df) >= n) {
		# Sample
		df <- sample_n(df, n)
		# Merge dataframe
		controls_del <- rbind(controls_del, df)
	} else {print(paste0("Less than ", n, " samples remaining for DEL ", i, "! No selection performed."))}

}
rm(df, age, sex, BMI, TDI, income, EA)

print(paste0("Number of retained deletion carriers: ", length(unique(controls_del$ID))))

# Save results
fwrite(controls_del_count, paste0("/data/case_control/controls_del_min_", n, "_count.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
fwrite(controls_del, paste0("/data/case_control/controls_del_min_", n,".txt.gz"), col.names = T, row.names = F, quote = F, sep = "\t")