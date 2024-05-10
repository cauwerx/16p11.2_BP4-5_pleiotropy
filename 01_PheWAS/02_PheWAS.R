########################################################
# Libraries
########################################################
library(dplyr)
library(data.table)
library(logistf)
library(MetBrewer)
library(ggplot2)
library(ggforestplot)
library(showtext)
font_add("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")



########################################################
# STEP 1: Load data
########################################################

# Load sample file (vector with identifiers, with header "eid")
samples <- as.data.frame(fread("samples_white_british_All.txt", col.names = "eid"))

# Load CNV file and subset samples (output from CNV calling pipeline, as described in Auwerx et al., 2022 (AJHG))
cnvs <- as.data.frame(fread("ukb_cnv_global.gz"))
cnvs$Sample_Name <- sub("_.*", "", cnvs$Sample_Name)
cnvs <- cnvs[cnvs$Sample_Name %in% samples$eid, ]

# Load continuous phenotype data (inverse normal transformed + covariate-corrected; as described in Auwerx et al., 2022 (AJHG); identifier + 1 column/trait); 
# Sex-specific traits originate from sex-specific data files and are appended to the main dataframe
pheno_All <- as.data.frame(fread("pheno_continuous_WB_INT_age_age2_sex_batch_PCs_All.txt.gz"))
pheno_M <- as.data.frame(fread("pheno_continuous_WB_INT_age_age2_batch_PCs_M.txt.gz"))
pheno_F <- as.data.frame(fread("pheno_continuous_WB_INT_age_age2_batch_PCs_F.txt.gz"))
pheno <- left_join(pheno_All, pheno_M[, names(pheno_M) %in% c("IID", "balding", "facial_hair")], by = "IID"); rm(pheno_All, pheno_M)
pheno <- left_join(pheno, pheno_F[, names(pheno_F) %in% c("IID", "menarche", "menopause", "birth_weight_first_child")], by = "IID"); rm(pheno_F)

# Load disease status data (identifier + 1 column/trait); 
# Sex-specific diseases originate from sex-specific data files and are appended to the main dataframe; 
# Data file is encoded so that cases are 1, controls are 0, and missing data NA
disease_All <- as.data.frame(fread("pheno_ICD10_All.txt.gz"))
disease_M <- as.data.frame(fread("pheno_ICD10_M.txt.gz"))
disease_F <- as.data.frame(fread("pheno_ICD10_F.txt.gz"))
disease <- left_join(disease_All, disease_M[, names(disease_M) %in% c("IID", "PC")], by = "IID"); rm(disease_All, disease_M)
disease <- left_join(disease, disease_F[, names(disease_F) %in% c("IID", "BC", "OC", "endometriosis", "menstruation")], by = "IID"); rm(disease_F)

# Load covariate data (identifier + 1 column/covariate, including age, sex, array, and PC1-40)
cov <- as.data.frame(fread("covariates.txt"))



########################################################
# STEP 2: Generate genotype table
########################################################

# Identify CNV carriers according to the following criteria; Distinguish deletion and duplication carriers;
# Start: 29.4-29.8 Mb; End: 30.05-30.4 Mb; |QS| > 0.5 
cnvs <- cnvs[which(cnvs$Chromosome == 16 & cnvs$Start_Position_bp >= 29400000 & cnvs$Start_Position_bp <= 29800000 & cnvs$End_Position_bp >= 30050000 & cnvs$End_Position_bp <= 30400000 & abs(cnvs$Quality_Score) >= 0.5), ]
dup_carriers <- cnvs[which(cnvs$Copy_Number > 2), "Sample_Name"]
del_carriers <- cnvs[which(cnvs$Copy_Number < 2), "Sample_Name"]

# Create an empty data frame to store CNV genotype encoding to 4 models
cnv_geno <- data.frame(IID = samples$eid, DUP = 0, DEL = 0, M = 0, U = 0)

# Duplication-only model encoding
cnv_geno[cnv_geno$IID %in% del_carriers, "DUP"] <- NA
cnv_geno[cnv_geno$IID %in% dup_carriers, "DUP"] <- 1

# Deletion-only model encoding
cnv_geno[cnv_geno$IID %in% dup_carriers, "DEL"] <- NA
cnv_geno[cnv_geno$IID %in% del_carriers, "DEL"] <- 1

# Mirror model encoding
cnv_geno[cnv_geno$IID %in% del_carriers, "M"] <- -1
cnv_geno[cnv_geno$IID %in% dup_carriers, "M"] <- 1

# U-shape model encoding
cnv_geno[cnv_geno$IID %in% c(dup_carriers, del_carriers), "U"] <- 1



########################################################
# STEP 3: Continuous traits PheWAS
########################################################

# Create a dataframe to store results
pheno_result <- data.frame(PHENO = rep(names(pheno)[names(pheno) != "IID"], each = 4),
                           MODEL = rep(c("DUP", "DEL", "M", "U"), length(names(pheno)[names(pheno) != "IID"])))

# Loop over phenotypes
for (p in unique(pheno_result$PHENO)) {
  
  print(paste0("Analyzing: ", p))
  
  # Loop over models
  for (m in c("DUP", "DEL", "M", "U")) {
    
      print(m)
    
      # Subset data (phenotype & model of interest)
      df <- left_join(na.omit(pheno[, c("IID", p)]), cnv_geno[, c("IID", m)], by = "IID")[2:3]
      colnames(df) <- c("PHENO", "CNV") 
      
      # Fit model
      fit <- lm(PHENO ~ CNV, data = df)
    
      # Fill in summary table
      pheno_result[which(pheno_result$PHENO == p & pheno_result$MODEL == m), "BETA"] <- summary(fit)$coeff[2,1]
      pheno_result[which(pheno_result$PHENO == p & pheno_result$MODEL == m), "SE"] <- summary(fit)$coeff[2,2]
      pheno_result[which(pheno_result$PHENO == p & pheno_result$MODEL == m), "P"] <- summary(fit)$coeff[2,4]
      pheno_result[which(pheno_result$PHENO == p & pheno_result$MODEL == m), "L95"] <- summary(fit)$coeff[2,1] - 1.96*summary(fit)$coeff[2,2]
      pheno_result[which(pheno_result$PHENO == p & pheno_result$MODEL == m), "U95"] <- summary(fit)$coeff[2,1] + 1.96*summary(fit)$coeff[2,2]
  
    }

}; rm(p, m, df, fit)

# Save data
fwrite(pheno_result, "/data/16p11.2_BP4-5_continuousPheWAS.txt", col.name = T, row.name = F, sep = "\t", quote = F)

# Identify the best model (i.e., yielding the lowest p-value) and select significant associations
top_pheno_result <- pheno_result %>% 
                    group_by(PHENO) %>% 
                    slice(which.min(P))
top_pheno_result_sig <- top_pheno_result[which(top_pheno_result$P <= 0.05/117), ]

# Number of associations remaining significant after accounting for the four association models
top_pheno_result_sig_model <- top_pheno_result[which(top_pheno_result$P <= 0.05/(4*117)), ]



########################################################
# STEP 4: Disease PheWAS
########################################################

# Create a dataframe to store results
diseases_result <- data.frame(PHENO = rep(names(disease)[names(disease) != "IID"], each = 4),
                           MODEL = rep(c("DUP", "DEL", "M", "U"), length(names(disease)[names(disease) != "IID"])))


# Loop over diseases
for (p in unique(diseases_result$PHENO)) {
  
  print(paste0("Analyzing: ", p))
  
  # Loop over models
  for (m in c("DUP", "DEL", "M", "U")) {
    
    print(m)
    
    # Subset data (disease & model of interest)
    df <- left_join(na.omit(disease[, c("IID", p)]), cnv_geno[, c("IID", m)], by = "IID")
    colnames(df) <- c("IID", "PHENO", "CNV") 
    
    # Add covariates (unlike continuous traits that were pre-corrected, we include age/sex/array/PC1-40 in the logistic regression model)
    df <- left_join(df, cov, by = "IID")
    
    # Fit model (Firth regression; do not include sex as a covariate for sex-specific diseases)
    if (p %in% c("PC","BC","OC","endometriosis","menstruation")) {
      fit_firth <- logistf(PHENO ~ ., data = df[-c(1,5)], plconf = 2, maxit = 100, maxstep = 10)
    } else {
      fit_firth <- logistf(PHENO ~ ., data = df[-c(1)], plconf = 2, maxit = 100, maxstep = 10)
    }
    
    # Fill in summary table
    diseases_result[which(diseases_result$PHENO == p & diseases_result$MODEL == m), "BETA"] <- fit_firth$coefficients[2]
    diseases_result[which(diseases_result$PHENO == p & diseases_result$MODEL == m), "SE"] <- sqrt(diag(vcov(fit_firth)))[2]
    diseases_result[which(diseases_result$PHENO == p & diseases_result$MODEL == m), "P"] <- fit_firth$prob[2]
    diseases_result[which(diseases_result$PHENO == p & diseases_result$MODEL == m), "L95"] <- fit_firth$ci.lower[2]
    diseases_result[which(diseases_result$PHENO == p & diseases_result$MODEL == m), "U95"] <- fit_firth$ci.upper[2]
    
  }
  
}
rm(p, m, df, fit_firth)

# Save data
fwrite(diseases_result, "/data/16p11.2_BP4-5_binaryPheWAS_firth.txt", col.name = T, row.name = F, sep = "\t", quote = F)

# Identify the best model (i.e., yielding the lowest p-value) and select significant associations
top_diseases_result <- diseases_result %>% 
                       group_by(PHENO) %>% 
                       slice(which.min(P))
top_diseases_result_sig <- top_diseases_result[which(top_diseases_result$P <= 0.05/117), ]

# Number of associations remaining significant after accounting for the four association models
top_diseases_result_sig_model <- top_diseases_result[which(top_diseases_result$P <= 0.05/(4*117)), ]



########################################################
# STEP 5: Plot the PheWAS
########################################################

# Combine results from continuous and binary PheWAS (DUP & DEL models) into a single dataframe for plotting
df <- rbind(pheno_result[pheno_result$MODEL %in% c("DUP", "DEL"), ], diseases_result[diseases_result$MODEL %in% c("DUP", "DEL"), ])

# Define color based on significance (different color for DUP vs DEL)
df$SIG <- 1
df[which(df$P < 0.05/(57+60) & df$MODEL == "DUP"), "SIG"] <- 2
df[which(df$P < 0.05/(57+60) & df$MODEL == "DEL"), "SIG"] <- 3
colors <- c("gray85",  "cornflowerblue", "indianred2")
df$SIG <- factor(df$SIG, levels = 1:3)

# Define shape based on model
df$MODEL <- factor(df$MODEL, levels = c("DEL", "DUP"))
shapes <- c(16, 15)

# Set category based on physiological system 
df[which(df$PHENO %in% c("IDA", "B12A", "AA", "anaemia", "sarcoidosis", "WBC_count", "reticulocyte_count", "RBC_count", "platelet_count", "neutrophil_count", "monocyte_count", "MCH", "lymphocyte_count", "eosinophil_count")), "CAT"] <- "hematologic"
df[which(df$PHENO %in% c("systolic_BP", "heart_rate", "diastolic_BP", "IHD", "aneurysm", "PE", "HTN_essential", "valves", "cardiomyopathies", "conduction", "stroke_ischemic", "stroke_hemorrhagic")), "CAT"] <- "cardiovascular"
df[which(df$PHENO %in% c("FVC", "pneumonia", "emphysema", "COPD", "asthma")), "CAT"] <- "pulmonary"
df[which(df$PHENO %in% c("GGT", "bilirubin", "AST", "ALT", "albumin", "total_protein", "ALP", "hepatic_fibrosis", "cholelithiasis", "bilirubin_disorder", "IBD", "celiac", "CRC")), "CAT"] <- "hepatic &\ngastrointestinal"
df[which(df$PHENO %in% c("urea", "urate", "cystatinC", "Cr", "AKI", "CKD", "KS", "PKD", "KC")), "CAT"] <- "renal"
df[which(df$PHENO %in% c("height", "GS", "vitamin_D", "phosphate", "Ca", "BMD", "IGF1", "osteoporosis", "mineral_disorder", "OA", "RA", "gout", "SLE", "hernia")), "CAT"] <- "musculoskeletal &\nconnective tissue"
df[which(df$PHENO %in% c("cornea", "cataract", "glaucoma")), "CAT"] <- "ophthalmic"
df[which(df$PHENO %in% c("psoriasis")), "CAT"] <- "dermatologic"
df[which(df$PHENO %in% c("BC", "OC", "PC", "endometriosis", "menstruation", "testosterone", "SHBG", "menopause", "menarche", "facial_hair", "birth_weight_first_child", "balding")), "CAT"] <- "sex-specific"
df[which(df$PHENO %in% c("SCZ", "bipolar", "depression", "neuroticism", "PD", "AD", "MS", "epilepsy", "headaches", "apnea", "sleep", "fluid_intelligence")), "CAT"] <- "neuropsychiatric"
df[which(df$PHENO %in% c("CRP", "hypothyroidism", "hyperthyroidism", "pituitary", "WHRadjBMI", "WHR", "weight", "body_fat_mass", "BMR", "BMI", "birthweight", "TG", "LPA", "LDL", "HDL", "cholesterol", "ApoA", "ApoB", "lipid", "glucose", "HbA1c", "T1D")), "CAT"] <- "metabolic &\nendocrine"
df$CAT <- factor(df$CAT, levels = c("neuropsychiatric", "ophthalmic", "dermatologic", "hematologic", "musculoskeletal &\nconnective tissue", "pulmonary", "cardiovascular", "metabolic &\nendocrine", "hepatic &\ngastrointestinal", "renal", "sex-specific"))

# Set label (file with two columns: PHENO = short phenotype name used in above analyses; LABEL = phenotype name to be printed on the graph)
labels <- as.data.frame(fread("phenotype_labels.txt"))
df <- left_join(df, labels, by = "PHENO")

# Order by increasing BETA (deletion effect)
order_pheno <- reshape(df[, c(10,2,3)], idvar = "LABEL", timevar = "MODEL", direction = "wide")
order_pheno$DIFF <- abs(order_pheno$BETA.DUP - order_pheno$BETA.DEL)
order_pheno <- order_pheno[order(order_pheno$DIFF, decreasing = F), ]
df$LABEL <- factor(df$LABEL, levels = order_pheno$LABEL)

# Indicate significant mirror and U-shape models
df_annot <- unique(df[, names(df) %in% c("PHENO", "CAT", "LABEL")])
df_annot$TEXT <- ""
df_annot[df_annot$PHENO %in% c(top_pheno_result_sig$PHENO[which(top_pheno_result_sig$MODEL == "M")], as.character(top_diseases_result_sig$PHENO[which(top_diseases_result_sig$MODEL == "M")])), "TEXT"] <- "M"
df_annot[df_annot$PHENO %in% c(top_pheno_result_sig$PHENO[which(top_pheno_result_sig$MODEL == "U")], as.character(top_diseases_result_sig$PHENO[which(top_diseases_result_sig$MODEL == "U")])), "TEXT"] <- "U"
df_annot$x_pos <- 4.25

# Plot 
p_pheWAS <- ggplot(data = df) +
            # Background
            geom_stripes(aes(y = LABEL), odd = "#F4F4F4", even = "white") +
            geom_vline(xintercept = 0, color = "#2B4050", linewidth = 0.8, lty = "dashed") +
            geom_vline(xintercept = c(-2,2,4), color = "#2B4050", linewidth = 0.2, lty = "dotted") +
            # Effect
            geom_effect(aes(x = BETA, y = LABEL, xmin = L95, xmax = U95, colour = SIG, shape = MODEL),
                        position = ggstance::position_dodgev(height = 0.5), size = 0.6) +
            scale_color_manual("Significance", values = colors, guide = "none") +
            scale_shape_manual("CNV", values = shapes, guide = "none") +
            # Facet
            facet_grid(CAT~., scales = "free", space = "free") +
            # Annotation
            geom_text(data = df_annot, aes(x = x_pos, y = LABEL, label = TEXT), size = 2, color = "#2B4050") + 
            # Axis
            xlab("Effect ± 95% CI") + ylab("") +
            coord_cartesian(xlim = c(-2.5, 4.23)) +
            # Layout
            theme_bw() +
            theme(text = element_text(family = "Arial"),
                  axis.text.x = element_text(size = 11),
                  axis.title.x = element_text(size = 12, margin = margin(t = 5, r = 0, b = 0, l = 0)),
                  axis.text.y = element_text(size = 9),
                  plot.title = element_text(size = 14, face = "bold", color = "#2B4050"),
                  strip.text.y = element_text(size = 9, color = "#2B4050", face = "bold",  angle = 0),
                  strip.background = element_rect(fill = "#F4F4F4", color = "#2B4050"))

# Save plot        
showtext_auto()   
p_pheWAS
ggsave("/plot/16p11.2_PheWAS.pdf", width = 6, height = 15)
showtext_auto(FALSE)