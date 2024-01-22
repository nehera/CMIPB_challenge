# This script describes missingness at the patient level for:
# each X matrix (titers, cell frequency, gene expression, cytokines) & Y (IgG_PT at 14 days)
# Author: Aidan
# Date Initialized: December 20, 2023

# Load libraries
# TODO fix issue installing tidyverse & mice
# library(tidyverse)
# library(mice)
library(lubridate)

# Load "processed" data provided by CMI-PB
p.data <- readRDS('~/processed_data/master_processed_training_data.RDS')
attach(p.data)

# Define baseline specimen IDs
subject_specimen_BL <- subject_specimen %>%
  filter(planned_day_relative_to_boost == 0)
specimen_id_BL <- subject_specimen_BL %>%
  pull(specimen_id)

# Extract & Transform batchCorrected_data
extract_transform_batchCorrected_data <- function(omic_list, specimen_id_BL) {
  # Reshape into data.frame with specimen.id & omic features as columns
  X <- omic_list$batchCorrected_data %>% t() %>% data.frame()
  X$specimen_id <- as.numeric(rownames(X))
  rownames(X) <- NULL
  # Filter down to BL specimens
  X <- X %>% filter(specimen_id %in% specimen_id_BL)
  return(X)
}

X_list <- lapply(p.data[2:length(p.data)], extract_transform_batchCorrected_data, specimen_id_BL)

# Checkout the BL missingness rates
n_subjects <- subject_specimen$subject_id %>% unique %>% length()

# n observations at BL in each dataset
lapply(X_list, nrow)

# proportion of observations missing at BL
lapply(X_list, function(X, n_subjects) { (n_subjects - nrow(X)) / n_subjects }, n_subjects)

# Collapse into one data.frame
X_combined <- subject_specimen_BL %>% 
  left_join(X_list$abtiter_wide) %>% 
  left_join(X_list$plasma_cytokine_concentrations) %>%
  left_join(X_list$pbmc_cell_frequency) %>%
  left_join(X_list$pbmc_gene_expression)

# Extract outcome Y for task 1
subject_specimen_14 <- subject_specimen %>%
  filter(planned_day_relative_to_boost == 14)
specimen_id_14 <- subject_specimen_14 %>%
  pull(specimen_id)

# Reshape into data.frame with specimen.id & omic features as columns
Y1 <- p.data$abtiter_wide$batchCorrected_data %>% t() %>% data.frame()
Y1$specimen_id <- as.numeric(rownames(Y1))
rownames(Y1) <- NULL
Y1.1 <- Y1 %>% filter(specimen_id %in% specimen_id_14) %>%
  left_join(subject_specimen) %>%
  select(subject_id, IgG_PT) %>%
  rename(IgG_PT_14 = IgG_PT)
Y1.BL <- Y1 %>% filter(specimen_id %in% specimen_id_BL) %>%
  left_join(subject_specimen) %>%
  select(subject_id, IgG_PT) %>%
  rename(IgG_PT_BL = IgG_PT)


Y1_tidy <- merge(Y1.BL, Y1.1) %>% 
  mutate(Y1.2 = (IgG_PT_14+0.3)/ (IgG_PT_BL+0.3)) %>%  ####Shift by constant to avoid negatives
  rename(Y1.1 = IgG_PT_14) %>%
  select(subject_id, Y1.1, Y1.2)


# Extract outcome Y for task 2
subject_specimen_1 <- subject_specimen %>%
  filter(planned_day_relative_to_boost == 1)
specimen_id_1 <- subject_specimen_1 %>%
  pull(specimen_id)
Y2 <- p.data$pbmc_cell_frequency$batchCorrected_data %>% t() %>% data.frame()
Y2$specimen_id <- as.numeric(rownames(Y2))
rownames(Y2) <- NULL
Y2.1 <- Y2 %>% filter(specimen_id %in% specimen_id_1) %>%
  left_join(subject_specimen) %>%
  select(subject_id, Monocytes) %>%
  rename(Monocytes_1 = Monocytes)
Y2.BL <- Y2 %>% filter(specimen_id %in% specimen_id_BL) %>%
  left_join(subject_specimen) %>%
  select(subject_id, Monocytes) %>%
  rename(Monocytes_BL = Monocytes)

Y2_tidy <- merge(Y2.BL, Y2.1) %>% 
  mutate(Y2.2 = Monocytes_1/ Monocytes_BL) %>%
  rename(Y2.1 = Monocytes_1) %>%
  select(subject_id, Y2.1, Y2.2)

# Extract outcome Y for task 3 (gene CCL3 has ID ENSG00000277632.1)
subject_specimen_3 <- subject_specimen %>%
  filter(planned_day_relative_to_boost == 3)
specimen_id_3 <- subject_specimen_3 %>%
  pull(specimen_id)
Y3 <- p.data$pbmc_gene_expression$batchCorrected_data %>% t() %>% data.frame()
Y3$specimen_id <- as.numeric(rownames(Y3))
rownames(Y3) <- NULL
Y3.1 <- Y3 %>% filter(specimen_id %in% specimen_id_3) %>%
  left_join(subject_specimen) %>%
  select(subject_id, ENSG00000277632.1) %>%
  rename(CCL3_3 = ENSG00000277632.1)
Y3.BL <- Y3 %>% filter(specimen_id %in% specimen_id_BL) %>%
  left_join(subject_specimen) %>%
  select(subject_id, ENSG00000277632.1) %>%
  rename(CCL3_BL = ENSG00000277632.1)

Y3_tidy <- merge(Y3.BL, Y3.1) %>% 
  mutate(Y3.2 = CCL3_3/ CCL3_BL) %>%
  rename(Y3.1 = CCL3_3) %>%
  select(subject_id, Y3.1, Y3.2)

# Combine and format data for imputation
data_combined <- left_join(left_join(left_join(X_combined, Y1_tidy),Y2_tidy),Y3_tidy) %>%
  select(-c(specimen_id, planned_day_relative_to_boost, specimen_type,
            visit, timepoint, dataset)) %>% # dataset is not meaningful for test data
  mutate(infancy_vac = factor(infancy_vac), 
         biological_sex = factor(biological_sex),
         ethnicity = factor(ethnicity),
         race = factor(race)) %>%
  mutate(date_of_boost = ymd(date_of_boost),
         year_of_birth = ymd(year_of_birth)) 

# For now, we perform mean imputation

# Columns with NA values
cols_with_NA <- colnames(data_combined)[colSums(is.na(data_combined)) > 0]

# Impute NA values with column means
for (col in cols_with_NA) {
  col_mean <- mean(data_combined[[col]], na.rm = TRUE)
  data_combined[[col]][is.na(data_combined[[col]])] <- col_mean
}

saveRDS(data_combined, "~/processed_data/data_imputed_EFL.RDS")

