###
predict.data=readRDS('~/processed_data/master_processed_prediction_data.RDS')
library(lubridate)
attach(predict.data)
# Define baseline specimen IDs
subject_specimen_BL <- subject_specimen %>%
  filter(planned_day_relative_to_boost == 0)
specimen_id_BL <- subject_specimen_BL %>%
  pull(specimen_id)

# Extract & Transform data
extract_transform_processed_data <- function(omic_list, specimen_id_BL) {
  # Reshape into data.frame with specimen.id & omic features as columns
  X <- omic_list$processed_similar_to_training %>% t() %>% data.frame()
  X$specimen_id <- as.numeric(rownames(X))
  rownames(X) <- NULL
  # Filter down to BL specimens
  X <- X %>% filter(specimen_id %in% specimen_id_BL)
  return(X)
}

X_list <- lapply(predict.data[2:length(predict.data)], extract_transform_processed_data, specimen_id_BL)

# Checkout the BL missingness rates
n_subjects <- subject_specimen$subject_id %>% unique %>% length()

# n observations at BL in each dataset
lapply(X_list, nrow)

# proportion of observations missing at BL
lapply(X_list, function(X, n_subjects) { (n_subjects - nrow(X)) / n_subjects }, n_subjects)

# Collapse into one data.frame
X_combined <- subject_specimen_BL %>% 
  left_join(X_list$abtiter) %>% 
  left_join(X_list$plasma_cytokine_concentrations) %>%
  left_join(X_list$pbmc_cell_frequency) %>%
  left_join(X_list$pbmc_gene_expression)


# Combine and format data for imputation
data_combined <- X_combined %>%
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
data_combined_predict=data_combined
saveRDS(data_combined_predict, "~/processed_data/data_imputed_predict.RDS")

