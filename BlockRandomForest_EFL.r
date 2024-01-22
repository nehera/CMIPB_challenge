#Perform block random forest for all 6 challenges, and predict on test set

data_imputed <- readRDS("~/processed_data/data_imputed_EFL.RDS")
Y1.1 <- data_imputed$Y1.1
Y1.2 <- data_imputed$Y1.2
Y2.1 <- data_imputed$Y2.1
Y2.2 <- data_imputed$Y2.2
Y3.1 <- data_imputed$Y3.1
Y3.2 <- data_imputed$Y3.2


### now, I am testing the blockForest method proposed in the following paper:
# https://doi-org.ezp2.lib.umn.edu/10.1186/s12859-019-2942-y
# install.packages("blockForest")
library(blockForest)
# for now, remove any non-numeric values. factors might be one-hot-encoded
# TODO make at least infancy_vac binary & year_of_birth numeric (relative to a baseline year)

numeric_columns <- data_imputed[sapply(data_imputed, is.numeric)] %>%
  select(-subject_id, -Y1.1, -Y1.2,-Y2.1,-Y2.2,-Y3.1,-Y3.2)

# Add days at boost
numeric_columns$days_at_boost <- difftime(data_imputed$date_of_boost,
                                          data_imputed$year_of_birth, units='days') %>% as.numeric
# Add infancy_vac as binary
numeric_columns$infancy_vac <- data_imputed$infancy_vac %>% as.numeric -1

data_imputed_test=readRDS('~/processed_data/data_imputed_predict.RDS')
numeric_columns_test <- data_imputed_test[sapply(data_imputed_test, is.numeric)] %>%
  select(-subject_id) 
# Add days at boost
numeric_columns_test$days_at_boost <- difftime(data_imputed_test$date_of_boost,
                                          data_imputed_test$year_of_birth, units='days') %>% as.numeric
# Add infancy_vac as binary
numeric_columns_test$infancy_vac <- data_imputed_test$infancy_vac %>% as.numeric -1
# Match=match(colnames(numeric_columns),colnames(numeric_columns_test))
# numeric_columns_test=data_imputed_test[,Match]

# there are clinical covariates and four omics blocks
clinical_covar <- c(1,8321,8322) # only actual_day_relative_to_boost numeric
columns_starting_with_ENSG <- grep("^ENSG", colnames(numeric_columns))
columns_starting_with_IgG <- grep("^IgG", colnames(numeric_columns))
protein_block <- 29:58 # plasma cytokine
immuno_block <- 59:78
blocks.list=list(clinical_covar,columns_starting_with_ENSG,
                 columns_starting_with_IgG,
                 protein_block, 
                 immuno_block)
#log transform Y3.1: 
Y3.1=log(Y3.1+1)

##fit model on full training data for all outcomes
blockfor_fit_1.1 <- blockfor(X = numeric_columns, y = Y1.1, blocks = blocks.list)
blockfor_fit_1.2 <- blockfor(X = numeric_columns, y = Y1.2, blocks = blocks.list)
blockfor_fit_2.1 <- blockfor(X = numeric_columns, y = Y2.1, blocks = blocks.list)
blockfor_fit_2.2 <- blockfor(X = numeric_columns, y = Y2.2, blocks = blocks.list)
blockfor_fit_3.1 <- blockfor(X = numeric_columns, y = Y3.1, blocks = blocks.list)
blockfor_fit_3.2 <- blockfor(X = numeric_columns, y = Y3.2, blocks = blocks.list)

#predict model on test data for all outcomes
library(blockForest)
pred.test1.1 <- predict(blockfor_fit_1.1$forest, data = as.data.frame(numeric_columns_test))$predictions
pred.test1.2 <- predict(blockfor_fit_1.2$forest, data = as.data.frame(numeric_columns_test))$predictions
pred.test2.1 <- predict(blockfor_fit_2.1$forest, data = as.data.frame(numeric_columns_test))$predictions
pred.test2.2 <- predict(blockfor_fit_2.2$forest, data = as.data.frame(numeric_columns_test))$predictions
pred.test3.1 <- predict(blockfor_fit_3.1$forest, data = as.data.frame(numeric_columns_test))$predictions
pred.test3.2 <- predict(blockfor_fit_3.2$forest, data = as.data.frame(numeric_columns_test))$predictions
#get orders
Ranking_Table=cbind(data_imputed_test$subject_id,order(pred.test1.1),order(pred.test1.2),order(pred.test2.1),order(pred.test2.2),order(pred.test3.1),order(pred.test3.2))
Ranking_Table
write.csv(Ranking_Table, file= "2ndChallengeSubmission.csv")





