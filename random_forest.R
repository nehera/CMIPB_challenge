# This script uses random forest regression to predict Y1 and Y2
# Author: Aidan
# Date Initialized: December 21, 2023

data_imputed <- readRDS("processed_data/data_imputed.RDS")
Y1 <- data_imputed$Y1
Y2 <- data_imputed$Y2
rf_predictors <- data_imputed %>%
  select(-subject_id, -Y1, -Y2)

# Install and load the randomForest package if you haven't already
# install.packages("randomForest")
library(randomForest)

set.seed(123)  # For reproducibility

# Fit the random forest regression model

# TODO test a grid of ntree & mtry
rf_model_1 <- randomForest(Y1 ~ ., data = rf_predictors, ntree = 500, mtry = 2)

rf_model_2 <- randomForest(Y2 ~ ., data = rf_predictors, ntree = 500, mtry = 2)

# Evaluate train set performance
predictions_1 <- predict(rf_model_1, newdata = rf_predictors)
spearman_corr_1_train <- cor(Y1, predictions_1, method = "spearman")

predictions_2 <- predict(rf_model_2, newdata = rf_predictors)
spearman_corr_2_train <- cor(Y2, predictions_2, method = "spearman")

# TODO Evaluate test set performance

# Describe variable importance
var_importance_1 <- importance(rf_model_1) %>% data.frame()
var_importance_1$var <- rownames(var_importance_1)
rownames(var_importance_1) <- NULL

# Order node purity values from largest to smallest
ordered_indices <- order(var_importance_1$IncNodePurity, decreasing = TRUE)
ordered_var_importance_1 <- var_importance_1[ordered_indices, ]

head(ordered_var_importance_1)

# TODO Describe variable importance to Y2

#test and training split
dim(Y1)
#half and half data split
train_ind=sample(c(1:96)[!is.na(Y1)],45)
test_ind=c(1:96)[-train_ind]
test_Y1=Y1[test_ind]
test_X=rf_predictors[test_ind,]
train_Y1=Y1[train_ind]
train_X=rf_predictors[train_ind,]

rf_predictors <- data_imputed %>%
  select(-subject_id, Y1, Y2)
rf_model_train=randomForest(train_Y1 ~ ., data = train_X, ntree = 500, mtry = 2)
predictions_train <- predict(rf_model_train, newdata = train_X)
spearman_corr_1_train <- cor(train_Y1, predictions_train, method = "spearman")
predictions_test <- predict(rf_model_train, newdata = test_X)
spearman_corr_1_test <- cor(test_Y1, predictions_test, method = "spearman",use="pairwise.complete.obs")

####5-fold cross-Validation
set.seed(2024)
Folds=createFolds(c(1:96),k=5)
save(Folds,file='~/Folds.rda')

##
Test_Spearman=c()
for(i in 1:5){
  test_Y1=Y1[Folds[[i]]]
  test_X=rf_predictors[Folds[[i]],]

  train_Y1=Y1[-Folds[[i]]]
  train_X=rf_predictors[-Folds[[i]],]

rf_model_train=randomForest(y=train_Y1,x = train_X, ntree = 500, mtry = 2)
predictions_train <- predict(rf_model_train, newdata = train_X)
spearman_corr_1_train <- cor(train_Y1, predictions_train, method = "spearman")
predictions_test <- predict(rf_model_train, newdata = test_X)
Test_Spearman[i] <- cor(test_Y1, predictions_test, method = "spearman",use="pairwise.complete.obs")
}
mean(Test_Spearman)


### now, I am testing the blockForest method proposed in the following paper:
# https://doi-org.ezp2.lib.umn.edu/10.1186/s12859-019-2942-y
# install.packages("blockForest")
library(blockForest)
# for now, remove any non-numeric values. factors might be one-hot-encoded
# TODO make at least infancy_vac binary & year_of_birth numeric (relative to a baseline year)
numeric_columns <- data_imputed[sapply(data_imputed, is.numeric)] %>%
  select(-subject_id, -Y1, -Y2)

# Add days at boost
numeric_columns$days_at_boost <- difftime(data_imputed$date_of_boost,
         data_imputed$year_of_birth, units='days') %>% as.numeric
# Add infancy_vac as binary
numeric_columns$infancy_vac <- data_imputed$infancy_vac %>% as.numeric -1

# there are clinical covariates and four omics blocks
clinical_covar <- 1 # only actual_day_relative_to_boost numeric
columns_starting_with_ENSG <- grep("^ENSG", colnames(numeric_columns))
columns_starting_with_IgG <- grep("^IgG", colnames(numeric_columns))
protein_block <- 29:58 # plasma cytokine
immuno_block <- 59:78

blockfor_fit_1 <- blockfor(X = numeric_columns, y = Y1, blocks = list(clinical_covar,
                             columns_starting_with_ENSG,
                             columns_starting_with_IgG,
                             protein_block, 
                             immuno_block))

pred <- predict(blockfor_fit_1$forest, data = as.data.frame(numeric_columns))

####5-fold cross-Validation
set.seed(2024)
Folds=createFolds(c(1:96),k=5)
save(Folds,file='~/Folds.rda')

##
Test_Spearman=c()
for(i in 1:5){
  test_Y1=Y1[Folds[[i]]]
  test_X=numeric_columns[Folds[[i]],]
  
  train_Y1=Y1[-Folds[[i]]]
  train_X=numeric_columns[-Folds[[i]],]
  
  model_train= blockfor(X = train_X, y = train_Y1, blocks = list(clinical_covar,
                                                                   columns_starting_with_ENSG,
                                                                   columns_starting_with_IgG,
                                                                   protein_block, 
                                                                   immuno_block))
  pred <- predict(model_train$forest, data = train_X)
  predictions_train <- pred$predictions
  spearman_corr_1_train <- cor(train_Y1, predictions_train, method = "spearman")
  pred <- predict(model_train$forest, data = test_X)
  predictions_test <- pred$predictions
  Test_Spearman[i] <- cor(test_Y1, predictions_test, method = "spearman",use="pairwise.complete.obs")
}
mean(Test_Spearman)

Test_Spearman=c()
for(i in 1:5){
  test_Y2=Y2[Folds[[i]]]
  test_X=numeric_columns[Folds[[i]],]
  
  train_Y2=Y2[-Folds[[i]]]
  train_X=numeric_columns[-Folds[[i]],]
  
  model_train= blockfor(X = train_X, y = train_Y2, blocks = list(clinical_covar,
                                                                 columns_starting_with_ENSG,
                                                                 columns_starting_with_IgG,
                                                                 protein_block, 
                                                                 immuno_block))
  pred <- predict(model_train$forest, data = train_X)
  predictions_train <- pred$predictions
  spearman_corr_1_train <- cor(train_Y2, predictions_train, method = "spearman")
  pred <- predict(model_train$forest, data = test_X)
  predictions_test <- pred$predictions
  Test_Spearman[i] <- cor(test_Y2, predictions_test, method = "spearman",use="pairwise.complete.obs")
}
mean(Test_Spearman)

# and comparing it to BIP from https://github.com/chekouo/BIPnet
# TODO I ran into an error installing devtools onto the server... 
# install.packages("devtools")
# library(devtools)
# install_github('chekouo/BIPnet')
# library(BIPnet)