# Clear existing data and graphics
rm(list=ls())
graphics.off()

# Load libraries -----
library(Hmisc)
library(tidyverse)
library(dplyr)
library(rms)
library(lubridate)
library(ResourceSelection)
library(caret)
library(pROC)
library(mice)
library(skimr)
library(ggplot2)
library(GGally)
library(ggbeeswarm)
library(crosstable)
library(flextable)

library(broom)
library(pander)
library(gridExtra)
library(grid)
library(sjPlot)
library(sjmisc)
library(sjlabelled)

# Load data -----
data <- read.csv("Y:/Research/Methodology/IPD_DVT/E_ResearchData/2_ResearchData/data_ipd_dvt_general.csv")

# Load imputed data ----
load("~/Desktop/imp.mlmi.RData")
data_imp <- complete(imp.mlmi, "long", include = FALSE) #TRUE still has some missing values otherwise FALSE
# write.csv(data_imp,"~/Desktop/Data1.csv", row.names = FALSE) #to get the data in csv format

# Getting to know the data
str(data_imp) 
skim(data_imp)

# Training and testing data set ----
data_train <- data_imp[data_imp$.imp == 1, ] #training
data_test <- data_imp[data_imp$.imp == 2, ] #testing

# Select variable for train data
data_train = data_train %>% select(.id, #patient id
                                   studyid, #study cluster
                                   age, 
                                   sex, #female is coded as 0
                                   ddimdich, #d-dimer dichotomic 
                                   malign, #cancer
                                   hist, #previous history of DVT
                                   altdiagn, #alternative diagnostic
                                   dvt) 
# Select variable for test data
data_test = data_test %>% select(.id, #patient id
                                 studyid, #study cluster
                                 age, 
                                 sex, #female is coded as 0
                                 ddimdich, #d-dimer dichotomic 
                                 malign, #cancer
                                 hist, #previous history of DVT
                                 altdiagn, #alternative diagnostic
                                 dvt) 

# Descriptive Table 1 ----
dt1 = crosstable(data_train, c(age, sex, ddimdich, malign, hist, altdiagn), by=dvt, total="both", 
                 #funs = "Mean(std)",
                 percent_pattern="{n} ({p_col})", percent_digits=1) %>%
  as_flextable()
dt1

# Export to word 
tf <- tempfile(fileext = ".docx")
save_as_docx("Table 1"= dt1, path = ("/Users/candaurn/Desktop/Absolute risk/Table_1.docx"))

# Sample size calculation for logistic regression ----
library(pmsampsize)
samplesize <- pmsampsize(type="b", cstatistic = 0.89, parameters = 6, prevalence = 0.18, seed=123)
samplesize #227

# Cross-validation methods ----
fitControl <- trainControl(method = "repeatedcv", # Cross-validation, default is bootstrap
                           number = 10,
                           repeats = 3,
                           classProbs = TRUE,
                           savePredictions = TRUE,
                           index = createFolds(data_train$dvt,5),
                           summaryFunction = twoClassSummary) # ROC, sens, spec under 50% cut-off

# Fit logistic regression ----
set.seed(825)
glm_mod <- train(dvt ~ age + sex + ddimdich  + hist + malign + altdiagn,
                 data = data_train,
                 method = "glm",
                 family = binomial (link = 'logit'),
                 trControl = fitControl,
                 metric = "ROC")

getTrainPerf(glm_mod) 
glm_mod$results # AUC, Sens, Spec with SD
glm_mod$finalModel
summary(glm_mod)

## Predictions with test data ----
glm.pred <- predict(glm_mod,
                    newdata = data_test)

## Probabilities with test data ----
glm.probs <- predict(glm_mod,
                     newdata = data_test,
                     type = "prob",  
                     se.fit = TRUE) 
head(glm.probs)

## Confusion matrix ----
confusionMatrix(data = glm.pred, 
                reference = data_test$dvt,
                positive = "Pos")

# Fit Random Forest ---- 
set.seed(825)
rf_mod.1 <- train(dvt ~ age + sex + ddimdich  + hist + malign + altdiagn, 
                  data = data_train,
                  method = "ranger",
                  trControl = fitControl,
                  #importance = TRUE,
                  metric = "ROC", # to optimize the model based on ROC
                  #ntree = 50,
)
rf_mod.1$results
rf_mod.1$finalModel

set.seed(825)
rf_mod.2 <- train(dvt ~ age + sex + ddimdich  + hist + malign + altdiagn, 
                  data = data_train,
                  method = "rf",
                  trControl = fitControl,
                  #importance = TRUE,
                  metric = "ROC", # to optimize the model based on ROC
                  #ntree = 50,
)
rf_mod.2$results
rf_mod.2$finalModel

## Predict on test data ----
ranger.pred <- predict(rf_mod.1,  # for confusion matrix
                       data_test,
                       type = "raw") 

ranger.probs <- predict(rf_mod.1,
                        data_test,
                        type = "prob")

rf.pred <- predict(rf_mod.2,
                   data_test,
                   type = "raw") 

rf.probs <- predict(rf_mod.2,
                    data_test,
                    type = "prob")

## Confusion matrix ----
confusionMatrix(data = ranger.pred,
                reference = data_test$dvt,
                #mode = "prec_recall"
                mode = "everything",
                positive = "Pos") 

confusionMatrix(data = rf.pred,
                reference = data_test$dvt,
                #mode = "prec_recall"
                mode = "everything",
                positive = "Pos") 

# Fit linear SVM ----
set.seed(825)
svm_mod.1 <- train(dvt ~ age + sex + ddimdich + hist + malign + altdiagn,
                   data = data_train,
                   method = "svmLinear", #kernlab
                   trControl = fitControl,
                   importance = TRUE)

svm_mod.1$results

svm_mod.2 <- train(dvt ~ age + sex + ddimdich + hist + malign + altdiagn,
                   data = data_train,
                   method = "svmLinear2", #e1071
                   trControl = fitControl,
                   importance = TRUE)
svm_mod.2$results

## Prediction on test data ----
kernlab.pred <- predict(svm_mod.1,
                        data_test,
                        type= "raw")

kernlab.probs <- predict(svm_mod.1, 
                         data_test,
                         type = "prob") 

e1071.pred <- predict(svm_mod.2,
                      data_test,
                      type= "raw")

e1071.probs <- predict(svm_mod.2, 
                       data_test,
                       type = "prob") 

## Confusion matrix ---- 
confusionMatrix(data = kernlab.pred, 
                reference = data_test$dvt,
                #mode = "prec_recall"
                mode = "everything",
                positive = "Pos") 

confusionMatrix(data = e1071.pred, 
                reference = data_test$dvt,
                #mode = "prec_recall"
                mode = "everything",
                positive = "Pos") 
