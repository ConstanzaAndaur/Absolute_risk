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