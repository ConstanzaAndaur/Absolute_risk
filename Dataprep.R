###### load data
load("/Users/ymaerziy/Networkshares/datamanagement/Research/Methodology/IPD_DVT/E_ResearchData/2_ResearchData/imp.mlmi.RData")

#TRUE still has some missing values otherwise FALSE
data_imp <- complete(imp.mlmi, "long", include = FALSE) 

# load packages
library(mice)
library(dplyr)
library(flextable)
library(crosstable)
library(caret)
library(skimr)
library(pROC)
library(purrr)
library(tidyr)
library(GGally)
 
# subset .imp=1
data_imp_subset <- subset(data_imp, .imp == 1)

# Check
table(data_imp_subset$.imp)
table(data_imp_subset$studyid)

# Splitting training and testing data set ----
set.seed(123)

train_ids <- data_imp_subset$studyid %in% c(1, 2, 3, 4, 5, 6, 7, 8)
train_data <- data_imp_subset %>% filter(train_ids)  # 1,2,3,4,5,6,7,8
#test_data  <- data_imp_subset %>% filter(!train_ids) # 9,10,11,12,13

table(train_data$studyid)

# Select variable for train data
train_data = train_data %>% select(.id, #patient id
                                   studyid, #study cluster
                                   age, 
                                   sex, #female is coded as 0
                                   ddimdich, #d-dimer dichotomic 
                                   malign, #cancer
                                   hist, #previous history of DVT
                                   altdiagn, #alternative diagnostic
                                   dvt) 

train_data <- train_data %>%
  dplyr::select(.id, studyid, age, sex, ddimdich, malign, hist, altdiagn, dvt)

# Select variable for test data
test_data = test_data %>% select(.id, #patient id
                                 studyid, #study cluster
                                 age, 
                                 sex, #female is coded as 0
                                 ddimdich, #d-dimer dichotomic 
                                 malign, #cancer
                                 hist, #previous history of DVT
                                 altdiagn, #alternative diagnostic
                                 dvt) 


# Descriptive Table 1 ----
dt1 = crosstable(train_data, c(age, sex, ddimdich, malign, hist, altdiagn), by=dvt, total="both", 
                 #funs = "Mean(std)",
                 percent_pattern="{n} ({p_col})", percent_digits=1) %>%
  as_flextable()
dt1


# Outcome ----
train_data$dvt = as.character(train_data$dvt)
train_data$dvt[train_data$dvt == "0"] = "Neg"
train_data$dvt[train_data$dvt == "1"] = "Pos"
train_data$dvt = factor(train_data$dvt, levels = c('Neg', 'Pos'))

train_data %>% 
  dplyr::count(dvt) #1146 cases of DVT / 4941 non-cases  1146/(1146+4941)=0.18 prevalence

#.  dvt    n
#1 Neg.  4941
#2 Pos   1146


# Sample size calculation for logistic regression ----
library(pmsampsize)
samplesize <- pmsampsize(type="b", cstatistic = 0.89, parameters = 6, prevalence = 0.18, seed=123)
samplesize #227  with 41 events (assuming an outcome prevalence = 0.18) and an EPP = 6.81 

samplesize <- pmsampsize(type="b", cstatistic = 0.81, parameters = 6, prevalence = 0.18, seed=123)
samplesize #266   with 48 events (assuming an outcome prevalence = 0.18) and an EPP = 7.98 

# Cox-Snell R-sq = 0.1818  
# Shrinkage 0.9
