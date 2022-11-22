# Absolute risk in DVT
# Constanza L. Andaur Navarro

# Clear existing data and graphics
rm(list=ls())
graphics.off()

# Load required packages -----
packages <- c("Hmisc", "tidyverse", "dplyr", "rms", "lubridate", "ResourceSelection", "caret", "pROC", 
              "mice", "skimr", "ggplot2", "GGally", "ggbeeswarm", "crosstable", "flextable", "broom", 
              "pander", "gridExtra", "grid", "sjPlot", "sjmisc", "sjlabelled", "data.table", "GGally")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Load libraries ----
lapply(packages, library, character.only = TRUE) 

# Load data -----
data <- read.csv("Y:/Research/Methodology/IPD_DVT/E_ResearchData/2_ResearchData/data_ipd_dvt_general.csv")

# Load imputed data ----
load("~/Desktop/imp.mlmi.RData")
data_imp <- complete(imp.mlmi, "long", include = FALSE) #TRUE still has some missing values otherwise FALSE
# write.csv(data_imp,"~/Desktop/Data1.csv", row.names = FALSE) #to get the data in csv format

# Getting to know the data
str(data_imp) 
skim(data_imp)

# Splitting training and testing data set ----
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

# Outcome ----
data_train$dvt = as.character(data_train$dvt)
data_train$dvt[data_train$dvt == "0"] = "Neg"
data_train$dvt[data_train$dvt == "1"] = "Pos"
data_train$dvt = factor(data_train$dvt, levels = c('Neg', 'Pos'))

data_train %>% 
  dplyr::count(dvt) #1864 cases of DVT / 8138 non-cases 0.18 prevalence

data_test$dvt = as.character(data_test$dvt)
data_test$dvt[data_test$dvt == "0"] = "Neg"
data_test$dvt[data_test$dvt == "1"] = "Pos"
data_test$dvt = factor(data_test$dvt, levels = c('Neg', 'Pos'))

# Cross-validation ----
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

# Evaluation ----

## AUC plot for all five models ----

rocobj1 <- plot.roc(as.numeric(glm_mod$trainingData$.outcome=='Pos'),
                    aggregate(Pos~rowIndex,glm_mod$pred,mean)[,'Pos'], 
                    print.auc=TRUE)

rocobj2 <- lines.roc(as.numeric(rf_mod.1$trainingData$.outcome=='Pos'),
                     aggregate(Pos~rowIndex,rf_mod.1$pred,mean)[,'Pos'])

rocobj3 <- lines.roc(as.numeric(rf_mod.2$trainingData$.outcome=='Pos'),
                     aggregate(Pos~rowIndex,rf_mod.2$pred,mean)[,'Pos'])

rocobj4 <- lines.roc(as.numeric(svm_mod.1$trainingData$.outcome=='Pos'),
                     aggregate(Pos~rowIndex,svm_mod.1$pred,mean)[,'Pos'])

rocobj5 <- lines.roc(as.numeric(svm_mod.2$trainingData$.outcome=='Pos'),
                     aggregate(Pos~rowIndex,svm_mod.2$pred,mean)[,'Pos'])

roc.list <- list(glm = rocobj1, rf1 = rocobj2, rf2 = rocobj3, 
                 svm1 = rocobj4, svm2 = rocobj5)
# extract AUC
roc.list %>% 
  map(~tibble(AUC = .x$auc)) %>% 
  bind_rows(.id = "name") -> data.auc

roc.list[["svm2"]][["ci"]] #to find CI

# generate labels
data.auc %>% 
  mutate(label_long=paste0(name,", AUC = ", paste(round(AUC,2))),
         label_AUC=paste0("AUC = ", paste(round(AUC,2)))) -> data.labels

# plot 
ggroc(roc.list, legacy.axes = TRUE) +
  ggtitle("Area Under the Curve") +
  scale_color_discrete(labels=data.labels$label_long) +
  labs(x = "1 - Specificity / False positive rate", y = "Sensitivity / True positive rate",
       colour = "Model") +
  geom_abline() +
  theme_light() +
  theme(
    legend.position = c(0.9, 0.1),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6), 
    legend.text = element_text(size = 9, colour = "black", face="bold"), 
    legend.title = element_text(colour="black", size=10, 
                                face="bold"), 
    legend.background = element_rect(fill="lightblue",
                                     size=0.5, linetype="solid", 
                                     colour ="darkblue"))
ggsave("plot_AUC.pdf", width = 20, height = 20, units = "cm")

## Density plot for frequency of observation ----

# Prepare data frame for density plot
graph_data <- cbind(glm.probs$Pos, ranger.probs$Pos, rf.probs$Pos, kernlab.probs$Pos, e1071.probs$Pos)
colnames(graph_data) <- c('glm','ranger', 'randomForest', 'kernlab', 'e1071')
graph_data <- as.data.frame(graph_data)
graph_data <- mutate(graph_data, subject = row_number())

graph_data_long <- gather(graph_data, "algorithm", "probability", glm:e1071, factor_key = TRUE) #wide to long

# plot
ggplot(graph_data_long, aes(probability, color= algorithm)) +
  geom_density(alpha=0.3,
               kernel = "rectangular") + #smoothing parameter
  ggtitle("Comparison of distribution of risk probabilities") +
  theme_light() +
  theme(
    legend.position = c(0.9, 0.7),
    legend.background = element_rect(fill="lightblue",
                                     size=0.5, linetype="solid", 
                                     colour ="darkblue"),
    legend.title = element_text(colour="black", size=10, 
                                face="bold")) 
ggsave("plot_density.pdf", width = 20, height = 20, units = "cm")

## Individual subject comparison ---- 

# Sequences of 10 number to obtain random subjects
floor(runif(10, min=1, max=10002)) # 4549  566 7086 4362 1647 4855   24  808 8883 1147

graph_data_2 = graph_data

# Prepare data
graph_data_small <-  graph_data_2[graph_data_2$subject  %in% c(4549, 566, 4362, 7086, 1647, 4855, 24, 808, 8883, 1174), ]

graph_data_small_long <- gather(graph_data_small, "algorithm", "probability", glm:e1071, factor_key = TRUE) #wide to long

label_names <- c(
  `24` = "ID 24",
  `566` = "ID 566",
  `808` = "ID 808",
  `1174` = "ID 1174",
  `1647` = "ID 1647",
  `4362` = "ID 4362",
  `4549` = "ID 4549",
  `4855` = "ID 4855",
  `7086` = "ID 7086",
  `8883` = "ID 8883"
)

# plot
ggplot(graph_data_small_long, aes(x=subject, y=probability)) + 
  geom_point(aes(colour = algorithm), size = 3) +
  facet_grid(~subject, labeller = as_labeller(label_names), scale = "free", switch = "both") + 
  scale_y_continuous(breaks=seq(0,0.4,0.05)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11, color = "black", face = "bold"),
    legend.position = "bottom", 
    legend.title = element_blank(),
    legend.background = element_rect(fill="lightblue",
                                     size=0.5, linetype="solid", 
                                     colour ="darkblue"),
    strip.text.x = element_text(
      size = 8, color = "black", face = "bold.italic")
  ) 
ggsave("plot_subject_comparison.pdf", width = 20, height = 20, units = "cm") 

## Matrix scatter plot ----

# prepare data
graph_plot <- graph_data[ ,1:5]

# Using your function
custom_range <- function(data, mapping, ...) { 
  ggplot(data = data, mapping = mapping, ...) + 
    geom_point(size=0.5) + # to resize dots
    scale_x_continuous(limits = c(0, 1)) + # to rescale lower axes
    scale_y_continuous(limits = c(0, 1)) 
}

ggpairs(
  graph_plot,
  upper = list(continuous = GGally::wrap(ggally_cor, stars = F)), # to remove stars
  lower = list(continuous = custom_range)) # to resize points and axis scale

ggsave("plot_scatter.pdf", width = 20, height = 20, units = "cm")