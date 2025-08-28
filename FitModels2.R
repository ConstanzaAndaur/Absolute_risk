# Cross-validation ----

n_cv_folds <- 10
n_cv_repetitions <- 1

fitControl <- trainControl(method = "repeatedcv", # Cross-validation, default is bootstrap
                           number = n_cv_folds,
                           repeats = n_cv_repetitions,
                           classProbs = TRUE,
                           savePredictions = TRUE,
                           index = createMultiFolds(train_data$dvt,n_cv_folds, n_cv_repetitions)
)


##################################################################################################
##################################################################################################
#1 Fit logistic regression ----
set.seed(825)
glm_mod <- train(dvt ~ age + sex + ddimdich  + hist + malign + altdiagn,
                 data = train_data,
                 method = "glm",
                 family = binomial (link = 'logit'),
                 trControl = fitControl,
                 metric = "ROC")

getTrainPerf(glm_mod) 
glm_mod$results # AUC, Sens, Spec with SD
glm_mod$finalModel

summary(glm_mod)
# Coefficients:
#                 Estimate  Std. Error  z value  Pr(>|z|)    
# (Intercept).  -2.840151   0.168170  -16.889    < 2e-16 ***
#  age          -0.004932   0.002283  -2.160     0.0307 *  
#  sex1          0.439916   0.076537   5.748    9.04e-09 ***
#  ddimdich1     2.606409   0.112015   23.268    < 2e-16 ***
#  hist1         0.749407   0.155067   4.833.    1.35e-06 ***
#  malign1       0.888034   0.107660   8.249     < 2e-16 ***
#  altdiagn1    -1.371382   0.086180   -15.913   < 2e-16 ***
  
## Predictions with test data ----

glm.probs <- predict(glm_mod,
                     newdata = train_data,
                     type = "prob",  
                     se.fit = TRUE) 
head(glm.probs)
hist(glm.probs$Pos)

# Set threshold to 0.02 and classify
threshold <- 0.02

glm.pred <- ifelse(glm.probs[, "Pos"] > threshold, "Pos", "Neg")
glm.pred <- factor(glm.pred, levels = c("Neg", "Pos"))


## Confusion matrix ----
confusionMatrix(data = glm.pred, 
                reference = train_data$dvt,
                positive = "Pos")



#  Sensitivity : 0.9843 Specificity : 0.3180 PPV : 0.2508  NPV : 0.9887 
##################################################################################################


##################################################################################################
#2 Fit Random Forest
tune_grid <- expand.grid(
  .mtry = c(2, 3, 4, 5, 6) # Adjust based on the number of predictors
)

set.seed(825)
rf_mod <- train(
  dvt ~ age + sex + ddimdich + hist + malign + altdiagn, 
  data = train_data,
  method = "rf",
  trControl = fitControl,
  metric = "ROC",
  tuneGrid = tune_grid,
  ntree = 1000
) 

rf_mod$results
rf_mod$finalModel

## Predict on test data ----

rf.probs <- predict(rf_mod,
                    train_data,
                    type = "prob")

hist(rf.probs$Pos, breaks = 100, ylim = c(0, 2500))
hist(rf_probs2$Pos, breaks = 100, ylim = c(0, 2500))
hist(rf.probs_new$Pos, breaks = 100, ylim = c(0, 2500))

# Set threshold to 0.02

threshold <- 0.02
rf.pred <- ifelse(rf.probs[, "Pos"] > threshold, "Pos", "Neg")
rf.pred <- factor(rf.pred, levels = c("Neg", "Pos"))

confusionMatrix(data = rf.pred, 
                reference = train_data$dvt, 
                mode = "everything", 
                positive = "Pos")

# Sensitivity : 0.8098, Specificity : 0.8106, PPV : 0.4979, NPV : 0.9484

# rf
# Extract variable importance
rf_importance <- varImp(rf_mod)

# Print the importance
print(rf_importance)

# rf variable importance
#           Overall
#ddimdich1  100.000
#altdiagn1  37.114
#age        34.024
#malign1     7.461
#sex1        1.564
#hist1       0.000

# Plot importance
plot(rf_importance)

##################################################################################################

##################################################################################################
#3 Fit SVM 
set.seed(825)

tune_grid = data.frame(C = 1/2^(3:7))

svm_mod <- train(dvt ~ age + sex + ddimdich + hist + malign + altdiagn,
                 data = train_data,
                 method = "svmLinear", #kernlab
                 trControl = fitControl,
                 tuneGrid = tune_grid)

svm_mod$results
svm_mod



## Prediction on test data ----

kernlab.probs <- predict(svm_mod, 
                         train_data,
                         type = "prob") 


hist(kernlab.probs$Pos)

# Set threshold to 0.02
threshold <- 0.02
svm.pred <- ifelse(kernlab.probs[, "Pos"] > threshold, "Pos", "Neg")
svm.pred <- factor(svm.pred, levels = c("Neg", "Pos"))

confusionMatrix(data = svm.pred, 
                reference = train_data$dvt, 
                mode = "everything", 
                positive = "Pos")


# Sensitivity : 1.000, Specificity : 0.000, PPV : 0.1883, NPV : NaN 

# Extract the coefficients of the SVM model
SVM_final_model <- svm_mod$finalModel

# Extract the weight vector (coefficients)
coefs <- t(SVM_final_model@coef[[1]]) %*% SVM_final_model@xmatrix[[1]]

# Convert to named vector (match to predictors)
coef_named <- setNames(as.vector(coefs), colnames(train_data[, c("age", "sex", "ddimdich", "hist", "malign", "altdiagn")]))

# Print coefficients
coef_named
#           age           sex      ddimdich          hist        malign      altdiagn 
# -4.779067e-05  8.897928e-05  1.084276e-05  7.924919e-06  2.049564e-05 -9.351424e-05 

abs(coef_named)
#          age          sex     ddimdich         hist       malign     altdiagn 
# 4.779067e-05 8.897928e-05 1.084276e-05 7.924919e-06 2.049564e-05 9.351424e-05 




##################################################################################################

##################################################################################################
#4 # Fit Ridge Model ----
set.seed(825)
ridge_mod <- train(dvt ~ age + sex + ddimdich + hist + malign + altdiagn,
                   data = train_data,
                   method = "glmnet",
                   trControl = fitControl,
                   tuneGrid = expand.grid(alpha = 0, # Alpha = 0 for Ridge
                                          lambda = seq(0.001, 0.1, by = 0.001)), # Adjust lambda range if needed
                   metric = "ROC") # Optimize model based on ROC

ridge_mod$results
ridge_mod$finalModel
ridge_mod$bestTune


## Predict on test data ----

ridge.probs <- predict(ridge_mod,
                       train_data,
                       type = "prob")

# Set threshold to 0.02
threshold <- 0.02
ridge.pred <- ifelse(ridge.probs[, "Pos"] > threshold, "Pos", "Neg")
ridge.pred <- factor(ridge.pred, levels = c("Neg", "Pos"))

confusionMatrix(data = ridge.pred, 
                reference = train_data$dvt, 
                mode = "everything", 
                positive = "Pos")

hist(ridge.probs$Pos)
# Sensitivity : 1.0000, Specificity : 0.0000, PPV : 0.1883, NPV : na 

# RR

# Get the final model object from the train object
ridge_model <- ridge_mod$finalModel

# Extract coefficients (for the best lambda value)
ridge_coefs <- coef(ridge_model, s = ridge_model$lambdaOpt)

# Print the coefficients
print(ridge_coefs)
# 7 x 1 sparse Matrix of class "dgCMatrix"
                  s1
#  (Intercept) -2.39584310 
# age         -0.00164548
# sex1         0.34630624
# ddimdich1    1.80638644
# hist1        0.56343851
# malign1      0.75402642
# altdiagn1   -1.05779130
                  
##################################################################################################
# neural network

library(caret)
library(nnet)

set.seed(825)
nn_mod <- train(dvt ~ age + sex + ddimdich + hist + malign + altdiagn,
                data = train_data,
                method = "nnet",
                trControl = fitControl,
                metric = "ROC",
                trace = FALSE)  # Suppresses output

nn_probs <- predict(nn_mod, newdata = train_data, type = "prob")  
head(nn_probs)  # Check probability predictions
nn_mod$results      # Shows model performance for each combination
nn_mod$bestTune     # Shows the best combination selected
#   size decay
#3    1   0.1


threshold <- 0.2

nn_pred <- ifelse(nn_probs[, "Pos"] > threshold, "Pos", "Neg")
nn_pred <- factor(nn_pred, levels = c("Neg", "Pos"))  # Ensure factor levels

confusionMatrix(data = nn_pred, 
                reference = train_data$dvt,
                positive = "Pos")


# Sensitivity : 0.8159, Specificity : 0.7472, PPV : 0.4281, NPV : 0.9459 


# nn
# Extract the weights from the trained neural network model
nn_weights <- coef(nn_mod$finalModel)

# Print the weights
print(nn_weights)

#     b->h1       i1->h1       i2->h1       i3->h1       i4->h1       i5->h1       i6->h1         b->o        h1->o 
# 0.483738873  0.004272649 -0.333799210 -1.751333937 -0.612318372 -0.717689994  0.989380509  1.292043744 -6.669148868 



# Step 1
nn_weights <- coef(nn_mod$finalModel)

# Step 2
predictor_names <- c("age", "sex", "ddimdich", "hist", "malign", "altdiagn")
nn_weight_names <- names(nn_weights)

# Step 3: Replace "i1", "i2", etc. with actual predictor names
for (i in seq_along(predictor_names)) {
  nn_weight_names <- gsub(paste0("i", i, "->"), paste0(predictor_names[i], "->"), nn_weight_names)
}

# Step 4: Assign the new names back
names(nn_weights) <- nn_weight_names

# Step 5: Convert to data frame for easier viewing
nn_weights_df <- data.frame(Connection = names(nn_weights),
                            Weight = nn_weights,
                            row.names = NULL)

# 
print(nn_weights_df)

#     Connection       Weight
#1        b->h1  0.483738873
#2      age->h1  0.004272649
#3      sex->h1 -0.333799210
#4 ddimdich->h1 -1.751333937
#5     hist->h1 -0.612318372
#6   malign->h1 -0.717689994
#7 altdiagn->h1  0.989380509
#8         b->o  1.292043744
#9        h1->o -6.669148868

# b->h1, b->o, and h1->o — are the bias and connection weights inside the neural network

# b->h1  Hidden layer intercept. Bias weight going to the hidden neuron (h1) — a constant added to its input before activation
# b->o   Output layer intercept. Weight from the hidden neuron to the output neuron
# h1->o  Influence of hidden neuron on output. Bias weight going to the output neuron — a constant added before output activation

##################################################################################################


##################################################################################################
# Evaluation ----

## AUC plot for all four models (train set) ----

rocobj1 <- plot.roc(as.numeric(glm_mod$trainingData$.outcome=='Pos'),
                    aggregate(Pos~rowIndex,glm_mod$pred,mean)[,'Pos'], 
                    print.auc=TRUE)

rocobj2 <- lines.roc(as.numeric(rf_mod$trainingData$.outcome=='Pos'),
                     aggregate(Pos~rowIndex,rf_mod$pred,mean)[,'Pos'])

rocobj3 <- lines.roc(as.numeric(svm_mod$trainingData$.outcome=='Pos'),
                     aggregate(Pos~rowIndex,svm_mod$pred,mean)[,'Pos'])

rocobj4 <- lines.roc(as.numeric(ridge_mod$trainingData$.outcome=='Pos'),
                     aggregate(Pos~rowIndex,ridge_mod$pred,mean)[,'Pos'])

rocobj5 <- lines.roc(as.numeric(nn_mod$trainingData$.outcome=='Pos'),
                     aggregate(Pos~rowIndex,nn_mod$pred,mean)[,'Pos'])


roc.list <- list(UR = rocobj1, 
                 RR = rocobj4,
                 RF = rocobj2, 
                 SVM = rocobj3, 
                 NN = rocobj5
)

ci.auc(rocobj1) # Logistic 95% CI: 0.8302-0.853 (DeLong)
ci.auc(rocobj4) # Ridge 95% CI: 0.8283-0.8513 (DeLong)
ci.auc(rocobj2) # RandomForest 95% CI: 0.7918-0.8187 (DeLong)
ci.auc(rocobj3) # 95% CI: 0.8061-0.8311 (DeLong)
ci.auc(rocobj5) # NeuralNetwork 95% CI: 0.8294-0.8525 (DeLong)

# extract AUC
data.auc <- roc.list %>% 
  map(~ tibble(AUC = as.numeric(.x$auc))) %>%  # Ensure AUC is numeric
  bind_rows(.id = "name")

data.auc
# name            AUC
#1 Logistic      0.842
#2 Ridge         0.840
#3 RandomForest  0.805
#4 SVM.kernlab   0.819
#5 NeuralNetwork 0.841


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


########################################################################