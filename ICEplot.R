library(iml)
library(ggplot2)

# Create predictor object
X <- train_data[, c("age", "sex", "ddimdich", "hist", "malign", "altdiagn")]  # Select features
y <- train_data$dvt  # Outcome variable


predictor_glm <- Predictor$new(glm_mod, data = X, y = y, type = "prob")
predictor_ridge <- Predictor$new(ridge_mod, data = X, y = y, type = "prob")
predictor_rf <- Predictor$new(rf_mod, data = X, y = y, type = "prob")
predictor_svm <- Predictor$new(svm_mod, data = X, y = y, type = "prob")
predictor_nn <- Predictor$new(nn_mod, data = X, y = y, type = "prob")

# Compute ICE curves for the feature "malign"
ice_glm1 <- FeatureEffect$new(predictor_glm, feature = "malign", method = "ice")
ice_glm2 <- FeatureEffect$new(predictor_glm, feature = "ddimdich", method = "ice")
ice_glm3 <- FeatureEffect$new(predictor_glm, feature = "altdiagn", method = "ice")
ice_glm4 <- FeatureEffect$new(predictor_glm, feature = "age", method = "ice")
ice_glm5 <- FeatureEffect$new(predictor_glm, feature = "sex", method = "ice")
ice_glm6 <- FeatureEffect$new(predictor_glm, feature = "hist", method = "ice")

# Plot the ICE curves
p1 <- plot(ice_glm1) + ggtitle("ICE Plot for malign (Logistic Regression)")
p2 <- plot(ice_glm2) + ggtitle("ICE Plot for ddimdich (Logistic Regression)")
p3 <- plot(ice_glm3) + ggtitle("ICE Plot for altdiagn (Logistic Regression)")
p4 <- plot(ice_glm4) + ggtitle("ICE Plot for age (Logistic Regression)")
p5 <- plot(ice_glm5) + ggtitle("ICE Plot for his (Logistic Regression)")
p6 <- plot(ice_glm6) + ggtitle("ICE Plot for sex (Logistic Regression)")

# Random Forest
ice_rf1 <- FeatureEffect$new(predictor_rf, feature = "malign", method = "ice")
ice_rf2 <- FeatureEffect$new(predictor_rf, feature = "ddimdich", method = "ice")
ice_rf3 <- FeatureEffect$new(predictor_rf, feature = "altdiagn", method = "ice")
ice_rf4 <- FeatureEffect$new(predictor_rf, feature = "age", method = "ice")
ice_rf5 <- FeatureEffect$new(predictor_rf, feature = "sex", method = "ice")
ice_rf6 <- FeatureEffect$new(predictor_rf, feature = "hist", method = "ice")

p7 <- plot(ice_rf1) + ggtitle("ICE Plot for malign (Random Forest)")
p8 <- plot(ice_rf2) + ggtitle("ICE Plot for ddimdich (Random Forest)")
p9 <- plot(ice_rf3) + ggtitle("ICE Plot for altdiagn (Random Forest)")
p10 <- plot(ice_rf4) + ggtitle("ICE Plot for age (Random Forest)")
p11 <- plot(ice_rf5) + ggtitle("ICE Plot for his (Random Forest)")
p12 <- plot(ice_rf6) + ggtitle("ICE Plot for sex (Random Forest)")


# SVM
ice_svm1 <- FeatureEffect$new(predictor_svm, feature = "malign", method = "ice")
ice_svm2 <- FeatureEffect$new(predictor_svm, feature = "ddimdich", method = "ice")
ice_svm3 <- FeatureEffect$new(predictor_svm, feature = "altdiagn", method = "ice")
ice_svm4 <- FeatureEffect$new(predictor_svm, feature = "age", method = "ice")
ice_svm5 <- FeatureEffect$new(predictor_svm, feature = "sex", method = "ice")
ice_svm6 <- FeatureEffect$new(predictor_svm, feature = "hist", method = "ice")

p13 <- plot(ice_svm1) + ggtitle("ICE Plot for malign (SVM)")
p14 <- plot(ice_svm2) + ggtitle("ICE Plot for ddimdich (SVM)")
p15 <- plot(ice_svm3) + ggtitle("ICE Plot for altdiagn (SVM)")
p16 <- plot(ice_svm4) + ggtitle("ICE Plot for age (SVM)")
p17 <- plot(ice_svm5) + ggtitle("ICE Plot for his (SVM)")
p18 <- plot(ice_svm6) + ggtitle("ICE Plot for sex (SVM)")


# Ridge Regression
ice_ridge1 <- FeatureEffect$new(predictor_ridge, feature = "malign", method = "ice")
ice_ridge2 <- FeatureEffect$new(predictor_ridge, feature = "ddimdich", method = "ice")
ice_ridge3 <- FeatureEffect$new(predictor_ridge, feature = "altdiagn", method = "ice")
ice_ridge4 <- FeatureEffect$new(predictor_ridge, feature = "age", method = "ice")
ice_ridge5 <- FeatureEffect$new(predictor_ridge, feature = "sex", method = "ice")
ice_ridge6 <- FeatureEffect$new(predictor_ridge, feature = "hist", method = "ice")

p19 <- plot(ice_ridge1) + ggtitle("ICE Plot for malign (Ridge Regression)")
p20 <- plot(ice_ridge2) + ggtitle("ICE Plot for ddimdich (Ridge Regression)")
p21 <- plot(ice_ridge3) + ggtitle("ICE Plot for altdiagn (Ridge Regression)")
p22 <- plot(ice_ridge4) + ggtitle("ICE Plot for age (Ridge Regression)")
p23 <- plot(ice_ridge5) + ggtitle("ICE Plot for his (Ridge Regression)")
p24 <- plot(ice_ridge6) + ggtitle("ICE Plot for sex (Ridge Regression)")
 
# Neural Network
ice_nn1 <- FeatureEffect$new(predictor_nn, feature = "malign", method = "ice")
ice_nn2 <- FeatureEffect$new(predictor_nn, feature = "ddimdich", method = "ice")
ice_nn3 <- FeatureEffect$new(predictor_nn, feature = "altdiagn", method = "ice")
ice_nn4 <- FeatureEffect$new(predictor_nn, feature = "age", method = "ice")
ice_nn5 <- FeatureEffect$new(predictor_nn, feature = "sex", method = "ice")
ice_nn6 <- FeatureEffect$new(predictor_nn, feature = "hist", method = "ice")

p25 <- plot(ice_nn1) + ggtitle("ICE Plot for malign (Neural Network)")
p26 <- plot(ice_nn2) + ggtitle("ICE Plot for ddimdich (Neural Network)")
p27 <- plot(ice_nn3) + ggtitle("ICE Plot for altdiagn (Neural Network)")
p28 <- plot(ice_nn4) + ggtitle("ICE Plot for age (Neural Network)")
p29 <- plot(ice_nn5) + ggtitle("ICE Plot for his (Neural Network)")
p30 <- plot(ice_nn6) + ggtitle("ICE Plot for sex (Neural Network)")


library(gridExtra)

#malign
grid.arrange(p1, p7, p13, p19, p25, ncol = 2)

#ddimdich
grid.arrange(p2, p8, p14, p20, p26, ncol = 2)

#altdiagn
grid.arrange(p3, p9, p15, p21, p27, ncol = 2)

#age
grid.arrange(p4, p10, p16, p22, p28, ncol = 2)

#his
grid.arrange(p5, p11, p17, p23, p29, ncol = 2)

#sex
grid.arrange(p6, p12, p18, p24, p30, ncol = 2)








# nn
# Extract the weights from the trained neural network model
nn_weights <- coef(nn_mod$finalModel)

# Print the weights
print(nn_weights)

#     b->h1       i1->h1       i2->h1       i3->h1       i4->h1       i5->h1       i6->h1         b->o        h1->o 
# 0.483738873  0.004272649 -0.333799210 -1.751333937 -0.612318372 -0.717689994  0.989380509  1.292043744 -6.669148868 



# rf
# Extract variable importance
rf_importance <- varImp(rf_mod)

# Print the importance
print(rf_importance)

# Plot importance
plot(rf_importance)


# SVM
# Get the SVM model object from the caret train object
svm_model <- svm_mod$finalModel

# Extract coefficients for the linear SVM model
svm_coefs <- svm_model@coef

# Print the coefficients
print(svm_coefs)


# Extract coefficients (weights) from the model
coefficients <- svm_model@coef


intercept <- svm_model@rho

predictor_names <- colnames(svm_model@xmatrix[[1]])
coef_df <- data.frame(Predictor = predictor_names, Coefficient = coefficients)


# LR




# RR

# Get the final model object from the train object
ridge_model <- ridge_mod$finalModel

# Extract coefficients (for the best lambda value)
ridge_coefs <- coef(ridge_model, s = ridge_model$lambdaOpt)

# Print the coefficients
print(ridge_coefs)
