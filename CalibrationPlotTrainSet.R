library(CalibrationCurves)

# Calibration Plot for the Train Set

# Obtain predicted probabilities from the training set

train_probs_glm <- predict(glm_mod, train_data, type = "prob")
train_probs_rg <- predict(ridge_mod, train_data, type = "prob")
train_probs_rf <- predict(rf_mod, train_data, type = "prob")
train_probs_svm <- predict(svm_mod, train_data, type = "prob")
train_probs_nn <- predict(nn_mod, train_data, type = "prob")


# Extract the true labels from the training set
yTrain <- ifelse(train_data$dvt == "Neg", 0, 1)

# Generate Calibration Plots for the Training Set

pred_probs2 <- train_probs_glm[, "Pos"]
valProbggplot(pred_probs2, yTrain)

valProbggplot(pred_probs2, yTrain, statloc = FALSE)

#valProbggplot(pred_probs2, yTrain, dostats = c ("Intercept", "Slope"))


pred_probs_rf3 <- train_probs_rf[, "Pos"]
pred_probs_SVM2 <- train_probs_svm[, "Pos"]
pred_probs_rg2 <- train_probs_rg[, "Pos"]
pred_probs_nn <- train_probs_nn[, "Pos"]

#valProbggplot(pred_probs_rf3, yTrain)
#valProbggplot(pred_probs_SVM2, yTrain)
#valProbggplot(pred_probs_rg2, yTrain)
#valProbggplot(pred_probs_nn, yTrain)

valProbggplot(pred_probs_SVM2, yTrain, statloc = FALSE)
valProbggplot(pred_probs_rg2, yTrain, statloc = FALSE)
valProbggplot(pred_probs_nn, yTrain, statloc = FALSE)




# RF
pred_probs_rf4 <- pmin(pmax(pred_probs_rf3, 1e-5), 1 - 1e-5)
#valProbggplot(pred_probs_rf4, yTrain)
valProbggplot(pred_probs_rf4, yTrain, statloc = FALSE)

 

#
# Re-plot

library(gridGraphics)
#library(gridExtra)
library(grid)


valProbggplot(pred_probs2, yTrain, statloc = FALSE)
grid1 <- grid.grab()

valProbggplot(pred_probs_rg2, yTrain, statloc = FALSE)
grid2 <- grid.grab()

valProbggplot(pred_probs_rf4, yTrain, statloc = FALSE)
grid3 <- grid.grab()

valProbggplot(pred_probs_SVM2, yTrain, statloc = FALSE)
grid4 <- grid.grab()

valProbggplot(pred_probs_nn, yTrain, statloc = FALSE)
grid5 <- grid.grab()


# Create labeled grobs with the model name above the plot
g1 <- arrangeGrob(grid1, top = textGrob("UR", gp = gpar(fontsize = 14, fontface = "bold")))
g2 <- arrangeGrob(grid2, top = textGrob("RR", gp = gpar(fontsize = 14, fontface = "bold")))
g3 <- arrangeGrob(grid3, top = textGrob("RF", gp = gpar(fontsize = 14, fontface = "bold")))
g4 <- arrangeGrob(grid4, top = textGrob("SVM", gp = gpar(fontsize = 14, fontface = "bold")))
g5 <- arrangeGrob(grid5, top = textGrob("NN", gp = gpar(fontsize = 14, fontface = "bold")))

# Arrange all labeled plots in a grid
grid.arrange(g1, g2, g3, g4, g5, nrow = 2, ncol = 3)

