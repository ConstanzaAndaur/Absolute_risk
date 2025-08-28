## Density plot for frequency of observation ----

# Prepare data frame for density plot
graph_data <- cbind(glm.probs$Pos, ridge.probs$Pos, rf.probs$Pos, kernlab.probs$Pos, nn_probs$Pos)
colnames(graph_data) <- c('UR', 'RR', 'RF', 'SVM', 'NN')
graph_data <- as.data.frame(graph_data)
graph_data <- mutate(graph_data, subject = row_number())

graph_data_long <- gather(graph_data, "algorithm", "probability", subject, factor_key = TRUE) #wide to long

library(tidyr)

graph_data_long <- graph_data %>%
  pivot_longer(
    cols = UR:NN,
    names_to = "algorithm",
    values_to = "probability"
  )

print(graph_data_long)

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





#
#################################################################################
## Individual subject comparison ---- 

# Generate 10 random numbers
set.seed(123) # For reproducibility
random_numbers <- sample(1:3915, 10)

# Print the random numbers
print(random_numbers)
#  [1] 2463 2511 2227  526  195 2986 1842 1142 3371 1253

graph_data_2 = graph_data

# Prepare data
graph_data_small <-  graph_data_2[graph_data_2$subject  %in% c(2463, 2511, 2227, 526, 195, 2986, 1842, 1142, 3371, 1253), ]

graph_data_small_long <- gather(graph_data_small, "algorithm", "probability", UR:NN, factor_key = TRUE) #wide to long

graph_data_small_long <- graph_data_small %>%
  pivot_longer(
    cols = UR:NN,
    names_to = "algorithm",
    values_to = "probability"
  )



label_names <- c(
  `195` = "ID 195",
  `526` = "ID 526",
  `1142` = "ID 1142",
  `1253` = "ID 1253",
  `1842` = "ID 1842",
  `2227` = "ID 2227",
  `2463` = "ID 2463",
  `2511` = "ID 2511",
  `2986` = "ID 2986",
  `3371` = "ID 3371"
)


# plot
ggplot(graph_data_small_long, aes(x=subject, y=probability)) + 
  geom_point(aes(colour = algorithm), size = 3) +
  facet_grid(~subject, labeller = as_labeller(label_names), scale = "free", switch = "both") + 
  scale_y_continuous(breaks=seq(0,1.0,0.05)) +
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




###################################################################################
#3
## Individual subject comparison ---- 

# Generate 10 random numbers
set.seed(123) # For reproducibility
random_numbers <- sample(1:3751, 10)

# Print the random numbers
print(random_numbers)
# [1] 2463 2511 2227  526 195 2986 1842 1142 3371 1253

graph_data_2 = graph_data

# Prepare data
graph_data_small <-  graph_data_2[graph_data_2$subject  %in% c(2463, 2511, 2227,  526, 195, 2986, 1842, 1142, 3371, 1253), ]

graph_data_small_long <- gather(graph_data_small, "algorithm", "probability", UR:NN, factor_key = TRUE) #wide to long

graph_data_small_long <- graph_data_small %>%
  pivot_longer(
    cols = UR:NN,
    names_to = "algorithm",
    values_to = "probability"
  )


label_names <- c(
  `195` = "ID 195",
  `526` = "ID 526",
  `1142` = "ID 1142",
  `1842` = "ID 1842",
  `2227` = "ID 2227",
  `2463` = "ID 2463",
  `2511` = "ID 2511",
  `2986` = "ID 2986",
  `3371` = "ID 3371",
  `1253` = "ID 1253",
  `5349` = "ID 5349"
)


# plot
ggplot(graph_data_small_long, aes(x=subject, y=probability)) + 
  geom_point(aes(colour = algorithm), size = 3) +
  facet_grid(~subject, labeller = as_labeller(label_names), scale = "free", switch = "both") + 
  scale_y_continuous(breaks=seq(0,0.6,0.05)) +
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


graph_data_small
#                UR         RR    RF        SVM         NN subject
# 195  0.011224804 0.02812270 0.000 0.10004679 0.01273403     195
# 526  0.141800737 0.15285210 0.002 0.11963595 0.13459500     526
# 1142 0.364056526 0.33234988 0.031 0.16455275 0.37311930    1142
# 1253 0.406857971 0.34586998 0.320 0.21506884 0.42074165    1253
# 1842 0.336416401 0.33601967 0.011 0.17397294 0.37141305    1842
# 2227 0.042039376 0.07647168 0.000 0.16443123 0.03876860    2227
# 2463 0.009795901 0.02689468 0.000 0.07959240 0.01154400    2463
# 2511 0.024742080 0.04922146 0.000 0.12078759 0.02554907    2511
# 2986 0.121851214 0.14546040 0.000 0.08992021 0.10784579    2986
# 3371 0.117180824 0.14362922 0.244 0.08350611 0.10183611    3371


library(dplyr)

train_data %>%
  filter(.id %in% c(2463, 2511, 2227,  526, 195, 2986, 1842, 1142, 3371, 1253)) %>%
  dplyr::select(.id, dvt)


#     .id dvt
#  1   195 Neg
#  2   526 Neg
#  3  1142 Pos
#  4  1253 Pos
#  5  1842 Neg
#  6  2227 Neg
#  7  2463 Neg
#  8  2511 Neg
#  9  2986 Neg
#  10 3371 Pos
