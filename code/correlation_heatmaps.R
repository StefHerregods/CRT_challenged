rm(list=ls())

library(lattice)
library(Rcpp)  
library(DEoptim)
library(gridExtra)

sourceCpp("code/ConfidenceBounds_scale4.cpp") 

text_settings <- list(
  par.settings = list(
    fontsize = list(text = 12, axis.text = 10),  # Adjusts font sizes
    fontfamily = "Helvetica"                     # Sets font family to Helvetica
  )
)
# binary confidence options variant --------------------------------------------


# Define parameter values
v <- 1
a <- 1
a_slope <- 0
ter <- 0
a2 <- seq(from = 1, to = 10, by = 0.5)
postdriftmod <- c(1, 3, 7)
ter2 <- 0
z <- 0.5  
ntrials <- 50000  
s <- 1  
dt <- 0.001 
a2_slope <- seq(from = 1, to = 10, by = 0.5)
starting_point_confidence <- 0.5

# Make predictions
df_correlations_spc1 <- NULL
df_correlations_spc2 <- NULL
df_correlations_spc3 <- NULL
for (i_a2 in a2){
  for (i_a2_slope in a2_slope){
    predictions_spc1 <- data.frame(DDM_confidence_bounds(v = v, a = a, a_slope = a_slope, ter = ter, a2 = i_a2, postdriftmod = postdriftmod[1], ter2 = ter2, a2_slope_upper = i_a2_slope, a2_slope_lower = i_a2_slope, starting_point_confidence = starting_point_confidence, z = z, ntrials = ntrials, s = s, dt = dt))
    names(predictions_spc1) <- c('rt', 'cor', 'rtconf', 'cj')
    df_correlations_spc1 <- rbind(df_correlations_spc1, c(postdriftmod[1], i_a2_slope, i_a2, cor(predictions_spc1$cj, predictions_spc1$rtconf)))
    predictions_spc2 <- data.frame(DDM_confidence_bounds(v = v, a = a, a_slope = a_slope, ter = ter, a2 = i_a2, postdriftmod = postdriftmod[2], ter2 = ter2, a2_slope_upper = i_a2_slope, a2_slope_lower = i_a2_slope, starting_point_confidence = starting_point_confidence, z = z, ntrials = ntrials, s = s, dt = dt))
    names(predictions_spc2) <- c('rt', 'cor', 'rtconf', 'cj')
    df_correlations_spc2 <- rbind(df_correlations_spc2, c(postdriftmod[2], i_a2_slope, i_a2, cor(predictions_spc2$cj, predictions_spc2$rtconf)))
    predictions_spc3 <- data.frame(DDM_confidence_bounds(v = v, a = a, a_slope = a_slope, ter = ter, a2 = i_a2, postdriftmod = postdriftmod[3], ter2 = ter2, a2_slope_upper = i_a2_slope, a2_slope_lower = i_a2_slope, starting_point_confidence = starting_point_confidence, z = z, ntrials = ntrials, s = s, dt = dt))
    names(predictions_spc3) <- c('rt', 'cor', 'rtconf', 'cj')
    df_correlations_spc3 <- rbind(df_correlations_spc3, c(postdriftmod[3], i_a2_slope, i_a2, cor(predictions_spc3$cj, predictions_spc3$rtconf)))
  }
  print(i_a2)
}  
df_correlations_spc1 <- data.frame(df_correlations_spc1)
df_correlations_spc2 <- data.frame(df_correlations_spc2)
df_correlations_spc3 <- data.frame(df_correlations_spc3)
names(df_correlations_spc1) <- c('v_ratio', 'a2_slope', 'a2', 'correlation')
names(df_correlations_spc2) <- c('v_ratio', 'a2_slope', 'a2', 'correlation')
names(df_correlations_spc3) <- c('v_ratio', 'a2_slope', 'a2', 'correlation')

# Transforming dataframe into a matrix
matrix_1 <- matrix(df_correlations_spc1$correlation, nrow = length(a2_slope), ncol = length(a2))
rownames(matrix_1) <- a2_slope
colnames(matrix_1) <- a2
matrix_2 <- matrix(df_correlations_spc2$correlation, nrow = length(a2_slope), ncol = length(a2))
rownames(matrix_2) <- a2_slope
colnames(matrix_2) <- a2
matrix_3 <- matrix(df_correlations_spc3$correlation, nrow = length(a2_slope), ncol = length(a2))
rownames(matrix_3) <- a2_slope
colnames(matrix_3) <- a2



# Define the breakpoints for the levels
breaks <- c(seq(-0.5, -0.2, length.out = 20), 
            seq(-0.199, 0.199, length.out = 60), 
            seq(0.2, 0.5, length.out = 20))




library(extrafont)
#font_import()
loadfonts(device = "win")

# Create the levelplot
p1 <- levelplot(matrix_1, 
          colorkey = FALSE,
          xlab = "CB urgency", 
          ylab = "CB separation", 
          at = breaks, 
          col.regions = colorRampPalette(c("#2A6A98", "#3B97DA", "#9ECAED", "#FFFFFF", "#FFCDA8", "#FE9C51", "#B16E39")),
          #main = "v_ratio = 1",
          scales = list(x = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10)), 
                        y = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10)),
                        tck = c(1,0)),
          par.settings = text_settings$par.settings)
p2 <- levelplot(matrix_2, 
          colorkey = FALSE,
          xlab = "CB urgency", 
          ylab = "CB separation", 
          at = breaks, 
          col.regions = colorRampPalette(c("#2A6A98", "#3B97DA", "#9ECAED", "#FFFFFF", "#FFCDA8", "#FE9C51", "#B16E39")),
          #main = "v_ratio = 3",
          scales = list(x = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10)), 
                        y = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10)),
                        tck = c(1,0)),
          par.settings = text_settings$par.settings)
p3 <- levelplot(matrix_3, 
          colorkey = list(space = "bottom", width = 0.7, cex.labels = 0.7),
          #colorkey = FALSE,
          xlab = "CB urgency", 
          ylab = "CB separation", 
          at = breaks, 
          col.regions = colorRampPalette(c("#2A6A98", "#3B97DA", "#9ECAED", "#FFFFFF", "#FFCDA8", "#FE9C51", "#B16E39")),
          #main = "v_ratio = 7",
          scales = list(x = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10), alternating = c(1,0)), 
                        y = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10)),
                        tck = c(1,0)),
          par.settings = text_settings$par.settings)

# Save plots
png("figures/heatmap_v1_1.png", width = 2.4, height = 2.4, units = "in", res = 1000)
print(p1)
dev.off()
png("figures/heatmap_v1_2.png", width = 2.4, height = 2.4, units = "in", res = 1000)
print(p2)
dev.off()
png("figures/heatmap_v1_3.png", width = 2.5, height = 2.85, units = "in", res = 1000)
print(p3)
dev.off()


# n-choice model variant -------------------------------------------------------


# Define parameter values
v <- 1
a <- 1
a_slope <- 0
ter <- 0
a2 <- 5
postdriftmod <- 3
ter2 <- 0
z <- 0.5  
ntrials <- 50000  
s <- 1  
dt <- 0.001 
a2_slope_upper <- seq(from = 1, to = 10, by = 0.5)
a2_slope_lower <- seq(from = 1, to = 10, by = 0.5)
starting_point_confidence <- c(0.4, 0.5, 0.6)

# Make predictions
df_correlations_spc1 <- NULL
df_correlations_spc2 <- NULL
df_correlations_spc3 <- NULL
for (i_a2_slope_upper in a2_slope_upper){
  for (i_a2_slope_lower in a2_slope_lower){
    predictions_spc1 <- data.frame(DDM_confidence_bounds(v = v, a = a, a_slope = a_slope, ter = ter, a2 = a2, postdriftmod = postdriftmod, ter2 = ter2, a2_slope_upper = i_a2_slope_upper, a2_slope_lower = i_a2_slope_lower, starting_point_confidence = starting_point_confidence[1], z = z, ntrials = ntrials, s = s, dt = dt))
    names(predictions_spc1) <- c('rt', 'cor', 'rtconf', 'cj')
    df_correlations_spc1 <- rbind(df_correlations_spc1, c(starting_point_confidence[1], i_a2_slope_upper, i_a2_slope_lower, cor(predictions_spc1$cj, predictions_spc1$rtconf)))
    predictions_spc2 <- data.frame(DDM_confidence_bounds(v = v, a = a, a_slope = a_slope, ter = ter, a2 = a2, postdriftmod = postdriftmod, ter2 = ter2, a2_slope_upper = i_a2_slope_upper, a2_slope_lower = i_a2_slope_lower, starting_point_confidence = starting_point_confidence[2], z = z, ntrials = ntrials, s = s, dt = dt))
    names(predictions_spc2) <- c('rt', 'cor', 'rtconf', 'cj')
    df_correlations_spc2 <- rbind(df_correlations_spc2, c(starting_point_confidence[2], i_a2_slope_upper, i_a2_slope_lower, cor(predictions_spc2$cj, predictions_spc2$rtconf)))
    predictions_spc3 <- data.frame(DDM_confidence_bounds(v = v, a = a, a_slope = a_slope, ter = ter, a2 = a2, postdriftmod = postdriftmod, ter2 = ter2, a2_slope_upper = i_a2_slope_upper, a2_slope_lower = i_a2_slope_lower, starting_point_confidence = starting_point_confidence[3], z = z, ntrials = ntrials, s = s, dt = dt))
    names(predictions_spc3) <- c('rt', 'cor', 'rtconf', 'cj')
    df_correlations_spc3 <- rbind(df_correlations_spc3, c(starting_point_confidence[3], i_a2_slope_upper, i_a2_slope_lower, cor(predictions_spc3$cj, predictions_spc3$rtconf)))
  }
  print(i_a2_slope_upper)
}  
df_correlations_spc1 <- data.frame(df_correlations_spc1)
df_correlations_spc2 <- data.frame(df_correlations_spc2)
df_correlations_spc3 <- data.frame(df_correlations_spc3)
names(df_correlations_spc1) <- c('starting_point_confidence', 'a2_slope_upper', 'a2_slope_lower', 'correlation')
names(df_correlations_spc2) <- c('starting_point_confidence', 'a2_slope_upper', 'a2_slope_lower', 'correlation')
names(df_correlations_spc3) <- c('starting_point_confidence', 'a2_slope_upper', 'a2_slope_lower', 'correlation')

# Transforming dataframe into a matrix
matrix_1 <- matrix(df_correlations_spc1$correlation, nrow = length(a2_slope_upper), ncol = length(a2_slope_lower))
rownames(matrix_1) <- a2_slope_upper
colnames(matrix_1) <- a2_slope_lower
matrix_2 <- matrix(df_correlations_spc2$correlation, nrow = length(a2_slope_upper), ncol = length(a2_slope_lower))
rownames(matrix_2) <- a2_slope_upper
colnames(matrix_2) <- a2_slope_lower
matrix_3 <- matrix(df_correlations_spc3$correlation, nrow = length(a2_slope_upper), ncol = length(a2_slope_lower))
rownames(matrix_3) <- a2_slope_upper
colnames(matrix_3) <- a2_slope_lower

# Define the breakpoints for the levels
breaks <- c(seq(-1, -0.5, length.out = 20), 
            seq(-0.499, 0.499, length.out = 60), 
            seq(0.5, 1, length.out = 20))
breaks <- c(seq(-0.5, -0.2, length.out = 20), 
            seq(-0.199, 0.199, length.out = 60), 
            seq(0.2, 0.5, length.out = 20))
breaks <- c(seq(-1, 1, length.out = 100))

# Creating the heatmap
p1 <- levelplot(matrix_1, 
          colorkey = FALSE,
          xlab = "Lower CB urgency", 
          ylab = "Upper CB urgency", 
          #main = "confidence starting point bias = 0.4",
          at = breaks, 
          col.regions = colorRampPalette(c("#2A6A98", "#3B97DA", "#9ECAED", "#FFFFFF", "#FFCDA8", "#FE9C51", "#B16E39")),
          scales = list(x = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10), alternating = c(1,0)), 
                        y = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10)),
                        tck = c(1,0)),
          par.settings = text_settings$par.settings)
p2 <- levelplot(matrix_2, 
          colorkey = FALSE,
          xlab = "Lower CB urgency", 
          ylab = "Upper CB urgency", 
          #main = "confidence starting point bias = 0.5",
          at = breaks, 
          col.regions = colorRampPalette(c("#2A6A98", "#3B97DA", "#9ECAED", "#FFFFFF", "#FFCDA8", "#FE9C51", "#B16E39")),
          scales = list(x = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10), alternating = c(1,0)), 
                        y = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10)),
                        tck = c(1,0)),
          par.settings = text_settings$par.settings)
p3 <- levelplot(matrix_3, 
          #colorkey = list(space = "bottom", width = 0.7, cex.labels = 0.7),
          colorkey = FALSE,
          xlab = "Lower CB urgency", 
          ylab = "Upper CB urgency", 
          #main = "confidence starting point bias = 0.6",
          at = breaks, 
          col.regions = colorRampPalette(c("#2A6A98", "#3B97DA", "#9ECAED", "#FFFFFF", "#FFCDA8", "#FE9C51", "#B16E39")),
          scales = list(x = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10), alternating = c(1,0)), 
                        y = list(at = c(3,7,11,15,19), labels = c(2,4,6,8,10)),
                        tck = c(1,0)),
          par.settings = text_settings$par.settings)

# Save plots
png("figures/heatmap_v2_1.png", width = 2.4, height = 2.4, units = "in", res = 1000)
print(p1)
dev.off()
png("figures/heatmap_v2_2.png", width = 2.4, height = 2.4, units = "in", res = 1000)
print(p2)
dev.off()
png("figures/heatmap_v2_3.png", width = 2.4, height = 2.4, units = "in", res = 1000) 
print(p3)
dev.off()
















