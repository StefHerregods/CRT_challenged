################################################################################
# Analysis of the estimated model parameters and their respective predictions  #
################################################################################

rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)


# Separate dataset analyses ----------------------------------------------------

# Choose a dataset to analyse
# [1] Bang et al. (2019) 
# [2]
dataset <- 1

# Load data
if (dataset == 1){
  df_obs <- read.csv('data/observations/data_Bang_2019_Exp2_preprocessed.csv')
  df_params <- read.csv('data/model_parameters/Bang2019/Bang2019_params.csv')
  df_pred <- read.csv('data/model_predictions/Bang2019_predictions.csv')
}


# Investigating the confidence RT vs. accuracy correlations
correlations_accuracy <- NULL
for (sub in unique(df_obs$sub)){
  correlations_accuracy <- rbind(correlations_accuracy, c(cor(df_obs[df_obs$sub==sub,]$rtconf, df_obs[df_obs$sub==sub,]$cor), cor(df_pred[df_pred$sub==sub,]$rtconf, df_pred[df_pred$sub==sub,]$cor)))
}
correlations_accuracy <- data.frame(correlations_accuracy)
names(correlations_accuracy) <- c('observed_correlations', 'predicted_correlations')
plot_accuracy <- ggplot(data = correlations_accuracy, aes(x = observed_correlations, y = predicted_correlations)) +
  geom_point(shape = 16, size = 3, alpha = 0.3, color = '#145369') +
  xlim(-1,1) +
  ylim(-1,1) +
  ggtitle('confidence RT ~ accuracy') +
  xlab('Observed correlation') +
  ylab('Simulated correlation') + 
  theme_classic() +
  geom_abline(slope = 1) +
  geom_abline(intercept = 0, slope = 0, color = 'darkred', size = 2, lty = 2)

# Investigating the confidence RT vs. confidence judgment correlations
correlations_cj <- NULL
for (sub in unique(df_obs$sub)){
  correlations_cj <- rbind(correlations_cj, c(cor(df_obs[df_obs$sub==sub,]$rtconf, df_obs[df_obs$sub==sub,]$cj), cor(df_pred[df_pred$sub==sub,]$rtconf, df_pred[df_pred$sub==sub,]$cj)))
}
correlations_cj <- data.frame(correlations_cj)
names(correlations_cj) <- c('observed_correlations', 'predicted_correlations')
plot_cj <- ggplot(data = correlations_cj, aes(x = observed_correlations, y = predicted_correlations)) +
  geom_point(shape = 16, size = 3, alpha = 0.3, color = '#145369') +
  xlim(-1,1) +
  ylim(-1,1) +
  ggtitle('confidence RT ~ CJ') +
  xlab('Observed correlation') +
  ylab('Simulated correlation') + 
  theme_classic() +
  geom_abline(slope = 1) +
  geom_abline(intercept = 0, slope = 0, color = 'darkred', size = 2, lty = 2)

# Combined plot
correlation_plots <- ggarrange(plot_accuracy, plot_cj)

# confidence RT for each possible confidence rating
ggplot() +
  stat_summary(data = df_obs, aes(x = cj, y = rtconf * 1000), fun = 'mean', geom = 'line', size = 1) +
  stat_summary(data = df_obs, aes(x = cj, y = rtconf * 1000), fun = 'mean', geom = 'point', size = 3) +
  
  stat_summary(data = df_pred, aes(x = cj, y = rtconf * 1000), fun = 'mean', geom = 'point', color = '#145369', size = 3) +
  xlab('Confidence') +
  ylab('cRT (ms)') +
  theme_classic()

