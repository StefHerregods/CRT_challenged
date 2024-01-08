################################################################################
# Analysis of the estimated model parameters and their respective predictions  #
# Code for several plots is adapted code from Bang et al. (2019)               #
################################################################################

rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)

# Settings
font_size = 18

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

# Add unique identifier
df_params$id <- row.names(df_params)
df_obs <- merge(df_params[, c('sub', 'experiment', 'id')], df_obs, by = c('sub', 'experiment'))
df_pred <- merge(df_params[, c('sub', 'experiment', 'id')], df_pred, by = c('sub', 'experiment'))

# Investigating the confidence RT vs. accuracy correlations
correlations_accuracy <- NULL
for (id in unique(df_obs$id)){
  correlations_accuracy <- rbind(correlations_accuracy, c(cor(df_obs[df_obs$id==id,]$rtconf, df_obs[df_obs$id==id,]$cor), cor(df_pred[df_pred$id==id,]$rtconf, df_pred[df_pred$id==id,]$cor)))
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
for (id in unique(df_obs$id)){
  correlations_cj <- rbind(correlations_cj, c(cor(df_obs[df_obs$id==id,]$rtconf, df_obs[df_obs$id==id,]$cj), cor(df_pred[df_pred$id==id,]$rtconf, df_pred[df_pred$id==id,]$cj)))
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
ggarrange(plot_accuracy, plot_cj)

# Confidence RT for each possible confidence rating
ggplot() +
  stat_summary(data = df_obs, aes(x = cj, y = rtconf * 1000), fun = 'mean', geom = 'line', size = 1) +
  stat_summary(data = df_obs, aes(x = cj, y = rtconf * 1000), fun = 'mean', geom = 'point', size = 3) +
  
  stat_summary(data = df_pred, aes(x = cj, y = rtconf * 1000), fun = 'mean', geom = 'point', color = '#145369', size = 3) +
  xlab('Confidence') +
  ylab('cRT (ms)') +
  theme_classic()


options(repr.plot.width = 4, repr.plot.height = 2.5, warn = -1)

data_count = df_obs %>%
  group_by(id) %>%
  summarise(cj_category = mean(cj_category), .groups = 'drop') %>%
  group_by(cj_category) %>%
  summarise(count = n(),.groups = 'drop') %>%
  mutate(proportion = count / sum(count))

cRT_cj_categorized = df_obs %>%
  mutate(cj_category = case_when(cj_category == 1 ~ 'Group 1',
                                 cj_category == 2 ~ 'Group 2',
                                 cj_category == 3 ~ 'Group 3',
                                 cj_category == 4 ~ 'Group 4'),
         rtconf = rtconf * 1000) %>%
  group_by(cj_category, cj) %>%
  summarise(mean = mean(rtconf), se = sd(rtconf) / sqrt(n()), .groups = 'drop') %>%
  group_by(cj_category) %>%
  mutate(min_RT_conf = min(mean), count = n()) %>%
  ungroup() %>%
  mutate(count = rep(data_count$count, each = 4)) %>%
  ggplot(aes(x = cj, y = mean)) +
  facet_grid(. ~ cj_category) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  #geom_line(data = df_pred, aes(x = cj, y = rtconf * 1000)) +
  geom_text(aes(x = 2.2, y = 670, label = paste('N = ', count, sep = '')), size = 4) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0, size = 0.6) +
  geom_hline(aes(yintercept = min_RT_conf), linetype = 2, size = 0.8) +
  scale_y_continuous(limits = c(290, 700)) +
  labs(y = 'cRT (ms)') +
  theme(text = element_text(size = font_size),
        plot.title = element_text(size = rel(1), hjust = 0.5),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        panel.spacing.x = unit(10, 'pt'),
        panel.background = element_rect(color = 'white', fill = 'white'),
        axis.line = element_line(colour = 'black', size = 1),
        strip.background = element_rect(color = 'white', fill = 'white'),
        strip.text = element_text(color = 'black', size = 14),
        legend.position = 'bottom')







