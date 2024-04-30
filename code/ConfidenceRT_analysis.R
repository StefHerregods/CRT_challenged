################################################################################
# Analysis of the estimated model parameters and their respective predictions  #
# Code for several plots is adapted code from Bang et al. (2019)               #
################################################################################

rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggpubr)
library(extrafont)

# Settings
font_size = 12

# Separate dataset analyses ----------------------------------------------------

# Choose a dataset to analyse
# [1] Bang et al. (2019) 
# [2] Haddara & Rahnev (2022) - Experiment 1
# [3] Haddara & Rahnev (2022) - Experiment 2
dataset <- 2

# Load data
if (dataset == 1){
  df_obs <- read.csv('data/observations/data_Bang_2019_Exp2_preprocessed.csv')
  df_params <- read.csv('data/model_parameters/Bang2019/Bang2019_params.csv')
  df_pred <- read.csv('data/model_predictions/Bang2019_predictions.csv')
} else if (dataset == 2){
  df_obs <- read.csv('data/observations/data_Haddara_2022_Expt1_preprocessed.csv')
  df_params <- read.csv('data/model_parameters/Haddara2022_Expt1/Haddara2022_Expt1_params.csv')
  df_pred <- read.csv('data/model_predictions/Haddara2022_Expt1_predictions.csv')
} else if (dataset == 3){
  df_obs <- read.csv('data/observations/data_Haddara_2022_Expt2_preprocessed.csv')
  df_params <- read.csv('data/model_parameters/Haddara2022_Expt2/Haddara2022_Expt2_params.csv')
  df_pred <- read.csv('data/model_predictions/Haddara2022_Expt2_predictions.csv')
} 

# Add unique identifier
df_params$id <- row.names(df_params)
df_obs <- merge(df_params[, c('sub', 'experiment', 'id')], df_obs, by = c('sub', 'experiment'))
df_pred <- merge(df_params[, c('sub', 'experiment', 'id')], df_pred, by = c('sub', 'experiment'))




# Investigating the confidence RT vs. confidence judgment correlations
correlations_cj <- NULL
for (id in unique(df_obs$id)){
  correlations_cj <- rbind(correlations_cj, c(cor(df_obs[df_obs$id==id,]$rtconf, df_obs[df_obs$id==id,]$cj, method = 'spearman'), cor(df_pred[df_pred$id==id,]$rtconf, df_pred[df_pred$id==id,]$cj, method = 'spearman')))
}
correlations_cj <- data.frame(correlations_cj)
names(correlations_cj) <- c('observed_correlations', 'predicted_correlations')
plot_cj_corr <- ggplot(data = correlations_cj, aes(x = observed_correlations, y = predicted_correlations)) +
  geom_point(shape = 16, size = 2, alpha = 0.3, color = '#145369') +
  xlim(-1,1) +
  ylim(-1,1) +
  xlab('Observed correlation (*r*)') +
  ylab('Simulated correlation (*r*)') + 
  theme_classic() +
  geom_abline(slope = 1, size = 0.7) +
  geom_abline(intercept = 0, slope = 0, color = 'darkred', size = 0.7, lty = 2) +
  theme(text=element_text(family="sans", size=font_size),
        axis.title.y = ggtext::element_markdown(),
        axis.title.x = ggtext::element_markdown(),
        axis.line = element_line(size = 0.7))

ggsave(filename = 'p2.png', plot = plot_cj_corr, dpi = 500)



plot_model <- function(v, a, ter, z, u_upper, u_lower, tmax = .5, dt = .001){
  # data
  time <- seq(0, tmax, dt)
  a_upper <- a - time*u_upper
  a_lower <- 0 + time*u_lower
  drift_rate <- a*z + (time*v)
  # plot
  plot(0, a*z, type = 'p', col='red', xlim = c(min(0,ter), tmax), ylim = c(0, a), xlab = "Time (s)", ylab = "Evidence", main = "Evidence Accumulation", bty="n")
  polygon(c(0, ter, ter, 0), c(0, 0, a, a), col = adjustcolor("grey", alpha = 0.3)) # non decision time
  lines(time+ter, a_upper); lines(time+ter, a_lower); lines(time+ter, drift_rate) # upper and lower bound and drift rate
  abline(h = a*z, lty = 2) # starting point line
}

id = 49
plot_model(v = df_params$v[df_params$id == id], a = df_params$a[df_params$id == id], ter = df_params$ter[df_params$id == id], z = df_params$starting_point_confidence[df_params$id == id], u_upper = df_params$a2_slope_upper[df_params$id == id], u_lower = df_params$a2_slope_lower[df_params$id == id])
#ggsave('Cj_correlations.png', width = 5.2, height = 5, units = 'cm', scale = 1)




# Combined plot




cor.test(correlations_cj$observed_correlations, correlations_cj$predicted_correlations, use="complete.obs", method = 'pearson')



t.test(atanh(correlations_cj$predicted_correlations[correlations_cj$observed_correlations<=0]),
       alternative = "less")
t.test(atanh(correlations_cj$predicted_correlations[correlations_cj$observed_correlations>0]),
       alternative = "greater")





# ------------------------------------------------------------------------------

# Compute mean values (grouped by individual, cj)
df_obs_mean = df_obs%>%
  group_by(id, cj)%>%
  summarise(cj_category = mean(cj_category),
            rtconf = mean(rtconf), .groups = 'drop')
df_pred_mean = df_pred%>%
  group_by(id, cj)%>%
  summarise(cj_category = mean(cj_category),
            rtconf = mean(rtconf), .groups = 'drop')


# Count participants per group
data_count_obs = df_obs %>%
  group_by(id) %>%
  summarise(cj_category = mean(cj_category), .groups = 'drop') %>%
  group_by(cj_category) %>%
  summarise(count = n(),.groups = 'drop') %>%
  mutate(proportion = count / sum(count))
data_count_pred = df_pred %>%
  group_by(id) %>%
  summarise(cj_category = mean(cj_category), .groups = 'drop') %>%
  group_by(cj_category) %>%
  summarise(count = n(),.groups = 'drop') %>%
  mutate(proportion = count / sum(count))


# Compute mean values (grouped by cj_category, cj)
cRT_cj_categorized_obs = df_obs_mean %>%
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
  mutate(count = rep(data_count_obs$count, each = 4))
cRT_cj_categorized_pred = df_pred_mean %>%
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
  mutate(count = rep(data_count_pred$count, each = 4))



cRT_cj_categorized_obs$source <- 'Observed'
cRT_cj_categorized_pred$source <- 'Simulated'


df_merged <- rbind(cRT_cj_categorized_pred, cRT_cj_categorized_obs)

color_scheme <- c("Observed" = "black", "Simulated" = "red")


ggplot(data = df_merged, aes(x = cj, y = mean, color = source)) +
  facet_grid(. ~ cj_category) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0, size = 0.6)+
  geom_line(size = 1) +
  geom_text(aes(x = 2.5, y = 870, label = paste('N = ', count, sep = '')), size = 4, fontface = 'bold') +
  scale_y_continuous(limits = c(50, 1100)) +
  labs(y = 'cRT (ms)', x = 'Confidence') +
  scale_color_manual(values = color_scheme) +
  theme(text = element_text(size = font_size),
        plot.title = element_text(size = rel(1), hjust = 0.5),
        axis.text = element_text(color = 'black'),
        axis.ticks = element_line(color = 'black'),
        panel.spacing.x = unit(10, 'pt'),
        panel.background = element_rect(color = 'white', fill = 'white'),
        axis.line = element_line(colour = 'black', size = 1),
        strip.background = element_rect(color = 'white', fill = 'white'),
        strip.text = element_text(color = 'black', size = 14),
        legend.position = 'none'
  )

ggsave('Bang_rtconf.png', device = 'png', units = 'cm', width = 10, height = 5)

