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

# Experimental datasets analyses ----------------------------------------------------

# Choose a dataset to analyse
# [1] Bang et al. (2019) 
# [2] Haddara & Rahnev (2022) - Experiment 1
# [3] Haddara & Rahnev (2022) - Experiment 2
dataset <- 3

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





library(extrafont)
#font_import()
loadfonts(device = "win")

figure_font <- 'Helvetica'

theme_set(
  theme_bw(base_size = 12) +
    theme(
      axis.title = element_text(family = figure_font),       
      axis.text = element_text(size = 10, family = figure_font),   # Set axis text size to 10
      legend.title = element_text(size = 12, family = figure_font),    
      legend.text = element_text(size = 12, family = figure_font),
      legend.background = element_blank(),
      axis.line = element_line(),
      axis.title.y = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(),
      plot.title = element_text(size = 12, family = figure_font, hjust = 0.5))
)


p1 <- ggplot(data = correlations_cj, aes(x = observed_correlations, y = predicted_correlations)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlim(-1,1) +
  ylim(-1,1) +
  xlab('Observed correlation (*r*)') +
  ylab('Simulated correlation (*r*)') + 
  geom_abline(slope = 1, size = 0.5) +
  geom_abline(intercept = 0, slope = 0, color = 'darkred', size = 0.7, lty = 2) +
  annotate("text", x = 1, y = -0.9, 
           label = bquote(italic(r) == 0.67 ~ "," ~ italic(p) < 0.05), 
           family = figure_font, hjust = 1)
p1_title <- p1 + ggtitle('Bang')

p2 <- ggplot(data = correlations_cj, aes(x = observed_correlations, y = predicted_correlations)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlim(-1,1) +
  ylim(-1,1) +
  xlab('Observed correlation (*r*)') +
  ylab('Simulated correlation (*r*)') + 
  geom_abline(slope = 1, size = 0.5) +
  geom_abline(intercept = 0, slope = 0, color = 'darkred', size = 0.7, lty = 2) +
  annotate("text", x = 1, y = -0.9, 
           label = bquote(italic(r) == 0.72 ~ "," ~ italic(p) < 0.05), 
           family = figure_font, hjust = 1)
p2_title <- p2 + ggtitle('Haddara1')

p3 <- ggplot(data = correlations_cj, aes(x = observed_correlations, y = predicted_correlations)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlim(-1,1) +
  ylim(-1,1) +
  xlab('Observed correlation (*r*)') +
  ylab('Simulated correlation (*r*)') + 
  geom_abline(slope = 1, size = 0.5) +
  geom_abline(intercept = 0, slope = 0, color = 'darkred', size = 0.7, lty = 2) +
  annotate("text", x = 1, y = -0.9, 
           label = bquote(italic(r) == 0.74 ~ "," ~ italic(p) < 0.05), 
           family = figure_font, hjust = 1)
p3_title <- p3 + ggtitle('Haddara2')


figure2 <- ggarrange(p1_title, p2_title, p3_title, labels = c('A', 'B', 'C'), font.label = list(size = 12, family = figure_font))




ggsave(filename = 'figure2.png', plot = figure2, dpi = 1000, width = 6, height = 6)







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


p3 <- ggplot(data = df_merged, aes(x = cj, y = mean, color = source)) +
  facet_grid(. ~ cj_category) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0, size = 0.6) +
  geom_line(size = 1) +
  geom_text(aes(x = 2.5, y = 870, label = paste('N = ', count, sep = '')), size = 3.5, family = "Helvetica", fontface = 'plain') +  # Explicitly set family here
  scale_y_continuous(limits = c(50, 1100)) +
  labs(y = 'cRT (ms)', x = 'Confidence') +
  scale_color_manual(values = color_scheme) +
  theme(
    text = element_text(size = 12, family = "Helvetica"),
    plot.title = element_text(size = 12, hjust = 0.5, family = "Helvetica"),
    axis.text = element_text(color = 'black', size = 10, family = "Helvetica"),
    axis.ticks = element_line(color = 'black'),
    panel.spacing.x = unit(10, 'pt'),
    strip.background = element_rect(color = 'black', fill = 'white'),
    strip.text = element_text(color = 'black', size = 12, family = "Helvetica"),
    legend.position = 'none'
  )

plot <- ggarrange(p1, p2, p3, ncol = 2, nrow = 2)

ggsave('rtconf.png', device = 'png', units = 'in', width = 7.5, height = 3.5, dpi = 1000)


# Parameter recovery analysis --------------------------------------------------


# Load data sets
params_pred <- read.csv("data/model_parameters/model_recovery/model_recovery_separateTer2_params.csv")
params_obs <- read.csv("data/model_parameters/model_recovery/parameters_recovery.csv")

# merge data sets
names(params_pred) <- paste0(names(params_pred), "_pred")
names(params_obs) <- paste0(names(params_obs), "_obs")
df <- cbind(params_pred, params_obs)






plot(params_pred$a, params_obs$a)
plot(params_pred$v, params_obs$v)
plot(params_pred$ter, params_obs$ter)
plot(params_pred$a2, params_obs$a2, xlab="", ylab="")
title(main = 'a2')
cor(params_pred$a2, params_obs$a2)

plot(params_pred$a2_slope_upper, params_obs$a2_slope_upper, xlab="", ylab="")
title(main = 'a2_slope_upper')

cor(params_pred$a2_slope_upper, params_obs$a2_slope_upper)

plot(params_pred$a2_slope_lower, params_obs$a2_slope_lower, xlab="", ylab="")
title(main = 'a2_slope_lower')

cor(params_pred$a2_slope_lower, params_obs$a2_slope_lower)

plot(params_pred$postdriftmod, params_obs$postdriftmod, xlab="", ylab="")
cor(params_pred$postdriftmod, params_obs$postdriftmod)

title(main = 'vratio')

plot(params_pred$ter2_1, params_obs$ter2_1)
title(main = 'ter_1')

plot(params_pred$ter2_2, params_obs$ter2_2)
title(main = 'ter_2')

plot(params_pred$ter2_3, params_obs$ter2_3)
title(main = 'ter_3')

plot(params_pred$ter2_4, params_obs$ter2_4)
title(main = 'ter_4')

par(mfrow=c(2,2))

names(params_pred)
names(params_obs)


cor(params_pred$ter2_4, params_obs$ter2_4)







cor(df$a_obs, df$a_pred)
cor(df$v_obs, df$v_pred)
cor(df$postdriftmod_obs, df$postdriftmod_pred)
cor(df$ter_obs, df$ter_pred)
cor(df$starting_point_confidence_obs, df$starting_point_confidence_pred)
cor(df$a2_obs, df$a2_pred)
cor(df$a2_slope_upper_obs, df$a2_slope_upper_pred)
cor(df$a2_slope_lower_obs, df$a2_slope_lower_pred)
cor(df$ter2_1_obs, df$ter2_1_pred)
cor(df$ter2_2_obs, df$ter2_2_pred)
cor(df$ter2_3_obs, df$ter2_3_pred)
cor(df$ter2_4_obs, df$ter2_4_pred)



theme_set(
  theme_bw(base_size = 11) +
    theme(
      axis.title = element_text(family = figure_font),       
      axis.text = element_text(size = 10, family = figure_font),   # Set axis text size to 10
      legend.title = element_text(size = 10, family = figure_font),    
      legend.text = element_text(size = 10, family = figure_font),
      legend.background = element_blank(),
      axis.line = element_line(),
      axis.title.y = ggtext::element_markdown(),
      axis.title.x = ggtext::element_markdown(),
      plot.title = element_text(size = 11, family = figure_font, hjust = 0.5))
)


a <- ggplot(data = df, aes(x = a_obs, y = a_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle('Decision boundary\nseparation') +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(0.5, 2, by = 0.5), limits = c(0.5, 2)) +
  scale_y_continuous(breaks = seq(0.5, 2, by = 0.5), limits = c(0.5, 2)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))
v <- ggplot(data = df, aes(x = v_obs, y = v_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle('Drift rate\n') +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(0, 2, by = 1), limits = c(0, 2.1)) +
  scale_y_continuous(breaks = seq(0, 2, by = 1), limits = c(0, 2.1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))
ter <- ggplot(data = df, aes(x = ter_obs, y = ter_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle('Non-decision time\n') +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.2), limits = c(0, 0.55)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.2), limits = c(0, 0.55)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))  
a2 <- ggplot(data = df, aes(x = a2_obs, y = a2_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle('Confidence boundary\nseparation') +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10.1)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10.1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))
a2_slope_upper <- ggplot(data = df, aes(x = a2_slope_upper_obs, y = a2_slope_upper_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle('Upper confidence\nboundary slope') +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(0, 7, by = 1), limits = c(0, 7.1)) +
  scale_y_continuous(breaks = seq(0, 7, by = 1), limits = c(0, 7.1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))
a2_slope_lower <- ggplot(data = df, aes(x = a2_slope_lower_obs, y = a2_slope_lower_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle('Lower confidence\nboundary slope') +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(0, 7, by = 1), limits = c(0, 7.1)) +
  scale_y_continuous(breaks = seq(0, 7, by = 1), limits = c(0, 7.1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) 
starting_point_confidence <- ggplot(data = df, aes(x = starting_point_confidence_obs, y = starting_point_confidence_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle('Post-decision starting\npoint bias') +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(0.2, 0.8, by = 0.2), limits = c(0.15, 0.85)) +
  scale_y_continuous(breaks = seq(0.2, 0.8, by = 0.2), limits = c(0.15, 0.85)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))
vratio <- ggplot(data = df, aes(x = postdriftmod_obs, y = postdriftmod_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle(expression(atop(italic('v') * '-ratio'), paste(''))) +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(0, 5, by = 1), limits = c(0, 5)) +
  scale_y_continuous(breaks = seq(0, 5, by = 1), limits = c(0, 5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))
ter2_1 <- ggplot(data = df, aes(x = ter2_1_obs, y = ter2_1_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle('Non-decision time\n(cj = 1)') +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(-2, 2, by = 1), limits = c(-2.2, 2.2)) +
  scale_y_continuous(breaks = seq(-2, 2, by = 1), limits = c(-2.2, 2.2)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11))  
ter2_2 <- ggplot(data = df, aes(x = ter2_2_obs, y = ter2_2_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle('Non-decision time\n(cj = 2)') +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(-2, 2, by = 1), limits = c(-2.2, 2.2)) +
  scale_y_continuous(breaks = seq(-2, 2, by = 1), limits = c(-2.2, 2.2)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) 
ter2_3 <- ggplot(data = df, aes(x = ter2_3_obs, y = ter2_3_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle('Non-decision time\n(cj = 3)') +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(-2, 2, by = 1), limits = c(-2.2, 2.2)) +
  scale_y_continuous(breaks = seq(-2, 2, by = 1), limits = c(-2.2, 2.2)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) 
ter2_4 <- ggplot(data = df, aes(x = ter2_4_obs, y = ter2_4_pred)) +
  geom_point(shape = 16, size = 2, alpha = 0.3) +
  xlab('Observed') +
  ylab('Simulated') + 
  ggtitle('Non-decision time\n(cj = 4)') +
  coord_fixed() +  
  scale_x_continuous(breaks = seq(-2, 2, by = 1), limits = c(-2.2, 2.2)) +
  scale_y_continuous(breaks = seq(-2, 2, by = 1), limits = c(-2.2, 2.2)) +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) 

recovery_plots <- ggarrange(a, v, ter, vratio, starting_point_confidence, a2, a2_slope_upper, a2_slope_lower, ter2_1, ter2_2, ter2_3, ter2_4, ncol = 4, nrow = 3, align = "hv")


ggsave('figures/ParamRecovery.png', recovery_plots, device = 'png', units = 'in', width = 7.5, height = 6.5, dpi = 1000)




