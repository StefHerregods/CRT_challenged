################################################################################
# Collapsing confidence bounds analysis on Haddara et al. (2022) data using    #
# the Herregods et al. (2023) model                                            #
################################################################################

rm(list=ls())

# Packages

library(dplyr)
library(Rcpp)  
library(DEoptim)

# Settings

overwrite <- F  # Overwrite already existing files?

z <- 0.5  # Starting point (accuracy-coded dataset -> 0.5)
ntrials <- 1000  # Number of trial simulations per observation
sigma <- 1  # Within-trial noise
dt <- 0.001  # Precision

itermax <- 1000  # Number of DeOptim iterations

make_predictions <- T  # Make model predictions after parameter estimation?
n_predictions <- 1000  # Number of model predictions for every participant 

# Load c++ script

sourceCpp("code/ConfidenceBounds_scale4.cpp") 

# Load and recode data; Remove catch trials

df <- read.csv('data/observations/data_Haddara_2022_Expt2.csv')
df_recoded <- df %>% 
  mutate(sub = Subj_idx,
         cj = Confidence,
         rt = RT_dec,
         rtconf = RT_conf,
         cor = as.numeric(Stimulus == Response)) %>%
  select(sub, cj, rt, rtconf, cor)
df_recoded$manipulation <- 1

# Removing outliers as described in Chen & Rahnev (2023): 
# more than 3 SD of deviation from the mean

df_obs <- df_recoded %>%
  group_by(sub, manipulation) %>%
  filter((rt < mean(rt) + 3 * sd(rt) & rt > mean(rt) - 3 * sd(rt)) & (rtconf < mean(rtconf) + 3 * sd(rtconf) & rtconf > mean(rtconf) - 3 * sd(rtconf))
  )
cat("Trials removed: ", nrow(df_recoded) - nrow(df_obs), ' (', round((nrow(df_recoded) - nrow(df_obs))/nrow(df_recoded),3), '%)', sep = "")

# Add max frequency cj categories for later analyses

df_preprocessed = df_obs %>%
  group_by(sub, manipulation, cj) %>%
  mutate(count = n(),
         experiment = manipulation) %>%
  ungroup() %>%
  group_by(sub, manipulation) %>%
  mutate(max_count = max(count),
         flag = as.numeric(count == max_count),
         cj_category = case_when(flag == 0 ~ 0,
                               flag == 1 ~ as.double(cj)),
         cj_category = max(cj_category)) %>%
  ungroup() %>%
  select(-flag, -max_count, -manipulation)
write.csv(df_preprocessed, "data/observations/data_Haddara_2022_Expt2_preprocessed.csv", row.names = F)

# Function to compute cost

chi_square_optim <- function(params, all_observations, returnFit){  
  
  # Reset chi-square
  
  chiSquare <- 0
  
  # Name parameters
  
  names(params) <- c('v', 'a', 'a_slope', 'ter', 'a2', 'starting_point_confidence', 'postdriftmod', 'a2_slope_upper', 'a2_slope_lower', 'ter2')
  
  observations <- all_observations
  
  # Generate predictions 
  
  predictions <- data.frame(DDM_confidence_bounds(v = params['v'], a = params['a'], a_slope = params['a_slope'], ter = params['ter'], z = z, ntrials = ntrials, s = sigma, dt = dt, a2 = params['a2'], starting_point_confidence = params['starting_point_confidence'], postdriftmod = params['postdriftmod'], a2_slope_upper = params['a2_slope_upper'], a2_slope_lower = params['a2_slope_lower'], ter2 = params['ter2']))
  names(predictions) <- c('rt', 'cor', 'rtconf', 'cj')
  
  # Separate predictions according to the response
  
  c_predicted <- predictions[predictions$cor == 1,]
  e_predicted <- predictions[predictions$cor == 0,]
  
  # Separate predictions according to the cj
  
  conf_predicted_1 <- predictions[predictions$cj == 1,]
  conf_predicted_2 <- predictions[predictions$cj == 2,]
  conf_predicted_3 <- predictions[predictions$cj == 3,]
  conf_predicted_4 <- predictions[predictions$cj == 4,]
  
  # RT data frame
  
  c_predicted_rt <- c_predicted$rt
  e_predicted_rt <- e_predicted$rt
  
  # RTconf data frame (1)
  
  conf_predicted_1_rtconf <- conf_predicted_1$rtconf
  conf_predicted_2_rtconf <- conf_predicted_2$rtconf
  conf_predicted_3_rtconf <- conf_predicted_3$rtconf
  conf_predicted_4_rtconf <- conf_predicted_4$rtconf
  
  # RTconf data frame (2)
  
  c_conf_predicted_rtconf <- c_predicted$rtconf
  e_conf_predicted_rtconf <- e_predicted$rtconf
  
  # If we are only simulating data: return the predictions
  
  if(returnFit==0){ 
    return(predictions)
    
    # If we are fitting the model: compare these predictions to the observations 
    
  }else{ 
    
    ### 1 - Decision RT comparison ###
    
    # Separate observations into correct and error trials
    
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]
    
    # Get the quantile RTs on the "observed data" for correct and error distributions separately (for quantiles .1, .3, .5, .7, .9)
    
    c_quantiles <- quantile(c_observed$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    e_quantiles <- quantile(e_observed$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    
    if (any(is.na(c_quantiles))) {
      c_quantiles <- rep(0,5)
    }
    if (any(is.na(e_quantiles))) {
      e_quantiles <- rep(0,5)
    }
    
    # To combine correct and incorrect trials, we scale the expected interquantile probability by the proportion of correct and incorrect respectively
    
    prop_obs_c <- dim(c_observed)[1] / dim(observations)[1]
    prop_obs_e <- dim(e_observed)[1] / dim(observations)[1]
    
    c_obs_proportion = prop_obs_c * c(.1, .2, .2, .2, .2, .1)
    e_obs_proportion = prop_obs_e * c(.1, .2, .2, .2, .2, .1)
    obs_props <- c(c_obs_proportion, e_obs_proportion)
    
    # Calculate proportion of responses that fall between the observed quantiles when applied to the predicted data 
    
    c_pred_proportion <- c(
      sum(c_predicted_rt <= c_quantiles[1]),
      sum(c_predicted_rt <= c_quantiles[2]) - sum(c_predicted_rt <= c_quantiles[1]),
      sum(c_predicted_rt <= c_quantiles[3]) - sum(c_predicted_rt <= c_quantiles[2]),
      sum(c_predicted_rt <= c_quantiles[4]) - sum(c_predicted_rt <= c_quantiles[3]),
      sum(c_predicted_rt <= c_quantiles[5]) - sum(c_predicted_rt <= c_quantiles[4]),
      sum(c_predicted_rt > c_quantiles[5])
    ) / dim(predictions)[1]
    
    e_pred_proportion <- c(
      sum(e_predicted_rt <= e_quantiles[1]),
      sum(e_predicted_rt <= e_quantiles[2]) - sum(e_predicted_rt <= e_quantiles[1]),
      sum(e_predicted_rt <= e_quantiles[3]) - sum(e_predicted_rt <= e_quantiles[2]),
      sum(e_predicted_rt <= e_quantiles[4]) - sum(e_predicted_rt <= e_quantiles[3]),
      sum(e_predicted_rt <= e_quantiles[5]) - sum(e_predicted_rt <= e_quantiles[4]),
      sum(e_predicted_rt > e_quantiles[5])
    ) / dim(predictions)[1]
    
    pred_props_rt <- c(c_pred_proportion, e_pred_proportion)
    
    # Avoid zeros in the the data (because of division by predictions for chi square statistic) -> set to small number
    
    pred_props_rt[pred_props_rt == 0] <- .0000001
    
    ### 2 - Confidence rating RT comparison (1) ###
    
    # Separate observations into confidence ratings
    
    conf_observed_1 <- observations[observations$cj == 1,]
    conf_observed_2 <- observations[observations$cj == 2,]
    conf_observed_3 <- observations[observations$cj == 3,]
    conf_observed_4 <- observations[observations$cj == 4,]
    
    # Get the quantile confidence RT's on the "observed data" for each distributions separately (for quantiles .1, .3, .5, .7, .9)
    
    conf_quantiles_1 <- quantile(conf_observed_1$rtconf, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    conf_quantiles_2 <- quantile(conf_observed_2$rtconf, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    conf_quantiles_3 <- quantile(conf_observed_3$rtconf, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    conf_quantiles_4 <- quantile(conf_observed_4$rtconf, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    
    if (any(is.na(conf_quantiles_1))) {
      conf_quantiles_1 <- rep(0,5)
    }
    if (any(is.na(conf_quantiles_2))) {
      conf_quantiles_2 <- rep(0,5)
    }
    if (any(is.na(conf_quantiles_3))) {
      conf_quantiles_3 <- rep(0,5)
    }
    if (any(is.na(conf_quantiles_4))) {
      conf_quantiles_4 <- rep(0,5)
    }
    
    # To combine correct and incorrect trials, we scale the expected interquantile probability by the proportion of correct and incorrect respectively
    
    prop_obs_conf_1 <- dim(conf_observed_1)[1] / dim(observations)[1]
    prop_obs_conf_2 <- dim(conf_observed_2)[1] / dim(observations)[1]
    prop_obs_conf_3 <- dim(conf_observed_3)[1] / dim(observations)[1]
    prop_obs_conf_4 <- dim(conf_observed_4)[1] / dim(observations)[1]
    
    conf_obs_1_proportion = prop_obs_conf_1 * c(.1, .2, .2, .2, .2, .1)
    conf_obs_2_proportion = prop_obs_conf_2 * c(.1, .2, .2, .2, .2, .1)
    conf_obs_3_proportion = prop_obs_conf_3 * c(.1, .2, .2, .2, .2, .1)
    conf_obs_4_proportion = prop_obs_conf_4 * c(.1, .2, .2, .2, .2, .1)
    
    conf_obs_props_1 <- c(conf_obs_1_proportion, conf_obs_2_proportion, conf_obs_3_proportion, conf_obs_4_proportion)
    
    # Calculate proportion of responses that fall between the observed quantiles when applied to the predicted data 
    
    conf_1_pred_proportion <- c(
      sum(conf_predicted_1_rtconf <= conf_quantiles_1[1]),
      sum(conf_predicted_1_rtconf <= conf_quantiles_1[2]) - sum(conf_predicted_1_rtconf <= conf_quantiles_1[1]),
      sum(conf_predicted_1_rtconf <= conf_quantiles_1[3]) - sum(conf_predicted_1_rtconf <= conf_quantiles_1[2]),
      sum(conf_predicted_1_rtconf <= conf_quantiles_1[4]) - sum(conf_predicted_1_rtconf <= conf_quantiles_1[3]),
      sum(conf_predicted_1_rtconf <= conf_quantiles_1[5]) - sum(conf_predicted_1_rtconf <= conf_quantiles_1[4]),
      sum(conf_predicted_1_rtconf > conf_quantiles_1[5])
    ) / dim(predictions)[1]
    
    conf_2_pred_proportion <- c(
      sum(conf_predicted_2_rtconf <= conf_quantiles_2[1]),
      sum(conf_predicted_2_rtconf <= conf_quantiles_2[2]) - sum(conf_predicted_2_rtconf <= conf_quantiles_2[1]),
      sum(conf_predicted_2_rtconf <= conf_quantiles_2[3]) - sum(conf_predicted_2_rtconf <= conf_quantiles_2[2]),
      sum(conf_predicted_2_rtconf <= conf_quantiles_2[4]) - sum(conf_predicted_2_rtconf <= conf_quantiles_2[3]),
      sum(conf_predicted_2_rtconf <= conf_quantiles_2[5]) - sum(conf_predicted_2_rtconf <= conf_quantiles_2[4]),
      sum(conf_predicted_2_rtconf > conf_quantiles_2[5])
    ) / dim(predictions)[1]
    
    conf_3_pred_proportion <- c(
      sum(conf_predicted_3_rtconf <= conf_quantiles_3[1]),
      sum(conf_predicted_3_rtconf <= conf_quantiles_3[2]) - sum(conf_predicted_3_rtconf <= conf_quantiles_3[1]),
      sum(conf_predicted_3_rtconf <= conf_quantiles_3[3]) - sum(conf_predicted_3_rtconf <= conf_quantiles_3[2]),
      sum(conf_predicted_3_rtconf <= conf_quantiles_3[4]) - sum(conf_predicted_3_rtconf <= conf_quantiles_3[3]),
      sum(conf_predicted_3_rtconf <= conf_quantiles_3[5]) - sum(conf_predicted_3_rtconf <= conf_quantiles_3[4]),
      sum(conf_predicted_3_rtconf > conf_quantiles_3[5])
    ) / dim(predictions)[1]
    
    conf_4_pred_proportion <- c(
      sum(conf_predicted_4_rtconf <= conf_quantiles_4[1]),
      sum(conf_predicted_4_rtconf <= conf_quantiles_4[2]) - sum(conf_predicted_4_rtconf <= conf_quantiles_4[1]),
      sum(conf_predicted_4_rtconf <= conf_quantiles_4[3]) - sum(conf_predicted_4_rtconf <= conf_quantiles_4[2]),
      sum(conf_predicted_4_rtconf <= conf_quantiles_4[4]) - sum(conf_predicted_4_rtconf <= conf_quantiles_4[3]),
      sum(conf_predicted_4_rtconf <= conf_quantiles_4[5]) - sum(conf_predicted_4_rtconf <= conf_quantiles_4[4]),
      sum(conf_predicted_4_rtconf > conf_quantiles_4[5])
    ) / dim(predictions)[1]
    
    pred_props_rtconf_1 <- c(conf_1_pred_proportion, conf_2_pred_proportion, conf_3_pred_proportion, conf_4_pred_proportion)
    
    # Avoid zeros in the the data (because of division by predictions for chi square statistic) -> set to small number
    
    pred_props_rtconf_1[pred_props_rtconf_1 == 0] <- .0000001
    
    ### 3 - Confidence rating RT comparison (2) ###
    
    # Separate observations into correct and wrong trials
    
    c_conf_observed <- observations[observations$cor == 1,]
    e_conf_observed <- observations[observations$cor == 0,]
    
    # Get the quantile confidence RTs on the "observed data" for correct and error distributions separately (for quantiles .1, .3, .5, .7, .9)
    
    c_conf_quantiles <- quantile(c_conf_observed$rtconf, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    e_conf_quantiles <- quantile(e_conf_observed$rtconf, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    
    if (any(is.na(c_conf_quantiles))) {
      c_conf_quantiles <- rep(0,5)
    }
    if (any(is.na(e_conf_quantiles))) {
      e_conf_quantiles <- rep(0,5)
    }
    
    # To combine correct and incorrect trials, we scale the expected interquantile probability by the proportion of correct and incorrect respectively
    
    prop_obs_c_conf <- dim(c_conf_observed)[1] / dim(observations)[1]
    prop_obs_e_conf <- dim(e_conf_observed)[1] / dim(observations)[1]
    
    c_conf_obs_proportion = prop_obs_c_conf * c(.1, .2, .2, .2, .2, .1)
    e_conf_obs_proportion = prop_obs_e_conf * c(.1, .2, .2, .2, .2, .1)
    conf_obs_props_2 <- c(c_conf_obs_proportion, e_conf_obs_proportion)
    
    # Calculate proportion of responses that fall between the observed quantiles when applied to the predicted data 
    
    c_conf_pred_proportion <- c(
      sum(c_conf_predicted_rtconf <= c_conf_quantiles[1]),
      sum(c_conf_predicted_rtconf <= c_conf_quantiles[2]) - sum(c_conf_predicted_rtconf <= c_conf_quantiles[1]),
      sum(c_conf_predicted_rtconf <= c_conf_quantiles[3]) - sum(c_conf_predicted_rtconf <= c_conf_quantiles[2]),
      sum(c_conf_predicted_rtconf <= c_conf_quantiles[4]) - sum(c_conf_predicted_rtconf <= c_conf_quantiles[3]),
      sum(c_conf_predicted_rtconf <= c_conf_quantiles[5]) - sum(c_conf_predicted_rtconf <= c_conf_quantiles[4]),
      sum(c_conf_predicted_rtconf > c_conf_quantiles[5])
    ) / dim(predictions)[1]
    
    e_conf_pred_proportion <- c(
      sum(e_conf_predicted_rtconf <= e_conf_quantiles[1]),
      sum(e_conf_predicted_rtconf <= e_conf_quantiles[2]) - sum(e_conf_predicted_rtconf <= e_conf_quantiles[1]),
      sum(e_conf_predicted_rtconf <= e_conf_quantiles[3]) - sum(e_conf_predicted_rtconf <= e_conf_quantiles[2]),
      sum(e_conf_predicted_rtconf <= e_conf_quantiles[4]) - sum(e_conf_predicted_rtconf <= e_conf_quantiles[3]),
      sum(e_conf_predicted_rtconf <= e_conf_quantiles[5]) - sum(e_conf_predicted_rtconf <= e_conf_quantiles[4]),
      sum(e_conf_predicted_rtconf > e_conf_quantiles[5])
    ) / dim(predictions)[1]
    
    pred_props_rtconf_2 <- c(c_conf_pred_proportion, e_conf_pred_proportion)
    
    # Avoid zeros in the the data (because of division by predictions for chi square statistic) -> set to small number
    
    pred_props_rtconf_2[pred_props_rtconf_2 == 0] <- .0000001
    
    ### 4 - Calculating cost
    
    # Combine the quantiles for RT's and confidence RT's
    obs_props <- c(obs_props, obs_props, conf_obs_props_1, conf_obs_props_2)
    pred_props <- c(pred_props_rt, pred_props_rt, pred_props_rtconf_1, pred_props_rtconf_2)
    
    # Calculate cost
    
    chiSquare_temp = sum( ( (obs_props - pred_props) ^ 2) )
    
    # Add cost
    
    chiSquare <- chiSquare + chiSquare_temp
    
  }
  
  # Return cost
  
  return(chiSquare)
  
}

# Load data

subs <- unique(df_obs$sub)
N_subs <- length(subs)
experiment_labels <- unique(df_obs$manipulation) 
N_experiments <- length(experiment_labels)

# Optimize (extended) DDM parameters 

for(i in 1:N_subs){  # For each participant separately
  
  print(paste('Running participant', subs[i], 'from', N_subs))
  tempAll <- subset(df_obs, sub == subs[i])
  
  for(c in 1:N_experiments){  # For each condition separately 
    
    tempDat <- subset(tempAll, manipulation == experiment_labels[c])
    tempDat <- tempDat[,c('rt', 'cor', 'cj', 'manipulation', 'rtconf')]
    
    # Load existing individual results if these already exist
    
    file_name <- paste0('data/model_parameters/Haddara2022_Expt2/Haddara2022_Expt2_results_sub_', subs[i], '.Rdata')
    if (overwrite == F & file.exists(file_name)){
      
      load(file_name)
      
      # Else, estimate parameters
      
    } else {
      
      # Optimization function
      
      optimal_params <- DEoptim(chi_square_optim,  # Function to optimize
                                # Possible values for 'v', 'a', 'a_slope', 'ter', 'a2', 'starting_point_confidence', 'postdriftmod', 'a2_slope_upper', 'a2_slope_lower', 'ter2'
                                lower = c(0, .5, 0,  0,     0.0001, 0, 0,  0.0001,  0.0001, -4),  
                                upper = c(3,  4, 0, 1.5,     15,     1, 15, 15,      15,      4),
                                all_observations = tempDat, returnFit = 1, control = c(itermax = itermax))
      
      results <- summary(optimal_params)
      
      # Save individual results
      
      save(results, file = file_name)
      
    }
  }
}

# Combine estimated parameters

df_params <- data.frame(matrix(ncol = 13, nrow = N_subs*length(experiment_labels)))
df_temp <- aggregate(cj_category ~ sub + experiment, data = df_preprocessed, FUN = mean)
colnames(df_params) <- c('sub', 'experiment', 'cj_category', 'v', 'a', 'a_slope', 'ter', 'a2', 'starting_point_confidence', 'postdriftmod', 'a2_slope_upper', 'a2_slope_lower', 'ter2')
row <- 1
for (i in (1:N_subs)){ 
  for(c in 1:N_experiments){
    file_name <- paste0('data/model_parameters/Haddara2022_Expt2/Haddara2022_Expt2_results_sub_', subs[i], '.Rdata')
    load(file_name)
    df_params[row,] <- c(subs[i], experiment_labels[c], df_temp$cj_category[df_temp$sub == subs[i] & df_temp$experiment == experiment_labels[c]], results$optim$bestmem[1], results$optim$bestmem[2], results$optim$bestmem[3], results$optim$bestmem[4], results$optim$bestmem[5], results$optim$bestmem[6], results$optim$bestmem[7], results$optim$bestmem[8], results$optim$bestmem[9], results$optim$bestmem[10])
    row <- row + 1
  }
}
df_params[3:13] <- lapply(df_params[3:13], as.numeric)
write.csv(df_params, file = 'data/model_parameters/Haddara2022_Expt2/Haddara2022_Expt2_params.csv', row.names = F)

# Make predictions for every set of estimated model parameters

if(make_predictions == TRUE){
  pb = txtProgressBar(min = 0, max = nrow(df_params), initial = 0, style = 3) 
  df_predictions <- NULL
  for (set in 1:nrow(df_params)){
    df_predictions_temp <- data.frame(rep(df_params[set,1], each = n_predictions),
                                      rep(df_params[set,2], each = n_predictions),
                                      DDM_confidence_bounds(v = df_params[set,]$v,
                                                            a = df_params[set,]$a,
                                                            a_slope = df_params[set,]$a_slope,
                                                            ter = df_params[set,]$ter,
                                                            a2 = df_params[set,]$a2,
                                                            starting_point_confidence = df_params[set,]$starting_point_confidence,
                                                            postdriftmod = df_params[set,]$postdriftmod,
                                                            a2_slope_upper = df_params[set,]$a2_slope_upper,
                                                            a2_slope_lower = df_params[set,]$a2_slope_lower,
                                                            ter2 = df_params[set,]$ter2,
                                                            s = sigma,
                                                            dt = dt,
                                                            ntrials = n_predictions,
                                                            z = z))
    names(df_predictions_temp) <- c('sub', 'experiment', 'rt', 'cor', 'rtconf', 'cj')
    df_predictions <- rbind(df_predictions, df_predictions_temp)
    setTxtProgressBar(pb, set)
  }
}
df_predictions <- merge(df_predictions, df_temp, by = c('sub', 'experiment'))
write.csv(df_predictions, file = 'data/model_predictions/Haddara2022_Expt2_predictions.csv', row.names = F)




