library(ggplot2)

params_pred <- read.csv("data/model_parameters/model_recovery/model_recovery_separateTer2_params.csv")
params_obs <- read.csv("data/model_parameters/model_recovery/parameters_recovery.csv")

params_pred$type <- 'pred'
params_obs$type <- 'obs'
params_obs$set <- 0
params_pred$X <- 0

df <- rbind(params_pred, params_obs)

ggplot(data = params_pred, )


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

