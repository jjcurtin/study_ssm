library(tidymodels)
library(tidyposterior)
library(data.table)

auc_vals <- read.csv("synthetic_performance_dataframe.csv")

# mod <- perf_mod(auc_vals, seed = 1,iter = 5000,
#                 hetero_var = FALSE, transform = tidyposterior::logit_trans,
#                 chains = 5, adapt_delta=0.99, refresh = 1)
mod <- perf_mod(auc_vals, seed = 1,iter = 5000,
                hetero_var = FALSE,
                chains = 5, adapt_delta=0.99, refresh = 1)


# auc_vals2 <- read.csv("synthetic_performance_dataframe2.csv")
# 
# # mod <- perf_mod(auc_vals, seed = 1,iter = 5000,
# #                 hetero_var = FALSE, transform = tidyposterior::logit_trans,
# #                 chains = 5, adapt_delta=0.99, refresh = 1)
# mod2 <- perf_mod(auc_vals2, seed = 1,iter = 5000,
#                 hetero_var = FALSE,
#                 chains = 5, adapt_delta=0.99, refresh = 1)
print(summary(mod), digits = 3)
# print(summary(mod2), digits = 3)
rstanarm::prior_summary(mod$stan)
preproc_diff <- contrast_models(mod, seed = 1, list_1=c("model_1","model_1"),list_2=c("model_2","model_3"))
summary(preproc_diff, seed = 5,size=0.01)
autoplot(preproc_diff)

extraction <- as.data.frame(mod$stan)
mod1 <- extraction$`(Intercept)`
mod2 <- mod1+extraction$modelmodel_2
mod3 <- mod1+extraction$modelmodel_3
df <- as.data.frame(cbind(mod1,mod2,mod3))
r1_func <- function(input){
  a = input[1]
  b = input[2]
  c = input[3]
  if(a>=b & a >= c){
    return(1)
  }else {
    return(0)
  }
}
df$mod1r1 <- apply(df,1,r1_func)
mean(df$mod1r1)

auc_post <- as.data.frame(tidy(mod)) # get object of class data.frame
ggplot(auc_post, aes(x = posterior, col = model, fill = model)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity", bins = 50) +
  geom_density(alpha=.2) + xlim(min = .7, max = .9) +
  geom_rug()

ggplot(df, aes(x = mod1, y = mod2)) + 
  geom_point(size=2)
# 
# 
# x <- seq(-2, 4, length=10000)
# library('boot')
# xp <- inv.logit(x)
# y <- dnorm(x, mean=1.4, sd=.77)
# plot(xp, y, type="l", lwd=1)
