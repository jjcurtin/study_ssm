library(tidymodels)
library(tidyposterior)
library(data.table)

# Select either auroc or aucprc
# all_auc_vals <- read.csv("auroc_vals.csv")
# export_path <- 'post_auroc_exports.csv'
all_auc_vals <- read.csv("auprc_vals.csv")
export_path <- 'post_auprc_exports.csv'

# Get widths and windows 
widths <- sort(unique(all_auc_vals$width))
windows <- sort(unique(all_auc_vals$window))
exports <- data.frame()

# Loop over windows and widths
for (window in windows){
  window_auc_vals <- all_auc_vals[all_auc_vals$window==window,]
  cat(paste(window,"\n"))
  for (width in widths){
    cat(paste("\t",width, "\n"))
    window_width_auc_vals <- window_auc_vals[window_auc_vals$width==width,]
    # Select relevant subset and rename columns for perf_mod
    aucs <- subset(window_width_auc_vals, select=-c(width,window)) |>
      rename(id=rs_id,id2=fold_id)
    # Run tidyposterior perf_mod function
    mod <- aucs |>
      perf_mod(
        formula = statistic ~ model + (1 | id / id2),
        family = gaussian,
        hetero_var = FALSE, transform = tidyposterior::logit_trans,
        iter = 6000, chains = 4, adapt_delta = 0.999,refresh=0
      )    # defaults: iter = 2000, chains = 4, adapt_delta = 0.8
    # will adjust to address warnings
    # Gather the aurocs for each method
    post_aurocs <- tidy(mod)
    
    # Credible intervals for mean model performance
    mean_perf_cri<-post_aurocs |>
      summary(prob=0.95)
    # Confirm sample count
    sample_count <- dim(post_aurocs)[1]/3
    # Calculate posterior probability that the ssm is rank 1
    pr_ssm_r1<-post_aurocs |> 
      mutate(sample = rep(1:sample_count, each = 3)) |> 
      pivot_wider(names_from = model,
                  values_from = posterior) |> 
      mutate(ssm_best = if_else(ssm > lr & ssm > xgb, 1, 0)) |> 
      pull(ssm_best) |> 
      mean()
    # Credible intervals for mean model performance differences
    mean_perf_diff <- mod |> contrast_models(list_1 = c("ssm", "ssm"),
                           list_2 = c("lr", "xgb"),
                           seed = 123) |>
      summary(prob = 0.95) |>
      select(contrast, probability, mean, lower, upper)
    # Gather information for exporting the performance samples (for use in other plotting)
    iter_exports <- copy(post_aurocs)
    iter_exports$window <- window
    iter_exports$width <- width
    exports <- rbind(exports,iter_exports)
    priors <- rstanarm::prior_summary(mod$stan) 
    # Report information to console
    cat(paste("\t\t Prob. SSM best: ",pr_ssm_r1,"\n"))
    cat(paste("\t\t Mean perf. CrI"))
    for (i in 1:dim(mean_perf_cri)[1]){
      cat(paste("\t\t ",mean_perf_cri[i,],"\n"))
    }
    cat(paste("\t\t, Mean perf. diff. CrI"))
    for (i in 1:dim(mean_perf_diff)[1]){
      cat(paste("\t\t ",mean_perf_diff[i,],"\n"))
    }
    
    cat(paste("\t\t Prior params \n"))
    cat(paste("\t\t intercept",priors$prior_intercept$dist,priors$prior_intercept$location,priors$prior_intercept$adjusted_scale,"\n"))
    cat(paste("\t\t coefs",priors$prior$dist,priors$prior$location,priors$prior$adjusted_scale,"\n"))
    cat(paste("\t\t sigma",priors$prior_aux$dist,priors$prior_aux$adjusted_scale,"\n"))
    
    # # Correlation export (for demonstrating the correlation in mean performance discussed in Appendix)
    # if(window==0 & width==15){
    #   corr_samples <- post_aurocs |> 
    #     mutate(sample = rep(1:sample_count, each = 3)) |> 
    #     pivot_wider(names_from = model,
    #                 values_from = posterior)
    #   write.csv(corr_samples,'corr_samples_0_15.csv',row.names=FALSE)
    # }
    # if(window==0 & width==30){
    #   corr_samples <- post_aurocs |> 
    #     mutate(sample = rep(1:sample_count, each = 3)) |> 
    #     pivot_wider(names_from = model,
    #                 values_from = posterior)
    #   write.csv(corr_samples,'corr_samples_0_30.csv',row.names=FALSE)
    # }
  }
}
# If desired, write out the collected posterior sample information
#write.csv(exports,export_path,row.names=FALSE)