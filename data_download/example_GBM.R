######################
### tcga.gbm.RData download from github
  load("./tcga.gbm.RData")
  surv_dat <- download_tcga[[3]]
  if(!("survival" %in% rownames(installed.packages())) ) {
    install.packages("survival")
  }# end fi
  require(survival)
  #idx without time information
  idx_time <- which(is.na(as.numeric(surv_dat$days_to_last_follow_up)) & is.na(as.numeric(surv_dat$days_to_death)))
  #idx without survival information
  idx_sur <- which(is.na(surv_dat$vital_status))
  # idx with wrong survival information
  time <- ifelse(surv_dat$vital_status == "Alive" | "alive", surv_dat$days_to_last_follow_up, surv_dat$days_to_death)
  censor <- ifelse(surv_dat$vital_status == "Alive" | "alive", 0, 1)
  idx_err <- which(is.na(time) | is.na(censor))
  # idx without follwoing up information
  idx_f <- which(time==0)
  idx <- unique(c(idx_time, idx_sur, idx_err, idx_f))
  s <- Surv(time, censor)[-idx]

  ######################
  ##deconvolution with bulk RNAseq data from TCGA database
  exp_mat <- download_tcga[[2]]
  exp_mat <- exp_mat[, -idx]
  ## deconvolute, resulting cell_composition
  num_clust <- 3 # for example clustering 3 groups
  group <- kmeans(cell_composition, num_clust, nstart=10000, iter.max=10000000)$cluster
  
  #################
  chisq_val <- survdiff(s ~ group)$chisq
  logrank <- 1 - pchisq(chisq_val, num_clust-1)
  
  library(survminer)

  sfit <- survfit(s~group, data=dat)
  ggsurvplot(
    sfit, risk.table = TRUE, ggtheme = theme_bw(),
    pval = TRUE, pval.coord = c(0, 0.03)
  )
  
  p <- ggsurvplot(sfit, conf.int=F, pval=TRUE)
  p
  
  