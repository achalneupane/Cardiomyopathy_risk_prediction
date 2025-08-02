## CAD: STJUDE
## this code is to use censored date to re-run all the analysis for CAD models
## repeat the process: prediction_models_cv, select _cv, post_EN_cox, evaluation_prediction 
## time_to_mi should be changed throughout the entire analysis 
rm(list = ls())
library(haven)
library(dplyr)
library(janitor)
library(lubridate)
library(glmnet)
library(stringi)
library(risksetROC)
library(timeROC)
library(RColorBrewer)
source("../utilities.R")

load("results/prediction_model_cox_cad_censordate.Rdata")


cols <- brewer.pal(3, "Dark2")

clinical_model_mi_coef <- coef(clinical_model_mi)
clinical_model_mi_covariates <- names(clinical_model_mi_coef)

prs_model_mi_coef <- coef(prs_model_mi)
prs_model_mi_covariates <- names(prs_model_mi_coef)

CCSSEric_model_mi_coef <- coef(CCSSEric_model_mi)
CCSSEric_model_mi_covariates <- names(CCSSEric_model_mi_coef)
#############################
# Report coefficients
############################
clinical_log_hr <- clinical_model_mi_coef
clinical_hr <- exp(clinical_log_hr)
clinical_coef_hr <- cbind(clinical_log_hr, clinical_hr)

prs_log_hr <- prs_model_mi_coef 
prs_hr <- exp(prs_log_hr)
prs_coef_hr <- cbind(prs_log_hr, prs_hr)

CCSSEric_log_hr <- CCSSEric_model_mi_coef
CCSSEric_hr <- exp(CCSSEric_log_hr)
CCSSEric_coef_hr <- cbind(CCSSEric_log_hr, CCSSEric_hr)


write.csv(round(clinical_coef_hr, 3), "results/clinical_model_mi_coef_censordate.csv")
write.csv(round(prs_coef_hr, 3), "results/prs_model_mi_coef_censordate.csv")
write.csv(round(CCSSEric_coef_hr, 3), "results/CCSSEric_model_mi_coef_censordate.csv")

## continue with the code
## C index and time dependent AUC


#####################
## Evaluate C-index 
#####################
eval_data_mi_df <- analysis_data_final_mi %>% filter(complete.cases(.))
eval_data_mi <- eval_data_mi_df %>% as.matrix


eval_data_mi_y <- with(eval_data_mi_df,
                       Surv(time = time_to_mi, event = mi))


clinical_model_mi_marker <- as.numeric(eval_data_mi[, clinical_model_mi_covariates] %*%
                                         clinical_model_mi_coef)

roc_clinical_mi_model <- risksetROC(Stime = eval_data_mi[, "time_to_mi"],
                                    status = eval_data_mi[, "mi"],
                                    marker = clinical_model_mi_marker,
                                    predict.time = 25,
                                    method = "Cox",lty = 2, col = "red")

auc_clinical_mi_model <- risksetAUC(Stime = eval_data_mi[, "time_to_mi"],
                                    status = eval_data_mi[, "mi"],
                                    marker = clinical_model_mi_marker,
                                    tmax = 25,
                                    method = "Cox",lty = 2, col = "red")

auc_clinical_mi_model$Cindex # 0.7863476



prs_model_mi_marker <- as.numeric(eval_data_mi[, prs_model_mi_covariates] %*%
                                    prs_model_mi_coef)

roc_prs_mi_model <- risksetROC(Stime = eval_data_mi[, "time_to_mi"],
                               status = eval_data_mi[, "mi"],
                               marker = prs_model_mi_marker,
                               predict.time = 25,
                               method = "Cox",lty = 2, col = "red")

auc_prs_mi_model <- risksetAUC(Stime = eval_data_mi[, "time_to_mi"],
                               status = eval_data_mi[, "mi"],
                               marker = prs_model_mi_marker,
                               tmax = 25,
                               method = "Cox",lty = 2, col = "red")

auc_prs_mi_model$Cindex #  0.788902


CCSSEric_model_mi_marker <- as.numeric(eval_data_mi[, CCSSEric_model_mi_covariates] %*%
                                         CCSSEric_model_mi_coef)

roc_CCSSEric_mi_model <- risksetROC(Stime = eval_data_mi[, "time_to_mi"],
                                    status = eval_data_mi[, "mi"],
                                    marker = CCSSEric_model_mi_marker,
                                    predict.time = 25,
                                    method = "Cox",lty = 2, col = "red")

auc_CCSSEric_mi_model <- risksetAUC(Stime = eval_data_mi[, "time_to_mi"],
                                    status = eval_data_mi[, "mi"],
                                    marker = CCSSEric_model_mi_marker,
                                    tmax = 25,
                                    method = "Cox",lty = 2, col = "red")

auc_CCSSEric_mi_model$Cindex #  0.7178797


##############################
## Evaluate Time-dependent AUC 
##############################
clinical_model_mi_ROC <- timeROC(T = eval_data_mi[, "time_to_mi"],
                                 delta = eval_data_mi[, "mi"],
                                 marker = clinical_model_mi_marker,
                                 # other_markers = eval_data[, eval_model_covariates],
                                 # weighting = "cox",
                                 cause = 1,
                                 times = c(10, 20, 25, 40, 45), iid = T) ## time should be less than 49.6, manually add 30

confint(clinical_model_mi_ROC) 

# Time-dependent-Roc curve estimated using IPCW  (n=4145, without competing risks). 
# Cases Survivors Censored AUC (%)   se
# t=10     4      3080     1067   78.52 7.46
# t=20    41      1762     2348   85.77 2.65
# t=25    71      1213     2867   85.57 2.30
# t=40   139       134     3878   71.13 3.56
# t=45   142        37     3972   69.64 5.25


# $CI_AUC
# 2.5% 97.5%
# t=10 63.89 93.14
# t=20 80.58 90.96
# t=25 81.07 90.07
# t=40 64.15 78.10
# t=45 59.35 79.93



prs_model_mi_ROC <- timeROC(T = eval_data_mi[, "time_to_mi"],
                            delta = eval_data_mi[, "mi"],
                            marker = prs_model_mi_marker,
                            #other_markers = eval_data[, selected_covariates],
                            cause = 1,
                            times =  c(10, 20, 25, 40, 45), 
                            weighting = "marginal",
                            iid = T) 

confint(prs_model_mi_ROC)

# Time-dependent-Roc curve estimated using IPCW  (n=4145, without competing risks). 
# Cases Survivors Censored AUC (%)    se
# t=10     4      3080     1067   77.04 9.14
# t=20    41      1762     2348   85.85 2.55
# t=25    71      1213     2867   86.48 2.23
# t=40   139       134     3878   73.18 3.40
# t=45   142        37     3972   72.64 5.39

# $CI_AUC
# 2.5% 97.5%
# t=10 59.13 94.95
# t=20 80.84 90.86
# t=25 82.12 90.85
# t=40 66.51 79.85
# t=45 62.07 83.21

CCSSEric_model_mi_ROC <- timeROC(T = eval_data_mi[, "time_to_mi"],
                                 delta = eval_data_mi[, "mi"],
                                 marker = CCSSEric_model_mi_marker,
                                 #other_markers = eval_data[, selected_covariates],
                                 cause = 1,
                                 times =  c(10, 20, 25, 40, 45), 
                                 weighting = "marginal",
                                 iid = T) 

confint(CCSSEric_model_mi_ROC)

# Time-dependent-Roc curve estimated using IPCW  (n=4145, without competing risks). 
# Cases Survivors Censored AUC (%)   se
# t=10     4      3080     1067   84.28 6.42
# t=20    41      1762     2348   78.69 3.56
# t=25    71      1213     2867   81.96 2.56
# t=40   139       134     3878   70.57 3.53
# t=45   142        37     3972   70.32 5.36


# $CI_AUC
# 2.5% 97.5%
# t=10 71.70 96.86
# t=20 71.72 85.67
# t=25 76.93 86.98
# t=40 63.66 77.49
# t=45 59.81 80.82



compare(prs_model_mi_ROC, clinical_model_mi_ROC)
# $p_values_AUC
# t=10       t=20       t=25       t=40       t=45 
# 0.3933268 0.8991319 0.0476425 0.1242113 0.1830747 

compare(prs_model_mi_ROC, CCSSEric_model_mi_ROC)
# $p_values_AUC
# t=10         t=20         t=25         t=40         t=45 
# 0.008614605 0.015170718 0.003632826 0.062284904 0.197669722 

compare(clinical_model_mi_ROC, CCSSEric_model_mi_ROC)
# t=10         t=20         t=25         t=40         t=45 
# 1.560097e-07 1.503246e-02 2.368769e-02 4.424106e-01 6.276049e-01  



#######################################
## Risk Stratification ################
#######################################
clinical_model_mi_cutoffs <- c(-Inf, quantile(clinical_model_mi_marker, c(0.15, 0.85)), Inf)
prs_model_mi_cutoffs <- c(-Inf, quantile(prs_model_mi_marker, c(0.15, 0.85)), Inf)
CCSSEric_model_mi_cutoffs <- c(-Inf, quantile(CCSSEric_model_mi_marker, c(0.15, 0.85)), Inf)

save(clinical_model_mi_cutoffs, prs_model_mi_cutoffs, CCSSEric_model_mi_cutoffs,
     file = "results/mi_risk_group_cutoffs.RData")
clinical_risk_group <- factor(cut(clinical_model_mi_marker, 
                                  clinical_model_mi_cutoffs,
                                  include.lowest = T, right = F),
                              labels = c("Low risk", "Medium risk", "High risk"))

table(clinical_risk_group)

prs_risk_group <- factor(cut(prs_model_mi_marker, 
                             prs_model_mi_cutoffs,
                             include.lowest = T, right = F),
                         labels = c("Low risk", "Medium risk", "High risk"))
table(prs_risk_group)

CCSSEric_risk_group <- factor(cut(CCSSEric_model_mi_marker, 
                                  CCSSEric_model_mi_cutoffs,
                                  include.lowest = T, right = F),
                              labels = c("Low risk", "Medium risk", "High risk"))

table(CCSSEric_risk_group)

km_data <- data.frame(fevent = eval_data_mi_df$mi,
                      ftime = eval_data_mi_df$time_to_mi,
                      clinical_risk_group = clinical_risk_group,
                      prs_risk_group = prs_risk_group,
                      CCSSEric_risk_group = CCSSEric_risk_group)
# 
# dysl_cuminc_clinical <- with(km_data,
#                              cuminc(ftime, fevent, clinical_risk_group))
# summary_dysl_cuminc_clinical <- get_cuminc_summary(rslt = dysl_cuminc_clinical, 
#                                               times = 25)
# 
# dysl_cuminc_prs <- with(km_data,
#                              cuminc(ftime, fevent, prs_risk_group))
# summary_dysl_cuminc_prs <- get_cuminc_summary(rslt = dysl_cuminc_prs, 
#                                               times = 25)

km_clinical <- survfit(Surv(ftime, fevent) ~ clinical_risk_group,
                       data = km_data)
km_prs <- survfit(Surv(ftime, fevent) ~ prs_risk_group,
                  data = km_data)
km_CCSSEric <- survfit(Surv(ftime, fevent) ~ CCSSEric_risk_group,
                       data = km_data)

clinical_cuminc <- get_survfit_summary(km_clinical, times = 25)
prs_cuminc <- get_survfit_summary(km_prs, times = 25)
CCSSEric_cuminc <- get_survfit_summary(km_CCSSEric, times = 25)


## Make plot of cumulative incidence
## get number at risk
clinical_ss <- summary(km_clinical, times = seq(0, 25, 5))
clinical_nrisk_l <- clinical_ss$n.risk[1:6]
clinical_nrisk_m <- clinical_ss$n.risk[7:12]
clinical_nrisk_h <- clinical_ss$n.risk[13:18]

prs_ss <- summary(km_prs, times = seq(0, 25, 5))
prs_nrisk_l <- prs_ss$n.risk[1:6]
prs_nrisk_m <- prs_ss$n.risk[7:12]
prs_nrisk_h <- prs_ss$n.risk[13:18]

CCSSEric_ss <- summary(km_CCSSEric, times = seq(0, 25, 5))
CCSSEric_nrisk_l <- CCSSEric_ss$n.risk[1:6]
CCSSEric_nrisk_m <- CCSSEric_ss$n.risk[7:12]
CCSSEric_nrisk_h <- CCSSEric_ss$n.risk[13:18]

time_grid <- seq(0, 25, 0.05)
ntimes <- length(time_grid)

png("plots/CAD_cumulative_incidence.png", width = 7, height = 6, units = "in", res = 120)
par(mar = c(5, 5, 1.5, 0.5))
plot(km_clinical, fun = "event", ylim = c(0, 0.2), xlim = c(0, 25), 
     xaxt = "n", yaxt = "n", cex.axis = 1.6, cex.lab = 2, lwd = 2,
     col = cols, lty = 2,
     ylab = "Cumulative incidence", xlab = "Time since cancer diagnosis (years)")
lines(km_prs, col = cols, lty = 1, lwd = 2, fun = "event")
lines(km_CCSSEric, col = cols, lty = 3, lwd = 2, fun = "event")
axis(1, at = seq(0, 25, 5), labels = seq(5, 30, 5), cex.axis = 1.6)
axis(2, at = seq(0, 0.2, 0.05), labels = seq(0, 0.2, 0.05), cex.axis = 1.6)
legend(x = 0, y = 0.2,
       legend = c("Low risk", "Moderate risk", "High risk"),
       fill = cols, bty = "n", xjust = 0, yjust = 1,
       cex = 1.2)
legend(x = 0, y = 0.16,
       legend = c("Clinical Model", "Clinical Model + PRS", 
                  "Chow et al. 2017"),
       lty = c(2, 1, 3), lwd = 2, bty = "n", xjust = 0, yjust = 1,
       cex = 1.2)
dev.off()

# 
# mtext(side = 1, text = "Number at risk", adj = -0.32, line = 1.5, cex = 0.6, font = 2)
# mtext(side = 1, text = "Clinical Model", adj = -0.32, line = 2.2, cex = 0.6, font = 2)
# mtext(side = 1, text = "Low risk", adj = -0.32, line = 2.9, cex = 0.6, font = 2)
# mtext(side = 1, text = "Moderate risk", adj = -0.32, line = 3.6, cex = 0.6, font = 2)
# mtext(side = 1, text = "High risk", adj = -0.32, line = 4.3, cex = 0.6, font = 2)
# 
# mtext(side = 1, text = "Clinical Model + PRS", adj = -0.32, line = 5.5, cex = 0.6, font = 2)
# mtext(side = 1, text = "Low risk", adj = -0.32, line = 6.2, cex = 0.6, font = 2)
# mtext(side = 1, text = "Moderate risk", adj = -0.32, line = 6.9, cex = 0.6, font = 2)
# mtext(side = 1, text = "High risk", adj = -0.32, line = 7.6, cex = 0.6, font = 2)
# 
# mtext(side = 1, text = "Chow et al.", adj = -0.32, line = 8.8, cex = 0.6, font = 2)
# mtext(side = 1, text = "Low risk", adj = -0.32, line = 9.5, cex = 0.6, font = 2)
# mtext(side = 1, text = "Moderate risk", adj = -0.32, line = 10.2, cex = 0.6, font = 2)
# mtext(side = 1, text = "High risk", adj = -0.32, line = 10.9, cex = 0.6, font = 2)
# 
# mtext(side = 1, text="Time (years)", line = 11.9, cex = 0.8, font = 1)
# dev.off()
