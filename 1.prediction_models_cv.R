rm(list = ls())
library(lubridate)
library(survival)
library(mgcv)
library(dplyr)
library(splines)
library(glmnet)
source("../utilities.R")

seed <- 2000
set.seed(seed)


#################################################
###########  MI models   ########################
#################################################
## replace lstvisitdt by censordate

mi_data <- readRDS("results/mi_complete_data.rds")

mi_data$time_to_mi <-
  with(mi_data, 
       ifelse(mi == 1, 
              yrdif(eligdt, mi_dt),
              yrdif(eligdt, censordate)))

analysis_data_mi <- mi_data %>%
  dplyr::transmute(mrn = mrn, 
                   agedx = agedx,
                   male = as.numeric(gender == "Male"),
                   ## white = as.numeric(racegrp == "White"), ## this is demographic variable, 3382 white vs others 
                   white = white, ## this is ancestry variable, 3287 white vs others 
                   agedx2 = as.numeric(agedx > 5 & agedx <= 10),
                   agedx3 = as.numeric(agedx > 10 & agedx <= 15),
                   agedx4 = as.numeric(agedx > 15),
                   cardiomyopathy = cardiomyopathy,
                   hypertension = hypertension, 
                   dyslipidemia = dyslipidemia,
                   diabetes = diabetes,
                   alkylating = aaclassic_5,
                   carboplatin = carboplatin_5,
                   cisplatin = cisplatin_5,
                   vinca = vinca_5,
                   segrt = as.numeric(maxsegrtdose > 0),
                   maxsegrtdose1 = as.numeric(maxsegrtdose > 0 & maxsegrtdose < 2000),
                   maxsegrtdose2 = as.numeric(maxsegrtdose >= 2000 & maxsegrtdose < 3000),
                   maxsegrtdose3 = as.numeric(maxsegrtdose >= 3000 & maxsegrtdose < 5000),
                   maxsegrtdose4 = as.numeric(maxsegrtdose >= 5000),
                   maxneckrtdose1 = as.numeric(maxneckrtdose > 0),
                   maxchestrtdose1 = as.numeric(maxchestrtdose > 0 &  maxchestrtdose < 500),
                   maxchestrtdose2 = as.numeric(maxchestrtdose >= 500 & maxchestrtdose < 1500),
                   maxchestrtdose3 = as.numeric(maxchestrtdose >= 1500 & maxchestrtdose < 3500),
                   maxchestrtdose4 = as.numeric(maxchestrtdose >= 3500),
                   score3725_z, score0710_z,
                   time_to_mi = time_to_mi,
                   mi
                   ## time_to_mi= pmin(time_to_mi, 40), 
                   ## mi = ifelse(time_to_mi > 40, 0, mi)  censor later
  )
## censoring at certain years after cancer diagnosis 

analysis_data_x <- analysis_data_mi %>%
  dplyr::select(male, agedx2:cisplatin, maxsegrtdose1:maxsegrtdose4,
                maxneckrtdose1:score0710_z)

analysis_data_y <- analysis_data_mi %>%
  dplyr::select(time_to_mi, mi)


analysis_data_final_mi <- cbind(mrn = analysis_data_mi$mrn,
                                analysis_data_x,
                                analysis_data_y)

## N = 4151
clinical_model_data_mi <- analysis_data_final_mi %>%
  filter(complete.cases(.)) %>% 
  dplyr::select(-mrn,
                -score3725_z, -score0710_z) ## no prs

clinical_model_mi_x <- clinical_model_data_mi %>%
  dplyr::select(-time_to_mi, -mi) %>%
  as.matrix

clinical_model_mi_y <- with(clinical_model_data_mi, 
                            Surv(time = time_to_mi,
                                 event = mi))

aa <- seq(0, 1, 0.05)
clinical_model_list <- vector("list", length = length(aa))

for (i in seq_along(aa)) {
  clinical_model_list[[i]] <- cv.glmnet(x = clinical_model_mi_x,
                                        y = clinical_model_mi_y,
                                        family = "cox",
                                        type.measure = "C",
                                        alpha = aa[i])
}

# prs_names <- c("score3725", "score0710")
prs_model_data_mi <- analysis_data_final_mi %>%
  filter(complete.cases(.)) ## 4145 


prs_model_mi_x <- prs_model_data_mi %>%
  dplyr::select(-mrn, -time_to_mi, -mi) %>%
  as.matrix

prs_model_mi_y <- with(prs_model_data_mi,
                       Surv(time = time_to_mi, event = mi))


aa <- seq(0, 1, 0.05)
prs_model_list <- vector("list", length = length(aa))
for (i in seq_along(aa)) {
  prs_model_list[[i]] <- cv.glmnet(x = prs_model_mi_x,
                                   y = prs_model_mi_y,
                                   family = "cox",
                                   type.measure = "C",
                                   alpha = aa[i])
  
}

clinical_coef <- prs_coef <- vector("list", 5)
clinical_order <- order(sapply(clinical_model_list, function(x) min(x$cvm)))
for (i in 1:5) {
  clinical_coef[[i]] <- coef(clinical_model_list[[clinical_order[i]]])
}

prs_order <- order(sapply(prs_model_list, function(x) min(x$cvm)))
for (i in 1:5) {
  prs_coef[[i]] <- coef(prs_model_list[[prs_order[i]]])
}


## display clinical_coef and prs_coef to select coefficients 

clinical_model_mi <- coxph(Surv(time_to_mi, mi) ~
                             male + ## white + agedx2 + agedx3 + 
                             agedx2 + agedx3 + agedx4 + 
                             # cardiomyopathy  + ## hypertension + 
                             # dyslipidemia + ## diabetes +
                             ## alkylating + carboplatin + 
                             cisplatin + 
                             #maxsegrtdose1 + maxsegrtdose2 +
                             #maxsegrtdose3 + maxsegrtdose4 +
                             #maxneckrtdose1 + 
                             maxchestrtdose1 + maxchestrtdose2 +
                             maxchestrtdose3 + maxchestrtdose4,
                           data = analysis_data_final_mi)

prs_model_mi <- coxph(Surv(time_to_mi, mi) ~
                        male + ## white + agedx2 + agedx3 + 
                        agedx2 + agedx3 + agedx4 + 
                        # cardiomyopathy  + ## hypertension + 
                        # dyslipidemia + ## diabetes +
                        ## alkylating + carboplatin + 
                        cisplatin + 
                        # maxsegrtdose1 + maxsegrtdose2 +
                        # maxsegrtdose3 + maxsegrtdose4 +
                        ## maxneckrtdose1 not selected 
                        maxchestrtdose1 + maxchestrtdose2 +
                        maxchestrtdose3 + maxchestrtdose4 +
                        score3725_z ## these two prs highly correlated, + score0710_s
                      ,
                      data = analysis_data_final_mi)

## we notice that censordate selected the exact same variables as lstvisitdt

CCSSEric_model_mi <- coxph(Surv(time_to_mi, mi) ~
                             male + agedx2 + agedx3 + agedx4 + 
                             maxchestrtdose1 + maxchestrtdose2 + 
                             maxchestrtdose3 + maxchestrtdose4,
                           data = analysis_data_final_mi)


save(analysis_data_final_mi, clinical_model_mi, prs_model_mi, CCSSEric_model_mi,
     file = sprintf("results/prediction_model_cox_cad_censordate.Rdata"))




## this code is to use censored date to re-run all the analysis for stroke models
## repeat the code: prediction_coefficients_calculation_stroke.R 
## time_to_stroke should be changed throughout the entire analysis

#################################################
###########  Stroke models   ########################
#################################################
stroke_data <- readRDS("results/stroke_complete_data.rds")

stroke_data$time_to_stroke <-
  with(stroke_data, 
       ifelse(stroke == 1, 
              yrdif(eligdt, stroke_dt),
              yrdif(eligdt, censordate)))

analysis_data_stroke <- stroke_data %>%
  dplyr::transmute(mrn = mrn, 
                   male = as.numeric(gender == "Male"),
                   ## white = as.numeric(racegrp == "White"), ## this is demographic variable, 3382 white vs others 
                   white = white, ## this is ancestry variable, 3287 white vs others 
                   agedx2 = as.numeric(agedx > 5 & agedx <= 10),
                   agedx3 = as.numeric(agedx > 10 & agedx <= 15),
                   agedx4 = as.numeric(agedx > 15),
                   cardiomyopathy = cardiomyopathy,
                   hypertension = hypertension, 
                   dyslipidemia = dyslipidemia,
                   diabetes = diabetes,
                   alkylating = aaclassic_5,
                   carboplatin = carboplatin_5,
                   cisplatin = cisplatin_5,
                   vinca = vinca_5, ## Eric's paper
                   anthra_dose1 = as.numeric(anthracyclines_dose_5 > 0 &
                                               anthracyclines_dose_5 < 100),
                   anthra_dose2 = as.numeric(anthracyclines_dose_5 >= 100 &
                                               anthracyclines_dose_5 < 250),
                   anthra_dose3 = as.numeric(anthracyclines_dose_5 >= 250),
                   segrt = as.numeric(maxsegrtdose > 0),
                   chestrt = as.numeric(maxchestrtdose > 0),
                   methotrexate = methotrexate_5,
                   maxsegrtdose1 = as.numeric(maxsegrtdose > 0 & maxsegrtdose < 2000),
                   maxsegrtdose2 = as.numeric(maxsegrtdose >= 2000 & maxsegrtdose < 3000),
                   maxsegrtdose3 = as.numeric(maxsegrtdose >= 3000 & maxsegrtdose < 5000),
                   maxsegrtdose4 = as.numeric(maxsegrtdose >= 5000),
                   maxneckrtdose1 = as.numeric(maxneckrtdose > 0),
                   maxchestrtdose1 = as.numeric(maxchestrtdose > 0 &  maxchestrtdose < 500),
                   maxchestrtdose2 = as.numeric(maxchestrtdose >= 500 & maxchestrtdose < 1500),
                   maxchestrtdose3 = as.numeric(maxchestrtdose >= 1500 & maxchestrtdose < 3500),
                   maxchestrtdose4 = as.numeric(maxchestrtdose >= 3500),
                   score2724_z,
                   time_to_stroke = time_to_stroke, stroke
                   ## time_to_stroke= pmin(time_to_stroke, 40), 
                   ## stroke = ifelse(time_to_stroke > 40, 0, stroke)  censor later
  )
## censoring at certain years after cancer diagnosis 

with(analysis_data_stroke, table(maxneckrtdose1, segrt))
with(analysis_data_stroke, table(chestrt, segrt))
## almost everyone with cranial radiation also have nonzero neck radiation


analysis_data_x <- analysis_data_stroke %>%
  dplyr::select(male, agedx2:cisplatin, vinca, maxsegrtdose1:maxsegrtdose4,
                maxchestrtdose1:score2724_z)

analysis_data_y <- analysis_data_stroke %>%
  dplyr::select(time_to_stroke, stroke)


analysis_data_final_stroke <- cbind(mrn = analysis_data_stroke$mrn,
                                    analysis_data_x,
                                    analysis_data_y)

## N = 4151
clinical_model_data_stroke <- analysis_data_final_stroke %>%
  filter(complete.cases(.)) %>% 
  dplyr::select(-mrn,
                -score2724_z) ## no prs

clinical_model_stroke_x <- clinical_model_data_stroke %>%
  dplyr::select(-time_to_stroke, -stroke) %>%
  as.matrix

clinical_model_stroke_y <- with(clinical_model_data_stroke, 
                                Surv(time = time_to_stroke,
                                     event = stroke))

aa <- seq(0, 1, 0.05)
clinical_model_list <- vector("list", length = length(aa))

for (i in seq_along(aa)) {
  clinical_model_list[[i]] <- cv.glmnet(x = clinical_model_stroke_x,
                                        y = clinical_model_stroke_y,
                                        family = "cox",
                                        type.measure = "C",
                                        alpha = aa[i])
}


prs_model_data_stroke <- analysis_data_final_stroke %>%
  filter(complete.cases(.)) ## 4151 


prs_model_stroke_x <- prs_model_data_stroke %>%
  dplyr::select(-mrn, -time_to_stroke, -stroke) %>%
  as.matrix


## minus four comorbidity
## prs_model_stroke_x2 <- prs_model_data_stroke %>%
##  dplyr::select(-mrn, -time_to_stroke, -stroke, -cardiomyopathy, -hypertension, -dyslipidemia, -diabetes) %>%
##  as.matrix

prs_model_stroke_y <- with(prs_model_data_stroke,
                           Surv(time = time_to_stroke, event = stroke))



aa <- seq(0, 1, 0.05)
prs_model_list <- vector("list", length = length(aa))
for (i in seq_along(aa)) {
  prs_model_list[[i]] <- cv.glmnet(x = prs_model_stroke_x,
                                   y = prs_model_stroke_y,
                                   family = "cox",
                                   type.measure = "C",
                                   alpha = aa[i])
  
}

## prs_model_list2 <- vector("list", length = length(aa))
## for (i in seq_along(aa)) {
##  prs_model_list2[[i]] <- cv.glmnet(x = prs_model_stroke_x2,
##                                   y = prs_model_stroke_y,
##                                   family = "cox",
##                                   type.measure = "C",
##                                   alpha = aa[i])
##}


clinical_coef <- prs_coef  <- vector("list", 5)
clinical_order <- order(sapply(clinical_model_list, function(x) min(x$cvm)))
for (i in 1:5) {
  clinical_coef[[i]] <- coef(clinical_model_list[[clinical_order[i]]])
}

prs_order <- order(sapply(prs_model_list, function(x) min(x$cvm)))
for (i in 1:5) {
  prs_coef[[i]] <- coef(prs_model_list[[prs_order[i]]])
}

## prs_order2 <- order(sapply(prs_model_list2, function(x) min(x$cvm)))
## for (i in 1:5) {
##  prs_coef2[[i]] <- coef(prs_model_list2[[prs_order2[i]]])
## }

clinical_coef 
prs_coef
## prs_coef2

## we notice that censordate selected the exact same variables as lstvisitdt
clinical_model_stroke <- coxph(Surv(time_to_stroke, stroke) ~
                                 dyslipidemia +
                                 maxsegrtdose1 + maxsegrtdose2 + maxsegrtdose3 + maxsegrtdose4 
                               #maxneckrtdose1
                               ,
                               data = analysis_data_final_stroke)

prs_model_stroke <- coxph(Surv(time_to_stroke, stroke) ~
                            #hypertension + 
                            dyslipidemia +
                            maxsegrtdose1 + maxsegrtdose2 + maxsegrtdose3 + maxsegrtdose4 +
                            ## maxneckrtdose1 + 
                            #maxchestrtdose1 +
                            score2724_z, ## force this variable in, suggested by PI
                          data = analysis_data_final_stroke)


## Eric's model in A1 of the manuscript, predictors in the standard model are 
## sex (not significant), age (not significant), cranial rt, chest rt
## neck was not listed in Eric's paper for stroke
## alkylating (significant), vinca (not significant), 

CCSSEric_model_stroke <- coxph(Surv(time_to_stroke, stroke) ~
                                 ## male +
                                 ## agedx4 + 
                                 ## vinca + alkylating + ## put this in Erik's stroke model 2
                                 maxsegrtdose1 + maxsegrtdose2  + maxsegrtdose3 + maxsegrtdose4 +
                                 maxchestrtdose3 + maxchestrtdose4
                                 ## maxneckrtdose1 + ## neck is selected 
                                 ,
                               data = analysis_data_final_stroke)

CCSSEric_model_stroke2 <- coxph(Surv(time_to_stroke, stroke) ~
                                  ## male +
                                  ## agedx4 + 
                                  vinca + alkylating + ## should we add this to the model, they are not significant
                                  maxsegrtdose1 + maxsegrtdose2  + maxsegrtdose3 + maxsegrtdose4 +
                                  maxchestrtdose3 + maxchestrtdose4,
                                data = analysis_data_final_stroke)

save(analysis_data_final_stroke, clinical_model_stroke, prs_model_stroke, CCSSEric_model_stroke, CCSSEric_model_stroke2,
     file = sprintf("results/prediction_model_cox_stroke_censordate.Rdata"))
