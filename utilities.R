impute_NA <- function(x, impute) {
  ifelse(is.na(x), impute, x)
}

## for binary variables, show counts and percentage
get_N_percentage <- function(x, digits = 1) {
  if (any(is.na(x))) {
    warning(sprintf("%s observation(s) missing", sum(is.na(x))))
  }
  paste0(sum(x, na.rm = T), " (", format(round(mean(x, na.rm = T) * 100, digits), nsmall = digits), "%)")
}

get_percentage_N <- function(x, digits = 1) {
  if (any(is.na(x))) {
    warning(sprintf("%s observation(s) missing", sum(is.na(x))))
  }
  paste0(format(round(mean(x, na.rm = T) * 100, digits), nsmall = digits),
         "% (", sum(x, na.rm = T), ")")
}

get_mean_sd <- function(x, digits = 2) {
  if (any(is.na(x))) {
    warning(sprintf("%s observation(s) missing", sum(is.na(x))))
  }
  paste0(format(round(mean(x, na.rm = T), digits), nsmall = digits), " (",
         format(round(sd(x, na.rm = T), digits), nsmall = digits), ")")
}

get_median_range <- function(x, digits = 2) {
  if (any(is.na(x))) {
    warning(sprintf("%s observation(s) missing", sum(is.na(x))))
  }
  paste0(format(round(median(x, na.rm = T), digits), nsmall = digits), " (",
         format(round(range(x, na.rm = T)[1], digits), nsmall = digits), ", ",
         format(round(range(x, na.rm = T)[2], digits), nsmall = digits), ")")
}

# p-values of Pearson Chi-squared test for group comparison
get_pval_pearson <- function(x, group_var, levs, alternative = "two.sided",
                             adjust = FALSE, digits = 3){
  x <- as.numeric(x)
  pvals <- pval_num <- rep(NA, choose(length(levs), 2))
  
  k <- 0
  for (i in 1:(length(levs) - 1)) {
    for (j in (i + 1):length(levs)) {
      k <- k + 1
      xi <- c(sum(x[group_var == levs[i]], na.rm = T), sum(x[group_var == levs[j]], na.rm = T))
      ni <- c(sum(group_var == levs[i] & !is.na(x)), sum(group_var == levs[j] & !is.na(x)))
      pval_num[k] <- prop.test(xi, ni, alternative = alternative)$p.value
      names(pval_num)[k] <- paste(levs[i], "vs.", levs[j])
    }
  }
  
  if (adjust) {
    pval_num <- p.adjust(pval_num, method = "bonferroni")
  }
  
  for (i in seq_along(pval_num)) {
    if (is.na(pval_num[i])) {
      pvals[i] <- "/"
    } else if (pval_num[i] < 10 ^ (-digits)) {
      pvals[i] <- paste0("<", 10 ^ (-digits))
    } else {
      pvals[i] <- format(round(pval_num[i], digits), nsmall = digits)
    }
  }
  names(pvals) <- names(pval_num)
  return(pvals)
}


## p-values of Pearson Chi-squared test for group comparison
get_pval_fisher <- function(x0, group_var0, levs, alternative = "two.sided",
                            adjust = FALSE, digits = 3, subset = rep(T, length(x0))){
  x <- x0[subset]
  group_var <- group_var0[subset]
  x <- as.numeric(x)
  pvals <- pval_num <- rep(NA, choose(length(levs), 2))
  
  k <- 0
  for (i in 1:(length(levs) - 1)) {
    for (j in (i + 1):length(levs)) {
      k <- k + 1
      
      cont_tab <- table(x[group_var %in% c(levs[i], levs[j])], 
                        group_var[group_var %in% c(levs[i], levs[j])])[, c(levs[i], levs[j])]
      xi <- x[group_var %in% c(levs[i], levs[j])]
      yi <- group_var[group_var %in% c(levs[i], levs[j])]
      
      result_ij <- try(fisher.test(x = xi, y = yi, alternative = alternative),
                       silent = T)
      if (class(result_ij) == "try-error") {
        pval_num[k] <- NA
      } else {
        pval_num[k] <- result_ij$p.value
      }
      names(pval_num)[k] <- paste(levs[i], "vs.", levs[j])
    }
  }
  
  if (adjust) {
    pval_num <- p.adjust(pval_num, method = "bonferroni")
  }
  
  for (i in seq_along(pval_num)) {
    if (is.na(pval_num[i])) {
      pvals[i] <- "/"
    } else if (pval_num[i] < 10 ^ (-digits)) {
      pvals[i] <- paste0("<", 10 ^ (-digits))
    } else {
      pvals[i] <- format(round(pval_num[i], digits), nsmall = digits)
    }
  }
  names(pvals) <- names(pval_num)
  return(pvals)
}

## p-values of Fisher's exact test for the association of two vectors
get_pval_fisher_2v <- function(x, y, alternative = "two.sided", digits = 3, ...){
  x0 <- x[complete.cases(cbind(x, y))]
  y0 <- y[complete.cases(cbind(x, y))]
  pval <- fisher.test(x = x0, y = y0, alternative = alternative, ...)$p.value
  
  
  
  if (is.na(pval)) {
    res <- "/"
  } else if (pval < 10 ^ (-digits)) {
    res <- paste0("<", 10 ^ (-digits))
  } else {
    res <- format(round(pval, digits), nsmall = digits)
  }
  
  res
}
## p-values of two-sample t-test for group comparison
get_pval_ttest <- function(x, group_var, levs, alternative = "two.sided",
                           adjust = FALSE, digits = 3){
  x <- as.numeric(x)
  pvals <- pval_num <- rep(NA, choose(length(levs), 2))
  
  k <- 0
  for (i in 1:(length(levs) - 1)) {
    for (j in (i + 1):length(levs)) {
      k <- k + 1
      x1i <- x[group_var == levs[i]]
      x2i <- x[group_var == levs[j]]
      pval_num[k] <- t.test(x1i, x2i, alternative = alternative)$p.value
      names(pval_num)[k] <- paste(levs[i], "vs.", levs[j])
    }
  }
  if (adjust) {
    pval_num <- p.adjust(pval_num, method = "bonferroni")
  }
  
  for (i in seq_along(pval_num)) {
    if (is.na(pval_num[i])) {
      pvals[i] <- "/"
    } else if (pval_num[i] < 10 ^ (-digits)) {
      pvals[i] <- paste0("<", 10 ^ (-digits))
    } else {
      pvals[i] <- format(round(pval_num[i], digits), nsmall = digits)
    }
  }
  names(pvals) <- names(pval_num)
  return(pvals)
}

## p-values of two-sample t-test for group comparison
get_pval_wilcox <- function(x, group_var, levs, alternative = "two.sided",
                            adjust = FALSE, digits = 3){
  x <- as.numeric(x)
  pvals <- pval_num <- rep(NA, choose(length(levs), 2))
  
  k <- 0
  for (i in 1:(length(levs) - 1)) {
    for (j in (i + 1):length(levs)) {
      k <- k + 1
      x1i <- x[group_var == levs[i]]
      x2i <- x[group_var == levs[j]]
      pval_num[k] <- wilcox.test(x1i, x2i, alternative = alternative)$p.value
      names(pval_num)[k] <- paste(levs[i], "vs.", levs[j])
    }
  }
  if (adjust) {
    pval_num <- p.adjust(pval_num, method = "bonferroni")
  }
  
  for (i in seq_along(pval_num)) {
    if (is.na(pval_num[i])) {
      pvals[i] <- "/"
    } else if (pval_num[i] < 10 ^ (-digits)) {
      pvals[i] <- paste0("<", 10 ^ (-digits))
    } else {
      pvals[i] <- format(round(pval_num[i], digits), nsmall = digits)
    }
  }
  names(pvals) <- names(pval_num)
  return(pvals)
}


## p-values of two-sample t-test for two-group comparison
get_pval_ttest_2g <- function(x, group_var, alternative = "two.sided", digits = 3){
  x0 <- as.numeric(x)[complete.cases(cbind(x, group_var))]
  group_var0 <- group_var[complete.cases(cbind(x, group_var))]
  
  levs <- unique(group_var0)
  
  x1 <- x0[group_var0 == levs[1]]
  x2 <- x0[group_var0 == levs[2]]
  pval_num <- t.test(x1, x2, alternative = alternative)$p.value
  
  if (is.na(pval_num)) {
    pval <- "/"
  } else if (pval_num < 10 ^ (-digits)) {
    pval <- paste0("<", 10 ^ (-digits))
  } else {
    pval <- format(round(pval_num, digits), nsmall = digits)
  }
  return(pval)
}

## p-values of two-sample wilcoxin-test for two-group comparison
get_pval_wilcox_2g <- function(x, group_var, alternative = "two.sided", digits = 3){
  x0 <- as.numeric(x)[complete.cases(cbind(x, group_var))]
  group_var0 <- group_var[complete.cases(cbind(x, group_var))]
  
  levs <- unique(group_var0)
  
  x1 <- x0[group_var0 == levs[1]]
  x2 <- x0[group_var0 == levs[2]]
  pval_num <- wilcox.test(x1, x2, alternative = alternative)$p.value
  
  if (is.na(pval_num)) {
    pval <- "/"
  } else if (pval_num < 10 ^ (-digits)) {
    pval <- paste0("<", 10 ^ (-digits))
  } else {
    pval <- format(round(pval_num, digits), nsmall = digits)
  }
  return(pval)
}


## get the summary statistics for a regression model
get_summary_reg <- function(reg, digits = 3, pval_digits = 3,
                            transformation = NULL, sandwich = T, ...) {
  if (sandwich) {
    reg_tab <- data.frame(est = coef(reg)[!is.na(coef(reg))],
                          ci_lb = lmtest::coefci(reg, vcov. = sandwich::vcovHC(reg, "HC0"))[, 1],
                          ci_ub = lmtest::coefci(reg, vcov. = sandwich::vcovHC(reg, "HC0"))[, 2],
                          pval = lmtest::coeftest(reg, vcov. = sandwich::vcovHC(reg, "HC0"))[, "Pr(>|z|)"])
  } else {
    reg_tab <- data.frame(est = coef(reg)[!is.na(coef(reg))],
                          ci_lb = confint.default(reg)[!is.na(coef(reg)), 1],
                          ci_ub = confint.default(reg)[!is.na(coef(reg)), 2],
                          pval = summary(reg)$coefficients[, "Pr(>|z|)"])
  }
  
  if (!is.null(transformation)) {
    reg_tab <- reg_tab %>%
      mutate(est = transformation(est),
             ci_lb = transformation(ci_lb),
             ci_ub = transformation(ci_ub))
  }
  
  result <- with(reg_tab, data.frame(
    ESTIMATE =  format(round(est, digits), digits),
    CI = paste0(format(round(ci_lb, digits), digits), ", ",
                format(round(ci_ub, digits), digits)),
    PVAL = format(round(pval, pval_digits), pval_digits)
  ))
  
  rownames(result) <- names(coef(reg)[!is.na(coef(reg))])
  return(result)
}

####################################################
## Get the summary of cumulative incidence function
###################################################
get_cuminc_summary <- function(rslt, times, confint.level = 0.95) {
  rslt_time <- timepoints(rslt, times)
  est_names <- names(rslt)[1:(length(rslt) - 1)]
  x <- rslt_time$est
  s <- sqrt(rslt_time$var)
  z <- qnorm(1 - (1 - confint.level) / 2)
  
  ci_lb <- x ^ exp(- z * s / (x * log(x)))
  ci_ub <- x ^ exp(z * s / (x * log(x)))
  
  
  return(list(est = x, lb = ci_lb, ub = ci_ub))
}

get_survfit_summary <- function(rslt, times, digits = 3) {
  ss <- summary(rslt, times = times)
  out_tab <- data.frame(group = rep(rownames(ss$table), each = length(times)),
                        times = rep(times, length(rownames(ss$table))),
                        cuminc = format(1 - round(ss$surv, digits), nsmall = digits),
                        cuminc_lb = format(1 - round(ss$upper, digits), nsmall = digits),
                        cuminc_ub = format(1 - round(ss$lower, digits), nsmall = digits))
  out_tab$summary <- with(out_tab, sprintf("%s (%s, %s)", cuminc, cuminc_lb, cuminc_ub))
  
  return(out_tab)
}

## calculate the year difference between startdate and stopdate
yrdif <- function(startdate, stopdate, type = 'actual'){
  ndate <- length(startdate)
  diff <- rep(NA, ndate)
  for (i in 1:ndate) {
    if (type == 'actual'){
      ys <- year(startdate[i])
      ye <- year(stopdate[i])
      
      if (is.na(ys) | is.na(ye)) {
        diff[i] <- NA
      } else {
        yslast <- as.Date(paste0(year(startdate[i]) + 1,'-01-', '01'))
        yelast <- as.Date(paste0(year(stopdate[i]) + 1,'-01-', '01'))
        
        if (ys == ye){
          
          if (leap_year(ys)) {
            diff[i] <- (stopdate[i] - startdate[i]) / 366
          }
          
          if (!leap_year(ys)){
            diff[i] <- (stopdate[i] - startdate[i]) / 365
          }
        }
        
        if (ys != ye){
          
          denom <- rep(1, length(ys:ye))
          
          if (leap_year(ys)) {denom[1] <- (yslast - startdate[i]) / 366}
          if (!leap_year(ys)) {denom[1] <- (yslast - startdate[i]) / 365}
          
          if (leap_year(ye)) {denom[length(ys:ye)] <- abs(1 - (yelast - stopdate[i]) / 366)}
          if (!leap_year(ye)) {denom[length(ys:ye)] <- abs(1 - (yelast - stopdate[i]) / 365)}
          
          diff[i] <- sum(denom)
          
        }
      }
    }
    
    if (type == '365'){
      diff[i] <- as.numeric(as.Date(stopdate[i]) - as.Date(startdate[i])) / 365
    }
  }
  return(diff)
}

## get race/ethnicity
get_race_ethnicity <- function(racegrp, hispanic) {
  race_ethnicity <- rep(NA, length(racegrp))
  
  for (i in 1:length(racegrp)) {
    racegrp_i <- racegrp[i]
    hispanic_i <- hispanic[i]
    if (racegrp_i == "White" & hispanic_i == "Non Hispanic/Latino") {
      race_ethnicity[i] <- "Non-Hispanic White"
    } else if (racegrp_i == "Black" & hispanic_i == "Non Hispanic/Latino") {
      race_ethnicity[i] <- "Non-Hispanic Black"
    } else if (hispanic_i == "Hispanic/Latino") {
      race_ethnicity[i] <- "Hispanic"
    } else if (racegrp_i == "Other") {
      race_ethnicity[i] <- "Other"
    } 
  }
  return(race_ethnicity)
}


## Using Slaughter's skinfold equation to get body fat percent for
## adolescents
get_pfat_from_skinfold <- function(sftric1, sftric2, sfscap1, sfscap2,
                                   age_sf, gender, racegrp) {
  nn <-  length(sftric1)
  pfat <- rep(NA, nn)
  
  sftric <- (sftric1 + sftric2) / 2
  sfscap <- (sfscap1 + sfscap2) / 2
  
  for (i in 1:nn) {
    if (sftric[i] + sfscap[i] > 35) {
      if (gender[i] == "Male") {
        pfat[i] <- 0.783 * (sftric[i] + sfscap[i]) + 1.6
      } else {
        pfat[i] <- 0.546 * (sftric[i] + sfscap[i]) + 9.7
      }
    } else {
      if (age_sf[i] < 9 & racegrp[i] == "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 1.7
      }
      
      if (age_sf[i] < 9 & racegrp[i] != "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 3.2
      }
      
      if (age_sf[i] >= 9 & age_sf[i] < 14 & 
          racegrp[i] == "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 3.4
      }
      
      if (age_sf[i] >= 9 & age_sf[i] < 14 & 
          racegrp[i] != "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 5.2
      }
      
      if (age_sf[i] >= 14 & age_sf[i] < 18 & 
          racegrp[i] == "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 5.5
      }
      
      if (age_sf[i] >= 14 & age_sf[i] < 18 & 
          racegrp[i] != "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 6.8
      }
      
      if (gender[i] == "Female") {
        pfat[i] <- 1.33 * (sftric[i] + sfscap[i]) - 
          0.013 * (sftric[i] + sfscap[i]) ^ 2 - 2.5
      }
    }
  }
  return(pfat)
}



## Using Slaughter's skinfold equation to get body fat percent for
## adolescents
get_pfat_from_skinfold_adolescent <- function(sftric1, sftric2, sfscap1, sfscap2,
                                              age_sf, gender, racegrp) {
  nn <-  length(sftric1)
  pfat <- rep(NA, nn)
  
  sftric <- (sftric1 + sftric2) / 2
  sfscap <- (sfscap1 + sfscap2) / 2
  
  for (i in 1:nn) {
    if (sftric[i] + sfscap[i] > 35) {
      if (gender[i] == "Male") {
        pfat[i] <- 0.783 * (sftric[i] + sfscap[i]) + 1.6
      } else {
        pfat[i] <- 0.546 * (sftric[i] + sfscap[i]) + 9.7
      }
    } else {
      if (age_sf[i] < 9 & racegrp[i] == "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 1.7
      }
      
      if (age_sf[i] < 9 & racegrp[i] != "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 3.2
      }
      
      if (age_sf[i] >= 9 & age_sf[i] < 14 & 
          racegrp[i] == "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 3.4
      }
      
      if (age_sf[i] >= 9 & age_sf[i] < 14 & 
          racegrp[i] != "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 5.2
      }
      
      if (age_sf[i] >= 14 & age_sf[i] < 18 & 
          racegrp[i] == "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 5.5
      }
      
      if (age_sf[i] >= 14 & age_sf[i] < 18 & 
          racegrp[i] != "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 6.8
      }
      
      if (gender[i] == "Female") {
        pfat[i] <- 1.33 * (sftric[i] + sfscap[i]) - 
          0.013 * (sftric[i] + sfscap[i]) ^ 2 - 2.5
      }
    }
  }
  return(pfat)
}



## Using Slaughter's skinfold equation to get body fat percent for
## adults

get_pfat_from_skinfold_adult <- function(sfchest, sfabdom, sfthigh_m,
                                         sftric, sfiliac, sfthigh_f,
                                         age_sf, gender, racegrp) {
  nn <-  length(sftric)
  
  sum_sf <- ifelse(gender == "Male", 
                   sfchest + sfabdom + sfthigh_m,
                   sftric + sfiliac + sfthigh_f)
  
  body_density <- ifelse(gender == "Male",
                         1.10938 - 0.0008267 * sum_sf + 0.0000016 * sum_sf ^ 2 - 0.0002574 * age_sf,
                         1.099421 - 0.0009929 * sum_sf + 0.0000023 * sum_sf ^ 2 - 0.0001392 * age_sf)
  
  
  pfat <- case_when(gender == "Male" & racegrp == "Black" ~ 4.86 / body_density - 4.39,
                    gender == "Female" & racegrp == "Black" ~ 4.85 / body_density - 4.39,
                    gender == "Male" & racegrp != "Black" ~ 4.95 / body_density - 4.5,
                    gender == "Female" & racegrp != "Black" ~ 5.01 / body_density - 4.57)
  
  return(pfat)
}

## among a list of visits, match an mrn and a date with a specific visit
## method can be one of "closest", "closest_prior", "closest_after",
## "within", "within_prior", "within_after"
## action_missing_date_a can be one of "NA", "earliest" or "latest"
match_visit <- function(id = NULL, id_a = NULL, id_b = NULL,
                        data_a, data_b,
                        date_a, date_b,
                        date_name = date_a,
                        var_b = NULL,
                        method = "closest", 
                        action_missing_date_a = "earliest",
                        within = 180){
  
  
  if (!is.null(id)) {
    id_a <- id_b <- id
  }
  
  
  if (is.null(var_b)) {
    data_b1 <- data_b
  } else {
    data_b1 <- data_b %>% 
      filter(c(!is.na(data_b[, var_b]))) 
  }
  
  data_b_names <- setdiff(names(data_b1), c(id_b, date_b))
  
  data_out <- data.frame(person_name = data_a[, id_a],
                         person_visit_date = data_a[, date_a]) 
  
  
  for (i in 1:nrow(data_a)) {
    id_ai <- unlist(data_a[i, id_a])
    date_ai <- unlist(data_a[i, date_a])
    data_bi <- data_b1[data_b1[, id_b] == id_ai, ]
    
    if (is.na(date_ai)) {
      if (nrow(data_bi) > 0) {
        if (action_missing_date_a == "earliest") {
          id_select <- which.min(data_bi[, date_b])
          data_out[i, data_b_names] <- data_bi[id_select, data_b_names]
        } else if (action_missing_date_a == "lastest") {
          id_select <- which.max(data_bi[, date_b])
          data_out[i, data_b_names] <- data_bi[id_select, data_b_names]
        } else if (action_missing_date_a == "na") {
          data_out[i, data_b_names] <- NA
        }
      }
    } else if (nrow(data_bi) > 0) {
      diff_date <- as.Date(date_ai) - 
        as.Date(unlist(data_bi[, date_b]))
      if (method == "closest") {
        id_select <- which.min(abs(diff_date))[1]
        data_out[i, data_b_names] <- data_bi[id_select, data_b_names]
      }
      
      if (method == "within") {
        elig_dates_n <- sum(abs(diff_date) <= within)
        if (elig_dates_n >= 1) {
          id_select <- which.min(abs(diff_date))[1]
          data_out[i, data_b_names] <- data_bi[id_select, data_b_names]
        }
      }
    }
  }
  return(data_out)
}



###################################################
## Derived Variables
###################################################
## get participants' SJLIFE status
get_status <- function(x) {
  y <- rep(NA, length(x))
  
  for (i in seq_along(x)) {
    if (x[i] == 3) {
      y[i] <- 1 ## Eligible survivor participants
    } else if (x[i] %in% c(8, 22)) {
      y[i] <- 2 ## SJLIFE ineligible
    } else if (x[i] %in% c(15, 18)) {
      y[i] <- 3 ## Not recruited
    } else if (x[i] %in% c(1, 2, 13, 21, 11)) {
      y[i] <- 4 ## Not visited yet
    } else if (x[i] %in% c(17, 20)) {
      y[i] <- 5 ## died prior to visit
    } else if (x[i] %in% c(4, 5, 24)) {
      y[i] <- 6 ## survey only
    } else if (x[i] %in% c(9, 12, 19, 23)) {
      y[i] <- 7 ## Refused
    } else if (x[i] %in% c(6, 7, 10, 14, 99)) {
      y[i] <- 8 ## Nonresponse
    }
  }
  
  return(y)
}



get_hsct <- function(sjlife_path = "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/",
                     mrn, dates, datesname = "date") {
  transplant <- read_sas(paste0(sjlife_path, "Clinical data/transplant.sas7bdat"))%>%
    rename_with(tolower) %>%
    arrange(mrn)
  
  data_out <- data.frame(mrn = mrn, dates = dates)
  names(data_out)[2] <- datesname
  data_out$hsc_2g <- NA
  
  hsc_tran <- transplant %>%
    dplyr::select(mrn, tpdt, tptype)
  
  for (i in 1:nrow(data_out)) {
    hsc_tran_i <- filter(hsc_tran,
                         mrn == data_out$mrn[i] &
                           tpdt <= data_out[, datesname][i])
    if (nrow(hsc_tran_i) > 1) {
      data_out$hsc_2g[i] <- 1
    } else {
      data_out$hsc_2g[i] <- 0
    }
  }
  
  return(data_out)
}


## adolescent physical activity from health habits survey
get_adolescent_pa <- function(sjlife_path = "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/") {
  adolescent_healthhabits <- read_sas(paste0(sjlife_path, "Survey data/adolescent_healthhabits.sas7bdat")) %>%
    rename_with(tolower)
  
  ## compute MET for each adolescent with MRN in the study population
  pa_adoles <- adolescent_healthhabits %>%
    dplyr::select(mrn, datecomp, percomp,vpa10, vpadays, vpamin,
                  mpa10, mpadays, mpamin, lpa10, lpadays, lpamin,
                  activity, play60, pa20days, patone,
                  screen, relation, survey) %>%
    mutate(wvpa = NA, wmpa = NA, mvpawk = NA, wmspa = NA, mets = NA)
  
  
  for (i in 1:nrow(pa_adoles)) {
    if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 1 & is.na(pa_adoles$vpadays[i])) {
      pa_adoles$vpadays[i] <- 1
    }
    
    if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 1 & is.na(pa_adoles$vpamin[i])) {
      pa_adoles$vpamin[i] <- 10
    }
    
    if (!is.na(pa_adoles$vpamin[i]) & pa_adoles$vpamin[i] > 360) {
      pa_adoles$vpamin[i] <- 360
    }
    
    if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 1) {
      pa_adoles$wvpa[i] = pa_adoles$vpadays[i] * pa_adoles$vpamin[i]
    }
    
    if (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2) {
      pa_adoles$wvpa[i] = 0
    }
    
    if (is.na(pa_adoles$wvpa[i])) {
      if ((is.na(pa_adoles$vpa10[i]) | (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2)) &
          (!is.na(pa_adoles$activity[i]) & pa_adoles$activity[i] == 2)) {
        pa_adoles$wvpa[i] <- 0
      } 
      
      if ((is.na(pa_adoles$vpa10[i])|(!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2)) &
          !(is.na(pa_adoles$play60[i])|(!is.na(pa_adoles$play60[i]) & pa_adoles$play60[i] == 1))) {
        pa_adoles$wvpa[i] <- 0
      }
      
      if ((is.na(pa_adoles$vpa10[i]) | (!is.na(pa_adoles$vpa10[i]) & pa_adoles$vpa10[i] == 2)) &
          !(is.na(pa_adoles$pa20days[i]) | (!is.na(pa_adoles$pa20days[i]) & pa_adoles$pa20days[i] == 1))) {
        pa_adoles$wvpa[i] <- (pa_adoles$pa20days[i] - 1) * 20
      }
    }
    
    if (is.na(pa_adoles$mpa10[i]) & 
        ((!is.na(pa_adoles$mpadays[i]))|(!is.na(pa_adoles$mpamin[i])))) {
      pa_adoles$mpa10[i] <- 1
    }
    
    if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 1 & is.na(pa_adoles$mpadays[i])) {
      pa_adoles$mpadays[i] <- 1
    }
    
    if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 1 & is.na(pa_adoles$mpamin[i])) {
      pa_adoles$mpamin[i] <- 10
    }
    
    if (!is.na(pa_adoles$mpamin[i]) & pa_adoles$mpamin[i] > 360) pa_adoles$mpamin[i] <- 360
    if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 1) {
      pa_adoles$wmpa[i] <- 
        pa_adoles$mpadays[i] * pa_adoles$mpamin[i]
    }
    if (!is.na(pa_adoles$mpa10[i]) & pa_adoles$mpa10[i] == 2) {
      pa_adoles$wmpa[i] <-0
    }
    
    if (is.na(pa_adoles$wmpa[i])) {
      if ((is.na(pa_adoles$mpa10[i]) | pa_adoles$mpa10[i] == 2) &
          (!is.na(pa_adoles$activity[i]) & pa_adoles$activity[i] == 2)) {
        pa_adoles$wmpa[i] <- 0
      } 
      
      if ((is.na(pa_adoles$mpa10[i]) | pa_adoles$mpa10[i] == 2) &
          !(is.na(pa_adoles$play60[i]) | pa_adoles$play60[i] == 1)) {
        pa_adoles$wmpa[i] <- (pa_adoles$play60[i] - 1) * 60
      } 
      
      if (!is.na(pa_adoles$wvpa[i])) pa_adoles$wmpa[i] <- 0
    }
    
    
    if (is.na(pa_adoles$wmpa[i]) & is.na(pa_adoles$wvpa[i])) {
      if (!is.na(pa_adoles$play60[i]) & pa_adoles$play60[i] == 1) {
        pa_adoles$wvpa[i] <- pa_adoles$wmpa[i] <- 0
      }
    }
    
    if (!is.na(pa_adoles$patone[i]) & pa_adoles$patone[i] > 1) {
      pa_adoles$wmspa[i] <- (pa_adoles$patone[i] - 1) * 20
    } else {
      pa_adoles$wmspa[i] <- 0
    }
    
    pa_adoles$mvpawk[i] <- pa_adoles$wmpa[i] + 2 * pa_adoles$wvpa[i] +
      2 * pa_adoles$wmspa[i]
    
    if (is.na(pa_adoles$wvpa[i])) pa_adoles$mvpawk[i] <- pa_adoles$wmpa[i]
    
    if (!is.na(pa_adoles$mvpawk[i]) & pa_adoles$mvpawk[i] > 2520) pa_adoles$mvpawk[i] <- 2520
    
    pa_adoles$mets[i] <- 3 * pa_adoles$mvpawk[i]
    
  }
  
  pa_adoles_final <- dplyr::select(pa_adoles, mrn, datecomp, relation, percomp, 
                                   survey, screen, wvpa, wmpa, mvpawk, mets)
  return(pa_adoles_final)
}

get_smoking_status <- function(evsm, cigmo, smnow) {
  smoker_str <- case_when(is.na(evsm) & is.na(cigmo) & is.na(smnow) ~ NA,
                          cigmo == 2 & evsm == 1 ~ "Former smoke",
                          cigmo == 1 ~ "Now smoke",
                          cigmo == 2 & evsm == 2 ~ "Never smoke",
                          .default = NA)
  smoker <- factor(smoker_str, 
                   levels = c("Now smoke", "Former smoke", "Never smoke"),
                   labels = c("Now smoke", "Former smoke", "Never smoke"))
  
  return(smoker)
}



## adult pa from health habits survey
get_adult_pa <- function(sjlife_path = "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/") {
  adult_healthhabits <- read_sas(paste0(sjlife_path, "Survey data/adult_healthhabits.sas7bdat")) %>%
    rename_with(tolower)
  
  pa_adults <- adult_healthhabits %>%
    dplyr::select(mrn, datecomp, relation, nopa,
                  vpa10, vpadays, vpamin, mpa10, mpadays, mpamin,
                  lpa10, lpadays, lpamin, ltpaw, wtlt, yoga, pa20) %>%
    mutate(wvpa = NA, wmpa = NA, mvpawk = NA,
           mets = NA, wlpa = NA)
  
  
  for (i in 1:nrow(pa_adults)) {
    if (pa_adults$vpa10[i] %in% 1 & is.na(pa_adults$vpadays[i])) {
      pa_adults$vpadays[i] <- 1
    }
    
    if (pa_adults$vpa10[i] %in% 1 & is.na(pa_adults$vpamin[i])) {
      pa_adults$vpamin[i] <- 10
    }
    
    if (!is.na(pa_adults$vpamin[i]) & pa_adults$vpamin[i] > 360) {
      pa_adults$vpamin[i] <- 360
    }
    
    if (pa_adults$vpa10[i] %in% 1) {
      pa_adults$wvpa[i] <- pa_adults$vpadays[i] * pa_adults$vpamin[i]
    }
    
    if (pa_adults$vpa10[i] %in% 2) {
      pa_adults$wvpa[i] <- 0
    }
    
    if (is.na(pa_adults$wvpa[i])) {
      if ((is.na(pa_adults$vpa10[i]) | (pa_adults$vpa10[i] %in% 2)) &
          (pa_adults$nopa[i] %in% 1) &
          !(is.na(pa_adults$pa20[i]) | (pa_adults$pa20[i] %in% 1))) {
        pa_adults$wvpa[i] <- 0
      } 
      
      if ((is.na(pa_adults$vpa10[i])|(pa_adults$vpa10[i] %in% 2)) &
          !(is.na(pa_adults$pa20[i])|(pa_adults$pa20[i] %in% 1))) {
        pa_adults$wvpa[i] <- (pa_adults$pa20[i] - 1) * 20
      }
      
      if ((is.na(pa_adults$vpa10[i]) | (pa_adults$vpa10[i] %in% 2)) &
          (pa_adults$nopa[i] %in% 2) &
          (is.na(pa_adults$pa20[i]) | (pa_adults$pa20[i] %in% 1))) {
        pa_adults$wvpa[i] <- 0
      }
    }
    
    if (is.na(pa_adults$mpa10[i]) & 
        ((!is.na(pa_adults$mpadays[i]))|(!is.na(pa_adults$mpamin[i])))) {
      pa_adults$mpa10[i] <- 1
    }
    
    if (pa_adults$mpa10[i] %in% 1 & is.na(pa_adults$mpadays[i])) {
      pa_adults$mpadays[i] <- 1
    }
    
    if (pa_adults$mpa10[i] %in% 1 & is.na(pa_adults$mpamin[i])) {
      pa_adults$mpamin[i] <- 10
    }
    
    if (!is.na(pa_adults$mpamin[i]) & pa_adults$mpamin[i] > 360) pa_adults$mpamin[i] <- 360
    if (pa_adults$mpa10[i] %in% 1) {
      pa_adults$wmpa[i] <- 
        pa_adults$mpadays[i] * pa_adults$mpamin[i]
    }
    if (pa_adults$mpa10[i] %in% 2) {
      pa_adults$wmpa[i] <- 0
    }
    
    if (is.na(pa_adults$wmpa[i])) {
      if ((is.na(pa_adults$mpa10[i]) | 
           (pa_adults$mpa10[i] %in% 2)) &
          (pa_adults$nopa[i] %in% 1)) {
        pa_adults$wmpa[i] <- 0
      } 
      
      if (!is.na(pa_adults$wvpa[i])) {
        pa_adults$wmpa[i] <- 0
      } 
      
    }
    
    pa_adults$mvpawk[i] <- pa_adults$wmpa[i] + 2 * pa_adults$wvpa[i] 
    if (!is.na(pa_adults$mvpawk[i]) & pa_adults$mvpawk[i] > 2520) pa_adults$mvpawk[i] <- 2520
    
    pa_adults$mets[i] <- 3 * pa_adults$mvpawk[i]
    
    
    if (is.na(pa_adults$lpa10[i]) & 
        ((!is.na(pa_adults$lpadays[i]))|(!is.na(pa_adults$lpamin[i])))) {
      pa_adults$lpa10[i] <- 1
    }
    
    if (pa_adults$lpa10[i] %in% 1 & is.na(pa_adults$lpadays[i])) {
      pa_adults$lpadays[i] <- 1
    }
    
    if (pa_adults$lpa10[i] %in% 1 & is.na(pa_adults$lpamin[i])) {
      pa_adults$lpamin[i] <- 10
    }
    
    if (!is.na(pa_adults$lpamin[i]) & pa_adults$lpamin[i] > 360) pa_adults$lpamin[i] <- 360
    if (pa_adults$lpa10[i] %in% 1) {
      pa_adults$wlpa[i] <- 
        pa_adults$lpadays[i] * pa_adults$lpamin[i]
    }
    if (pa_adults$lpa10[i] %in% 2) {
      pa_adults$wlpa[i] <- 0
    }
    
    if (is.na(pa_adults$wlpa[i])) {
      if ((is.na(pa_adults$lpa10[i]) | 
           (pa_adults$lpa10[i] %in% 2)) &
          (pa_adults$nopa[i] %in% 1)) {
        pa_adults$wlpa[i] <- 0
      } 
      
      if (!is.na(pa_adults$wmpa[i])) {
        pa_adults$wlpa[i] <- 0
      } 
      
    }
  }
  
  pa_adults_final <- dplyr::select(pa_adults, mrn, datecomp, relation, 
                                   wmpa, wvpa, mvpawk, mets, wlpa, 
                                   mvpawk)
  return(pa_adults_final)
}

get_offspring <- function(stlife_path = "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/",
                          mrn_list) {
  adult_home <- read_sas(paste0(sjlife_path, "Survey data/adult_home.sas7bdat")) %>%
    rename_with(tolower) %>%
    filter(mrn %in% mrn_list)
  
  everpreg <- adult_home$everpreg
  n_preg <- adult_home$n_preg
  n_offspring <- rep(NA, length(everpreg))
  pregout1 <- adult_home$pregout1
  pregout2 <- adult_home$pregout2
  pregout3 <- adult_home$pregout3
  pregout4 <- adult_home$pregout4
  pregout5 <- adult_home$pregout5
  pregout6 <- adult_home$pregout6
  pregout7 <- adult_home$pregout7
  pregout8 <- adult_home$pregout8
  
  for (i in 1:length(n_offspring)) {
    if (everpreg[i] %in% 2) {
      n_preg[i] <- 0
      n_offspring[i] <- 0
    } else if (n_preg[i] %in% 1) {
      n_offspring[i] <- sum(pregout1[i] == 1)
    } else if (n_preg[i] %in% 2) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) 
    } else if (n_preg[i] %in% 3) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) +
        as.numeric(pregout3[i] == 1)
    } else if (n_preg[i] %in% 4) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) +
        as.numeric(pregout3[i] == 1) +
        as.numeric(pregout4[i] == 1)
    } else if (n_preg[i] %in% 5) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) +
        as.numeric(pregout3[i] == 1) +
        as.numeric(pregout4[i] == 1) + 
        as.numeric(pregout5[i] == 1)
    } else if (n_preg[i] %in% 6) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) +
        as.numeric(pregout3[i] == 1) +
        as.numeric(pregout4[i] == 1) + 
        as.numeric(pregout5[i] == 1) +
        as.numeric(pregout6[i] == 1)
    } else if (n_preg[i] %in% 7) {
      n_offspring[i] <- as.numeric(pregout1[i] == 1) + 
        as.numeric(pregout2[i] == 1) +
        as.numeric(pregout3[i] == 1) +
        as.numeric(pregout4[i] == 1) + 
        as.numeric(pregout5[i] == 1) +
        as.numeric(pregout6[i] == 1) +
        as.numeric(pregout7[i] == 1)
    } else if (is.na(n_preg[i]) | n_preg[i] >= 8) {
      n_offspring[i] <- sum(pregout1[i] == 1, pregout2[i] == 1,
                            pregout3[i] == 1, pregout4[i] == 1,
                            pregout5[i] == 1, pregout6[i] == 1,
                            pregout7[i] == 1, pregout8[i] == 1,
                            na.rm = T)
    }
  }
  
  return(data.frame(mrn = adult_home$mrn,
                    datecomp = adult_home$datecomp,
                    relation = adult_home$relation,
                    everpreg = everpreg,
                    n_preg = n_preg,
                    n_offspring = n_offspring))
}
## method: closest, most recent, or past max

get_ctcae_grade <- function(stlife_path = "Z:/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/",
                            mrn, dates, conditions, 
                            method = "history", within = 180) {
  ctcaegrades <- read_sas(paste0(sjlife_path, "Event data/ctcaegrades.sas7bdat")) %>% 
    rename_with(tolower) %>%
    arrange(mrn) %>%
    filter(!grade %in% c(NA, -9))
  
  dat1 <- data.frame(mrn = mrn, dates = dates) %>%
    mutate(mrn = as.numeric(mrn))
  
  dat2 <- dat1
  
  for (cc in 1:length(conditions)) {
    sub_grades <- ctcaegrades %>% 
      filter(condition == conditions[cc]) %>%
      dplyr::select(mrn, gradedt, grade) %>%
      mutate(mrn = as.numeric(mrn))
    
    grade_dat <- dat1 %>% 
      inner_join(sub_grades, by = "mrn", relationship = "many-to-many") %>%
      filter(gradedt <= dates + within) %>%
      arrange(mrn) %>%
      group_by(mrn, dates) %>%
      summarize(grade = max(grade)) %>%
      ungroup() %>%
      select(mrn, dates, grade) 
    
    names(grade_dat)[3] <- conditions[cc] %>%
      tolower() %>%
      stringr::str_replace_all(" ", "_")
    
    dat2 <- left_join(dat2, grade_dat, by = c("mrn", "dates"))
  }
  
  return(dat2)
}

#########################################################################
## Get body fat percentage from skinfold data using Slaughter's equation
########################################################################
## Using Slaughter's skinfold equation to get body fat percent
get_pfat_from_skinfold <- function(sftric1, sftric2, sfscap1, sfscap2,
                                   age_sf, gender, racegrp) {
  nn <-  length(sftric1)
  pfat <- rep(NA, nn)
  
  sftric <- (sftric1 + sftric2) / 2
  sfscap <- (sfscap1 + sfscap2) / 2
  
  for (i in 1:nn) {
    if (sftric[i] + sfscap[i] > 35) {
      if (gender[i] == "Male") {
        pfat[i] <- 0.783 * (sftric[i] + sfscap[i]) + 1.6
      } else {
        pfat[i] <- 0.546 * (sftric[i] + sfscap[i]) + 9.7
      }
    } else {
      if (age_sf[i] < 9 & racegrp[i] == "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 1.7
      }
      
      if (age_sf[i] < 9 & racegrp[i] != "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 3.2
      }
      
      if (age_sf[i] >= 9 & age_sf[i] < 14 & 
          racegrp[i] == "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 3.4
      }
      
      if (age_sf[i] >= 9 & age_sf[i] < 14 & 
          racegrp[i] != "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 5.2
      }
      
      if (age_sf[i] >= 14 & age_sf[i] < 18 & 
          racegrp[i] == "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 5.5
      }
      
      if (age_sf[i] >= 14 & age_sf[i] < 18 & 
          racegrp[i] != "White" & gender[i] == "Male") {
        pfat[i] <- 1.12 * (sftric[i] + sfscap[i]) - 
          0.008 * (sftric[i] + sfscap[i]) ^ 2 - 6.8
      }
      
      if (gender[i] == "Female") {
        pfat[i] <- 1.33 * (sftric[i] + sfscap[i]) - 
          0.013 * (sftric[i] + sfscap[i]) ^ 2 - 2.5
      }
    }
  }
  return(pfat)
}


#############
## Get SF6D
############
get_sf6d <- function(vigorous, moderate, bathe, phkind, epless, painmuch, painintf,
                     nervous, energy, depres, social2) {
  sf3 <- case_when(vigorous == 1 ~ 2, vigorous == 2 ~ 1, vigorous == 3 ~ 3)
  sf4 <- case_when(moderate == 1 ~ 2, moderate == 2 ~ 1, moderate == 3 ~ 3)
  sf12 <- case_when(bathe == 1 ~ 2, bathe == 2 ~ 1, bathe == 3 ~ 3)
  sf21 <- painmuch; sf22 <- painintf; sf24 <- nervous;
  sf27 <- energy; sf28 <- depres; sf32 <- social2
  
  ## code questions 15 and 18 to have 2 response 
  ## categories instead of 5 
  sf15 <- case_when(phkind < 5 ~ 1,
                    phkind == 5 ~ 2,
                    is.na(phkind) ~ NA)
  
  sf18 <- case_when(epless < 5 ~ 1,
                    epless == 5 ~ 2,
                    is.na(epless) ~ NA)
  
  sf6d <- rep(NA, length(sf3))
  
  
  SFPhys <- case_when(sf3 == 3 & sf4 == 3 & sf12 == 3 ~ 1,
                      sf3 == 2 & sf4 == 3 & sf12 == 3 ~ 2,
                      sf3 == 1 & sf4 == 3 & sf12 == 3 ~ 2,
                      sf4 == 2 & sf12 == 3 ~ 3,
                      sf4 == 1 & sf12 == 3 ~ 4,
                      sf12 == 2 ~ 5,
                      sf12 == 1 ~ 6,
                      .default = NA)
  
  SFRole <- case_when(sf15 == 2 & sf18 == 2 ~ 1,
                      sf15 == 1 & sf18 == 2 ~ 2,
                      sf18 == 1 ~ 3,
                      .default = NA)
  
  SFSocial <- 6 - sf32
  
  SFPain <- case_when(sf21 == 1 & sf22 == 1 ~ 1,
                      sf21 > 1 & sf22 == 1 ~ 2,
                      sf22 > 1 ~ sf22 + 1,
                      .default = NA)
  
  SFMental <- case_when(sf24 < 1 | sf24 > 5 | sf28 < 1 | sf28 > 5 ~ NA,
                        sf24 == 1 | sf28 == 1 ~ 5,
                        sf24 == 2 | sf28 == 2 ~ 4,
                        sf24 == 3 | sf28 == 3 ~ 3,
                        sf24 == 4 | sf28 == 4 ~ 2,
                        sf24 == 5 | sf28 == 5 ~ 1,
                        .default = NA)
  
  
  SFVital <- sf27
  
  # create MOST category if any dimention is at its worst state 
  most <- ifelse(SFPhys %in% 4:6 | SFRole %in% 3:4 |
                   SFSocial %in% 4:5 | SFPain %in% 5:6 |
                   SFMental %in% 4:5 | SFVital %in% 4:5, 1, 0)
  
  # assign decriments based on level from consistent in table 4 and calcuate score
  pf1 <- case_match(SFPhys,
                    1 ~ 0,
                    2 ~ -0.035,
                    3 ~ -0.035,
                    4 ~ -0.044,
                    5 ~ -0.056,
                    6 ~ -0.117,
                    .default = NA)
  
  rl1 <- case_when(SFRole > 1 ~ -0.053,
                   SFRole == 1 ~ 0,
                   .default = NA)
  
  sc1 <- case_match(SFSocial,
                    1 ~ 0,
                    2 ~ -0.057,
                    3 ~ -0.059,
                    4 ~ -0.072,
                    5 ~ -0.087,
                    .default = NA)
  
  pn1 <- case_match(SFPain,
                    1 ~ 0,
                    2 ~ -0.042,
                    3 ~ -0.042,
                    4 ~ -0.065,
                    5 ~ -0.102,
                    6 ~ -0.171,
                    .default = NA)
  
  mh1 <- case_match(SFMental,
                    1 ~ 0,
                    2 ~ -0.042, 
                    3 ~ -0.042,
                    4 ~ -0.100,
                    5 ~ -0.118,
                    .default = NA)
  
  v1 <- case_match(SFVital,
                   1 ~ 0,
                   2 ~ -0.071,
                   3 ~ -0.071,
                   4 ~ -0.071,
                   5 ~ -0.092,
                   .default = NA)
  
  mst1 <- most * -0.061
  
  SF6D <- 1 + pf1 + rl1 + sc1 + pn1 + mh1 + v1 + mst1
  
  return(SF6D)
}




## categorize education

get_education <- function(mrn, grade, gradespe) {
  
  education <- rep(NA, length(grade))
  
  for (i in 1:length(grade)) {
    grade_i <- grade[i]
    gradespe_i <- gradespe[i]
    
    if (grade_i %in% c(6, 7)) {
      education[i] <- 3
    } else if (grade_i %in% c(4, 5)) {
      education[i] <- 2 
    } else if (grade_i %in% 3) {
      education[i] <- 1
    } else if (grade_i %in% 1:2) {
      education[i] <- 0
    } else if (grade_i %in% 9:10) {
      education[i] <- NA
    } else if (grade_i %in% 8 &
               gradespe_i %in% c("11th",
                                 "12th grader -",
                                 "HOME SCHOOLED",
                                 "I have one more GED test to take the I will received my GED",
                                 "I went to Special Education from 1989 to 1996",
                                 "Transitioned through high school special education. Did not earn diploma.",
                                 "currently in 12th grade",
                                 "Still in high school",
                                 "Currently enrolled in 11th grade public high school.",
                                 "I went to Special Education from 1989 to 1996",
                                 "still attending high school, just completed 10th grade",
                                 "12th grader -",
                                 "Still in high school",
                                 "currently in 12th grade",
                                 "Still in High school - 12th grade",
                                 "still in high school",
                                 "11th"
               )) {
      education[i] <- 0
    } else if (grade_i %in% 8 &
               gradespe_i %in% c("11th then got my G.E.D. in- 2014","12th grader -",
                                 "9-12 years (high school with high school diploma)",
                                 "GED",
                                 "GED AGS",
                                 "Graduated through 12th grade but didn't get diploma only a certificate. I'm trying to get my GED.",
                                 "I got my GED",
                                 "I graduate high school in May",
                                 "entering college",
                                 "main streamed Kindergarten through Grade 12.",
                                 "USMC",
                                 "GED",
                                 "9-12 years (high school with high school diploma)",
                                 "GED AGS",
                                 "Graduated through 12th grade but didn't get diploma only a certificate. I'm trying to get my GED.",
                                 "I graduate high school in May",
                                 "I finished high school witha special diploma. I attempted to further my education on two occassions. Once in the automotive field and once in massage therapy. I could not complete either course because of inability to take written exams."
               )) { 
      
      education[i] <- 1
      
      
    } else if (grade_i %in% 8 &
               gradespe_i %in% c("(2 years college) Associates Degree","1 yr of college",
                                 "1/2 years of Schooling in Medical Administration","2 year Degree Associate of Science Business Administration",
                                 "2 years of college","2 yr college","2nd year in college","4th year in college","AA Degree",
                                 "AFTER HIGH SCHOOL ATTENDED COMMUNTY COLLEDGE WITH NO DEGREE. ATTENDED A 4 YEAR ELECTRICAL APPRENTICESHIP. RECEIVED JOURNEYMAN WIREMAN CERTIFICATE",
                                 "Adult College Handicap H.S.D","Adult Education class in Accounting","Assoc. of Arts","Associate Degree in Business Administration from Belhaven",
                                 "Associate degree","Associate graduate Computer Information System","Associate's Degree","Associates","Associates Criminal Justice","Associates Degree",
                                 "Associates Degree + extra hours toward bachelors","Associates Degree 2 years","Associates Degree in Electronics and Computer Technology",
                                 "Associates degree","Associates degree in Business","Associates in Applied Science","Associates in Electronic Engineering ITT Tech.",
                                 "Associates in General Education","Associates in Science Degree","Associates of Arts degree A.A.","Auto-Body tech.","Aviationi Mechanics, Tech School",
                                 "CPI Heavy equipment operating","Certified EKG/ECG Tech.","Community College",
                                 "Completed 2 year acting program at Stella Adler Academy of Acting and Theater Junior @ Grand Canyon University",
                                 "Completed Vocational training.","Cosmetology School","Cosmetology school","Cosmetology school and ASU","Currently a college junior",
                                 "Currently a junior in college","Currently college freshman","Currently enrolled in college","Currently in college",
                                 "Currently working on degree","Degree in Early Childhood Education","Dental Assistant/Trade school PIMA","Dental Assisting school",
                                 "Diesel Tech School","Drafting & Design degree Civil Tech degree 2 yr","Dropped out of high school GED 1981 Associate's Degree 1988",
                                 "Finishing CDA @ ASV Child Development Associate","Flight school, Commercial Pilot","Grad Cosm., Grad Med Asst.",
                                 "Graduated from a vo-tech school for machinist school","Graduated high school 1 1/2 years community college","Graduating from college May 2019",
                                 "I am going into my junior year of college next year","I am in college right now.","I completed 9/ 1/2 yers. of Airconditiong, Heating Repair. In 1979",
                                 "I took full course Cosmetology school. Graduated and licensed and I have been professionally trained to be a cake decorator by Wilton Academy graduates teacher.",
                                 "I went to Barbering School","I went to Concorde Career College. Trained to be a certified nursing assistant.",
                                 "I went to Delta Technical College got a degree in Dental Assistant and Radiology License",
                                 "I went to North Carolina Master's Commission. If I could afford the whole term, I could of became a pastor/youth pastor...",
                                 "I went to college for 1 year and my son got hit by a car so I didn't go back until I was 40. I completed a course to get my CMA and wanted to get my RN. I went back at age 40 but had to drop out because of my health problems.",
                                 "In college now","In college now.","Jr College","Meat processing","Medical Assistant, Phlebotomy Tech, I.V. Tech certified","Medical Asst. Training",
                                 "Paramedic Firefighters Course","Passed CNA test. Continuing on to nursing program.","Several fitness certs, real estate classes & license",
                                 "Some collage after dropout got my Ged","Some college and training after highschool.","Some college, culinary school, sommmelier school, medical coding school",
                                 "Start Nursing School June 1 at UMMC.","Still in college","Tech college grad","Technical School","Technical college Assoc. Applied Science","Trade school",
                                 "VOC - Diesel Engines - Boat & Locomotives","Vocational school","Vocational school industrical mechanic","Will complete college degree in Dec. 2016",
                                 "Will graduate this May","Working on Bachelors-Mechanical Engineering","currently enrolled at NWCC","currently in college","graduated from Lincoln Tech in Nashville",
                                 "graduated trade school","graduated vocational college","had some technical school traning","have a Master's degree in Cosmetology","in college","nursing school",
                                 "philosyphy 130/logics Heredity & evolution Astronomy. Anthropology misc....","some course in local college after high school","some graduate course",
                                 "still attending college","trade school","LMT","Current student","I have 2 career degrees.",
                                 "Jr College","AFTER HIGH SCHOOL ATTENDED COMMUNTY COLLEDGE WITH NO DEGREE. ATTENDED A 4 YEAR ELECTRICAL APPRENTICESHIP. RECEIVED JOURNEYMAN WIREMAN CERTIFICATE",
                                 "nursing school","Associates degree","Associates Degree in Electronics and Computer Technology",
                                 "Associate's Degree","Dropped out of high school GED 1981 Associate's Degree 1988","Adult College Handicap H.S.D",
                                 "2 year Degree Associate of Science Business Administration","some graduate course","Diesel Tech School","Vocational school","Dental Assistant/Trade school PIMA",
                                 "Completed Vocational training.","Flight school, Commercial Pilot","nursing school","Associate degree",
                                 "Associates in Electronic Engineering ITT Tech.","have a Master's degree in Cosmetology",
                                 "graduated trade school","Medical Assistant, Phlebotomy Tech, I.V. Tech certified",
                                 "Medical Asst. Training","Trade school","1 yr of college","Some college and training after highschool.",
                                 "I am in college right now.","Currently a junior in college","Associates of Arts degree A.A.",
                                 "Associates Degree 2 years","4th year in college","Associates in Science Degree","Some collage after dropout got my Ged","2nd year in college",
                                 "graduated vocational college","Currently enrolled in college","in college","Technical college Assoc. Applied Science","Currently in college",
                                 "Will complete college degree in Dec. 2016","Currently in college","Certified EKG/ECG Tech.",
                                 "Cosmetology School","Associates in Electronic Engineering ITT Tech.","Dropped out of high school GED 1981 Associate's Degree 1988", "USMC"
               )) {
      education[i] <- 2
    } else if (grade_i %in% 8 &
               gradespe_i %in% c("1 yr post grad","2 degrees - (a) mgmt (b) counseling","3x college degrees 1x postgraduate",
                                 "After completing my bachelor's in Health Sciences, I specialized in Radiologic Technology & Radiation Therapy.","BS - Nursing",
                                 "BS degree","Bachelor's in advertising","Currently enrolled in graduate classes","Currently in Pharmacy school","Currently in a master's program",
                                 "Currently in graduate school will have master's on 12-15-12","Currently in school to obtain master's degree","Doctorate-PharmD",
                                 "Ed. D. currently in progress","Graduate college in May 2016","Have 2 master's degrees 31 hrs above each masters called a rank I",
                                 "I am in my 3rd year of Pharmacy school","I have a Masters and will finish Ed. Admin Specialist in April 2015",
                                 "I have my master's degree-sorry I marked college grad.","In MA Ed program","Juirs Doctorate","Law School","M.D. completion on 05/19/2018",
                                 "M.S.; M.D.","MD","Master's Degree","Masters Degree","Masters Degree in School Administration","Masters FNP","Medical Degree","Medical Doctor",
                                 "Nurse","PhD Candidate","Phlebotomy school- graduated nursing school- graduated and passed boards","Registered Nurse (Associate)",
                                 "Some Masters level work","Working on masters. Will graduate May 2017","currently in graduate school","graduate with PharmD in 5/17",
                                 "halfway done with master's","l am a physician assistant","master's degree","starts med school Feb 2019","teaching licensure 1 year of master finishing next year",
                                 "two bachelor degrees","Specialized certification",
                                 "I have a Masters and will finish Ed. Admin Specialist in April 2015",
                                 "Masters FNP","Masters Degree in School Administration","Juirs Doctorate",
                                 "Master's Degree","Medical Degree","master's degree","Law School",
                                 "After completing my bachelor's in Health Sciences, I specialized in Radiologic Technology & Radiation Therapy.",
                                 "Currently enrolled in graduate classes","In MA Ed program","Ed. D. currently in progress",
                                 "Bachelor's in advertising","l am a physician assistant","currently in graduate school",
                                 "Medical Doctor","teaching licensure 1 year of master finishing next year",
                                 "Graduate college in May 2016","Doctorate-PharmD","I am in my 3rd year of Pharmacy school",
                                 "Currently in Pharmacy school")
    ) {
      education[i] <- 3
    } else if (grepl("jack of all trades", gradespe_i)) {
      education[i] <- 0
    } else if (is.na(grade_i) &
               gradespe_i %in% c("10th grade","Special Ed","Special Ed.","graduated Special Ed","Certificate","7 grade",
                                 "No","Special Diploma-Mentally Handicap","still attending high school","Haley is a senior",
                                 "Community Works job support offered by school district until age 21."
               )
    ) {
      education[i] <- 0
    } else if (is.na(grade_i) &
               gradespe_i %in% c("Some certificates working at hospital","Meat processing",
                                 "Manicurist school and medical assistant","Associate Degree",
                                 "Phlebotomy school- graduated nursing school- graduated and passed boards",
                                 "LPN","Still in college")) {
      education[i] <- 2
    } else if (mrn[i] == 12230) {
      education[i] <- 0
    } else {
      education[i] <- NA
    }
    
    if (grade_i %in% 8 &
        grepl("GED 1981|Electronic", gradespe_i)){
      education[i] <- 1
    }
  }
  
  education <- factor(education, levels = c(0:3), 
                      labels = c("Below high school",
                                 "High school/GED",
                                 "Training after high school/Some college",
                                 "College graduate/post-graduate"))
  
  return(education)
}
