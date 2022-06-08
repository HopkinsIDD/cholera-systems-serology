#' @title Get binary variable
#' 
#' @param indat data to pull from
#' @param group group number (1-5) as described in comments
#' @param cat category number (1-3)
#' 
#' @return vector with binary variable
get_binary_var <- function(indat, group, cat) {
  
  # testing group 1: Uninfected (no bacteria in stool) vs. Infected (bacteria in stool)
  if (group == 1) {
    indat[, binary_var := ifelse(get(paste0('category', cat)) == 'uninfected', 0, 1)]
  }
  
  # testing group 2: Symptomatic (symptoms) vs. Asymptomatic/Uninfected (No symptoms, regardless of bacteria in stool)
  if (group == 2) {
    indat[, binary_var := ifelse(get(paste0('category', cat)) == 'symptomatic', 1, 0)]
  }
  
  # testing group 3: Symptomatic (symptoms) vs. Uninfected (no bacteria in stool)
  if (group == 3) {
    indat[, binary_var := ifelse(get(paste0('category', cat)) == 'symptomatic', 1, 0)]
    indat[get(paste0('category', cat)) == 'asymptomatic', binary_var := NA]
  }
  
  # testing group 4: Asymptomatic (no symptoms and bacteria in stool) vs. Uninfected (no bacteria in stool)
  if (group == 4) {
    indat[, binary_var := ifelse(get(paste0('category', cat)) == 'asymptomatic', 1, 0)]
    indat[get(paste0('category', cat)) == 'symptomatic', binary_var := NA]
  }
  
  # testing group 5: Symptomatic (symptoms) vs. Asymptomatic (no symptoms and bacteria in stool)
  if (group == 5) {
    indat[, binary_var := ifelse(get(paste0('category', cat)) == 'symptomatic', 1, 0)]
    indat[get(paste0('category', cat)) == 'uninfected', binary_var := NA]
  }
  
  # end function
  return(indat$binary_var)
}


#' @title Get adjusted odds ratio using logistic regression
#' @description Logistic regression adjusting for age
#' 
#' @param indat full data table with biomarkers and outcome information
#' @param biomarker biomarker variable to use
#' @param outcome_group outcome group to use
#' @param outcome_cat outcome category to use
#' @param adjusted_values T/F whether to use adjusted values, default to TRUE
#' 
#' @return vector of odds ratio, confidence interval, and p value
glm_odds_ratio <- function(indat, biomarker, 
                           outcome_group, outcome_cat,
                           adjusted_values = TRUE) {
  
  # subset to biomarker
  indat <- indat[variable == biomarker]
  
  # set which version of the value to use
  value_col <- ifelse(adjusted_values, 'value_adj', 'value')
  
  # get binary variable for outcome
  binary_values <- get_binary_var(indat = indat, 
                                  group = outcome_group, 
                                  cat = outcome_cat)
  indat[, outcome := binary_values]
  
  # fit logistic regression
  fit <- glm(outcome ~ get(value_col) + yoa, 
             data = indat, family = binomial(link = 'logit'))
  
  # extract values
  coefs <- summary(fit)$coefficients
  ci <- suppressMessages(confint(fit))
  vals <- c(coefs[2,1], ci[2,1], ci[2,2])
  
  # convert to OR
  vals <- exp(vals)
  
  # return with p value
  return(c(vals, coefs[2,4]))
}


#' @title Get adjusted odds ratio using GEE
#' @description GEE adjusting for age and family/study
#' 
#' @param indat full data table with biomarkers and outcome information
#' @param biomarker biomarker variable to use
#' @param outcome_group outcome group to use
#' @param outcome_cat outcome category to use
#' @param adjusted_values T/F whether to use adjusted values, default to TRUE
#' 
#' @return vector of odds ratio and confidence interval
gee_odds_ratio <- function(indat, biomarker,
                           outcome_group, outcome_cat,
                           adjusted_values = TRUE) {
  
  # subset to biomarker
  indat <- indat[variable == biomarker]
  
  # set which version of the value to use
  value_col <- ifelse(adjusted_values, 'value_adj', 'value')
  
  # get binary variable for outcome
  binary_values <- get_binary_var(indat = indat, 
                                  group = outcome_group, 
                                  cat = outcome_cat)
  indat[, outcome := binary_values]
  
  # fit gee
  # http://www.personal.soton.ac.uk/dab1f10/MixedModels/Lecture4.pdf
  fit <- gee(outcome ~ get(value_col) + yoa, 
             id = as.factor(hh), data = indat, 
             family = binomial(link = 'logit'),
             silent = TRUE,
             corstr = 'unstructured')
  
  # extract values
  coefs <- summary(fit)$coefficients
  rob_se <- coefs[2, 4]
  ci <- coef(fit)['get(value_col)'] + c(-1,1)*rob_se*qnorm(0.975)
  vals <- exp(c(coefs[2,1], ci[1], ci[2]))
  
  # assign a dummy p value
  vals[4] <- ifelse(sum(vals < 1) == 3 | sum(vals > 1) == 3, 0, 1)
  
  # end function
  return(vals)
}


#' @title Get odds ratio table with biomarkers, ranked by effect size if desired
#' 
#' @param indat full data table with biomarkers and outcome information
#' @param vars biomarker variables to use
#' @param category outcome category to use (1-3)
#' @param outcome_group outcome group to use (1-5)
#' @param ranked T/F whether to rank odds ratios by effect size
#' @param model_type whether to calculate the OR using gee or glm ('gee' or 'glm')
#' @param adjusted_values whether to run adjusted values (TRUE), unadjusted (FALSE), or both (c(TRUE, FALSE))
#' @param formatted_values_only whether to return formatted values for a table (TRUE), or include raw estimates (FALSE)
#' @param return_pretty_table whether or not to return formatted table (if TRUE, returns kable, if FALSE returns data frame)
#' 
#' @return Kable table with ranked odds ratios, highlighting p < 0.05 in blue
odds_ratio_table <- function(indat, 
                             vars = unique(indat$variable),
                             category = 1,
                             outcome_group = 1,
                             ranked = FALSE,
                             model_type = 'glm',
                             adjusted_values = FALSE,
                             formatted_values_only = TRUE,
                             return_pretty_table = FALSE) {
  
  # table to write to
  tab <- data.table(
    biomarker = rep(vars, 2),
    adjusted_data = c(rep(FALSE, length(vars)), rep(TRUE, length(vars))),
    OR = '',
    p = ''
  )
  
  # loop over biomarkers
  for (v in vars) {
    print(v)
    # loop over adjusted/unadjusted data
    for (a in adjusted_values) {
      # get values
      if (model_type == 'glm') {
        tmp <- glm_odds_ratio(indat, biomarker = v, 
                              outcome_group = outcome_group,
                              outcome_cat = category,
                              adjusted_values = a)
      } else if (model_type == 'gee') {
        tmp <- suppressMessages( # don't need to print all the starting values
          gee_odds_ratio(indat, biomarker = v, 
                         outcome_group = outcome_group,
                         outcome_cat = category,
                         adjusted_values = a
                         )
        )
      } else {
        stop('Model type must be either "glm" or "gee"')
      }
      # format
      tmp <- round(tmp, 3)
      # add summary values to table
      tab[biomarker == v & adjusted_data == a, 
          OR := paste0(tmp[1], ' (', tmp[2], '-', tmp[3], ')')]
      tab[biomarker == v & adjusted_data == a, 
          p := tmp[4]]
      # add estimates to table
      tab[biomarker == v & adjusted_data == a, 
          OR_mean := tmp[1]]
      tab[biomarker == v & adjusted_data == a, 
          OR_lower := tmp[2]]
      tab[biomarker == v & adjusted_data == a, 
          OR_upper := tmp[3]]
    }
  }
  
  # return estimates if we don't want the formatted values
  if (formatted_values_only == FALSE) {
    
    return(tab[!is.na(OR_mean)])
    
  } else {
    # make pretty table
    tab2 <- cbind(tab[adjusted_data == FALSE],
                  tab[adjusted_data == TRUE])
    tab2 <- tab2[, c(1,3,4,7,8)]
    tab2[, biomarker := gsub('total|_Lx|_new_titers|_PS', ' ',
                             biomarker)]
    tab2[, biomarker := gsub('_', ' ', biomarker)]
    tab2[, biomarker := gsub('AUC', '', biomarker)]
    
    # order by effect size, if applicable
    if (ranked) setorderv(tab2, 'OR')
    setnames(tab2, c('Biomarker', 'Odds Ratio', 'p value', 'OR', 'p'))
    
    # make a pretty table
    if (return_pretty_table) {
      tab2 %>%
        # set up conditional formatting
        mutate_all(~cell_spec(.x, background = ifelse(.x < 0.05, '#80CDC1', 
                                                      ifelse(.x == 'Ogvctiter2 ', '#CC99CC', 
                                                             'white')))) %>%
        # make table with informative title and headers
        kable(escape = F, font = 12) %>%
        # # simple striped, bordered table
        kable_styling(bootstrap_options = c('bordered'), 
                      position = 'left', full_width = FALSE) %>%
        # rows are black
        row_spec(1:nrow(tab2), color='black') %>%
        # add header
        add_header_above(c(' ' = 1, 'Unadjusted plate data' = 2,
                           'Adjusted plate data' = 2)) %>%
        # make it scroll if it's large
        scroll_box(width = '100%', height = '500px')
    
      } else {
        
      # or return data frame
      return(tab2)
        
      }
    
  }
  
}


#' @title Perform GLM-PCA and unspervised feature selection
#' @description Selecting for two features and plotting the 2 components
#' 
#' @param abs_dt full data table with biomarkers only
#' @param cat vector of categories to color on plot
#' @param color_vals values for colors by category
#' @param return_data if TRUE, return GLM-PCA data in addition to plot
#' 
#' @return plot of principal components and pca data if desired
plot_glmpca <- function(abs_dt, cat, color_vals,
                        return_data = FALSE) {
  
  # add groups to data for plotting 
  abs_pca <- copy(abs_dt)
  abs_pca[, Outcome := cat]
  abs_pca <- na.omit(abs_pca)
  
  # report the amount that we're dropping
  message('Dropping ', (nrow(abs_dt)-nrow(abs_pca)), '/', nrow(abs_dt), ' (', 
          round((nrow(abs_dt)-nrow(abs_pca))/nrow(abs_dt)*100, 1),
          '%) of participants due to missing values in at least 1 biomarker')
  
  # transpose data so that observations are columns and features are rows
  abs_tr <- data.table::transpose(abs_pca[, -c('Outcome')])
  
  # select 2 features
  glmfit <- glmpca(abs_tr, 2, fam = 'nb')
  
  # get principal components
  dt_glmpca <- glmfit$factors
  dt_glmpca$Outcome <- abs_pca$Outcome
  
  # plot
  ggplot(dt_glmpca, aes(x = dim1, y = dim2, 
                        colour = Outcome)) +
    scale_color_manual(values = color_vals) +
    geom_point(size = 1) +
    theme_dark() -> p1
  print(p1)
  
  # end function
  if (return_data) return(glmfit)
}


#' @title Simple wrapper to fit different GEE models
#' 
#' @param fit_formula character string of formula to fit
#' @param indat data table with covariates and outcome variables
#' 
#' @return model fit object
fit_gee <- function(fit_formula, indat) {
  
  # fit model
  fit <- gee(as.formula(fit_formula), 
             id = as.factor(hh), data = indat, 
             family = binomial(link = 'logit'),
             silent = TRUE,
             corstr = 'unstructured')
  
  # return
  return(fit)
}


#' @title Extract GEE coefficients for a given covariate
#' 
#' @param variable variable to extract coefficients for
#' @param fit GEE model fit object
#' 
#' @return vector of exponentiated estimate, CI, and whether or not CI crosses 1
gee_coef <- function(variable, fit) {
  
  # extract coefficients
  coefs <- summary(fit)$coefficients
  rob_se <- coefs[which(rownames(coefs) == variable), 'Robust S.E.']
  var_est <- coefs[which(rownames(coefs) == variable), 'Estimate']
  ci <- var_est + c(-1,1)*rob_se*qnorm(0.975)
  vals <- exp(c(var_est, ci[1], ci[2]))
  
  # assign a dummy p value
  vals[4] <- ifelse(sum(vals < 1) == 3 | sum(vals > 1) == 3, 0, 1)
  vals <- round(vals, 2)
  
  # return
  return(vals)
}


#' @title Run random forest analysis
#' 
#' @param fit_dt data frame include all predictors to include, as well as outcome variable with col name "outcome"
#' @param variation_explained whether or not to print variation explained using RF regression
#' @param select_method selection method to pick biomarkers, options: 'rf', 'rf_regression', 'crf', 'crf_regression', 'plsda'
#' @param num_top_markers number of "top" markers to plot on the graph
#' @param print_rank_table whether or not to print all biomarkers ranked in a table
#' @param use_weights whether or not to use weights to account for imbalanced sample sizes when running random
#' 
#' @return plot and table of results
run_rf_analysis <- function(fit_dt,
                            variation_explained = FALSE,
                            select_method = 'crf',
                            num_top_markers = 20,
                            print_rank_table = FALSE,
                            use_weights = TRUE) {
  
  # set weights if using them
  if (use_weights) {
    wt_neg <- nrow(fit_dt[outcome==0])/nrow(fit_dt[outcome==1])
    wt_pos <- 1
  } else {
    wt_neg <- 1
    wt_pos <- 1
  }
  fit_dt[, class_wt := ifelse(outcome==1, wt_neg, wt_pos)]
  class_wt <- fit_dt$class_wt
  fit_dt$class_wt <- NULL
  
  # fit regression
  if (variation_explained) {
    rf <- ranger(as.numeric(outcome) ~ ., 
                 data = fit_dt,
                 num.trees = 1000,
                 replace = FALSE)
    message(round(rf$r.squared*100, 1), '% of variation in susceptibility can be explained by these variables (based on random forest regression model).')
  }

  # conditional random forest to select biomarkers
  if (select_method == 'crf') {
    crf <- cforest(as.factor(outcome)~.,
                   data=fit_dt,
                   control=cforest_unbiased(ntree=1000,mtry=ceiling(sqrt(ncol(fit_dt)-1))),
                   weights=class_wt)
    vimp <- permimp(crf, conditional = TRUE, asParty = TRUE, progressBar = FALSE)
    dt_vimp <- data.table(Biomarker=names(vimp$values),
                          `Conditional Importance`=as.numeric(vimp$values))
    setorderv(dt_vimp, 'Conditional Importance', order = -1)
  }
  
  # conditional random forest to select biomarkers
  if (select_method == 'crf_regression') {
    crf <- cforest(as.numeric(outcome)~.,
                   data=fit_dt,
                   control=cforest_unbiased(ntree=1000,mtry=ceiling(sqrt(ncol(fit_dt)-1))),
                   weights=class_wt)
    vimp <- permimp(crf, conditional = TRUE, asParty = TRUE, progressBar = FALSE)
    dt_vimp <- data.table(Biomarker=names(vimp$values),
                          `Conditional Importance`=as.numeric(vimp$values))
    setorderv(dt_vimp, 'Conditional Importance', order = -1)
  }
  
  # random forest to select using permutation 
  if (select_method == 'rf') {
    vrf <- ranger(as.factor(outcome) ~ ., 
                  data = fit_dt,
                  num.trees = 1000,
                  importance = 'permutation',
                  class.weights = c('0'=wt_pos, '1'=wt_neg),
                  replace = FALSE)
    vimp <- vrf$variable.importance
    dt_vimp <- data.table(Biomarker=names(vimp),
                          Importance=as.numeric(vimp))
    setorderv(dt_vimp, 'Importance', order = -1)
  }
  
  # random forest to select using permutation 
  if (select_method == 'rf_regression') {
    vrf <- ranger(as.numeric(outcome) ~ ., 
                  data = fit_dt,
                  num.trees = 1000,
                  importance = 'permutation',
                  class.weights = c('0'=wt_pos, '1'=wt_neg),
                  replace = FALSE)
    vimp <- vrf$variable.importance
    dt_vimp <- data.table(Biomarker=names(vimp),
                          Importance=as.numeric(vimp))
    setorderv(dt_vimp, 'Importance', order = -1)
  }
  
  # PLS-DA to select by VIP values
  if (select_method == 'plsda') {
    # center/scaling strongly advised for PLS-DA
    plsda_dt <- copy(fit_dt)
    cols <- names(plsda_dt[, -c('outcome')])
    plsda_dt[, (cols) := lapply(.SD, scale), .SDcols = cols]
    # run fit
    model <- plsda(plsda_dt[, -c('outcome')],
                   as.logical(plsda_dt$outcome), cv = 1, 
                   classname = 'protected')
    vimp <- vipscores(model)
    dt_vimp <- data.table(Biomarker=names(vimp[, 1]),
                          `VIP scores`=as.numeric(vimp[, 1]))
    setorderv(dt_vimp, 'VIP scores', order = -1)
  }
  
  # data table with importance values
  labels <- sapply(dt_vimp$Biomarker, format_biomarker)
  dt_vimp[, Biomarker := labels]
  
  # plot
  title_start <- names(dt_vimp)[2]
  dt_plot <- dt_vimp[1:num_top_markers, ]
  dt_plot[, Biomarker := factor(Biomarker, levels = rev(dt_plot$Biomarker))]
  ggplot(dt_plot, aes(x = Biomarker, y = get(title_start))) +
    geom_bar(stat = 'identity') + coord_flip() + xlab('') + ylab(title_start) +
    theme_minimal() + ggtitle(paste0(title_start, ' of top ', num_top_markers, ' variables/biomarkers')) -> p1
  print(p1)
  
  # full table
  if (print_rank_table) {
    tab <- cbind(data.table(Rank = 1:nrow(dt_vimp)), dt_vimp)
    tab %>%
      kable(font = 12, caption = paste0(title_start, ' for all variables/biomarkers')) %>%
      # # simple striped, bordered table
      kable_styling(bootstrap_options = c('bordered'), 
                    position = 'left', full_width = FALSE) %>%
      # rows are black
      row_spec(1:nrow(tab), color='black') %>%
      # make it scroll if it's large
      scroll_box(width = '100%', height = '500px')
  }
  
}


#' @title Run conditional random forest with leave-one-out cross-validation
#' 
#' @param variables either list of biomarkers, or 'topX' to select top X, or 'all' to include all biomarkers
#' @param indat input data table in long format
#' @param outcome_name name of the column to use as the outcome
#' @param use_weights whether or not to use weights to account for imbalanced sample sizes when running random
#' 
#' @return data table with detailed results
loo_crf <- function(variables, indat, outcome_name,
                    include_age = TRUE,
                    use_weights = TRUE,
                    quietly = TRUE) {
  
  # report variable
  if (quietly == FALSE) print(variables)
  
  # add outcome
  indat[, outcome := get(outcome_name)]
  
  # if applicable, run crf to select top X biomarkers using cRF (not including age)
  if (variables %like% 'top') {
    
    # cast data
    if(include_age) {
      if (!'yoa'%in%names(indat)) setnames(indat, 'age', 'yoa')
      dt_roc <- dcast(indat, sample + outcome + yoa ~ variable, value.var = 'value')
    } else {
      dt_roc <- dcast(indat, sample + outcome ~ variable, value.var = 'value')
    }
    dt_roc <- na.omit(dt_roc[, -c('sample')])
    
    # set weights if using them
    if (use_weights) {
      wt_neg <- nrow(dt_roc[outcome==0])/nrow(dt_roc[outcome==1])
      wt_pos <- 1
    } else {
      wt_neg <- 1
      wt_pos <- 1
    }
    dt_roc[, class_wt := ifelse(outcome==1, wt_neg, wt_pos)]
    class_wt <- dt_roc$class_wt
    dt_roc$class_wt <- NULL
    
    # run crf on full dataset to select top 3, excluding age from vimp
    crf <- cforest(as.factor(outcome)~.,
                   data=dt_roc,
                   control=cforest_unbiased(ntree=1000,mtry=ceiling(sqrt(ncol(dt_roc)-1))),
                   weights=class_wt)
    vimp <- permimp(crf, conditional = TRUE, asParty = TRUE, progressBar = FALSE)
    if (include_age) {
      vimp <- vimp$values[-which(names(vimp$values)=='yoa')]
    } else {
      vimp <- vimp$values
    }
    
    # subset
    variables <- names(sort(vimp, decreasing = TRUE))[1:as.numeric(gsub('top', '', variables))]
    indat <- indat[variable %in% variables]
    
    # otherwise, unless keeping all, subset to the variables indicated
  } else if (variables != 'all') {
    
    # subset
    indat <- indat[variable %in% variables]
  }
  
  # set up data for LOO validation
  if (include_age) {
    if (!'yoa'%in%names(indat)) setnames(indat, 'age', 'yoa')
    dt_roc <- dcast(indat, sample + outcome + yoa ~ variable, value.var = 'value')
  } else {
    dt_roc <- dcast(indat, sample + outcome ~ variable, value.var = 'value')
  }
  dt_roc <- na.omit(dt_roc)
  
  # unique ids
  unique_ids <- dt_roc$sample
  
  # true values
  true_vals <- dt_roc$outcome
  
  # data frame to save results to
  res <- data.table(sample = unique_ids,
                    outcome_label = outcome_name,
                    outcome_value = true_vals,
                    prediction = '',
                    prob0 = '', prob1 = '',
                    biomarkers = ifelse(length(variables > 1), 
                                        paste0(variables, collapse= ', '), variables))
  
  # loop over unique IDs and run LOO validation
  for (i in unique_ids) {
    
    # split train/test sets
    dt_train <- dt_roc[!(sample %in% i)]
    dt_test <- dt_roc[(sample %in% i)]
    
    # set weights if using them
    if (use_weights) {
      wt_neg <- nrow(dt_train[outcome==0])/nrow(dt_train[outcome==1])
      wt_pos <- 1
    } else {
      wt_neg <- 1
      wt_pos <- 1
    }
    dt_train[, class_wt := ifelse(outcome==1, wt_neg, wt_pos)]
    class_wt <- dt_train$class_wt
    dt_train$class_wt <- NULL
    
    # run conditional random forest on training dataset
    crf <- cforest(as.factor(outcome)~.,
                   data=dt_train[, -c('sample')],
                   control=cforest_unbiased(ntree=1000,mtry=ceiling(sqrt(ncol(dt_train)-1))),
                   weights=class_wt)
    
    # predict with testing dataset
    crf_pred <- predict(crf, newdata=dt_test[, -c('sample')])
    crf_prob <- predict(crf, newdata=dt_test[, -c('sample')], type = 'prob')
    
    # save values
    res[sample==i, prediction := crf_pred]
    res[sample==i, prob0 := as.numeric(crf_prob[[1]][1])]
    res[sample==i, prob1 := as.numeric(crf_prob[[1]][2])]
  }
  
  # return results
  return(res)
  
}


#' @title Run classification using an ensemble of penalized logistic regression, BRT, and SVM
#' 
#' @param variables either list of biomarkers, or 'topX' to select top X, or 'all' to include all biomarkers
#' @param data input data table in long format
#' @param outcome_name name of the column to use as the outcome
#' @param use_weights whether or not to use weights to account for imbalanced sample sizes when running BRT
#' 
#' @return data table with detailed results, including ensemble and individual models
loo_ensemble <- function(variables, data, outcome_name,
                         include_age = TRUE,
                         use_weights = TRUE,
                         quietly = TRUE) {
  
  # report variable
  if (quietly == FALSE) print(variables)
  
  # add outcome
  indat <- copy(data)
  indat[, outcome := get(outcome_name)]
  
  # if applicable, run crf to select top X biomarkers using cRF (not including age)
  if (variables %like% 'top') {
    
    # cast data
    if(include_age) {
      if (!'yoa'%in%names(indat)) setnames(indat, 'age', 'yoa')
      dt_roc <- dcast(indat, sample + outcome + yoa ~ variable, value.var = 'value')
    } else {
      dt_roc <- dcast(indat, sample + outcome ~ variable, value.var = 'value')
    }
    dt_roc <- na.omit(dt_roc[, -c('sample')])
    
    # set weights if using them
    if (use_weights) {
      wt_neg <- nrow(dt_roc[outcome==0])/nrow(dt_roc[outcome==1])
      wt_pos <- 1
    } else {
      wt_neg <- 1
      wt_pos <- 1
    }
    dt_roc[, class_wt := ifelse(outcome==1, wt_neg, wt_pos)]
    class_wt <- dt_roc$class_wt
    dt_roc$class_wt <- NULL
    
    # run crf on full dataset to select top 3, excluding age from vimp
    crf <- cforest(as.factor(outcome)~.,
                   data=dt_roc,
                   control=cforest_unbiased(ntree=1000,mtry=ceiling(sqrt(ncol(dt_roc)-1))),
                   weights=class_wt)
    vimp <- permimp(crf, conditional = TRUE, asParty = TRUE, progressBar = FALSE)
    if (include_age) {
      vimp <- vimp$values[-which(names(vimp$values)=='yoa')]
    } else {
      vimp <- vimp$values
    }
    
    # subset
    variables <- names(sort(vimp, decreasing = TRUE))[1:as.numeric(gsub('top', '', variables))]
    indat <- indat[variable %in% variables]
    
    # otherwise, unless keeping all, subset to the variables indicated
  } else if (variables != 'all') {
    
    # subset
    indat <- indat[variable %in% variables]
  }
  
  # set up data for LOO validation
  if (include_age) {
    if (!'yoa'%in%names(indat)) setnames(indat, 'age', 'yoa')
    dt_roc <- dcast(indat, sample + outcome + yoa ~ variable, value.var = 'value')
  } else {
    dt_roc <- dcast(indat, sample + outcome ~ variable, value.var = 'value')
  }
  dt_roc <- na.omit(dt_roc)
  
  # unique ids
  unique_ids <- dt_roc$sample
  
  # true values
  true_vals <- dt_roc$outcome
  
  # data frame to save results to
  res <- data.table(sample = unique_ids,
                    outcome_label = outcome_name,
                    outcome_value = true_vals,
                    prob1 = '',
                    coef_rf = '', prob1_rf = '',
                    coef_plr = '', prob1_plr = '',
                    coef_svm = '', prob1_svm = '',
                    biomarkers = ifelse(length(variables > 1), 
                                        paste0(variables, collapse= ', '), variables))
  
  # loop over unique IDs to run LOO validation
  for (i in unique_ids) {
    
    if (quietly == FALSE) print(i)
    
    # split train/test sets
    dt_train <- dt_roc[!(sample %in% i)]
    dt_test <- dt_roc[(sample %in% i)]
    
    # set weights if using them
    if (use_weights) {
      wt_neg <- nrow(dt_train[outcome==0])/nrow(dt_train[outcome==1])
      wt_pos <- 1
    } else {
      wt_neg <- 1
      wt_pos <- 1
    }
    dt_train[, class_wt := ifelse(outcome==1, wt_neg, wt_pos)]
    class_wt <- dt_train$class_wt
    dt_train$class_wt <- NULL
    
    # super learner
    #TODO: reduce the number of variables?
    sl = SuperLearner(Y = dt_train$outcome, 
                      X = dt_train[, -c('sample', 'outcome')], 
                      family = binomial(),
                      SL.library = c('SL.glmnet', 'SL.svm', 'SL.ranger'),
                      obsWeights = class_wt)

    # print time elapsed
    if (quietly == FALSE) print(sl$times$everything[3])
        
    # extract coefficients and predictions
    coef <- sl$coef
    sl_pred <- predict(sl, dt_test[, -c('sample', 'outcome')], onlySL = TRUE)
    sl_lib <- sl_pred$library.predict
        
    # add to results table
    res[sample == i, prob1 := sl_pred$pred]
    res[sample == i, coef_rf := coef[3]] 
    res[sample == i, prob1_rf := sl_lib[,3]]
    res[sample == i, coef_plr := coef[1]]
    res[sample == i, prob1_plr := sl_lib[,1]]
    res[sample == i, coef_svm := coef[2]]
    res[sample == i, prob1_svm := sl_lib[,2]]

  }

  # return results
  return(res)
  
}


#' @title Extract isotype from raw biomarker label
#' 
#' @param label biomarker label in the raw dataset
#' 
#' @return data table with detailed results, including ensemble and individual models
get_isotype <- function(label) {
  
  # return functional response for those that aren't isotype specific
  if (label %in% c('vibinab', 'vibogaw', 'Ogvctiter2_new_titers', 
                   'ADNP_PS', 'ADCD', 'ADCP', 'ADCP_PS')) {
    
    # vibriocidals
    if (label %in% c('vibinab', 'vibogaw')) newlab <- 'Vibriocidal'
    if (label == 'Ogvctiter2_new_titers') newlab <- 'Vibriocidal (new)'
    
    # other functional responses
    if (label %in% c('ADNP_PS', 'ADCD', 'ADCP', 'ADCP_PS')) newlab <- gsub('_PS', '', label)
  
  # pull first 3 letters and any subsequent numbers for luminex markers  
  } else {
    # isotype
    iso <- substr(label, 1, 3)
    # subtype
    subtype <- str_extract(label, '[[:digit:]]+')
    # combine
    newlab <- ifelse(!is.na(subtype), paste0(iso, subtype), iso)
  }
  
  # end function
  return(newlab)
}


#' @title Extract antigen from raw biomarker label
#' 
#' @param label biomarker label in the raw dataset
#' 
#' @return data table with detailed results, including ensemble and individual models
get_antigen <- function(label) {
  
  # return functional response for those that aren't isotype specific
  if (label %in% c('vibinab', 'vibogaw', 'Ogvctiter2_new_titers', 'vib4fold',
                   'ADNP_PS', 'ADCD', 'ADCP', 'ADCP_PS', 'vib_fold')) {
    ag <- 'Functional responses'
    
  # manually assign antigens
  } else {
    if (grepl('CtxB', label)) ag <- 'CtxB'
    if (grepl('Sial', label)) ag <- 'Sialidase'
    if (grepl('TcpA', label)) ag <- 'TcpA'
    if (grepl('CTH', label)) ag <- 'CT-HT'
    if (grepl('OSP', label) & grepl('In', label)) ag <- 'Inaba-OSP'
    if (grepl('OSP', label) & grepl('Og', label)) ag <- 'Ogawa-OSP'
  }
  
  # end function
  return(ag)
}


#' @title Format biomarker label into something readable
#' 
#' @param label biomarker label in the raw dataset
#' 
#' @return data table with detailed results, including ensemble and individual models
format_biomarker <- function(label) {
  
  # functional responses and age
  if (label %in% c('vibinab', 'vibogaw', 'Ogvctiter2_new_titers', 'vib_fold', 'vib4fold',
                   'ADNP_PS', 'ADCD', 'ADCP', 'ADCP_PS', 'yoa', 'age')) {
    
    # vibriocidals
    if (label %in% c('vibinab', 'vibogaw')) newlab <- 'Vibriocidal'
    if (label %in% c('vib_fold')) newlab <- 'Vibriocidal fold increase'
    if (label %in% c('vib4fold')) newlab <- 'Vibriocidal 4-fold increase'
    if (label == 'Ogvctiter2_new_titers') newlab <- 'Vibriocidal (new)'
    
    # other functional responses
    if (label %in% c('ADNP_PS', 'ADCD', 'ADCP', 'ADCP_PS')) newlab <- gsub('_PS', '', label)
  
    # age
    if (label %in% c('yoa', 'age')) newlab <- 'age'
    
  # luminex biomarkers  
  } else {
    newlab <- paste0(get_antigen(label), ' ', get_isotype(label))
  }
  
  # end function
  return(newlab)
}


##' Pair-wise spearman correlations between biomarkers
##' 
##' @param indat Data frame with biomarker data
##' @param return_values Whether or not to return tables with corr and CI
##' 
run_corr_analysis <- function(indat, return_values = FALSE) {
  
  # rename
  names(indat) <- sapply(names(indat), format_biomarker)
  
  # correlation matrix
  cor_mat <- cor(indat, 
                 method = 'spearman',
                 use = 'pairwise.complete.obs')
  
  # p values
  testRes <- cor.mtest(indat, conf.level = 0.95)
  
  # return values, if needed
  if (return_values) {
    print(cor_mat)
    print(testRes)
  }
  
  # correlation plot
  ggcorrplot(cor_mat, lab = FALSE, ggtheme=ggplot2::theme_bw,
             colors = brewer.pal(n = 9, name = 'PRGn')[c(1,5,9)],
             tl.col = 'black', tl.srt = 45, insig = 'pch',
             pch.cex = 2, tl.cex = 8,
             p.mat = testRes$p, sig.level = 0.05) -> p1
  
  # center title
  p1 <- p1 + 
    theme(plot.title = element_text(hjust = 0.5))
  
  # return plot
  return(p1)
}

#' @title Run conditional random forest with leave-one-out cross-validation, training with one dataset and predicting with the other
#' 
#' @param variables either list of biomarkers, or 'topX' to select top X, or 'all' to include all biomarkers
#' @param train_dat input training data table in long format
#' @param test_dat input testing data table in long format
#' @param outcome_name name of the column to use as the outcome
#' @param use_weights whether or not to use weights to account for imbalanced sample sizes when running random
#' 
#' @return data table with detailed results
cross_pred <- function(variables, 
                       train_dat, 
                       test_dat,
                       train_outcome,
                       test_outcome,
                       include_age = TRUE,
                       use_weights = TRUE,
                       quietly = TRUE) {
  
  # report variable
  if (quietly == FALSE) print(variables)
  
  # add outcome
  train_dat[, outcome := get(train_outcome)]
  test_dat[, outcome := get(test_outcome)]
  
  # subset
  train_dat <- train_dat[variable_clean %in% variables]
  test_dat <- test_dat[variable_clean %in% variables]
  
  # set up data for LOO validation
  if (include_age) {
    if (!'yoa'%in%names(train_dat)) setnames(train_dat, 'age', 'yoa')
    if (!'yoa'%in%names(test_dat)) setnames(test_dat, 'age', 'yoa')
    dt_train <- dcast(train_dat, sample + outcome + yoa ~ variable_clean, value.var = 'value')
    dt_test <-  dcast(test_dat, sample + outcome + yoa ~ variable_clean, value.var = 'value')
  } else {
    dt_train <- dcast(train_dat, sample + outcome ~ variable_clean, value.var = 'value')
    dt_test <- dcast(test_dat, sample + outcome ~ variable_clean, value.var = 'value')
  }
  dt_train <- na.omit(dt_train)
  dt_test <- na.omit(dt_test)
  
  # unique ids for loo validation
  unique_ids <- dt_test$sample
  
  # true values
  true_vals <- dt_test$outcome
  
  # data frame to save results to
  res <- data.table(sample = unique_ids,
                    outcome_label = test_outcome,
                    outcome_value = true_vals,
                    prediction = '',
                    prob0 = '', prob1 = '',
                    biomarkers = ifelse(length(variables > 1), 
                                        paste0(variables, collapse= ', '), variables))
  
  # set weights if using them
  if (use_weights) {
    wt_neg <- nrow(dt_train[outcome==0])/nrow(dt_train[outcome==1])
    wt_pos <- 1
  } else {
    wt_neg <- 1
    wt_pos <- 1
  }
  dt_train[, class_wt := ifelse(outcome==1, wt_neg, wt_pos)]
  class_wt <- dt_train$class_wt
  dt_train$class_wt <- NULL
  
  # run conditional random forest on training dataset
  crf <- cforest(as.factor(outcome)~.,
                 data=dt_train[, -c('sample')],
                 control=cforest_unbiased(ntree=1000,mtry=ceiling(sqrt(ncol(dt_train)-1))),
                 weights=class_wt)
  
  # loop over unique IDs and run LOO validation
  for (i in unique_ids) {
    
    # subset
    test_row <- dt_test[sample == i]
    
    # predict with testing dataset
    crf_pred <- predict(crf, newdata=test_row[, -c('sample')])
    crf_prob <- predict(crf, newdata=test_row[, -c('sample')], type = 'prob')
    
    # save values
    res[sample==i, prediction := crf_pred]
    res[sample==i, prob0 := as.numeric(crf_prob[[1]][1])]
    res[sample==i, prob1 := as.numeric(crf_prob[[1]][2])]
  }
  
  # return results
  return(res)
  
}

