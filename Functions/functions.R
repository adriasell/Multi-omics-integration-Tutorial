bootstrap_figure <- function(bootstrap,
                             type="loadings",
                             block,
                             stratified_by=NULL,
                             comp=1,
                             color= "#000000" ){
  
  ## boostrap = boostrap data ($stats), or list of boostrap stats for stratified plots
  ## type = "weights" or "loadings"
  ## block = name of the block, as used in the model
  ## stratified = T if you want stratified plots (combine multiple forest plots)
  defined_block = block
  defined_type = type
  defined_comp = comp
  
  # Not stratified
  if (is.null(stratified_by)){
    data_bootstrap <- bootstrap %>%
    dplyr::rename(
      variable = .data$var,
      beta = .data$estimate,
      `CI 2.5` = .data$lower_bound,
      `CI 97.5` = .data$upper_bound
    ) %>%
    dplyr::filter(block == defined_block & comp == defined_comp & type == defined_type) %>%
    dplyr::mutate(variable = forcats::fct_reorder(variable, abs(beta)))
    
    
    boostrap_plot <- forest_plot(data_bootstrap, color = color)
    boostrap_plot
    
    # Stratified
  } else {
    data_bootstrap <- bootstrap %>%
      filter(block == defined_block, comp == defined_comp, type == defined_type) %>%
      mutate(var = fct_reorder(var, abs(mean)))
    
    group <- stratified_by
    boostrap_plot <- ggplot(data_bootstrap, aes(var, mean, fill= Group)) +
      # range from CI25% to CI75%, colored by group (stratification)
      geom_pointrange(
        aes(ymin = lower_bound, ymax = upper_bound, color = color,
            shape=Group),
        position = position_dodge(0.7),size=.5)+
      # vertical line at y=0
      geom_hline(yintercept=0, lty=2)+
      scale_colour_paletteer_d("ggthemes::wsj_red_green") +
      # flip the coordinates
      coord_flip() +
      # Parameters for x and y axis
      # Label of y axis
      ylab(paste0(type," for comp",comp," of block ",block)) +
      theme_bw()
    
    
    boostrap_plot
    
    
  }
  return(boostrap_plot)
}


#-------------------------------------------------------------------------------

forest_plot= function(results, color){
  # Purpose : from the multivariate model, plot for each 
  #               exposure (y-axis): beta and CI (x-axis)
  # Inputs: results should be in the format
  # - variable: all exposures, ploted on the x-axis
  # - label: label of the exposure
  # - beta: estimated coefficient of the regression
  # - sd: standard deviation of the estimation
  #   and if color=T: if is a risk or a protective factor.
  # Output : forest plot
  
  # Chose the colors to use: red for risk factors and green for protective factor
  # (gray for insignificant associations)
  
  results <- mutate(results,
                    beta=as.numeric(beta),
                    IC_low=as.numeric(`CI 2.5`),
                    IC_high=as.numeric(`CI 97.5`))
  
  # Plot of estimated coefficient (y) in function of the exposure (x)
  p <- ggplot(results, aes(variable, beta)) +
    # Point at beta + range from CI25 (beta-1.96*sd) to CI75 (beta+1.96*sd)
    geom_pointrange(
      aes(x=variable,ymin = IC_low, ymax = IC_high),
      position = position_dodge(0.5),size=1, color = color)+
    # Vertical line at y=0
    geom_hline(yintercept=0, lty=2,color="black")+
    coord_flip() +
    # Limits of the y axis
    ylim(min(as.numeric(results$`CI 2.5`))-0.02,max(as.numeric(results$`CI 97.5`))+0.02)+ 
    # Label of y axis
    ylab("Beta (CI 95%)")+
    xlab("")+
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 24),
      axis.text.x = element_text(size = 22),
      axis.text.y = element_text(size = 22) 
    )  
  return(p)
}


#Cross validation function
cross_validation_single_outcome <- function(X_test, response, rgcca_res, Nfold = 5, n_run = 10) {
  
  # Function to run one repetition of cross-validation
  cross_val_run <- function(rep) {
    id_cv <- sample(1:Nfold, nrow(X_test[[1]]), replace = TRUE)
    
    quality_cv_list <- list()  # store results for each fold
    
    for (i in 1:Nfold) {
      # Prepare test blocks
      X_test_i <- lapply(X_test, function(block) block[id_cv == i, , drop = FALSE])
      
      # RGCCA prediction
      pred_quality <- rgcca_predict(rgcca_res, blocks_test = X_test_i, prediction_model = "lm")
      
      # Combine latent variables from projection
      latent_variables_test <- purrr::reduce(pred_quality$projection, cbind)
      colnames(latent_variables_test) <- paste0("Comp", seq_len(ncol(latent_variables_test)))
      
      # Correlation with outcome
      cor_values <- cor(latent_variables_test, X_test_i[[response]])[,1]
      
      # Combine latent variables and outcome
      data_r2 <- cbind(latent_variables_test, X_test_i[[response]]) %>% as.data.frame()
      
      # Get actual outcome column name
      outcome_name <- colnames(X_test_i[[response]])
      
      # Compute R2 for each latent variable
      r2_values <- sapply(colnames(latent_variables_test), function(comp) {
        summary(lm(as.formula(paste0(outcome_name, " ~ ", comp)), data = data_r2))$r.squared
      })
      
      # Combine R2 and correlation
      fold_result <- c(r2_values, cor_values)
      names(fold_result) <- c(paste0("R2_", colnames(latent_variables_test)), paste0("cor_", colnames(latent_variables_test)))
      
      quality_cv_list[[i]] <- fold_result
    }
    
    # Combine folds into a matrix
    do.call(rbind, quality_cv_list)
  }
  
  # Run multiple repetitions
  res <- lapply(1:n_run, cross_val_run)
  res_combined <- do.call(rbind, res)
  
  # Convert to long format for plotting
  res_long <- res_combined %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Fold") %>%
    tidyr::pivot_longer(cols = -Fold, names_to = "Indicator", values_to = "value") %>%
    mutate(group = ifelse(grepl("cor_", Indicator), "Correlation with outcome", "R2 of components"))
  
  # Compute medians for labels
  medians <- res_long %>%
    group_by(Indicator, group) %>%
    summarise(md = median(value), .groups = "drop") %>%
    mutate(md = round(md, 3))
  
  # Plot
  final_plot <- ggplot(res_long, aes(x = Indicator, y = value)) +
    geom_boxplot() +
    ggforce::facet_row(vars(group), scales = 'free', space = 'free') +
    ggrepel::geom_text_repel(data = medians, aes(x = Indicator, y = abs(md), label = abs(md)), direction = "y") +
    xlab("Cross-validation folds") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  return(list(plot = final_plot, median_quality = medians, quality_all = res_long))
}
