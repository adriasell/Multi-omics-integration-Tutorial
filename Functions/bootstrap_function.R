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
      rename(variable=var,
             beta = estimate,
             `CI 2.5`=lower_bound,
             `CI 97.5` = upper_bound) %>%
      filter(block == defined_block & comp==defined_comp & type == defined_type)%>%
      mutate(variable = fct_reorder(variable, abs(beta)))
    
    
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

