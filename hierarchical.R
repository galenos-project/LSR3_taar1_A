
hierarchical_subgroup_analysis <- function(df, experiment_type, outcome, rho_value) {
  # Filter data according to your requirements
  df_filtered <- df %>% 
    filter(SortLabel == experiment_type) %>% 
    filter(outcome_type == outcome) %>%  
    filter(!CDI_level2 == 'DAT KO') %>%
    filter(!is.na(SMDv))
  
  df_filtered <- df_filtered %>% mutate(effect_id = row_number()) # add effect_id column
  
  #calculate variance-covariance matrix of the sampling errors for dependent effect sizes
  factors_to_consider <- c("Strain", "StudyId", "ExperimentID_I")
  
  # Create the random effects formula
  random_formula <- create_formula(factors_to_consider, df_filtered)
  
  VCVM_SMD <- vcalc(vi = SMDv,
                    cluster = StudyId, 
                    subgroup= ExperimentID_I,
                    obs=effect_id,
                    data = df_filtered, 
                    rho = rho_value) 
  
  ## ML model on df_filtered without subgroup
  global_analysis <- rma.mv(yi = SMD,
                                 V = VCVM_SMD,
                                 random = random_formula, # nested levels
                                 test = "t", # use t- and F-tests for making inferences
                                 data = df_filtered,
                                 rho = rho_value,
                                 dfs="contain", # improve degree of freedom estimation for t- and F-distributions
                                 control=list(optimizer="nlminb"))
  
  k_subgroups <- df_filtered %>%
    count() %>%
    pull(n)
  
  
  global_analysis_plotdata <- data.frame('Global estimate', k_subgroups, global_analysis$beta, global_analysis$se, global_analysis$pval, global_analysis$ci.lb, global_analysis$ci.ub)
  colnames(global_analysis_plotdata) <- c('subgroup', "k", "SMD", "se","p", "ci_l", "ci_u") #, "pi.lb", "pi.ub")
  global_analysis_plotdata$symbol <- 15
  global_analysis_plotdata$size <- (1/global_analysis_plotdata$se)
  global_analysis_plotdata$summary <- 'TRUE'
  global_analysis_plotdata$fontfaace <- "bold"
  global_analysis_plotdata$fontsize <- 7
  #global_analysis_plotdata=rbind(global_analysis_plotdata, c("Overall estimate",  overall_estimate_rma$k, overall_estimate_rma$beta, overall_estimate_rma$se, overall_estimate_rma$pval, overall_estimate_rma$ci.lb,overall_estimate_rma$ci.ub, 18,1,TRUE,"bold",5)) #overall_estimate_rma_predict$pi.lb, overall_estimate_rma_predict$pi.ub))
  
  
  
  
  rownames(global_analysis_plotdata) <- 1:nrow(global_analysis_plotdata)
  global_analysis_plotdata$k <- as.numeric(global_analysis_plotdata$k)
  global_analysis_plotdata$SMD <- as.numeric(global_analysis_plotdata$SMD)
  global_analysis_plotdata$se <- as.numeric(global_analysis_plotdata$se)
  global_analysis_plotdata$ci_l <- as.numeric(global_analysis_plotdata$ci_l)
  global_analysis_plotdata$ci_u <- as.numeric(global_analysis_plotdata$ci_u)
  global_analysis_plotdata$p <- as.numeric(global_analysis_plotdata$p)
  global_analysis_plotdata$symbol <- as.numeric(global_analysis_plotdata$symbol)
  global_analysis_plotdata$size <- as.numeric(global_analysis_plotdata$size)
  global_analysis_plotdata$summary <- as.character(global_analysis_plotdata$summary)
  global_analysis_plotdata$fontsize <- as.numeric(global_analysis_plotdata$fontsize)
  
  global_analysis_plotdata$d1 <- (global_analysis_plotdata$SMD - global_analysis_plotdata$ci_l)/1.92
  global_analysis_plotdata$d2 <- (global_analysis_plotdata$ci_u - global_analysis_plotdata$SMD)/1.92
  global_analysis_plotdata$CDI_level2 <- ''
  
  
  # Get unique CDI_level1 values
  cdi_levels <- unique(df_filtered$CDI_level2)
  
  # Initialize lists to store results
  all_plotdata <- list()
  all_analyses <- list()
  
  all_plotdata[['global']] <- global_analysis_plotdata
  all_analyses[['global']] <- global_analysis
  
  # Process each CDI_level1
  for (cdi1 in cdi_levels) {
    df_cdi1 <- df_filtered %>% filter(CDI_level2 == cdi1)
    
    # Check if we have CDI_level3 subgroups
    if (length(unique(df_cdi1$CDI_level3)) > 1) {
      # Perform subgroup analysis on CDI_level3 within this CDI_level1
      res <- subgroup_analysis(df_cdi1, experiment_type, outcome, "CDI_level3", rho_value)
      
      if (!is.null(res)) {
        # Add CDI_level2 information
        res$plotdata$CDI_level2 <- cdi1
        
        # Rename the moderator column to "subgroup" for consistency
        # The column name will be "CDI_level3" from the subgroup_analysis
        res$plotdata <- res$plotdata %>%
          rename(subgroup = CDI_level3)
        
        # Store results
        all_plotdata[[cdi1]] <- res$plotdata
        all_analyses[[cdi1]] <- res$analysis
      }
    } else {
      # If only one CDI_level3, create a simple overall analysis for this CDI_level1
      # This is a simplified version of your subgroup_analysis without moderator
      # (Implementation would go here)
    }
  }
  
  # Combine all plotdata
  combined_plotdata <- bind_rows(all_plotdata)
  
  combined_plotdata <- combined_plotdata %>%
               mutate(subgroup = case_when(subgroup == 'Overall estimate' ~ 'Subgroup estimate',
                                           subgroup == 'Global estimate' ~ 'Overall estimate',
                       TRUE ~ subgroup))
  lt <- nrow(combined_plotdata)
  combined_plotdata <- combined_plotdata[c(2:lt,1),]
  
  return(list(plotdata = combined_plotdata, analyses = all_analyses))
}

forest_hierarchical <- function(plotdata, outcome) {
  # Check if plotdata is empty
  if (nrow(plotdata) == 0) {
    stop("plotdata is empty - cannot create forest plot")
  }
  
  plotdata <- plotdata %>%
    mutate(size = case_when(
      str_detect(subgroup, "Subgroup") ~ 10,
      TRUE ~ `size`  # keep original value if no match
    ))
  
  # Prepare data for plotting - create the display labels
  plotdata_prepared <- plotdata %>%
    mutate(subgroup = case_when(
      str_detect(subgroup, "Subgroup estimate") ~ CDI_level2,
                     TRUE ~ subgroup))
  plotdata_prepared$order <- 1
  
  # Reorder factors for proper display
  plotdata_prepared <- plotdata_prepared %>%
    mutate(display_label = factor(display_label, levels = unique(display_label)))
  
  # Create a copy with the moderator column renamed for the forest plot function
  plotdata_for_forest <- plotdata_prepared %>%
    mutate(subgroup = display_label) %>%
    select(-display_label)  # Remove the temporary column
  
  # Call the forest plot function
  
  p <- forest_subgroup_ml(plotdata_prepared, "moderator", outcome, "CDI Level")
  
  return(p)
}

# Usage example:
# hierarchical_results <- hierarchical_subgroup_analysis(df, "your_experiment_type", "Locomotor Activity", 0.8)
# forest_plot <- forest_hierarchical(hierarchical_results$plotdata, "Locomotor Activity")
# print(forest_plot)