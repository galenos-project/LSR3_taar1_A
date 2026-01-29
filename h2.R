hierarchical_subgroup_analysis <- function(df, experiment_type, outcome, rho_value) {
  # Debug: Print input parameters
  cat("Input - experiment_type:", experiment_type, "outcome:", outcome, "\n")
  
  df_filtered <- df %>% 
    filter(SortLabel == experiment_type) %>% 
    filter(outcome_type == outcome) %>%  
    filter(!is.na(SMDv))
  
  # Debug: Check filtering results
  cat("Records after filtering:", nrow(df_filtered), "\n")
  if (nrow(df_filtered) == 0) return(list(plotdata = tibble(), analyses = list()))
  
  cdi_levels <- unique(df_filtered$CDI_level2)
  cat("CDI_level2 groups found:", toString(cdi_levels), "\n")
  
  all_plotdata <- list()
  all_analyses <- list()
  
  for (cdi1 in cdi_levels) {
    df_cdi1 <- df_filtered %>% filter(CDI_level2 == cdi1)
    cat("\nProcessing", cdi1, "- Records:", nrow(df_cdi1), "\n")
    
    if (length(unique(df_cdi1$CDI_level3)) > 1) {
      cat("Multiple CDI_level3 subgroups found\n")
      res <- subgroup_analysis(df_cdi1, experiment_type, outcome, "CDI_level3", rho_value)
      
      if (!is.null(res) && !is.null(res$plotdata) && nrow(res$plotdata) > 0) {
        cat("Successfully processed", cdi1, "\n")
        res$plotdata$CDI_level2 <- cdi1
        res$plotdata <- res$plotdata %>% rename(subgroup = CDI_level3)
        all_plotdata[[cdi1]] <- res$plotdata
        all_analyses[[cdi1]] <- res$analysis
      } else {
        cat("Subgroup analysis failed or returned empty for", cdi1, "\n")
      }
    } else {
      cat("Insufficient subgroups in", cdi1, "\n")
    }
  }
  
  # Check final results
  if (length(all_plotdata) == 0) {
    cat("WARNING: No plotdata generated!\n")
  } else {
    cat("Generated plotdata for", length(all_plotdata), "groups\n")
  }
  
  return(list(plotdata = bind_rows(all_plotdata), analyses = all_analyses))
}
