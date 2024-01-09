filter_experiment_outcome_type <- function(df, experiment_type, outcome) {
  df_by_outcome <- df %>%
    filter(outcome_type == outcome) %>%
    filter(!is.na(SMD) | !is.na(SMDv))
  
  df_by_experiment_outcome <- df_by_outcome %>%
    filter(SortLabel == experiment_type)
  
  return(df_by_experiment_outcome)
}

run_ML_SMD <- function(df, experiment, outcome, rho_value) {

  df<-filter_experiment_outcome_type(df, experiment, outcome)

  df<-df %>% 
    filter(!is.na(SMDv))
  
  # Check if data comes from 2 or more unique studies
  if (n_distinct(df$StudyId) > 1) {
  
  df <- df %>% mutate(effect_id = row_number()) # add effect_id column
  
  #calculate variance-covariance matrix of the sampling errors for dependent effect sizes
  
  VCVM_SMD <- vcalc(vi = SMDv,
                    cluster = StudyId, 
                    subgroup= ExperimentID_I,
                    obs=effect_id,
                    data = df, 
                    rho = rho_value)
  
  SMD_ML <- rma.mv(yi = SMD,
                   V = VCVM_SMD,
                   random = ~1 | Strain / StudyId / ExperimentID_I, # nested levels
                   test = "t", # use t- and F-tests for making inferences
                   data = df,
                   dfs="contain",
                   control=list(optimizer="nlm")
  )
  
  #if (length(unique(df$StudyId)) > 1) {
    #SMD_ML <- robust(SMD_ML, cluster = StudyId, clubSandwich = FALSE)
  #}
  
  cat("Meta analysis summary:\n")
  print(summary(SMD_ML))
  
  cat("\n-------------------------\n")
  
  cat("Prediction Interval:\n")
  
  pred_interval <- predict(SMD_ML)
  
  print(pred_interval)
  
  return(SMD_ML)
}}

forest_metafor <- function(model, experiment_type, outcome_title){ #outcome title is what you want outcome to be written as, it doesn't have to match outcome type
  if(!is.null(model)){
  
  
  lower_x <- floor((min(model[["yi"]])-mean(model[["vi"]])) - 1)
  upper_x <- ceiling((max(model[["yi"]])+mean(model[["vi"]])) + 1)
  summary_x <- model[["beta"]]
  
  at_values <- seq(floor(lower_x / 5) * 5, ceiling(upper_x / 5) * 5, by = 5)
  
         forest_plot <- if(experiment_type == "TvC"){
                               forest(model,
                                      xlim=c((lower_x-4), (upper_x+2)),
                                      mlab="SMD [95% C.I.]", 
                                      alim=c((lower_x-4), (upper_x+2)),
                                      slab=paste(word(Authors, 1), Year, Strain),
                                      at = at_values,
                                      col = c("darkred","darkred"),
                                      addfit = TRUE,
                                      addpred = TRUE,
                                      annotate = TRUE,
                                      header = "Study and Strain",
                                      order=StudyId,
                                      xlab = "",
                                      ilab = cbind(ARRIVEScore, SMD),
                                      ilab.xpos = c(-3, -2),
                                      cex = 0.6, 
                                      cex.axis = 1.0, 
                                      cex.lab = 1.2,
                                      efac = c(1,1,3))
           text(c(-3,-2), model$k+2, c("Reporting\n completeness", "SMD"), cex=0.75, font=2)
         } else {
                               forest(model,
                                      xlim=c((lower_x-4), (upper_x+2)),
                                      mlab="SMD [95% C.I.]",
                                      alim=c((lower_x-4), (upper_x+2)),
                                      slab=paste(word(Authors, 1), Year, Strain),
                                      at = at_values,
                                      col = c("darkred","darkred"),
                                      addfit = TRUE,
                                      addpred = TRUE,
                                      annotate = TRUE,
                                      header = "Study and Strain",
                                      order=StudyId,
                                      xlab = "", 
                                      ilab = cbind(ARRIVEScore, SMD),
                                      ilab.xpos = c(-3, -2),
                                      cex = 0.6, 
                                      cex.axis = 1.0, 
                                      cex.lab = 1.2,
                                      efac = c(1,1,3))
           text(c(-3,-2), model$k+2, c("Reporting\n completeness", "SMD"), cex=0.75, font=2)
         }
         
cixlower <- model[["ci.lb"]]
cixhigher <- model[["ci.ub"]]


  #mtext(outcome_title, side = 1, line = 3, cex = 1.2, font = 2)
  
  if (experiment_type == "TvA") {
    mtext("Favours conventional \nantipsychotic", side = 1, line = 3, at = (lower_x*0.6), cex = 1.1, col = "red", font = 1)
    mtext("Favours TAAR1 \nagonist", side = 1, line = 3, at = (upper_x), cex = 1.1, col = "darkgreen", font = 1)
    #addpoly(model, row = 0.25, cex = 0.4, col = "darkred", mlab = "SMD", annotate = FALSE, xvals = c(cixlower, cixhigher))
    mtext(paste0("SMD: ", round(model$beta, 2), " (", round(model$ci.lb, 2), " to ", round(model$ci.ub, 2), ")"), side = 3, line = -1, cex = 1, font = 2)
    title(paste0("TAAR1 agonists effect on ", outcome_title, " compared with\nconventional antipsychotic in psychosis (SMD)"))
    
  } else if (experiment_type == "AvC") {
    mtext("Favours control", side = 1, line = 3, at = (lower_x*0.6), cex = 1.1, col = "red", font = 1)

    mtext("Favours conventional \nantipsychotic", side = 1, line = 3, at = (upper_x), cex = 1.1, col = "darkgreen", font = 1)

    #addpoly(model, row = 0.25, cex = 0.4, col = "darkred", mlab = "SMD", annotate = FALSE, xvals = c(cixlower, cixhigher))    
    mtext(paste0("SMD: ", round(model$beta, 2), " (", round(model$ci.lb, 2), " to ", round(model$ci.ub, 2), ")"), side = 3, line = -1, cex = 1, font = 2)
    title(paste0("Conventional antipsychotic effect on ", outcome_title, "\nin psychosis (SMD)"))
    
  } else {  
    
    mtext("Favours control", side = 1, line = 3, at = (lower_x*0.7), cex = 1.2, col = "red", font = 1)
    mtext("Favours TAAR1 agonist", side = 1, line = 3, at = (upper_x*0.4), cex = 1.2, col = "darkgreen", font = 1)
  }
}
}

plot_subgroup_analysis <- function(df, experiment_type, outcome, moderator, rho_value) {
  
  # Ensure the moderator is a character string for later conversion to symbol
  moderator <- as.character(moderator)
  
  df2 <- df %>% 
    filter(SortLabel == experiment_type) %>% 
    filter(outcome_type == outcome) %>%  
    filter(!is.na(SMDv)) %>%
    filter(!is.na(!!sym(moderator))) # Filter out NA values in moderator column
  
  # Convert character to factor if necessary
  if (is.character(df2[[moderator]])) {
    df2[[moderator]] <- factor(df2[[moderator]])}
  
  # Add a check for the number of levels in the moderator variable
  if (length(levels(df2[[moderator]])) <= 1) {
    message("In this iteration of the review, there was insufficient data to perform subgroup analysis for this variable (data for one subgroup only)")
    return(NULL)
  }
  
  if ((n_distinct(df$StudyId) > 2) & (n_distinct(df$ExperimentID_I) >10)) {
  #df2$RoBScore <- as.numeric(df2$RoBScore)
  #df2$RoBScore <- factor(df2$RoBScore, levels = c(0, 1, 2))
  
  
  #df2<-df2 %>% 
    #filter(SMD>-6) %>% 
    #filter(SMD<6) # delete missing values and some weirdly large values, like -15 and 16
  
  df2 <- df2 %>% mutate(effect_id = row_number()) # add effect_id column
  
  #calculate variance-covariance matrix of the sampling errors for dependent effect sizes
  
  VCVM_SMD <- vcalc(vi = SMDv,
                    cluster = StudyId, 
                    subgroup= ExperimentID_I,
                    obs=effect_id,
                    data = df2, 
                    rho = rho_value) 
  
  # ML model on df2 with subgroup
  subgroup_analysis <- rma.mv(
    yi = SMD,
    V = VCVM_SMD,
    random = ~1 | Strain / StudyId / ExperimentID_I,
    data = df2,
    mods = as.formula(paste("~", moderator, "-1")),
    method = 'REML',
    test = "t",
    dfs = "contain"
  )
  
  #subgroup_analysis_predict <- predict(subgroup_analysis)
  
  ## ML model on df2 without subgroup
  overall_estimate_rma <- rma.mv(yi = SMD,
                   V = VCVM_SMD,
                   random = ~1 | Strain / StudyId / ExperimentID_I, # nested levels
                   test = "t", # use t- and F-tests for making inferences
                   data = df2,
                   dfs="contain", # improve degree of freedom estimation for t- and F-distributions
                   control=list(optimizer="nlm"))

  #overall_estimate_rma_predict <- predict(overall_estimate_rma)
  
  k_subgroups <- df2 %>%
    group_by(df2[[moderator]]) %>%
    count() %>%
    pull(n)
  
  
  subgroup_analysis_plotdata <- data.frame(levels(df2[[moderator]]), k_subgroups, subgroup_analysis$beta, subgroup_analysis$se) #subgroup_analysis_predict$pi.lb, subgroup_analysis_predict$pi.ub)
  colnames(subgroup_analysis_plotdata) <- c(moderator, "k", "SMD", "se") #, "pi.lb", "pi.ub")
  subgroup_analysis_plotdata=rbind(subgroup_analysis_plotdata, c("Overall estimate",  overall_estimate_rma$k, overall_estimate_rma$beta, overall_estimate_rma$se)) #overall_estimate_rma_predict$pi.lb, overall_estimate_rma_predict$pi.ub))
  
  rownames(subgroup_analysis_plotdata) <- 1:nrow(subgroup_analysis_plotdata)
  subgroup_analysis_plotdata$k <- as.numeric(subgroup_analysis_plotdata$k)
  subgroup_analysis_plotdata$SMD <- as.numeric(subgroup_analysis_plotdata$SMD)
  subgroup_analysis_plotdata$se <- as.numeric(subgroup_analysis_plotdata$se)
  
  overall_estimate_index <-dim(subgroup_analysis_plotdata)[1]
  
  if (moderator == "ARRIVEScoreCat") {
    sorted_data <- subgroup_analysis_plotdata[-overall_estimate_index, ]
    sorted_data <- sorted_data[order(sorted_data[[moderator]]), ]
  } else {
    sorted_data <- subgroup_analysis_plotdata[-overall_estimate_index, ]
  }
  
  options(digits=3)
  
  meta.all = metagen(TE = sorted_data$SMD, 
                     seTE = sorted_data$se, 
                     studlab = sorted_data[[moderator]], 
                     data = sorted_data, 
                     sm = "SMD", 
                     common = F)
  meta.all$TE.random <- subgroup_analysis_plotdata$SMD[overall_estimate_index]
  meta.all$seTE.random <- subgroup_analysis_plotdata$se[overall_estimate_index]

  
  if (moderator == "ARRIVEScoreCat") {
    
    # forest() call without sortvar
    x <- forest(meta.all,
                xlab="SMD",
                smlab=outcome,
                just="right",
                addrow=F,
                overall=T,
                overall.hetstat =F,
                print.pval.Q=F,
                col.square="black",
                col.by="black",
                fill.equi="aliceblue",
                leftcols = c(moderator, "k"),
                leftlabs = c(moderator, "Number of experiments")
    )
  } else {
    # forest() call with sortvar=seTE
    x <- forest(meta.all,
                xlab="SMD",
                smlab=outcome,
                just="right",
                addrow=F,
                overall=T,
                overall.hetstat =F,
                print.pval.Q=F,
                col.square="black",
                sortvar=seTE,
                col.by="black",
                fill.equi="aliceblue",
                leftcols = c(moderator, "k"),
                leftlabs = c(moderator, "Number of experiments")
    )
  }
  
  
  return(list(
    subgroup_analysis = subgroup_analysis,
    subgroup_rma_summary = subgroup_analysis_plotdata))
}}

metaregression_analysis <- function(df, experiment_type, outcome, moderator, rho_value) {
  
  # Ensure the moderator is a character string for evaluation in sym() function (can't convert numerics to symbol)
  moderator <- as.character(moderator)
  
  df2 <- df %>% 
    filter(SortLabel == experiment_type) %>% 
    filter(outcome_type == outcome) %>%  
    filter(!is.na(SMDv)) %>%
    filter(!is.na(!!sym(moderator))) # Filter out NA values in moderator column
  
  # Convert moderator back to numeric
  df2[[moderator]] <- as.numeric(df2[[moderator]])
  
  #df2<-df2 %>% 
    #filter(SMD>-6) %>% 
    #filter(SMD<6) # delete missing values and some weirdly large values, like -15 and 16
  
  df2 <- df2 %>% mutate(effect_id = row_number()) # add effect_id column
  
  #calculate variance-covariance matrix of the sampling errors for dependent effect sizes
  
  VCVM_SMD <- vcalc(vi = SMDv,
                    cluster = StudyId, 
                    subgroup= ExperimentID_I,
                    obs=effect_id,
                    data = df2, 
                    rho = rho_value) 
  
  # Metaregression
  
  metaregression <- rma.mv(
    yi = SMD,
    V = VCVM_SMD,
    random = ~1 | Strain / StudyId / ExperimentID_I,
    data = df2,
    mods = as.formula(paste("~", moderator)),
    method = 'REML',
    test = "t",
    dfs = "contain"
  )
  
  metaregression_summary <- summary(metaregression)
  
  x <- bubble_plot(metaregression, 
              group = "StudyId",
              mod = moderator, 
              xlab = moderator, 
              ylab = "SMD", 
              legend.pos = "none", 
              k = TRUE) 
  
  return(list(
    metaregression_summary = metaregression_summary,
    regression_plot = x))
}

metaregression_analysis_by_drug <- function(df, experiment_type, outcome, drug_name, moderator, rho_value) {
  
  # Ensure the moderator is a character string for evaluation in sym() function (can't convert numerics to symbol)
  moderator <- as.character(moderator)
  
  df2 <- df %>% 
    filter(SortLabel == experiment_type) %>% 
    filter(outcome_type == outcome) %>% 
    filter(DrugName == drug_name) %>% #THIS IS ONLY CHANGE
    filter(!is.na(SMDv)) %>%
    filter(!is.na(!!sym(moderator))) # Filter out NA values in moderator column
  
  if ((n_distinct(df2$StudyId) > 2) & (n_distinct(df2) >10)) {
    
  
  
  # Convert moderator back to numeric
  df2[[moderator]] <- as.numeric(df2[[moderator]])
  
  #df2<-df2 %>% 
  #filter(SMD>-6) %>% 
  #filter(SMD<6) # delete missing values and some weirdly large values, like -15 and 16
  
  df2 <- df2 %>% mutate(effect_id = row_number()) # add effect_id column
  
  #calculate variance-covariance matrix of the sampling errors for dependent effect sizes
  
  VCVM_SMD <- vcalc(vi = SMDv,
                    cluster = StudyId, 
                    subgroup= ExperimentID_I,
                    obs=effect_id,
                    data = df2, 
                    rho = rho_value) 
  
  # Metaregression
  
  metaregression <- rma.mv(
    yi = SMD,
    V = VCVM_SMD,
    random = ~1 | Strain / StudyId / ExperimentID_I,
    data = df2,
    mods = as.formula(paste("~", moderator)),
    method = 'REML',
    test = "t",
    dfs = "contain"
  )
  
  return(metaregression)
  }}
  
metaregression_plot_by_drug <- function(x, df, experiment_type, outcome, moderator, drug_name) {
  
  moderator <- as.character(moderator)
  
  df2 <- df %>% 
    filter(SortLabel == experiment_type) %>% 
    filter(outcome_type == outcome) %>% 
    filter(DrugName == drug_name) %>% #THIS IS ONLY CHANGE
    filter(!is.na(SMDv)) %>%
    filter(!is.na(!!sym(moderator))) # Filter out NA values in moderator column
  
  if ((n_distinct(df2$StudyId) > 2) & (n_distinct(df2) >10)) {
  plot <- bubble_plot(x,
                   group = "StudyId",
                   mod = moderator, 
                   xlab = paste0("Dose of ", drug_name, " (mg/kg)"),
                   ylab = "SMD", 
                   legend.pos = "none", 
                   k = TRUE) 
  
  return(plot)
  }
}


###########################################
### risk of bias visualisation function ###
###########################################

SyRCLE_RoB_summary <- function(df, experiment_type, outcome) {
  
  df <- df %>% 
    filter(SortLabel == experiment_type) %>% 
    filter(outcome_type == outcome) 


RoB <- unique(df[,c(4,6,12,25:58)])

#change studyId to Author, year
RoB$StudyId <- toupper(paste0(str_extract(RoB$Authors,"\\b\\w+\\b"),', ',RoB$Year))

# fix >1 publication per first author in a year

# Assuming your data frame is named RoB and the column is named StudyId
unique_study_ids <- unique(RoB$StudyId)
suffix_list <- character(length = nrow(RoB))

for (study_id in unique_study_ids) {
  indices <- RoB$StudyId == study_id
  if (sum(indices) > 1) {
    suffix_list[indices] <- letters[seq_along(suffix_list[indices])]
  }
}

RoB$suffix <- suffix_list

# Add the suffix to the original column
RoB$StudyId <- paste(RoB$StudyId, RoB$suffix, sep = "")

# Remove the 'suffix' column if you no longer need it
RoB <- select(RoB, -suffix)
RoB <- RoB[order(RoB$StudyId),]


#extract Syrcle RoB scores
SyRCLE <- RoB[,c(1,5:14)]

#Change "yes' to 'low' and 'No' to 'high'
SyRCLE <- mutate_all(SyRCLE, list(~ ifelse(. == 'Yes', 'Low', .)))
SyRCLE <- mutate_all(SyRCLE, list(~ ifelse(. == 'No', 'High', .)))


colnames(SyRCLE) <- c('Study ID','Allocation sequence','Baseline similarity','Concealment of allocation sequence','Random housing','Caregivers blinded','Random selection for outcome assessment','Blinded outcome assessor','Incomplete data reporting addressed','Free from selective outcome reporting','Free of other risks of bias')
RoB_summary <- rob_summary(data <- SyRCLE, tool = "Generic", weighted = FALSE, overall = FALSE)

return(RoB_summary)
}

SyRCLE_RoB_traffic <- function(df, experiment_type, outcome) {
  
  df <- df %>% 
    filter(SortLabel == experiment_type) %>% 
    filter(outcome_type == outcome) 
  

  RoB <- unique(df[,c(4,6,12,25:58)])

  #change studyId to Author, year
  RoB$StudyId <- toupper(paste0(str_extract(RoB$Authors,"\\b\\w+\\b"),', ',RoB$Year))
  
  # fix >1 publication per first author in a year
  
  # Assuming your data frame is named RoB and the column is named StudyId
  unique_study_ids <- unique(RoB$StudyId)
  suffix_list <- character(length = nrow(RoB))
  
  for (study_id in unique_study_ids) {
    indices <- RoB$StudyId == study_id
    if (sum(indices) > 1) {
      suffix_list[indices] <- letters[seq_along(suffix_list[indices])]
    }
  }
  
  RoB$suffix <- suffix_list
  
  # Add the suffix to the original column
  RoB$StudyId <- paste(RoB$StudyId, RoB$suffix, sep = "")
  
  # Remove the 'suffix' column if you no longer need it
  RoB <- select(RoB, -suffix)
  RoB <- RoB[order(RoB$StudyId),]
  
  
  #extract Syrcle RoB scores
  SyRCLE <- RoB[,c(1,5:14)]
  
  #Change "yes' to 'low' and 'No' to 'high'
  SyRCLE <- mutate_all(SyRCLE, list(~ ifelse(. == 'Yes', 'Low', .)))
  SyRCLE <- mutate_all(SyRCLE, list(~ ifelse(. == 'No', 'High', .)))
  
  
  colnames(SyRCLE) <- c('Study ID','Allocation sequence','Baseline similarity','Concealment of allocation sequence','Random housing','Caregivers blinded','Random selection for outcome assessment','Blinded outcome assessor','Incomplete data reporting addressed','Free from selective outcome reporting','Free of other risks of bias')
  RoB_TL <- rob_traffic_light(data <- SyRCLE, tool = "Generic", psize = 10, overall = FALSE)
  
  return(RoB_TL)
}


ARRIVE_summary <- function(df, experiment_type, outcome) {
  
  df <- df %>% 
    filter(SortLabel == experiment_type) %>% 
    filter(outcome_type == outcome)
  

  RoB <- unique(df[,c(4,6,12,25:58)])

  #change studyId to Author, year
  RoB$StudyId <- toupper(paste0(str_extract(RoB$Authors,"\\b\\w+\\b"),', ',RoB$Year))
  
  # fix >1 publication per first author in a year
  
  # Assuming your data frame is named RoB and the column is named StudyId
  unique_study_ids <- unique(RoB$StudyId)
  suffix_list <- character(length = nrow(RoB))
  
  for (study_id in unique_study_ids) {
    indices <- RoB$StudyId == study_id
    if (sum(indices) > 1) {
      suffix_list[indices] <- letters[seq_along(suffix_list[indices])]
    }
  }
  
  RoB$suffix <- suffix_list
  
  # Add the suffix to the original column
  RoB$StudyId <- paste(RoB$StudyId, RoB$suffix, sep = "")
  
  # Remove the 'suffix' column if you no longer need it
  RoB <- select(RoB, -suffix)
  RoB <- RoB[order(RoB$StudyId),]

  #extract ARRIVE reporting scores
ARRIVE <- RoB[,c(1,15:37)]

#Change "yes' to 'low' and 'No' to 'high'
ARRIVE <- mutate_all(ARRIVE, list(~ ifelse(. == 'Yes', 'Low', .)))
ARRIVE <- mutate_all(ARRIVE, list(~ ifelse(. == 'No', 'High', .)))
#but ethics NAs to justification into ethics low risk
ARRIVE <- mutate_all(ARRIVE, list(~ ifelse(. == 'NA (ethical approval declared)', 'Low', .)))
#combine desc stats and variance with ES and CI
ARRIVE <- ARRIVE %>%
  mutate(Data_reporting = ifelse(ARRIVE$`(ARRIVE) Are desc stats for each exp group provided with measure of variability?_I` == 'Low' | ARRIVE$`(ARRIVE) Is the effect size and confidence interval provided?_I` == 'Low', 'Low', 'High'))
ARRIVE <- ARRIVE[,c(1:17,25,20:24)]


colnames(ARRIVE) <- c('Study','Groups clearly defined','Experimental unit defined','Exact number of experimental units','Sample size justification',
                      'Inclusion and exclusion criteria given','Any exclusions reported','Randomisation for any experiments','Blinding to group allocation',
                      'Details of what was measured','Statistical approach for each outcome','Assessment of whether data met statistical assumptions',
                      'All species specified','Animal sex specified','Age, weight or developmental stage specified','Timing and frequency of proceedures described',
                      'Any acclimitisation described','Data with variance, or Effect size and CI','Ethical approval with approval number',
                      'Ethical approval with or without approval number','Conflicts of interest statement','Funding sources','Description of any role of funder')
Rep_summary <- rob_summary(data <- ARRIVE, tool = "Generic", weighted = FALSE, overall = FALSE)

return(Rep_summary)
}

ARRIVE_traffic <- function(df, experiment_type, outcome) {
  
  dfa <- subset(df, df$SortLabel == experiment_type)
  dfb <- subset(dfa, dfa$outcome_type == outcome)


  RoB <- unique(dfb[,c(4,6,12,25:58)])

  #change studyId to Author, year
  RoB$StudyId <- toupper(paste0(str_extract(RoB$Authors,"\\b\\w+\\b"),', ',RoB$Year))
  
  # fix >1 publication per first author in a year
  
  # Assuming your data frame is named RoB and the column is named StudyId
  unique_study_ids <- unique(RoB$StudyId)
  suffix_list <- character(length = nrow(RoB))
  
  for (study_id in unique_study_ids) {
    indices <- RoB$StudyId == study_id
    if (sum(indices) > 1) {
      suffix_list[indices] <- letters[seq_along(suffix_list[indices])]
    }
  }
  
  RoB$suffix <- suffix_list
  
  # Add the suffix to the original column
  RoB$StudyId <- paste(RoB$StudyId, RoB$suffix, sep = "")
  
  # Remove the 'suffix' column if you no longer need it
  RoB <- select(RoB, -suffix)
  RoB <- RoB[order(RoB$StudyId),]
  
  #extract ARRIVE reporting scores
  ARRIVE <- RoB[,c(1,15:37)]
  
  #Change "yes' to 'low' and 'No' to 'high'
  ARRIVE <- mutate_all(ARRIVE, list(~ ifelse(. == 'Yes', 'Low', .)))
  ARRIVE <- mutate_all(ARRIVE, list(~ ifelse(. == 'No', 'High', .)))
  #but ethics NAs to justification into ethics low risk
  ARRIVE <- mutate_all(ARRIVE, list(~ ifelse(. == 'NA (ethical approval declared)', 'Low', .)))
  #combine desc stats and variance with ES and CI
  ARRIVE <- ARRIVE %>%
    mutate(Data_reporting = ifelse(ARRIVE$`(ARRIVE) Are desc stats for each exp group provided with measure of variability?_I` == 'Low' | ARRIVE$`(ARRIVE) Is the effect size and confidence interval provided?_I` == 'Low', 'Low', 'High'))
  ARRIVE <- ARRIVE[,c(1:17,25,20:24)]
  
  
  colnames(ARRIVE) <- c('Study','Groups clearly defined','Experimental unit defined','Exact number of experimental units','Sample size justification',
                        'Inclusion and exclusion criteria given','Any exclusions reported','Randomisation for any experiments','Blinding to group allocation',
                        'Details of what was measured','Statistical approach for each outcome','Assessment of whether data met statistical assumptions',
                        'All species specified','Animal sex specified','Age, weight or developmental stage specified','Timing and frequency of proceedures described',
                        'Any acclimitisation described','Data with variance, or Effect size and CI','Ethical approval with approval number',
                        'Ethical approval with or without approval number','Conflicts of interest statement','Funding sources','Description of any role of funder')
  Rep_TL <- rob_traffic_light(data = ARRIVE, tool = "Generic", psize = 10, overall = FALSE)
  
  return(Rep_TL)
}

run_sse_NMD <- function(df, rho_value = 0.5) {
  
  #  df<-filter_experiment_outcome_type(df, experiment, outcome)
  
  df<-df %>% 
    filter(!is.na(NMDv)) %>%
    filter(outcome_type == "Locomotor activity") %>%
    filter(SortLabel == "TvC")
  
  df <- df %>% mutate(effect_id = row_number()) # add effect_id column
  df$NMDSE <- sqrt(df$NMDv)
  
  #calculate variance-covariance matrix of the sampling errors for dependent effect sizes
  
  
  VCVM_NMD <- vcalc(vi = NMDv,
                    cluster = StudyId, 
                    subgroup= ExperimentID_I,
                    obs=effect_id,
                    data = df, 
                    rho = rho_value)
  
  NMD_sse <- rma.mv(yi = NMD,
                    V = VCVM_NMD,
                    random = ~1 | Strain / StudyId / ExperimentID_I, # nested levels
                    mods = ~ NMDSE, # sampling error (squart root of sampling variance SMDV);
                    test = "t", # use t- and F-tests for making inferences
                    data = df,
                    dfs="contain",
                    control=list(optimizer="nlm")
  )
  
  return(NMD_sse)
}

run_sse_plot <- function(df, rho_value = 0.5) {
  
  #  df<-filter_experiment_outcome_type(df, experiment, outcome)
  
  df<-df %>% 
    filter(!is.na(NMDv)) %>%
    filter(outcome_type == "Locomotor activity") %>%
    filter(SortLabel == "TvC")
  
  df <- df %>% mutate(effect_id = row_number()) # add effect_id column
  df$NMDSE <- sqrt(df$NMDv)
  
  #calculate variance-covariance matrix of the sampling errors for dependent effect sizes
  
  
  VCVM_NMD <- vcalc(vi = NMDv,
                    cluster = StudyId, 
                    subgroup= ExperimentID_I,
                    obs=effect_id,
                    data = df, 
                    rho = rho_value)
  
  NMD_sse <- rma.mv(yi = NMD,
                    V = VCVM_NMD,
                    random = ~1 | Strain / StudyId / ExperimentID_I, # nested levels
                    mods = ~ NMDSE, # sampling error (squart root of sampling variance SMDV);
                    test = "t", # use t- and F-tests for making inferences
                    data = df,
                    dfs="contain",
                    control=list(optimizer="nlm")
  )
  
  plot <- bubble_plot(NMD_sse, mod = "NMDSE", group = "StudyId", xlab = "SE of NMD estimate", ylab = "NMD estimate")
  return(plot)
}

############################################################
run_ML_NMD <- function(df, experiment, outcome, rho_value) {
  
  df<-filter_experiment_outcome_type(df, experiment, outcome) 
  
  #df<-df %>% 
    #filter(!is.na(NMDv)) %>% 
    #filter(NMD>-600) %>% 
    #filter(NMD<600) # delete missing values and some weirdly large values, like -15 and 16
  
  df <- df %>% mutate(effect_id = row_number()) # add effect_id column
  
  #calculate variance-covariance matrix of the sampling errors for dependent effect sizes
  
  VCVM_NMD <- vcalc(vi = NMDv,
                    cluster = StudyId, 
                    subgroup= ExperimentID_I,
                    obs=effect_id,
                    data = df, 
                    rho = rho_value) 
  
  NMD_ML <- rma.mv(yi = NMD,
                   V = VCVM_NMD,
                   random = ~1 | Strain / StudyId / ExperimentID_I, # nested levels
                   test = "t", # use t- and F-tests for making inferences
                   data = df,
                   dfs="contain",
                   control=list(optimizer="nlm")
  )
  
  #if (length(unique(df$StudyId)) > 1) {
  #NMD_ML <- robust(NMD_ML, cluster = StudyId, clubSandwich = FALSE)
  #}
  
  cat("Meta analysis summary:\n")
  print(summary(NMD_ML))
  
  cat("\n-------------------------\n")
  
  cat("Prediction Interval:\n")
  
  pred_interval <- predict(NMD_ML)
  
  print(pred_interval)
  
  return(NMD_ML)
}

forest_metafor_NMD <- function(model, outcome){
  
  lower_x <- min(model[["yi"]])-mean(model[["vi"]])
  upper_x <- max(model[["yi"]])+mean(model[["vi"]])
  summary_x <- model[["beta"]]
  
  at_values <- seq(floor(lower_x / 5) * 5, ceiling(upper_x / 5) * 5, by = 5)
  
  pred_interval <- predict(model)
  
  ifelse(model[["k"]] > 25, 
         forest_plot <- forest(model, 
                               xlim=c(lower_x-20, upper_x+20),
                               mlab="NMD",
                               slab=NA,
                               alim=c(lower_x-20, upper_x+20),
                               at = at_values,
                               col = c("darkred","darkred"),
                               addfit = TRUE,
                               addpred = TRUE,
                               annotate = FALSE,
                               order="obs",
                               xlab = "", 
                               cex = 1.0, 
                               cex.axis = 1.0, 
                               cex.lab = 1.2,
                               efac = c(1,1,3)), 
         forest_plot <- forest(model, 
                               xlim=c(lower_x-20, upper_x+20),
                               mlab="NMD",
                               slab=NA,
                               alim=c(lower_x-20, upper_x+20),
                               at = at_values,
                               col = c("darkred","darkred"),
                               addfit = TRUE,
                               addpred = TRUE,
                               annotate = TRUE,
                               order="obs",
                               xlab = ""))
  
  #mtext(outcome, side = 1, line = 3, cex = 1.2, font = 2)
  mtext("Favours control", side = 1, line = 3, at = (lower_x + 30), cex = 1.2, col = "red", font = 1)
  mtext("Favours TAAR1 agonist", side = 1, line = 3, at = (upper_x - 40), cex = 1.2, col = "darkgreen", font = 1)
  mtext(paste0("NMD: ", round(model$beta, 2), " (", round(model$ci.lb, 2), " to ", round(model$ci.ub, 2), ")"), side = 2, line = 3, cex = 1.2, font = 2)
  title(paste0("TAAR1 agonist effect on ", outcome, " in psychosis (NMD)"))
  
}


###### Function to check number of levels in moderator variables for an experiment type (SortLabel) and outcome for this iteration #######

check_moderator_levels <- function(df, experiment, outcome) {
  
  moderator_vars <- c("Sex", "CategoryDiseaseInduction", "InterventionAdministrationRoute", 
                      "ProphylacticOrTherapeutic", "TreatmentDurationCategory", 
                      "DrugName", "Efficacy", "Selectivity")

  
  df2 <- filter_experiment_outcome_type(df, experiment, outcome)
  single_level_mods <- character()  # Initialize an empty character vector to store moderators
  
  for (moderator in moderator_vars) {
    # convert to factor if not already
    if (!is.factor(df2[[moderator]])) {
      df2[[moderator]] <- factor(df2[[moderator]])
    }
    
    # check number of levels
    moderator_levels <- df2 %>%
      group_by(!!sym(moderator)) %>%
      summarise(n = n_distinct(StudyId)) %>%
      ungroup() %>%
      summarise(moderator_levels = n_distinct(!!sym(moderator))) %>%
      pull(moderator_levels)
    
    if (moderator_levels < 2) {
      single_level_mods <- c(single_level_mods, moderator)
    }
  }
  
  # rename moderator variables in list to match names in the inline text
  
  single_level_mods1 <- single_level_mods %>% 
    str_replace("TreatmentDurationCategory", "Duration of treatment period") %>% 
    str_replace("InterventionAdministrationRoute", "Route of intervention administration") %>%
    str_replace("ProphylacticOrTherapeutic", "Prophylactic or therapeutic intervention") %>%
    str_replace("CategoryDiseaseInduction", "Disease induction method") %>%
    str_replace("DrugName", "Intervention admnistered (drug)") %>% 
    tolower()
  
  if (length(single_level_mods1) > 1) {
    x <- paste(head(single_level_mods1, -1), collapse = ", ")
    x <- paste(x, "and", tail(single_level_mods1, 1))
  } else {
    x <- single_level_mods1
  }
  
  return(x)  # Return the formatted string
}

run_sse_SMD <- function(df, rho_value = 0.5) {
  
  #  df<-filter_experiment_outcome_type(df, experiment, outcome)
  
  df<-df %>% 
    filter(!is.na(SMDv)) %>%
    filter(outcome_type == "Locomotor activity") %>%
    filter(SortLabel == "TvC")
  
  df <- df %>% mutate(effect_id = row_number()) # add effect_id column
  df$SMDN <- 1/sqrt(as.numeric(df$NumberOfAnimals_C) + as.numeric(df$NumberOfAnimals_I))
  
  #calculate variance-covariance matrix of the sampling errors for dependent effect sizes
  
  
  VCVM_SMD <- vcalc(vi = SMDv,
                    cluster = StudyId, 
                    subgroup= ExperimentID_I,
                    obs=effect_id,
                    data = df, 
                    rho = rho_value)
  
  SMD_sse <- rma.mv(yi = SMD,
                    V = VCVM_SMD,
                    random = ~1 | Strain / StudyId / ExperimentID_I, # nested levels
                    mods = ~ SMDN, # sampling error (squart root of N);
                    test = "t", # use t- and F-tests for making inferences
                    data = df,
                    dfs="contain",
                    control=list(optimizer="nlm")
  )
  
  return(SMD_sse)
}

run_sse_plot_SMD <- function(df, rho_value = 0.5) {
  
  #  df<-filter_experiment_outcome_type(df, experiment, outcome)
  
  df<-df %>% 
    filter(!is.na(SMDv)) %>%
    filter(outcome_type == "Locomotor activity") %>%
    filter(SortLabel == "TvC")
  
  df <- df %>% mutate(effect_id = row_number()) # add effect_id column
  df$SMDN <- 1/sqrt(as.numeric(df$NumberOfAnimals_C) + as.numeric(df$NumberOfAnimals_I))
  
  #calculate variance-covariance matrix of the sampling errors for dependent effect sizes
  
  
  VCVM_SMD <- vcalc(vi = SMDv,
                    cluster = StudyId, 
                    subgroup= ExperimentID_I,
                    obs=effect_id,
                    data = df, 
                    rho = rho_value)
  
  SMD_sse <- rma.mv(yi = SMD,
                    V = VCVM_SMD,
                    random = ~1 | Strain / StudyId / ExperimentID_I, # nested levels
                    mods = ~ SMDN, # sampling error (squart root of N);
                    test = "t", # use t- and F-tests for making inferences
                    data = df,
                    dfs="contain",
                    control=list(optimizer="nlm")
  )
  
  plot <- bubble_plot(SMD_sse, mod = "SMDN", group = "StudyId", xlab = "1/SQRT(N) associated with SMD estimate", ylab = "SMD estimate", legend.pos = "none")
  return(plot)
}





