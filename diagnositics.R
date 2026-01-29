data <- LSR3_clean_data_2026_01_29
data_TvC <- subset(data, data$SortLabel == "TvC")


data_TVC_Loc <- subset(data_TvC, data_TvC$outcome_type == "Prepulse inhibition")
nr <- nrow(data_TVC_Loc)
data_T_L_A <- data.frame(matrix(nrow=nr))


data_T_L_A$ID <- data_TVC_Loc$StudyId 
data_T_L_A$AS <- data_TVC_Loc$ARRIVEScore
data_T_L_A$Sy <- data_TVC_Loc$RoBTF
data_T_L_A <- data_T_L_A[,c(2,3,4)]
d_t_l_a <- unique(data_T_L_A)

med <- median(d_t_l_a$AS)
