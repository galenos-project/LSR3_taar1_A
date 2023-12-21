library(dplyr)
library(stringr)
library(tibble)
library(tidyverse)
library(tidyr)
library(purrr)
library(metafor)
#setwd("Francesca_analysis")
source("wrangling/wrangling_functions.R", local = TRUE)

LSR <- 'LSR3'

# Import SyRF outcome data
LSR3_SyRFOutcomes <- read_csv("data/Quantitative_data_-_2023_12_18_-_c494b1ae-4cf4-4618-b91a-e69a2b815bdd_-_Investigators_Unblinded.csv")

#clean ; from TiAb etc
LSR3_SyRFOutcomes$Title <- gsub(";", ":", LSR3_SyRFOutcomes$Title)
LSR3_SyRFOutcomes$Abstract <- gsub(";", ":", LSR3_SyRFOutcomes$Abstract)
LSR3_SyRFOutcomes$Authors <- gsub(";", ":", LSR3_SyRFOutcomes$Authors)
LSR3_SyRFOutcomes$PublicationName <- gsub(";", ":", LSR3_SyRFOutcomes$PublicationName)
LSR3_SyRFOutcomes$AlternateName <- gsub(";", ":", LSR3_SyRFOutcomes$AlternateName)
LSR3_SyRFOutcomes$Url <- gsub(";", ":", LSR3_SyRFOutcomes$Url)
LSR3_SyRFOutcomes$AuthorAddress <- gsub(";", ":", LSR3_SyRFOutcomes$AuthorAddress)
LSR3_SyRFOutcomes$Doi <- gsub(";", ":", LSR3_SyRFOutcomes$Doi)
LSR3_SyRFOutcomes$Keywords <- gsub(";", ":", LSR3_SyRFOutcomes$Keywords)
LSR3_SyRFOutcomes$CustomId <- gsub(";", ":", LSR3_SyRFOutcomes$CustomId)
LSR3_SyRFOutcomes$Year<- gsub(";", ":", LSR3_SyRFOutcomes$Year)
LSR3_SyRFOutcomes$ReferenceType <- gsub(";", ":", LSR3_SyRFOutcomes$ReferenceType)
LSR3_SyRFOutcomes$PdfRelativePath <- gsub(";", ":", LSR3_SyRFOutcomes$PdfRelativePath)


# Filter for reconciled studies (and rename columns for consistency with shiny outcomes/remove SYRF columnns)
LSR3_reconciled <- LSR3_SyRFOutcomes %>% 
  filter(`Is this a reconciliation?` == TRUE) %>%  #for PTSD `is this a reconciliation?`
  rename(ModelID = `DiseaseModelId(s)`, 
         ExperimentID = ExperimentId, 
         InterventionLabel = `TreatmentLabel(s)`,
         InterventionID = `TreatmentId(s)`,
         OutcomeResult = OutcomeAverage,
         `Dose of pharmacological induction (dose)` = `Dose of pharmacological induction...61`,
         `Dose of pharmacological induction (units)` = `Dose of pharmacological induction...62`)

## first pass for removing duplicate reconciliations
  recent_reconciledID <- LSR3_reconciled %>%
  arrange(desc(DateTimeDataEntered)) %>%
  group_by(StudyId) %>%
  slice(1) %>%
  ungroup()

# Filter all records that have a StudyIdStr and AnnotatorIdStr that is in the reconciled_study_annotator_pairs, to get all reconciled data
#reconciled_records_unique <- LSR3_reconciled %>%
#  semi_join(recent_reconciledID, by = "StudyId")

## Some studies were split due to their complexity, with RoB/ Arrive only entered once
## the ARRIVE/RoB data for the second split needs overwritten with the values from the first
#if this is the case, identify the studyIds concerned, recort in ll48 and 49, and change ll52-59 accordingly

## Match split studies [Get pairs from bibliographic download on SyRF]
## Find position of the start and end of ARRIVE ROB columns to store column names between the start and end index
col_names <- names(LSR3_reconciled)
start_index <- match("Title", col_names)
end_index <- match("Is any role of the funder in the design/analysis/reporting of study described?", col_names)
ARRIVEROB_columns <- col_names[start_index:end_index] 





## 1. Overwrite appendix ARRIVE/ROB with main ARRIVE/ROB for both pairs
## PAIR 1: fb8ed201-f663-48db-aae1-8ee88b355abd - main, dd9dddac-739b-412f-9fac-6f36e0a21494 - appendix
r_to_r <- subset(LSR3_reconciled, LSR3_reconciled$StudyId == 'dd9dddac-739b-412f-9fac-6f36e0a21494')
source_row <- subset(LSR3_reconciled, LSR3_reconciled$StudyId == 'fb8ed201-f663-48db-aae1-8ee88b355abd') %>% slice(1)
for(i in start_index: end_index) {
  for(j in 1:nrow(r_to_r)) {
    r_to_r[j,i] <- source_row[1,i]
  }
}
residual_rows <- subset(LSR3_reconciled, !LSR3_reconciled$StudyId %in% r_to_r$StudyId)
LSR3_reconciled <- rbind(residual_rows, r_to_r)

## PAIR 2: 4f31dcd6-e041-4882-acc8-dbf3ccfd2368 - main, 283c541f-6f14-462b-9f3a-800b69a3f44c - appendix
r_to_r <- subset(LSR3_reconciled, LSR3_reconciled$StudyId == '283c541f-6f14-462b-9f3a-800b69a3f44c')
source_row <- subset(LSR3_reconciled, LSR3_reconciled$StudyId == '4f31dcd6-e041-4882-acc8-dbf3ccfd2368') %>% slice(1)
for(i in start_index: end_index) {
  for(j in 1:nrow(r_to_r)) {
    r_to_r[j,i] <- source_row[1,i]
  }
}
residual_rows <- subset(LSR3_reconciled, !LSR3_reconciled$StudyId %in% r_to_r$StudyId)
LSR3_reconciled <- rbind(residual_rows, r_to_r)

## 2. Replace StudyId of appendix studies with the StudyId of corresponding main paper
LSR3_reconciled <- LSR3_reconciled %>% 
  mutate(StudyId = ifelse(StudyId == 'dd9dddac-739b-412f-9fac-6f36e0a21494', 'fb8ed201-f663-48db-aae1-8ee88b355abd', StudyId)) %>% 
  mutate(StudyId = ifelse(StudyId == '283c541f-6f14-462b-9f3a-800b69a3f44c', '4f31dcd6-e041-4882-acc8-dbf3ccfd2368', StudyId))

## 3. Further step to remove observations from an accidental dual-reconciliation: choose the most recent reconciliation (unlikely to be necessary)
# Choose the most recent reconciliation 
recent_reconciledID <- LSR3_reconciled %>%
  arrange(desc(DateTimeDataEntered)) %>%
  group_by(StudyId) %>%
  slice(1) %>%
  ungroup()

# Filter all records that have a StudyIdStr and AnnotatorIdStr that is in the reconciled_study_annotator_pairs to get final set of reconciled data
reconciled_records_unique <- LSR3_reconciled %>%
  semi_join(recent_reconciledID, by = "StudyId")

# Rearrange rows to be grouped by studies and then experiment within studies. Reorder columns for readability
reconciled_records <- reconciled_records_unique %>%
  arrange(StudyId, ExperimentID, CohortId, OutcomeId, InterventionID) %>% 
  rename(`Dose of positive control treatment?` = `Dose of control treatment?`) %>% 
  relocate(ExperimentLabel, .after = ExperimentID) %>% 
  relocate(CohortLabel, .after = CohortId) %>% 
  relocate(OutcomeLabel, .after = OutcomeId) %>% 
  relocate(c(OutcomeLabel, OutcomeId), .after = ExperimentLabel) %>% 
  relocate(c(InterventionLabel), .after = InterventionID)

# this segmented removed by MM 151223
# Remove ARRIVE/ROB columns and put into separate dataframe (-> reconciled_record_ROB). Join later
# Now two dataframes, one which contains information that is the same for every observation within a study (reconciled_studyconstants), and one which contains information that differs across every observation (reconciled_studyvaried)
# reonciled_record_ROB <- reconciled_records %>% 
#  select(StudyId, ExperimentID, OutcomeId, any_of(ARRIVEROB_columns))
#reconciled_records_noROB <- reconciled_records %>% 
#  select(-any_of(ARRIVEROB_columns))

savename <- paste0('reconciled_records_',Sys.Date(),'.csv')
write.csv(reconciled_records, savename)

# Split any column which contains multiple responses (separated by ';') into seperate columns
reconciled_records <- split_columns(reconciled_records)

# Label types of cohorts
# Aim = Make each experiment one comparison (with sham where present, control and intervention)
## Differentiate between positive and negative controls and TAAR1KO's

reconciled_cohort_label <- reconciled_records %>%
  mutate(
    Treatment1Type = case_when(
      str_detect(`InterventionLabel[1]`, regex("veh|con", ignore_case = TRUE)) ~ "Negative control",
      str_detect(`InterventionLabel[1]`, "TAAR1 KO") ~ "TAAR1KO",
      str_detect(`InterventionLabel[1]`, regex("cloza|risp|ari|olan|olz", ignore_case = TRUE)) & !str_detect(`InterventionLabel[1]`, regex("ro", ignore_case = TRUE)) ~ "Positive control",
      is.na(`InterventionLabel[1]`) ~ NA_character_,
      TRUE ~ "Intervention"
    ),
    Treatment2Type = case_when(
      str_detect(`InterventionLabel[2]`, regex("veh|con", ignore_case = TRUE)) ~ "Negative control",
      str_detect(`InterventionLabel[2]`, "TAAR1 KO") ~ "TAAR1KO",
      str_detect(`InterventionLabel[2]`, regex("cloza|risp|ari|olan|olz", ignore_case = TRUE)) & !str_detect(`InterventionLabel[2]`, regex("ro", ignore_case = TRUE)) ~ "Positive control",
      is.na(`InterventionLabel[2]`) ~ NA_character_,
      TRUE ~ "Intervention"
    )
  )

## sort the blank disease models (which are equivalent to sham)
reconciled_cohort_label <- reconciled_cohort_label %>%
  mutate(
    IsDiseaseModelControl = case_when(
      IsDiseaseModelControl == 'TRUE' ~ TRUE,
      IsDiseaseModelControl == 'FALSE' ~ FALSE,
      is.na(IsDiseaseModelControl) ~ TRUE,
      TRUE ~ as.logical(NA)  # default case if none of the conditions are met
    )
  )



## Make CohortType column
# Combination interventions are interventions where currently licensed antipsychotic is an intervention

reconciled_cohort_type <- reconciled_cohort_label %>%
  mutate(CohortType = case_when(
    IsDiseaseModelControl == TRUE & 
      ((Treatment1Type == "TAAR1KO" & Treatment2Type == "Negative control") |
         (Treatment2Type == "TAAR1KO" & Treatment1Type == "Negative control")) ~ "Sham for TAAR1KO comb.",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "TAAR1KO" & Treatment2Type == "Negative control") |
         (Treatment2Type == "TAAR1KO" & Treatment1Type == "Negative control")) ~ "Negative control for TAAR1KO comb.",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "TAAR1KO" & Treatment2Type == "Positive control") |
         (Treatment2Type == "TAAR1KO" & Treatment1Type == "Positive control")) ~ "Positive control for TAAR1KO comb.",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "TAAR1KO" & Treatment2Type == "Intervention") |
         (Treatment2Type == "TAAR1KO" & Treatment1Type == "Intervention")) ~ "Intervention and TAAR1KO comb.",
    
    (IsDiseaseModelControl == TRUE | is.na(IsDiseaseModelControl)) & 
      (Treatment1Type %in% c("TAAR1KO", "Negative control", NA) & 
         Treatment2Type %in% c("TAAR1KO", "Negative control", NA)) ~ "Sham for simple intervention",
    
    (IsDiseaseModelControl == TRUE | is.na(IsDiseaseModelControl)) & 
      (Treatment1Type == "Positive control" | Treatment2Type == "Positive control") ~ "Sham for combination intervention",
    
    (IsDiseaseModelControl == TRUE | is.na(IsDiseaseModelControl)) & 
      (Treatment1Type == "Positive control" & (Treatment2Type != "Intervention" | is.na(Treatment2Type)) | 
         (Treatment2Type == "Positive control" & (Treatment1Type != "Intervention" | is.na(Treatment1Type)))) ~ "Positive control treated sham",
    
    (IsDiseaseModelControl == TRUE | is.na(IsDiseaseModelControl)) & 
      (Treatment1Type == "Intervention" | Treatment2Type == "Intervention") ~ "Intervention treated sham",
    
    IsDiseaseModelControl == FALSE & 
      (Treatment1Type == "Positive control" & (Treatment2Type != "Intervention" | is.na(Treatment2Type)) | 
         (Treatment2Type == "Positive control" & (Treatment1Type != "Intervention" | is.na(Treatment1Type)))) ~ "Positive control",
    
    IsDiseaseModelControl == FALSE & 
      (Treatment1Type != "Intervention" | is.na(Treatment1Type)) & 
      (Treatment2Type != "Intervention" | is.na(Treatment2Type)) ~ "Negative control for simple intervention",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "Intervention" & Treatment2Type == "Positive control") | 
         (Treatment2Type == "Intervention" & Treatment1Type == "Positive control")) ~ "Combination intervention",
    
    IsDiseaseModelControl == FALSE & 
      (Treatment1Type == "Intervention" & is.na(Treatment2Type)) | 
      (Treatment2Type == "Intervention" & is.na(Treatment1Type)) ~ "Simple intervention"
  )) %>% 
  relocate(CohortType, .after = `InterventionID[2]`) %>% 
  relocate(c(Treatment1Type, Treatment2Type), .after = `InterventionLabel[2]`) %>% 
  relocate(IsDiseaseModelControl, .after = ModelID)


## Label how positive controls are being used
# If for an observation with the value "Positive control" in the CohortType column, there is an observation with the CohortType "Combination intervention" with the same combination of values in the StudyId, ExperimentID, OutcomeId, TimeInMinute variables, then the value in the CohortType column needs changed to "Control for combination intervention". Similarly, if for an observation with the value "Positive control treated sham" in the CohortType column, there is an observation with the CohortType "Combination intervention" with the same combination of values in the StudyId, ExperimentID, OutcomeId, TimeInMinute variables, then the value in the CohortType column needs changed to "Sham for combination intervention"

reconciled_cohort_type <- positive_control_cohort_type(reconciled_cohort_type)

############## FIX THIS SECTION IF EXCLUDING SOME POSITIVE CONTROLS #####################
#only interested in positive controls that serve as a comparison to an intervention

#test <- reconciled_cohort_type %>%  
 # filter(CohortType != "Control for combination intervention") %>% # remove positive controls that aren't for a simple intervention
  # filter(CohortType != "Sham for combination intervention")

#test <- reconciled_cohort_type %>%  
  #group_by(GroupID) %>%   
  #filter(CohortType == "Positive control" & any(CohortType == "Simple intervention")) %>%
  #ungroup()

#########################################################################################


## Identify experiments that are just 'Simple' comparison, those which are 'Combination' comparisons and those which have 'TAAR1KO' involved
reconciled_comparison_type <- reconciled_cohort_type %>%
  mutate(ExperimentType = case_when(
    str_detect(CohortType, "combination") | str_detect(CohortType, "Combination") ~ "Antipsychotic combination",
    str_detect(CohortType, "TAAR1")  ~ "TAAR1KO combination",
    str_detect(CohortType, "simple") | str_detect(CohortType, "Simple") ~ "Simple intervention"
  ))

reconciled_cohort_role <- reconciled_comparison_type %>% 
  mutate(RoleOfCohort = case_when(
    str_detect(CohortType, "Sham") | str_detect(CohortType, "sham") ~ "S", 
    str_detect(CohortType, "Control") | str_detect(CohortType, "control") ~ "C",
    str_detect(CohortType, "Intervention") | str_detect(CohortType, "intervention") ~ "I"
  ))

## Add placeholder sham rows if not present: MM not sure why you would do this
#  reconciled_with_shams <- add_missing_sham(reconciled_cohort_role)

## Wrangle wide so each observations is a single comparison
data <- reconciled_cohort_role
# remove DA knockouts
data_rem <- subset(data, data$`DiseaseModelLabel(s)` == 'DAT +/-')
data <- anti_join(data, data_rem)
#remove data reported from subgroups
data_rem <- data %>% filter(str_detect(`DiseaseModelLabel(s)`, "SUBGROUP") | str_detect(CohortLabel, 'SUBGROUP'))
data <- anti_join(data, data_rem)

groups <- as.data.frame(unique(data$GroupID))
groups$`unique(data$GroupID)` <- as.character(groups$`unique(data$GroupID)`)

###run 1  - moving relevant control and sham to another column###

group <- groups[1,1]
cohset <- subset(data, data$GroupID == group)

chhrows <- nrow(cohset)

lab_n <- match('CohortLabel', names(data))
n_n <- match('NumberOfAnimals', names(data))
m_n <- match('OutcomeResult', names(data))
v_n <- match('OutcomeError', names(data))

for (j in 1:chhrows) {
  if (cohset$CohortType[j] == 'Negative control for simple intervention') {
    cohset$Control_label <- cohset[j, lab_n]
    cohset$Control_n <- cohset[j, n_n]
    cohset$Control_m <- cohset[j, m_n]
    cohset$Control_v <- cohset[j, v_n]
  }
}

for (j in 1:chhrows) {
  if(cohset$CohortType[j] == 'Sham for simple intervention'){
    cohset$Sham_label <- cohset[j,lab_n]
    cohset$Sham_n <- cohset[j,n_n]
    cohset$Sham_m <- cohset[j,m_n]
    cohset$Sham_v <- cohset[j,v_n]
  }
}
for (j in 1:chhrows) {
  if(cohset$CohortType[j] == 'Control for combination intervention'){
    cohset$CombC_label <- cohset[j,lab_n]
    cohset$CombC_n <- cohset[j,n_n]
    cohset$CombC_m <- cohset[j,m_n]
    cohset$CombC_v <- cohset[j,v_n]
  }
}
for (j in 1:chhrows) {
  if(cohset$CohortType[j] == 'Sham for TAAR1KO comb.'){
    cohset$ShamTAAR_label <- cohset[j,lab_n]
    cohset$ShamTAAR_n <- cohset[j,n_n]
    cohset$ShamTAAR_m <- cohset[j,m_n]
    cohset$ShamTAAR_v <- cohset[j,v_n]
  }
}
for (j in 1:chhrows) {
  if(cohset$CohortType[j] == 'Negative control for TAAR1KO comb.'){
    cohset$CTAAR_label <- cohset[j,lab_n]
    cohset$CTAAR_n <- cohset[j,n_n]
    cohset$CTAAR_m <- cohset[j,m_n]
    cohset$CTAAR_v <- cohset[j,v_n]
  }
}

data_org <- cohset

### now run the rest###
for(i in 2:nrow(groups)){
  group <- groups[i,1]
  cohset <- subset(data, data$GroupID == group)
  #chhrows <- nrow(cohset)
  # Assuming chhrows is the number of rows in cohset
  chhrows <- nrow(cohset)
  
  for (j in 1:chhrows) {
    if (cohset$CohortType[j] == 'Negative control for simple intervention') {
      cohset$Control_label <- cohset[j, lab_n]
      cohset$Control_n <- cohset[j, n_n]
      cohset$Control_m <- cohset[j, m_n]
      cohset$Control_v <- cohset[j, v_n]
    }
  }
  
  for (j in 1:chhrows) {
    if(cohset$CohortType[j] == 'Sham for simple intervention'){
      cohset$Sham_label <- cohset[j,lab_n]
      cohset$Sham_n <- cohset[j,n_n]
      cohset$Sham_m <- cohset[j,m_n]
      cohset$Sham_v <- cohset[j,v_n]
    }
  }
  for (j in 1:chhrows) {
    if(cohset$CohortType[j] == 'Control for combination intervention'){
      cohset$CombC_label <- cohset[j,lab_n]
      cohset$CombC_n <- cohset[j,n_n]
      cohset$CombC_m <- cohset[j,m_n]
      cohset$CombC_v <- cohset[j,v_n]
    }
  }
  for (j in 1:chhrows) {
    if(cohset$CohortType[j] == 'Sham for TAAR1KO comb.'){
      cohset$ShamTAAR_label <- cohset[j,lab_n]
      cohset$ShamTAAR_n <- cohset[j,n_n]
      cohset$ShamTAAR_m <- cohset[j,m_n]
      cohset$ShamTAAR_v <- cohset[j,v_n]
    }
  }
  for (j in 1:chhrows) {
    if(cohset$CohortType[j] == 'Negative control for TAAR1KO comb.'){
      cohset$CTAAR_label <- cohset[j,lab_n]
      cohset$CTAAR_n <- cohset[j,n_n]
      cohset$CTAAR_m <- cohset[j,m_n]
      cohset$CTAAR_v <- cohset[j,v_n]
    }
  }
  
  data_org <- bind_rows(data_org, cohset)  
}

data_SI <- subset(data_org, data_org$CohortType== 'Simple intervention')
data_CI <- subset(data_org, data_org$CohortType== 'Combination intervention')
data_TI <- subset(data_org, data_org$CohortType== 'Intervention and TAAR1KO comb.')

### get standard placement - Single I ###
data_SI$F_C_L <- data_SI$Control_label$CohortLabel
data_SI$F_C_n <- data_SI$Control_n$NumberOfAnimals
data_SI$F_C_m <- data_SI$Control_m$OutcomeResult
data_SI$F_C_v <- data_SI$Control_v$OutcomeError

data_SI$F_T_L <- data_SI$CohortLabel
data_SI$F_T_n <- data_SI$NumberOfAnimals
data_SI$F_T_m <- data_SI$OutcomeResult
data_SI$F_T_v <- data_SI$OutcomeError

data_SI$F_S_L <- data_SI$Sham_label$CohortLabel
data_SI$F_S_n <- data_SI$Sham_n$NumberOfAnimals
data_SI$F_S_m <- data_SI$Sham_m$OutcomeResult
data_SI$F_S_v <- data_SI$Sham_v$OutcomeError

### get standard placement - Combined I ###
data_CI$F_C_L <- data_CI$CombC_label$CohortLabel
data_CI$F_C_n <- data_CI$CombC_n$NumberOfAnimals
data_CI$F_C_m <- data_CI$CombC_m$OutcomeResult
data_CI$F_C_v <- data_CI$CombC_v$OutcomeError

data_CI$F_T_L <- data_CI$CohortLabel
data_CI$F_T_n <- data_CI$NumberOfAnimals
data_CI$F_T_m <- data_CI$OutcomeResult
data_CI$F_T_v <- data_CI$OutcomeError

data_CI$F_S_L <- data_CI$Sham_label$CohortLabel
data_CI$F_S_n <- data_CI$Sham_n$NumberOfAnimals
data_CI$F_S_m <- data_CI$Sham_m$OutcomeResult
data_CI$F_S_v <- data_CI$Sham_v$OutcomeError

### get standard placement - TAAR I ###
data_TI$F_C_L <- data_TI$CTAAR_label$CohortLabel
data_TI$F_C_n <- data_TI$CTAAR_n$NumberOfAnimals
data_TI$F_C_m <- data_TI$CTAAR_m$OutcomeResult
data_TI$F_C_v <- data_TI$CTAAR_v$OutcomeError

data_TI$F_T_L <- data_TI$CohortLabel
data_TI$F_T_n <- data_TI$NumberOfAnimals
data_TI$F_T_m <- data_TI$OutcomeResult
data_TI$F_T_v <- data_TI$OutcomeError

data_TI$F_S_L <- data_TI$ShamTAAR_label$CohortLabel
data_TI$F_S_n <- data_TI$ShamTAAR_n$NumberOfAnimals
data_TI$F_S_m <- data_TI$ShamTAAR_m$OutcomeResult
data_TI$F_S_v <- data_TI$ShamTAAR_v$OutcomeError

start_flag <- match('Control_label', names(data_CI))
end_flag <- match('CTAAR_v', names(data_CI))

data_CI_F <- data_CI[,-c(start_flag:end_flag)]
data_SI_F <- data_SI[,-c(start_flag:end_flag)]
data_TI_F <- data_TI[,-c(start_flag:end_flag)]

data_all_F <- rbind(data_CI_F, data_SI_F, data_TI_F)

data_all_F <- data_all_F%>%
  mutate(drugname1 = case_when(
    grepl("RO5263397", data_all_F$`InterventionLabel[1]`) ~ "RO5263397",
    grepl("olanzepine", data_all_F$`InterventionLabel[1]`) ~ "olanzapine",
    grepl("SEP-363856", data_all_F$`InterventionLabel[1]`) ~ "SEP-363856",
    grepl("SEP-383856", data_all_F$`InterventionLabel[1]`) ~ "SEP-363856",
    grepl("SEP-856", data_all_F$`InterventionLabel[1]`) ~ "SEP-363856",
    grepl("SEP", data_all_F$`InterventionLabel[1]`) ~ "SEP",
    grepl("RO5203648", data_all_F$`InterventionLabel[1]`) ~ "RO5203648",
    grepl("LK000764", data_all_F$`InterventionLabel[1]`) ~ "LK000764",
    grepl("RO5256390", data_all_F$`InterventionLabel[1]`) ~ "RO5256390",
    grepl("SEP-856", data_all_F$`InterventionLabel[1]`) ~ "SEP-856",
    grepl("RO5073012", data_all_F$`InterventionLabel[1]`) ~ "RO5073012",
    grepl("Compound 50B", data_all_F$`InterventionLabel[1]`) ~ "Compound 50B",
    grepl("Compound 50A", data_all_F$`InterventionLabel[1]`) ~ "Compound 50A",
    grepl("RO5166017", data_all_F$`InterventionLabel[1]`) ~ "RO5166017",
    grepl("risperidone", data_all_F$`InterventionLabel[1]`) ~ "risperidone",
    grepl("Olanzapine", data_all_F$`InterventionLabel[1]`) ~ "olanzapine",
    grepl("OLZ", data_all_F$`InterventionLabel[1]`) ~ "olanzapine",
    grepl("TAAR1 KO", data_all_F$`InterventionLabel[1]`) ~ "TAAR1 KO",
    grepl("AP163", data_all_F$`InterventionLabel[1]`) ~ "AP163",
    TRUE ~ "Other"
  ))
data_all_F <- data_all_F%>%
  mutate(drugname2 = case_when(
    grepl("RO5263397", data_all_F$`InterventionLabel[2]`) ~ "RO5263397",
    grepl("olanzepine", data_all_F$`InterventionLabel[2]`) ~ "olanzapine",
    grepl("SEP-363856", data_all_F$`InterventionLabel[2]`) ~ "SEP-363856",
    grepl("SEP-383856", data_all_F$`InterventionLabel[1]`) ~ "SEP-363856",
    grepl("SEP-856", data_all_F$`InterventionLabel[2]`) ~ "SEP-363856",
    grepl("SEP", data_all_F$`InterventionLabel[2]`) ~ "SEP",
    grepl("RO5203648", data_all_F$`InterventionLabel[2]`) ~ "RO5203648",
    grepl("LK000764", data_all_F$`InterventionLabel[2]`) ~ "LK000764",
    grepl("RO5256390", data_all_F$`InterventionLabel[2]`) ~ "RO5256390",
    grepl("SEP-856", data_all_F$`InterventionLabel[2]`) ~ "SEP-856",
    grepl("RO5073012", data_all_F$`InterventionLabel[2]`) ~ "RO5073012",
    grepl("Compound 50B", data_all_F$`InterventionLabel[2]`) ~ "Compound 50B",
    grepl("Compound 50A", data_all_F$`InterventionLabel[2]`) ~ "Compound 50A",
    grepl("RO5166017", data_all_F$`InterventionLabel[2]`) ~ "RO5166017",
    grepl("risperidone", data_all_F$`InterventionLabel[2]`) ~ "risperidone",
    grepl("Olanzapine", data_all_F$`InterventionLabel[2]`) ~ "olanzapine",
    grepl("OLZ", data_all_F$`InterventionLabel[2]`) ~ "olanzapine",
    grepl("TAAR1 KO", data_all_F$`InterventionLabel[2]`) ~ "TAAR1 KO",
    grepl("AP163", data_all_F$`InterventionLabel[2]`) ~ "AP163",
    TRUE ~ "Other"
  ))

savename_SI <- paste0(LSR,'dataSI_',Sys.Date(),'.csv')
savename_CI <- paste0(LSR,'dataCI_',Sys.Date(),'.csv')
savename_TI <- paste0(LSR,'dataTI_',Sys.Date(),'.csv')
savename_all <- paste0(LSR,'data_all_',Sys.Date(),'.csv')

write_csv(data_CI_F, savename_CI) #move to LSR3data in RMA folder 
write_csv(data_SI_F, savename_SI)
write_csv(data_TI_F, savename_TI)
write_csv(data_all_F, savename_all)


# Data are ready to calculate effect sizes with numerical data for one comparison (test/control) per line

# Note: current code doesn't have cohort level questions split into treatment and control annotations per line
# Not an issue currently since control and treatment cohorts have had the same characteristics (sex, strain, etc) so far, but probably the subject of future edits (e.g. using pivot_wider function)

#####################################################################################################################################################
# Calculate effect size for simple interventions

## Read in data and edit outcome label 
#pass though better
# dataall <- read_csv("dataall.csv")
dataall <- data_all_F


dataall <- dataall %>% 
  mutate(OutcomeType = case_when(str_detect(OutcomeLabel, regex("lma|locomotor|oft|horizontal|distance|vertical|climbing", ignore_case = TRUE)) ~ "Locomotor activity",
                                 str_detect(OutcomeLabel, "PPI") ~ "Prepulse inhibition", 
                                 TRUE ~ `Type of outcome?`)) %>% 
  relocate(OutcomeType, .after = OutcomeLabel)

outcomeFrequencies <- dataall %>% group_by(OutcomeType, OutcomeLabel, `Type of outcome?`) %>% count()

savename_of <- paste0(LSR,'_',Sys.Date(),'.csv')

write.csv(outcomeFrequencies, savename_of)
## 1. Calculate SD for all comparisons 
#F = final
#C/T/S = control/treatment/sham
#L/n/m/v = cohort label/n in cohort/mean/variance

dataall <- dataall %>% 
  mutate(F_C_v.SD = case_when(ErrorType == "IQR" ~ (F_C_v/1.35), 
                              ErrorType == "SD" ~ F_C_v, 
                              ErrorType == "SEM" ~ sqrt(F_C_n)*F_C_v)) %>%  
  relocate(F_C_v.SD, .after = F_C_v) %>% 
  mutate(F_T_v.SD = case_when(ErrorType == "IQR" ~ (F_T_v/1.35), 
                              ErrorType == "SD" ~ F_T_v, 
                              ErrorType == "SEM" ~ sqrt(F_T_n)*F_T_v)) %>%  
  relocate(F_T_v.SD, .after = F_T_v) %>% 
  mutate(F_S_v.SD = case_when(ErrorType == "IQR" ~ (F_S_v/1.35), 
                              ErrorType == "SD" ~ F_S_v, 
                              ErrorType == "SEM" ~ sqrt(F_S_n)*F_S_v)) %>%  
  relocate(F_S_v.SD, .after = F_S_v)

# Check number of comparisons where NMD can be calculated

dataall <- dataall %>%
  rowwise() %>%
  mutate(`NMD_possible` = all(!is.na(F_S_L)))

## 2. Calculate effect sizes (SMD for all, NMD where possible)

### Calculate true n for control groups (n'c)

# Step 1: Number of groups served by control group
F_C_L_frequencies <- dataall %>%
  group_by(StudyId, OutcomeId, F_C_L) %>%
  summarise(Frequency_FCL = n()) %>% 
  ungroup()

# Step 2: Join the frequencies back to the original dataframe
dataall <- dataall %>%
  left_join(F_C_L_frequencies)

# Step 3: Divide F_C_n by the frequency count
dataall <- dataall %>%
  mutate(F_C_n_true = F_C_n / Frequency_FCL)

### SMD

#Hedges g to account for small sample sizes (default for SMD when using the escalc() function) - Hedge’s g (statistically corrects for variance that may be introduced when sample sizes are small (Larry V. Hedges 1981))
# m1i = control, m2i = rx 

# Hedge's g effect size
SMD_data_all.nottrue <- escalc(
  measure = "SMD", 
  m1i = dataall$F_C_m, 
  m2i = dataall$F_T_m, 
  sd1i = dataall$F_C_v.SD, 
  sd2i = dataall$F_T_v.SD, 
  n1i = dataall$F_C_n, 
  n2i = dataall$F_T_n, 
  data = dataall) %>% 
  select(yi, vi)

dataall$SMD <- SMD_data_all.nottrue$yi
dataall$SMDv <- SMD_data_all.nottrue$vi

#escalc (m1 - m2) = (control - treatment)
# so if greater is better, then *-1

#dataall$SMD_true <- SMD_data_all.true$yi #Row 112, SMD calculated, but SMD_true not calculated
#dataall$SMDv_true <- SMD_data_all.true$vi

### NMD

# Assume that treatments are closer to shams than controls are 
# So C-S > T-S 


dataall <- dataall %>% 
  mutate(`NMD` = 100*(((F_C_m - F_S_m) - (F_T_m - F_S_m))/(F_C_m - F_S_m))) %>% 
  mutate(`NMD_SDc*` = 100*((F_C_v.SD/(F_C_m - F_S_m)))) %>% 
  mutate(`NMD_SDrx*` = 100*((F_T_v.SD/(F_C_m - F_S_m)))) %>% 
  mutate(NMDv = sqrt(((`NMD_SDc*`)^2/F_C_n) + ((`NMD_SDrx*`)^2/F_T_n))) 

dataall.direction <- dataall %>% 
  mutate(SMD = if_else((GreaterIsWorse == "FALSE"), -1*SMD, SMD)) %>% 
  rename(CategoryDiseaseInduction = `Category of disease model induction method:`)

diagnostic <- dataall.direction %>% 
  select(F_C_m, F_T_m, F_S_m, SMD, NMD, GreaterIsWorse) %>% 
  mutate(CbiggerT = if_else(F_C_m > F_T_m, "Yes", "No"))

# nicely named columns for subgroup analysis

df <- dataall.direction

# Correct "SEP" : All SEP VALUES FOR DRUG NAMES ARE SEP-363856. Checked each paper on 14.12.23. StudyID's: eba6e60f, c064173a, 84b834
df <- df %>% 
  mutate(drugname1 = str_replace_all(drugname1, "SEP(?!-363856)", "SEP-363856"))

### Duration of treatment (categorical)

# Create variable standardised to Weeks
df <- df %>% 
  mutate(`Duration of treatment[1]` = as.numeric(`Duration of treatment[1]`)) %>% 
  mutate(DurationOfTreatmentWeeks = case_when(`Unit of measurement for treatment duration[1]` == "Days" ~ `Duration of treatment[1]`/7, 
                                              `Unit of measurement for treatment duration[1]` == "Months" ~ `Duration of treatment[1]`/4.345, 
                                              `Unit of measurement for treatment duration[1]` == "Weeks" ~ `Duration of treatment[1]`, 
                                              `Unit of measurement for treatment duration[1]` == "Single dose" ~ 1/7))

# Create categorical variable for Duration of treatment. Grouped into up to a week, between a week and 4 weeks, more than 4 weeks
df <- df %>% 
  mutate(TreatmentDurationCategory = case_when(DurationOfTreatmentWeeks <= 1 ~ "Less than 1 week", 
                                               DurationOfTreatmentWeeks > 1 & DurationOfTreatmentWeeks < 4 ~ "Between 1-4 weeks", 
                                               DurationOfTreatmentWeeks >= 4 ~ "More than 4 weeks"))


### Prophylactic or therapeutic

df <- df %>% 
  mutate(ProphylacticOrTherapeutic = case_when(`Timing of treatment administration:[1]` == "After disease model induction" ~ "Therapeutic", 
                                               TRUE ~ "Prophylactic"))

## Give original variables better names for analysis in new columns - this is all for simple 

df <- df %>% 
  mutate(Species = `Species of animals?`, 
         Strain = `Animal strain?`,
         Sex = `Sex of animals?`, 
         DrugName = drugname1,  #change for combination
         InterventionAdministrationRoute = `Treatment administration route:[1]`, #change for combination
         DoseOfIntervention_mgkg = `Dose of treatment used:[1]`)  #FOR LSR3 THESE ARE ALL IN UNIT mg/kg SO NO CONVERSION NEEDED

## Categorise by outcome type - requires checking with each iteration
df <-  df %>%
  mutate(outcome_type = case_when(
    OutcomeType == "Locomotor activity" ~ "Locomotor activity",
    OutcomeType == "Prepulse inhibition" ~ "Prepulse inhibition",
    OutcomeType == "Neurobiological outcome (mechanistic) e.g. dopaminergic/seritonergic/glutamatergic signalling" ~ "Neurobiological outcome",
    str_detect(OutcomeLabel, regex("MWM|reward|soon|NOR|attention|soon", ignore_case = TRUE)) ~ "Cognition",
    str_detect(OutcomeLabel, regex("social", ignore_case = TRUE)) ~ "Social interaction",
    str_detect(OutcomeLabel, regex("stereotyp", ignore_case = TRUE)) ~ "Stereotypy",
    str_detect(OutcomeLabel, regex("movements", ignore_case = TRUE)) ~ "Locomotor activity",
    TRUE ~ "Other"
  )) 

## Add drug characteristics - Info taken from Taar1_drugs_spiros_18.12.23.xlsx
df <- df %>% 
  mutate(pE50 = case_when(DrugName == "LK000764" ~ 8.40, 
                          DrugName == "RO5073012" ~ 7.64, 
                          DrugName == "RO5166017" ~ 7.23,
                          DrugName == "RO5203648" ~ 7.52,
                          DrugName == "RO5256390" ~ 7.74,
                          DrugName == "RO5263397" ~ 7.77,
                          DrugName == "Ulotaront" | DrugName == "SEP-363856" ~ 6.85,
                          DrugName == "Ralmitaront" | DrugName == "RO6889450" ~ 7.23,
                          DrugName == "AP163" ~ 6.95,
                          DrugName == "Compound 50A" ~ 5.20,
                          DrugName == "Compound 50B" ~ 6.39,
                          TRUE ~ NA_real_)) %>% 
  mutate(Efficacy = case_when(DrugName == "LK000764" ~ "TAAR1 full agonist", 
                              DrugName == "RO5073012" ~ "TAAR1 partial agonist", 
                              DrugName == "RO5166017" ~ "TAAR1 full agonist",
                              DrugName == "RO5203648" ~ "TAAR1 partial agonist",
                              DrugName == "RO5256390" ~ "TAAR1 full agonist",
                              DrugName == "RO5263397" ~ "TAAR1 partial agonist",
                              DrugName == "Ulotaront" | DrugName == "SEP-363856" ~ "TAAR1 full agonist",
                              DrugName == "Ralmitaront" | DrugName == "RO6889450"~ "TAAR1 partial agonist",
                              DrugName == "AP163" ~ "TAAR1 full agonist",
                              DrugName == "Compound 50A" ~ "TAAR1 partial agonist",
                              DrugName == "Compound 50B" ~ "TAAR1 partial agonist", 
                              TRUE ~ "NA")) %>% 
  mutate(Selectivity = case_when(DrugName == "LK000764" ~ "Unclear", 
                                 DrugName == "RO5073012" ~ "High", 
                                 DrugName == "RO5166017" ~ "High",
                                 DrugName == "RO5203648" ~ "High",
                                 DrugName == "RO5256390" ~ "High",
                                 DrugName == "RO5263397" ~ "High",
                                 DrugName == "Ulotaront" | DrugName == "SEP-363856" ~ "Low (5HT1A partial agonism)",
                                 DrugName == "Ralmitaront" | DrugName == "RO6889450" ~ "Unclear",
                                 DrugName == "AP163" ~ "Unclear",
                                 DrugName == "Compound 50A" ~ "Unclear",
                                 DrugName == "Compound 50B" ~ "Low (5HT1A partial agonism)", 
                                 TRUE ~ "NA")) %>% 
  mutate(MolarMass = case_when(DrugName == "LK000764" ~ 299.1058, 
                               DrugName == "RO5073012" ~ 249.742, 
                               DrugName == "RO5166017" ~ 219.288,
                               DrugName == "RO5203648" ~ 231.08,
                               DrugName == "RO5256390" ~ 218.29,
                               DrugName == "RO5263397" ~ 194.21,
                               DrugName == "Ulotaront" | DrugName == "SEP-363856" ~ 183.27,
                               DrugName == "Ralmitaront" | DrugName == "RO6889450" ~ 314.38,
                               DrugName == "AP163" ~ 282.1368,
                               DrugName == "Compound 50A" ~ 170.3,
                               DrugName == "Compound 50B" ~ 170.3, 
                               TRUE ~ NA_real_)) %>% 
  mutate(EC50mM = case_when(DrugName == "LK000764" ~ 0.004, 
                            DrugName == "RO5073012" ~ 0.023, 
                            DrugName == "RO5166017" ~ 0.059,
                            DrugName == "RO5203648" ~ 0.03,
                            DrugName == "RO5256390" ~ 0.018,
                            DrugName == "RO5263397" ~ 0.017,
                            DrugName == "Ulotaront" | DrugName == "SEP-363856" ~ 0.14,
                            DrugName == "Ralmitaront" | DrugName == "RO6889450" ~ 0.059,
                            DrugName == "AP163" ~ 0.112,
                            DrugName == "Compound 50A" ~ 6.25,
                            DrugName == "Compound 50B" ~ 0.405, 
                            TRUE ~ NA_real_)) 


# Calculate standardised dose
df <- df %>% 
  mutate(DoseOfIntervention_mgkg = as.numeric(DoseOfIntervention_mgkg)) %>% 
  mutate(StandardisedDose = log(DoseOfIntervention_mgkg)/(MolarMass*EC50mM*0.001)) 

# SAVE FILE
savefile_output <- paste0(LSR,'_','clean_data_',Sys.Date(),'.csv')
write.csv(df, savefile_output, row.names = FALSE)


