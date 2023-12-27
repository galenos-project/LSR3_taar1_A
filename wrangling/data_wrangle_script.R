library(dplyr)
library(stringr)
library(tibble)
library(tidyverse)
library(tidyr)
library(purrr)
library(metafor)

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
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "TAAR1KO" & (Treatment2Type == "Negative control"| is.na(Treatment2Type))) |
       (Treatment2Type == "TAAR1KO" & (Treatment1Type == "Negative control"| is.na(Treatment1Type)))) ~ "Negative control for TAAR1KO",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "Negative control"| is.na(Treatment1Type)) & (Treatment2Type == "Negative control"| is.na(Treatment2Type))) ~ "Negative control",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "Positive control" & (Treatment2Type == "Negative control"| is.na(Treatment2Type))) |
         (Treatment2Type == "Positive control" & (Treatment1Type == "Negative control"| is.na(Treatment1Type)))) ~ "Positive control",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "Intervention" & (Treatment2Type == "Negative control"| is.na(Treatment2Type))) |
         (Treatment2Type == "Intervention" & (Treatment1Type == "Negative control"| is.na(Treatment1Type)))) ~ "Simple intervention",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "TAAR1KO" & Treatment2Type == "Positive control") |
         (Treatment2Type == "TAAR1KO" & Treatment1Type == "Positive control")) ~ "Positive control TAAR1KO",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "TAAR1KO" & Treatment2Type == "Intervention") |
         (Treatment2Type == "TAAR1KO" & Treatment1Type == "Intervention")) ~ "Intervention for TAAR1KO",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "TAAR1KO" & Treatment2Type == "Positive control") |
         (Treatment2Type == "TAAR1KO" & Treatment1Type == "Positive control")) ~ "Positive control for TAAR1KO",
    
    IsDiseaseModelControl == FALSE & 
      ((Treatment1Type == "Intervention" & Treatment2Type == "Positive control") | 
         (Treatment2Type == "Intervention" & Treatment1Type == "Positive control")) ~ "Combination intervention",

    IsDiseaseModelControl == TRUE & 
      ((Treatment1Type == "TAAR1KO" & (Treatment2Type == "Negative control"| is.na(Treatment2Type))) |
         (Treatment2Type == "TAAR1KO" & (Treatment1Type == "Negative control"| is.na(Treatment1Type)))) ~ "Sham TAAR1KO",
    
    IsDiseaseModelControl == TRUE & 
      ((Treatment1Type == "Negative control"| is.na(Treatment1Type)) & (Treatment2Type == "Negative control"| is.na(Treatment2Type))) ~ "Sham",
        
    IsDiseaseModelControl == TRUE & 
          ((Treatment1Type == "Positive control" & (Treatment2Type == "Negative control"| is.na(Treatment2Type))) |
             (Treatment2Type == "Positive control" & (Treatment1Type == "Negative control"| is.na(Treatment1Type)))) ~ "Sham Positive control",
        
    IsDiseaseModelControl == TRUE & 
          ((Treatment1Type == "Intervention" & (Treatment2Type == "Negative control"| is.na(Treatment2Type))) |
             (Treatment2Type == "Intervention" & (Treatment1Type == "Negative control"| is.na(Treatment1Type)))) ~ "Sham intervention",
        
    IsDiseaseModelControl == TRUE & 
          ((Treatment1Type == "TAAR1KO" & Treatment2Type == "Positive control") |
             (Treatment2Type == "TAAR1KO" & Treatment1Type == "Positive control")) ~ "Sham TAAR1KO Positive control",
        
    IsDiseaseModelControl == TRUE & 
          ((Treatment1Type == "TAAR1KO" & Treatment2Type == "Intervention") |
             (Treatment2Type == "TAAR1KO" & Treatment1Type == "Intervention")) ~ "Sham TAAR1KO intervention",
        
    IsDiseaseModelControl == TRUE & 
          ((Treatment1Type == "TAAR1KO" & Treatment2Type == "Positive control") |
             (Treatment2Type == "TAAR1KO" & Treatment1Type == "Positive control")) ~ "Sham TAAR1KO Positive control intervention for TAAR1KO comb",
        
    IsDiseaseModelControl == TRUE & 
          ((Treatment1Type == "Intervention" & Treatment2Type == "Positive control") | 
             (Treatment2Type == "Intervention" & Treatment1Type == "Positive control")) ~ "Sham Positive control intervention"
      )) %>% 
  relocate(CohortType, .after = `InterventionID[2]`) %>% 
  relocate(c(Treatment1Type, Treatment2Type), .after = `InterventionLabel[2]`) %>% 
  relocate(IsDiseaseModelControl, .after = ModelID)

reconciled_cohort_type <- reconciled_cohort_type %>%
  mutate(GroupID = interaction(StudyId, ExperimentID, OutcomeId, TimeInMinute))


## Label how positive controls are being used in positive control v TAAR1 Ag experiments
# If for an observation with the value "Positive control" in the CohortType column, there is an observation 
# with the CohortType "Combination intervention" with the same combination of values in the StudyId, ExperimentID, 
# OutcomeId, TimeInMinute variables, then the value in the CohortType column needs changed to "Control for 
# combination intervention". Similarly, if for an observation with the value "Positive control treated sham" 
# in the CohortType column, there is an observation with the CohortType "Combination intervention" with 
# the same combination of values in the StudyId, ExperimentID, OutcomeId, TimeInMinute variables, then the 
# value in the CohortType column needs changed to "Sham for combination intervention"

## Identify experiments that are just 'Simple' comparison, those which are 'Combination' comparisons and those which have 'TAAR1KO' involved
reconciled_comparison_type <- reconciled_cohort_type %>%
  mutate(ExperimentType = case_when(
    str_detect(CohortType, "combination") | str_detect(CohortType, "Combination") ~ "Antipsychotic combination",
    str_detect(CohortType, "TAAR1")  ~ "TAAR1KO background",
    str_detect(CohortType, "simple") | str_detect(CohortType, "Simple") ~ "Simple intervention",
    str_detect(CohortType, "positive") ~ "Antipsychotic drug"
  ))

reconciled_cohort_role <- reconciled_comparison_type %>% 
  mutate(RoleOfCohort = case_when(
    str_detect(CohortType, "Sham") | str_detect(CohortType, "sham") ~ "S", 
    str_detect(CohortType, "Control") | str_detect(CohortType, "control") ~ "C",
    str_detect(CohortType, "Intervention") | str_detect(CohortType, "intervention") ~ "I"
  ))


## Wrangle wide so each observations is a single comparison
data <- reconciled_cohort_role

####following section can be removed when SyRF records have been corrected
dfdelta <- data
timecol <- match('TimeInMinute', names(dfdelta))

groupid_col <- match('GroupId', names(dfdelta))
dfdelta[325,138] <- 1
dfdelta[280,138] <- 1
dfdelta[277,138] <- 1
data <- dfdelta
data <- data %>%
  mutate(GroupID = interaction(StudyId, ExperimentID, OutcomeId, TimeInMinute))

#### SyRF correction section ends


# remove DA knockouts
data_rem <- subset(data, data$`DiseaseModelLabel(s)` == 'DAT +/-')
data <- anti_join(data, data_rem)
#remove data reported from subgroups
data_rem <- data %>% filter(str_detect(`DiseaseModelLabel(s)`, "SUBGROUP") | str_detect(CohortLabel, 'SUBGROUP'))
data <- anti_join(data, data_rem)
data <- data %>%
  mutate(`DiseaseModelLabel(s)` = ifelse(`DiseaseModelLabel(s)` == "DAT -/-", "DAT KO", `DiseaseModelLabel(s)`))

lab_n <- match('CohortLabel', names(data))
n_n <- match('NumberOfAnimals', names(data))
m_n <- match('OutcomeResult', names(data))
v_n <- match('OutcomeError', names(data))

###for each cohort, label sham and control groups, along with TAAR1KO -ve control Contol for combination interventions (== positive control)
groups <- as.data.frame(unique(data$GroupID))
groups$`unique(data$GroupID)` <- as.character(groups$`unique(data$GroupID)`)

for (i in 1:nrow(groups)) {
  group <- groups[i, 1]
  cohset <- subset(data, GroupID == group)
  
  cohset <- process_cohort(cohset, 'Negative control', 'Control')
  cohset <- process_cohort(cohset, 'Sham', 'Sham')
  cohset <- process_cohort(cohset, 'Sham intervention', 'ShInt')
  cohset <- process_cohort(cohset, 'Sham TAAR1KO', 'ShamTAAR')
  cohset <- process_cohort(cohset, 'Negative control for TAAR1KO', 'CTAAR')
  
  # For the first group, update 'data_org' directly
  if (i == 1) {
    data_org <- cohset
  } else {
    data_org <- bind_rows(data_org, cohset)
  }
}

### identify studies which have positive control v TAAR1 Ag with new label comp_exp
data$comp_exp <- "FALSE"

for (i in 1:nrow(groups)) {
  group <- groups[i, 1]
  cohset <- subset(data, GroupID == group)
  if (any(cohset$CohortType == "Positive control")) {
    data[data$GroupID == group, "comp_exp"] <- "TRUE"
  }
}

### select Groups which include a positive control
datapc <- subset(data, data$comp_exp == "TRUE")

groupspc <- as.data.frame(unique(datapc$GroupID))
groupspc$`unique(datapc$GroupID)` <- as.character(groupspc$`unique(datapc$GroupID)`)

for (i in 1:nrow(groups)) {
  grouppc <- groupspc[i, 1]
  cohset <- subset(datapc, GroupID == grouppc)
  
  cohset <- process_cohort(cohset, "Positive control", "Control")
  cohset <- process_cohort(cohset, "Sham", "ShamC")
  
  if (i == 1) {
    data_org_p <- cohset
  } else {
    data_org_p <- bind_rows(data_org_p, cohset)
  }
}

data_org_p1 <- subset(data_org_p, data_org_p$CohortType == "Simple intervention")
data_org_p1$ExperimentType <- "Head to head"

data_org <- bind_rows(data_org, data_org_p1)

data_SI <- subset(data_org, (data_org$CohortType== 'Simple intervention' & !data_org$ExperimentType == "Head to head"))
data_CI <- subset(data_org, data_org$CohortType== 'Combination intervention')
data_TI <- subset(data_org, data_org$CohortType== 'Intervention for TAAR1KO')
data_PC <- subset(data_org, data_org$ExperimentType == "Head to head")

### get standard placement - Single I ###
data_SI$F_C_L <- data_SI$Control_label
data_SI$F_C_n <- data_SI$Control_n
data_SI$F_C_m <- data_SI$Control_m
data_SI$F_C_v <- data_SI$Control_v

data_SI$F_T_L <- data_SI$CohortLabel
data_SI$F_T_n <- data_SI$NumberOfAnimals
data_SI$F_T_m <- data_SI$OutcomeResult
data_SI$F_T_v <- data_SI$OutcomeError

data_SI$F_S_L <- data_SI$Sham_label
data_SI$F_S_n <- data_SI$Sham_n
data_SI$F_S_m <- data_SI$Sham_m
data_SI$F_S_v <- data_SI$Sham_v

### get standard placement - Combined I ###
data_CI$F_C_L <- data_CI$Control_label
data_CI$F_C_n <- data_CI$Control_n
data_CI$F_C_m <- data_CI$Control_m
data_CI$F_C_v <- data_CI$Control_v

data_CI$F_T_L <- data_CI$CohortLabel
data_CI$F_T_n <- data_CI$NumberOfAnimals
data_CI$F_T_m <- data_CI$OutcomeResult
data_CI$F_T_v <- data_CI$OutcomeError

data_CI$F_S_L <- data_CI$Sham_label
data_CI$F_S_n <- data_CI$Sham_n
data_CI$F_S_m <- data_CI$Sham_m
data_CI$F_S_v <- data_CI$Sham_v

### get standard placement - TAAR I ###
data_TI$F_C_L <- data_TI$CTAAR_label
data_TI$F_C_n <- data_TI$CTAAR_n
data_TI$F_C_m <- data_TI$CTAAR_m
data_TI$F_C_v <- data_TI$CTAAR_v

data_TI$F_T_L <- data_TI$CohortLabel
data_TI$F_T_n <- data_TI$NumberOfAnimals
data_TI$F_T_m <- data_TI$OutcomeResult
data_TI$F_T_v <- data_TI$OutcomeError

data_TI$F_S_L <- data_TI$ShamTAAR_label
data_TI$F_S_n <- data_TI$ShamTAAR_n
data_TI$F_S_m <- data_TI$ShamTAAR_m
data_TI$F_S_v <- data_TI$ShamTAAR_v

### get standard placement - PC ###
data_PC$F_C_L <- data_PC$Control_label
data_PC$F_C_n <- data_PC$Control_n
data_PC$F_C_m <- data_PC$Control_m
data_PC$F_C_v <- data_PC$Control_v

data_PC$F_T_L <- data_PC$CohortLabel
data_PC$F_T_n <- data_PC$NumberOfAnimals
data_PC$F_T_m <- data_PC$OutcomeResult
data_PC$F_T_v <- data_PC$OutcomeError

data_PC$F_S_L <- data_PC$ShamC_label
data_PC$F_S_n <- data_PC$ShamC_n
data_PC$F_S_m <- data_PC$ShamC_m
data_PC$F_S_v <- data_PC$ShamC_v

start_flag <- match('Control_label', names(data_CI))
end_flag <- match('CTAAR_v', names(data_CI))

data_CI_F <- data_CI[,-c(start_flag:end_flag)]
data_SI_F <- data_SI[,-c(start_flag:end_flag)]
data_TI_F <- data_TI[,-c(start_flag:end_flag)]

data_all_F <- bind_rows(data_CI, data_SI, data_TI, data_PC)

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
savename_PC <- paste0(LSR,'dataPC_',Sys.Date(),'.csv')
savename_all <- paste0(LSR,'data_all_',Sys.Date(),'.csv')

#write_csv(data_CI, savename_CI) #move to LSR3data in RMA folder 
#write_csv(data_SI, savename_SI)
#write_csv(data_TI, savename_TI)
#write_csv(data_PC, savename_PC)
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
                              DrugName == "RO5263397" ~ "TAAR1 full agonist",
                              DrugName == "Ulotaront" | DrugName == "SEP-363856" ~ "TAAR1 partial agonist",
                              DrugName == "Ralmitaront" | DrugName == "RO6889450"~ "TAAR1 partial agonist",
                              DrugName == "AP163" ~ "TAAR1 full agonist",
                              DrugName == "Compound 50A" ~ "TAAR1 partial agonist",
                              DrugName == "Compound 50B" ~ "TAAR1 partial agonist", 
                              TRUE ~ "NA")) %>% 
  mutate(Selectivity = case_when(DrugName == "LK000764" ~ "Unclear", 
                                 DrugName == "RO5073012" ~ "High", 
                                 DrugName == "RO5166017" ~ "High",
                                 DrugName == "RO5203648" ~ "High",
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

# Calculate standardised dose for overall dose-response meta-regression
df <- df %>% 
  mutate(DoseOfIntervention_mgkg = as.numeric(DoseOfIntervention_mgkg)) %>% 
  mutate(StandardisedDose = (log(DoseOfIntervention_mgkg))/((MolarMass*1000)*(EC50mM/1000000))) 


###### For RoB subgroup analysis ######
df <- df %>%
  rename(`(RoB) Were caregivers/investigator blinded to which intervention each animal received?` = `Were caregivers/investigator blinded to which intervention each animal received?`)

# Calculate overall RoB score

df <- df %>%
  rowwise() %>%
  mutate(RoBScore = sum(c_across(contains("RoB")) == "Yes", na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(RoBScore = paste0(RoBScore, " criteria met"))

##### For reporting quality subgroup analysis #####
df <- df %>% 
  rename(`(ARRIVE) Is any role of the funder in the design/analysis/reporting of the study described?` = `Is any role of the funder in the design/analysis/reporting of study described?`)

# Calculate overall ARRIVE score
df <- df %>%
  mutate(ARRIVEScore = rowSums(across(contains("ARRIVE"), ~ (.x == "Yes") | (.x == "NA (ethical approval declared)")), na.rm = TRUE)) %>% 
  mutate(ARRIVEScoreCat = case_when(ARRIVEScore <= 3 ~ "A: < 3 criteria met",
                                    ARRIVEScore > 3 & ARRIVEScore <= 7 ~ "B: 4-7 criteria met",
                                    ARRIVEScore > 7 & ARRIVEScore <= 11 ~ "C: 8-11 criteria met",
                                    ARRIVEScore > 11 & ARRIVEScore <= 15 ~ "D: 12-15 criteria met",
                                    ARRIVEScore > 15 & ARRIVEScore <= 19 ~ "E: 16-19 criteria met",
                                    ARRIVEScore > 19 ~ "F: > 20 criteria met")) %>% 
  mutate(ARRIVEScoreCat = as.factor(ARRIVEScoreCat))
                                  

# SAVE FILE
savefile_output <- paste0(LSR,'_','clean_data_',Sys.Date(),'.csv')
write.csv(df, savefile_output, row.names = FALSE)

