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

###### Tidying and cleaning the data ######
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

    IsDiseaseModelControl == TRUE ~ "Sham"
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
    str_detect(CohortType, "combination") | str_detect(CohortType, "Combination") ~ "Antipsychotic combination (TA v A)",
    str_detect(CohortType, "TAAR1")  ~ "TAAR1KO",
    str_detect(CohortType, "simple") | str_detect(CohortType, "Simple") ~ "Simple intervention (T v Cont)",
    str_detect(CohortType, "Positive") | str_detect(CohortType, "positive") ~ "Antipsychotic drug (A v Cont)"
  ))

reconciled_cohort_role <- reconciled_comparison_type %>% 
  mutate(RoleOfCohort = case_when(
    str_detect(CohortType, "Sham") | str_detect(CohortType, "sham") ~ "S", 
    str_detect(CohortType, "Control") | str_detect(CohortType, "control") ~ "C",
    str_detect(CohortType, "Intervention") | str_detect(CohortType, "intervention") ~ "I",
    str_detect(CohortType, "Positive") | str_detect(CohortType, "positive") ~ "P"
  ))


## Wrangle wide so each observations is a single comparison
data <- reconciled_cohort_role

####following section can be removed when SyRF records have been corrected
dfdelta <- data
dfdelta[339,138] <- 1
dfdelta[291,138] <- 1
dfdelta[294,138] <- 1
data <- dfdelta
data <- data %>%
  mutate(GroupID = interaction(StudyId, ExperimentID, OutcomeId, TimeInMinute))

#### SyRF correction section ends


##### Remove DA knockouts and data reported from subgroups #####
data_rem <- subset(data, data$`DiseaseModelLabel(s)` == 'DAT +/-')
data <- anti_join(data, data_rem)
#remove data reported from subgroups
data_rem <- data %>% filter(str_detect(`DiseaseModelLabel(s)`, "SUBGROUP") | str_detect(CohortLabel, 'SUBGROUP'))
data <- anti_join(data, data_rem)
data <- data %>%
  mutate(`DiseaseModelLabel(s)` = ifelse(`DiseaseModelLabel(s)` == "DAT -/-", "DAT KO", `DiseaseModelLabel(s)`))
data <- data %>% mutate_all(trimws)

##### Extract treatment names and doses ####
data <- data%>%
  mutate(drugname1 = case_when(
    grepl("RO5263397", data$`InterventionLabel[1]`) ~ "RO5263397",
    grepl("olanzepine", data$`InterventionLabel[1]`) ~ "olanzapine",
    grepl("SEP-363856", data$`InterventionLabel[1]`) ~ "SEP-363856",
    grepl("SEP-383856", data$`InterventionLabel[1]`) ~ "SEP-363856",
    grepl("SEP-856", data$`InterventionLabel[1]`) ~ "SEP-363856",
    grepl("SEP", data$`InterventionLabel[1]`) ~ "SEP",
    grepl("RO5203648", data$`InterventionLabel[1]`) ~ "RO5203648",
    grepl("LK000764", data$`InterventionLabel[1]`) ~ "LK000764",
    grepl("RO5256390", data$`InterventionLabel[1]`) ~ "RO5256390",
    grepl("SEP-856", data$`InterventionLabel[1]`) ~ "SEP-856",
    grepl("RO5073012", data$`InterventionLabel[1]`) ~ "RO5073012",
    grepl("Compound 50B", data$`InterventionLabel[1]`) ~ "Compound 50B",
    grepl("Compound 50A", data$`InterventionLabel[1]`) ~ "Compound 50A",
    grepl("RO5166017", data$`InterventionLabel[1]`) ~ "RO5166017",
    grepl("risperidone", data$`InterventionLabel[1]`) ~ "risperidone",
    grepl("Risperidone", data$`InterventionLabel[1]`) ~ "risperidone",
    grepl("clozapine", data$`InterventionLabel[1]`) ~ "clozapine",
    grepl("Clozapine", data$`InterventionLabel[1]`) ~ "clozapine",
    grepl("Aripiprazole", data$`InterventionLabel[1]`) ~ "aripiprazole",
    grepl("aripiprazole", data$`InterventionLabel[1]`) ~ "aripiprazole",
    grepl("Olanzapine", data$`InterventionLabel[1]`) ~ "olanzapine",
    grepl("OLZ", data$`InterventionLabel[1]`) ~ "olanzapine",
    grepl("TAAR1 KO", data$`InterventionLabel[1]`) ~ "TAAR1 KO",
    grepl("AP163", data$`InterventionLabel[1]`) ~ "AP163",
    TRUE ~ "Other"
  ))
data <- data%>%
  mutate(drugname2 = case_when(
    grepl("RO5263397", data$`InterventionLabel[2]`) ~ "RO5263397",
    grepl("olanzepine", data$`InterventionLabel[2]`) ~ "olanzapine",
    grepl("SEP-363856", data$`InterventionLabel[2]`) ~ "SEP-363856",
    grepl("SEP-383856", data$`InterventionLabel[1]`) ~ "SEP-363856",
    grepl("SEP-856", data$`InterventionLabel[2]`) ~ "SEP-363856",
    grepl("SEP", data$`InterventionLabel[2]`) ~ "SEP",
    grepl("RO5203648", data$`InterventionLabel[2]`) ~ "RO5203648",
    grepl("LK000764", data$`InterventionLabel[2]`) ~ "LK000764",
    grepl("RO5256390", data$`InterventionLabel[2]`) ~ "RO5256390",
    grepl("SEP-856", data$`InterventionLabel[2]`) ~ "SEP-856",
    grepl("RO5073012", data$`InterventionLabel[2]`) ~ "RO5073012",
    grepl("Compound 50B", data$`InterventionLabel[2]`) ~ "Compound 50B",
    grepl("Compound 50A", data$`InterventionLabel[2]`) ~ "Compound 50A",
    grepl("RO5166017", data$`InterventionLabel[2]`) ~ "RO5166017",
    grepl("risperidone", data$`InterventionLabel[2]`) ~ "risperidone",
    grepl("Risperidone", data$`InterventionLabel[2]`) ~ "risperidone",
    grepl("clozapine", data$`InterventionLabel[2]`) ~ "clozapine",
    grepl("Clozapine", data$`InterventionLabel[2]`) ~ "clozapine",
    grepl("Aripiprazole", data$`InterventionLabel[2]`) ~ "aripiprazole",
    grepl("aripiprazole", data$`InterventionLabel[2]`) ~ "aripiprazole",
    grepl("Olanzapine", data$`InterventionLabel[2]`) ~ "olanzapine",
    grepl("OLZ", data$`InterventionLabel[2]`) ~ "olanzapine",
    grepl("TAAR1 KO", data$`InterventionLabel[2]`) ~ "TAAR1 KO",
    grepl("AP163", data$`InterventionLabel[2]`) ~ "AP163",
    TRUE ~ "Other"
  ))

##### Get names of each cohort as drug, dose, unit ####
data <- data %>% 
  mutate(TreatmentLabel1 = case_when(
    grepl("Intervention", data$Treatment1Type) ~ paste0(drugname1, ", ", `Dose of treatment used:[1]`, ' ', `Measurement unit of treatment dose:[1]`),
    (grepl("Positive control", data$Treatment1Type) & is.na(`Dose of positive control treatment?[1]`)) ~ paste0(drugname1, ", ", `Dose of treatment used:[1]`),
    (grepl("Positive control", data$Treatment1Type) & !is.na(`Dose of positive control treatment?[1]`)) ~ paste0(drugname1, ", ", `Dose of positive control treatment?[1]`)
  ))

data <- data %>% 
  mutate(TreatmentLabel2 = case_when(
    grepl("Intervention", data$Treatment2Type) ~ paste0(drugname2, ", ", `Dose of treatment used:[2]`, ' ', `Measurement unit of treatment dose:[2]`),
    (grepl("Positive control", data$Treatment2Type) & is.na(`Dose of positive control treatment?[2]`)) ~ paste0(drugname2, ", ", `Dose of treatment used:[2]`),
    (grepl("Positive control", data$Treatment2Type) & !is.na(`Dose of positive control treatment?[2]`)) ~ paste0(drugname2, ", ", `Dose of positive control treatment?[2]`)
  ))

data$TreatmentLabel1 <- gsub("miligrams \\(mg\\) per kg", "mg/kg", data$TreatmentLabel1)
data$TreatmentLabel2 <- gsub("miligrams \\(mg\\) per kg", "mg/kg", data$TreatmentLabel2)
data$TreatmentLabel1 <- gsub("Other, NA NA", "", data$TreatmentLabel1)
data$TreatmentLabel2 <- gsub("Other, NA NA", "", data$TreatmentLabel2)

data <- data %>% 
  mutate(
    TreatmentLabel = case_when(
      !is.na(TreatmentLabel1) & !is.na(TreatmentLabel2) ~ paste(TreatmentLabel1, "&", TreatmentLabel2),
      !is.na(TreatmentLabel1) ~ TreatmentLabel1,
      !is.na(TreatmentLabel2) ~ TreatmentLabel2,
      TRUE ~ NA_character_
    )
  )


## may not be required
lab_n <- match('CohortLabel', names(data))
int_n <- match('InterventionLabel[1]', names(data))
n_n <- match('NumberOfAnimals', names(data))
m_n <- match('OutcomeResult', names(data))
v_n <- match('OutcomeError', names(data))

###for each cohort, label sham and control groups, along with TAAR1KO -ve control Contol for combination interventions (== positive control)
### first, we need to work out the attributes of each group - wgat types of comparisons will they allow?
# Assuming "df" is your data frame with columns "groupId" and "cohortType"

##### Count occurrences for each cohort type within each group #####
group_characteristics <- data %>%
  group_by(GroupID, CohortType) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = CohortType, values_from = count, values_fill = 0) %>%
  ungroup()
group_characteristics$total <- rowSums(select(group_characteristics, -GroupID))
group_characteristics$GroupID <- as.character(group_characteristics$GroupID)
data$GroupID <- as.character(data$GroupID)


##### Wrangling experiment type and outcome data #####

##### experimenting from here ####

cols2retain <- c('CohortId','CohortLabel','NumberOfAnimals','TreatmentLabel1',
                 'TreatmentLabel2','TreatmentLabel','CohortType','Treatment1Type','Treatment2Type',
                 'OutcomeResult','OutcomeError','GroupID','RoleOfCohort')
##### 3.1 T v C  with sham ####
## using same approach as for 3.3 and 3.4 to get consistency in final df ##
data_TvCs <- data.frame()

for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  cohset <- cohset %>% select(all_of(cols2retain))
  
  # Create combinations of interventions and positive controls
  combinations <- expand.grid(
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Simple intervention"]),
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Negative control"]),
    Sham = unique(cohset$CohortLabel [cohset$CohortType == "Sham"])
  )
  
  # Merge combinations with the original data
  cohset <- combinations %>%
    left_join(cohset %>% filter(CohortType == "Simple intervention"), by = c("Intervention" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Negative control"), by = c("Control" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Sham"), by = c("Sham" = "CohortLabel"))
  data_TvCs <- bind_rows(data_TvCs, cohset)
}
data_TvCs$Label <- data_TvCs$TreatmentLabel.x
data_TvCs$SortLabel <- "TvC"

##### 3.1 T v C  without sham ####

data_TvC <- data.frame()

for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  cohset <- cohset %>% select(all_of(cols2retain))
  
  # Create combinations of interventions and positive controls
  combinations <- expand.grid(
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Simple intervention"]),
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Negative control"])
  )
  
  # Merge combinations with the original data
  cohset <- combinations %>%
    left_join(cohset %>% filter(CohortType == "Simple intervention"), by = c("Intervention" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Negative control"), by = c("Control" = "CohortLabel"))
  data_TvC <- bind_rows(data_TvC, cohset)
}
data_TvC$Label <- data_TvC$TreatmentLabel.x
data_TvC$SortLabel <- "TvC"

##### 3.2 A v C with sham #####
data_AvCs <- data.frame()

for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  cohset <- cohset %>% select(all_of(cols2retain))
  
  # Create combinations of interventions and positive controls
  combinations <- expand.grid(
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Positive control"]),
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Negative control"]),
    Sham = unique(cohset$CohortLabel [cohset$CohortType == "Sham"])
  )
  
  # Merge combinations with the original data
  cohset <- combinations %>%
    left_join(cohset %>% filter(CohortType == "Positive control"), by = c("Intervention" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Negative control"), by = c("Control" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Sham"), by = c("Sham" = "CohortLabel"))
  
  data_AvCs <- bind_rows(data_AvCs, cohset)
}
data_AvCs$Label <- data_AvCs$TreatmentLabel.x
data_AvCs$SortLabel <- "AvC"

##### 3.2 A v C without sham #####
data_AvC <- data.frame()

for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  cohset <- cohset %>% select(all_of(cols2retain))
  
  # Create combinations of interventions and positive controls
  combinations <- expand.grid(
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Positive control"]),
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Negative control"])
  )
  
  # Merge combinations with the original data
  cohset <- combinations %>%
    left_join(cohset %>% filter(CohortType == "Positive control"), by = c("Intervention" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Negative control"), by = c("Control" = "CohortLabel"))
  
  data_AvC <- bind_rows(data_AvC, cohset)
}
data_AvC$Label <- data_AvC$TreatmentLabel.x
data_AvC$SortLabel <- "AvC"

##### 3.3 T v A - with sham - possibility of multiple control (A) conditions #####
data_TvAs <- data.frame()

for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  cohset <- cohset %>% select(all_of(cols2retain))

  # Create combinations of interventions and positive controls
  combinations <- expand.grid(
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Simple intervention"]),
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Positive control"]),
    Sham = unique(cohset$CohortLabel [cohset$CohortType == "Sham"])
  )
  
  # Merge combinations with the original data
  cohset <- combinations %>%
    left_join(cohset %>% filter(CohortType == "Simple intervention"), by = c("Intervention" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Positive control"), by = c("Control" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Sham"), by = c("Sham" = "CohortLabel"))
  
  data_TvAs <- bind_rows(data_TvAs, cohset)
  }
  data_TvAs$Label <- paste0(data_TvAs$TreatmentLabel.x, " v. ", data_TvAs$TreatmentLabel.y)
  data_TvAs$SortLabel <- "TvA"
  
  ##### 3.3 T v A without sham - possibility of multiple control (A) conditions #####
  data_TvA <- data.frame()
  
  for (i in 1:nrow(group_characteristics)) {
    group <- group_characteristics[i, 1]
    cohset <- data[data$GroupID == as.character(group), ]  
    cohset <- cohset %>% select(all_of(cols2retain))
    
    # Create combinations of interventions and positive controls
    combinations <- expand.grid(
      Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Simple intervention"]),
      Control = unique(cohset$CohortLabel [cohset$CohortType == "Positive control"])
    )
    
    # Merge combinations with the original data
    cohset <- combinations %>%
      left_join(cohset %>% filter(CohortType == "Simple intervention"), by = c("Intervention" = "CohortLabel")) %>%
      left_join(cohset %>% filter(CohortType == "Positive control"), by = c("Control" = "CohortLabel"))
    
    data_TvA <- bind_rows(data_TvA, cohset)
  }
  data_TvA$Label <- paste0(data_TvA$TreatmentLabel.x, " v. ", data_TvA$TreatmentLabel.y)
  data_TvA$SortLabel <- "TvA"
  
  
  ##### 3.4 TA v A - with sham - possibility of multiple control (A) conditions #####
  data_TAvAs <- data.frame()

for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  cohset <- cohset %>% select(all_of(cols2retain))
  
  # Create combinations of interventions and positive controls
  combinations <- expand.grid(
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Combination intervention"]),
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Positive control"]),
    Sham = unique(cohset$CohortLabel [cohset$CohortType == "Sham"])
  )
  
  # Merge combinations with the original data
  cohset <- combinations %>%
    left_join(cohset %>% filter(CohortType == "Combination intervention"), by = c("Intervention" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Positive control"), by = c("Control" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Sham"), by = c("Sham" = "CohortLabel"))
  
  data_TAvAs <- bind_rows(data_TAvAs, cohset)
}

diag1 <- subset(data_TAvAs, data_TAvAs$TreatmentLabel1.x == data_TAvAs$TreatmentLabel1.y)
diag2 <- subset(data_TAvAs, data_TAvAs$TreatmentLabel2.x == data_TAvAs$TreatmentLabel1.y)
diag3 <- subset(data_TAvAs, data_TAvAs$TreatmentLable1.x == data_TAvAs$TreatmentLabel2.y)
diag4 <- subset(data_TAvAs, data_TAvAs$TreatmentLabel2.x == data_TAvAs$TreatmentLabel2.y)
data_TAvAs <- rbind(diag1, diag2, diag3, diag4)
data_TAvAs$Label <- paste0(data_TAvAs$TreatmentLabel.x, " v. ", data_TAvAs$TreatmentLabel.y)
data_TAvAs$SortLabel <- "TAvA"
### ends

##### 3.4 TA v A - without sham - possibility of multiple control (A) conditions #####
data_TAvA <- data.frame()

for (i in 1:nrow(group_characteristics)) {
  group <- group_characteristics[i, 1]
  cohset <- data[data$GroupID == as.character(group), ]  
  cohset <- cohset %>% select(all_of(cols2retain))
  
  # Create combinations of interventions and positive controls
  combinations <- expand.grid(
    Intervention = unique(cohset$CohortLabel [cohset$CohortType == "Combination intervention"]),
    Control = unique(cohset$CohortLabel [cohset$CohortType == "Positive control"])
  )
  
  # Merge combinations with the original data
  cohset <- combinations %>%
    left_join(cohset %>% filter(CohortType == "Combination intervention"), by = c("Intervention" = "CohortLabel")) %>%
    left_join(cohset %>% filter(CohortType == "Positive control"), by = c("Control" = "CohortLabel"))
  
  data_TAvA <- bind_rows(data_TAvA, cohset)
}

diag1 <- subset(data_TAvA, data_TAvA$TreatmentLabel1.x == data_TAvA$TreatmentLabel1.y)
diag2 <- subset(data_TAvA, data_TAvA$TreatmentLabel2.x == data_TAvA$TreatmentLabel1.y)
diag3 <- subset(data_TAvA, data_TAvA$TreatmentLable1.x == data_TAvA$TreatmentLabel2.y)
diag4 <- subset(data_TAvA, data_TAvA$TreatmentLabel2.x == data_TAvA$TreatmentLabel2.y)
data_TAvA <- rbind(diag1, diag2, diag3, diag4)
data_TAvA$Label <- paste0(data_TAvA$TreatmentLabel.x, " v. ", data_TAvA$TreatmentLabel.y)
data_TAvA$SortLabel <- "TAvA"
### ends

##### combine to a single df and merge back to main #####
data_sub <- subset(data_TvC, !data_TvC$GroupID.x %in% data_TvCs$GroupID.x)
data_TvCc <- bind_rows(data_TvCs, data_sub)

data_sub <- subset(data_AvC, !data_AvC$GroupID.x %in% data_AvCs$GroupID.x)
data_AvCc <- bind_rows(data_AvCs, data_sub)

data_sub <- subset(data_TvA, !data_TvA$GroupID.x %in% data_TvAs$GroupID.x)
data_TvAc <- bind_rows(data_TvAs, data_sub)

data_sub <- subset(data_TAvA, !data_TAvA$GroupID.x %in% data_TAvAs$GroupID.x)
data_TAvAc <- bind_rows(data_TAvAs, data_sub)

data2 <- rbind(data_TvAc, data_AvCc, data_TvAc, data_TAvAc)

GroupIdlist <- unique(data2$GroupID.x)

data$OutcomeResult.x <-data$OutcomeResult
data$GroupID.x <- data$GroupID
merged_data1 <- merge(data, data_TvCc, by = c("GroupID.x","OutcomeResult.x" ), all.x = TRUE)
merged_data2 <- merge(data, data_AvCc, by = c("GroupID.x","OutcomeResult.x" ), all.x = TRUE)
merged_data3 <- merge(data, data_TvAc, by = c("GroupID.x","OutcomeResult.x" ), all.x = TRUE)
merged_data4 <- merge(data, data_TAvAc, by = c("GroupID.x","OutcomeResult.x" ), all.x = TRUE)

data1 <- rbind(merged_data1,merged_data2,merged_data3,merged_data4)
data1 <- subset(data1, !is.na(data1$SortLabel))

##### get standard placement - Single I #####
data1$F_C_L <- data1$Control
data1$F_C_n <- as.numeric(data1$NumberOfAnimals.y)
data1$F_C_m <- as.numeric(data1$OutcomeResult.y)
data1$F_C_v <- as.numeric(data1$OutcomeError.y)

data1$F_T_L <- data1$Intervention
data1$F_T_n <- as.numeric(data1$NumberOfAnimals.x.1)
data1$F_T_m <- as.numeric(data1$OutcomeResult.x.1)
data1$F_T_v <- as.numeric(data1$OutcomeError.x.1)

data1$F_S_L <- data1$Sham
data1$F_S_n <- as.numeric(data1$NumberOfAnimals.y.1)
data1$F_S_m <- as.numeric(data1$OutcomeResult.y.1)
data1$F_S_v <- as.numeric(data1$OutcomeError.y.1)

data_all_F <- data1

savename_all <- paste0(LSR,'data_all_',Sys.Date(),'.csv')

write_csv(data_all_F, savename_all)

###### Calculate effect size for simple interventions #####

# Data are ready to calculate effect sizes with numerical data for one comparison (test/control) per line

# Note: current code doesn't have cohort level questions split into treatment and control annotations per line
# Not an issue currently since control and treatment cohorts have had the same characteristics (sex, strain, etc) so far, but probably the subject of future edits (e.g. using pivot_wider function)



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

# Replace unit of measurements for drugs with abbreviations 
df <- df %>% 
  mutate(`Measurement unit of treatment dose:[1]` = str_replace_all(`Measurement unit of treatment dose:[1]`, "miligrams \\(mg\\) per kg", "mg/kg"),
         `Measurement unit of treatment dose:[2]` = str_replace_all(`Measurement unit of treatment dose:[2]`, "miligrams \\(mg\\) per kg", "mg/kg")) 

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
                          DrugName == "RO5256390" ~ 7.74,
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
                              DrugName == "RO5263397" ~ "TAAR1 partial agonist",
                              DrugName == "RO5256390" ~ "TAAR1 full agonist",
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
                                 DrugName == "RO5263397" ~ "High",
                                 DrugName == "RO5256390" ~ "High",
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
                                    ARRIVEScore > 19 ~ "F: > 20 criteria met")) 



# SAVE FILE
savefile_output <- paste0(LSR,'_','clean_data_',Sys.Date(),'.csv')
write.csv(df, savefile_output, row.names = FALSE)

