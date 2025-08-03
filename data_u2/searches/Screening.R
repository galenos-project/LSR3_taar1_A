### Data preparation

# Import SyRF screening decisions
tiab_screening_LSR3 <- read_csv("data_u2/Screening_data_-_2025_07_30_-_Long_format_-_c494b1ae-4cf4-4618-b91a-e69a2b815bdd_-_Investigators_Unblinded.csv")
ft_screening_LSR3 <- read_csv("data_u2/Annotation_data_-_2025_08_01_-_Long_format_-_c494b1ae-4cf4-4618-b91a-e69a2b815bdd_-_Investigators_Unblinded.csv")

# Get bibliographic details of studies included/excluded during tiab and ft screening
tiab_screening_LSR3 <- tiab_screening_LSR3 %>% 
  select(StudyId, Title, Authors, PublicationName, Year, Doi, SystematicSearchName, ScreeningStatus) %>%
  filter(!SystematicSearchName %in% c("dup", "dup2")) %>%
  distinct()

ft_screening_LSR3 <- ft_screening_LSR3 %>%
  select(StudyId, Title, Authors, PublicationName, Year, Doi, SystematicSearchName, Question, Answer, Comments) %>%
  filter(Question == "should have been excluded") %>%
  filter(!SystematicSearchName %in% c("dup", "dup2")) %>%
  group_by(StudyId) %>%
  arrange(is.na(Comments), .by_group = TRUE) %>%
  slice(1) %>%
  distinct()


### Data for all searches (original and updates)

# Included at tiab
tiab_screening_LSR3_included <- tiab_screening_LSR3 %>%
  filter(ScreeningStatus == "Included")

# Included after FT screening
ft_screening_LSR3_included <- ft_screening_LSR3 %>%
  filter(Answer == "False")

# Excluded based on FT screening
ft_screening_LSR3_excluded <- ft_screening_LSR3 %>%
  filter(Answer == "True")
  

### Data for update 1 (u1)

# Included at tiab
tiab_screening_u1_included <- tiab_screening_LSR3_included %>%
  filter(SystematicSearchName == "150224_SOLES_update")

# Included after FT screening
ft_screening_u1_included <- ft_screening_LSR3_included %>%
  filter(SystematicSearchName == "150224_SOLES_update")

# Excluded based on FT screening
ft_screening_u1_excluded <- ft_screening_LSR3_excluded %>%
  filter(SystematicSearchName == "150224_SOLES_update")

### Data for update 2 (u2)

# Included at tiab
tiab_screening_u2_included <- tiab_screening_LSR3_included %>%
  filter(SystematicSearchName == "LSR3_update_08112024")

# Included after FT screening
ft_screening_u2_included <- ft_screening_LSR3_included %>%
  filter(SystematicSearchName == "LSR3_update_08112024")

# Excluded based on FT screening
ft_screening_u2_excluded <- ft_screening_LSR3_excluded %>%
  filter(SystematicSearchName == "LSR3_update_08112024")


### Write .csv files for studies with different screening decisions

# all studies (original and updates)
write.csv(ft_screening_LSR3_included, "data_u2/searches/LSR3_all_included.csv", row.names = FALSE)
write.csv(ft_screening_LSR3_excluded, "data_u2/searches/LSR3_all_excluded.csv", row.names = FALSE)

# update 1
write.csv(ft_screening_u1_included, "data_u2/searches/LSR3_u1_included.csv", row.names = FALSE)
write.csv(ft_screening_u1_excluded, "data_u2/searches/LSR3_u1_excluded.csv", row.names = FALSE)

# update 2
write.csv(ft_screening_u2_included, "data_u2/searches/LSR3_u2_included.csv", row.names = FALSE)
write.csv(ft_screening_u2_excluded, "data_u2/searches/LSR3_u2_excluded.csv", row.names = FALSE)




