LSR3 documents for TAAR1 agonists in psychosis

- Raw data is found in the /data folder
- /wrangling contains files required for data wrangling (data_wrangle_script.R which calls functions in wrangling_functions.R)
- Assuming the raw data and files requred for wrangling are downloaded with the markdown document, wrangling and meta-analysis can be performed entirely within the R markdown document
- 2 input files are called from the main RMA script and this should be amended for later downloads:
    1. Qualitative ... ##date## project id ## csv : a SyRF system generated file
    2. LSR3_prisma_##date##.csv : a user generated file describing the fate of citations identifed in systematic review
- the code is specific to the LSR3 project, based on the structure of that project in SyRF; but should be adaptable for other SyRF projects

nb the code requires orchaRd 2.0, available from https://daniel1noble.github.io/orchaRd/, rather than the current (Dec 2023) CRAN version

- This project is licensed under the terms of the Creative Commons Attribution 4.0 International license (CC-BY 4.0) (https://creativecommons.org/licenses/by/4.0/). 
