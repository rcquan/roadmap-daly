#############################
# Ryan Quan
# Columbia University
# GRAPH
# DOHMH Roadmap, Piece 1
# 2015-04-11
# rcq2102@columbia.edu
#
# The following script runs
# the DALY estimations for 2013 NYC data
#############################

source("R/preprocess.R")
source("R/calculate_daly.R")
source("R/plot_daly.R")
dir.create("results")
dir.create("results/sensitivity_analysis")
dir.create("data")
## Reading in the Data -----------------------------------------------------
url <- "http://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_USA_GBD_2010_RESULTS_1990_2010_BY_CAUSE_Y2013M08D29.CSV"
cause <- readData(url) %>%
    preprocessGBD()

url <- "http://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_USA_GBD_2010_RESULTS_1990_2010_BY_RISK_UPDATED_Y2013M11D21.CSV"
risk <- readData(url) %>%
    rename(cause_name = risk_name) %>%
    preprocessGBD()

mortality <- read.csv("data/2013_nyc_mortality.csv", stringsAsFactors=FALSE)
population <- read.csv("data/2013_nyc_population.csv", stringsAsFactors=FALSE)
prevalence <- read.csv("data/2013_nyc_prevalence.csv", stringsAsFactors=FALSE)

nationalRates <- rbind(cause, risk) %>%
    ungroup() %>%
    mutate(sex = ifelse(sex == "Females", "Female", "Male")) %>%
    mutate(cause_name = ifelse(cause_name == "Road injury", "Motor vehicle accidents", cause_name),
           cause_name = ifelse(cause_name == "Trachea, bronchus, and lung cancers", "Lung cancer", cause_name)) %>% 
    arrange(cause_name)
write.csv(nationalRates, "results/national_yldyll_rates.csv", row.names=FALSE)

## Preprocessing -----------------------------------------------
nycYLL <- calculateYLL(mortality)
write.csv(nycYLL, "results/nyc_yll_by_age_sex_race.csv", row.names=FALSE)

nycYLL %<>%
    group_by(cause_name, sex, ageGroup) %>%
    summarize(yll = sum(yll))
write.csv(nycYLL, "results/nyc_yll_by_age_sex.csv", row.names=FALSE)

nycYLD <- calculatePrevalenceYLD(prevalence)
write.csv(nycYLD, "results/nyc_yld_by_age_sex.csv", row.names=FALSE)

nycYLD %<>%
    group_by(cause_name, sex) %>%
    summarize(yld = sum(yld, na.rm=TRUE),
              yld_upper = sum(yld_upper, na.rm=TRUE),
              yld_lower = sum(yld_lower, na.rm=TRUE))
write.csv(nycYLD, "results/nyc_yld_by_sex.csv", row.names=FALSE)

## DALY Estimation -------------------------------------

### Michaud YLD Approach

## This section contains an implementation of the Michaud approach 
## described in the above methods section. We first create a search index 
## containing all the disease conditions of interest. 

## create a search index
disease <- unique(c(nycYLL$cause_name, nycYLD$cause_name))
drug <- "Cannabis"
mental <- c("Major depressive disorder", "Anxiety", "Bipolar")
index <- unique(c(disease, drug, mental))

## This search index is then fed through the `calculateDALY` workhorse function to 
## estimate DALYs for each disease condition. The result is a `data.frame` object containing the 
## following columns: `cause_name`, `sex`, `yll`, `yld`, `yld_upper`, `yld_lower`, `daly`, 
## `daly_upper`, `daly_lower`.

michaudDALY <- lapply(index, calculateDALY, population, nycYLL=nycYLL, nationalRates=nationalRates)
michaudDALY <- do.call(rbind.fill, michaudDALY)
write.csv(michaudDALY, "results/nyc_daly_michaud.csv", row.names=FALSE)

### Prevalence-Based YLD Approach --------------------------------------

## Similar to the section, we implement the prevalence-based YLD approach 
## here using the same search index.

prevalenceDALY <- lapply(index, calculateDALY, population, nycYLL=nycYLL, nycYLD=nycYLD)
prevalenceDALY <- do.call(rbind.fill, prevalenceDALY)
prevalenceDALY$cause_name <- sapply(prevalenceDALY$cause_name, renameDiseaseLabel, USE.NAMES=FALSE)
write.csv(prevalenceDALY, "results/nyc_daly_prevalence.csv", row.names=FALSE)

### Sensitivity Analysis ---------------------------------------------

sensitivity <- read.csv("data/2005_nyc_mortality.csv", stringsAsFactors=FALSE)

nycYLLSensitivity <- calculateYLL(sensitivity)
write.csv(nycYLLSensitivity, "results/sensitivity_analysis/nyc_yll_by_age_sex_race.csv", row.names=FALSE)

nycYLLSensitivity %<>%
    group_by(cause_name, sex, ageGroup) %>%
    summarize(yll = sum(yll))
write.csv(nycYLLSensitivity, "results/sensitivity_analysis/nyc_yll_by_age_sex.csv", row.names=FALSE)

disease <- unique(c(nycYLL$cause_name, nycYLD$cause_name))
drug <- "Cannabis"
mental <- c("Major depressive disorder", "Anxiety", "Bipolar")
index <- unique(c(disease, drug, mental))

michaudDALYSensitivity <- lapply(index, calculateDALY, population, nycYLL=nycYLL, nationalRates=nationalRates)
michaudDALYSensitivity <- do.call(rbind.fill, michaudDALYSensitivity)
write.csv(michaudDALYSensitivity, "results/sensitivity_analysis/nyc_daly_michaud.csv", row.names=FALSE)

## END SCRIPT -------------------------------------------------------------
