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

## Reading in the Data -----------------------------------------------------

# To make our analysis reproducible, we download the 2010 Global Burden of Disease data 
# straight from the source using the `readData()` function.

url1 <- "http://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_USA_GBD_2010_RESULTS_1990_2010_BY_CAUSE_Y2013M08D29.CSV"
cause <- readData(url1) %>%
    preprocessGBD()

url2 <- "http://ghdx.healthdata.org/sites/default/files/record-attached-files/IHME_USA_GBD_2010_RESULTS_1990_2010_BY_RISK_UPDATED_Y2013M11D21.CSV"
risk <- readData(url2) %>%
    rename(cause_name = risk_name) %>%
    preprocessGBD()

## Next, we read in the mortality, population, and prevalence data provided by NYCDOHMH.
mortality <- read.csv("data/2013_nyc_mortality.csv", stringsAsFactors=FALSE)
population <- read.csv("data/2013_nyc_population.csv", stringsAsFactors=FALSE)
prevalence <- read.csv("data/2013_nyc_prevalence.csv", stringsAsFactors=FALSE)


## Data Preparation --------------------------------------------------------

## We pre-process the national YLD/YLL rates by substituting values for `cause_name` 
## in order to match the indices of the other datasets. This will allow us to merge 
## datasets using `cause_name` as the key. We also write out the resulting dataset for inspection.

nationalRates <- rbind(cause, risk) %>%
    ungroup() %>%
    mutate(sex = ifelse(sex == "Females", "Female", "Male")) %>%
    mutate(cause_name = ifelse(cause_name == "Road injury", "Motor vehicle accidents", cause_name),
           cause_name = ifelse(cause_name == "Trachea, bronchus, and lung cancers", "Lung cancer", cause_name)) %>% 
    arrange(cause_name)
write.csv(nationalRates, "results/national_yldyll_rates.csv")

## Next, we pre-process the NYC mortality and calculate the YLLs for each disease by age, sex, and race. 
## For the analysis, we only use YLLs stratified by age and sex.

nycYLL <- calculateYLL(mortality)
write.csv(nycYLL, "results/nyc_yll_by_age_sex_race.csv")

nycYLL %<>%
    group_by(cause_name, sex, ageGroup) %>%
    summarize(yll = sum(yll))
write.csv(nycYLL, "results/nyc_yll_by_age_sex.csv")

## We calculate YLDs for each condition using NYC prevalence data, which also contains the associated disability 
## weights for each disease. To capture the level of uncertainty around disability weights, 
## we include the upper and lower bounds of the resulting YLDs in the output.

nycYLD <- calculatePrevalenceYLD(prevalence)
write.csv(nycYLD, "results/nyc_yld_by_age_sex.csv")

nycYLD %<>%
    group_by(cause_name, sex) %>%
    summarize(yld = sum(yld, na.rm=TRUE),
              yld_upper = sum(yld_upper, na.rm=TRUE),
              yld_lower = sum(yld_lower, na.rm=TRUE))
write.csv(nycYLD, "results/nyc_yld_by_sex.csv")

## DALY Estimation ------------------------------------------------------------------------

### Michaud YLD Approach

## create a search index
disease <- unique(c(nycYLL$cause_name, nycYLD$cause_name))
drug <- c("Amphetamine", "Heroin", "Cocaine", "Cannabis")
mental <- c("Major depressive disorder", "Anxiety", "Bipolar")
index <- unique(c(disease, drug, mental))

## This search index is then fed through the `calculateDALY` workhorse function to estimate DALYs 
## for each disease condition. The result is a `data.frame` object containing the following columns: `
## cause_name`, `sex`, `yll`, `yld`, `yld_upper`, `yld_lower`, `daly`, `daly_upper`, `daly_lower`.

michaudDALY <- lapply(index, calculateDALY, population, nycYLL=nycYLL, nationalRates=nationalRates)
michaudDALY <- do.call(rbind.fill, michaudDALY)
write.csv(michaudDALY, "results/nyc_daly_michaud.csv")

### Prevalence-Based YLD Approach

prevalenceDALY <- lapply(index, calculateDALY, population, nycYLL=nycYLL, nycYLD=nycYLD)
prevalenceDALY <- do.call(rbind.fill, prevalenceDALY)
write.csv(prevalenceDALY, "results/nyc_daly_prevalence.csv")

## Output --------------------------------------------------

## Michaud YLD Approach ----------------------------------

### 2013 NYC DALY Estimates, Total

michaudTotal <- segmentDALY(michaudDALY, strata="total")
michaudTotal

plotDALY(michaudTotal, "Leading Causes of DALYs, NYC 2013")
plotDALY(michaudTotal, "Leading Causes of DALYs, NYC 2013", stackedBar=TRUE)


### 2013 NYC DALY Estimates, Male

michaudMale <- segmentDALY(michaudDALY, strata="male")
michaudMale

plotDALY(michaudMale, "Leading Causes of DALYs in Males, NYC 2013")
plotDALY(michaudMale, "Leading Causes of DALYs in Males, NYC 2013", stackedBar=TRUE)

### 2013 NYC DALY Estimates, Female

michaudFemale <- segmentDALY(michaudDALY, strata="female")
michaudFemale

plotDALY(michaudFemale, "Leading Causes of DALYs in Females, NYC 2013")
plotDALY(michaudFemale, "Leading Causes of DALYs in Females, NYC 2013", stackedBar=TRUE)

## Prevalence-Based YLD Approach ---------------------------------------

### 2013 NYC DALY Estimates, Total


prevalenceTotal <- segmentDALY(prevalenceDALY, strata="total")
prevalenceTotal

plotDALY(prevalenceTotal, "Leading Causes of DALYs, NYC 2013")
plotDALY(prevalenceTotal, "Leading Causes of DALYs, NYC 2013", stackedBar=TRUE)


### 2013 NYC DALY Estimates, Male


prevalenceMale <- segmentDALY(prevalenceDALY, strata="male")
prevalenceMale

plotDALY(prevalenceMale, "Leading Causes of DALYs in Males, NYC 2013")
plotDALY(prevalenceMale, "Leading Causes of DALYs in Males, NYC 2013", stackedBar=TRUE)

### 2013 NYC DALY Estimates, Female

prevalenceFemale <- segmentDALY(prevalenceDALY, strata="female")
prevalenceFemale

plotDALY(prevalenceFemale, "Leading Causes of DALYs in Females, NYC 2013")
plotDALY(prevalenceFemale, "Leading Causes of DALYs in Females, NYC 2013", stackedBar=TRUE)


## Michaud YLDs vs. Prevalence-Based YLDs: Side-by-Side Comparison -------------------------

### Total
multiplot(plotDALY(michaudTotal, "Michaud YLDs"), plotDALY(prevalenceTotal, "Prevalence-Based YLDs"))

### Male
multiplot(plotDALY(michaudMale, "Michaud YLDs"), plotDALY(prevalenceMale, "Prevalence-Based YLDs"))

### Female
multiplot(plotDALY(michaudFemale, "Michaud YLDs"), plotDALY(prevalenceFemale, "Prevalence-Based YLDs"))

## Disease Conditions with Small Sample Sizes
prevalence[prevalence$small_sample == "yes", c("cause_name", "sequlae", "sex", "age")]

## END SCRIPT -------------------------------------------------------------
