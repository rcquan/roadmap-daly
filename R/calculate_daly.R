#############################
# Ryan Quan
# Columbia University
# GRAPH
# DOHMH Roadmap, Piece 1
# 2015-04-11
# rcq2102@columbia.edu
#
# The following script defines methods
# to calculate DALYs 
#############################

library("plyr")
library("dplyr")
library("reshape2")
library("magrittr")

# ------------------------------------------
calculateDALY <- function(diseaseName, population, nycYLL, nycYLD=NULL, nationalRates=NULL) {
    ## workhorse function to calculate DALY scores for specified disease using either
    ## prevalence-based YLD estimates or the Michaud approach using national YLD/YLL rates
    diseaseYLL <- subsetDataByDisease(diseaseName, nycYLL)
    if (!is.null(nycYLD) & !is.null(nationalRates)) {
        stop("You cannot provide values to both nycYLD and nationalRates parameters.")
        
    } else if (!is.null(nycYLD)) {
        nycYLD <- subsetDataByDisease(diseaseName, nycYLD)
        dalys <- calculatePrevalenceDALY(diseaseName, nycYLL, nycYLD)
        return(dalys)
        
    } else if (!is.null(nationalRates)) {
        ## subset datasets for specified disease
        diseaseRates <- subsetDataByDisease(diseaseName, nationalRates)
        ## if disease not found in gbdData, return YLL data as DALYs
        if (nrow(diseaseRates) == 0) {
            dalys <- diseaseYLL %>%
                group_by(cause_name, sex) %>%
                summarize(yll = sum(yll),
                          daly = sum(yll))
            return(dalys)
        }
        ## compute national YLD:YLL ratio and join to NYC YLL and population data by age, sex
        dalys <- diseaseRates %>% 
            ## compute national YLD:YLL ratio
            mutate(yldyll_ratio_mean = yld_nm_mean / yll_nm_mean,
                   yldyll_ratio_upper = yld_nm_upper / yll_nm_mean,
                   yldyll_ratio_lower = yld_nm_lower / yll_nm_mean) %>%
            # join tables
            join(population, by=c("ageGroup", "sex")) %>%
            join(diseaseYLL, by=c("cause_name", "ageGroup", "sex")) %>%
            ## estimate YLDs using Michaud logic
            mutate(yld = calculateMichaudYLD(yldyll_ratio_mean, yldyll_ratio_mean, yld_rt_mean, population, yll),
                   yld_upper = calculateMichaudYLD(yldyll_ratio_mean, yldyll_ratio_upper, yld_rt_upper, population, yll),
                   yld_lower = calculateMichaudYLD(yldyll_ratio_mean, yldyll_ratio_lower, yld_rt_lower, population, yll)) %>%
            ## collapse age groups
            group_by(cause_name, sex) %>%
            summarise_each(funs(sum(., na.rm=TRUE)), -c(cause_name, sex, ageGroup)) %>%
            ## calculate DALY estimates with lower and upper bounds
            mutate(daly = yll + yld,
                   daly_upper = yll + yld_upper,
                   daly_lower = yll + yld_lower) %>%
            select(cause_name, sex, yll, yld, yld_upper, yld_lower, daly, daly_upper, daly_lower)
        return(dalys)
    }
}
# ------------------------------------------
calculatePrevalenceDALY <- function(diseaseName, nycYLL, nycYLD) {
    ## calculates DALYs using prevalence-based YLDs from the 2010 GBD study
    ## Args:
    ##      diseaseName: chr. The disease of interest.
    ##      nycYLL: data.frame. New York City YLL estimates
    ##      nycYLD: data.frame. New York City YLD estimates 
    ## Returns:
    ##      dalys: data.frame. New York City DALY estimates
    diseaseYLL <- subsetDataByDisease(diseaseName, nycYLL)
    nycYLD <- subsetDataByDisease(diseaseName, nycYLD)
    dalys <- diseaseYLL %>%
        group_by(cause_name, sex) %>%
        summarize(yll = sum(yll)) %>%
        join(nycYLD, c("cause_name", "sex"), type = "right") %>%
        ungroup() %>%
        mutate(daly = ifelse(is.na(yll), 0 + yld, yll + yld),
               daly_upper = ifelse(is.na(yll), 0 + yld_upper, yll + yld_upper),
               daly_lower = ifelse(is.na(yll), 0 + yld_lower, yll + yld_lower))
    return(dalys)
}
# ------------------------------------------
getDiseaseIndex <- function(diseaseName, data) {
    ## searches disease index and returns indices of the first match
    ## Args:
    ##      diseaseName: string vector denoting diseases of interest
    ##      data: data.frame to be searched
    ## Returns:
    ##      indices of the first string match
    index <- grep(diseaseName, data$cause_name)
    pattern <- unique(data$cause_name[index])[1]
    return(which(data$cause_name == pattern))
}

# ------------------------------------------
subsetDataByDisease <- function(diseaseName, data) {
    ## subsets data frame from first string match
    index <- getDiseaseIndex(diseaseName, data)
    return(data[index, ])
}

# ------------------------------------------
calculateMichaudYLD <- function(checkRatio, yldyllRatio, nationalYLD, nycPop, nycYLL) {
    ## calculates YLDs based on the 2006 Michaud study
    ## Args:
    ##      checkRatio: numeric. National YLD:YLL ratio to check if > 10 or < 10
    ##      yldyllRatio: numeric. National YLD:YLL ratio to evaluate
    ##      nationalYLD: numeric. National YLD rate
    ##      nycPop: numeric. NYC Population
    ##      nycYLL: numeric. NYC YLL
    ## Returns:
    ##      nycYLD: New York City YLD estimate
    nycYLDLogic <- (checkRatio >= 10 | is.na(checkRatio) | is.infinite(checkRatio) | is.na(nycYLL))
    nycYLD <- ifelse(nycYLDLogic, nationalYLD * (nycPop / 100000), yldyllRatio * nycYLL)
    return(nycYLD)
}

# ------------------------------------------
calculatePrevalenceYLD <- function(nycPrevalence) {
    ## calculates prevalence-based YLD estimates from 2010 GBD Study
    ## Args:
    ##      nycPrevalence: data.frame. NYC prevalence data with associated disability weights
    ## Returns:
    ##      nycYLD: data.frame. NYC YLD estimates.
    nycYLD <- nycPrevalence %>%
        mutate(yld = prevalence * dependence_rate * dw_estimate,
               yld_upper = prevalence * dependence_rate * dw_upper,
               yld_lower = prevalence * dependence_rate * dw_lower)
    return(nycYLD)
}

# ------------------------------------------
calculateYLL <- function(mortalityData) {
    ## calculates YLLs from mortality data
    nycYLL <- mortalityData %>%
        mutate(le = sle - mean_age,
               yll = mortality * (1 - exp((-0.03 * le))) / 0.03)
    return(nycYLL)
}

# ------------------------------------------
segmentDALY <- function(dalyObj, strata) {
    ## helper function to subset DALY data
    if (strata == "total") {
        dalyObj %>% group_by(cause_name) %>% summarise_each(funs(sum), -c(sex)) %>% arrange(desc(daly))
    } else if (strata == "male") {
        dalyObj %>% filter(sex == "Male") %>% arrange(desc(daly))
    } else if (strata == "female") {
        dalyObj %>% filter(sex == "Female") %>% arrange(desc(daly))
    }
}
