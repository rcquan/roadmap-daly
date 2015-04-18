#############################
# Ryan Quan
# Columbia University
# GRAPH
# DOHMH Roadmap, Piece 1
# 2015-04-11
# rcq2102@columbia.edu
#
# The following script defines methods
# to pre-process data required for DALY estimation
#############################

library("plyr")
library("dplyr")
library("reshape2")
library("magrittr")

# ------------------------------------------
readData <- function(url) {
    ## Reads CSV data from input URL string
    filename <- tail(unlist(strsplit(url, "/")), 1)
    filepath <- paste("data", "/", filename, sep="")
    if (!file.exists(filepath)) {
        download.file(url, filename, method="curl")
    }
    data <- read.csv(filepath, stringsAsFactors=FALSE)
    return(data)
}

# ------------------------------------------
assignAgeGroup <- function(ageVar) {
    ## logic for childhood, teenage, young adult, adult, and later-life age groups
    if (ageVar %in% c("Under 5 years", "5-14 years")) {
        return("00-14")
    } else if (ageVar %in% c("15-19 years", "20-24 years")) {
        return("15-24")
    } else if (ageVar %in% c("25-29 years", "30-34 years", "35-39 years", "40-44 years")) {
        return("25-44")
    } else if (ageVar %in% c("45-49 years", "50-54 years", "55-59 years", "60-64 years")) {
        return("45-64")
    } else if (ageVar %in% c("65-69 years", "70+ years")) {
        return("65+")
    } else {
        return("")
    }
}
# ------------------------------------------
addAgeGroup <- function(data, ageVar="age_name") {
    ## replaces age grouping in current data.frame to childhood, teenage, YA, adult, later-life
    ## Args:
    ##      data: data.frame object
    ##      ageVar: string denoting the column of ages to be replaced
    ## Returns:
    ##      data: data.frame object with new age groupings
    ageGroup <- vector(length=nrow(data))
    for (i in 1:nrow(data)) {
        ageGroup[i] <- assignAgeGroup(as.vector(data[i, ageVar]))
    }
    data$ageGroup <- ageGroup
    return(data)
}
# ------------------------------------------
preprocessGBD <- function(data) {
    ## extracts YLD and YLL rates from 2010 Global Burden of Disase data
    ## Args:
    ##      data: GBD dataset downloaded from the web
    ## Returns:
    ##      data: a pre-processed 2010 GBD dataset
    data %<>%
        ## filter out unnecessary variables
        select(-c(pc_mean, pc_upper, pc_lower)) %>%
        filter(year == 2010) %>%
        filter(sex %in% c("Females", "Males")) %>%
        ## extract only YLD and YLL rates
        filter(measure %in% c("yll", "yld")) %>%
        ## create long-form dataset
        melt(measure.vars=c("nm_mean", "nm_upper", "nm_lower", "rt_mean", "rt_upper", "rt_lower")) %>%
        ## create wide-form dataset with national YLD/YLL rates
        dcast(cause_name + age_name + sex ~ measure + variable, value.var="value") %>%
        ## age group manipulations
        addAgeGroup("age_name") %>%
        filter(ageGroup != "") %>%
        select(-age_name) %>%
        ## averaging YLD/YLL rates with respect to new age groupings
        group_by(cause_name, sex, ageGroup) %>%
        summarise_each(funs(mean))
    return(data)
}
