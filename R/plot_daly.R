#############################
# Ryan Quan
# Columbia University
# GRAPH
# DOHMH Roadmap, Piece 1
# 2015-04-11
# rcq2102@columbia.edu
#
# The following script defines methods
# to plot DALYs.
#############################

library("ggplot2")
library("dplyr")
library("reshape2")
library("grid")
library("scales")
library('stringr')

## -------------------------------------
segmentDALY <- function(dalyObj, strata) {
    ## helper function to subset DALY data
    if (strata == "total") {
        dalyObj %>% group_by(cause_name) %>% 
            summarise_each(funs(sum(., na.rm=TRUE)), -c(sex)) %>% 
            arrange(desc(daly)) %>% 
            as.data.frame()
    } else if (strata == "male") {
        dalyObj %>% filter(sex == "Male") %>% arrange(desc(daly)) %>%
            select(-sex)
    } else if (strata == "female") {
        dalyObj %>% filter(sex == "Female") %>% arrange(desc(daly)) %>%
            select(-sex)
    }
}
## -------------------------------------

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}
## -------------------------------------
plotDALY <- function(data, title, stackedBar=FALSE) {
    ## plot function for DALY object
    if (stackedBar) {
        meltedData <- melt(data, id.vars="cause_name", measure.vars=c("yll", "yld"), value.name="daly")
        ggplot(meltedData, aes(x=reorder(cause_name, daly, FUN=sum, na.rm=TRUE), y=daly, fill=variable)) + 
            geom_bar(stat="identity") +
            ggtitle(title) +
            ylab("Disability-Adjusted Life Years (DALYs)") + xlab("Causes") +
            scale_y_continuous(breaks=seq(0, max(data$daly_upper, na.rm=TRUE), by=100000), labels=comma) +
            scale_fill_brewer() + 
            coord_flip() +
            theme_bw()
    } else {
        limits <- aes(ymin=daly_lower, ymax=daly_upper)
        ggplot(data, aes(x=reorder(cause_name, daly), y=daly)) + 
            geom_pointrange(limits) + 
            ggtitle(title) +
            ylab("Disability-Adjusted Life Years (DALYs)") + xlab("Causes") +
            scale_y_continuous(breaks=seq(0, max(data$daly_upper, na.rm=TRUE), by=100000), labels=comma) +
            coord_flip() +
            theme_bw()
    }
}
