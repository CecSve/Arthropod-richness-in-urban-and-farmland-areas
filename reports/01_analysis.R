#### Load required libraries ################################################################
library(cowplot) # for visuals
library(ggplot2) # for visuals
library(ggpubr) # for visuals
library(scales) # for visuals
library(ggpmisc) # ggplot extensions, for visuals
library(grid) # for visuals
library(gridExtra) # for visuals
library(car) # Companion to Applied Regression
library(lme4) # Fit linear and generalized linear mixed-effects models
library(lmerTest) # Provides p-values in type I, II or III anova and summary tables for lmer model fits (cf. lme4) 
library(nlme) # Fit and compare Gaussian linear and nonlinear mixed-effects models
library(effects) # Graphical and tabular effect displays, e.g., of interactions, for various statistical models with linear predictors
library(MuMIn)#AIC, R2
library(multcomp) # generalized linear hypothesis, to compare between effects of land cover
library(sjPlot) # to plot nice tables
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(tidyverse) # for data wrangling
library(lubridate) # for date and time wrangling

library(ggcorrplot)
library(GGally)

#### Set colour scheme ################################################################

#landuseCols <- c("#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "darkgrey") # colour friendly, ordered by land cover 

landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest")
#landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest", "Unspecified") # if including unspecified/other category

### Load data ###################################################

# load Danish data
allInsects <- read_rds("data/data_richness_BINs.RDS")

# change the five land covers to be 0-100 instead of 0-1
allInsects_trans_landcovers <- allInsects
allInsects[,33:60] <- allInsects[,33:60]*100

#### Correlation checks #################

##### Richness variables ############

resFolders <- "data"

#read in data
df <- readRDS(paste(resFolders,"data_richness_BINs.RDS", sep="/"))

#subset to response metrics
df <- df %>%
  dplyr::select(contains("richness"),
                contains("shannon"),
                contains("Biomass"),
                contains("n_reads")) %>%
  select(-est_richness_model, -biomassUncertainty, -est_richness_lci, -est_richness_uci)

colnames(df)
colnames(df) <- c("Rarefied richness", "Shannon diversity", "Observed richness", "Est. richness", "Biomass", "# of reads")

#sep for Year and Time_band eventually
ggpairs(df)

#get and plot correlation matrix
corr <- round(cor(df, use="pairwise.complete.obs"), 2)

png(file="plots/richness_variable_correlations.png",
    width=1500, height=1200, res=300)

ggcorrplot(corr, hc.order = TRUE, type = "lower",lab = TRUE,
           outline.col = "white",
           ggtheme = ggplot2::theme_minimal,
           colors = c("#6D9EC1", "white", "#E46726"))
dev.off()

##### Land cover variables ############

# correlation plot for 1000 m buffer
someInsects <- allInsects[,c(85,54:55,57,59:60)]
colnames(someInsects)
colnames(someInsects) <- c("Stops", "Farmland", "Forest", "Grassland", "Urban", "Wetland")

ggpairs(someInsects)

#get and plot correlation matrix
corr <- round(cor(someInsects, use="pairwise.complete.obs"), 2)

png(file="plots/landcover_variable_correlations.png",
    width=1500, height=1200, res=300)

ggcorrplot(corr, hc.order = TRUE, type = "lower",lab = TRUE,
           outline.col = "white",
           ggtheme = ggplot2::theme_minimal,
           colors = c("#6D9EC1", "white", "#E46726"))
dev.off()

#### Linear Mixed Effects Model #################