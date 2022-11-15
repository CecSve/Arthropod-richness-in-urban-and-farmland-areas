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
library(SciViews) # for ln() function

library(ggcorrplot)
library(GGally)

#### Set colour scheme ################################################################

landuseCols <- c("#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "darkgrey") # colour friendly, ordered by land cover 

landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest")
#landuseOrder <- c("Urban","Farmland","Grassland","Wetland","Forest", "Unspecified") # if including unspecified/other category

### Load data ###################################################

# load Danish data
allInsects <- read_rds("data/data_richness_BINs.RDS")
asvs <- read_rds("data/asvs_BIN_combo.RDS")

# change the five land covers to be 0-100 instead of 0-1
allInsects_trans_landcovers <- allInsects
allInsects[,33:60] <- allInsects[,33:60]*100

#### Calculate evenness and other variables ########

allInsects$evenness <- allInsects$richness_rarefied_shannon/ln(allInsects$obs_richness)

allInsects$sDiversity_1000 <- scale(allInsects$Diversity_1000)
allInsects$sUrban_1000 <- scale(allInsects$Urban_1000)
allInsects$sAgriculture_1000 <- scale(allInsects$Agriculture_1000)
allInsects$scnumberTime <- scale(allInsects$cnumberTime)
allInsects$scTL <- scale(allInsects$cTL)
allInsects$scyDay <- scale(allInsects$cyDay)
allInsects$meanDiv <- mean(allInsects$Diversity_1000)
allInsects$sdDiv <- sd(allInsects$Diversity_1000)

min(allInsects$richness_est)
max(allInsects$richness_est)

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
someInsects <- allInsects[,c(85,54:55,57,59:60, 86)]
colnames(someInsects)
colnames(someInsects) <- c("Stops", "Farmland", "Forest", "Grassland", "Urban", "Wetland", "Complexity")

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

###### richness ~ land cover ####

model_lc <- lmer(richness_est ~ 
                     Urban_1000 +
                     Agriculture_1000 + 
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects) 
tab_model(model_lc)
summary(model_lc)
anova(model_lc)
shapiro.test(resid(model_lc))
qqnorm(resid(model_lc)) # plot of the normal theoretical quantiles against the exponential data quantiles
qqline(resid(model_lc)) # the residuals should fall along this line if they are normally distributed
AICc(model_lc)
#check variance inflation factor
vif(model_lc)

tab_model(model_lc, collapse.ci = F, show.intercept = F, pred.labels = c("Urban (1000 m)", "Farmland (1000 m)", "Time band: midday vs. evening", "Year", "Potential stops", "Day of year", "Time within evening", "Time within midday"), file = "plots/full_model_richness.html", digits = 2) #, file = "plots/full_model_richness.html"

# a quick plot
sjPlot::plot_model(model_lc)

# calculate % change urban cover
-0.93/mean(allInsects$richness_est)*100 # -0.7846992

###### richness ~ land cover x heterogeneity ####

# add three levels of heterogeneity
tertile_limits <- quantile(allInsects$Diversity_1000, seq(0, 1, 1/3), na.rm = TRUE)
allInsects$landscape_complexity <-
  cut(
    allInsects$Diversity_1000,
    tertile_limits,
    c('Low', 'Medium', 'High'),
    include.lowest = TRUE
  )


# focus on urban cover as this is the one with an effect
model_lc_h <- lmer(richness_est ~ 
                     (Urban_1000*Diversity_1000) +
                     #(Agriculture_1000*Diversity_1000) + 
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects) 

summary(model_lc_h)
anova(model_lc_h)
shapiro.test(resid(model_lc_h))
qqnorm(resid(model_lc_h)) # plot of the normal theoretical quantiles against the exponential data quantiles
qqline(resid(model_lc_h)) # the residuals should fall along this line if they are normally distributed
AICc(model_lc_h)
#check variance inflation factor
vif(model_lc_h)

tab_model(model_lc_h)
tab_model(model_lc_h, collapse.ci = F, show.intercept = F, pred.labels = c("Urban (1000 m)", "Heterogeneity (1000 m)", "Time band: midday vs. evening", "Year", "Potential stops", "Day of year", "Urban x Heterogeneity", "Time within midday", "Time within evening"), file = "plots/model_richness_landcover_heterogeneity.html", digits = 2) #, file = "plots/full_model_richness.html"

# a quick plot
sjPlot::plot_model(model_lc_h)

# calculate % change at scaled urban cover
-0.84/mean(allInsects$richness_est)*100 # -13.65208
15.36/mean(allInsects$richness_est)*100 # 5

# calculate with the scaled variables and manually insert into table
model_lc_h_scaled <- lmer(richness_est ~ 
                          (sUrban_1000*sDiversity_1000) +
                          #(Agriculture_1000*Diversity_1000) + 
                          Time_band + 
                          Year +
                          Time_band:scnumberTime + scTL + scyDay + 
                          (1|RouteID_JB) + (1|PID), data=allInsects)  

tab_model(model_lc_h_scaled)

# "Biomass (mg)", "mean DNA conc. (ng/ul)",

###### biomass ~ land cover ####

#full and final model - the effect of each land cover on estimated flying insect biomass. Notice the different buffer sizes from biomass!
model_b_lc <- lmer(log(totalBiomass_mg+1) ~ 
                     Urban_1000 +
                     Agriculture_1000 + 
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects) # data=subset(allInsects, c(Open.uncultivated.land_1000 < 20, Wetland_50 < 40))

tab_model(model_b_lc)
summary(model_b_lc)
AICc(model_b_lc)
#check variance inflation factor
vif(model_b_lc)

tab_model(model_b_lc, collapse.ci = F, show.intercept = F, pred.labels = c("Urban (1000 m)", "Farmland (1000 m)", "Time band: midday vs. evening", "Year", "Potential stops", "Day of year", "Time within midday", "Time within evening"), digits = 2, file = "plots/model_biomass_lc.html")

# a quick plot
sjPlot::plot_model(model_b_lc)

# scaled explanatory variables
model_b_lc_scaled <- lmer(log(totalBiomass_mg+1) ~ 
                     sUrban_1000 +
                     sAgriculture_1000 + 
                     Time_band + 
                     Year +
                     Time_band:scnumberTime + scTL + scyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects) 

tab_model(model_b_lc_scaled)

AICc(model_b_lc_scaled)
#check variance inflation factor
vif(model_b_lc_scaled)

# a quick plot
sjPlot::plot_model(model_b_lc_scaled)

# -46% biomass with 1% increase in urban cover (?!)

###### biomass ~ land cover x heterogeneity ####

# focus also only on urban cover 
model_b_lc_h <- lmer(log(totalBiomass_mg+1) ~ 
                     (Urban_1000*Diversity_1000) +
                     #(Agriculture_1000*Diversity_1000) + 
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects) 

tab_model(model_b_lc_h)
summary(model_b_lc_h)
anova(model_b_lc_h)
shapiro.test(resid(model_b_lc_h))
qqnorm(resid(model_b_lc_h)) # plot of the normal theoretical quantiles against the exponential data quantiles
qqline(resid(model_b_lc_h)) # the residuals should fall along this line if they are normally distributed
AICc(model_b_lc_h)
#check variance inflation factor
vif(model_b_lc_h)

tab_model(model_b_lc_h, collapse.ci = F, show.intercept = F, pred.labels = c("Urban (1000 m)", "Heterogeneity (1000 m)", "Time band: midday vs. evening", "Year", "Potential stops", "Day of year", "Urban x Heterogeneity", "Time within midday", "Time within evening"), file = "plots/model_biomass_landcover_heterogeneity.html", digits = 2) #, file = "plots/full_model_richness.html"

# a quick plot
sjPlot::plot_model(model_b_lc_h)

# with scaled variables

model_b_lc_h_scaled <- lmer(log(totalBiomass_mg+1) ~ 
                            (sUrban_1000*sDiversity_1000) +
                            #sAgriculture_1000 + 
                            Time_band + 
                            Year +
                            Time_band:scnumberTime + scTL + scyDay + 
                            (1|RouteID_JB) + (1|PID), data=allInsects)

tab_model(model_b_lc_h_scaled)

###### evenness ~ land cover ####

# Species evenness (J') is constrained between 0 and 1. The less evenness in communities between the species (and the presence of a dominant species), the lower J' is. And vice versa. 

# evenness
hist(allInsects$evenness)

model_e_lc <- lmer(evenness ~ 
                     Urban_1000 +
                     Agriculture_1000 + 
                     Time_band + 
                     Year +
                     Time_band:cnumberTime + cTL + cyDay + 
                     (1|RouteID_JB) + (1|PID), data=allInsects) # data=subset(allInsects, c(Open.uncultivated.land_1000 < 20, Wetland_50 < 40))

summary(model_e_lc)
AICc(model_e_lc)
#check variance inflation factor
vif(model_e_lc)
tab_model(model_e_lc, digits = 3)
tab_model(model_e_lc, collapse.ci = F, show.intercept = F, pred.labels = c("Urban (1000 m)", "Farmland (1000 m)", "Time band: midday vs. evening", "Year", "Potential stops", "Day of year", "Time within midday", "Time within evening"), digits = 3, file = "plots/model_evenness_lc.html")

# a quick plot
sjPlot::plot_model(model_e_lc)

# calculate % change at scaled urban cover
0.001/mean(allInsects$evenness, na.rm = T)*100 # 0.1833027 - one NA value in there

###### evenness ~ land cover x heterogeneity ####

# focus also only on urban cover 
model_e_lc_h <- lmer(evenness ~ 
                       (Urban_1000*Diversity_1000) +
                       #(Agriculture_1000*Diversity_1000) + 
                       Time_band + 
                       Year +
                       Time_band:cnumberTime + cTL + cyDay + 
                       (1|RouteID_JB) + (1|PID), data=allInsects) 

tab_model(model_e_lc_h, digits = 3, show.intercept = F)
summary(model_e_lc_h)
anova(model_e_lc_h)
shapiro.test(resid(model_e_lc_h))
qqnorm(resid(model_e_lc_h)) # plot of the normal theoretical quantiles against the exponential data quantiles
qqline(resid(model_e_lc_h)) # the residuals should fall along this line if they are normally distributed
AICc(model_e_lc_h)
#check variance inflation factor
vif(model_e_lc_h)

tab_model(model_e_lc_h, collapse.ci = F, show.intercept = F, pred.labels = c("Urban (1000 m)", "Heterogeneity (1000 m)", "Time band: midday vs. evening", "Year", "Potential stops", "Day of year", "Urban x Heterogeneity", "Time within midday", "Time within evening"), file = "plots/model_evenness_landcover_heterogeneity.html", digits = 3) #, file = "plots/full_model_richness.html"

# a quick plot
sjPlot::plot_model(model_e_lc_h)

# calculate % change at scaled urban cover
0.001/mean(allInsects$evenness, na.rm = T)*100 # 0.1833027 - one NA value in there

# scaled variables
model_e_lc_h_scaled <- lmer(evenness ~ 
                       (sUrban_1000*sDiversity_1000) +
                       #(Agriculture_1000*Diversity_1000) + 
                       Time_band +
                       Year +
                       Time_band:scnumberTime + scTL + scyDay +
                       (1 | RouteID_JB) + (1 | PID),
                     data = allInsects)

tab_model(model_e_lc_h_scaled, digits = 3, show.intercept = F)

### percent change between biomass and richness ####
(-0.78--0.03)/-0.03 # percent
(-0.93--0.03)/-0.03 # estimates

#### Effect plots for land cover effect ##########################

###### richness ~ land cover ##########

# extract effects
gls1.alleffects <- allEffects(model_lc)
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(model_lc)
plot(eall.lm1, lines=list(multiline=TRUE))
plot(predictorEffects(model_lc, ~ Agriculture_1000 + Urban_1000 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

# make data frames for each land cover
temp <- effectdata$Agriculture_1000
temp$landcover <- "Agriculture_1000"
farm <- temp %>% 
  dplyr::rename(
    propcover = Agriculture_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# urban
temp <- effectdata$Urban_1000
temp$landcover <- "Urban_1000"
urb <- temp %>% 
  dplyr::rename(
    propcover = Urban_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Timeband
temp <- effectdata$`Time_band:cnumberTime`
#temp$landcover <- "Time_band:cnumberTime"
time <- temp

# combine the data frames
r_lc <- rbind(urb, farm)

# Visualization
r_lc_effectplot <- r_lc %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_1000"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban cover",
      "Farmland cover"
    )
  ) +guides(colour=guide_legend(ncol=1, byrow = T)) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size =14),
    legend.position = "",
    legend.spacing.x = unit(1.0, 'cm'),
    legend.key.height = unit(1, 'cm')
  ) + scale_x_continuous(
    limits = c(0, 100),
    labels = function(x)
      paste0(x, "%")) + geom_ribbon(
        aes(
          ymin = fit-se,
          ymax = fit+se,
          group = landcover
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        #x = "Land cover extent",
        y = "Estimated richness",
        #subtitle = "A",
        colour = "Land cover type"
      ) + scale_fill_manual(values = landuseCols)

#cowplot::save_plot("plots/Fig_DK_effect_richness_landcover.tiff", r_lc_effectplot, base_width = 8, base_height = 5, dpi = 800)

#cowplot::save_plot("plots/Fig_DK_effect_richness_landcover.png", effectplot, base_width = 8, base_height = 5, dpi = 800)

###### evenness ~ land cover ##########

# extract effects
gls1.alleffects <- allEffects(model_e_lc)
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(model_e_lc)
plot(eall.lm1, lines=list(multiline=TRUE))
plot(predictorEffects(model_e_lc, ~ Agriculture_1000 + Urban_1000 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

# make data frames for each land cover
temp <- effectdata$Agriculture_1000
temp$landcover <- "Agriculture_1000"
farm <- temp %>% 
  dplyr::rename(
    propcover = Agriculture_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# urban
temp <- effectdata$Urban_1000
temp$landcover <- "Urban_1000"
urb <- temp %>% 
  dplyr::rename(
    propcover = Urban_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Timeband
temp <- effectdata$`Time_band:cnumberTime`
#temp$landcover <- "Time_band:cnumberTime"
time <- temp

# combine the data frames
e_lc <- rbind(urb, farm)

# Visualization
e_lc_effectplot <- e_lc %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_1000"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban cover",
      "Farmland cover"
    )
  ) +guides(colour=guide_legend(ncol=1, byrow = T)) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size =14),
    legend.position = "",
    legend.spacing.x = unit(1.0, 'cm'),
    legend.key.height = unit(1, 'cm')
  ) + scale_x_continuous(
    limits = c(0, 100),
    labels = function(x)
      paste0(x, "%")) + geom_ribbon(
        aes(
          ymin = fit-se,
          ymax = fit+se,
          group = landcover
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        #x = "Land cover extent",
        y = "Evenness",
        #subtitle = "B",
        colour = "Land cover type"
      ) + scale_fill_manual(values = landuseCols)

###### biomass ~ land cover ##########

# extract effects
gls1.alleffects <- allEffects(model_b_lc)
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(model_b_lc)
plot(eall.lm1, lines=list(multiline=TRUE))
plot(predictorEffects(model_b_lc, ~ Agriculture_1000 + Urban_1000 + cnumberTime, residuals = T), partial.residuals=list(smooth=TRUE, span=0.50, lty = "dashed"))

# make data frames for each land cover
temp <- effectdata$Agriculture_1000
temp$landcover <- "Agriculture_1000"
farm <- temp %>% 
  dplyr::rename(
    propcover = Agriculture_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# urban
temp <- effectdata$Urban_1000
temp$landcover <- "Urban_1000"
urb <- temp %>% 
  dplyr::rename(
    propcover = Urban_1000
  )%>% dplyr::select(landcover, propcover, fit, se, lower, upper)

# Timeband
temp <- effectdata$`Time_band:cnumberTime`
#temp$landcover <- "Time_band:cnumberTime"
time <- temp

# combine the data frames
b_lc <- rbind(urb, farm)

# Visualization
b_lc_effectplot <- b_lc %>% mutate(
  landcover = fct_relevel(
    landcover,
    "Urban_1000",
    "Agriculture_1000"
  )
) %>% ggplot(aes(x = propcover, y = fit, fill = landcover)) +
  geom_line(aes(color = landcover), size = 2) +
  scale_color_manual(
    values = landuseCols,
    labels = c(
      "Urban cover",
      "Farmland cover"
    )
  ) +guides(colour=guide_legend(nrow=1, byrow = T)) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size =10),
    legend.position = "",
    legend.spacing.x = unit(0.5, 'cm'),
    legend.key.height = unit(0.5, 'cm')
  ) + scale_x_continuous(
    limits = c(0, 100),
    labels = function(x)
      paste0(x, "%")) + geom_ribbon(
        aes(
          ymin = fit-se,
          ymax = fit+se,
          group = landcover
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        #x = "Land cover extent",
        y = "Estimated log(biomass mg+1)",
        #subtitle = "C",
        colour = "Land cover type"
      ) + scale_fill_manual(values = landuseCols)

##### arrange plots #####

effect_fig <-
  ggarrange(
    r_lc_effectplot + rremove("xlab"),
    ggarrange(
      e_lc_effectplot + rremove("xlab"),
      b_lc_effectplot + rremove("xlab"),
      ncol = 2,
      labels = c("B", "C")
    ),
    nrow = 2,
    labels = "A"
  ) 

fig_effect <- annotate_figure(effect_fig, bottom = textGrob("Land cover extent", gp = gpar(cex = 1.3)))

#cowplot::save_plot("plots/Fig_DK_effect_biomass_landcover.tiff", effectplot, base_width = 8, base_height = 5, dpi = 800)

cowplot::save_plot("plots/landcover_effectplot.png", fig_effect, base_width = 12, base_height = 8, dpi = 800)

#### Effect plots for land cover x complexity effect ##########################

###### richness ~ land cover ##########

#landcover_complex_col <- c("#944872", "#E0D75C", "#E084B7", "#6ED5E0", "#508D94")
urban_col <- c("#6B164F", "#AB247E", "#F734B6")

# extract effects
gls1.alleffects <- allEffects(model_lc_h, xlevels=list(Diversity_1000 = c(allInsects$meanDiv - allInsects$sdDiv, allInsects$meanDiv, allInsects$meanDiv + allInsects$sdDiv)))
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

plotdata <- effectdata$`Urban_1000:Diversity_1000`

# Visualization
r_lc_h_effectplot <- plotdata %>%
  ggplot(aes(x = Urban_1000, y = fit, group = Diversity_1000, fill = factor(Diversity_1000))) +
  geom_line(aes(color = factor(Diversity_1000)), size = 2) +
  scale_color_manual(values = urban_col, labels = c(
    "low complexity",
    "medium complexity", "high complexity")) + guides(colour=guide_legend(ncol=3, byrow = T)) + theme_minimal_grid() + theme(
    plot.subtitle = element_text(size = 20, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size =14),
    legend.position = "",
    legend.spacing.x = unit(1.0, 'cm'),
    legend.key.height = unit(1, 'cm'),
    legend.direction = "horizontal"
  ) + scale_x_continuous(
    limits = c(0, 100),
    labels = function(x)
      paste0(x, "%")) + geom_ribbon(
        aes(
          ymin = fit-se,
          ymax = fit+se,
          group = Diversity_1000
        ),
        linetype = 2,
        alpha = 0.2,
        show.legend = F
      ) + labs(
        #x = "Land cover extent",
        y = "Estimated richness",
        #subtitle = "A",
        colour = "Land cover type"
      ) + scale_fill_manual(values = urban_col)


#cowplot::save_plot("plots/Fig_DK_effect_richness_landcover.tiff", r_lc_effectplot, base_width = 8, base_height = 5, dpi = 800)

#cowplot::save_plot("plots/Fig_DK_effect_richness_landcover.png", effectplot, base_width = 8, base_height = 5, dpi = 800)

###### evenness ~ land cover ##########

# extract effects
gls1.alleffects <- allEffects(model_e_lc_h, xlevels=list(Diversity_1000 = c(allInsects$meanDiv - allInsects$sdDiv, allInsects$meanDiv, allInsects$meanDiv + allInsects$sdDiv)))
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

plotdata <- effectdata$`Urban_1000:Diversity_1000`

# Visualization
e_lc_h_effectplot <- plotdata %>%
  ggplot(aes(x = Urban_1000, y = fit, group = Diversity_1000, fill = factor(Diversity_1000))) +
  geom_line(aes(color = factor(Diversity_1000)), size = 2) +
  scale_color_manual(values = urban_col, labels = c(
    "low complexity",
    "medium complexity", "high complexity")) + guides(colour=guide_legend(ncol=3, byrow = T)) + theme_minimal_grid() + theme(
      plot.subtitle = element_text(size = 20, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size =14),
      legend.position = "",
      legend.spacing.x = unit(1.0, 'cm'),
      legend.key.height = unit(1, 'cm'),
      legend.direction = "horizontal"
    ) + scale_x_continuous(
      limits = c(0, 100),
      labels = function(x)
        paste0(x, "%")) + geom_ribbon(
          aes(
            ymin = fit-se,
            ymax = fit+se,
            group = Diversity_1000
          ),
          linetype = 2,
          alpha = 0.2,
          show.legend = F
        ) + labs(
          #x = "Land cover extent",
          y = "Evenness",
          #subtitle = "A",
          colour = "Land cover type"
        ) + scale_fill_manual(values = urban_col)

###### biomass ~ land cover ##########

# extract effects
gls1.alleffects <- allEffects(model_b_lc_h, xlevels=list(Diversity_1000 = c(allInsects$meanDiv - allInsects$sdDiv, allInsects$meanDiv, allInsects$meanDiv + allInsects$sdDiv)))
effectdata <- as.data.frame(gls1.alleffects, row.names=NULL, optional=TRUE)

eall.lm1 <- predictorEffects(model_b_lc_h)

plotdata <- effectdata$`Urban_1000:Diversity_1000`

# Visualization
b_lc_h_effectplot <- plotdata %>%
  ggplot(aes(x = Urban_1000, y = fit, group = Diversity_1000, fill = factor(Diversity_1000))) +
  geom_line(aes(color = factor(Diversity_1000)), size = 2) +
  scale_color_manual(values = urban_col, labels = c(
    "low complexity",
    "medium complexity", "high complexity")) + guides(colour=guide_legend(ncol=3, byrow = T)) + theme_minimal_grid() + theme(
      plot.subtitle = element_text(size = 20, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size =14),
      legend.position = "",
      legend.spacing.x = unit(1.0, 'cm'),
      legend.key.height = unit(1, 'cm')
    ) + scale_x_continuous(
      limits = c(0, 100),
      labels = function(x)
        paste0(x, "%")) + geom_ribbon(
          aes(
            ymin = fit-se,
            ymax = fit+se,
            group = Diversity_1000
          ),
          linetype = 2,
          alpha = 0.2,
          show.legend = F
        ) + labs(
          #x = "Land cover extent",
          y = "Estimated log(biomass mg+1)",
          #subtitle = "A",
          colour = "Land cover type"
        ) + scale_fill_manual(values = urban_col)

##### arrange plots #####

effect_fig <-
  ggarrange(
    r_lc_h_effectplot + rremove("xlab"),
    ggarrange(
      e_lc_h_effectplot + rremove("xlab"),
      b_lc_h_effectplot + rremove("xlab"),
      ncol = 2,
      labels = c("B", "C")
    ),
    nrow = 2,
    labels = "A",
    legend = "top",
    common.legend = T
  ) 

fig_effect <- annotate_figure(effect_fig, bottom = textGrob("Land cover extent", gp = gpar(cex = 1.3)))

#cowplot::save_plot("plots/Fig_DK_effect_biomass_landcover.tiff", effectplot, base_width = 8, base_height = 5, dpi = 800)

cowplot::save_plot("plots/landcover_heterogeneity_effectplot.png", fig_effect, base_width = 12, base_height = 8, dpi = 800)
