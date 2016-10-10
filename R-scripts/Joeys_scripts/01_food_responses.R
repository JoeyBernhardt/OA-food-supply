# OA-food-supply
# Joey Bernhardt
# Sept 10 2016

## To do: first, just do comparison with most conservative treatment, or closest to delta ph of 0.5

# load packages -----------------------------------------------------------

library(readr)
library(tidyr)
library(dplyr)
library(metafor)
library(Formula)
library(ggplot2)

# load data ---------------------------------------------------------------

growth <- read_csv("data-raw/growth_data.csv")


# clean data --------------------------------------------------------------

growth_clean <- growth %>% 
	filter(Food.supply != "Med")


# create effect sizes -----------------------------------------------------

effect_sizes <- escalc(Mean_Ambient/SD_Ambient ~ factor(Food.supply) | factor(Paper_no),data = growth_clean, weights = N_Ambient, measure = "ROM",slab = Paper_No)
effects_summary <- summary(effect_sizes)
effects_summary$study_number <- rownames(effects_summary)

ggplot(data = effects_summary, aes(y = yi, x = study_number)) + geom_point(size = 2) +
	geom_hline(yintercept = 0) +
	geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.1) + ggtitle("Growth responses") +
	ylab("log response ratio") + xlab("study number")


effect_sizes2 <- escalc(Mean_Elevated/SD_Elevated ~ factor(Food.supply) | factor(Paper_no),data = growth_clean, weights = N_Elevated, measure = "ROM",slab = Paper_No)
effects_summary2 <- summary(effect_sizes2)
effects_summary2$study_number <- rownames(effects_summary2)

ggplot(data = effects_summary2, aes(y = yi, x = study_number)) + geom_point(size = 2) +
	geom_hline(yintercept = 0) +
	geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.1) + ggtitle("Growth responses") +
	ylab("log response ratio") + xlab("study number")


## merge low and high CO2

effects_summary2$CO2_treatment <- "elevated CO2"
effects_summary$CO2_treatment <- "ambient CO2"

effects_all <- bind_rows(effects_summary, effects_summary2)


## plot all the of the food effects together



effects_all %>% 
	arrange(desc(yi)) %>% 
ggplot(data = ., aes(y = yi, x = study_number, group = CO2_treatment, color = CO2_treatment)) + geom_point(stat = "identity", size = 4) +
	geom_hline(yintercept = 0) +
	geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.1) + ggtitle("Growth responses, food effect at high and low CO2") +
	ylab("LnRR [growth at high food / low food]") + xlab("study number") + theme_bw()

## is the effect of food explained by CO2 treatment? No.
food_effect <- lm(yi ~ CO2_treatment, data = effects_all)

summary(food_effect)
