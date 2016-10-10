# Going from wide to long and back!
# Joey Bernhardt 
# Oct 1 2016 
# Main goal: Rework the calcification 
# data frame so that food supply is in columns to make it easier to analyze with
# the metafor package

# load packages -----------------------------------------------------------

library(tidyverse)
library(janitor)

# read in data ------------------------------------------------------------

calci <- read_csv("data-raw/calcification_multiple_CO2.csv") %>% 
	remove_empty_cols() %>% 
	clean_names() %>% 
	rename(extra_notes = x32)
growth <- read_csv("data-raw/growth_multiple_CO2.csv") %>% 
	remove_empty_cols() %>% 
	clean_names() %>%
	rename(extra_notes = x33)




# select out unnecessary columns---------------------------------------------------


growth_unique <- growth %>% 
	select(5:27) %>% 
	unite(uniqueID, author, co_treatment_coarse, unit, remove = FALSE)


# create a metadata table -------------------------------------------------

growth_metadata <- growth_unique %>% 
select(1:14)

# reshape the data.frame to plop food supply into cols --------------------


data_long <- calci_top %>% 
	gather(response_variable, value, 16:25) %>% 
	select(response_variable, value, everything()) %>%
	spread(Food.supply, value) %>% 
	select(High, Low, everything()) %>%
	separate(response_variable, c("statistic", "CO2_level"), sep = "_") %>% 
	rename(Highfood = High,
				 Lowfood = Low)

# now try with growth and mulitple CO2 levels -----------------------------


growth_long <- growth_unique %>% 
	gather(response_variable, value, 15:24) %>% 
	select(response_variable, value, everything()) %>%
	spread(food_supply, value) %>% 
	select(High, Low, everything()) %>%
	separate(response_variable, c("statistic", "CO2_level"), sep = "_") %>% 
	rename(Highfood = High,
				 Lowfood = Low) %>% 
	group_by(statistic, CO2_level, uniqueID) %>% 
	summarise(low_summary = mean(Lowfood, na.rm = TRUE),
						high_summary = mean(Highfood, na.rm = TRUE))





# get high food cols into wide format -------------------------------------


Highfood <- growth_long %>% 
	select(high_summary, uniqueID, CO2_level, statistic) %>% 
	# filter(!is.na(Highfood)) %>%  # get of NA's associated with the gather
	spread(statistic, high_summary) %>%
	rename(high_food_95CI = x95ci,
				 high_food_mean = mean,
				 high_food_N = n,
				 high_food_SE = se,
				 high_food_SD = sd) %>%
	separate(uniqueID, c("Author", "CO2_magnitude", "Unit"), sep = "_", remove = FALSE)


# get low food cols into wide format --------------------------------------

Lowfood <- growth_long %>% 
	select(low_summary, uniqueID, CO2_level, statistic) %>% 
	# filter(!is.na(Lowfood)) %>% # get of NA's associated with the gather
	spread(statistic, low_summary) %>%
	rename(low_food_95CI = x95ci,
				 low_food_mean = mean,
				 low_food_N = n,
				 low_food_SE = se,
				 low_food_SD = sd) %>%
	separate(uniqueID, c("Author", "CO2_magnitude", "Unit"), sep = "_", remove = FALSE)

glimpse(Lowfood)
glimpse(Highfood)

# join them all together! --------------------------------------------------


all <- left_join(Lowfood, Highfood)

all_2 <- left_join(all, growth_metadata) %>% 
	distinct(uniqueID, .keep_all = TRUE) %>% ## get rid of duplicated rows after join with metadata
	select(-food_supply) %>%
	select(1:2, 4:7, 9:12, 3, 8, everything()) %>% 
	select(-unit)


# write out as csv :) -----------------------------------------------------


write_csv(all_2, "data-processed/growth_by_foodsupply_allCO2_levels.csv")
 

# now onto calcification --------------------------------------------------

# select out unnecessary columns---------------------------------------------------


calci_unique <- calci %>% 
	select(5:27) %>% 
	unite(uniqueID, author, co_treatment_coarse, unit, remove = FALSE)


# create a metadata table -------------------------------------------------

calci_metadata <- calci_unique %>% 
	select(1:14)

# reshape the data.frame to plop food supply into cols --------------------

calci_long <- calci_unique %>% 
	gather(response_variable, value, 15:24) %>% 
	select(response_variable, value, everything()) %>%
	spread(food_supply, value) %>% 
	select(High, Low, everything()) %>%
	separate(response_variable, c("statistic", "CO2_level"), sep = "_") %>% 
	rename(Highfood = High,
				 Lowfood = Low) %>% 
	group_by(statistic, CO2_level, uniqueID) %>% 
	summarise(low_summary = mean(Lowfood, na.rm = TRUE),
						high_summary = mean(Highfood, na.rm = TRUE))





# get high food cols into wide format -------------------------------------


Highfood_calc <- calci_long %>% 
	select(high_summary, uniqueID, CO2_level, statistic) %>% 
	# filter(!is.na(Highfood)) %>%  # get of NA's associated with the gather
	spread(statistic, high_summary) %>%
	rename(high_food_95CI = x95ci,
				 high_food_mean = mean,
				 high_food_N = n,
				 high_food_SE = se,
				 high_food_SD = sd) %>%
	separate(uniqueID, c("Author", "CO2_magnitude", "Unit"), sep = "_", remove = FALSE)


# get low food cols into wide format --------------------------------------

Lowfood_calc <- calci_long %>% 
	select(low_summary, uniqueID, CO2_level, statistic) %>% 
	# filter(!is.na(Lowfood)) %>% # get of NA's associated with the gather
	spread(statistic, low_summary) %>%
	rename(low_food_95CI = x95ci,
				 low_food_mean = mean,
				 low_food_N = n,
				 low_food_SE = se,
				 low_food_SD = sd) %>%
	separate(uniqueID, c("Author", "CO2_magnitude", "Unit"), sep = "_", remove = FALSE)



# join them all together! --------------------------------------------------


all_calc <- left_join(Lowfood_calc, Highfood_calc)

all_2_calc <- left_join(all_calc, calci_metadata) %>% 
	distinct(uniqueID, .keep_all = TRUE) %>% ## get rid of duplicated rows after join with metadata
	select(-food_supply) %>%
	select(1:2, 4:7, 9:12, 3, 8, everything()) %>% 
	select(-unit)


# write out as csv :) -----------------------------------------------------


write_csv(all_2_calc, "data-processed/calci_by_foodsupply_allCO2_levels.csv")

