# Data wrangling


# load packages -----------------------------------------------------------

library(tidyverse)

# read in data ------------------------------------------------------------

calci <- read_csv("data-raw/calcification_CO2.csv")

calci$CO2_elevated[calci$CO2_elevated == "1094"] <- "1065" ## assuming this discrepancy in CO2 values is a typo, so I'm fixing it

calci_top <- calci %>%
	select(1:9, 11:30, 10) %>%
	group_by(Author, Food.supply) %>% 
	top_n(1)

calci_top_metadata <- calci_top %>% 
											select(1:15, 30)

data_long <- calci_top %>% 
	gather(response_variable, value, 16:25) %>% 
	select(response_variable, value, everything()) %>%
	spread(Food.supply, value) %>% 
	select(High, Low, everything()) %>%
	separate(response_variable, c("statistic", "CO2_level"), sep = "_") %>% 
	rename(Highfood = High,
				 Lowfood = Low)

Highfood <- data_long %>% 
	select(Highfood, Author, CO2_level, statistic) %>% 
	filter(!is.na(Highfood)) %>% 
	spread(statistic, Highfood) %>%
	rename(high_food_95CI = `95CI`,
				 high_food_mean = Mean,
				 high_food_N = N,
				 high_food_SE = SE,
				 high_food_SD = SD)

Lowfood <- data_long %>% 
	select(Lowfood, Author, CO2_level, statistic) %>% 
	filter(!is.na(Lowfood)) %>% 
	spread(statistic, Lowfood) %>%
	rename(low_food_95CI = `95CI`,
				 low_food_mean = Mean,
				 low_food_N = N,
				 low_food_SE = SE,
				 low_food_SD = SD)

all <- left_join(Lowfood, Highfood)
all_2 <- left_join(all, calci_top_metadata) %>%
	distinct(low_food_mean, .keep_all = TRUE) %>%
	select(-Food.supply) %>% 
	select(1:2, 4:7, 9:12, 3, 8, everything())

write_csv(all_2, "calcification_by_foodsupply.csv")

