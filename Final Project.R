# Set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Import libraries
require(tidyverse)

data = read_csv('data/REC - provided by Phil.csv')
census = read_csv('data/census_postcode.csv')

# Data Cleaning
agg.data = data %>%
  group_by(POSTCODE) %>%
  summarise(count=sum(COUNT), KW=sum(KW)) 

census$POSTCODE = parse_number(census$POA_CODE_2016)

data.merge = merge(agg.data, census, by = 'POSTCODE' )

# Relationship of Count and KW
ggplot(data.merge, aes(KW,count)) + geom_point()

# Building linear model
model = lm(count ~ Total_Persons_Persons + Median_age_persons + Median_mortgage_repay_monthly +
             Median_tot_hhd_inc_weekly + Average_household_size, data.merge)
summary(model)

