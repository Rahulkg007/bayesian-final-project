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
data.merge= data.merge %>% select(-c(KW,POA_CODE_2016,POSTCODE))


# Building linear model
model = lm(count ~ ., data.merge)
summary(model)

# Hi Kob
# Time to commit

