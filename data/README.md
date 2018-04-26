# MS_stops.csv

**Mississippi Traffic Stops by Police 2013-2016**

Data from https://openpolicing.stanford.edu

For the purpose of this tutorial I made the following changes:

- Selected Black, White, and unknown (empty string: "") "driver_race" only.
- Kept only the following columns "id", "state", "stop_date", "county_name", "county_fips", "police_department", "driver_gender", "driver_age_raw", "driver_race",  "violation_raw", "officer_id"
- Renamed "driver_age_raw" to "driver_birthdate"
- Set all empty strings "" to NA for:
    - "driver_gender"
    - "driver_race"
    - "driver_birthdate"
    - "officer_id"
- Set "driver_birthdate" to NA if it is later or equal than "stop_date"
- Calculated age of driver at the time stopped and added a column with driver_age (in years)
- "driver_gender": Renamed M and F to male and female
- added violations, recoded, based on violations_raw:
[1] Seat belt not used properly as required >> Seat belt
[2] Careless driving >> Carless driving
[3] Speeding - Regulated or posted speed limit and actual speed >> Speeding
[4] ?? >> Other or unknown
[5] Failure to obey sign or traffic control device >> Carless driving
[6] Driving while license suspended  >> License/Permit/Insurance
[7] Failure to maintain required liability insurance >> License/Permit/Insurance
[8] Other (non-mapped) >> Other or unknown
[9] Expired or no non-commercial driver license or permit >> License/Permit/Insurance
[10] Child or youth restraint not used properly as required >> Seat belt
[11] Improper turn >> Carless driving
[12] Driving wrong way >> Carless driving
[13] Operating without equipment as required by law >>  Breaks/Lights/etc
[14] Speeding  >> Speeding
[15] Following too closely >> Carless driving
[16] Reckless driving >> Carless driving
[17] Improper passing >> Carless driving
[18] Failure to yield right of way (FTY ROW) >> Carless driving
[19] Failure to comply with financial responsibility law >> License/Permit/Insurance

```r
traf <- read.csv("../R-intro/data/MS_trafficstops_bw.csv", na.strings=c("","NA"), stringsAsFactors = F)
   traf$driver_birthdate[ymd(traf$driver_birthdate) >= ymd(traf$stop_date)] <- NA
   traf$driver_age <- round((ymd(traf$stop_date) - ymd(traf$driver_birthdate))/365)
   traf$driver_gender[traf$driver_gender == "F"] <- "female"
   traf$driver_gender[traf$driver_gender == "M"] <- "male"
   
write.csv(traf, "data/MS_stops.csv", row.names=F)
## recoding factor levels manually...
```

## UPDATE: Added two more fields 
```r
traf %>% 
mutate(wk_day = wday(stop_date, label = TRUE, y_day = yday(stop_date))) %>%
write.csv("data/MS_stops.csv", row.names = F)
```
# MS_county_stops.csv

**Mississippi Traffic Stops by Police 2013-2016, aggregated per county with ACS population counts and ratios**

Using estimated values of the 5 year average of the 2011-2015 American Community Survey (ACS) for: B02001. Race: Black Black or African American alone and White alone.

```r
library(tidycensus)
#census_api_key("XXXXXX")

black_pop <- get_acs(geography = "county", variable = "B02001_003E", state="MS")
white_pop <- get_acs(geography = "county", variable = "B02001_002E", state="MS")

MS_pop <- data.frame(FIPS = as.numeric(black_pop$GEOID), black_pop = black_pop$estimate, white_pop = white_pop$estimate, total_pop = total_pop$estimate)

#write.csv(MS_bw, "data/MS_acs2015_bw.csv", row.names = F)`
```

Then:

```r
traf %>% 
     filter(!is.na(driver_race)) %>% 
     count(county_name, county_fips, driver_race) %>% 
     spread(driver_race, n, fill = 0, sep = "_") %>%  
     left_join(MS_bw, by = c("county_fips" = "FIPS")) %>% 
     mutate(pct_black_stopped = driver_race_Black/black_pop,
         pct_white_stopped = driver_race_White/white_pop) %>% 
     write.csv(file = "data_output/MS_county_stops.csv", row.names = FALSE)`
```

## UPDATE: Added two more fields 

```r
wb_delta <- pct_white_stopped - pct_black_stopped  
bias <- ifelse(wb_delta < 0, "black bias", "white bias")
```