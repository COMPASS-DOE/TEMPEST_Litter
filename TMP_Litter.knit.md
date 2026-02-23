---
title: "Alice's TEMPEST litter analysis"
author: "AES"
date: "23 February, 2026"
output: html_document
---



## Read in Data


``` r
message("Reading litter data...")
```

```
## Reading litter data...
```

``` r
litter <- read_csv("TEMPEST litter data - Masses.csv", 
                   skip = 1, show_col_types = FALSE) 
```

```
## New names:
## • `` -> `...17`
```

```
## Warning: One or more parsing issues, call `problems()` on your data frame for details,
## e.g.:
##   dat <- vroom(...)
##   problems(dat)
```

``` r
## Convert date collected column into an actual date format,
## and isolate columns we care about
litter %>%
  filter(!is.na(Date_collected)) %>% 
  mutate(Date_collected = mdy(Date_collected),
         Year = year(Date_collected)) %>% 
  select(Plot, Trap, Date_collected, Year, 
         Leaf_ACRU, Leaf_FAGR, Leaf_LITU, Leaf_Other, All_Leaf) ->
  litter

message("\tdata read ok! rows = ", nrow(litter))
```

```
## 	data read ok! rows = 1980
```

``` r
message("Reading climate data...")
```

```
## Reading climate data...
```

``` r
clim <- rast("ERA5_data/8b6cbe648672d22c1a733be76b08f899.grib")
climdat_era <- terra::extract(clim, data.frame(lon = -76.551, lat = 38.875))
```

```
## Warning: [extract] transforming vector data to the CRS of the raster
```

``` r
climdat_era %>% 
  pivot_longer(-ID) %>% 
  mutate(time = terra::time(clim),
         Year = year(time)) %>% 
  group_by(Year, name) %>% 
  summarise(value = mean(value), .groups = "drop") %>% 
  pivot_wider(names_from = "name") ->
  climdat

names(climdat) <- c("Year", "Tair2m_C", "Srad_Ws_m2", "Precip_cm")
# Convert precipitation from average m to total cm
climdat$Precip_cm <- climdat$Precip_cm * 365 * 100
# Convert Tair from K to degC
climdat$Tair2m_C <- climdat$Tair2m_C - 273.1
```

## Calculate/fill in All-Leaf column


``` r
## if else means; condition is na is all leaf; if true add species and put in all leaf; if false use what is in all leaf
litter %>% 
  mutate(All_Leaf = if_else(condition = is.na(All_Leaf),
                            true = Leaf_ACRU + Leaf_FAGR + Leaf_LITU + Leaf_Other,
                            false = All_Leaf)) ->
  litter
```

## Turn these collection masses into daily fluxes

* Step 1: Calculate how long it's been since the last collection
* Step 2: Divide leaf litter biomass by the time computed in Step 1


``` r
# for each plot and trap we want to compute the days

litter %>%
  arrange(Date_collected) %>%
  group_by(Plot, Trap) %>% 
  # compute the difference between collection dates. For a series of
  # size n, diff() returns n-1 intervals, so we add one NA at the front
  # because the first date has no elapsed time (since no previous date)
  mutate(elapsed_days = c(NA, diff(Date_collected))) %>% 
  ungroup() %>%
  arrange(Plot, Trap, Date_collected) %>% 
  ## finally compute flux: litter per day
  mutate(All_Leaf_perday = All_Leaf / elapsed_days) ->
  litter
```

## Step 3: Interpolate rates between collection dates

### Control A 2019 example


``` r
## to interpolate between two dates, we need to create empty rows
## for example, if we collected on 10/1 and 10/10 right now we only have those rows in our data set
## so we need to create rows from 10/2 . . . 10/9 and then interpolate between the 10/1 flux and 10/10 flux to fill in those rows
## we will use complete() to generate the non-measured intervening dates
## and then we will use zoo::na.approx() to interpolate

# Date example
test <- tibble(day = mdy(c("10/20/2021", "10/31/2021")),
               values = c(60, 97))
test %>% 
  complete(day = seq(min(day), max(day), by = "day")) %>% 
  mutate(values = na.approx(values))
```

```
## # A tibble: 12 × 2
##    day        values
##    <date>      <dbl>
##  1 2021-10-20   60  
##  2 2021-10-21   63.4
##  3 2021-10-22   66.7
##  4 2021-10-23   70.1
##  5 2021-10-24   73.5
##  6 2021-10-25   76.8
##  7 2021-10-26   80.2
##  8 2021-10-27   83.5
##  9 2021-10-28   86.9
## 10 2021-10-29   90.3
## 11 2021-10-30   93.6
## 12 2021-10-31   97
```

``` r
# Examine: just the control 'A' trap in 2019
litter %>% 
  filter(Plot == "C", Trap == "A", year(Date_collected) == 2019) %>% 
  select(Plot, Trap, Date_collected, All_Leaf) ->
  litter_C_A_example

# Generate a sequence between first and last collection dates
days_wanted <- seq(from = min(litter_C_A_example$Date_collected), 
                   to = max(litter_C_A_example$Date_collected), 
                   by = "day")

# Just as an EXAMPLE Generate for the entire year NOT IDEAL as
# this hard-codes dates which we don't want to do
## days_wanted <- seq(from = ymd ("2019-01-01"), 
                  # to = ymd ("2019-12-31"), 
                  # by = "day")

# We want all days for all combinations of plot and trap

litter_C_A_example %>%
  complete(Plot, Trap, Date_collected = days_wanted) %>%
  mutate(All_Leaf = na.approx(All_Leaf, na.rm = FALSE)) ->
  litter_C_A_2019

# Visualize our example control A 2019 data

ggplot(litter_C_A_2019, aes(x = Date_collected, y = All_Leaf)) +
  geom_point()
```

<img src="TMP_Litter_files/figure-html/control-a-2019-example-1.png" width="672" />

### Control trap A over all years


``` r
## Control trap A over all years

litter %>% 
  filter(Plot == "C", Trap == "A") %>% 
  select(Plot, Trap, Date_collected, Year, All_Leaf) ->
  litter_C_A_allyears

# Compute the dates we want to complete()
days_wanted <- seq(from = min(litter_C_A_allyears$Date_collected, na.rm = TRUE), 
                   to = max(litter_C_A_allyears$Date_collected, na.rm = TRUE), 
                   by = "day")

# Estimate daily fluxes for all days of all years
litter_C_A_allyears %>%
  # complete the dataset, inserting missing dates
  complete(Plot, Trap, Date_collected = days_wanted) %>%
  # interpolate
  mutate(All_Leaf = na.approx(All_Leaf, na.rm = FALSE),
         # re-compute this, since we inserted new dates
         Year = year(Date_collected)) %>% 
  # compute cumulative flux for each year
  group_by(Year) %>% 
  mutate(All_Leaf_Cumulative = cumsum(All_Leaf)) ->
  litter_C_A_daily

ggplot(litter_C_A_daily, aes(x = Date_collected, y = All_Leaf)) + 
  geom_point(na.rm = TRUE) +
  ylab("Leaf Litter Flux (g/day)") +
  xlab("Date") +
  geom_vline(xintercept = ymd(paste0(2020:2025, "-01-01")), 
             color = "blue", linetype = 2)
```

<img src="TMP_Litter_files/figure-html/control-A-example-1.png" width="672" />

``` r
ggplot(litter_C_A_daily, aes(x = yday(Date_collected),
                             y = All_Leaf_Cumulative, 
                             color = factor(Year))) +
  geom_point(na.rm = TRUE) 
```

<img src="TMP_Litter_files/figure-html/control-A-example-2.png" width="672" />

**This shows we have a problem.** In some years the first collection isn't
until well into the year, and as a result the interpolated litterfall
continues throughout winter (i.e., the first part of the year, even
though presumably is almost all occurred late in the previous year).

First and last collection dates by year:


``` r
litter %>% 
  group_by(Year = year(Date_collected)) %>% 
  summarise(first = min(Date_collected), 
            last = max(Date_collected)) %>% 
  knitr::kable()
```



| Year|first      |last       |
|----:|:----------|:----------|
| 2018|2018-11-08 |2018-12-19 |
| 2019|2019-04-16 |2019-12-20 |
| 2020|2020-06-19 |2020-12-30 |
| 2021|2021-03-09 |2021-12-10 |
| 2022|2022-01-25 |2022-10-28 |
| 2023|2023-01-27 |2023-12-14 |
| 2024|2024-09-06 |2024-12-03 |
| 2025|2025-04-25 |2025-11-18 |
| 2026|2026-01-06 |2026-01-06 |

## Account for lack of collections right at the year changeover


``` r
PREV_YEAR_FRAC <- 0.95
```

**Fix for this issue:** for each year, we find the first collection and 
assume that 0.95 of that actually occurred by 12/31 the 
previous year.

So, we insert 'fake' collections on 12/31 (previous year), as well as a
fake zero collection on 1/1 of the current year. Finally, adjust the 
first collection to 0.05 of its value.

This is a PITA! But this seems the best way to deal with it reproducibly.


``` r
# Flag the first collections for each year, plot, trap
litter %>% 
  group_by(Year, Plot, Trap) %>% 
  # compare date collected with min date collected and assign it to "first"
  mutate(first = Date_collected == min(Date_collected)) %>% 
  ungroup() %>% 
  # we don't do this for 2018 data, the first year
  mutate(first = if_else(Year == 2018, FALSE, first)) ->
  litter

# Create artificial collections on 12/31 of the previous year
litter %>% 
  filter(first) %>% 
  mutate(Leaf_ACRU = PREV_YEAR_FRAC * Leaf_ACRU,
         Leaf_FAGR = PREV_YEAR_FRAC * Leaf_FAGR,
         Leaf_LITU = PREV_YEAR_FRAC * Leaf_LITU,
         Leaf_Other = PREV_YEAR_FRAC * Leaf_Other,
         All_Leaf = PREV_YEAR_FRAC * All_Leaf,
         Date_collected = ymd(paste0(Year - 1, "-12-31")),
         Year = year(Date_collected)) ->
  inferred_end_of_year

# Create artificial collections on 1/1 of the current year
litter %>% 
  filter(first) %>% 
  mutate(Leaf_ACRU = 0,
         Leaf_FAGR = 0,
         Leaf_LITU = 0,
         Leaf_Other = 0,
         All_Leaf = 0,
         Date_collected = ymd(paste0(Year, "-01-01"))) ->
  inferred_start_of_year

# Adjust the first collections to 5% of their value
litter %>% 
  filter(first) %>% 
  mutate(Leaf_ACRU = (1 - PREV_YEAR_FRAC) * Leaf_ACRU,
         Leaf_FAGR = (1 - PREV_YEAR_FRAC) * Leaf_FAGR,
         Leaf_LITU = (1 - PREV_YEAR_FRAC) * Leaf_LITU,
         Leaf_Other = (1 - PREV_YEAR_FRAC) * Leaf_Other,
         All_Leaf = (1 - PREV_YEAR_FRAC) * All_Leaf) ->
  new_firsts

# Combine these together and re-compute things
litter %>% 
  # get rid of original first numbers/data "!" b/c we've recalculated these
  filter(!first) %>% 
  # add rows; take litter dataset, drop away first data collection and add on 
    # new data calculations, 12/31, & 1/1 dates
  bind_rows(new_firsts, inferred_end_of_year, inferred_start_of_year) %>% 
  # re-do the elapsed_days calculation
  arrange(Date_collected) %>%
  group_by(Plot, Trap) %>% 
  mutate(elapsed_days = c(NA, diff(Date_collected))) %>% 
  ungroup() %>%
  mutate(All_Leaf_perday = All_Leaf / elapsed_days) ->
  litter_adjusted
```

Whew! Re-do the Control A example above:


``` r
litter_adjusted %>% 
  filter(Plot == "C", Trap == "A") %>% 
  select(Plot, Trap, Date_collected, Year, All_Leaf) ->
  litter_C_A_allyears

# Compute the dates we want to complete()
days_wanted <- seq(from = min(litter_C_A_allyears$Date_collected, na.rm = TRUE), 
                   to = max(litter_C_A_allyears$Date_collected, na.rm = TRUE), 
                   by = "day")

# Estimate daily fluxes for all days of all years
litter_C_A_allyears %>%
  # complete the dataset, inserting missing dates
  complete(Plot, Trap, Date_collected = days_wanted) %>%
  # interpolate
  mutate(All_Leaf = na.approx(All_Leaf, na.rm = FALSE),
         # re-compute this, since we inserted new dates
         Year = year(Date_collected)) %>% 
  # compute cumulative flux for each year
  group_by(Year) %>% 
  mutate(All_Leaf_Cumulative = cumsum(All_Leaf)) ->
  litter_C_A_daily

ggplot(litter_C_A_daily, aes(x = Date_collected, y = All_Leaf)) + 
  geom_point(na.rm = TRUE) +
  ylab("Leaf Litter Flux (g/day)") +
  xlab("Date") +
  geom_vline(xintercept = ymd(paste0(2020:2025, "-01-01")), 
             color = "blue", linetype = 2)
```

<img src="TMP_Litter_files/figure-html/unnamed-chunk-2-1.png" width="672" />

``` r
ggplot(litter_C_A_daily, aes(x = yday(Date_collected),
                             y = All_Leaf_Cumulative, 
                             color = factor(Year))) +
  geom_point(na.rm = TRUE) 
```

<img src="TMP_Litter_files/figure-html/unnamed-chunk-2-2.png" width="672" />

## Full dataset over time


``` r
## Next steps:
# 1. extend litter control A example to plot total litter biomass by year; 
# take litter adjust and compute total flux for each year

litter_adjusted %>%
  arrange(Date_collected) %>%
  group_by(Plot, Trap) %>% 
  # compute the difference between collection dates. For a series of
  # size n, diff() returns n-1 intervals, so we add one NA at the front
  # because the first date has no elapsed time (since no previous date)
  mutate(elapsed_days = c(NA, diff(Date_collected))) %>% 
  ungroup() %>%
  arrange(Plot, Trap, Date_collected) %>% 
  ## finally compute flux: litter per day
  mutate(All_Leaf_perday = All_Leaf / elapsed_days) ->
  litter_adj_perday

# generating the sequence of dates from min day collected to max
days_wanted_adj <- seq(from = min(litter_adj_perday$Date_collected), 
                   to = max(litter_adj_perday$Date_collected), 
                   by = "day")

litter_adj_perday %>%
  complete(Plot, Trap, Date_collected = days_wanted_adj) %>%
  mutate(All_Leaf_perday = na.approx(All_Leaf_perday, na.rm = FALSE)) ->
  litter_adj_perday_complete

# ALICE: I have simplified things above and changed it so we're not
# filtering for control A anymore; we want to process and visualize
# ALL the data.

# 2. plot

ggplot(litter_adj_perday_complete, aes(x = Date_collected, y = All_Leaf_perday, color = Plot)) +
  geom_point() +
  facet_grid(Plot ~ .)
```

```
## Warning: Removed 2629 rows containing missing values or values outside the scale range
## (`geom_point()`).
```

<img src="TMP_Litter_files/figure-html/unnamed-chunk-3-1.png" width="672" />

``` r
# Annual flux by plot and trap
litter_adj_perday_complete %>%
  mutate(Year = year(Date_collected)) %>%
# group by plot, trap, year
  group_by(Plot, Trap, Year) %>%
# sum up the daily fluxes to an annual values
  summarise(All_Leaf_peryear = sum(All_Leaf_perday), .groups = "drop") ->
  litter_annual_fluxes

# boxplot
ggplot(litter_annual_fluxes, aes(x = Year, y = All_Leaf_peryear, color = Plot, group = paste(Year, Plot))) +
  geom_boxplot(na.rm = TRUE)
```

<img src="TMP_Litter_files/figure-html/unnamed-chunk-3-2.png" width="672" />

``` r
# play around with different visualizations

# histogram

ggplot(litter_annual_fluxes, aes(x = Year, y = All_Leaf_peryear, color = Plot, group = paste(Year, Plot))) +
  scale_color_manual(values = c("green", "blue", "red")) +
  geom_boxplot(na.rm = TRUE)
```

<img src="TMP_Litter_files/figure-html/unnamed-chunk-3-3.png" width="672" />




``` r
##LITU's
litter_adjusted %>%
  arrange(Date_collected) %>%
  group_by(Plot, Trap) %>% 
  # compute the difference between collection dates. For a series of
  # size n, diff() returns n-1 intervals, so we add one NA at the front
  # because the first date has no elapsed time (since no previous date)
  mutate(elapsed_days = c(NA, diff(Date_collected))) %>% 
  ungroup() %>%
  arrange(Plot, Trap, Date_collected) %>% 
  ## finally compute flux: litter per day
  mutate(Leaf_LITU_perday = Leaf_LITU / elapsed_days) ->
  litter_adj_LITU_perday


# generating the sequence of dates from min day collected to max
days_wanted_adj <- seq(from = min(litter_adj_LITU_perday$Date_collected), 
                   to = max(litter_adj_LITU_perday$Date_collected), 
                   by = "day")

litter_adj_LITU_perday %>%
  complete(Plot, Trap, Date_collected = days_wanted_adj) %>%
  mutate(LITU_perday = na.approx(Leaf_LITU_perday, na.rm = FALSE)) ->
  litter_adj_LITU_perday_complete

ggplot(litter_adj_LITU_perday_complete, aes(x = Date_collected, y = Leaf_LITU_perday, color = Plot)) +
  geom_point(na.rm = TRUE) +
  scale_color_manual(values = c("green", "blue", "red")) +
  facet_grid(Plot ~ .)
```

<img src="TMP_Litter_files/figure-html/unnamed-chunk-4-1.png" width="672" />

