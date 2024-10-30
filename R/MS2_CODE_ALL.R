# Long-term increases in species richness and fish community reorganization in a mid-atlantic estuary 
## Haley Oleynik
## August 2024 

# Load libraries -------------
# use ggplot theme_light()
require(tidyverse)
require(mice)
require(Kendall)
require(trend)
require(stats)
require(boot)
require(vegan)
require(reshape2)
require(rMR) 
require(patchwork)
require(zoo)
require(data.table)
require(ggplot2)
require(readr)
require(ggridges)

# Colorblind palette for all plots ----------- 
colorblind_palette <- c(
  "#E69F00",  # Orange
  "#56B4E9",  # Sky Blue
  "#009E73",  # Bluish Green
  "#F0E442",  # Yellow
  "#0072B2",  # Blue
  "#D55E00",  # Vermilion
  "#CC79A7"   # Reddish Purple
)
# Read Data ---------
# small trawl 
smalltrawl <- read.csv("/Users/haleyoleynik/Documents/Thesis/Data/SmallTrawl_DATA.csv")
groups <- read.csv("/Users/haleyoleynik/Documents/Thesis/Data/SpeciesGroups2.csv")

smalltrawl <- smalltrawl %>%
  left_join(groups, by = c("Common.name" = "CommonName"))

New.Data <- read.csv("/Users/haleyoleynik/Documents/Thesis/Data/2019TrawlData.csv")
New.Data$BottomDO_sat <- DO.unit.convert(New.Data$BottomDO, DO.units.in = "mg/L", DO.units.out = "pct", bar.units.in = "atm", bar.press = 1.013253, temp.C = New.Data$BottomTemp, salinity = New.Data$BottomSal, salinity.units = "pp.thou")

New.Data$SurfaceDO_sat <- DO.unit.convert(New.Data$SurfaceDO, DO.units.in = "mg/L", DO.units.out = "pct", bar.units.in = "atm", bar.press = 1.013253, temp.C = New.Data$SurfaceTemp, salinity = New.Data$SurfaceSal, salinity.units = "pp.thou")

bayregions <- read_csv("/Users/haleyoleynik/Documents/Thesis/Data/Map Stations/Bay_Regions.csv")

# new fish info with everything 
fishinfo.30 <- readr::read_csv("/Users/haleyoleynik/Documents/Thesis/Data/Fish Base/30-foot fish info v2.csv")
fishinfo.16 <- readr::read_csv("/Users/haleyoleynik/Documents/Thesis/Data/Fish Base/16-foot fish info v2.csv")

# Trends -----------
## 30-foot survey -------------
# from MS2_Trends.R 

### Plot env. variables through time ------------
temp.TS <- New.Data %>%
  group_by(Month, Year) %>%
  summarise(temp=mean(BottomTemp,na.rm=TRUE),
            sal = mean(BottomSal,na.rm=TRUE),
            DO = mean(BottomDO, na.rm=TRUE),
            DO_sat = mean(BottomDO_sat, na.rm=TRUE)) %>%
  filter(Year > 1989) %>%
  mutate(Date = as.Date(paste(Year, Month, "1", sep="-"))) %>%
  ungroup() %>%
  select(-c(Month,Year)) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month")) 

# do not use !! 
p1 <- ggplot(temp.TS, aes(x=Date, y = temp)) +
  geom_line(color = "black", alpha = 0.7) +
  geom_smooth(method = "gam", se = T, size = 1, col = "red") + 
  geom_hline(yintercept = mean(temp.TS$temp, na.rm=T), linetype = 'dashed') +
  xlab("") +
  ylab("Temperature (C)") + 
  theme_light() +
  theme(text = element_text(size=15),axis.text.x = element_blank())

p2 <- ggplot(temp.TS, aes(x=Date, y = sal)) +
  geom_line(color = "black", alpha = 0.7) +
  geom_smooth(method = "gam", se = T, size = 1, col = "red") + 
  geom_hline(yintercept = mean(temp.TS$sal, na.rm=T), linetype = 'dashed') +
  xlab("") +
  ylab("Salinity (ppt)") + 
  theme_light() +
  theme(text = element_text(size=15), axis.text.x = element_blank())

p3 <- ggplot(temp.TS, aes(x=Date, y = DO_sat)) +
  geom_line(color = "black", alpha = 0.7) +
  geom_smooth(method = "gam", se = T, size = 1, col = "red") + 
  geom_hline(yintercept = mean(temp.TS$DO_sat, na.rm=T), linetype = 'dashed') +
  xlab("Year") +
  ylab("Dissolved Oxygen (%)") + 
  theme_light() +
  theme(text = element_text(size=15))

p1 / p2 / p3

### Temperature by month -------
temp <- New.Data %>%
  filter(Month > 2) %>%
  group_by(Year, Month) %>%
  summarise(temp = mean(BottomTemp, na.rm=T))

temp <- filter(temp, Month < 12)

#trim years
temp <- temp[temp$Year %in% c(1990:2019), ]

temp$Month <- as.factor(temp$Month)

ggplot(temp, aes(x=Year, y=temp, color = Month)) +
  geom_line() + 
  #geom_smooth(method = "loess", size = 1) +
  ylab("Temperature (C)") + 
  theme_light() +
  theme(text = element_text(size=30))

### Model temperature through time (to determine % increase) --------
model <- glm(temp ~ Year + Month, data = temp)
summary(model)

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["Year"] # 0.04357815 -- means increased 0.04 per year on average 

# Calculate the total number of years
n_years <- max(temp$Year) - min(temp$Year)

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year 
start_year <- min(temp$Year)
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall # 23.05986% increase overall 

### Model DO through time (to determine % increase) --------
do <- New.Data %>%
  filter(Year %in% c(1990:2019)) %>%
  filter(Month %in% c(2:11)) %>%
  group_by(Year, Month) %>%
  summarise(DO = mean(BottomDO_sat, na.rm=T))

model <- glm(DO ~ Year + Month, data = do)
summary(model)

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["Year"] # -0.1756533 -- means decreased 0.18 per year on average 
coef_year

# Calculate the total number of years
n_years <- max(do$Year) - min(do$Year)

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year 
start_year <- min(do$Year)
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall # 4.95337% decrease overall 


### ALL Richness through time ------------------
rich.TS <- New.Data %>%
  group_by(Month, Year, SpeciesCode) %>%
  summarise(number=sum(CPUE,na.rm=TRUE)) %>%
  filter(Year > 1989) %>%
  mutate(Date = as.Date(paste(Year, Month, "1", sep="-"))) %>%
  ungroup() %>%
  select(-c(Month,Year)) %>%
  group_by(Date) %>% 
  tally() %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month")) %>%
  rename(richness = n) 

rich <- ggplot(rich.TS, aes(x=Date, y = richness)) +
  geom_line(color = "#E69F00", alpha = 0.7) +
  geom_smooth(method = "lm", se = T, size = 1, col = "#E69F00") + 
  geom_hline(yintercept = mean(rich.TS$richness, na.rm=T), linetype = 'dashed') +
  labs(x = "", y="Species Richness", tag = "a)") +
  theme_light()

### Model richness through time (to determine % increase) ----- 
try <- rich.TS %>%
  mutate(month = month(Date)) %>%
  mutate(year = year(Date))

model <- glm(richness ~ year + month, data = try)
summary(model)

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["year"] # 0.3439364 -- means increased .34 per year on average 

# Calculate the total number of years
n_years <- max(try$year) - min(try$year)

# Get the coefficient of year
coef_year <- coef(model)["year"]

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year (e.g., 1990)
start_year <- 1990
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall # 56.73357% increase overall 

### FISH Richness through time ------------------
fish.rich.TS <- New.Data %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Month, Year, SpeciesCode) %>%
  summarise(number=sum(CPUE,na.rm=TRUE)) %>%
  mutate(Date = as.Date(paste(Year, Month, "1", sep="-"))) %>%
  ungroup() %>%
  select(-c(Month,Year)) %>%
  group_by(Date) %>% 
  tally() %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month")) %>%
  rename(richness = n) 

fish.rich <- ggplot(fish.rich.TS, aes(x=Date, y = richness)) +
  geom_line(color = "#E69F00", alpha = 0.7) +
  geom_smooth(method = "lm", se = T, size = 1, col = "#E69F00") + 
  geom_hline(yintercept = mean(fish.rich.TS$richness, na.rm=T), linetype = 'dashed') +
  labs(x = "", y="Species Richness", tag = "a)") +
  theme_light()

### Model richness through time (to determine % increase) ----- 
try <- fish.rich.TS %>%
  mutate(month = month(Date)) %>%
  mutate(year = year(Date))

model <- glm(richness ~ year + month, data = try)
summary(model)

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["year"] # 0.3439364 -- means increased .34 per year on average 

# Calculate the total number of years
n_years <- max(try$year) - min(try$year)

# Get the coefficient of year
coef_year <- coef(model)["year"]

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year (e.g., 1990)
start_year <- 1966
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall # 50.96816% increase overall 


### ALL Diversity through time -----------------------

div.TS <- New.Data %>%
  group_by(Month, Year, SpeciesCode) %>%
  summarise(Number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  filter(Year %in% 1990:2019) %>%
  mutate(Month.Year = paste(Year, Month, "1", sep="-"),
         Date = as.Date(Month.Year)) %>%
  select(-c(Month.Year,Year,Month)) %>% 
  pivot_wider(names_from = SpeciesCode, values_from = Number) %>%
  tibble::column_to_rownames('Date') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
shann <- diversity(div.TS, index="shannon")
shann <- as.data.frame(shann)
shann$Date = row.names(shann)

#join to make a diversity table, fill in missing date values with NA 
div.TS <- shann %>%
  #left_join(simp, by = "Date") %>%
  mutate(Date = as.Date(Date)) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month"))

# plot 
div <- ggplot(div.TS, aes(x=Date, y = shann)) +
  geom_line(color = "#E69F00", alpha = 0.7) +
  geom_smooth(method = "lm", se = T, size = 1, col = "#E69F00") + 
  geom_hline(yintercept = mean(div.TS$shann, na.rm=T), linetype = 'dashed') +
  labs(x="Year",y="Species Diversity") +
  theme_light() 

rich / div

### model diversity through time ---------
try <- div.TS %>%
  mutate(month = month(Date)) %>%
  mutate(year = year(Date))

model <- glm(shann ~ year + month, data = try)
summary(model)

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["year"] # 0.3439364 -- means increased .34 per year on average 

# Calculate the total number of years
n_years <- max(try$year) - min(try$year)

# Get the coefficient of year
coef_year <- coef(model)["year"]

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year (e.g., 1990)
start_year <- 1990
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall # 11.73% increase overall 

### FISH Diversity through time -----------------------

fish.div.TS <- New.Data %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Month, Year, SpeciesCode) %>%
  summarise(Number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(Month.Year = paste(Year, Month, "1", sep="-"),
         Date = as.Date(Month.Year)) %>%
  select(-c(Month.Year,Year,Month)) %>% 
  pivot_wider(names_from = SpeciesCode, values_from = Number) %>%
  tibble::column_to_rownames('Date') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
shann <- diversity(fish.div.TS, index="shannon")
shann <- as.data.frame(shann)
shann$Date = row.names(shann)

#join to make a diversity table, fill in missing date values with NA 
fish.div.TS <- shann %>%
  #left_join(simp, by = "Date") %>%
  mutate(Date = as.Date(Date)) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month"))

# plot 
fish.div <- ggplot(fish.div.TS, aes(x=Date, y = shann)) +
  geom_line(color = "#E69F00", alpha = 0.7) +
  geom_smooth(method = "lm", se = T, size = 1, col = "#E69F00") + 
  geom_hline(yintercept = mean(fish.div.TS$shann, na.rm=T), linetype = 'dashed') +
  labs(x="Year",y="Species Diversity") +
  theme_light() 

fish.rich / fish.div

### model diversity through time ---------
try <- fish.div.TS %>%
  mutate(month = month(Date)) %>%
  mutate(year = year(Date))

model <- glm(shann ~ year + month, data = try)
summary(model)

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["year"] # 0.3439364 -- means increased .34 per year on average 

# Calculate the total number of years
n_years <- max(try$year) - min(try$year)

# Get the coefficient of year
coef_year <- coef(model)["year"]

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year (e.g., 1990)
start_year <- 1966
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall # 8.195311% increase overall

### Rate of introduction / departure ------------------------------
# filter by fish because too hard to tell when inverts started being monitored 
heatmap <- New.Data %>%
  group_by(Year, CommonName) %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  summarise(number=sum(CPUE,na.rm=TRUE))

#take only 1990-present 
#heatmap <- heatmap[heatmap$Year %in% c(1990:2019), ]

firstyear <- heatmap %>%
  group_by(CommonName) %>%
  summarise(Year = min(Year))

#plot 
ggplot(firstyear, aes(x=Year, y=reorder(CommonName, Year))) +
  geom_point() +
  theme(text = element_text(size=10)) +
  ylab("Species")

#filter out 1990, then plot again 
firstyear <- filter(firstyear, Year > 1966)

# LAST YEAR 
lastyear <- heatmap %>%
  group_by(CommonName) %>%
  summarise(Year = max(Year))

#filter out 1990, then plot again 
lastyear <- filter(lastyear, Year < 2018)

#plot 
ggplot(lastyear, aes(x=Year, y=reorder(CommonName, Year))) +
  geom_point() +
  theme(text = element_text(size=15)) +
  ylab("Species")

# with just fish 
introduction.rate <- firstyear %>%
  group_by(Year) %>%
  tally()

departure.rate <- lastyear %>%
  group_by(Year) %>%
  tally()

# adjust this based on how many years to take out at beginning and end 
introduction.rate <- filter(introduction.rate, Year > 1970)
departure.rate <- filter(departure.rate, Year < 2015)

mean(introduction.rate$n) #1.827586 (-2 yrs) , 2.1
sd(introduction.rate$n)  #1.28366

mean(departure.rate$n) #1.583333 (-2 yrs) , 1.8
sd(departure.rate$n) #0.8297022

## 16-foot survey--------------------------------------

### ALL Diversity through time --------------------------------
small.div.TS <- smalltrawl %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Month, Year4, SpecCode) %>%
  summarise(Number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  filter(Year4 %in% 1990:2019) %>%
  filter(!is.na(SpecCode)) %>%
  mutate(Month.Year = paste(Year4, Month, "1", sep="-"),
         Date = as.Date(Month.Year)) %>%
  select(-c(Month.Year,Year4,Month)) %>% 
  pivot_wider(names_from = SpecCode, values_from = Number) %>%
  tibble::column_to_rownames('Date') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
shann <- diversity(small.div.TS, index="shannon")
shann <- as.data.frame(shann)
shann$Date = row.names(shann)

#join to make a diversity table
small.div.TS <- shann %>%
  #left_join(simp, by = "Date") %>%
  mutate(Date = as.Date(Date)) %>%
  rename(small.shann = shann) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month"))

### FISH Diversity through time --------------------------------
fish.small.div.TS <- smalltrawl %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Month, Year4, SpecCode) %>%
  summarise(Number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  #filter(Year4 %in% 1990:2019) %>%
  filter(!is.na(SpecCode)) %>%
  mutate(Month.Year = paste(Year4, Month, "1", sep="-"),
         Date = as.Date(Month.Year)) %>%
  select(-c(Month.Year,Year4,Month)) %>% 
  pivot_wider(names_from = SpecCode, values_from = Number) %>%
  tibble::column_to_rownames('Date') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
shann <- diversity(fish.small.div.TS, index="shannon")
shann <- as.data.frame(shann)
shann$Date = row.names(shann)

#join to make a diversity table
fish.small.div.TS <- shann %>%
  #left_join(simp, by = "Date") %>%
  mutate(Date = as.Date(Date)) %>%
  rename(small.shann = shann) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month"))

### ALL Richness through time----------------------------------------------
small.rich.TS <- smalltrawl %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Month, Year4, SpecCode) %>%
  summarise(number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  filter(Year4 %in% 1990:2019) %>%
  filter(!is.na(SpecCode)) %>%
  mutate(Date = as.Date(paste(Year4, Month, "1", sep="-"))) %>%
  select(-c(Month,Year4)) %>%
  group_by(Date) %>% 
  tally() %>%
  rename(small.richness = n) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month")) # fill in missing dates

### FISH Richness through time----------------------------------------------
fish.small.rich.TS <- smalltrawl %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Month, Year4, SpecCode) %>%
  summarise(number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  #filter(Year4 %in% 1990:2019) %>%
  filter(!is.na(SpecCode)) %>%
  mutate(Date = as.Date(paste(Year4, Month, "1", sep="-"))) %>%
  select(-c(Month,Year4)) %>%
  group_by(Date) %>% 
  tally() %>%
  rename(small.richness = n) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month")) # fill in missing dates

### Plot ---------- 
small.rich <- ggplot(small.rich.TS, aes(x=Date, y = small.richness)) +
  geom_line(color = "#56B4E9", alpha = 0.7) +
  geom_smooth(method = "loess", se = T, size = 1, col = "#56B4E9") + 
  geom_hline(yintercept = mean(small.rich.TS$small.richness, na.rm=T), linetype = 'dashed') +
  labs(x="",y="Species Richness", tag = "b)") +
  theme_light()

small.div <- ggplot(small.div.TS, aes(x=Date, y = small.shann)) +
  geom_line(color = "#56B4E9", alpha = 0.7) +
  geom_smooth(method = "loess", se = T, size = 1, col = "#56B4E9") + 
  geom_hline(yintercept = mean(small.div.TS$small.shann, na.rm=T), linetype = 'dashed') +
  labs(x="Year",y="Species Diversity") +
  theme_light() 

small.rich.fish <- ggplot(fish.small.rich.TS, aes(x=Date, y = small.richness)) +
  geom_line(color = "#56B4E9", alpha = 0.7) +
  geom_smooth(method = "loess", se = T, size = 1, col = "#56B4E9") + 
  geom_hline(yintercept = mean(fish.small.rich.TS$small.richness, na.rm=T), linetype = 'dashed') +
  labs(x="",y="Species Richness", tag = "b)") +
  theme_light()

small.div.fish <- ggplot(fish.small.div.TS, aes(x=Date, y = small.shann)) +
  geom_line(color = "#56B4E9", alpha = 0.7) +
  geom_smooth(method = "loess", se = T, size = 1, col = "#56B4E9") + 
  geom_hline(yintercept = mean(fish.small.div.TS$small.shann, na.rm=T), linetype = 'dashed') +
  labs(x="Year",y="Species Diversity") +
  theme_light() 

small.rich / small.div
small.rich.fish / small.div.fish

# FIGURE 1 - richness and diversity trends -----------------------
(fish.rich / fish.div) | (small.rich.fish / small.div.fish)

ggsave("Figures/FIGURE_1.png", 
       dpi=300, height=7, width=10, units='in')

ggsave("Figures/FIGURE_1b.png", 
       dpi=300, height=7, width=10, units='in')

### Mann Kendall tests -------
# 30 foot 
#make a time series object
rich.TSobj = ts(rich.TS$richness, frequency=12, start=c(1990,4))

#plot time series
plot(rich.TSobj)

#The null hypothesis for this test is that there is no monotonic trend in the series
MK = MannKendall(rich.TSobj)
summary(MK)  #tau = 0.144, 2-sided pvalue =0.00038302

#seasonal man-kendall- use this for my data 
SMK = SeasonalMannKendall(rich.TSobj)
summary(SMK) #tau = 0.186, 2-sided pvalue =1.484e-05

#simpson 
simp <- diversity(div.TS, index="simpson") 
simp <- as.data.frame(simp)
simp$Date = row.names(simp)

#make a time series object
div.TSobj.shann = ts(div.TS$shann, frequency=12, start=c(1990,4))
div.TSobj.simp = ts(div.TS$simp, frequency=12, start=c(1990,4))

#plot time series
plot(div.TSobj.shann)
plot(div.TSobj.simp)

#seasonal man-kendall- use this for my data 
SMK = SeasonalMannKendall(div.TSobj.shann)
summary(SMK) 

SMK = SeasonalMannKendall(div.TSobj.simp)
summary(SMK)

# 16- foot 
#make a time series object
div.TSobj.shann = ts(div.TS$shann, frequency=12, start=c(1990,4))
#div.TSobj.simp = ts(div.TS$simp, frequency=12, start=c(1990,4))

#plot time series
plot(div.TSobj.shann)
#plot(div.TSobj.simp)

# #simpson 
# simp <- diversity(div.TS, index="simpson") 
# simp <- as.data.frame(simp)
# simp$Date = row.names(simp)


#seasonal man-kendall- use this for my data 
SMK = SeasonalMannKendall(div.TSobj.shann)
summary(SMK) 

SMK = SeasonalMannKendall(div.TSobj.simp)
summary(SMK)

#make a time series object
rich.TSobj = ts(small.rich.TS$richness, frequency=12, start=c(1990,4))

#plot time series
plot(rich.TSobj)

#The null hypothesis for this test is that there is no monotonic trend in the series
MK = MannKendall(rich.TSobj)
summary(MK)  #tau = 0.144, 2-sided pvalue =0.00038302

#seasonal man-kendall- use this for my data 
SMK = SeasonalMannKendall(rich.TSobj)
summary(SMK) #tau = 0.186, 2-sided pvalue =1.484e-05

# Regional Trends ------------
## 30-foot regional richness and diversity through time ----------------
# see MS2_regionaltrends.R

# wrangle
bayregions2 <- bayregions %>%
  select(-c(midbay.simp,lowerbay.simp,upperbay.simp)) %>%
  mutate(across(c(midbay.shann,upperbay.shann,lowerbay.shann,midbay.richness,upperbay.richness,lowerbay.richness), as.numeric)) %>%
  mutate(Date = as.Date(Date)) %>%
  pivot_longer(cols = 2:7, names_to = "metric", values_to = "value") %>%
  mutate(type = case_when(
    metric == "midbay.shann" ~ "Diversity",
    metric == "lowerbay.shann" ~ "Diversity",
    metric == "upperbay.shann" ~ "Diversity",
    TRUE ~ "Richness"  # Default category if none of the above match
  ))

# plot the richness bayregions 
p1 <- bayregions2 %>%
  as.data.frame() %>%
  filter(type == "Richness") %>%
  ggplot(aes(x=Date,y=value, color = metric)) +
  geom_line(alpha=0.5) +
  geom_smooth(method="lm") +
  labs(x="Year",y="Species Richness", tag = "a)") +
  theme_light() +
  scale_color_manual(values=  c("#009E73","#F0E442","#56B4E9"),
                     name = "",
                     labels = c("Lower Bay", "Mid Bay", "Upper Bay")) 

# plot the diversity bayregions 
p2 <- bayregions2 %>%
  as.data.frame() %>%
  filter(type == "Diversity") %>%
  ggplot(aes(x=Date,y=value, color = metric)) +
  geom_line(alpha=0.5) +
  geom_smooth(method="lm") +
  labs(x="Year",y="Species Diversity") +
  theme_light() +
  scale_color_manual(values=  c("#009E73","#F0E442","#56B4E9"),
                     name = "",
                     labels = c("Lower Bay", "Mid Bay", "Upper Bay")) 

p1 / p2

### model regional trends through time ---------
# swap out 'lowerbay.richness' and 'lowerbay.shann' to look at richness and diversity
regional.glm <- bayregions2 %>%
  filter(metric == "lowerbay.shann") %>%
  mutate(month = month(Date)) %>%
  mutate(year = year(Date)) %>%
  select(-c(metric,type))

model <- glm(value ~ year + month, data = regional.glm)
summary(model)

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["year"] # 0.1881718 -- means increased XX per year on average 

# Calculate the total number of years
n_years <- max(regional.glm$year) - min(regional.glm$year)

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year 
start_year <- min(regional.glm$year)
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall # 11.73% increase overall 


## 16-foot Regional Richness ------------
small.regional.rich <- smalltrawl %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Month, Year4, Station, SpecCode) %>%
  summarise(number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  filter(Year4 %in% 1990:2019) %>%
  filter(!is.na(SpecCode)) %>%
  mutate(Date = as.Date(paste(Year4, Month, "1", sep="-"))) %>%
  select(-c(Month,Year4)) 

# upper bay 
sm.upper.rich <- small.regional.rich %>% 
  filter(Station %in% c(94,93,95,92,96,91,80,81)) %>%
  group_by(Date) %>% 
  tally() %>%
  rename(sm.upper.rich = n) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month"))

# mid bay 
sm.mid.rich <- small.regional.rich %>% 
  filter(Station %in% c(83,85,86,88,87,8,7,10,11,13,99,14,16,20,22,26,25,32)) %>%
  group_by(Date) %>% 
  tally() %>%
  rename(sm.mid.rich = n) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month"))

sm.lower.rich <- small.regional.rich %>% 
  filter(Station %in% c(17,21,33,40,41,90,49,48,55,56,61,62,63,64)) %>%
  group_by(Date) %>% 
  tally() %>%
  rename(sm.lower.rich = n) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month"))

## 16-foot Regional Diversity  --------------------------------
small.regional.div <- smalltrawl %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Month, Year4, Station, SpecCode) %>%
  summarise(Number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  filter(Year4 %in% 1990:2019) %>%
  filter(!is.na(SpecCode)) %>%
  mutate(Date = as.Date(paste(Year4, Month, "1", sep="-"))) %>%
  select(-c(Year4,Month)) 

# upper bay 
sm.upper.div <- small.regional.div %>% 
  filter(Station %in% c(94,93,95,92,96,91,80,81)) %>%
  group_by(Date, SpecCode) %>%
  summarise(Number = sum(Number, na.rm=T)) %>%
  pivot_wider(names_from = SpecCode, values_from = Number) %>%
  tibble::column_to_rownames('Date') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
sm.upper.div <- diversity(sm.upper.div, index="shannon")

sm.upper.div <- as.data.frame(sm.upper.div) %>%
  rownames_to_column(var = "Date")

# mid bay 
sm.mid.div <- small.regional.div %>% 
  filter(Station %in% c(83,85,86,88,87,8,7,10,11,13,99,14,16,20,22,26,25,32)) %>%
  group_by(Date, SpecCode) %>%
  summarise(Number = sum(Number, na.rm=T)) %>%
  pivot_wider(names_from = SpecCode, values_from = Number) %>%
  tibble::column_to_rownames('Date') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
sm.mid.div <- diversity(sm.mid.div, index="shannon")

sm.mid.div <- as.data.frame(sm.mid.div) %>%
  rownames_to_column(var = "Date")

# lower bay 
sm.lower.div <- small.regional.div %>% 
  filter(Station %in% c(17,21,33,40,41,90,49,48,55,56,61,62,63,64)) %>%
  group_by(Date, SpecCode) %>%
  summarise(Number = sum(Number, na.rm=T)) %>%
  pivot_wider(names_from = SpecCode, values_from = Number) %>%
  tibble::column_to_rownames('Date') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
sm.lower.div <- diversity(sm.lower.div, index="shannon")

sm.lower.div <- as.data.frame(sm.lower.div) %>%
  rownames_to_column(var = "Date")

# richness plot 
sm.rich.region.plot <- sm.upper.rich %>%
  left_join(sm.mid.rich, by = "Date") %>%
  left_join(sm.lower.rich, by = "Date") %>%
  mutate(Date = as.Date(Date)) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month")) %>%
  rename("Upper Bay" = sm.upper.rich, "Mid Bay" = sm.mid.rich, "Lower Bay" = sm.lower.rich) %>%
  pivot_longer(cols = 2:4, names_to = "Region", values_to = "Richness") %>%
  ggplot(aes(x=Date,y=Richness, color = Region)) +
  geom_line(alpha=0.5) +
  geom_smooth(method = "lm") +
  labs(x="Year",y="Species Richness", tag = "b)") +
  theme_light() +
  scale_color_manual(values=  c("#F0E442","#56B4E9","#E69F00"),
                     name = "",
                     labels = c("Mid Bay", "Upper Bay", "River")) 

# diversity 
sm.div.region.plot <- sm.upper.div %>%
  left_join(sm.mid.div, by = "Date") %>%
  left_join(sm.lower.div, by = "Date") %>%
  mutate(Date = as.Date(Date)) %>%
  complete(Date = seq.Date(min(Date), max(Date), by="month")) %>%
  rename("Upper Bay" = sm.upper.div, "Mid Bay" = sm.mid.div, "Lower Bay" = sm.lower.div) %>%
  pivot_longer(cols = 2:4, names_to = "Region", values_to = "Diversity") %>%
  ggplot(aes(x=Date,y=Diversity, color = Region)) +
    geom_line(alpha=0.5) +
    geom_smooth(method = "lm") +
    labs(x="Year",y="Species Diversity") +
    theme_light() +
  scale_color_manual(values=  c("#F0E442","#56B4E9","#E69F00"),
                     name = "",
                     labels = c("Mid Bay", "Upper Bay", "River")) 

sm.rich.region.plot / sm.div.region.plot

# FIGURE 4 - regional trends  -----
(p1 / p2) | (sm.rich.region.plot / sm.div.region.plot)

ggsave("Figures/FIGURE_4.png", 
       dpi=300, height=7, width=12, units='in')

### Rate of intro / departure ------------------------------------
smalltrawl <- smalltrawl %>%
  rename(Year = Year4,
         CommonName = Common.name) %>%
  filter(smalltrawl, ! is.na(CommonName)) %>%
  left_join(fishinfo.16, by = "CommonName")

heatmap <- smalltrawl %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Year, CommonName) %>%
  summarise(number=sum(CPUE,na.rm=TRUE))

# FIRST YEAR
firstyear <- heatmap %>%
  group_by(CommonName) %>%
  summarise(Year = min(Year))

ggplot(firstyear, aes(x=Year, y=reorder(CommonName, Year))) +
  geom_point() +
  theme(text = element_text(size=10)) +
  ylab("Species")

# LAST YEAR 
lastyear <- heatmap %>%
  group_by(CommonName) %>%
  summarise(Year = max(Year))

# rates 
introduction.rate <- firstyear %>%
  group_by(Year) %>%
  tally()

departure.rate <- lastyear %>%
  group_by(Year) %>%
  tally()

# cut out first few years, last few years 
introduction.rate <- filter(introduction.rate, Year > 1984)
departure.rate <- filter(departure.rate, Year < 2013)

# calc rates 
mean(introduction.rate$n) #3.5 (1995-)
sd(introduction.rate$n)  #2.382534

mean(departure.rate$n) #2.5 (-2015)
sd(departure.rate$n) #1.414214

# Seasonal Trends -----------
## Env. Trends -----------
# From SeasonalTemperatureChanges.R
seasonal <- New.Data %>%
  mutate(season = case_when(
    Month %in% 3:5 ~ "Spring",
    Month %in% 6:8 ~ "Summer",
    Month %in% 9:11 ~ "Fall",
    Month %in% c(1,2,12) ~ "Winter",
  )) %>%
  group_by(Year, season) %>%
  dplyr::summarize(bottom.temp = mean(BottomTemp, na.rm = TRUE),
            surface.temp = mean(SurfaceTemp, na.rm=T),
            bottom.DO = mean(BottomDO_sat, na.rm=T),
            surface.DO = mean(SurfaceDO_sat,na.rm=T)) %>%
  mutate(across(where(is.numeric), ~ na_if(., NaN)))

p1 <- seasonal %>%
  filter(season %in% c("Spring","Summer","Fall")) %>%
ggplot(aes(x=Year,y=bottom.temp, color = season)) +
  geom_line(alpha=0.5) +
  geom_smooth(method="lm") +
  labs(x="Year",y="Bottom Temperature (C)", tag = "a)") +
  theme_light() +
  scale_color_manual(values=colorblind_palette, name = "",labels = c("Fall", "Spring", "Summer"))

p2 <- seasonal %>%
  filter(season %in% c("Spring","Summer","Fall")) %>%
  ggplot(aes(x=Year,y=surface.temp, color = season)) +
  geom_line(alpha=0.5) +
  geom_smooth(method="lm") +
  labs(x="Year",y="Surface Temperature (C)", tag = "b)") +
  theme_light() +
  scale_color_manual(values=colorblind_palette, name = "",labels = c("Fall", "Spring", "Summer"))

p3 <- seasonal %>%
  filter(season %in% c("Spring","Summer","Fall")) %>%
  ggplot(aes(x=Year,y=bottom.DO, color = season)) +
  geom_line(alpha=0.5) +
  geom_smooth(method="lm") +
  labs(x="Year",y="Bottom Dissolved Oxygen", tag = "c)") +
  theme_light() +
  scale_color_manual(values=colorblind_palette, name = "", labels = c("Fall", "Spring", "Summer"))

p4 <- seasonal %>%
  filter(season %in% c("Spring","Summer","Fall")) %>%
  ggplot(aes(x=Year,y=surface.DO, color = season)) +
  geom_line(alpha=0.5) +
  geom_smooth(method="lm") +
  labs(x="Year",y="Surface Dissolved Oxygen", tag="d)") +
  theme_light() +
  scale_color_manual(values=colorblind_palette, name = "",labels = c("Fall", "Spring", "Summer"))

# FIGURE 2 - env. trends ----------

(p1 / p2) | (p3 / p4)

ggsave("Figures/FIGURE_2.png", 
       dpi=300, height=7, width=10, units='in')

### Calculate rate -----------------------
# using regression line 
fall <- seasonal %>%
  filter(season == "Fall")

summer <- seasonal %>%
  filter(season == "Summer")

spring <- seasonal %>%
  filter(season == "Spring")

# run regressions 

#summer
ssum <- lm(surface.temp ~ Year, data = summer)
summary(ssum) # yes

bsum <- lm(bottom.temp ~ Year, data = summer)
summary(bsum) #yes 

#fall 
bfall <- lm(bottom.temp ~ Year, data = fall)
summary(bfall) #no 

sfall <- lm(surface.temp ~ Year, data = fall)
summary(sfall) #yes #no 

#spring
bspring <- lm(bottom.temp ~ Year, data = spring)
summary(bspring) #no 

sspring <- lm(surface.temp ~ Year, data = spring)
summary(sspring) #no

# slopes 
ssum$coefficients #0.0427988
bsum$coefficients #0.02935573

sspring$coefficients #-0.009953095 
bspring$coefficients #0.01875844 

sfall$coefficients #0.06051394
bfall$coefficients #0.01625886 

#slope = rise/run = temp / year 
2019-1980
0.0427988 * 39

# surface = 39 years

#summer surface 
0.0427988*39 #1.669153

#spring surface
-0.009953095*39 #-0.3881707

#fall surface
0.06051394*39 #2.360044

# bottom = 52 years 

#summer bottom
0.03899307*52 #2.02764

#spring bottom
0.01875844*52 #0.9754389

#fall bottom 
0.01973277*52 #1.026104

## Seasonal Richness Trends -----------
### 30-foot survey -----------
rich.seasonal <- New.Data %>%
  filter(Year %in% 1990:2019) %>%
  mutate(season = case_when(
    Month %in% 3:5 ~ "Spring",
    Month %in% 6:8 ~ "Summer",
    Month %in% 9:11 ~ "Fall",
    Month %in% c(1,2,12) ~ "Winter",
  )) %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(season, Year, SpeciesCode) %>%
  summarise(number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  group_by(Year, season) %>% 
  tally() %>%
  rename(richness = n)

big.richness.seasonal <- rich.seasonal %>%
  filter(season %in% c("Spring","Summer","Fall")) %>%
  ggplot(aes(x=Year,y=richness, color = season)) +
  geom_line(size=0.5) +
  geom_smooth(method="lm") +
  geom_point(alpha=0.5)+
  labs(x="",y="Species Richness", tag = "a)") +
  theme_light() +
  scale_color_manual(values=colorblind_palette, name = "",labels = c("Fall", "Spring", "Summer"))

# using regression line 
fall <- rich.seasonal %>%
  filter(season == "Fall")

summer <- rich.seasonal %>%
  filter(season == "Summer")

spring <- rich.seasonal %>%
  filter(season == "Spring")

# run regressions 

#summer
summer.rich <- lm(richness ~ Year, data = summer)
summary(summer.rich) 

#fall 
fall.rich <- lm(richness ~ Year, data = fall)
summary(fall.rich) 

#spring
spring.rich <- lm(richness ~ Year, data = spring)
summary(spring.rich) 


# slopes 
summer.rich$coefficients #0.6137931
fall.rich$coefficients #0.8160178
spring.rich$coefficients #0.4669633

### 16-foot survey -----------
small.rich.seasonal <- smalltrawl %>%
  rename(Year = Year4) %>%
  filter(Year %in% 1990:2019) %>%
  mutate(season = case_when(
    Month %in% 3:5 ~ "Spring",
    Month %in% 6:8 ~ "Summer",
    Month %in% 9:11 ~ "Fall",
    Month %in% c(1,2,12) ~ "Winter",
  )) %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(season, Year, SpecCode) %>%
  summarise(number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  group_by(Year, season) %>% 
  tally() %>%
  rename(richness = n)

small.richness.seasonal <- small.rich.seasonal %>%
  filter(season %in% c("Spring","Summer","Fall")) %>%
  ggplot(aes(x=Year,y=richness, color = season)) +
  geom_line(size=0.5) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm") +
  labs(x="",y="Species Richness", tag = "b)") +
  theme_light() +
  scale_color_manual(values=colorblind_palette, name = "",labels = c("Fall", "Spring", "Summer"))

### FIGURE - seasonal richness -----------
big.richness.seasonal / small.richness.seasonal

## Seasonal Diversity ------------
### 30-foot seasonal diversity ------------

# Summer
big.summer.div <- New.Data %>%
  filter(Year %in% 1990:2019) %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  mutate(season = case_when(
    Month %in% 3:5 ~ "Spring",
    Month %in% 6:8 ~ "Summer",
    Month %in% 9:11 ~ "Fall",
    Month %in% c(1,2,12) ~ "Winter",
  )) %>%
  mutate(period = paste(season, Year, sep="-")) %>%
  filter(season == "Summer") %>%
  group_by(period,SpeciesCode) %>%
  summarise(Number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  filter(!is.na(SpeciesCode)) %>%
  pivot_wider(names_from = SpeciesCode, values_from = Number) %>%
  tibble::column_to_rownames('period') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
summer.shann <- diversity(big.summer.div, index="shannon")
summer.shann <- as.data.frame(summer.shann)
summer.shann$period = row.names(summer.shann)

# Fall
big.fall.div <- New.Data %>%
  filter(Year %in% 1990:2019) %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  mutate(season = case_when(
    Month %in% 3:5 ~ "Spring",
    Month %in% 6:8 ~ "Summer",
    Month %in% 9:11 ~ "Fall",
    Month %in% c(1,2,12) ~ "Winter",
  )) %>%
  mutate(period = paste(season, Year, sep="-")) %>%
  filter(season == "Fall") %>%
  group_by(period,SpeciesCode) %>%
  summarise(Number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  filter(!is.na(SpeciesCode)) %>%
  pivot_wider(names_from = SpeciesCode, values_from = Number) %>%
  tibble::column_to_rownames('period') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
fall.shann <- diversity(big.fall.div, index="shannon")
fall.shann <- as.data.frame(fall.shann)
fall.shann$period = row.names(fall.shann)

# Spring
big.spring.div <- New.Data %>%
  filter(Year %in% 1990:2019) %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  mutate(season = case_when(
    Month %in% 3:5 ~ "Spring",
    Month %in% 6:8 ~ "Summer",
    Month %in% 9:11 ~ "Fall",
    Month %in% c(1,2,12) ~ "Winter",
  )) %>%
  mutate(period = paste(season, Year, sep="-")) %>%
  filter(season == "Spring") %>%
  group_by(period,SpeciesCode) %>%
  summarise(Number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  filter(!is.na(SpeciesCode)) %>%
  pivot_wider(names_from = SpeciesCode, values_from = Number) %>%
  tibble::column_to_rownames('period') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
spring.shann <- diversity(big.spring.div, index="shannon")
spring.shann <- as.data.frame(spring.shann)
spring.shann$period = row.names(spring.shann)

#join to make a diversity table
big.seasonal.div <- summer.shann %>%
  separate(period, into = c("season", "year"), sep = "-") %>%
  select(-season) %>%
  left_join(fall.shann %>% separate(period, into = c("season", "year"), sep = "-"), by = "year") %>%
  select(-season) %>%
  left_join(spring.shann %>% separate(period, into = c("season", "year"), sep = "-"), by = "year")%>%
  select(-season) %>%
  rename("Fall" = fall.shann, "Summer" = summer.shann, "Spring" = spring.shann) %>%
  pivot_longer(cols = c(1,3,4),names_to = "Season", values_to = "Diversity") %>%
  mutate(year = as.numeric(year))

# plot 
big.seasonal.div %>%
  ggplot(aes(x=year,y=Diversity, color = Season)) +
  geom_line(size=0.5) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm") +
  labs(x="Year",y="Species Diversity", tag = "b)") +
  theme_light() +
  scale_color_manual(values=colorblind_palette, name = "",labels = c("Fall", "Spring", "Summer"))

### 16-foot seasonal diversity ------------

# Summer
small.summer.div <- smalltrawl %>%
  filter(Year4 %in% 1990:2019) %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  mutate(season = case_when(
    Month %in% 3:5 ~ "Spring",
    Month %in% 6:8 ~ "Summer",
    Month %in% 9:11 ~ "Fall",
    Month %in% c(1,2,12) ~ "Winter",
  )) %>%
  mutate(period = paste(season, Year4, sep="-")) %>%
  filter(season == "Summer") %>%
  group_by(period,SpecCode) %>%
  summarise(Number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  filter(!is.na(SpecCode)) %>%
  pivot_wider(names_from = SpecCode, values_from = Number) %>%
  tibble::column_to_rownames('period') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
summer.shann <- diversity(small.summer.div, index="shannon")
summer.shann <- as.data.frame(summer.shann)
summer.shann$period = row.names(summer.shann)

# Fall
small.fall.div <- smalltrawl %>%
  filter(Year4 %in% 1990:2019) %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  mutate(season = case_when(
    Month %in% 3:5 ~ "Spring",
    Month %in% 6:8 ~ "Summer",
    Month %in% 9:11 ~ "Fall",
    Month %in% c(1,2,12) ~ "Winter",
  )) %>%
  mutate(period = paste(season, Year4, sep="-")) %>%
  filter(season == "Fall") %>%
  group_by(period,SpecCode) %>%
  summarise(Number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  filter(!is.na(SpecCode)) %>%
  pivot_wider(names_from = SpecCode, values_from = Number) %>%
  tibble::column_to_rownames('period') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
fall.shann <- diversity(small.fall.div, index="shannon")
fall.shann <- as.data.frame(fall.shann)
fall.shann$period = row.names(fall.shann)

# Spring
small.spring.div <- smalltrawl %>%
  filter(Year4 %in% 1990:2019) %>%
  #filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  mutate(season = case_when(
    Month %in% 3:5 ~ "Spring",
    Month %in% 6:8 ~ "Summer",
    Month %in% 9:11 ~ "Fall",
    Month %in% c(1,2,12) ~ "Winter",
  )) %>%
  mutate(period = paste(season, Year4, sep="-")) %>%
  filter(season == "Spring") %>%
  group_by(period,SpecCode) %>%
  summarise(Number=sum(CPUE,na.rm=TRUE)) %>%
  ungroup() %>%
  filter(!is.na(SpecCode)) %>%
  pivot_wider(names_from = SpecCode, values_from = Number) %>%
  tibble::column_to_rownames('period') %>%
  mutate(across(everything(), ~ replace_na(., 0)))

#Shannon-Wiener Function 
spring.shann <- diversity(small.spring.div, index="shannon")
spring.shann <- as.data.frame(spring.shann)
spring.shann$period = row.names(spring.shann)

#join to make a diversity table
small.seasonal.div <- summer.shann %>%
  separate(period, into = c("season", "year"), sep = "-") %>%
  select(-season) %>%
  left_join(fall.shann %>% separate(period, into = c("season", "year"), sep = "-"), by = "year") %>%
  select(-season) %>%
  left_join(spring.shann %>% separate(period, into = c("season", "year"), sep = "-"), by = "year")%>%
  select(-season) %>%
  rename("Fall" = fall.shann, "Summer" = summer.shann, "Spring" = spring.shann) %>%
  pivot_longer(cols = c(1,3,4),names_to = "Season", values_to = "Diversity") %>%
  mutate(year = as.numeric(year))

# plot 
small.seasonal.div %>%
  ggplot(aes(x=year,y=Diversity, color = Season)) +
  geom_line(size=0.5) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm") +
  labs(x="Year",y="Species Diversity", tag = "b)") +
  theme_light() +
  scale_color_manual(values=colorblind_palette, name = "",labels = c("Fall", "Spring", "Summer"))

## FIGURE - seasonal diversity 
big.div.s <- big.seasonal.div %>%
  ggplot(aes(x=year,y=Diversity, color = Season)) +
  geom_line(size=0.5) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm") +
  labs(x="Year",y="Species Diversity", tag = "") +
  theme_light() +
  scale_color_manual(values=colorblind_palette, name = "",labels = c("Fall", "Spring", "Summer"))

small.div.s <- small.seasonal.div %>%
  ggplot(aes(x=year,y=Diversity, color = Season)) +
  geom_line(size=0.5) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm") +
  labs(x="Year",y="Species Diversity", tag = "") +
  theme_light() +
  scale_color_manual(values=colorblind_palette, name = "",labels = c("Fall", "Spring", "Summer"))

# FIGURE 3 - seasonal - rich & div -------------
(big.richness.seasonal / big.div.s) | (small.richness.seasonal/ small.div.s)

ggsave("Figures/FIGURE_3.png", 
       dpi=300, height=7, width=12, units='in')

# Species accumulation curve -------------------------------------
# from MS2_species_accumulation.R 
## 30-foot survey ----------
richness = NULL

for(i in unique(New.Data$Year)) {
  
  trawl <- New.Data %>%
    group_by(Trawl, SpeciesCode) %>%
    filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
    filter(Year == i) %>%
    summarise(number = sum(CPUE, na.rm = TRUE))
  
  trawl <- na.omit(trawl)
  
  trawl <- dcast(trawl, Trawl ~ SpeciesCode, value.var = 'number')
  rownames(trawl) <- trawl$Trawl
  trawl$Trawl = NULL
  trawl[is.na(trawl)] <- 0
  
  accurve <- specaccum(trawl, method="random", permutations=100)
  
  accuum <- cbind(accurve$sites, accurve$richness, i)
  
  richness <- rbind(accuum, richness)
  
  print(i)
}

richness <- richness %>%
  as.data.frame() %>%
  rename(Year = i) 

### Slope loop ------------------------------------------
slope = NULL

#make dataframe- station, number, lat, long 
for(i in unique(richness$Year)) {
  
  trawl <- filter(richness, Year == i)
  
  lm <- lm(trawl$V2 ~ log(trawl$V1))
  
  est <- summary(lm)$coefficients[2,1]
  
  se <- coef(summary(lm))[2,2]
  
  estimate <- as.data.frame(cbind(est, se, i))
  
  slope <- rbind(estimate, slope)
  
  print(i)
}

# make dataframe of slopes for each year 
slope <- as.data.frame(slope)

### Fit weighted quadratic curve ------------------------------------

fit2<-lm(est~poly(i,2,raw=TRUE), weights = 1/se, data = slope)

summary(fit2)

fit2$coefficient[1]
fit2$coefficient[2]
fit2$coefficient[3]

quadratic = fit2$coefficient[3]*slope$i^2 + fit2$coefficient[2]*slope$i + fit2$coefficient[1]
quadratic

# write the dataframes 
slope.df <- slope %>%
  cbind(quadratic) 

#write_csv(slope.df, "/Users/haleyoleynik/Documents/Thesis/Data/MS 2/30-foot-sp-accum-slope.csv") 
# write_csv(richness, "/Users/haleyoleynik/Documents/Thesis/Data/MS 2/30-foot-sp-accum-richness.csv") 
# FIGURE Sp. accum
big.richness <- read_csv("/Users/haleyoleynik/Documents/Thesis/Data/MS 2/30-foot-sp-accum-richness.csv")
big.slope <- read_csv("/Users/haleyoleynik/Documents/Thesis/Data/MS 2/30-foot-sp-accum-slope.csv")

big.SA.plot <- ggplot(big.richness) + 
  geom_point(aes(x=V1, y=V2, color = Year), size = 3) +
  labs(x = "Samples", y="Number of Species", tag = "a)") +
  theme_light() + 
  scale_color_gradient2(midpoint = 1990, low = "#009E73", high = "#E69F00", mid = "#F0E442")

big.slope.plot <- ggplot(big.slope) + 
  geom_point(aes(x=i, y=est), color = "#56B4E9") +
  geom_errorbar(aes(x=i, y=est,ymin=est-se, ymax=est+se), color = "#56B4E9") +
  geom_line(aes(i, quadratic), col="#56B4E9", size = 1) +
  geom_line(aes(x=i, y=est), alpha = 0.3, size = 1, color = "#56B4E9") +
  xlab("Year") +
  ylab("Slope") +
  theme_light() +
  scale_x_continuous(breaks = c(1970, 1980, 1990, 2000, 2010, 2019))

big.SA.plot / big.slope.plot

## 16-foot survey ----------
small.accu <- smalltrawl %>%
  rename(Year = Year4,
         CommonName = Common.name) %>%
  filter(! is.na(CommonName)) %>%
  filter(Year > 1979)

richness = NULL

for(i in unique(small.accu$Year)) {
  trawl <- small.accu %>%
    group_by(Trawl, SpecCode) %>%
    filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
    filter(Year == i) %>%
    summarise(number = sum(CPUE, na.rm = TRUE))
  trawl <- na.omit(trawl)
  trawl <- dcast(trawl, Trawl ~ SpecCode, value.var = 'number')
  rownames(trawl) <- trawl$Trawl
  trawl$Trawl = NULL
  trawl[is.na(trawl)] <- 0
  accurve <- specaccum(trawl, method="random", permutations=100)
  accuum <- cbind(accurve$sites, accurve$richness, i)
  richness <- rbind(accuum, richness)
  print(i)
}

richness <- as.data.frame(richness)
richness$Year <- richness$i
richness$i = NULL

### Slope loop ------------------------------------------
#make a loop for this 
slope = NULL

#make dataframe- station, number, lat, long 
for(i in unique(richness$Year)) {
  trawl <- filter(richness, Year == i)
  lm <- lm(trawl$V2 ~ log(trawl$V1))
  est <- summary(lm)$coefficients[2,1]
  se <- coef(summary(lm))[2,2]
  estimate <- as.data.frame(cbind(est, se, i))
  slope <- rbind(estimate, slope)
  print(i)
}

# make dataframe of slopes for each year 
slope <- as.data.frame(slope)

plot(slope$i, slope$est, type = "l")

#plot slopes over time 
ggplot(slope, aes(x=i, y=est)) + 
  geom_point(col = "darkblue") +
  geom_errorbar(aes(ymin=est-se, ymax=est+se), col = "darkblue") +
  xlab("Year") +
  ylab("Slope") +
  theme(text = element_text(size=30))

### Fit weighted quadratic curve ------------------------------------

fit2<-lm(est~poly(i,2,raw=TRUE), weights = 1/se, data = slope)

summary(fit2)

fit2$coefficient[1]
fit2$coefficient[2]
fit2$coefficient[3]

quadratic = fit2$coefficient[3]*slope$i^2 + fit2$coefficient[2]*slope$i + fit2$coefficient[1]
quadratic


### FIGURE 16-foot sp. accum
small.SA <- ggplot(richness) + 
  geom_point(aes(x=V1, y=V2, color = Year), size = 3) +
  labs(x = "Samples", y="Number of Species", tag = "b)") +
  theme_light() + 
  scale_color_gradient2(midpoint = 2000, low = "#009E73", high = "#E69F00", mid = "#F0E442") 

small.slope <- ggplot(slope) + 
  geom_point(aes(x=i, y=est), col = "#56B4E9") +
  geom_errorbar(aes(x=i, y=est,ymin=est-se, ymax=est+se), col = "#56B4E9") +
  geom_line(aes(i, quadratic), col="#56B4E9", size = 1) +
  geom_line(aes(x=i, y=est), alpha = 0.3, size = 1, col="#56B4E9") +
  xlab("Year") +
  ylab("Slope") +
  theme_light() +
  scale_x_continuous(breaks = c(1970, 1980, 1990, 2000, 2010, 2018))

small.SA / small.slope

# FIGURE 5 - sp. accumulation ----------------------
(big.SA.plot / big.slope.plot) | (small.SA / small.slope)

ggsave("Figures/FIGURE_5.png", 
       dpi=300, height=7, width=10, units='in')

# taken from Month_PresenceAbsence_Th.R

# Occurrences through time --------------------------------

## 30-foot survey --------
## First half occurrences
sp.early <- New.Data %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Year, Month, CommonName) %>%
  summarise(number=sum(CPUE,na.rm=TRUE)) %>%
  filter(Year < 1990) %>%
  mutate(Month.Year = as.Date(paste(Year, Month,"1", sep="-"))) %>%
  ungroup() %>%
  select(-c(Month,Year)) %>%
  pivot_wider(names_from = Month.Year, values_from = number) %>%
  mutate(across(2:113, ~ replace_na(., 0))) %>%
  mutate(across(2:113, ~ ifelse(. > 0, 1, .))) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(total.early = sum(c_across(2:113))) %>%
  select(CommonName, total.early) %>%
  mutate(prop.early = total.early/112) %>% 
  filter(prop.early > 0.05) %>%
  filter(prop.early < 0.95)

# Second half occurrences
sp.late <- New.Data %>%
  filter(Year > 1989) %>% # 1990-2019
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Year, Month, CommonName) %>%
  summarise(number=sum(CPUE,na.rm=TRUE)) %>% 
  mutate(Month.Year = as.Date(paste(Year, Month,"1", sep="-"))) %>%
  ungroup() %>%
  select(-c(Month,Year)) %>%
  pivot_wider(names_from = Month.Year, values_from = number) %>%
  mutate(across(2:288, ~ replace_na(., 0))) %>%
  mutate(across(2:288, ~ ifelse(. > 0, 1, .))) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(total.late = sum(c_across(2:288))) %>%
  select(CommonName, total.late) %>%
  mutate(prop.late = total.late/287) %>% 
  filter(prop.late > 0.05) %>%
  filter(prop.late < 0.95)

# Combine 
sp.all <- sp.early %>%
  right_join(sp.late, by = "CommonName") %>%
  mutate(diff = prop.late - prop.early)  %>% # to sort by differences 
  select(CommonName, diff, prop.early,prop.late) %>%
  pivot_longer(cols = 3:4, names_to = "variable", values_to = "value")

## FIGURE - 30-foot occurrences ------------- 
big.oc <- ggplot(sp.all, aes(fill=variable, y=value, x=reorder(CommonName,diff))) + 
  geom_bar(position="dodge", stat="identity") + 
  coord_flip() + 
  labs(x="Species",y="Proportion of Months", tag = "a)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="",
                    labels=c("First Half", "Second Half"),
                    values = colorblind_palette)

## 16-foot survey -----------------------
smalltrawl <- smalltrawl %>%
  rename(Year = Year4,
         CommonName = Common.name) %>%
  filter(! is.na(CommonName)) %>%
  left_join(fishinfo.16, by = "CommonName")

small.early <- smalltrawl %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Year, Month, CommonName) %>%
  summarise(number=sum(CPUE,na.rm=TRUE)) %>%
  filter(Year < 1990) %>% # 1966-1984
  mutate(Month.Year = as.Date(paste(Year, Month,"1", sep="-"))) %>%
  ungroup() %>%
  select(-c(Month,Year)) %>%
  pivot_wider(names_from = Month.Year, values_from = number) %>%
  mutate(across(2:69, ~ replace_na(., 0))) %>%
  mutate(across(2:69, ~ ifelse(. > 0, 1, .))) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(total.early = sum(c_across(2:69))) %>%
  select(CommonName, total.early) %>%
  mutate(prop.early = total.early/68) %>% 
  filter(prop.early > 0.05) %>%
  filter(prop.early < 0.95)

# Second half
small.late <- smalltrawl %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Year, Month, CommonName) %>%
  summarise(number=sum(CPUE,na.rm=TRUE)) %>%
  filter(Year > 1989) %>%
  mutate(Month.Year = as.Date(paste(Year, Month,"1", sep="-"))) %>%
  ungroup() %>%
  select(-c(Month,Year)) %>%
  pivot_wider(names_from = Month.Year, values_from = number) %>%
  mutate(across(2:195, ~ replace_na(., 0))) %>%
  mutate(across(2:195, ~ ifelse(. > 0, 1, .))) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(total.late = sum(c_across(2:195))) %>%
  select(CommonName, total.late) %>%
  mutate(prop.late = total.late/194) %>% 
  filter(prop.late > 0.05) %>%
  filter(prop.late < 0.95)

# Combine 
small.all <- small.early %>%
  right_join(small.late, by = "CommonName") %>%
  mutate(diff = prop.late - prop.early)  %>% # to sort by differences 
  select(CommonName, diff, prop.early,prop.late) %>%
  pivot_longer(cols = 3:4, names_to = "variable", values_to = "value")

## FIGURE - 16-foot occurrences ------------- 
sm.oc <- ggplot(small.all, aes(fill=variable, y=value, x=reorder(CommonName,diff))) + 
  geom_bar(position="dodge", stat="identity") + 
  coord_flip() + 
  labs(x="Species",y="Proportion of months", tag = "b)") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(name="",
                    labels=c("First Half", "Second Half"),
                    values = colorblind_palette)

# FIGURE 6 - occurences ------
big.oc | sm.oc

ggsave("Figures/FIGURE_6.png", 
       dpi=300, height=7, width=10, units='in')

# Fish Base Metrics -------------------------
# look at MS2_community_metrics and MS2_length_freq 
# MS2_occurences

# big trawl 
fish <- New.Data %>% 
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>%
  group_by(Trawl, Year, Month, CommonName) %>%
  summarise(number = sum(CPUE, na.rm=T)) %>%
  left_join(fishinfo.30, by = "CommonName")

# small trawl 
small.fish <- smalltrawl %>% 
  group_by(Trawl, Year4, Month, Common.name) %>%
  summarise(number = sum(CPUE, na.rm=T)) %>%
  left_join(fishinfo.16, by = c("Common.name" = "CommonName")) %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish")

### 30-foot Environment -------------------------------------------------------------
try <- fish %>%
  group_by(Year, CommonName, Environment) %>%
  summarise(number = sum(number, na.rm=T)) %>%
  ungroup() %>%
  filter(!is.na(Environment)) %>%
  group_by(Year, Environment) %>%
  summarise(number = sum(number,na.rm=T)) %>%
  pivot_wider(values_from = number, names_from = Environment) %>%
  rowwise() %>%
  mutate(across(benthopelagic:`reef-associated`, ~ . / sum(c_across(benthopelagic:`reef-associated`)))) %>%
  ungroup() %>%
  pivot_longer(cols = 2:5, values_to = "number", names_to = "Environment")

# stacked bar 
p1 <- ggplot(data = try, aes(x = Year, y = number, fill = factor(Environment, levels=c("demersal", "pelagic", "benthopelagic", "reef-associated")))) + 
  geom_bar(stat = 'identity') + 
  ylab("Relative CPUE") +
  guides(fill=guide_legend(title="Environment")) +
  scale_fill_manual(values=colorblind_palette,name = "Environment", labels = c("Demersal", "Pelagic", "Benthopelagic", "Reef-Associated")) +
  theme_light() +
  theme(text = element_text(size=10))

p2 <- try %>%
  ggplot(aes(x=Year, y = number, color = factor(Environment, levels=c("demersal", "pelagic", "benthopelagic", "reef-associated")))) +
  geom_line() +
  geom_smooth(method = "lm") +
  xlab("Year") +
  ylab("Relative CPUE") +
  facet_wrap(vars(Environment),nrow=2, scales="free_y") +
  scale_color_manual(values = colorblind_palette, name = "Environment", labels = c("Demersal", "Pelagic", "Benthopelagic", "Reef-Associated")) +
  theme_light() +
  theme(text = element_text(size=10), strip.background = element_blank())

p1 | p2

### Model environment -------
env <- fish %>%
  group_by(Year, CommonName, Environment) %>%
  summarise(number = sum(number, na.rm=T)) %>%
  ungroup() %>%
  filter(!is.na(Environment)) %>%
  group_by(Year, Environment) %>%
  summarise(number = sum(number,na.rm=T)) %>%
  pivot_wider(values_from = number, names_from = Environment) %>%
  rowwise() %>%
  mutate(across(benthopelagic:`reef-associated`, ~ . / sum(c_across(benthopelagic:`reef-associated`)))) %>%
  ungroup()

model <- glm(`reef-associated` ~ Year, data = env)
summary(model)

### 16-foot Environment -------------------------------------------------------------
try <- small.fish %>%
  group_by(Year4, Common.name, Environment) %>%
  summarise(number = sum(number, na.rm=T)) %>%
  ungroup() %>%
  filter(!is.na(Environment)) %>%
  filter(Environment %in% c("demersal", "pelagic")) %>%
  group_by(Year4, Environment) %>%
  summarise(number = sum(number,na.rm=T)) %>%
  pivot_wider(values_from = number, names_from = Environment) %>%
  mutate(across(everything(), ~ replace_na(., 0))) %>%
  rowwise() %>%
  mutate(across(demersal:pelagic, ~ . / sum(c_across(demersal:pelagic), na.rm = TRUE))) %>%
  ungroup() %>%
  pivot_longer(cols = 2:3, values_to = "number", names_to = "Environment")

# stacked bar 
p1 <- ggplot(data = try, aes(x = Year4, y = number, fill = factor(Environment, levels=c("demersal", "pelagic")))) + 
  geom_bar(stat = 'identity') + 
  ylab("Relative CPUE") +
  guides(fill=guide_legend(title="Environment")) +
  scale_fill_manual(values=colorblind_palette,name = "Environment", labels = c("Demersal", "Pelagic")) +
  theme_light() +
  theme(text = element_text(size=10))

p2 <- try %>%
  ggplot(aes(x=Year4, y = number, color = factor(Environment, levels=c("demersal", "pelagic")))) +
  geom_line() +
  geom_smooth(method = "lm") +
  xlab("Year") +
  ylab("Relative CPUE") +
  facet_wrap(vars(Environment),nrow=2, scales="free_y") +
  scale_color_manual(values = colorblind_palette, name = "Environment", labels = c("Demersal", "Pelagic")) +
  theme_light() +
  theme(text = element_text(size=10), strip.background = element_blank())

p1 | p2


### 30-foot Feeding Habit ----------------------------------------------------------
fish2 <- fish %>%
  group_by(Year, CommonName, Feeding.Habit2) %>%
  summarise(number = sum(number, na.rm=T)) %>%
  ungroup() %>%
  mutate(Feeding.Habit2 = na_if(Feeding.Habit2, "")) %>%
  filter(!is.na(Feeding.Habit2)) %>%
  group_by(Year, Feeding.Habit2) %>%
  summarise(number = sum(number,na.rm=T)) %>%
  pivot_wider(values_from = number, names_from = Feeding.Habit2) %>%
  mutate(across(everything(), ~ replace_na(., 0))) %>%
  rowwise() %>%
  mutate(across(`grazing on aquatic plants`:`filtering plankton`, ~ . / sum(c_across(`grazing on aquatic plants`:`filtering plankton`))))  %>%
  ungroup() %>%
  pivot_longer(cols = `grazing on aquatic plants`:`filtering plankton`, values_to = "number", names_to = "feeding.habit")

# stacked bar 
f1 <- fish2 %>%
  filter(feeding.habit != "grazing on aquatic plants") %>%
ggplot(aes(x = Year, y = number, fill = factor(feeding.habit, levels=c("hunting macrofauna (predator)", "filtering plankton", "selective plankton feeding", "variable")))) + 
  geom_bar(stat = 'identity') + 
  ylab("Relative CPUE") +
  scale_fill_manual(values = colorblind_palette,name = "Feeding Habit", labels = c("Predator", "Filter Feeder", "Planktivore", "Variable")) +
  theme_light() +
  theme(text = element_text(size=10))

f2 <- fish2 %>%
  filter(feeding.habit != "grazing on aquatic plants") %>%
ggplot(aes(x=Year, y = number, color = feeding.habit)) +
  geom_line() +
  geom_smooth(method = "lm") +
  xlab("Year") +
  ylab("Relative CPUE") +
  facet_wrap(vars(feeding.habit),nrow=2, scales="free_y") +
  scale_color_manual(values = colorblind_palette,name = "Feeding Habit", labels = c("Predator", "Filter Feeder", "Planktivore", "Variable")) +
  theme_light() +
  theme(text = element_text(size=10), strip.background = element_blank())

f1 | f2

### 16-foot Feeding Habit ----------------------------------------------------------
fish2 <- small.fish %>%
  group_by(Year4, Common.name, Feeding.Habit) %>%
  summarise(number = sum(number, na.rm=T)) %>%
  ungroup() %>%
  mutate(Feeding.Habit = na_if(Feeding.Habit, "")) %>%
  filter(Feeding.Habit %in% c("predator", "variable")) %>%
  filter(!is.na(Feeding.Habit)) %>%
  group_by(Year4, Feeding.Habit) %>%
  summarise(number = sum(number,na.rm=T)) %>%
  pivot_wider(values_from = number, names_from = Feeding.Habit) %>%
  mutate(across(everything(), ~ replace_na(., 0))) %>%
  rowwise() %>%
  mutate(across(predator:variable, ~ . / sum(c_across(predator:variable))))  %>%
  ungroup() %>%
  pivot_longer(cols = predator:variable, values_to = "number", names_to = "feeding.habit")

# stacked bar 
f1 <- fish2 %>%
  ggplot(aes(x = Year4, y = number, fill = factor(feeding.habit, levels=c("predator", "variable")))) + 
  geom_bar(stat = 'identity') + 
  ylab("Relative CPUE") +
  scale_fill_manual(values = colorblind_palette,name = "Feeding Habit", labels = c("Predator", "Variable")) +
  theme_light() +
  theme(text = element_text(size=10))

f2 <- fish2 %>%
  ggplot(aes(x=Year4, y = number, color = feeding.habit)) +
  geom_line() +
  geom_smooth(method = "lm") +
  xlab("Year") +
  ylab("Relative CPUE") +
  facet_wrap(vars(feeding.habit),nrow=2, scales="free_y") +
  scale_color_manual(values = colorblind_palette,name = "Feeding Habit", labels = c("Predator", "Variable")) +
  theme_light() +
  theme(text = element_text(size=10), strip.background = element_blank())

f1 | f2

# FIGURE 8 & 8b - env. and feeding habits -----------
(p1 | p2) / (f1 | f2)

ggsave("Figures/FIGURE_8.png", 
       dpi=300, height=7, width=12, units='in')
ggsave("Figures/FIGURE_8b.png", 
       dpi=300, height=7, width=12, units='in')

### 30-foot Trophic Level -----------------------------
TL <- fish %>% 
  group_by(Year, CommonName) %>% 
  summarise(TL = mean(Trophic.Level, na.rm=T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(mean.TL = mean(TL, na.rm=T),
            #sd.TL = sd(TL, na.rm=T),
            se.TL = sd(TL, na.rm=T) / sqrt(n())) 

### Model trophic level ----- 
model <- glm(TL ~ Year, data = TL)
summary(model) # not signif. 

### 30-foot Latitude -----------------------------------------
#range <- read.csv("/Users/haleyoleynik/Documents/Thesis/Data/New_SpeciesRanges.csv")
#range <- range[,-1]

# USE NEW Spreadsheets!!!!!!!!!!!!!!!1
latitude <-  fish %>% 
  group_by(Year, CommonName) %>% 
  summarise(upperlat= mean(upperlat, na.rm=T),
            lowerlat= mean(lowerlat, na.rm=T))  %>%
  mutate(midlat = upperlat - ((upperlat-lowerlat)/2)) %>%
  ungroup() %>%
  group_by(Year) %>% 
  summarise(mean.lowerlat = mean(lowerlat,na.rm=T),
            #sd.lowerlat = sd(lowerlat, na.rm=T),
            se.lowerlat = sd(lowerlat, na.rm=T) / sqrt(n()),
            mean.midlat = mean(midlat,na.rm=T),
            #sd.midlat = sd(midlat, na.rm=T),
            se.midlat = sd(midlat, na.rm=T) / sqrt(n()),
            mean.upperlat = mean(upperlat,na.rm=T),
            #sd.upperlat = sd(upperlat, na.rm=T),
            se.upperlat = sd(upperlat, na.rm=T) / sqrt(n()))

### Model latitude ----- 
model <- glm(mean.midlat ~ Year, data = latitude)
summary(model) # signif

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["Year"] # -0.06566188 - means decreased .34 per year on average 

# Calculate the total number of years
n_years <- max(latitude$Year) - min(latitude$Year)

# Get the coefficient of year
coef_year <- coef(model)["Year"]

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year (e.g., 1990)
start_year <- min(latitude$Year)
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall # -12.26295% decrease overall 

# lowerlat
model <- glm(mean.lowerlat ~ Year, data = latitude)
summary(model) # signif

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["Year"] # -0.1305788  -- means increased .34 per year on average 

# Calculate the total number of years
n_years <- max(latitude$Year) - min(latitude$Year)

# Get the coefficient of year
coef_year <- coef(model)["Year"]

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year (e.g., 1990)
start_year <- min(latitude$Year)
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall #-88.1169 % decrease overall 

### 30-foot Preferred Temperature  ----------------------------
# try confidence intervals -- calculate SD (or SE)
temp <- fish %>% 
  group_by(Year, CommonName) %>% 
  summarise(temp = mean(Temp.mean, na.rm=T),
            maxtemp = mean(Temp.max, na.rm=T),
            mintemp = mean(Temp.min, na.rm=T)) %>%
  ungroup() %>%
  group_by(Year) %>% 
  filter(!is.nan(temp)) %>%
  summarise(mean.temp = mean(temp, na.rm=T),
            #sd.temp = sd(temp, na.rm=T),
            se.temp = sd(temp, na.rm=T) / sqrt(n()),
            mean.maxtemp = mean(maxtemp, na.rm=T),
            #sd.maxtemp = sd(maxtemp, na.rm=T),
            se.maxtemp = sd(maxtemp, na.rm=T) / sqrt(n()),
            mean.mintemp = mean(mintemp, na.rm=T),
            #sd.mintemp = sd(mintemp, na.rm=T),
            se.mintemp = sd(mintemp, na.rm=T) / sqrt(n()), # Calculate standard deviation
            .groups = "drop")

### Model preferred temp ----- 
model <- glm(mean.temp ~ Year, data = temp)
summary(model) # signif

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["Year"] # 0.3439364 -- means increased .34 per year on average 

# Calculate the total number of years
n_years <- max(temp$Year) - min(temp$Year)

# Get the coefficient of year
coef_year <- coef(model)["Year"]

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year (e.g., 1990)
start_year <- min(temp$Year)
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall # 10.68071% increase overall 

temp <- fish %>% 
  group_by(Year, CommonName) %>% 
  summarise(temp = mean(Temp.mean, na.rm=T),
            maxtemp = mean(Temp.max, na.rm=T),
            mintemp = mean(Temp.min, na.rm=T)) %>%
  ungroup() %>%
  group_by(Year) %>% 
  filter(!is.nan(temp)) %>%
  summarise(mean.temp = mean(temp, na.rm=T),
            #sd.temp = sd(temp, na.rm=T),
            se.temp = sd(temp, na.rm=T) / sqrt(n()),
            mean.maxtemp = mean(maxtemp, na.rm=T),
            #sd.maxtemp = sd(maxtemp, na.rm=T),
            se.maxtemp = sd(maxtemp, na.rm=T) / sqrt(n()),
            mean.mintemp = mean(mintemp, na.rm=T),
            #sd.mintemp = sd(mintemp, na.rm=T),
            se.mintemp = sd(mintemp, na.rm=T) / sqrt(n()), # Calculate standard deviation
            .groups = "drop")

model <- glm(mean.mintemp ~ Year, data = temp)
summary(model) # signif

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["Year"] # 0.3439364 -- means increased .34 per year on average 

# Calculate the total number of years
n_years <- max(temp$Year) - min(temp$Year)

# Get the coefficient of year
coef_year <- coef(model)["Year"]

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year (e.g., 1990)
start_year <- min(temp$Year)
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall 

# mean temp 10.68071% increase overall 
# maxtemp 5.697029% increase overall 
# mintemp 13.50025% increase overall 

### 30-foot length ---------------
## From MS2_Length_Analysis.R
dat <- select(New.Data, SpeciesCode, CommonName, Group) %>%
  distinct() %>%
  rename(SpecCode = SpeciesCode)

# lengths data from MS2_Length_Freq.R code 
lengths <- read.csv("/Users/haleyoleynik/Documents/Thesis/Data/Length_Data.csv")

# separate out 
len <- separate(lengths, col = Trawl, into = c("Year", "Month", "Day", "Time", "Station"), sep = "-") %>%
  mutate(Year = as.numeric(Year)) %>%
  right_join(dat, lengths, by = "SpecCode") %>%
  filter(length < 999) %>%
  filter(!is.na(CommonName)) %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>% # JUST FISH 
  uncount(freq) %>% # repeat rows by frequency
  group_by(Year, CommonName) %>% 
  summarise(length = mean(length, na.rm=T)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(mean.length = mean(length,na.rm=T),
            #sd.length = sd(length,na.rm=T),
            se.length = sd(length,na.rm=T)/sqrt(n()))

### Model length ----- 
model <- glm(mean.length ~ Year, data = len)
summary(model) # signif

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["Year"] # -0.932753 -- means increased X per year on average 

# Calculate the total number of years
n_years <- max(len$Year) - min(len$Year)

# Get the coefficient of year
coef_year <- coef(model)["Year"]

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year (e.g., 1990)
start_year <- min(len$Year)
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent change
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall 

# -19.20922% decrease  overall 

#### Length distribution ------------
# Model distribution through time 
dat <- select(New.Data, SpeciesCode, CommonName, Group) %>%
  distinct() %>%
  rename(SpecCode = SpeciesCode)

# lengths data from MS2_Length_Freq.R code 
lengths <- read.csv("/Users/haleyoleynik/Documents/Thesis/Data/Length_Data.csv")

# separate out 
length.dist <- separate(lengths, col = Trawl, into = c("Year", "Month", "Day", "Time", "Station"), sep = "-") %>%
  mutate(Year = as.numeric(Year)) %>%
  right_join(dat, lengths, by = "SpecCode") %>%
  filter(length < 999) %>%
  filter(!is.na(CommonName)) %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>% # JUST FISH 
  uncount(freq) %>% # repeat rows by frequency
  group_by(Year, CommonName) %>% 
  summarise(length = mean(length, na.rm=T)) %>%
  ungroup() 

# RAW LENGTH
p1 <- ggplot(length.dist, aes(x = length, y = as.factor(Year))) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) +
  labs(
    #title = "Distribution of Lengths by Year",
    x = "Length (mm)",
    y = "Year"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# LOG LENGTH
p2 <- ggplot(length.dist, aes(x = log(length), y = as.factor(Year))) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) +
  labs(
    #title = "Distribution of Lengths by Year",
    x = "log(Length)",
    y = "Year"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

p1 | p2 

ggsave("Figures/length_distributions.png", 
       dpi=300, height=10, width=7, units='in')

# Calculate mean log(length) for each year
mean_lengths <- length.dist %>%
  group_by(Year) %>%
  summarise(mean_log_length = mean(log(length), na.rm = TRUE))

# Plot with a single point for the mean length in each year
ggplot(length.dist, aes(x = log(length), y = as.factor(Year))) +
  geom_density_ridges(alpha = 0.7, scale = 1.2) +
  geom_point(data = mean_lengths, aes(x = mean_log_length, y = as.factor(Year)), 
             color = "red", size = 2, shape = 18) + # Using a small dot or shape for clarity
  labs(
    title = "Distribution of Lengths by Year",
    x = "Length (mm)",
    y = "Year"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


x %>%
  group_by(Year) %>%
  summarise(mean.length = mean(length,na.rm=T),
            #sd.length = sd(length,na.rm=T),
            se.length = sd(length,na.rm=T)/sqrt(n()))

### 16-foot Trophic Level -----------------------------
small.TL <- small.fish %>% 
  rename(Year = Year4, CommonName = Common.name) %>% 
  group_by(Year, CommonName) %>% 
  summarise(TL = mean(Trophic.Level, na.rm=T)) %>% 
  ungroup() %>%
  group_by(Year) %>% 
  summarise(mean.TL = mean(TL, na.rm=T),
            #sd.TL = sd(TL, na.rm=T),
            se.TL = sd(TL, na.rm=T) / sqrt(n())) 

### Model trophic level ----- 
model <- glm(mean.TL ~ Year, data = small.TL)
summary(model) # not signif. 

### 16-foot Latitude ----------------------------------------
small.latitude <-  small.fish %>% 
  rename(Year = Year4, CommonName = Common.name) %>% 
  group_by(Year, CommonName) %>% 
  summarise(upperlat= mean(upperlat, na.rm=T),
            lowerlat= mean(lowerlat, na.rm=T))  %>%
  mutate(midlat = upperlat - ((upperlat-lowerlat)/2)) %>%
  ungroup() %>%
  group_by(Year) %>% 
  summarise(mean.lowerlat = mean(lowerlat,na.rm=T),
            #sd.lowerlat = sd(lowerlat, na.rm=T),
            se.lowerlat = sd(lowerlat, na.rm=T) / sqrt(n()),
            mean.midlat = mean(midlat,na.rm=T),
            #sd.midlat = sd(midlat, na.rm=T),
            se.midlat = sd(midlat, na.rm=T) / sqrt(n()),
            mean.upperlat = mean(upperlat,na.rm=T),
            #sd.upperlat = sd(upperlat, na.rm=T),
            se.upperlat = sd(upperlat, na.rm=T) / sqrt(n()))

### Model Latitude ----- 
model <- glm(mean.upperlat ~ Year, data = small.latitude)
summary(model) 

# Get the coefficient for the year (the long-term trend)
coef_year <- coef(model)["Year"] # -0.06566188 - means decreased .34 per year on average 

# Calculate the total number of years
n_years <- max(small.latitude$Year) - min(small.latitude$Year)

# Get the coefficient of year
coef_year <- coef(model)["Year"]

# Calculate the total temperature change over the time series
total_change <- coef_year * n_years

# Get the intercept (for predicted starting temperature)
intercept <- coef(model)["(Intercept)"]

# Estimate the predicted temperature in the starting year (e.g., 1990)
start_year <- min(small.latitude$Year)
predicted_temp_start <- intercept + coef_year * start_year

# Calculate the overall percent increase
percent_increase_overall <- ((predicted_temp_start + total_change) / predicted_temp_start - 1) * 100

percent_increase_overall # -12.26295% decrease overall 

### 16-foot Preferred Temperature  ----------------------------
# weight by abundance? 
small.temp <- small.fish %>% 
  rename(Year = Year4, CommonName = Common.name) %>% 
  group_by(Year, CommonName) %>% 
  summarise(temp = mean(Temp.mean, na.rm=T),
            maxtemp = mean(Temp.max, na.rm=T),
            mintemp = mean(Temp.min, na.rm=T)) %>%
  ungroup() %>%
  group_by(Year) %>% 
  summarise(mean.temp = mean(temp, na.rm=T),
            #sd.temp = sd(temp, na.rm=T),
            se.temp = sd(temp, na.rm=T) / sqrt(n()),
            mean.maxtemp = mean(maxtemp, na.rm=T),
            #sd.maxtemp = sd(maxtemp, na.rm=T),
            se.maxtemp = sd(maxtemp, na.rm=T) / sqrt(n()),
            mean.mintemp = mean(mintemp, na.rm=T),
            #sd.mintemp = sd(mintemp, na.rm=T),
            se.mintemp = sd(mintemp, na.rm=T) / sqrt(n()), # Calculate standard deviation
            .groups = "drop")

### Model Preferred temp ----- 
model <- glm(mean.temp ~ Year, data = small.temp)
summary(model) 


## 16-foot Lengths ------------------
smalltrawl.2 <- smalltrawl %>%
  rename(Year = Year4,
         CommonName = Common.name) %>%
  filter(! is.na(CommonName))

# read length data 
small.lengths <- read_csv("/Users/haleyoleynik/Documents/Thesis/Data/Length_Data_small-trawl.csv")

small.len <- smalltrawl.2 %>%
  select(SpecCode, CommonName, Group) %>%
  distinct() %>%
  na.omit() %>%
  mutate(SpecCode = as.double(SpecCode)) %>%
  right_join(small.lengths, by = "SpecCode") %>%
  filter(Length < 999) %>%
  filter(!is.na(CommonName)) %>%
  filter(Group == "Fin Fish" | Group == "Cartilaginous Fish") %>% # JUST FISH 
  uncount(freq) %>% # repeat rows by frequency
  group_by(Year, CommonName) %>% 
  summarise(length = mean(Length, na.rm=T)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(mean.length = mean(length,na.rm=T),
            #sd.length = sd(length,na.rm=T),
            se.length = sd(length,na.rm=T)/sqrt(n()))

### Model length ----- 
model <- glm(mean.length ~ Year, data = small.len)
summary(model) 

# FIGURE 7 - 30-foot fishbase metrics -------------
all <- temp %>%
  left_join(TL, by = "Year") %>%
  left_join(latitude, by = "Year") %>%
  left_join(len, by = c("Year" = "Year")) %>%
  pivot_longer(cols = 2:17, names_to = "Metric", values_to = "Value") %>%
  mutate(Stat = case_when(
    Metric %in% c("mean.temp", "mean.maxtemp","mean.mintemp","mean.TL","mean.lowerlat","mean.midlat","mean.upperlat","mean.length") ~ "Mean",
    #Metric %in% c("sd.temp", "sd.maxtemp","sd.mintemp","sd.TL","sd.lowerlat","sd.midlat","sd.upperlat","sd.length") ~ "SD",
    Metric %in% c("se.temp", "se.maxtemp","se.mintemp","se.TL","se.lowerlat","se.midlat","se.upperlat","se.length") ~ "SE"
  ),
  Group = case_when(
    Metric %in% c("mean.temp", "mean.maxtemp","mean.mintemp",
                  #"sd.temp", "sd.maxtemp","sd.mintemp",
                  "se.temp", "se.maxtemp","se.mintemp") ~ "Temperature",
    Metric %in% c("mean.TL", 
                  #"sd.TL", 
                  "se.TL") ~ "Trophic Level",
    Metric %in% c("mean.lowerlat","mean.midlat","mean.upperlat",
                  #"sd.lowerlat","sd.midlat","sd.upperlat", 
                  "se.lowerlat","se.midlat","se.upperlat") ~ "Latitude",
    Metric %in% c("mean.length", 
                  #"sd.length", 
                  "se.length") ~ "Length"
  )) %>%
  mutate(Metric = sub("mean\\.|se\\.", "", Metric)) %>%
  pivot_wider(names_from = "Stat", values_from = "Value")

# plot 
all %>%
  ggplot(aes(x = Year, y = Mean, color = Group)) +
  geom_line() +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), 
                width = 0.2) +
  facet_wrap(Group ~ Metric, scales = "free_y", nrow = 2) +
  scale_color_manual(values = colorblind_palette) + 
  theme_light() 

ggsave("Figures/FIGURE_7.png", 
       dpi=300, height=7, width=14, units='in')

# FIGURE 7b - 16-foot fishbase metrics -------------

small.all <- small.temp %>%
  left_join(small.TL, by = "Year") %>%
  left_join(small.latitude, by = "Year") %>%
  left_join(small.len, by = c("Year" = "Year")) %>%
  pivot_longer(cols = 2:17, names_to = "Metric", values_to = "Value") %>%
  mutate(Stat = case_when(
    Metric %in% c("mean.temp", "mean.maxtemp","mean.mintemp","mean.TL","mean.lowerlat","mean.midlat","mean.upperlat","mean.length") ~ "Mean",
    #Metric %in% c("sd.temp", "sd.maxtemp","sd.mintemp","sd.TL","sd.lowerlat","sd.midlat","sd.upperlat","sd.length") ~ "SD",
    Metric %in% c("se.temp", "se.maxtemp","se.mintemp","se.TL","se.lowerlat","se.midlat","se.upperlat","se.length") ~ "SE"
  ),
  Group = case_when(
    Metric %in% c("mean.temp", "mean.maxtemp","mean.mintemp",
                  #"sd.temp", "sd.maxtemp","sd.mintemp",
                  "se.temp", "se.maxtemp","se.mintemp") ~ "Temperature",
    Metric %in% c("mean.TL", 
                  #"sd.TL", 
                  "se.TL") ~ "Trophic Level",
    Metric %in% c("mean.lowerlat","mean.midlat","mean.upperlat",
                  #"sd.lowerlat","sd.midlat","sd.upperlat", 
                  "se.lowerlat","se.midlat","se.upperlat") ~ "Latitude",
    Metric %in% c("mean.length", 
                  #"sd.length", 
                  "se.length") ~ "Length"
  )) %>%
  mutate(Metric = sub("mean\\.|se\\.", "", Metric)) %>%
  pivot_wider(names_from = "Stat", values_from = "Value")

# plot 
small.all %>%
  ggplot(aes(x = Year, y = Mean, color = Group)) +
  geom_line() +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), 
                width = 0.2) +
  facet_wrap(Group ~ Metric, scales = "free_y", nrow = 2) +
  scale_color_manual(values = colorblind_palette) + 
  theme_light() 

ggsave("Figures/FIGURE_7b.png", 
       dpi=300, height=7, width=14, units='in')

# GAM -------------------------
# written data !! read this and jump to time series plots
#data <- read.csv("/Users/haleyoleynik/Documents/Thesis/Data/MS 2/ALL_GLM_DATA_v2")
data <- read.csv("/Users/haleyoleynik/Documents/Thesis/Data/MS 2/ALL_GLM_DATA_v3")
data$Date <- as.Date(data$Date)
fssi <- read_csv("/Users/haleyoleynik/Documents/Thesis/Data/MS 2/Atlantic_FSSI_timeseries.csv")

climate <- read_csv("/Users/haleyoleynik/Documents/Thesis/Data/Climate Environment/Climate_Indices.csv")
climate$date <- as.Date(climate$date)


# combine 
new.data <- fssi %>%
  dplyr::select(Date, FSSI = score) %>%
  right_join(data, by = "Date") %>%
  right_join(climate, by = c('Date' = 'date')) %>%
  arrange(Date)

# FIGURE - Plot timeseries and relationships

p1 <- ggplot(new.data, aes(Date,richness)) +
  geom_line() +
  labs(x = "Year", y = "Richness") +
  theme_light()

p2 <- ggplot(new.data, aes(Date,FSSI)) +
  geom_line(col = "#E69F00") +
  labs(x = "Year", y = "FSSI") +
  theme_light()

p3 <- ggplot(new.data, aes(Date,temp)) +
  geom_line(col = "#56B4E9") +
  labs(x = "Year", y = "Temp (C)") +
  theme_light()

p4 <- ggplot(new.data, aes(Date,DO)) +
  geom_line(col = "#009E73") +
  labs(x = "Year", y = "Dissolved Oxygen") +
  theme_light()

p5 <- ggplot(new.data, aes(Date,all.trips.extrapolated)) +
  geom_line(col = "#CC79A7") +
  labs(x = "Year", y = "Trips") +
  theme_light()

p6 <- ggplot(new.data, aes(Date,River.DO)) +
  geom_line(col = "#D55E00") +
  labs(x = "Year", y = "River Dissolved Oxygen") +
  theme_light()


(p1 | p2 | p3) / (p4 | p5 | plot_spacer())

# relationships 
t1 <- ggplot(new.data, aes(FSSI,richness)) +
  geom_point(col = "#E69F00") +
  geom_smooth(method = "lm", col = "#E69F00",se=F) +
  labs(x = "FSSI", y = "Richness") +
  theme_light()

t2 <- ggplot(new.data, aes(temp,richness)) +
  geom_point(col = "#56B4E9") +
  geom_smooth(method = "lm", col = "#56B4E9",se=F) +
  labs(x = "Temp (C)", y = "Richness") +
  theme_light()

t3 <- ggplot(new.data, aes(DO,richness)) +
  geom_point(col = "#009E73") +
  geom_smooth(method = "lm", col = "#009E73",se=F) +
  labs(x = "Dissolved Oxygen", y = "Richness") +
  theme_light()

t4 <- ggplot(new.data, aes(all.trips.extrapolated, richness)) +
  geom_point(col = "#CC79A7") +
  geom_smooth(method = "lm", col = "#CC79A7",se=F) +
  labs(x = "Trips", y = "Richness") +
  theme_light()

t5 <- ggplot(new.data, aes(River.DO, richness)) +
  geom_point(col = "#D55E00") +
  geom_smooth(method = "lm", col = "#D55E00",se=F) +
  labs(x = "River Dissolved Oxygen", y = "Richness") +
  theme_light()

(p1 | p2 | p3) / (p4 | p5 | plot_spacer())
(t1 | t2) / (t3 | t4 )


# FIGURE  9 - GAM ------ 
(p1 ) /
  (p2 | t1) /
  (p3 | t2) /
  (p4 | t3) /
  (p5 | t4) /
  (p6 | t5)

ggsave("Figures/FIGURE_9.png", 
       dpi=300, height=12, width=7, units='in')

# clean up data pt. 2 - take out FMP amendments 
data <- new.data %>% 
  dplyr::select(Date, Month, Year, richness, all.trips.extrapolated, temp, DO, FSSI)

## Collinearity / Concurvity -----
#define multiple linear regression model
model <- lm(richness ~ Month + temp + FSSI + all.trips.extrapolated + DO , data=data)

#calculate the VIF for each predictor variable in the model
car::vif(model)

# could include both DOs, river do as a metric of freshwater quality, and do as a measure of bay quality, but river DO is strongly correlated (R2=0.55) to temperature, high VIF scores when River DO is included. 

## GLM ----------------------------------
# try GLM one by one 
GLM.temp <- glm(richness ~ temp + Month + Year, family = poisson, data = new.data)
GLM.DO <- glm(richness ~ DO + Month, family = poisson, data = new.data)
GLM.fssi<- glm(richness ~ FSSI + Month, family = poisson, data = new.data)
GLM.all.trips.ex <- glm(richness ~ all.trips.extrapolated + Month, family = poisson, data = new.data)

summary(GLM.temp)
summary(GLM.DO)
summary(GLM.fssi)
summary(GLM.all.trips.ex)

GLM <- glm(richness ~ temp + all.trips.extrapolated + DO + Month + Year, family = poisson, data = new.data)
summary(GLM)
glm_model <- glm(richness ~ temp + Year + factor(Month), data = new.data)
summary(glm_model)

## GAM --------------
# plot richness over time 
data %>%
  ggplot(aes(x=Date, y=richness)) +
  geom_point() +
  geom_smooth()

# just look at temp data 
temp <- new.data %>%
  select(c(richness,temp,Year,Month)) %>%
  na.omit() %>%
  filter(Year > 1989)

length(unique(temp$Month))
length(unique(temp$Year))
table(temp$Month)

# shape of month
temp %>%
  group_by(Month) %>%
  summarise(temp = mean(temp,na.rm=T)) %>%
  ggplot(aes(x=Month, y=temp)) +
  geom_line()

# shape of richness
temp %>%
  group_by(Month) %>%
  summarise(richness = mean(richness,na.rm=T)) %>%
  ggplot(aes(x=Month, y=richness)) +
  geom_line()

## TEMP sensitivities ---------------------- 
GAM.temp <- gam(richness ~ temp + Year + s(Month, bs = "cc", k = 10), family = poisson, data = temp)
summary(GAM.temp) # deviance explained = 46.1
coef(GAM.temp)
plot(GAM.temp)

GAM.temp2 <- gam(richness ~ s(temp) + s(Year) + s(Month, bs = "cc", k = 10), family = poisson, data = temp)
summary(GAM.temp2) # 62.9
coef(GAM.temp2)
plot(GAM.temp2)

GAM.temp2 <- gam(richness ~ s(temp) + Year + s(Month, bs = "cc", k = 10), family = poisson, data = temp)
summary(GAM.temp2) # 52
coef(GAM.temp2)
plot(GAM.temp2)

GAM.temp2 <- gam(richness ~ temp + s(Year) + s(Month, bs = "cc", k = 10), family = poisson, data = temp)
summary(GAM.temp2) # 57.5
coef(GAM.temp2)
plot(GAM.temp2)

# all covariates ****** no significance 

GAM1.temp <- gam(richness ~ temp + Year + s(Month, bs = "cc", k = 10), family = poisson, data = new.data)
summary(GAM1.temp) # deviance explained = 44.3
coef(GAM1.temp)
plot(GAM1.temp)

GAM1.trips <- gam(richness ~ all.trips.extrapolated + Year + s(Month, bs = "cc", k = 10), family = poisson, data = new.data)
summary(GAM1.trips) # deviance explained = 41.6 -- not significant! 

GAM1.DO <- gam(richness ~ DO + Year + s(Month, bs = "cc", k = 10), family = poisson, data = new.data)
summary(GAM1.DO) # deviance explained = 43.7 -- not significant! 

ggplot(new.data, aes(x=Date, y=River.DO)) +
  geom_line()

GAM1.RiverDO <- gam(richness ~ River.DO + Year + s(Month, bs = "cc", k = 10), family = poisson, data = new.data)
summary(GAM1.RiverDO) # deviance explained = 43.7 -- not significant! 

GAM1.temp <- gam(richness ~ temp + Year + s(Month, bs = "cc", k = 10), family = poisson, data = new.data)
summary(GAM1.temp)

# model all together 
GAM1.all <- gam(richness ~ temp + all.trips.extrapolated + DO + Year + s(Month, bs = "cc", k = 10), family = poisson, data = new.data)
summary(GAM1.all)

# sensitivity -- temp & year as a smooth term
GAM2.temp <- gam(richness ~ s(temp) + s(Year) + s(Month, bs = "cc", k = 10), family = poisson, data = new.data)
summary(GAM2.temp) 
coef(GAM2.temp)
plot(GAM2.temp)

# R-sq.(adj) =  0.354   Deviance explained = 49.2%
# AIC - 241.1226


# without FSSI 
data.2 <- new.data %>% 
  select(-FSSI) %>%
  na.omit()

GAM2 <- gam(richness ~ s(Year, bs = "cr", k = 27) + s(Month, bs = "cc", k = 10) + s(temp) + s(all.trips.extrapolated) + s(DO), family = poisson, data = data.2, na.action = "na.fail")

summary(GAM2)

# check individual 
GAM.DO <- gam(richness ~ s(Year, bs = "cr", k = 27) + s(Month, bs = "cc", k = 10) + s(temp), family = poisson, data = data.2, na.action = "na.fail")

summary(GAM.DO)

# try DO, FMP, rec effort as a linear term 
GAM1.2 <- gam(richness ~ s(Year, bs = "cr", k = 33) + s(Month, bs = "cc", k = 10) + s(temp) + FMP + all.trips.extrapolated + DO , family = poisson, data = data, na.action = "na.fail")

# try year as a factor 
GAM1.5 <- gam(richness ~ factor(Year) + s(Month, bs = "cc", k = 10) + s(temp) + s(FMP) + s(all.trips.extrapolated) + s(DO) , family = poisson, data = data, na.action = "na.fail")

summary(GAM1)
summary(GAM1.2)
summary(GAM1.5)

plot(GAM1.5)
gam.fit(GAM1)

plot(GAM1, rug = TRUE, residuals = TRUE,
     pch = 1, cex = 1)

#diagnostics 
gam.check(GAM1, old.style	= T)
qq.gam(GAM1)
AIC(GAM1) #851.6154
# checks concurity, anaolg to collinearity in a linear model 
concurvity(GAM1, full = TRUE)
concurvity(GAM1, full = F)

# Method 2 ********************************************************
# Model model month with a sine and cosine term

data$costerm = cos(2*pi*(data$Month)/10)
data$sinterm = sin(2*pi*(data$Month)/10)

data %>%
  ggplot(aes(x=Year)) +
  geom_line(aes(y=costerm)) +
  geom_line(aes(y=sinterm), col = "red") +
  theme_classic()

GAM2 <- gam(richness ~ costerm + sinterm + s(temp) + s(FMP) + s(all.trips.extrapolated) + s(DO), family = poisson, data = data, na.action = "na.fail")

GAM2.5 <- gam(richness ~ costerm + sinterm + s(FMP) + s(temp) + s(all.trips.extrapolated) + s(DO), family = poisson, data = data, na.action = "na.fail")


summary(GAM2)
summary(GAM2.5)

plot(GAM2)
gam.fit(GAM2)

#diagnostics 
gam.check(GAM2)  
AIC(GAM2) #AIC: 857.4236

# Method 3 ****************************************************************
# make Year and Month factors 
GAM.poi <- gam(richness ~ factor(Year) + factor(Month) + s(FMP) + s(temp) + s(all.trips) + s(DO), family = poisson, data = data, na.action = "na.fail")

summary(GAM.poi)
plot(GAM.poi)

## Environment ONLY model  ---------------------
# Method 1 *********************************************************
# Model month cubic cyclic spline

GAM.env <- gam(richness ~ s(Year, bs = "cr", k = 33) + s(Month, bs = "cc", k = 10) + s(temp) + s(DO) , family = poisson, data = data, na.action = "na.fail")

summary(GAM.env)
residuals.gam(GAM.env)

# add residuals to data 
df <- data %>%
  cbind(residuals.gam(GAM.env)) %>%
  rename("Residuals" = "residuals.gam(GAM.env)")

# plot 
df %>%
  ggplot(aes(x=Date, y=Residuals)) +
  geom_point() +
  geom_smooth() +
  theme_classic()

# check cross correlation 
ccf(df$all.trips.extrapolated, df$Residuals, lag.max = 48)
ccf(df$all.trips.extrapolated, df$richness, lag.max = 48)

length(unique(df$Year))

## Model selection -----------------------------
# Step wise model selection via null space penalization https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/gam.selection.html

# GAM 1 
mods1 <- dredge(GAM1,rank="AIC")  
top.mods1 <- get.models(mods1,cumsum(weight) <=.95)
mods.av1 <- model.avg(top.mods1,rank="AIC")
summary(mods.av1)
top.mods1 # more detail on specific model fits--top at the script output is the AIC best model.
summary(model.avg(get.models(mods1,cumsum(weight) <=0.95),rank="AIC")) # summary on a model.avg object gives variable importance
View(mods1)

# lowest aic = month + Year + Temp 
GAM1.3 <- gam(richness ~ s(Year, bs = "cr", k = 34) + s(Month, bs = "cc", k = 10) + s(temp), family = poisson, data = data, na.action = "na.fail")

summary(GAM1.3)
AIC(GAM1.3)

# GAM 2 
mods2 <- dredge(GAM2,rank="AIC") 
top.mods2 <- get.models(mods2,cumsum(weight) <=.95)
mods.av2 <- model.avg(top.mods2,rank="AIC")
summary(mods.av2)
top.mods2 # more detail on specific model fits--top at the script output is the AIC best model.
summary(model.avg(get.models(mods2,cumsum(weight) <=0.95),rank="AIC")) # summary on a model.avg object gives variable importance
View(mods2)

# lowest AIC = month + Year + Temp 
GAM1.2 <- gam(richness ~ s(Year, bs = "cr", k = 34) + s(Month, bs = "cc", k = 10) + s(temp), family = poisson, data = data, na.action = "na.fail")

summary(GAM1.2)
AIC(GAM1.2)

# Species Turnover ----------------------
# https://search.r-project.org/CRAN/refmans/codyn/html/turnover.html







