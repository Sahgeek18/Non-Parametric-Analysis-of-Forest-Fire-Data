library(tseries)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stats)
library(evd)
library(fitdistrplus)
library(corrplot)

df <- read.csv('forestfires.csv')
attach(df)
colnames(df)

summary(df)
sample(df, size = 5)
random_sample <- df[sample(1:nrow(df), 7, replace=FALSE), ]
random_sample
#--------------------------------------------------------------------------------
#------------------------Runs test for randomness--------------------------------
#--------------------------------------------------------------------------------

# 1. Temperature (high vs low)
# H0: High and low temperature days occur randomly.
# H1: High and low temperature days do not occur randomly.
temp_binary <- ifelse(temp > median(temp, na.rm = TRUE), 1, 0)
runs.test(as.factor(temp_binary))
# p-value = 1.66e-09 : Reject

#------------------------------------
# 2. Relative Humidity (high vs low)
# H0: High and low RH days occur randomly.
# H1: High and low RH days do not occur randomly.
RH_binary <- ifelse(RH > median(RH, na.rm = TRUE), 1, 0)
runs.test(as.factor(RH_binary))
# p-value = 0.004124 : Reject

#------------------------------------
# 3. Wind speed (high vs low)
# H0: High and low wind days occur randomly.
# H1: High and low wind days do not occur randomly.
wind_binary <- ifelse(wind > median(wind, na.rm = TRUE), 1, 0)
runs.test(as.factor(wind_binary))
#p-value = 0.007003 : Reject

#------------------------------------
# 4. Burned area (high vs low)
# H0: High and low burned area days occur randomly.
# H1: High and low burned area days do not occur randomly.
area_binary <- ifelse(df$area > median(df$area, na.rm = TRUE), 1, 0)
runs.test(as.factor(area_binary))
# p-value < 2.2e-16 : Reject

#------------------------------------
# 5. Rain occurrence
# H0: Rain/no-rain days occur randomly.
# H1: Rain/no-rain days do not occur randomly.
rain_binary <- ifelse(rain > 0, 1, 0)
runs.test(as.factor(rain_binary))
# p-value < 2.2e-16 : Reject


## Conclusion: All p-values are significantly less than 0.05. Thus the sequence
## of data are not random, there is a systematic pattern or trend or clustering 
## in the variable over time


#------------------------------------
# Visual representation

plot_data <- df %>%
  dplyr::select(temp, RH, wind, area, rain) %>%
  mutate(obs = 1:n()) %>%
  tidyr::pivot_longer(cols = c(temp, RH, wind, area, rain),
                      names_to = "Variable",
                      values_to = "Value")

ggplot(plot_data, aes(x = obs, y = Value, color = Variable)) +
  geom_line(size = 0.9) +
  facet_wrap(~Variable, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Forest Fire Dataset Variables Over Observations",
    x = "Observation Index",
    y = "Value"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )


#--------------------------------------------------------------------------------
#------------------------Chi sq GoF Test-----------------------------------------
#--------------------------------------------------------------------------------

## H₀: Forest fires are uniformly distributed across months.
## H₁: Forest fires are not uniformly distributed across months.
# Frequency table
month_table <- table(month)

# Chi-square goodness-of-fit test (expected uniform distribution)
chisq.test(month_table)
# p-value < 2.2e-16 : Reject


## Conclusion: p-value is significantly less than 0.05, thus The number 
## of fires varies significantly by month, indicating a seasonal pattern 
## in fire occurrences.

# Bar plot of observed frequencies
ggplot(data.frame(month = names(month_table), freq = as.numeric(month_table)),
       aes(x = month, y = freq, fill = month)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal(base_size = 13) +
  labs(title = "Distribution of Fires Across Months",
       x = "Month", y = "Frequency") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))

#--------------------------------------------------------------------------------
#----------------------------KS GoF Test-----------------------------------------
#--------------------------------------------------------------------------------

## H₀: The burnt area follows a shifted exponential distribution with parameters 
## x0 and lambda
## H₁: The burnt area does not follow a shifted exponential distribution

# Keep only positive areas
area_pos <- df$area[df$area > 0]

# Estimate shift (minimum observed area)
x0 <- min(area_pos)

# Estimate rate parameter for shifted exponential
lambda_hat <- 1 / (mean(area_pos) - x0)

# Define shifted exponential CDF
pexp_shifted <- function(x) {
  ifelse(x < x0, 0, pexp(x - x0, rate = lambda_hat))
}

# Perform one-sample Kolmogorov–Smirnov test
ks.test(area_pos, pexp_shifted)
# p-value < 2.2e-16 : Reject

# Visual representation
ggplot(data.frame(area_pos), aes(x = area_pos)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightgreen", color = "black") +
  stat_function(fun = function(x) dexp(x - x0, rate = lambda_hat),
                color = "red", size = 1.2) +
  labs(title = "Shifted Exponential Fit to Burnt Area",
       x = "Burnt Area (ha)", y = "Density") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#------------------------

# 1. Temperature
## H₀: Temperature follows a normal distribution.
## H₁: Temperature does not follow a normal distribution.

ks.test(temp, "pnorm", mean(temp, na.rm = TRUE), sd(temp, na.rm = TRUE))
# p-value = 0.1392 : Accept

# Visualization
ggplot(df, aes(x = temp)) +
  geom_histogram(aes(y = ..density..), bins = 25, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm,
                args = list(mean = mean(temp, na.rm = TRUE),
                            sd = sd(temp, na.rm = TRUE)),
                color = "red", size = 1.2) +
  theme_minimal(base_size = 13) +
  labs(title = "K–S Test: Temperature vs Normal Distribution",
       x = "Temperature (°C)", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#-------------------------------

# 2. Relative Humidity

## H₀: Relative humidity follows a normal distribution.
## H₁: Relative humidity does not follow a normal distribution.

ks.test(RH, "pnorm", mean(RH, na.rm = TRUE), sd(RH, na.rm = TRUE))
# p-value = 1.725e-05 : Reject

# Visualization
ggplot(df, aes(x = RH)) +
  geom_histogram(aes(y = ..density..), bins = 25, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm,
                args = list(mean = mean(RH, na.rm = TRUE),
                            sd = sd(RH, na.rm = TRUE)),
                color = "red", size = 1.2) +
  theme_minimal(base_size = 13) +
  labs(title = "K–S Test: RH vs Normal Distribution",
       x = "RH", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#------------------------

# 3. Wind Speed

## H₀: Wind speed follows a normal distribution.
## H₁: Wind speed does not follow a normal distribution.

ks.test(wind, "pnorm", mean(wind, na.rm = TRUE), sd(wind, na.rm = TRUE))
#p-value = 6.485e-05 : Reject

# Visualization
ggplot(df, aes(x = wind)) +
  geom_histogram(aes(y = ..density..), bins = 25, fill = "skyblue", color = "black") +
  stat_function(fun = dnorm,
                args = list(mean = mean(wind, na.rm = TRUE),
                            sd = sd(wind, na.rm = TRUE)),
                color = "red", size = 1.2) +
  theme_minimal(base_size = 13) +
  labs(title = "K–S Test: wind vs Normal Distribution",
       x = "wind", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



## H₀: The minimum relative humidity values follow a Gumbel (Extreme Value Type
## I) distribution.
## H₁: The minimum relative humidity values do not follow a Gumbel distribution.
# Step 1: Extract monthly minimum RH
RH_monthly_min <- df %>%
  group_by(month) %>%
  summarise(min_RH = min(RH, na.rm = TRUE))

# Step 2: Fit Gumbel distribution (shape = 0 for Gumbel)
fit_gumbel_RH <- fgev(RH_monthly_min$min_RH, shape = 0)
summary(fit_gumbel_RH)

# Step 3: Extract estimated parameters
mu_RH <- fit_gumbel_RH$estimate["loc"]
sigma_RH <- fit_gumbel_RH$estimate["scale"]

# Step 4: Define Gumbel CDF for K–S test
pgumbel <- function(x, mu, sigma) exp(-exp(-(x - mu) / sigma))

# Step 5: Perform K–S test
ks.test(RH_monthly_min$min_RH, pgumbel, mu_RH, sigma_RH)
# p-value = 0.5962 : Accept

ggplot(RH_monthly_min, aes(x = min_RH)) +
  geom_histogram(aes(y = ..density..),
                 bins = 10, fill = "lightblue", color = "black") +
  stat_function(fun = function(x) dgev(x, loc = mu_RH, scale = sigma_RH, shape = 0),
                color = "red", size = 1.2) +
  labs(title = "Extreme Value (Gumbel) Fit — Monthly Minimum Relative Humidity",
       x = "Minimum RH (%)",
       y = "Density") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


## H₀: The maxima of burnt areas follow an Extreme Value (Gumbel) distribution.
## H₁: The maxima do not follow an Extreme Value (Gumbel) distribution.

area_monthly_max <- df %>%
  group_by(month) %>%
  summarise(max_area = max(area, na.rm = TRUE))


# Fit Gumbel distribution to monthly maxima
fit_gumbel_area <- fgev(area_monthly_max$max_area, shape = 0)  # Gumbel = shape = 0
summary(fit_gumbel_area)

# Extract parameters
mu <- fit_gumbel_area$estimate["loc"]
sigma <- fit_gumbel_area$estimate["scale"]

# Define CDF for Gumbel
pgumbel <- function(x, mu, sigma) exp(-exp(-(x - mu) / sigma))

ks.test(area_monthly_max$max_area, pgumbel, mu, sigma)
# p-value = 0.04348

ggplot(data.frame(area_monthly_max), aes(x = max_area)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "orange", color = "black") +
  stat_function(fun = function(x) dgev(x, loc = mu, scale = sigma, shape = 0),
                color = "red", size = 1.2) +
  labs(title = "Extreme Value Fit (Gumbel) - Monthly Max Burnt Area",
       x = "Max Burnt Area (ha)", y = "Density") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


#===================================================================
#                  Two Sample NON-PARAMETRIC ANALYSIS PROJECT
#===================================================================

df['fire'] <- ifelse(df$area == 0, 'No', 'Yes')

sept_data <- df[df$month == 'sep', ]
aug_data <- df[df$month == 'aug', ]


#===================================================================
#                         LOCATION TESTS
#===================================================================

#---------------- Temperature vs Fire Occurrence -------------------

plot(density(df[df$fire == 'Yes', 'temp']), ylim = c(0, 0.08),
     col = 'green', lwd = 2,
     main = "Density Comparison: Temperature (Fire vs No Fire)",
     xlab = "Temperature")
lines(density(df[df$fire == 'No', 'temp']), col = 'pink', lwd = 2)
legend("topright", legend = c("Fire", "No Fire"),
       col = c("green", "pink"), lwd = 2)


# H0: Median temperature is equal for days with and without burned area
# H1: Median temperature is higher on days with burned area

wilcox.test(df[df$fire == 'Yes', 'temp'], df[df$fire == 'No', 'temp'],
            alternative = 'greater')

# p-value = 0.007353 < 0.05 ⇒ Reject H0.
# Conclusion: Days with burned forest area have significantly higher median temperatures.
# Implication: Elevated temperatures increase forest fire risk and burned area extent.




#---------------- DMC vs Fire Occurrence ----------------------------

plot(density(df[df$fire == 'Yes', 'DMC']),
     col = 'green', lwd = 2,
     main = "Density Comparison: DMC (Fire vs No Fire)",
     xlab = "DMC")
lines(density(df[df$fire == 'No', 'DMC']), col = 'pink', lwd = 2)
legend("topright", legend = c("Fire", "No Fire"),
       col = c("green", "pink"), lwd = 2)

# H0: Median DMC (Duff Moisture Code) is equal for days with and without burned area
# H1: Median DMC is higher on days with burned area

wilcox.test(df[df$fire == 'Yes', 'DMC'], df[df$fire == 'No', 'DMC'],
            alternative = 'greater')

# p-value = 0.03675 < 0.05 ⇒ Reject H0.
# Conclusion: Days with burned area exhibit significantly higher median DMC values.
# Implication: Drier duff layer conditions (higher DMC) increase vulnerability to forest fires.




#---------------- DC vs Month (August vs September) -----------------

plot(density(aug_data$DC), ylim = c(0, 0.02),
     col = 'green', lwd = 2,
     main = "Density Comparison: DC (August vs September)",
     xlab = "DC")
lines(density(sept_data$DC), col = 'pink', lwd = 2)
legend("topright", legend = c("August", "September"),
       col = c("green", "pink"), lwd = 2)

# H0: Median DC (Drought Code) in September is less than or equal to August
# H1: Median DC in September is greater than in August

wilcox.test(aug_data$DC, sept_data$DC, alternative = 'less')

# p-value < 2.2e-16 ⇒ Reject H0.
# Conclusion: September has significantly higher median DC values than August.
# Implication: Progressive seasonal drying leads to increased deep layer drought by September, elevating fire danger.



#---------------- Wind vs Month (August vs September) ---------------

plot(density(aug_data$wind), ylim = c(0, 0.3),
     col = 'green', lwd = 2,
     main = "Density Comparison: Wind (August vs September)",
     xlab = "Wind Speed")
lines(density(sept_data$wind), col = 'pink', lwd = 2)
legend("topright", legend = c("August", "September"),
       col = c("green", "pink"), lwd = 2)

# H0: Median wind speed in August is less than or equal to September
# H1: Median wind speed in August is greater than in September

wilcox.test(aug_data$wind, sept_data$wind, alternative = 'greater')

# p-value = 0.0004959 < 0.05 ⇒ Reject H0.
# Conclusion: August exhibits significantly higher median wind speeds than September.
# Implication: Stronger August winds can rapidly spread fires despite lower drought conditions.



#---------------- Temperature vs Month (August vs September) --------

plot(density(aug_data$temp), ylim = c(0, 0.3),
     col = 'green', lwd = 2,
     main = "Density Comparison: Temperature (August vs September)",
     xlab = "Temperature (°C)")
lines(density(sept_data$temp), col = 'pink', lwd = 2)
legend("topright", legend = c("August", "September"),
       col = c("green", "pink"), lwd = 2)

# H0: Median temperature in August is less than or equal to September
# H1: Median temperature in August is greater than in September

wilcox.test(aug_data$temp, sept_data$temp, alternative = 'greater')

# p-value = 1.308e-05 < 0.05 ⇒ Reject H0.
# Conclusion: August has significantly higher median temperatures than September.
# Implication: Combined high August temperatures and wind, with September drought, create distinct fire risk profiles for each month.



#===================================================================
#                         SCALE TESTS
#===================================================================

#---------------- DC vs Fire Occurrence -----------------------------

plot(density(df[df$fire == 'Yes', 'DC']),
     col = 'green', lwd = 2,
     main = "Density Comparison: DC (Fire vs No Fire)",
     xlab = "DC")
lines(density(df[df$fire == 'No', 'DC']), col = 'pink', lwd = 2)
legend("topright", legend = c("Fire", "No Fire"),
       col = c("green", "pink"), lwd = 2)

# H0: Variance of DC is greater or equal for burned area days versus no burned area
# H1: Variance of DC is less for burned area days

mood.test(df[df$fire == 'Yes', 'DC'], df[df$fire == 'No', 'DC'],
          alternative = 'less')

# p-value = 0.3451 > 0.05 ⇒ Fail to reject H0.
# Conclusion: DC variability shows no significant difference between days with and without burned area.
# Implication: Drought Code fluctuations are similar regardless of fire outcome; level matters more than variability.



#---------------- FFMC vs Month (August vs September) ---------------

plot(density(aug_data$FFMC), ylim = c(0, 0.3),
     col = 'green', lwd = 2,
     main = "Density Comparison: FFMC (August vs September)",
     xlab = "FFMC")
lines(density(sept_data$FFMC), col = 'pink', lwd = 2)
legend("topright", legend = c("August", "September"),
       col = c("green", "pink"), lwd = 2)

# H0: Variance of FFMC in August is less than or equal to September
# H1: Variance of FFMC in August is greater than in September

mood.test(aug_data$FFMC, sept_data$FFMC, alternative = 'greater')

# p-value = 1.173e-05 < 0.05 ⇒ Reject H0.
# Conclusion: FFMC shows significantly greater variability in August compared to September.
# Implication: Fine fuel moisture fluctuates more in August, creating unpredictable ignition conditions.



#---------------- ISI vs Month (August vs September) ----------------

plot(density(aug_data$ISI), ylim = c(0, 0.2),
     col = 'green', lwd = 2,
     main = "Density Comparison: ISI (August vs September)",
     xlab = "ISI")
lines(density(sept_data$ISI), col = 'pink', lwd = 2)
legend("topright", legend = c("August", "September"),
       col = c("green", "pink"), lwd = 2)

# H0: Variance of ISI in August is less than or equal to September
# H1: Variance of ISI in August is greater than in September

mood.test(aug_data$ISI, sept_data$ISI, alternative = 'greater')

# p-value = 0.001421 < 0.05 ⇒ Reject H0.
# Conclusion: ISI exhibits significantly higher variability in August than September.
# Implication: Fire spread potential is more erratic in August due to fluctuating wind and moisture interactions.



#===================================================================
#                     DISTRIBUTION TESTS
#===================================================================

#---------------- FFMC vs Fire Occurrence ---------------------------

plot(density(df[df$fire == 'Yes', 'FFMC']), ylim = c(0, 0.25),
     col = 'green', lwd = 2,
     main = "Density Comparison: FFMC (Fire vs No Fire)",
     xlab = "FFMC")
lines(density(df[df$fire == 'No', 'FFMC']), col = 'pink', lwd = 2)
legend("topleft", legend = c("Fire", "No Fire"),
       col = c("green", "pink"), lwd = 2)

# H0: The distribution of FFMC is identical for days with and without burned area
# H1: The distributions differ

ks.test(df[df$fire == 'Yes', 'FFMC'], df[df$fire == 'No', 'FFMC'])

# p-value = 0.7826 > 0.05 ⇒ Fail to reject H0.
# Conclusion: FFMC distributions are statistically indistinguishable between fire and no-fire days.
# Implication: Fine fuel moisture alone does not differentiate burned from non-burned days in this dataset.



#---------------- RH vs Fire Occurrence -----------------------------

plot(density(df[df$fire == 'Yes', 'RH']),
     col = 'green', lwd = 2,
     main = "Density Comparison: RH (Fire vs No Fire)",
     xlab = "RH")
lines(density(df[df$fire == 'No', 'RH']), col = 'pink', lwd = 2)
legend("topright", legend = c("Fire", "No Fire"),
       col = c("green", "pink"), lwd = 2)

# H0: The distribution of relative humidity is identical for days with and without burned area
# H1: The distributions differ

ks.test(df[df$fire == 'Yes', 'RH'], df[df$fire == 'No', 'RH'])

# p-value = 0.788 > 0.05 ⇒ Fail to reject H0.
# Conclusion: Relative humidity distributions show no significant difference between fire conditions.
# Implication: RH distribution is consistent across fire events; other moisture indices may be more discriminative.



#---------------- ISI vs Fire Occurrence ----------------------------

plot(density(df[df$fire == 'Yes', 'ISI']),
     col = 'green', lwd = 2,
     main = "Density Comparison: ISI (Fire vs No Fire)",
     xlab = "ISI")
lines(density(df[df$fire == 'No', 'ISI']), col = 'pink', lwd = 2)
legend("topright", legend = c("Fire", "No Fire"),
       col = c("green", "pink"), lwd = 2)

# H0: The distribution of ISI is identical for days with and without burned area
# H1: The distributions differ

ks.test(df[df$fire == 'Yes', 'ISI'], df[df$fire == 'No', 'ISI'])

# p-value = 0.7327 > 0.05 ⇒ Fail to reject H0.
# Conclusion: ISI distributions are statistically similar between fire and no-fire conditions.
# Implication: Initial spread index distribution does not vary with fire occurrence in this forest region.



#---------------- Wind vs Fire Occurrence ---------------------------

plot(density(df[df$fire == 'Yes', 'wind']),
     col = 'green', lwd = 2,
     main = "Density Comparison: Wind (Fire vs No Fire)",
     xlab = "Wind")
lines(density(df[df$fire == 'No', 'wind']), col = 'pink', lwd = 2)
legend("topright", legend = c("Fire", "No Fire"),
       col = c("green", "pink"), lwd = 2)

# H0: The distribution of wind speed is identical for days with and without burned area
# H1: The distributions differ

ks.test(df[df$fire == 'Yes', 'wind'], df[df$fire == 'No', 'wind'])

# p-value = 0.9269 > 0.05 ⇒ Fail to reject H0.
# Conclusion: Wind speed distributions are statistically equivalent across fire conditions.
# Implication: Wind distribution alone does not distinguish days with burned area from those without.



#---------------- RH vs Month (August vs September) -----------------

plot(density(aug_data$RH), ylim = c(0, 0.05),
     col = 'green', lwd = 2,
     main = "Density Comparison: RH (August vs September)",
     xlab = "Relative Humidity (%)")
lines(density(sept_data$RH), col = 'pink', lwd = 2)
legend("topright", legend = c("August", "September"),
       col = c("green", "pink"), lwd = 2)

# H0: The distribution of relative humidity is identical in August and September
# H1: The distributions differ

ks.test(aug_data$RH, sept_data$RH)

# p-value = 0.3249 > 0.05 ⇒ Fail to reject H0.
# Conclusion: Relative humidity distributions are statistically similar between the two months.
# Implication: Despite different drought conditions, RH patterns remain stable across August and September.



#---------------- Temperature vs Month (August vs September) --------

plot(density(aug_data$temp),
     col = 'green', lwd = 2, ylim = c(0, 0.3),
     main = "Distribution Comparison: Temperature (Aug vs Sept)",
     xlab = "Temperature")
lines(density(sept_data$temp), col = 'pink', lwd = 2)
legend("topright", legend = c("August", "September"),
       col = c("green", "pink"), lwd = 2)

# H0: Temperature distributions are identical in August and September
# H1: August temperatures are stochastically lower than September

ks.test(aug_data$temp, sept_data$temp, alternative = 'less')

# p-value = 0.0001669 < 0.05 ⇒ Reject H0.
# Conclusion: The temperature distribution in August is stochastically lower than September.
# Implication: Despite higher median in August, the overall temperature distribution differs, with September showing a rightward shift in some ranges.



#===================================================================
#                     CORRELATION TESTS
#===================================================================

data_z <- df[df$area != 0, ]
X <- c("FFMC", "DMC", "DC", "ISI", "temp", "RH", "wind", "rain")


# H0: No monotonic relationship exists between weather variables and burned area extent
# H1: A significant monotonic relationship exists

for (x in X) {
  print(paste(x, "vs area"))
  r <- cor.test(data_z[[x]], data_z$area, method = "spearman")
  print(r)
  if (r$p.value < 0.05) {
    print("Reject H0: Significant correlation.")
  } else {
    print("Fail to reject H0: No significant correlation.")
  }
}



# Conclusion: All p-values > 0.05 ⇒ Fail to reject H0.
# Implication: When fires do occur, individual weather indices show weak monotonic relationships with burned area size, suggesting complex multivariate interactions.


#===================================================================
#                   GROUP COMPARISON TESTS
#===================================================================

ggplot(df, aes(x = day, y = area, fill = day)) +
  geom_boxplot(outlier.color = "red", alpha = 0.7) +
  labs(title = "Burned Area vs Day of the Week",
       subtitle = "Kruskal-Wallis test on median burned area",
       x = "Day of the Week", y = "Burned Area (ha)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")


# H0: Median burned area is equal across all days of the week
# H1: At least one weekday has a different median burned area

kruskal.test(area ~ day, data = df)

# p-value > 0.05 ⇒ Fail to reject H0.
# Conclusion: No significant difference in median burned area across weekdays.
# Implication: Forest fire severity shows no weekly temporal pattern; human activity schedules do not influence fire outcomes.


