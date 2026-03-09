library(readr)
library(readxl)
library(lubridate)
library(tidyverse)
library(dplyr)
library(stringi) # or library(stringr)

dir <- "/Users/johanna/Uni/ws_25_26/ecohydrology/ecohydro_project/data/Water stable isotope analysis"

fname_stand <- "+ 2024 cep in-house  standards.xlsx"
fname_template <- "2026 Ecohydrology Template water stable isotopes.csv"
fname_data <- "20260309 Water stable isotope data.csv"

data <- read_csv(file.path(dir, fname_data))
standards <- read_excel(file.path(dir, fname_stand), skip = 1)
template <- read_csv(file.path(dir, fname_template))

colnames(data) <- gsub(" ", "_", colnames(data))

# Correct time stamp
data$Time_Code <- ymd_hms(data$Time_Code)

cols <- data %>%
  select(where(is.character)) %>%
  colnames()

data$Identifier_2[which(data$Identifier_2==data$Identifier_2[55])] <- "spuelen"

# Remove blank spaces (using stringi for safety)
for(i in cols) {
  data[[i]] <- stringi::stri_replace_all_fixed(data[[i]], " ", "", vectorize_all = FALSE)
}

# Fix entries
data[data$Identifier_1=="SStandardheavy", "Identifier_1"] <- "Standardheavy"
data[data$Identifier_1=="Standardheavy2", "Identifier_1"] <- "Superheavy"

data <- data %>%
  mutate(
    d.18_16.Mean = as.numeric(`d(18_16)Mean`),
    d.D_H.Mean = as.numeric(`d(D_H)Mean`)
  )

# Only keep the last 3
ids <- unique(data$Identifier_1)[-c(1:3)] #removing standards
# Filter the data and keep only the last 3 injections per ID
data_short <- data %>%
  filter(Identifier_1 %in% ids) %>%
  group_by(Identifier_1) %>%
  arrange(Identifier_1, desc(Inj_Nr)) %>%  # Sort by ID and descending injection number
  slice_head(n = 3) %>%                     # Keep only the top 3 (last 3 injections)
  ungroup()

# Summarize
split_date <- ymd_hms("2026-03-04 00:00:00 UTC")
data_short_1 <- data_short %>%
  filter(Time_Code <= split_date) %>%
  mutate(
    d.18_16.Mean = as.numeric(d.18_16.Mean),
    d.D_H.Mean = as.numeric(d.D_H.Mean)
  ) %>%
  group_by(Identifier_1) %>%
  summarise(
    d.18_16.Mean.Mean = mean(d.18_16.Mean, na.rm = TRUE),
    `d(18_16)Mean.Sd` = sd(d.18_16.Mean, na.rm = TRUE),
    `d(D_H)Mean.Mean` = mean(d.D_H.Mean, na.rm = TRUE),
    `d(D_H)Mean.Sd` = sd(d.D_H.Mean, na.rm = TRUE),
    .groups = "drop"
  )
data_short_2 <- data_short %>%
  filter(Time_Code > split_date) %>%
  mutate(
    d.18_16.Mean = as.numeric(d.18_16.Mean),
    d.D_H.Mean = as.numeric(d.D_H.Mean)
  ) %>%
  group_by(Identifier_1) %>%
  summarise(
    d.18_16.Mean.Mean = mean(d.18_16.Mean, na.rm = TRUE),
    `d(18_16)Mean.Sd` = sd(d.18_16.Mean, na.rm = TRUE),
    `d(D_H)Mean.Mean` = mean(d.D_H.Mean, na.rm = TRUE),
    `d(D_H)Mean.Sd` = sd(d.D_H.Mean, na.rm = TRUE),
    .groups = "drop"
  )
data_short <- data_short %>%
  mutate(
    d.18_16.Mean = as.numeric(d.18_16.Mean),
    d.D_H.Mean = as.numeric(d.D_H.Mean)
  ) %>%
  group_by(Identifier_1) %>%
  summarise(
    d.18_16.Mean.Mean = mean(d.18_16.Mean, na.rm = TRUE),
    `d(18_16)Mean.Sd` = sd(d.18_16.Mean, na.rm = TRUE),
    `d(D_H)Mean.Mean` = mean(d.D_H.Mean, na.rm = TRUE),
    `d(D_H)Mean.Sd` = sd(d.D_H.Mean, na.rm = TRUE),
    .groups = "drop"
  )

ids_standards <- unique(data$Identifier_1)[1:3]
data_standards <- data %>%
  filter(Identifier_1 %in% ids_standards) %>%
  group_by(Port) %>%
  arrange(Port, desc(Inj_Nr)) %>%  # Sort by ID and descending injection number
  slice_head(n = 3) %>%                     # Keep only the top 3 (last 3 injections)
  ungroup()

# Identify jumps
test1 <- data_standards %>% filter(Identifier_1 == "Standardmedium")
plot(test1$Time_Code, test1$d.D_H.Mean)
test2 <- data_standards %>% filter(Identifier_1 == "Standardheavy")
plot(test2$Time_Code, test2$d.D_H.Mean)
test3 <- data_standards %>% filter(Identifier_1 == "Superheavy")
plot(test3$Time_Code, test3$d.D_H.Mean)

# Fit model for calibration
standards <- data.frame(standard= c("Standardmedium", "Standardheavy", "Superheavy"),
                        d.18_16= c(-9.05, -4.91, -0.18),
                        d.D_H= c(-63.63, -10.18, 54.05))

# kalibriert alle laufenden Standards mit dem ersten Standardblock

# berechnet dann die Differenz der Werte zu den Soll Werten (18O Dif und 2 H Dif). Dann könnt ihr diese Differenz gegen die Injektionnummer oder besser Zeit bei euch auftragen und dann könnt ihr sehen, wie die Werte wegdriften. Wenn das so aussieht, wie bei 2H in dem Beispiel, dann könnt ihr für x die Zeit/Injektionsnummer einsetzen und diesen Wert dann auf euren standardisiertes Ergebnis draufrechnen (Tab Summary auf der rechten Seite bei 18Ocor und 2Hcor).

Ich würden den letzten Standardblock hinten dran auch mal wie normale Proben behandeln, dann könnt ihr schauen, wie dieser nach der Driftkorrektur aussieht.




# ------ Drift correction

# First attempt

first_block_end_date <- ymd_hms("2026-03-02 23:59:00 UTC")
last_block_start_date <- ymd_hms("2026-03-08 00:00:00 UTC")

first_calib <- data_standards %>% filter(Time_Code <= first_block_end_date)
second_calib <- data_standards %>% filter(Time_Code >= last_block_start_date)

test_data_1 <- full_join(first_calib, standards, by= join_by(Identifier_1== standard))
test_data_2 <- full_join(second_calib, standards, by= join_by(Identifier_1== standard))
test_data <- full_join(data_standards, standards, by= join_by(Identifier_1== standard))

# First Block
lm_d.18_16_1 <- lm(d.18_16 ~ d.18_16.Mean, test_data_1)
data_short_1$d.18_16.corr <- unname(
  predict(lm_d.18_16_1, newdata = data.frame("d.18_16.Mean"= data_short_1$d.18_16.Mean.Mean))
)

plot(d.18_16~d.18_16.Mean, test_data_1, col=as.factor(Identifier_1), pch=16,
     xlab="measured Delta18", ylab="Standard Delta18")
curve(lm_d.18_16_1$coefficients[2]*x+ lm_d.18_16_1$coefficients[1], add=T)

# Second Block
lm_d.18_16_2 <- lm(d.18_16 ~ d.18_16.Mean, test_data_2)
data_short_2$d.18_16.corr <- unname(
  predict(lm_d.18_16_2, newdata = data.frame("d.18_16.Mean"= data_short_2$d.18_16.Mean.Mean))
)

plot(d.18_16~d.18_16.Mean, test_data_2, col=as.factor(Identifier_1), pch=16,
     xlab="measured Delta18 (2nd Block)", ylab="Standard Delta18 (2nd Block)")
curve(lm_d.18_16_2$coefficients[2]*x+ lm_d.18_16_2$coefficients[1], add=T)

# -- Second attempt
lm_d.18_16 <- lm(d.18_16 ~ d.18_16.Mean + Time_Code, test_data)
data_short$d.18_16.corr <- unname(
  predict(lm_d.18_16, newdata = data.frame(
    "d.18_16.Mean"= data_short$d.18_16.Mean.Mean,

  ))
)





# Plot Standards
plot(data_short_1$Time_Code, data_short_1$d.18_16.Mean.Mean)