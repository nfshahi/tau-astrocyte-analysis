# Load required packages
install.packages(c("readxl", "dplyr", "ARTool", "ggplot2", "writexl"))
library(readxl)
library(dplyr)
library(ARTool)
library(ggplot2)
library(writexl)
library(stringr)


# Load the Excel data
file_path <- "path/to/your/data.xlsx"
raw_data <- read_excel(file_path, sheet = "Structured")


# Set variable types
raw_data$Condition <- as.factor(raw_data$Condition)
raw_data$Time <- as.factor(raw_data$Time)
raw_data$Tau_Puncta <- as.numeric(raw_data$Tau_Puncta)


# Fit the ART model with Well_ID included as a random effect
art_model <- art(
  Tau_Puncta ~ Condition * Time + (1 | Well_ID),
  data = raw_data
)


# Save ANOVA results
anova_results <- anova(art_model)
write_xlsx(anova_results, "ANOVA_results.xlsx")


# Extract post-hoc contrasts for the Condition × Time interaction
posthoc_df <- art.con(
  art_model,
  "Condition:Time",
  adjust = "none"
) %>%
  as.data.frame()


# Conditions included in the analysis
target_conditions <- c(
  "Tau_Veh",
  "Tau_LRP1",
  "S/T P-Tau_Veh",
  "S/T P-Tau_LRP1"
)


# Keep pairwise comparisons between conditions at the same time point
posthoc_filtered <- posthoc_df %>%
  mutate(
    parts = str_match(
      contrast,
      "\\(?\\s*([A-Za-z0-9 /_\\-]+)\\s*,\\s*(\\d+)\\s*\\)?\\s*-\\s*\\(?\\s*([A-Za-z0-9 /_\\-]+)\\s*,\\s*(\\d+)\\s*\\)?"
    ),
    cond1 = str_trim(parts[, 2]),
    time1 = as.numeric(parts[, 3]),
    cond2 = str_trim(parts[, 4]),
    time2 = as.numeric(parts[, 5])
  ) %>%
  filter(
    !is.na(time1) & !is.na(time2),
    cond1 %in% target_conditions,
    cond2 %in% target_conditions,
    cond1 != cond2,
    time1 == time2
  ) %>%
  mutate(
    comparison = paste0(cond1, " vs ", cond2),
    FDR = p.adjust(p.value, method = "fdr"),
    bonferroni = p.adjust(p.value, method = "bonferroni")
  ) %>%
  select(
    time1,
    comparison,
    contrast,
    p.value,
    FDR,
    bonferroni
  ) %>%
  arrange(time1, comparison)


# Save filtered post-hoc results
write_xlsx(
  posthoc_filtered,
  "Posthoc_4_conditions_same_timepoint.xlsx"
)


# Print the first few post-hoc results
print(head(posthoc_filtered))


# Calculate mean and SEM for each condition and time point
summary_data <- raw_data %>%
  group_by(Condition, Time) %>%
  summarise(
    mean_tau = mean(Tau_Puncta, na.rm = TRUE),
    sem_tau = sd(Tau_Puncta, na.rm = TRUE) / sqrt(n())
  )


# Plot mean tau puncta ± SEM
ggplot(
  summary_data,
  aes(
    x = as.numeric(Time),
    y = mean_tau,
    color = Condition,
    group = Condition
  )
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = mean_tau - sem_tau,
      ymax = mean_tau + sem_tau
    ),
    width = 0.3
  ) +
  labs(
    x = "Time (hours)",
    y = "Mean Tau Puncta ± SEM",
    title = "Tau Puncta over Time by Condition"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.title = element_blank()
  )


# Plot without grid lines
ggplot(
  summary_data,
  aes(
    x = as.numeric(Time),
    y = mean_tau,
    color = Condition,
    group = Condition
  )
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = mean_tau - sem_tau,
      ymax = mean_tau + sem_tau
    ),
    width = 0.3
  ) +
  labs(
    x = "Time (hours)",
    y = "Mean Tau Puncta ± SEM",
    title = "Tau Puncta over Time by Condition"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    legend.title = element_blank()
  )


# Prepare data for the final plot
summary_data <- raw_data %>%
  mutate(Time_num = as.numeric(as.character(Time))) %>%
  group_by(Condition, Time_num) %>%
  summarise(
    mean_tau = mean(Tau_Puncta, na.rm = TRUE),
    sem_tau = sd(Tau_Puncta, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


# Get colors from the Royal1 palette
royal_colors <- wes_palette("Royal1", n = 5, type = "continuous")


# Assign colors to the conditions
custom_colors <- c(
  "Tau_Veh" = royal_colors[1],
  "Tau_LRP1" = royal_colors[1],
  "S/T P-Tau_Veh" = royal_colors[2],
  "S/T P-Tau_LRP1" = royal_colors[2],
  "No Tau" = royal_colors[4]
)


# Assign plotting symbols to the conditions
custom_shapes <- c(
  "Tau_Veh" = 16,
  "Tau_LRP1" = 0,
  "S/T P-Tau_Veh" = 16,
  "S/T P-Tau_LRP1" = 0,
  "No Tau" = 16
)


# Plot with condition-specific colors and symbols
ggplot(
  summary_data,
  aes(
    x = Time_num,
    y = mean_tau,
    color = Condition,
    group = Condition,
    shape = Condition
  )
) +
  geom_line(size = 0.5) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = mean_tau - sem_tau,
      ymax = mean_tau + sem_tau
    ),
    width = 0.3
  ) +
  scale_color_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  scale_x_continuous(
    limits = c(0, 25),
    breaks = seq(0, 25, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 15000),
    breaks = seq(0, 15000, by = 3000),
    expand = c(0, 0)
  ) +
  labs(
    x = "Time (Hours)",
    y = "Number of Tau Puncta",
    title = "Tau Puncta over Time by Condition"
  ) +
  theme_classic()


# Plot with SEM shown as a shaded area
ggplot(
  summary_data,
  aes(
    x = Time_num,
    y = mean_tau,
    color = Condition,
    group = Condition,
    shape = Condition
  )
) +
  geom_line(size = 0.5) +
  geom_point(size = 3) +
  geom_ribbon(
    aes(
      ymin = mean_tau - sem_tau,
      ymax = mean_tau + sem_tau,
      fill = Condition
    ),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_shape_manual(values = custom_shapes) +
  scale_x_continuous(
    limits = c(0, 25),
    breaks = seq(0, 25, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 15000),
    breaks = seq(0, 15000, by = 3000),
    expand = c(0, 0)
  ) +
  labs(
    x = "Time (Hours)",
    y = "Tau Puncta",
    title = "Tau Puncta over Time by Condition"
  ) +
  theme_classic()
