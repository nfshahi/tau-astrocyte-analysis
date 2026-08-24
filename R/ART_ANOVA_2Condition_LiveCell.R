# Load required packages
library(readxl)
library(dplyr)
library(ARTool)
library(ggplot2)
library(writexl)
library(stringr)
library(wesanderson)

# Set variables to the appropriate data types
raw_data$Condition <- as.factor(raw_data$Condition)
raw_data$Time <- as.factor(raw_data$Time)
raw_data$Tau_Puncta <- as.numeric(raw_data$Tau_Puncta)

# Fit the ART model with Condition, Time, and their interaction.
# Well_ID is included as a random intercept to account for repeated
# measurements from the same well.
art_model <- art(Tau_Puncta ~ Condition * Time + (1|Well_ID), data = raw_data)

# Extract the ANOVA results
anova_results <- anova(art_model)
print(anova_results)

# Save the ANOVA results
anova_df <- as.data.frame(anova_results)
write_xlsx(anova_df, "ART_ANOVA_mixed_results.xlsx")

# Calculate post hoc contrasts for the Condition × Time interaction
posthoc_df <- art.con(art_model, "Condition:Time", adjust = "none") %>%
  as.data.frame()

# Keep contrasts comparing Control and Heparin conditions
posthoc_filtered <- posthoc_df %>%
  filter(grepl("Control", contrast) & grepl("Heparin", contrast))

# Extract the two time points from each contrast and retain comparisons
# between conditions at the same time point
posthoc_filtered <- posthoc_filtered %>%
  mutate(
    times = str_extract_all(contrast, ",\\s*(\\d+)"),
    time1 = as.numeric(sapply(times, function(x) gsub(",", "", x[1]))),
    time2 = as.numeric(sapply(times, function(x) gsub(",", "", x[2])))
  ) %>%
  filter(time1 == time2) %>%
  select(-times)

# Apply FDR and Bonferroni corrections to the post hoc p-values
posthoc_filtered$FDR <- p.adjust(posthoc_filtered$p.value, method = "fdr")
posthoc_filtered$bonferroni <- p.adjust(
  posthoc_filtered$p.value,
  method = "bonferroni"
)

# View the filtered post hoc results
print(posthoc_filtered)

# Save the post hoc results
write_xlsx(posthoc_filtered, "Tau_posthoc_results.xlsx")

# Calculate mean and SEM for each condition and time point
summary_data <- raw_data %>%
  group_by(Condition, Time) %>%
  summarise(
    mean_tau = mean(Tau_Puncta, na.rm = TRUE),
    sem_tau = sd(Tau_Puncta, na.rm = TRUE) / sqrt(n())
  )

# Plot mean tau puncta over time with SEM
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

# Plot the same data using a classic theme without grid lines
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

# Convert time to numeric for plotting
summary_data$Time <- as.numeric(as.character(summary_data$Time))

summary_data <- raw_data %>%
  mutate(Time_num = as.numeric(as.character(Time))) %>%
  group_by(Condition, Time_num) %>%
  summarise(
    mean_tau = mean(Tau_Puncta, na.rm = TRUE),
    sem_tau = sd(Tau_Puncta, na.rm = TRUE) / sqrt(n())
  )

# Set the colors used for the plot
palette_colors <- wes_palette("Royal1")

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
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = mean_tau - sem_tau,
      ymax = mean_tau + sem_tau
    ),
    width = 0.3
  ) +
  scale_color_manual(
    values = c(
      "Control" = palette_colors[1],
      "Heparin" = palette_colors[1]
    )
  ) +
  scale_shape_manual(
    values = c(
      "Control" = 16,
      "Heparin" = 16
    )
  ) +
  scale_x_continuous(
    limits = c(0, 25),
    breaks = seq(0, 25, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 5000),
    breaks = seq(0, 5000, by = 1000),
    expand = c(0, 0)
  ) +
  labs(
    x = "Time (Hours)",
    y = "Number of Tau Puncta",
    title = "Tau Puncta over Time by Condition"
  ) +
  theme_classic()

# Assign colors to the two conditions
royal_colors <- wes_palette("Royal1")

custom_colors <- c(
  "Control" = royal_colors[1],
  "Heparin" = royal_colors[1]
)

# Recalculate summary statistics for the final plot
summary_data <- raw_data %>%
  mutate(Time_num = as.numeric(as.character(Time))) %>%
  group_by(Condition, Time_num) %>%
  summarise(
    mean_tau = mean(Tau_Puncta, na.rm = TRUE),
    sem_tau = sd(Tau_Puncta, na.rm = TRUE) / sqrt(n())
  )

# Plot mean tau puncta with SEM ribbons
ggplot(
  summary_data,
  aes(
    x = Time_num,
    y = mean_tau,
    color = Condition,
    fill = Condition,
    group = Condition,
    shape = Condition
  )
) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_ribbon(
    aes(
      ymin = mean_tau - sem_tau,
      ymax = mean_tau + sem_tau
    ),
    alpha = 0.25,
    color = NA
  ) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_shape_manual(
    values = c(
      "Control" = 16,
      "Heparin" = 2
    )
  ) +
  scale_x_continuous(
    limits = c(0, 25),
    breaks = seq(0, 25, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 4000),
    breaks = seq(0, 4000, by = 1000),
    expand = c(0, 0)
  ) +
  labs(
    x = "Time (Hours)",
    y = "Number of Tau Puncta",
    title = "Tau Puncta over Time by Condition"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 18),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank()
  )
