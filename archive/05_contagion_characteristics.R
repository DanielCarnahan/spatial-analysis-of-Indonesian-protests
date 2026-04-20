# Contagion Characteristics Analysis
# What makes protests more likely to trigger other protests?
# Author: Daniel Carnahan
# Date: 2025-10-27

library(dplyr)
library(ggplot2)
library(tidyr)

# Load data
cat("Loading protest data...\n")
protests <- readRDS("protests_prepared.rds")
protest_stpp <- readRDS("protest_stpp.rds")

cat("Loaded", nrow(protests), "events\n\n")

# ====================================================================
# 1. IDENTIFY AVAILABLE CHARACTERISTICS (MARKS)
# ====================================================================

cat("=== AVAILABLE PROTEST CHARACTERISTICS ===\n\n")

# Key characteristics we can analyze
characteristics <- list(
  "Violence Level" = c("is_violent", "is_peaceful"),
  "Fatalities" = "fatalities",
  "State Intervention" = "state_intervention",
  "Event Type" = "event_type",
  "Sub-Event Type" = "sub_event_type",
  "Actor 1" = "actor1",
  "Actor 2" = "actor2",
  "Interaction Type" = "interaction",
  "Civilian Targeting" = "civilian_targeting",
  "Region" = "admin1"
)

cat("Available characteristics for contagion analysis:\n")
for(i in seq_along(characteristics)) {
  cat(sprintf("  %d. %s\n", i, names(characteristics)[i]))
}

cat("\n")

# ====================================================================
# 2. HYPOTHESIS GENERATION
# ====================================================================

cat("=== CONTAGION HYPOTHESES ===\n\n")

hypotheses <- data.frame(
  H = paste0("H", 1:10),
  Characteristic = c(
    "Violence",
    "Fatalities",
    "State repression",
    "Civilian targeting",
    "Event type",
    "Actor type",
    "Geographic proximity",
    "Temporal recency",
    "Urban location",
    "Prior local activity"
  ),
  Hypothesis = c(
    "Violent protests trigger MORE subsequent protests (anger/mobilization)",
    "Protests with fatalities trigger MORE protests (martyrdom effect)",
    "State intervention triggers MORE protests (backlash effect)",
    "Civilian-targeting events trigger MORE protests (outrage)",
    "Riots are MORE contagious than peaceful protests",
    "Certain actors (students, workers) are MORE mobilizing",
    "Nearby protests are MORE likely to be triggered (diffusion)",
    "Recent protests are MORE likely to trigger (memory decay)",
    "Urban protests are MORE contagious (visibility/networks)",
    "Areas with history trigger MORE easily (activation)"
  ),
  Alternative = c(
    "Violent protests trigger FEWER (fear/deterrence)",
    "Fatalities SUPPRESS protests (fear)",
    "Repression DETERS future protests",
    "Targeting causes fear, REDUCES protests",
    "Peaceful protests MORE contagious (safety)",
    "Elite actors LESS mobilizing (distance)",
    "Distance does not matter (media diffusion)",
    "Memory is long-lasting (no decay)",
    "Rural protests equally contagious (local grievances)",
    "No path dependence"
  ),
  stringsAsFactors = FALSE
)

print(hypotheses)
cat("\n")

# ====================================================================
# 3. EXPLORATORY: TEMPORAL PATTERNS BY CHARACTERISTIC
# ====================================================================

cat("=== EXPLORATORY: Do different protest types show different temporal patterns? ===\n\n")

# Function to measure "triggering potential"
# Simple approach: After an event, count how many events occur nearby in next N days

measure_triggering <- function(data, focal_condition, spatial_thresh = 0.5, temporal_window = 30) {

  # Events matching focal condition
  focal_events <- data %>% filter(eval(parse(text = focal_condition)))

  cat(sprintf("  Focal events (%s): %d\n", focal_condition, nrow(focal_events)))

  if(nrow(focal_events) == 0) return(NULL)

  # For each focal event, count subsequent events nearby
  triggered_counts <- numeric(nrow(focal_events))

  # Sample if too many (for speed)
  if(nrow(focal_events) > 500) {
    sample_idx <- sample(1:nrow(focal_events), 500)
    focal_events <- focal_events[sample_idx, ]
    triggered_counts <- numeric(500)
  }

  for(i in 1:nrow(focal_events)) {
    focal_x <- focal_events$longitude[i]
    focal_y <- focal_events$latitude[i]
    focal_t <- focal_events$days_since_start[i]

    # Find subsequent events within spatio-temporal window
    subsequent <- data %>%
      filter(
        days_since_start > focal_t,
        days_since_start <= focal_t + temporal_window,
        sqrt((longitude - focal_x)^2 + (latitude - focal_y)^2) < spatial_thresh
      )

    triggered_counts[i] <- nrow(subsequent)
  }

  return(data.frame(
    condition = focal_condition,
    mean_triggered = mean(triggered_counts),
    median_triggered = median(triggered_counts),
    sd_triggered = sd(triggered_counts),
    n_focal = nrow(focal_events)
  ))
}

cat("\n--- Measuring Triggering Potential ---\n")
cat("(For each event type, count how many protests occur within 55km in next 30 days)\n\n")

# Test different conditions
conditions_to_test <- c(
  "is_violent == TRUE",
  "is_peaceful == TRUE",
  "fatalities > 0",
  "fatalities == 0",
  "state_intervention == TRUE",
  "state_intervention == FALSE",
  "event_type == 'Riots'",
  "event_type == 'Protests'"
)

results_list <- list()

for(cond in conditions_to_test) {
  cat(sprintf("\nTesting: %s\n", cond))
  result <- measure_triggering(protests, cond, spatial_thresh = 0.5, temporal_window = 30)
  if(!is.null(result)) {
    results_list[[length(results_list) + 1]] <- result
  }
}

triggering_comparison <- bind_rows(results_list)

cat("\n=== TRIGGERING POTENTIAL COMPARISON ===\n")
print(triggering_comparison)
cat("\n")

# Save results
write.csv(triggering_comparison, "triggering_comparison.csv", row.names = FALSE)
cat("Saved: triggering_comparison.csv\n\n")

# ====================================================================
# 4. VISUALIZATION: Triggering by Characteristic
# ====================================================================

cat("=== Creating visualizations ===\n")

# Plot 1: Triggering potential by violence
png("plots/09_triggering_by_violence.png", width = 1000, height = 600, res = 120)

violence_comparison <- triggering_comparison %>%
  filter(grepl("violent|peaceful", condition, ignore.case = TRUE))

if(nrow(violence_comparison) > 0) {
  violence_comparison <- violence_comparison %>%
    mutate(
      type = case_when(
        grepl("violent == TRUE", condition) ~ "Violent",
        grepl("peaceful == TRUE", condition) ~ "Peaceful"
      )
    )

  ggplot(violence_comparison, aes(x = type, y = mean_triggered, fill = type)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = mean_triggered - sd_triggered/sqrt(n_focal),
                      ymax = mean_triggered + sd_triggered/sqrt(n_focal)),
                  width = 0.2) +
    scale_fill_manual(values = c("Violent" = "#d62728", "Peaceful" = "#2ca02c")) +
    labs(
      title = "Triggering Potential: Violent vs Peaceful Protests",
      subtitle = "Mean number of subsequent protests within 55km in next 30 days",
      x = "Protest Type",
      y = "Mean Triggered Events",
      fill = "Type"
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          text = element_text(size = 12))
}

dev.off()
cat("Saved: plots/09_triggering_by_violence.png\n")

# Plot 2: Triggering by fatalities
png("plots/10_triggering_by_fatalities.png", width = 1000, height = 600, res = 120)

fatality_comparison <- triggering_comparison %>%
  filter(grepl("fatalities", condition))

if(nrow(fatality_comparison) > 0) {
  fatality_comparison <- fatality_comparison %>%
    mutate(
      type = case_when(
        grepl("fatalities > 0", condition) ~ "With Fatalities",
        grepl("fatalities == 0", condition) ~ "No Fatalities"
      )
    )

  ggplot(fatality_comparison, aes(x = type, y = mean_triggered, fill = type)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = mean_triggered - sd_triggered/sqrt(n_focal),
                      ymax = mean_triggered + sd_triggered/sqrt(n_focal)),
                  width = 0.2) +
    scale_fill_manual(values = c("With Fatalities" = "#ff7f0e",
                                  "No Fatalities" = "#1f77b4")) +
    labs(
      title = "Triggering Potential: Protests With vs Without Fatalities",
      subtitle = "Mean number of subsequent protests within 55km in next 30 days",
      x = "Fatality Status",
      y = "Mean Triggered Events",
      fill = "Type"
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          text = element_text(size = 12))
}

dev.off()
cat("Saved: plots/10_triggering_by_fatalities.png\n")

# Plot 3: Triggering by state intervention
png("plots/11_triggering_by_state_intervention.png", width = 1000, height = 600, res = 120)

state_comparison <- triggering_comparison %>%
  filter(grepl("state_intervention", condition))

if(nrow(state_comparison) > 0) {
  state_comparison <- state_comparison %>%
    mutate(
      type = case_when(
        grepl("TRUE", condition) ~ "With State Intervention",
        grepl("FALSE", condition) ~ "No State Intervention"
      )
    )

  ggplot(state_comparison, aes(x = type, y = mean_triggered, fill = type)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = mean_triggered - sd_triggered/sqrt(n_focal),
                      ymax = mean_triggered + sd_triggered/sqrt(n_focal)),
                  width = 0.2) +
    scale_fill_manual(values = c("With State Intervention" = "#9467bd",
                                  "No State Intervention" = "#8c564b")) +
    labs(
      title = "Triggering Potential: State Intervention Effect",
      subtitle = "Mean number of subsequent protests within 55km in next 30 days",
      x = "State Response",
      y = "Mean Triggered Events",
      fill = "Type"
    ) +
    theme_minimal() +
    theme(legend.position = "none",
          text = element_text(size = 12))
}

dev.off()
cat("Saved: plots/11_triggering_by_state_intervention.png\n\n")

# ====================================================================
# 5. ACTOR ANALYSIS
# ====================================================================

cat("=== ACTOR CONTAGIOUSNESS ANALYSIS ===\n\n")

# Which actors are most associated with subsequent protests?
actor_counts <- protests %>%
  count(actor1, sort = TRUE) %>%
  head(15)

cat("Top 15 actors by frequency:\n")
print(actor_counts)
cat("\n")

# For top actors, measure triggering
top_actors <- actor_counts$actor1[1:5]

actor_triggering <- data.frame()

for(actor in top_actors) {
  condition <- sprintf("actor1 == '%s'", actor)
  result <- measure_triggering(protests, condition,
                               spatial_thresh = 0.5,
                               temporal_window = 30)
  if(!is.null(result)) {
    result$actor <- actor
    actor_triggering <- bind_rows(actor_triggering, result)
  }
}

cat("\n=== TRIGGERING BY ACTOR ===\n")
print(actor_triggering %>% select(actor, mean_triggered, n_focal))
cat("\n")

# Save
write.csv(actor_triggering, "actor_triggering.csv", row.names = FALSE)
cat("Saved: actor_triggering.csv\n\n")

# ====================================================================
# 6. REGIONAL CONTAGIOUSNESS
# ====================================================================

cat("=== REGIONAL CONTAGIOUSNESS ===\n\n")

# Which provinces show strongest local contagion?
top_provinces <- protests %>%
  count(admin1, sort = TRUE) %>%
  head(10) %>%
  pull(admin1)

regional_triggering <- data.frame()

for(province in top_provinces) {
  condition <- sprintf("admin1 == '%s'", province)
  result <- measure_triggering(protests, condition,
                               spatial_thresh = 0.5,
                               temporal_window = 30)
  if(!is.null(result)) {
    result$province <- province
    regional_triggering <- bind_rows(regional_triggering, result)
  }
}

cat("=== TRIGGERING BY PROVINCE (Top 10) ===\n")
print(regional_triggering %>%
      select(province, mean_triggered, n_focal) %>%
      arrange(desc(mean_triggered)))
cat("\n")

# Save
write.csv(regional_triggering, "regional_triggering.csv", row.names = FALSE)
cat("Saved: regional_triggering.csv\n")

# ====================================================================
# 7. SUMMARY OF FINDINGS
# ====================================================================

cat("\n=== PRELIMINARY FINDINGS ON CONTAGION CHARACTERISTICS ===\n\n")

# Compare violent vs peaceful
violent_trigger <- triggering_comparison %>%
  filter(condition == "is_violent == TRUE") %>%
  pull(mean_triggered)

peaceful_trigger <- triggering_comparison %>%
  filter(condition == "is_peaceful == TRUE") %>%
  pull(mean_triggered)

cat("Violence Effect:\n")
if(length(violent_trigger) > 0 & length(peaceful_trigger) > 0) {
  cat(sprintf("  Violent protests trigger: %.2f events on average\n", violent_trigger))
  cat(sprintf("  Peaceful protests trigger: %.2f events on average\n", peaceful_trigger))
  cat(sprintf("  Ratio (Violent/Peaceful): %.2f\n", violent_trigger/peaceful_trigger))

  if(violent_trigger > peaceful_trigger) {
    cat("  → Violent protests appear MORE contagious (mobilization hypothesis)\n")
  } else {
    cat("  → Peaceful protests appear MORE contagious (safety hypothesis)\n")
  }
}
cat("\n")

# Compare fatalities
fatal_trigger <- triggering_comparison %>%
  filter(condition == "fatalities > 0") %>%
  pull(mean_triggered)

nonfatal_trigger <- triggering_comparison %>%
  filter(condition == "fatalities == 0") %>%
  pull(mean_triggered)

cat("Fatality Effect:\n")
if(length(fatal_trigger) > 0 & length(nonfatal_trigger) > 0) {
  cat(sprintf("  Protests with fatalities trigger: %.2f events\n", fatal_trigger))
  cat(sprintf("  Protests without fatalities trigger: %.2f events\n", nonfatal_trigger))
  cat(sprintf("  Ratio (Fatal/Non-fatal): %.2f\n", fatal_trigger/nonfatal_trigger))

  if(fatal_trigger > nonfatal_trigger) {
    cat("  → Fatalities INCREASE contagion (martyrdom effect)\n")
  } else {
    cat("  → Fatalities DECREASE contagion (deterrence effect)\n")
  }
}
cat("\n")

# Compare state intervention
state_yes <- triggering_comparison %>%
  filter(condition == "state_intervention == TRUE") %>%
  pull(mean_triggered)

state_no <- triggering_comparison %>%
  filter(condition == "state_intervention == FALSE") %>%
  pull(mean_triggered)

cat("State Intervention Effect:\n")
if(length(state_yes) > 0 & length(state_no) > 0) {
  cat(sprintf("  With state intervention: %.2f events triggered\n", state_yes))
  cat(sprintf("  Without state intervention: %.2f events triggered\n", state_no))
  cat(sprintf("  Ratio (With/Without): %.2f\n", state_yes/state_no))

  if(state_yes > state_no) {
    cat("  → State intervention INCREASES contagion (backlash effect)\n")
  } else {
    cat("  → State intervention DECREASES contagion (deterrence effect)\n")
  }
}
cat("\n")

cat("INTERPRETATION:\n")
cat("These are descriptive statistics. For causal inference, we need:\n")
cat("  1. Self-exciting point process models (Hawkes/ETAS)\n")
cat("  2. Control for confounders (location, time, etc.)\n")
cat("  3. Statistical significance tests\n")
cat("  4. Conditional intensity modeling\n\n")

cat("=== NEXT STEP: Fit marked point process models ===\n")
cat("We will model: λ(s,t|H_t) where triggering depends on event marks\n\n")

cat("=== COMPLETE ===\n")
