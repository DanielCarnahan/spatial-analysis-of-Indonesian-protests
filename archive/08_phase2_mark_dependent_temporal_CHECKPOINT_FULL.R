# Phase 2: Mark-Dependent Temporal Hawkes Models - CHECKPOINT FULL DATASET
# Test hypotheses about violence, fatalities, and state intervention
# This version uses ALL 16,467 events with checkpoint saving
# Author: Daniel Carnahan
# Date: 2025-10-28

library(dplyr)
library(ggplot2)

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  PHASE 2: MARK-DEPENDENT TEMPORAL HAWKES (FULL DATASET)     ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat("CHECKPOINT MODE:\n")
cat("  - Saves progress after each model fit\n")
cat("  - Can resume if interrupted\n")
cat("  - Safe for computer sleep/closure\n\n")

cat("FULL DATASET:\n")
cat("  - Using ALL 16,467 events (100% of data)\n")
cat("  - Max 50 iterations per model (relaxed convergence)\n")
cat("  - Expected time: 30-60 hours total (3-8 hours per model)\n")
cat("  - Higher statistical power than sample versions\n\n")

cat("HYPOTHESES TO TEST:\n")
cat("  H1: Violence effect - peaceful vs violent triggering\n")
cat("  H2: Fatality effect - deaths suppress/amplify contagion\n")
cat("  H3: State intervention - repression deters/amplifies\n\n")

# ====================================================================
# CHECKPOINT FUNCTIONS
# ====================================================================

checkpoint_dir <- "checkpoints_phase2_full"
if(!dir.exists(checkpoint_dir)) {
  dir.create(checkpoint_dir)
  cat(sprintf("✓ Created checkpoint directory: %s\n", checkpoint_dir))
}

save_model_checkpoint <- function(model_fit, model_id) {
  checkpoint_file <- file.path(checkpoint_dir, sprintf("model_%s.rds", model_id))
  saveRDS(model_fit, checkpoint_file)
  cat(sprintf("✓ Checkpoint saved: %s\n\n", checkpoint_file))
}

load_model_checkpoint <- function(model_id) {
  checkpoint_file <- file.path(checkpoint_dir, sprintf("model_%s.rds", model_id))
  if(file.exists(checkpoint_file)) {
    return(readRDS(checkpoint_file))
  }
  return(NULL)
}

model_is_complete <- function(model_id) {
  checkpoint_file <- file.path(checkpoint_dir, sprintf("model_%s.rds", model_id))
  return(file.exists(checkpoint_file))
}

# ====================================================================
# 1. LOAD DATA
# ====================================================================

cat("\n=== LOADING DATA ===\n")
protests <- readRDS("protests_prepared.rds")

# Prepare marks (characteristics)
marks_data <- protests %>%
  select(
    event_id = event_id_cnty,
    time = days_since_start,
    is_violent,
    is_peaceful,
    fatalities,
    state_intervention,
    longitude,
    latitude
  ) %>%
  arrange(time)

cat("Loaded", nrow(marks_data), "events\n")
cat("Time range:", min(marks_data$time), "to", max(marks_data$time), "days\n\n")

# Summary of marks
cat("Mark Statistics:\n")
cat(sprintf("  Violent events: %d (%.1f%%)\n",
            sum(marks_data$is_violent),
            100*mean(marks_data$is_violent)))
cat(sprintf("  Peaceful events: %d (%.1f%%)\n",
            sum(marks_data$is_peaceful),
            100*mean(marks_data$is_peaceful)))
cat(sprintf("  Events with fatalities: %d (%.1f%%)\n",
            sum(marks_data$fatalities > 0),
            100*mean(marks_data$fatalities > 0)))
cat(sprintf("  Events with state intervention: %d (%.1f%%)\n",
            sum(marks_data$state_intervention),
            100*mean(marks_data$state_intervention)))
cat("\n")

# ====================================================================
# 2. MARK-DEPENDENT HAWKES MODEL IMPLEMENTATION
# ====================================================================

cat("=== IMPLEMENTING MARK-DEPENDENT HAWKES ===\n\n")

# Model: λ(t | H_t, marks) = μ + Σ α(marks_i) · exp(-β(t - t_i))
# where: α(marks) = exp(β₀ + β₁·violent + β₂·fatalities + β₃·intervention)

mark_hawkes_functions <- function() {

  # Conditional intensity with mark-dependent excitation
  lambda_mark <- function(idx, times, marks, params) {

    mu <- params[["mu"]]
    beta_0 <- params[["beta_0"]]
    beta_violence <- params[["beta_violence"]]
    beta_fatalities <- params[["beta_fatalities"]]
    beta_state <- params[["beta_state"]]
    decay <- params[["decay"]]

    # Background
    intensity <- mu

    # Triggered component
    if(idx > 1) {
      for(j in 1:(idx-1)) {

        # Mark-dependent excitation
        alpha_j <- exp(beta_0 +
                      beta_violence * marks$is_violent[j] +
                      beta_fatalities * marks$fatalities[j] +
                      beta_state * marks$state_intervention[j])

        # Temporal kernel
        tau <- times[idx] - times[j]
        if(tau > 0) {
          intensity <- intensity + alpha_j * exp(-decay * tau)
        }
      }
    }

    return(intensity)
  }

  # Log-likelihood
  loglik_mark <- function(params, times, marks) {

    # Parameter constraints
    if(params[["mu"]] <= 0 || params[["decay"]] <= 0) return(-Inf)

    n <- length(times)
    ll <- 0

    # Sum log-intensities at event times
    for(i in 1:n) {
      lambda_i <- lambda_mark(i, times, marks, params)
      ll <- ll + log(max(lambda_i, 1e-10))

      # Progress indicator
      if(i %% 1000 == 0) {
        cat(sprintf("    Event %d/%d (%.1f%%)\r", i, n, 100*i/n))
      }
    }

    # Compensator (integral of intensity)
    T_max <- max(times)
    compensator <- params[["mu"]] * T_max

    # Triggered component contribution
    for(i in 1:n) {
      if(times[i] < T_max) {

        # Mark-dependent alpha for event i
        alpha_i <- exp(params[["beta_0"]] +
                      params[["beta_violence"]] * marks$is_violent[i] +
                      params[["beta_fatalities"]] * marks$fatalities[i] +
                      params[["beta_state"]] * marks$state_intervention[i])

        # Integral of exponential kernel
        tau_remain <- T_max - times[i]
        compensator <- compensator +
          (alpha_i / params[["decay"]]) * (1 - exp(-params[["decay"]] * tau_remain))
      }
    }

    ll_final <- ll - compensator

    cat(sprintf("\n    LL: %.2f (events: %.2f, comp: %.2f)\n", ll_final, ll, -compensator))

    return(ll_final)
  }

  return(list(
    lambda_mark = lambda_mark,
    loglik_mark = loglik_mark
  ))
}

hawkes_funcs <- mark_hawkes_functions()

# ====================================================================
# 3. PREPARE FULL DATASET (NO SAMPLING)
# ====================================================================

cat("=== PREPARING FULL DATASET ===\n")
cat("Using ALL events for maximum statistical power\n\n")

# Use full dataset
times_full <- marks_data$time
marks_full <- marks_data

cat(sprintf("Using full dataset: %d events (100%% of data)\n\n", nrow(marks_full)))

# ====================================================================
# 4. OPTIMIZATION FUNCTION
# ====================================================================

# Optimization function
fit_mark_hawkes_optim <- function(times, marks, include_violence = FALSE,
                                  include_fatalities = FALSE,
                                  include_state = FALSE,
                                  model_name = "Model") {

  cat(sprintf("\n--- Fitting %s ---\n", model_name))
  start_time <- Sys.time()

  # Initial parameters
  # From previous analysis, we know approximate values
  params_init <- list(
    mu = 0.1,
    beta_0 = log(0.1),
    beta_violence = 0,
    beta_fatalities = 0,
    beta_state = 0,
    decay = 0.2
  )

  # Which parameters to estimate?
  if(!include_violence) params_init$beta_violence <- 0
  if(!include_fatalities) params_init$beta_fatalities <- 0
  if(!include_state) params_init$beta_state <- 0

  # Parameter vector for optimization
  param_names <- c("mu", "beta_0", "decay")
  param_vec <- c(params_init$mu, params_init$beta_0, params_init$decay)

  if(include_violence) {
    param_names <- c(param_names, "beta_violence")
    param_vec <- c(param_vec, params_init$beta_violence)
  }
  if(include_fatalities) {
    param_names <- c(param_names, "beta_fatalities")
    param_vec <- c(param_vec, params_init$beta_fatalities)
  }
  if(include_state) {
    param_names <- c(param_names, "beta_state")
    param_vec <- c(param_vec, params_init$beta_state)
  }

  names(param_vec) <- param_names

  # Objective function for optim
  obj_fun <- function(par) {

    # Reconstruct full parameter list
    params_full <- list(
      mu = par[["mu"]],
      beta_0 = par[["beta_0"]],
      beta_violence = if(include_violence) par[["beta_violence"]] else 0,
      beta_fatalities = if(include_fatalities) par[["beta_fatalities"]] else 0,
      beta_state = if(include_state) par[["beta_state"]] else 0,
      decay = par[["decay"]]
    )

    # Return negative log-likelihood for minimization
    -hawkes_funcs$loglik_mark(params_full, times, marks)
  }

  # Optimize
  cat("  Starting optimization (max 50 iterations, 3-8 hours per model)...\n")

  fit <- optim(
    par = param_vec,
    fn = obj_fun,
    method = "L-BFGS-B",
    lower = c(mu = 1e-6, beta_0 = -5, decay = 0.01,
              beta_violence = -3, beta_fatalities = -3, beta_state = -3)[param_names],
    upper = c(mu = 10, beta_0 = 5, decay = 2,
              beta_violence = 3, beta_fatalities = 3, beta_state = 3)[param_names],
    control = list(maxit = 50, factr = 1e10, trace = 1, REPORT = 10)
  )

  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")

  cat(sprintf("\n  ✓ Optimization complete in %.2f minutes\n", as.numeric(runtime)))
  cat(sprintf("  Log-likelihood: %.2f\n", -fit$value))
  cat(sprintf("  Convergence: %d\n", fit$convergence))

  # Extract results
  results <- list(
    model_name = model_name,
    params = fit$par,
    loglik = -fit$value,
    convergence = fit$convergence,
    n_params = length(fit$par),
    AIC = 2*length(fit$par) - 2*(-fit$value),
    BIC = length(fit$par)*log(length(times)) - 2*(-fit$value),
    runtime_mins = as.numeric(runtime)
  )

  return(results)
}

# ====================================================================
# 5. FIT MODELS WITH CHECKPOINTS
# ====================================================================

cat("\n╔══════════════════════════════════════════╗\n")
cat("║      FITTING MODELS WITH CHECKPOINTS    ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# Define all models to fit
models_to_fit <- list(
  list(id = "M0", name = "M0: Baseline (no marks)",
       violence = FALSE, fatalities = FALSE, state = FALSE),
  list(id = "M1", name = "M1: Violence effect",
       violence = TRUE, fatalities = FALSE, state = FALSE),
  list(id = "M2", name = "M2: Fatality effect",
       violence = FALSE, fatalities = TRUE, state = FALSE),
  list(id = "M3", name = "M3: State intervention",
       violence = FALSE, fatalities = FALSE, state = TRUE),
  list(id = "M4", name = "M4: Violence + Fatalities",
       violence = TRUE, fatalities = TRUE, state = FALSE),
  list(id = "M5", name = "M5: Violence + State",
       violence = TRUE, fatalities = FALSE, state = TRUE),
  list(id = "M6", name = "M6: Fatalities + State",
       violence = FALSE, fatalities = TRUE, state = TRUE),
  list(id = "M7", name = "M7: Full model (all marks)",
       violence = TRUE, fatalities = TRUE, state = TRUE)
)

# Check which models are already complete
completed_models <- sapply(models_to_fit, function(m) model_is_complete(m$id))

if(any(completed_models)) {
  cat("RESUMING FROM CHECKPOINT:\n")
  cat(sprintf("  Already completed: %s\n",
              paste(sapply(models_to_fit[completed_models], function(m) m$id), collapse = ", ")))
  cat(sprintf("  Remaining: %s\n\n",
              paste(sapply(models_to_fit[!completed_models], function(m) m$id), collapse = ", ")))
} else {
  cat("Starting fresh (no checkpoints found)\n\n")
}

# Fit each model (or load from checkpoint)
all_fits <- list()

for(i in 1:length(models_to_fit)) {
  model_spec <- models_to_fit[[i]]

  if(model_is_complete(model_spec$id)) {
    cat(sprintf("⏩ %s: Loading from checkpoint\n", model_spec$name))
    all_fits[[model_spec$id]] <- load_model_checkpoint(model_spec$id)
  } else {
    # Fit the model
    fit_result <- fit_mark_hawkes_optim(
      times_full, marks_full,
      include_violence = model_spec$violence,
      include_fatalities = model_spec$fatalities,
      include_state = model_spec$state,
      model_name = model_spec$name
    )

    # Save checkpoint immediately
    all_fits[[model_spec$id]] <- fit_result
    save_model_checkpoint(fit_result, model_spec$id)
  }
}

cat("\n✓ All models complete!\n\n")

# ====================================================================
# 6. MODEL COMPARISON
# ====================================================================

cat("\n╔══════════════════════════════════════════╗\n")
cat("║         MODEL COMPARISON                ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# Compile results
comparison_df <- data.frame(
  Model = sapply(all_fits, function(x) x$model_name),
  LogLik = sapply(all_fits, function(x) x$loglik),
  nParams = sapply(all_fits, function(x) x$n_params),
  AIC = sapply(all_fits, function(x) x$AIC),
  BIC = sapply(all_fits, function(x) x$BIC)
) %>%
  arrange(AIC)

print(comparison_df)

# Best model
best_aic <- comparison_df$Model[1]
cat(sprintf("\n✓ Best model by AIC: %s\n\n", best_aic))

# Save comparison
write.csv(comparison_df, "model_comparison_phase2_full.csv", row.names = FALSE)
cat("✓ Saved: model_comparison_phase2_full.csv\n\n")

# ====================================================================
# 7. LIKELIHOOD RATIO TESTS
# ====================================================================

cat("╔══════════════════════════════════════════╗\n")
cat("║      LIKELIHOOD RATIO TESTS             ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# Function for LR test
lr_test <- function(fit_null, fit_alt, hypothesis_name) {

  LR_stat <- 2 * (fit_alt$loglik - fit_null$loglik)
  df <- fit_alt$n_params - fit_null$n_params
  p_value <- 1 - pchisq(LR_stat, df)

  cat(sprintf("%s:\n", hypothesis_name))
  cat(sprintf("  LR statistic: %.4f\n", LR_stat))
  cat(sprintf("  df: %d\n", df))
  cat(sprintf("  p-value: %.6f", p_value))

  if(p_value < 0.001) {
    cat(" ***\n")
  } else if(p_value < 0.01) {
    cat(" **\n")
  } else if(p_value < 0.05) {
    cat(" *\n")
  } else {
    cat(" (not significant)\n")
  }

  return(data.frame(
    Hypothesis = hypothesis_name,
    LR = LR_stat,
    df = df,
    p_value = p_value,
    significant = p_value < 0.05
  ))
}

# Test each hypothesis
cat("\nH1: Violence Effect\n")
cat("--------------------\n")
test_H1 <- lr_test(all_fits$M0, all_fits$M1, "H1: Violence effect")

cat("\nH2: Fatality Effect\n")
cat("--------------------\n")
test_H2 <- lr_test(all_fits$M0, all_fits$M2, "H2: Fatality effect")

cat("\nH3: State Intervention Effect\n")
cat("------------------------------\n")
test_H3 <- lr_test(all_fits$M0, all_fits$M3, "H3: State intervention")

cat("\nFull Model vs Baseline\n")
cat("----------------------\n")
test_Full <- lr_test(all_fits$M0, all_fits$M7, "Full model vs baseline")

# Combine tests
lr_tests_df <- bind_rows(test_H1, test_H2, test_H3, test_Full)
write.csv(lr_tests_df, "likelihood_ratio_tests_phase2_full.csv", row.names = FALSE)
cat("\n✓ Saved: likelihood_ratio_tests_phase2_full.csv\n\n")

# ====================================================================
# 8. PARAMETER INTERPRETATION
# ====================================================================

cat("╔══════════════════════════════════════════╗\n")
cat("║     PARAMETER INTERPRETATION            ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# Extract parameters from full model
best_fit <- all_fits$M7

cat("Full Model Parameters:\n")
cat("----------------------\n")
print(best_fit$params)
cat("\n")

# Interpret beta coefficients
interpret_beta <- function(beta, name) {

  multiplicative_effect <- exp(beta)

  cat(sprintf("%s:\n", name))
  cat(sprintf("  β = %.4f\n", beta))
  cat(sprintf("  exp(β) = %.4f\n", multiplicative_effect))

  if(multiplicative_effect > 1.1) {
    percent_increase <- (multiplicative_effect - 1) * 100
    cat(sprintf("  → INCREASES triggering by %.1f%%\n", percent_increase))
  } else if(multiplicative_effect < 0.9) {
    percent_decrease <- (1 - multiplicative_effect) * 100
    cat(sprintf("  → DECREASES triggering by %.1f%%\n", percent_decrease))
  } else {
    cat(sprintf("  → Minimal effect (near neutral)\n"))
  }

  cat("\n")
}

# Interpret each coefficient
if("beta_violence" %in% names(best_fit$params)) {
  interpret_beta(best_fit$params[["beta_violence"]], "Violence Effect")
}

if("beta_fatalities" %in% names(best_fit$params)) {
  interpret_beta(best_fit$params[["beta_fatalities"]], "Per-Fatality Effect")
}

if("beta_state" %in% names(best_fit$params)) {
  interpret_beta(best_fit$params[["beta_state"]], "State Intervention Effect")
}

# Save all model objects
saveRDS(all_fits, "mark_hawkes_all_models_full.rds")
cat("✓ Saved: mark_hawkes_all_models_full.rds\n\n")

# ====================================================================
# 9. VISUALIZATIONS
# ====================================================================

cat("╔══════════════════════════════════════════╗\n")
cat("║          VISUALIZATIONS                 ║\n")
cat("╚══════════════════════════════════════════╝\n\n")

# Create plots directory if needed
if(!dir.exists("plots")) {
  dir.create("plots")
}

# Plot 1: Model comparison (AIC)
png("plots/21_model_comparison_aic_full.png", width = 1200, height = 600, res = 120)

comparison_plot <- comparison_df %>%
  mutate(Model_short = gsub("M[0-9]: ", "", Model))

ggplot(comparison_plot, aes(x = reorder(Model_short, AIC), y = AIC)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_text(aes(label = sprintf("%.0f", AIC)), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "Model Comparison: Akaike Information Criterion (AIC) - FULL DATASET",
       subtitle = "Lower AIC = Better model fit | n = 16,467 events",
       x = "Model",
       y = "AIC") +
  theme_minimal() +
  theme(text = element_text(size = 11))

dev.off()
cat("✓ Saved: plots/21_model_comparison_aic_full.png\n")

# Plot 2: Parameter effects (forest plot)
png("plots/22_parameter_effects_forest_full.png", width = 1000, height = 600, res = 120)

# Extract from full model
params_full <- best_fit$params

effects_plot_data <- data.frame(
  Parameter = character(),
  Beta = numeric(),
  Effect = numeric(),
  stringsAsFactors = FALSE
)

if("beta_violence" %in% names(params_full)) {
  effects_plot_data <- rbind(effects_plot_data,
                            data.frame(Parameter = "Violent (vs Peaceful)",
                                      Beta = params_full[["beta_violence"]],
                                      Effect = exp(params_full[["beta_violence"]])))
}

if("beta_fatalities" %in% names(params_full)) {
  effects_plot_data <- rbind(effects_plot_data,
                            data.frame(Parameter = "Per Fatality",
                                      Beta = params_full[["beta_fatalities"]],
                                      Effect = exp(params_full[["beta_fatalities"]])))
}

if("beta_state" %in% names(params_full)) {
  effects_plot_data <- rbind(effects_plot_data,
                            data.frame(Parameter = "State Intervention",
                                      Beta = params_full[["beta_state"]],
                                      Effect = exp(params_full[["beta_state"]])))
}

if(nrow(effects_plot_data) > 0) {
  ggplot(effects_plot_data, aes(x = Parameter, y = Effect)) +
    geom_point(size = 5, color = "darkblue") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
    geom_text(aes(label = sprintf("%.3f", Effect)), vjust = -1) +
    coord_flip() +
    labs(title = "Mark Effects on Protest Triggering Strength - FULL DATASET",
         subtitle = "Multiplicative effect on excitation parameter α | n = 16,467 events\nEffect > 1: Increases triggering; Effect < 1: Decreases triggering",
         x = "Event Characteristic",
         y = "Multiplicative Effect (exp(β))") +
    theme_minimal() +
    theme(text = element_text(size = 12))
}

dev.off()
cat("✓ Saved: plots/22_parameter_effects_forest_full.png\n\n")

# ====================================================================
# 10. SUMMARY
# ====================================================================

cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║      PHASE 2 COMPLETE: MARK-DEPENDENT TEMPORAL (FULL)       ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat("OUTPUTS:\n")
cat("  1. mark_hawkes_all_models_full.rds - All fitted models\n")
cat("  2. model_comparison_phase2_full.csv - AIC/BIC comparison\n")
cat("  3. likelihood_ratio_tests_phase2_full.csv - Hypothesis tests\n")
cat("  4. plots/21_model_comparison_aic_full.png\n")
cat("  5. plots/22_parameter_effects_forest_full.png\n")
cat("  6. checkpoints_phase2_full/ - Individual model checkpoints\n\n")

cat("KEY FINDINGS:\n")
cat("  Best model:", best_aic, "\n")
cat("  Statistical significance:\n")
print(lr_tests_df %>% select(Hypothesis, p_value, significant))

cat("\n\nNEXT STEP:\n")
cat("  → Phase 3: Add spatial component to mark-dependent model\n")
cat("  → Full spatial-temporal ETAS with marks\n\n")

cat("NOTE: To clean up checkpoints after successful completion:\n")
cat("  Run: unlink('checkpoints_phase2_full', recursive = TRUE)\n\n")

cat("=== END PHASE 2 FULL DATASET ===\n")
