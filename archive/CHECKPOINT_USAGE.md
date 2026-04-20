# Phase 2 Checkpoint Version - Usage Guide

## Overview

The checkpoint version (`08_phase2_mark_dependent_temporal_CHECKPOINT.R`) allows you to safely run Phase 2 analysis even if your computer sleeps or closes.

## Key Features

1. **Saves after each model**: Each of the 8 models is saved immediately after fitting
2. **Resume from interruption**: If interrupted, re-running the script will skip completed models
3. **Safe checkpoints**: All checkpoints stored in `checkpoints_phase2/` directory

## How to Use

### First Run

```r
Rscript 08_phase2_mark_dependent_temporal_CHECKPOINT.R
```

The script will:
- Create `checkpoints_phase2/` directory
- Fit models M0 through M7 sequentially
- Save each model as `checkpoints_phase2/model_M0.rds`, etc.
- Generate final outputs when all models complete

### If Interrupted

Simply run the same command again:

```r
Rscript 08_phase2_mark_dependent_temporal_CHECKPOINT.R
```

The script will:
- Check `checkpoints_phase2/` for completed models
- Print: "RESUMING FROM CHECKPOINT: Already completed: M0, M1, M2..."
- Skip to first incomplete model
- Continue from where it left off

### After Successful Completion

You'll have:
- `mark_hawkes_all_models.rds` - All 8 fitted models
- `model_comparison_phase2.csv` - AIC/BIC comparison table
- `likelihood_ratio_tests_phase2.csv` - Hypothesis test results
- `plots/21_model_comparison_aic.png` - Model comparison plot
- `plots/22_parameter_effects_forest.png` - Parameter effects plot
- `checkpoints_phase2/` - 8 individual checkpoint files

### Clean Up (Optional)

After successful completion, you can remove checkpoints:

```r
# In R:
unlink('checkpoints_phase2', recursive = TRUE)

# Or in terminal:
rm -rf checkpoints_phase2
```

## Example Progress

### Initial Run
```
╔══════════════════════════════════════════╗
║      FITTING MODELS WITH CHECKPOINTS    ║
╚══════════════════════════════════════════╝

Starting fresh (no checkpoints found)

--- Fitting M0: Baseline (no marks) ---
  Starting optimization (this may take 10-30 minutes)...
  [... optimization progress ...]
  ✓ Optimization complete in 15.32 minutes
✓ Checkpoint saved: checkpoints_phase2/model_M0.rds

--- Fitting M1: Violence effect ---
  [... interrupted here ...]
```

### Resumed Run
```
╔══════════════════════════════════════════╗
║      FITTING MODELS WITH CHECKPOINTS    ║
╚══════════════════════════════════════════╝

RESUMING FROM CHECKPOINT:
  Already completed: M0
  Remaining: M1, M2, M3, M4, M5, M6, M7

⏩ M0: Baseline (no marks): Loading from checkpoint

--- Fitting M1: Violence effect ---
  [... continues from where it left off ...]
```

## Differences from Original Version

| Feature | Original | Checkpoint Version |
|---------|----------|-------------------|
| Computer sleep | ❌ Loses all progress | ✅ Resumes from last model |
| Partial results | ❌ None if interrupted | ✅ Saved after each model |
| Re-run behavior | Starts from scratch | Skips completed models |
| Storage | Single final file | Checkpoints + final file |

## Technical Details

Each checkpoint file contains:
- `model_name`: e.g., "M0: Baseline (no marks)"
- `params`: Named vector of fitted parameters
- `loglik`: Log-likelihood value
- `convergence`: Optimization convergence code
- `n_params`: Number of parameters
- `AIC`: Akaike Information Criterion
- `BIC`: Bayesian Information Criterion
- `runtime_mins`: Time taken to fit this model

## Troubleshooting

**Problem**: Script says "Already completed: M0, M1, ..." but I want to re-fit
**Solution**: Delete the checkpoint directory first:
```bash
rm -rf checkpoints_phase2
```

**Problem**: Optimization fails for one model
**Solution**: Delete just that model's checkpoint and re-run:
```bash
rm checkpoints_phase2/model_M3.rds
Rscript 08_phase2_mark_dependent_temporal_CHECKPOINT.R
```

**Problem**: Want to check checkpoint contents
**Solution**: In R:
```r
checkpoint <- readRDS("checkpoints_phase2/model_M0.rds")
str(checkpoint)
```
