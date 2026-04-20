# Recommended Visualizations for RMarkdown Report

## FIGURE 1: Study Area and Data
**Panel A: Spatial Distribution Map**
- Map of Indonesia with all 16,467 protest locations
- Color by type: peaceful (blue), violent (red)
- Size by state intervention
- Shows spatial coverage and clustering

**Panel B: Temporal Pattern**
- Time series of daily protest counts (2017-2021)
- Stacked by type (peaceful/violent)
- Shows temporal clustering and trends

## FIGURE 2: Phase 2 - Mark-Dependent Temporal Triggering
**Panel A: Excitation Coefficients (Bar Chart)**
- 4 bars: Peaceful-nostate, Violent-nostate, Peaceful-state, Violent-state
- Y-axis: Alpha (excitation strength)
- Error bars if available
- Shows peaceful > violent

**Panel B: Branching Ratios (Bar Chart)**
- Same 4 categories
- Y-axis: Expected offspring per event
- Horizontal line at 1.0 (criticality threshold)
- Shows only peaceful-nostate > 1.0

**Panel C: Model Comparison (Table)**
- Models M0-M7
- LogLik, AIC, BIC, nParams
- Highlight best model (M5)
- p-values for LR tests

## FIGURE 3: Phase 3 - Spatial Contagion Evidence
**Panel A: Parameter Convergence Across Configurations**
- Grid showing all tested configurations:
  - Rows: Power-law 500km, Gaussian 500km, Gaussian 100km, Gaussian 1000km
  - Columns: β₀, decay, σ (or σ_violent, σ_peaceful)
- Color code: RED if at boundary, GREEN if interior
- Shows ALL configurations hit boundaries

**Panel B: Spatial Parameters (Table)**
- Rows: Different configurations
- Columns: σ_violent, σ_peaceful, ω (if power-law)
- Lower/Upper bounds shown
- Highlight boundary convergence

**Panel C: Model Fit Comparison**
- Phase 2 (temporal only) vs Phase 3 (spatial-temporal)
- Show that adding spatial components doesn't improve fit
- AIC/BIC comparison

## FIGURE 4: Convergence Diagnostics (Supplementary)
**Panel A: Iteration Counts**
- Histogram of function evaluations across all starts
- Shows most converge in <30 iterations
- Supports "true result" not "insufficient iterations"

**Panel B: Multi-start Consistency**
- Log-likelihood values from 4 random starts per model
- Shows all successful starts reach same optimum
- Demonstrates robustness

## FIGURE 5: Interpretation Diagram (Conceptual)
**Schematic showing:**
- Temporal triggering: ✓ (arrows in time)
- Spatial triggering: ✗ (no arrows in space)
- Mark effects: Peaceful > Violent
- Common national triggers: ✓ (shared μ)

## TABLE 1: Phase 2 Final Parameters (M5 - Best Model)
| Parameter | Estimate | Interpretation |
|-----------|----------|----------------|
| μ | 0.272 | Background rate (events/day) |
| β₀ | -3.48 | Base excitation (log scale) |
| β_violence | -0.305 | Violence reduces triggering |
| β_state | -0.284 | State intervention reduces triggering |
| decay | 0.030 | Memory time ≈ 33 days |

## TABLE 2: Phase 3 Robustness Tests
| Configuration | σ (km) | At Boundary? | LogLik | Conclusion |
|---------------|--------|--------------|--------|------------|
| Power-law 500km | 0.1 | ✓ | -8629 | No spatial |
| Gaussian 500km | 1.0 | ✓ | -8629 | No spatial |
| Gaussian 100km | 1.0 | ✓ | -8629 | No spatial |
| Gaussian 1000km | 1.0 | ✓ | -8629 | No spatial |

