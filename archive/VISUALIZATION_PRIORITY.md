# Visualization Priority Guide

## If you can only make 3 figures, make these:

### FIGURE 1: The Story in One Image (CRITICAL)
**Branching Ratios by Event Type**
- Bar chart: 4 bars (Peaceful-nostate, Violent-nostate, Peaceful-state, Violent-state)
- Y-axis: Expected offspring per event
- Horizontal line at 1.0 (criticality)
- Shows: Only peaceful protests are self-sustaining
- **Impact**: Single most important finding about temporal contagion

### FIGURE 2: Spatial Null Result (CRITICAL)
**Heatmap of Boundary Convergence**
- Rows: All tested configurations (Power-law 500km, Gaussian 500km/100km/1000km)
- Columns: σ, β₀, decay
- Colors: Red = boundary, Green = interior
- Shows: ALL red → no spatial contagion
- **Impact**: Visually demonstrates robustness of null finding

### FIGURE 3: Map + Timeline (CONTEXT)
**Study Area and Temporal Pattern**
- Panel A: Map of Indonesia with protest locations (peaceful=blue, violent=red)
- Panel B: Daily time series showing temporal clustering
- **Impact**: Sets up the problem - lots of events, clear temporal patterns

---

## If you have time for 2 more:

### FIGURE 4: Model Selection
**Phase 2 Model Comparison (Bar Chart)**
- X: Models M0-M7
- Y: AIC or BIC
- Lower is better
- Highlight M5 as winner
- **Impact**: Justifies using M5 for interpretation

### FIGURE 5: Convergence Diagnostics
**Iteration Count Histogram**
- Shows most models converge in <30 iterations
- Addresses "insufficient iterations" concern
- **Impact**: Builds confidence in results

---

## Tables to Include:

### TABLE 1: Phase 2 Parameters (M5)
```
Parameter | Estimate | Interpretation
----------|----------|----------------
μ         | 0.272    | Background rate (events/day)
β₀        | -3.48    | Base excitation
β_violence| -0.305   | Violence REDUCES triggering
β_state   | -0.284   | State intervention REDUCES triggering
decay     | 0.030    | Memory ≈ 33 days
```

### TABLE 2: Spatial Parameters Across All Configs
```
Configuration        | σ (km) | At Bound? | LogLik
---------------------|--------|-----------|--------
Power-law 500km      | 0.1    | ✓         | 8629
Gaussian 500km       | 1.0    | ✓         | 8629
Gaussian 100km       | 1.0    | ✓         | 8629
Gaussian 1000km      | 1.0    | ✓         | 8629
```
**Message**: All configurations → same result

### TABLE 3: Likelihood Ratio Tests
```
Hypothesis               | LR    | p-value  | Sig?
-------------------------|-------|----------|------
Violence effect          | 24.4  | <0.001   | ***
State intervention       | 26.2  | <0.001   | ***
Fatality effect          | 2.1   | 0.146    | ns
```

---

## Visual Design Principles:

1. **Color coding**:
   - Peaceful protests: Blue (#3498db)
   - Violent protests: Red (#e74c3c)
   - State intervention: Orange (#f39c12)
   - Boundary convergence: Red
   - Interior convergence: Green

2. **Annotations**:
   - Always label critical thresholds (e.g., branching ratio = 1.0)
   - Add brief interpretations as subtitles
   - Use arrows/text to highlight key findings

3. **Consistency**:
   - Use same color scheme across all figures
   - Same category order (Peaceful → Violent → State)
   - Same fonts and sizes

---

## What NOT to include:

❌ **Skip these unless specifically asked:**
- Raw parameter traces from optimization
- Individual model diagnostics (M1, M2, M3, M4, M6, M7)
- Distance matrix visualizations
- Individual checkpoint files
- Detailed code (use code folding in RMarkdown)

---

## One-Slide Summary Version:

If you need a single PowerPoint slide:

**Title**: "Temporal but Not Spatial: Indonesian Protest Contagion"

**Layout** (2x2 grid):
1. **Top-left**: Map showing no spatial pattern
2. **Top-right**: Time series showing temporal clustering
3. **Bottom-left**: Branching ratios (peaceful > violent)
4. **Bottom-right**: Boundary convergence heatmap (all red)

**Key message box**:
- ✓ Temporal contagion (p<0.001)
- ✗ No spatial spread
- 🕊️ Peaceful > Violent (1.01 vs 0.75 offspring)

---

## Data Availability Statement:

Include in report:
```
Data: ACLED (Armed Conflict Location & Event Data Project)
Period: 2017-2021
Events: 16,467 protests in Indonesia
Variables: Location, date, violence, state intervention
```

---

## Key Messages for Different Audiences:

**For academics:**
- Robust null finding for spatial contagion
- Mark-dependent temporal excitation
- Methodological contribution: multi-scale robustness tests

**For policymakers:**
- Protests respond to national triggers, not local spread
- Geographic containment strategies won't work
- Address root causes (economic, political grievances)
- Peaceful movements more sustainable than violent ones

**For journalists:**
- Protests cluster in time (waves), not space
- Violence backfires - reduces follow-on mobilization
- Indonesia's protests are about national issues, not local contagion
