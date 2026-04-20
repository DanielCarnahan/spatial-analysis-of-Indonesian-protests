# Statistical Methods for Contagion Analysis
## Detailed Explanation of Models and Tests

**Author:** Daniel Carnahan
**Date:** October 27, 2025

---

## Overview

We used **two complementary approaches** to analyze protest contagion:

1. **Descriptive "Triggering Potential" Analysis** (Exploratory)
2. **Hawkes Self-Exciting Point Process Models** (Formal Statistical Models)

Let me explain each in detail, what they tell us, and their limitations.

---

## Method 1: Descriptive Triggering Analysis

### What We Did

For each protest event, we counted **how many subsequent protests occurred within:**
- **Spatial radius:** 0.5 degrees (~55 km)
- **Time window:** 30 days after the event

Then we compared the **mean number of triggered events** across different event types.

### Mathematical Definition

For event `i` occurring at location `(x_i, y_i)` and time `t_i`:

```
Triggered(i) = Count of events j where:
  1. t_j > t_i  (j occurs AFTER i)
  2. t_j ≤ t_i + 30 days  (within time window)
  3. distance(i, j) < 55 km  (nearby in space)
```

Then for each characteristic (e.g., "violent"), we compute:

```
Mean_Triggered = Average of Triggered(i) for all events with that characteristic
```

### Example

```
Peaceful protests: Mean = 8.77 events triggered
Violent protests:  Mean = 6.59 events triggered
```

### What This Tells Us

This is a **descriptive statistic** that measures **association**, not causation. It tells us:

✅ **What we CAN say:**
- "On average, peaceful protests are followed by more nearby protests in the next 30 days than violent protests"
- "There is a correlation between peaceful protest type and higher subsequent activity"
- "Spatial-temporal clustering is stronger around peaceful events"

❌ **What we CANNOT say:**
- "Peaceful protests CAUSE more protests" (causal claim)
- "The difference is statistically significant" (no p-value yet)
- "This accounts for confounders" (no controls)

### Limitations

**1. Confounding:**
- Peaceful protests might occur in areas that were already going to have more protests anyway
- Time trends: If peaceful protests increased over time, they'd naturally be followed by more events
- Location effects: Peaceful protests in Jakarta vs violent in Papua

**2. Selection bias:**
- What makes a protest peaceful vs violent?
- Are they fundamentally different populations?

**3. No statistical inference:**
- No confidence intervals
- No hypothesis tests
- No p-values

**4. Arbitrary thresholds:**
- Why 55 km? Why 30 days?
- Results might change with different cutoffs

### What Would Make This Stronger?

1. **Bootstrap confidence intervals:**
```r
# Estimate uncertainty in mean differences
boot_result <- boot(data, function(data, indices) {
  sample_data <- data[indices,]
  mean_peaceful <- mean(triggered[is_peaceful])
  mean_violent <- mean(triggered[is_violent])
  return(mean_peaceful - mean_violent)
}, R = 1000)

# Get 95% CI
boot.ci(boot_result)
```

2. **Permutation test:**
```r
# Test if difference is significant
observed_diff <- mean_peaceful - mean_violent

# Shuffle labels and recompute 1000 times
null_distribution <- replicate(1000, {
  shuffled_labels <- sample(is_peaceful)
  mean(triggered[shuffled_labels]) - mean(triggered[!shuffled_labels])
})

p_value <- mean(abs(null_distribution) >= abs(observed_diff))
```

3. **Matched case-control:**
- For each violent protest, find a peaceful protest at same time/location
- Compare only these matched pairs
- Controls for confounding

---

## Method 2: Hawkes Self-Exciting Point Process

### What Is a Hawkes Process?

A **Hawkes process** is a statistical model where **events trigger more events**. It's widely used for:
- Earthquake aftershocks (original application)
- Financial market trades
- Social media cascades
- Conflict events

### The Model

The **conditional intensity function** describes the instantaneous rate of events:

```
λ(t) = μ + Σ_{t_i < t} α · exp(-β(t - t_i))
```

Where:
- **λ(t)** = instantaneous event rate at time t
- **μ** = background rate (constant baseline)
- **α** = excitation parameter (how much each event triggers)
- **β** = decay rate (how fast influence fades)
- **Σ** sums over all past events t_i < t

### Intuition

Think of each protest as:
1. **Background component:** Random "spontaneous" protests at rate μ
2. **Triggered component:** Each past protest increases the rate temporarily

The rate boost from a past event decays exponentially: `α · exp(-β · time_elapsed)`

### Key Parameters

**1. Background Rate (μ)**
- Events per day that occur "spontaneously"
- Captures baseline protest risk
- Higher μ = more endemic protest activity

**2. Excitation (α)**
- Immediate increase in rate after an event
- Higher α = stronger triggering
- If α = 0, no contagion (just Poisson process)

**3. Decay Rate (β)**
- How fast influence fades
- Units: per day
- 1/β = mean influence duration
- Higher β = faster decay

**4. Branching Ratio (α/β)**
- **MOST IMPORTANT PARAMETER**
- Average number of events directly triggered by one event
- If α/β < 1: **Subcritical** (stable, self-limiting)
- If α/β ≥ 1: **Supercritical** (explosive, cascading)

### Example: Jakarta

```
μ = 0.10 events/day
α = 0.10
β = 0.23 per day
```

**Interpretation:**
1. **Background:** 0.10 protests/day = 1 protest every 10 days spontaneously
2. **Excitation:** Each protest temporarily increases rate by 0.10
3. **Decay:** Influence halves every 3 days (ln(2)/β)
4. **Branching:** α/β = 0.44 → each protest triggers 0.44 more on average
5. **Stability:** Subcritical → stable system

### How We Estimated Parameters

We used **Maximum Likelihood Estimation (MLE)**, which finds parameters that make the observed data most probable.

**Log-Likelihood Function:**

```
L(μ, α, β) = Σ log(λ(t_i)) - ∫₀ᵀ λ(t) dt
```

Where:
- First term: Sum of log-intensities at each event time (rewards fitting events)
- Second term: Integral of intensity (penalizes high overall rate)

**Optimization:**

We used a **grid search** (simplified approach):
- Try many combinations of (μ, α, β)
- Compute log-likelihood for each
- Choose parameters with highest likelihood

**Note:** This is computationally feasible but suboptimal. Better approaches:
- Gradient-based optimization (BFGS, L-BFGS)
- Specialized Hawkes packages
- Bayesian MCMC

### What This Tells Us

✅ **What we CAN say:**

1. **Temporal self-excitation exists:**
   - If α > 0, there's evidence of triggering
   - Past events increase future rates

2. **Quantify influence duration:**
   - 1/β tells us how long influence lasts
   - Jakarta: ~4.4 days

3. **System stability:**
   - Branching ratio < 1 → won't explode
   - α/β = 0.44 → 44% of protests triggered, 56% spontaneous

4. **Compare across types:**
   - Different branching ratios → different contagion strength

❌ **What we CANNOT say:**

1. **Not causal identification:**
   - Model assumes past → future causation
   - But could be common causes (e.g., election drives all protests)

2. **Temporal only:**
   - Ignores spatial heterogeneity
   - All locations treated equally
   - Jakarta and Papua pooled together

3. **Stationarity assumption:**
   - Assumes μ, α, β constant over time
   - But we know protests increased 2015→2024
   - Violates model assumptions

4. **No covariates:**
   - Can't test "peaceful vs violent" within model
   - Need mark-dependent Hawkes

### Statistical Significance

We did NOT report p-values because:
1. MLE standard errors require Hessian matrix (second derivatives)
2. Our grid search doesn't provide this
3. Would need bootstrap or asymptotic theory

**What we SHOULD do:**

```r
# Likelihood Ratio Test
# H0: α = 0 (no excitation)
# H1: α > 0 (excitation exists)

L0 <- loglik(mu = mu_mle, alpha = 0, beta = beta_mle)
L1 <- loglik(mu = mu_mle, alpha = alpha_mle, beta = beta_mle)

LR_stat <- 2 * (L1 - L0)
# LR_stat ~ χ²(1) under H0
p_value <- 1 - pchisq(LR_stat, df = 1)
```

### Limitations of Our Hawkes Analysis

**1. Coarse estimation:**
- Grid search with 10³ combinations
- True optimum might be between grid points
- No uncertainty quantification

**2. No spatial component:**
- Current model: all past events trigger equally
- Reality: nearby events trigger more
- Need spatial-temporal Hawkes (ETAS)

**3. No mark-dependence:**
- We fit separate models for peaceful/violent
- But can't test if peaceful → violent cross-triggering
- Need multivariate Hawkes

**4. Supercritical results for peaceful/violent:**
- Branching = 1.21 (explosive)
- Likely due to:
  - Spatial heterogeneity (Jakarta inflates rate)
  - Non-stationarity (time trends)
  - Model misspecification

**5. Small sample for fatal protests:**
- Only 152 events
- High variance in estimates
- Results less reliable

---

## Comparing the Two Methods

| Aspect | Descriptive Triggering | Hawkes Process |
|--------|------------------------|----------------|
| **Type** | Exploratory | Formal statistical model |
| **Assumptions** | Minimal | Self-exciting, stationarity |
| **Output** | Mean triggered events | Parameters (μ, α, β) |
| **Causation** | No | No (but stronger inference) |
| **Inference** | None (no p-values) | Possible (we didn't do) |
| **Confounders** | Uncontrolled | Partially (time dependence) |
| **Interpretability** | High (simple counts) | Moderate (need math) |

### Why Both?

1. **Descriptive → Exploratory:** Shows patterns, generates hypotheses
2. **Hawkes → Confirmatory:** Tests temporal dependence rigorously

They **complement** each other:
- Descriptive: "Peaceful protests followed by more events"
- Hawkes: "Past protests increase future rates (α > 0)"

---

## What We Actually Concluded

Let me be precise about our claims:

### Strong Claims (Supported by Data)

✅ **"Peaceful protests are followed by more nearby protests"**
- Evidence: Mean 8.77 vs 6.59
- Type: Descriptive association
- Caveats: Not causal, unadjusted

✅ **"Temporal self-excitation exists"**
- Evidence: Hawkes α > 0
- Type: Statistical model fit
- Caveats: No significance test, could be spatial heterogeneity

✅ **"Jakarta shows extreme clustering"**
- Evidence: 33.62 vs 8.77 national average
- Type: Descriptive statistic
- Caveats: Obvious (more protests in capital)

### Moderate Claims (Partially Supported)

⚠️ **"Peaceful protests are MORE contagious"**
- Evidence: Higher triggering in both methods
- Type: Suggestive pattern
- Caveats:
  - No significance test
  - Could be confounded by location/time
  - Hawkes models show similar parameters

⚠️ **"Fatalities SUPPRESS contagion"**
- Evidence: Strong difference (3.73 vs 9.48)
- Type: Association
- Caveats:
  - Selection bias (what causes fatalities?)
  - Fatal protests might be in different contexts
  - Small sample (152 events)

### Weak Claims (Speculative)

🔶 **"Safety hypothesis supported over anger hypothesis"**
- Evidence: Pattern consistent with risk-aversion
- Type: Theoretical interpretation
- Caveats:
  - Alternative explanations exist
  - No direct test of mechanism
  - Need individual-level data

🔶 **"Repression deters (no backlash effect)"**
- Evidence: Lower triggering with state intervention
- Type: Association
- Caveats:
  - MAJOR selection bias (police intervene when protests expected)
  - Reverse causation possible
  - Need instrumental variable or natural experiment

---

## What Would Strengthen Our Conclusions?

### 1. Statistical Inference

**Add p-values and confidence intervals:**

```r
# Bootstrap for triggering comparison
boot_ci <- boot(data, function(d, i) {
  mean(d$triggered[d$is_peaceful][i]) -
  mean(d$triggered[d$is_violent][i])
}, R = 10000)

# If CI excludes 0 → significant difference
```

### 2. Spatial-Temporal Models (ETAS)

**Include spatial decay:**

```
λ(s, t) = μ(s) + Σ α · g(s - s_i) · h(t - t_i)

where:
  g(r) = spatial kernel (e.g., 1/(r² + d)^q)
  h(τ) = temporal kernel (e.g., exp(-β·τ))
```

This controls for spatial heterogeneity (Jakarta vs Papua).

### 3. Mark-Dependent Hawkes

**Allow triggering to depend on event type:**

```
λ_peaceful(t) = μ₁ + Σ α₁₁·φ(t-t_i) [if i peaceful] +
                      α₁₂·φ(t-t_i) [if i violent]

λ_violent(t)  = μ₂ + Σ α₂₁·φ(t-t_i) [if i peaceful] +
                      α₂₂·φ(t-t_i) [if i violent]
```

**α matrix:**
- α₁₁ = peaceful → peaceful
- α₁₂ = violent → peaceful
- α₂₁ = peaceful → violent
- α₂₂ = violent → violent

Can test: "Does α₁₁ > α₂₂?" (peaceful more self-exciting?)

### 4. Regression Adjustment

**Control for confounders:**

```r
# Conditional intensity with covariates
λ(s, t | X) = exp(β₀ + β₁·pop_density + β₂·is_jakarta +
                   β₃·year + ...) + Σ α·φ(s-s_i, t-t_i)
```

### 5. Causal Identification

**For causal claims, need:**

**Option A: Natural Experiment**
- Find exogenous variation in protest type
- E.g., weather on protest day affects violence
- Instrument: weather → violence → subsequent protests

**Option B: Regression Discontinuity**
- E.g., permit threshold for police intervention
- Compare protests just above/below threshold

**Option C: Difference-in-Differences**
- Policy change affecting protest type
- Compare regions with/without policy, before/after

---

## Summary: What Our Models Tell Us

### Descriptive Triggering Analysis

**Does:**
- Measures spatial-temporal clustering by event type
- Shows which events are followed by more activity
- Generates hypotheses

**Doesn't:**
- Test causation
- Control for confounders
- Provide statistical significance

**Conclusion Strength:** **Suggestive patterns, not causal**

### Hawkes Process Models

**Does:**
- Formally models temporal self-excitation
- Quantifies triggering parameters
- Decomposes background vs triggered activity

**Doesn't:**
- Account for spatial heterogeneity (our implementation)
- Test mark-dependence (our implementation)
- Establish causation (structural assumption)

**Conclusion Strength:** **Strong evidence of temporal dependence, but spatial/mark effects confounded**

---

## Honest Assessment

### What We Can Confidently Say

1. ✅ **Strong spatial-temporal clustering exists** (highly significant in quadrat test)
2. ✅ **Temporal self-excitation detected** (Hawkes α > 0)
3. ✅ **Peaceful protests followed by more activity** (descriptive fact)
4. ✅ **Jakarta exceptional** (obvious from data)
5. ✅ **Influence decays in days, not months** (Hawkes 1/β ≈ 4-12 days)

### What We Cannot Confidently Say

1. ❌ **Peaceful protests CAUSE more contagion** (confounding not addressed)
2. ❌ **Repression deters** (severe selection bias)
3. ❌ **Safety hypothesis proven** (interpretation, not test)
4. ❌ **Exact magnitudes** (no confidence intervals)

### What's Needed for Stronger Claims

1. **Statistical inference** (p-values, CIs)
2. **Spatial-temporal models** (ETAS)
3. **Covariate adjustment** (regression controls)
4. **Causal identification** (natural experiment, IV, RDD)
5. **Robustness checks** (different thresholds, time periods)

---

## Recommended Next Steps

### Immediate (Strengthen Current Analysis)

1. **Compute confidence intervals** for triggering differences
2. **Proper MLE** for Hawkes (not grid search)
3. **Likelihood ratio tests** for significance
4. **Robustness:** Try different spatial/temporal thresholds

### Medium-Term (Better Models)

1. **Fit ETAS models** with `stpp` or `etasFLP` package
2. **Mark-dependent Hawkes** with `ppmlasso`
3. **Regression adjustment** for confounders
4. **Model diagnostics** (residual analysis, Q-Q plots)

### Long-Term (Causal Inference)

1. **Identify instruments** for protest characteristics
2. **Search for natural experiments**
3. **Match on observables** (propensity score)
4. **Collect individual data** (mechanisms)

---

## References for Further Reading

### Point Process Methods

- Baddeley, A., Rubak, E., & Turner, R. (2015). *Spatial Point Patterns: Methodology and Applications with R*. CRC Press.
  - Chapter 9: Testing for spatial patterns
  - Chapter 10-11: Point process models

- Daley, D. J., & Vere-Jones, D. (2003). *An Introduction to the Theory of Point Processes*. Springer.
  - Chapter 7: Self-exciting processes

### Hawkes Processes

- Hawkes, A. G. (1971). "Spectra of some self-exciting and mutually exciting point processes." *Biometrika*, 58(1), 83-90.
  - Original paper

- Reinhart, A. (2018). "A review of self-exciting spatio-temporal point processes and their applications." *Statistical Science*, 33(3), 299-318.
  - Excellent overview

### Applications to Conflict

- Zhukov, Y. M. (2012). "Roads and the diffusion of insurgent violence." *American Journal of Political Science*, 56(2), 449-469.
  - Spatial point processes for conflict

- White, G., Ruggeri, F., & Porter, M. D. (2020). "Modeling insurgent activity using self-exciting processes." *Journal of Quantitative Criminology*, 36, 483-515.
  - Hawkes for insurgency

### Causal Inference

- Imbens, G. W., & Rubin, D. B. (2015). *Causal Inference for Statistics, Social, and Biomedical Sciences*. Cambridge.
  - Gold standard textbook

- Cunningham, S. (2021). *Causal Inference: The Mixtape*. Yale University Press.
  - Accessible intro with code

---

**End of Methods Explanation**

*The key takeaway: Our analysis provides strong evidence of patterns and associations, but causal claims require additional modeling and identification strategies.*
