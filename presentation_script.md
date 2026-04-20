# Presentation Script: Do Protests Really Spread?

**Target length:** 25-30 minutes

---

## Slide 1: Title

> Thank you for having me. Today I'm going to talk about protest diffusion—specifically, whether protests actually spread from place to place, or whether the clustering we observe is just regions responding to the same underlying conditions.
>
> I'll be presenting evidence from Indonesia using a spatial Hawkes process, which lets me separate these two explanations and quantify how far protest influence travels.

---

## Slide 2: Research Question

> The central question is deceptively simple: do protests really spread? We've all seen protest waves where demonstrations seem to cascade across cities and regions. But there's a puzzle here that the literature has struggled with.
>
> When we observe protests clustering in space and time, it could reflect true contagion—protests actually causing other protests through media coverage, activist networks, or demonstration effects. Or it could just be common shocks—regions responding independently to the same national events like elections or price spikes.
>
> So I'm asking two specific sub-questions. First, do protests elsewhere predict more protests here? And second, if they do, does geography matter—is the spread local or national?

---

## Slide 3: Prior Work on Protest Diffusion

> There's a rich qualitative literature on this. Tarrow's work on protest cycles, Beissinger on revolutionary cascades, McAdam and Soule on tactical diffusion—these scholars have documented how protest tactics and mobilization spread across movements and regions.
>
> The proposed mechanisms include media coverage, organizational networks, and demonstration effects where success in one place emboldens activists elsewhere.
>
> But the quantitative evidence is more limited. The early event-history models from Spilerman and Myers were important first steps, but they have limitations. They can't easily separate local versus national spread, and they can't quantify the rate of diffusion in an interpretable way.

---

## Slide 4: The Identification Problem

> This brings me to the core identification challenge. Protests cluster in space and time—but why?
>
> Standard approaches struggle to distinguish three explanations. First, common shocks: a national event like an election or a fuel price increase triggers protests everywhere simultaneously. Second, correlated fundamentals: regions with similar demographics or economic conditions protest for similar structural reasons. Third, true contagion: protests in one place actually cause protests elsewhere.
>
> These are observationally similar but theoretically distinct. If we see Jakarta and Surabaya both protest on the same week, we need to know whether one caused the other or whether both were responding to the same national news.

---

## Slide 5: Approach — Discrete Spatial Hawkes Process

> My approach uses a Hawkes process, which comes from seismology—it was originally developed to model earthquake aftershocks.
>
> The key insight is that Hawkes processes decompose event intensity into two components. There's a background rate driven by structural factors—this captures why some places protest more than others regardless of what's happening elsewhere. And there's a triggered component—events caused by prior events.
>
> Crucially, I can incorporate spatial structure. By comparing different spatial kernels, I can test whether spread is local—with nearby regions mattering more—or national, with all regions weighted equally. The parameter theta captures the characteristic decay distance.
>
> This approach builds on Hawkes's original 1971 paper, adapted for discrete count data following work by Mohler and Schoenberg.

---

## Slide 6: Scope — Cross-Region Contagion

> Before getting into the model details, let me clarify two simplifying choices.
>
> First, on event definition: I aggregate protests that occur within 7 days at the same location into a single event. This avoids double-counting what's really just continuation of the same protest, and focuses our analysis on distinct mobilization episodes.
>
> Second, I focus specifically on cross-region contagion—whether protests elsewhere predict protests here. This is a cleaner test because it controls for local persistence, and it directly tests spatial diffusion.
>
> There is a limitation to acknowledge: this design cannot study within-region persistence. I can't answer whether protests in Jakarta lead to more protests in Jakarta. That's a different question requiring a different identification strategy.

---

## Slide 7: Model Specification

> Here's the formal model. The baseline model M0 is an inhomogeneous Poisson with just covariates—population, poverty rate, CPI, and year fixed effects. This captures the background rate.
>
> The Hawkes model M1 adds the cross-region contagion term. The key addition is alpha times this cross-region exposure measure.
>
> The exposure term sums over all other regions, weighting each by the spatial kernel, and counts how many protests occurred there over the past 30 days. So if alpha is positive and significant, that means protests elsewhere predict more protests here, controlling for structural factors.

---

## Slide 8: Spatial Kernel Specification

> The spatial kernel determines how influence decays with distance. I use an exponential decay kernel, where the weight between regions r and s declines exponentially with the distance between them.
>
> The intuition is that nearby regions matter more than distant ones. Theta is the characteristic decay distance in kilometers. At distance theta, influence drops to about 37 percent of its maximum. At twice theta, it's down to 14 percent. At three times theta, just 5 percent.
>
> The key test is whether theta is finite or infinite. If finite, spread is local—nearby regions matter more. If infinite, all regions are weighted equally, which would suggest national-level spread rather than geographic diffusion.

---

## Slide 9: Estimation Strategy

> For estimation, I use quasi-Poisson GLM with a log link. This accommodates overdispersion—the estimated dispersion parameter is about 1.05, so there's mild overdispersion but nothing severe. The dispersion-adjusted standard errors handle this without needing clustered inference.
>
> For theta, I use profile likelihood. I fix theta at a grid of values from 50 kilometers up to infinity, estimate the other parameters conditional on each theta, and select the theta that minimizes deviance.
>
> I then conduct two hypothesis tests. First, I test whether alpha equals zero—that's testing whether contagion exists at all. Second, I test whether theta equals infinity—that's testing whether distance decay matters or whether spread is uniform across the country.

---

## Slide 10: Data — ACLED Indonesia

> The data come from ACLED—the Armed Conflict Location and Event Data Project. I focus on Indonesia from January 2015 through October 2024.
>
> The sample includes about 12,500 protest events across 434 districts—these are the kabupaten and kota administrative units. This gives me 1.56 million district-day observations.
>
> For covariates, I include log population and poverty rate at the district level, CPI at the national monthly level to capture economic conditions, and year fixed effects.

---

## Slide 11: Spatial Distribution

> This map shows where protests occur across Indonesia. Each dot is a protest event.
>
> You can immediately see the concentration in Java—that's the densely populated central region. The eastern islands are much sparser. This spatial heterogeneity is exactly why we need to control for population and other structural factors in the background rate.

---

## Slide 12: Temporal Distribution

> The time series shows daily protest counts over the study period.
>
> Notice the clustering and temporal bursts—there are periods of relative quiet followed by spikes of activity. This pattern is consistent with self-exciting dynamics, where events beget more events. But of course, it could also reflect common shocks. The model will help us distinguish these.

---

## Slide 13: Test 1 — Does Cross-District Contagion Exist?

> Now for the results. The first test asks whether cross-district contagion exists at all.
>
> The table compares the baseline model M0 with the Hawkes model M1. Adding the contagion term reduces deviance from about 108,700 to 108,500—a reduction of about 200 points for one additional parameter.
>
> The F-test for this nested comparison gives an F-statistic of 185 with a p-value below 10 to the minus 41. So there's overwhelming evidence that the contagion term improves model fit. Protests elsewhere do predict more protests here.

---

## Slide 14: Coefficient Estimates

> Here's the full coefficient table. The covariates behave as expected: larger populations see more protests, higher CPI is associated with more protests, and interestingly, higher poverty rates are associated with fewer protests—perhaps reflecting capacity constraints on mobilization.
>
> The key coefficient is alpha, the cross-district contagion effect, estimated at 0.32 with a standard error of 0.023. Highly significant.
>
> In multiplicative terms, exp of 0.32 is about 1.38, so a one-unit increase in cross-district exposure multiplies the expected protest rate by 38 percent.

---

## Slide 15: Magnitude — Branching Ratio and Endogenous Fraction

> Beyond statistical significance, we want to know the magnitude. How much does contagion actually matter?
>
> The branching ratio tells us how many additional protests each protest triggers elsewhere. The estimate is 0.11—so each protest triggers about one-tenth of an additional protest in other districts.
>
> Importantly, this is subcritical—the branching ratio is less than one—meaning the process is stable, not explosive. Protests don't cascade infinitely.
>
> The endogenous fraction tells us what share of protests are due to contagion versus background factors. About 14 percent of protests are attributable to cross-district contagion; the remaining 86 percent are driven by exogenous structural factors.
>
> So contagion is real and statistically significant, but most protest activity is still driven by local conditions rather than spatial diffusion.

---

## Slide 16: Test 2 — Does Distance Matter?

> The second test asks whether distance matters—is spread local or national?
>
> This figure shows the profile likelihood over theta. The x-axis is the decay distance in kilometers on a log scale, with infinity on the right.
>
> The minimum deviance occurs at theta equals 100 kilometers, marked by the red dot. This is our optimal estimate of the characteristic decay distance.

---

## Slide 17: Test 2 — Formal Hypothesis Test

> To formally test whether distance decay matters, I compare theta equals 100 kilometers to theta equals infinity.
>
> At 100 km, deviance is 108,518. At infinity, it's 108,591—about 73 points higher.
>
> The likelihood ratio test gives a chi-squared statistic of 69 with p less than 0.001. So we reject uniform national weighting. Distance decay is significant—spread is local, not national.
>
> Notice also that the alpha coefficient changes substantially: 0.32 at 100 km versus 1.02 at infinity. The model with local spread fits better and gives a more modest contagion effect.

---

## Slide 18: Conclusion

> Let me summarize the three main findings.
>
> First, recent protests elsewhere do predict protests today. The contagion effect is highly significant with an F-statistic of 185.
>
> Second, this effect is local, not national. The optimal decay distance is 100 kilometers, and we can reject uniform national spread. Protests in nearby districts matter more than protests far away.
>
> Third, the effect is modest but significant. The branching ratio is 0.11, meaning about 14 percent of protests are triggered by cross-district contagion. Most protest activity is still driven by structural factors.
>
> The methodological contribution is a discrete spatial-temporal Hawkes model that can distinguish local from national diffusion and quantify the rate of spread.

---

## Slide 19: Next Steps — Violent vs. Non-Violent Events

> Looking ahead, I'm interested in whether violent and non-violent events have different mobilization potential.
>
> The approach would be to separate events by type—protests versus riots or violence—and estimate cross-type contagion effects. Do non-violent protests trigger more protests? Do violent events trigger more violence, or do they actually suppress mobilization by raising the costs of participation?
>
> This matters for the literature on tactical choice in social movements. It would provide quantitative evidence on whether violence escalates or dampens collective action.

---

## Slide 20: Thank You

> Thank you. I'm happy to take questions.

---

## Appendix Notes

*For appendix slides, use if questions arise:*

- **Temporal Lag Structure**: Effect is concentrated in Week 1 (lags 1-7 days); fades rapidly after that.
- **Covariate Effects**: Population positive, CPI positive, Poverty negative.
- **Profile Likelihood Detail**: Full grid of theta values tested, showing minimum at 100 km.
