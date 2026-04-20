#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double hawkes_basic_negloglik_cpp(
    NumericVector times,           // Event times
    NumericVector log_pop,         // Log population for each event
    NumericVector poverty_decimal, // Poverty rate (decimal) for each event
    IntegerVector year,            // Year for each event (2015-2024)
    NumericVector district_year_exposure, // Exposure time for each district-year
    NumericVector district_year_log_pop,  // Log pop for each district-year
    NumericVector district_year_poverty,  // Poverty for each district-year
    IntegerVector district_year_year,     // Year for district-year (2015-2024)
    NumericVector params,          // Parameter vector (14 params)
    double temporal_cutoff         // Temporal cutoff in days
) {

    int n = times.size();
    int n_dy = district_year_exposure.size();

    // Extract parameters (14 total: 3 background base + 9 year effects + 2 triggering)
    double beta_0_bg = params[0];
    double gamma_raw = params[1];
    double delta_raw = params[2];
    // Year effects (9 parameters for 2016-2024)
    double beta_2016_raw = params[3];
    double beta_2017_raw = params[4];
    double beta_2018_raw = params[5];
    double beta_2019_raw = params[6];
    double beta_2020_raw = params[7];
    double beta_2021_raw = params[8];
    double beta_2022_raw = params[9];
    double beta_2023_raw = params[10];
    double beta_2024_raw = params[11];
    // Triggering parameters (NO MARK EFFECTS - constant α)
    double beta_0_trig_raw = params[12];
    double decay_raw = params[13];

    // Transform raw (unconstrained) parameters to bounded parameters using sigmoid
    // Background rate parameters (RELAXED BOUNDS):
    double gamma = 0.001 + 4.999 / (1.0 + exp(-gamma_raw));       // [0.001, 5] (was [0.01, 2])
    double delta = -20.0 + 40.0 / (1.0 + exp(-delta_raw));        // [-20, 20]

    // Year effects (9 parameters) (RELAXED BOUNDS):
    double beta_2016 = -5.0 + 10.0 / (1.0 + exp(-beta_2016_raw)); // [-5, 5] (was [-3, 3])
    double beta_2017 = -5.0 + 10.0 / (1.0 + exp(-beta_2017_raw));
    double beta_2018 = -5.0 + 10.0 / (1.0 + exp(-beta_2018_raw));
    double beta_2019 = -5.0 + 10.0 / (1.0 + exp(-beta_2019_raw));
    double beta_2020 = -5.0 + 10.0 / (1.0 + exp(-beta_2020_raw));
    double beta_2021 = -5.0 + 10.0 / (1.0 + exp(-beta_2021_raw));
    double beta_2022 = -5.0 + 10.0 / (1.0 + exp(-beta_2022_raw));
    double beta_2023 = -5.0 + 10.0 / (1.0 + exp(-beta_2023_raw));
    double beta_2024 = -5.0 + 10.0 / (1.0 + exp(-beta_2024_raw));

    // Triggering parameters (constant α, no mark dependence) (RELAXED BOUNDS):
    double beta_0_trig = -15.0 + 20.0 / (1.0 + exp(-beta_0_trig_raw));    // [-15, 5] (was [-10, 10])
    double decay = 0.01 + 49.99 / (1.0 + exp(-decay_raw));                 // [0.01, 50] (was [0.001, 10])

    // Constant alpha (exp(beta_0_trig) - same for all events)
    double alpha = exp(beta_0_trig);

    // Initialize log-likelihood
    double loglik = 0.0;

    // ===== Part 1: Sum of log intensities at event times =====
    for (int i = 0; i < n; i++) {
        // Lookup year effect based on discrete year (2015 is baseline)
        double beta_year = 0.0;  // 2015 baseline
        if (year[i] == 2016) beta_year = beta_2016;
        else if (year[i] == 2017) beta_year = beta_2017;
        else if (year[i] == 2018) beta_year = beta_2018;
        else if (year[i] == 2019) beta_year = beta_2019;
        else if (year[i] == 2020) beta_year = beta_2020;
        else if (year[i] == 2021) beta_year = beta_2021;
        else if (year[i] == 2022) beta_year = beta_2022;
        else if (year[i] == 2023) beta_year = beta_2023;
        else if (year[i] == 2024) beta_year = beta_2024;

        // Background rate with year fixed effects
        double mu = exp(beta_0_bg + gamma * log_pop[i] +
                       delta * poverty_decimal[i] + beta_year);

        // Triggering intensity from past events (CONSTANT ALPHA - no mark dependence)
        double trigger_sum = 0.0;
        for (int j = 0; j < i; j++) {
            double dt = times[i] - times[j];
            if (dt <= temporal_cutoff && dt > 0) {
                // Alpha is constant (no mark effects)
                // Kernel
                double kernel = exp(-decay * dt);
                trigger_sum += alpha * kernel;
            }
        }

        // Total intensity
        double lambda = mu + trigger_sum;

        // Add to log-likelihood (with safeguard)
        if (lambda > 0) {
            loglik += log(lambda);
        } else {
            // Penalty for invalid parameters
            return 1e10;
        }
    }

    // ===== Part 2: Compensator (integral term) =====

    // Background compensator: sum over district-years
    double background_compensator = 0.0;
    for (int k = 0; k < n_dy; k++) {
        // Lookup year effect based on discrete year (2015 is baseline)
        double beta_year_k = 0.0;  // 2015 baseline
        if (district_year_year[k] == 2016) beta_year_k = beta_2016;
        else if (district_year_year[k] == 2017) beta_year_k = beta_2017;
        else if (district_year_year[k] == 2018) beta_year_k = beta_2018;
        else if (district_year_year[k] == 2019) beta_year_k = beta_2019;
        else if (district_year_year[k] == 2020) beta_year_k = beta_2020;
        else if (district_year_year[k] == 2021) beta_year_k = beta_2021;
        else if (district_year_year[k] == 2022) beta_year_k = beta_2022;
        else if (district_year_year[k] == 2023) beta_year_k = beta_2023;
        else if (district_year_year[k] == 2024) beta_year_k = beta_2024;

        // Background rate with year fixed effects
        double mu_k = exp(beta_0_bg + gamma * district_year_log_pop[k] +
                         delta * district_year_poverty[k] + beta_year_k);

        background_compensator += mu_k * district_year_exposure[k];
    }

    // Triggering compensator: sum over all event pairs (CONSTANT ALPHA)
    double trigger_compensator = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double dt = times[i] - times[j];
            if (dt <= temporal_cutoff && dt > 0) {
                // Alpha is constant (no mark dependence)

                // Integral of exp(-decay * s) from 0 to dt
                double integral = (1.0 - exp(-decay * dt)) / decay;

                trigger_compensator += alpha * integral;
            }
        }
    }

    // Total compensator
    double compensator = background_compensator + trigger_compensator;

    // Log-likelihood = sum(log lambda) - compensator
    double neg_loglik = -(loglik - compensator);

    return neg_loglik;
}
