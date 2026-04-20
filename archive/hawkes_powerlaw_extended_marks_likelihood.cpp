#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double hawkes_powerlaw_extended_marks_negloglik_cpp(
    NumericVector times,
    NumericVector log_pop,
    NumericVector poverty_decimal,
    IntegerVector year,
    IntegerVector has_fatalities,
    IntegerVector is_violent,
    IntegerVector is_student,
    IntegerVector is_labor,
    IntegerVector is_papua,
    NumericVector district_year_exposure,
    NumericVector district_year_log_pop,
    NumericVector district_year_poverty,
    IntegerVector district_year_year,
    NumericVector params,
    double temporal_cutoff
) {

    int n = times.size();
    int n_dy = district_year_exposure.size();

    // Extract parameters (20 total: 12 background + 8 triggering)
    // Background parameters (12)
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

    // Triggering parameters (8): baseline + 5 marks + 2 power-law kernel
    double beta_0_trig_raw = params[12];
    double beta_fatal_raw = params[13];
    double beta_violent_raw = params[14];        // Violence main effect
    double beta_student_raw = params[15];
    double beta_labor_raw = params[16];
    double beta_papua_raw = params[17];
    // Power-law kernel parameters
    double alpha_raw = params[18];
    double c_raw = params[19];

    // Transform raw (unconstrained) parameters to bounded parameters using sigmoid
    // Background rate parameters:
    double gamma = 0.01 + 1.99 / (1.0 + exp(-gamma_raw));         // [0.01, 2]
    double delta = -20.0 + 40.0 / (1.0 + exp(-delta_raw));        // [-20, 20]

    // Year effects (9 parameters):
    double beta_2016 = -3.0 + 6.0 / (1.0 + exp(-beta_2016_raw));  // [-3, 3]
    double beta_2017 = -3.0 + 6.0 / (1.0 + exp(-beta_2017_raw));
    double beta_2018 = -3.0 + 6.0 / (1.0 + exp(-beta_2018_raw));
    double beta_2019 = -3.0 + 6.0 / (1.0 + exp(-beta_2019_raw));
    double beta_2020 = -3.0 + 6.0 / (1.0 + exp(-beta_2020_raw));
    double beta_2021 = -3.0 + 6.0 / (1.0 + exp(-beta_2021_raw));
    double beta_2022 = -3.0 + 6.0 / (1.0 + exp(-beta_2022_raw));
    double beta_2023 = -3.0 + 6.0 / (1.0 + exp(-beta_2023_raw));
    double beta_2024 = -3.0 + 6.0 / (1.0 + exp(-beta_2024_raw));

    // Triggering parameters:
    double beta_0_trig = -10.0 + 20.0 / (1.0 + exp(-beta_0_trig_raw));           // [-10, 10]
    double beta_fatal = -5.0 + 10.0 / (1.0 + exp(-beta_fatal_raw));              // [-5, 5]
    double beta_violent = -5.0 + 10.0 / (1.0 + exp(-beta_violent_raw));          // [-5, 5]
    double beta_student = -5.0 + 10.0 / (1.0 + exp(-beta_student_raw));          // [-5, 5]
    double beta_labor = -5.0 + 10.0 / (1.0 + exp(-beta_labor_raw));              // [-5, 5]
    double beta_papua = -5.0 + 10.0 / (1.0 + exp(-beta_papua_raw));              // [-5, 5]

    // Power-law kernel parameters:
    double alpha = 1.05 + 2.45 / (1.0 + exp(-alpha_raw));                 // [1.05, 3.5]
    double c = 0.01 + 0.99 / (1.0 + exp(-c_raw));                         // [0.01, 1.0]

    // Safety checks for alpha and c
    if (alpha <= 1.0 || alpha >= 3.6) {
        return 1e10;
    }
    if (c <= 0.005 || c >= 1.5) {
        return 1e10;
    }

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

        // Triggering intensity from past events
        double trigger_sum = 0.0;
        for (int j = 0; j < i; j++) {
            double dt = times[i] - times[j];
            if (dt <= temporal_cutoff && dt > 0) {
                // Alpha effect depends on 5 marks of triggering event j
                // fatality, violence, 3 actor types
                double alpha_effect = exp(beta_0_trig +
                                  beta_fatal * has_fatalities[j] +
                                  beta_violent * is_violent[j] +
                                  beta_student * is_student[j] +
                                  beta_labor * is_labor[j] +
                                  beta_papua * is_papua[j]);

                // POWER-LAW KERNEL: g(dt) = (dt + c)^(-alpha)
                double kernel = pow(dt + c, -alpha);

                // Check for numerical issues
                if (isnan(kernel) || isinf(kernel)) {
                    return 1e10;
                }

                trigger_sum += alpha_effect * kernel;
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

    // Triggering compensator: sum over all event pairs
    double trigger_compensator = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            double dt = times[i] - times[j];
            if (dt <= temporal_cutoff && dt > 0) {
                // Alpha effect for event j (5 mark effects)
                double alpha_effect = exp(beta_0_trig +
                                  beta_fatal * has_fatalities[j] +
                                  beta_violent * is_violent[j] +
                                  beta_student * is_student[j] +
                                  beta_labor * is_labor[j] +
                                  beta_papua * is_papua[j]);

                // POWER-LAW INTEGRAL:
                // ∫₀^dt (s + c)^(-α) ds = [c^(1-α) - (dt+c)^(1-α)] / (α - 1)
                double c_power = pow(c, 1.0 - alpha);
                double dt_c_power = pow(dt + c, 1.0 - alpha);
                double integral = (c_power - dt_c_power) / (alpha - 1.0);

                // Check for numerical issues
                if (isnan(integral) || isinf(integral) || integral < 0) {
                    return 1e10;
                }

                trigger_compensator += alpha_effect * integral;
            }
        }
    }

    // Check for unreasonable compensator values
    if (trigger_compensator < 0 || trigger_compensator > 1e8) {
        return 1e10;
    }

    // Total compensator
    double compensator = background_compensator + trigger_compensator;

    // Log-likelihood = sum(log lambda) - compensator
    double neg_loglik = -(loglik - compensator);

    return neg_loglik;
}
