# Load packages
library(pbapply)
library(ivreg)
library(lmtest)
library(sandwich)
source('/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run/B_solve_coef.R')
source('/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run/B_calculate_intersections.R')
source('/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run/B_find_intervals.R')
source('/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run/B_AR_algo2.R')


gen_data_dose_response <- function(N, rho = 0.95, gamma = 0) {
  # Generate standard normal variables
  W1i <- rnorm(N, 0, 1)
  W2i <- rnorm(N, 0, 1)
  Xi <- rnorm(N, 0, 1)
  
  # Generate potential doses
  W_tilde_1i <- abs(Xi) * W1i  # |Xi|W1i
  d0 <- W_tilde_1i
  d1 <- gamma + d0
  
  # Generate potential outcomes (NO TREATMENT EFFECT)
  # Yi(di(1)) = Yi(di(0)) = ρW1i + sqrt(1-ρ²)W2i
  y0 <- rho * W1i + sqrt(1 - rho^2) * W2i
  y1 <- y0  # No treatment effect!
  
  # Create output table
  output <- data.frame(
    indices = 1:N,
    y1 = y1,
    y0 = y0,
    d1 = d1,
    d0 = d0
  )
  
  return(output)
}

# Parameters
N <- 100
N1 <- 50  # Treated
N0 <- 50  # Control
rho <- 0.95
gamma <- 0
n_permutations <- 1000
n_runs <- 4  # Run 4 times to check variability

# Storage for results
runtime_results <- data.frame(
  run = 1:n_runs,
  runtime_seconds = NA,
  CI_lower = NA,
  CI_upper = NA
)

for (run in 1:n_runs) {
  
  cat("\n=== RUN", run, "===\n")
  
  # Generate data
  potential_outcomes <- gen_data_dose_response(N, rho, gamma)
  
  # Random assignment
  set.seed(123 + run)  # Different seed each run
  assignment <- c(rep(1, N1), rep(0, N0))
  assignment <- sample(assignment)  # Randomize
  
  # Observed data
  y_obs <- ifelse(assignment == 1, 
                  potential_outcomes$y1, 
                  potential_outcomes$y0)
  d_obs <- ifelse(assignment == 1, 
                  potential_outcomes$d1, 
                  potential_outcomes$d0)
  
  # Create data table for AR algorithm
  data_table <- data.frame(
    Y_observed = y_obs,
    D_observed = d_obs,
    assignment = assignment
  )
  
  # Generate permutations
  zsim <- pbsapply(1:n_permutations, gen_assignment_CR_index, 
                   N1 = N1, N0 = N0)
  
  # TIME THE ALGORITHM
  start_time <- proc.time()
  
  CI <- AR_algo2_custom(
    data_table = data_table,
    N1 = N1,
    N0 = N0,
    zsim = zsim,
    tol = 1e-8,
    alpha = 0.95
  )
  
  elapsed <- proc.time() - start_time
  
  # Store results
  runtime_results$runtime_seconds[run] <- elapsed[3]
  runtime_results$CI_lower[run] <- CI[[1]][1]
  runtime_results$CI_upper[run] <- CI[[1]][2]
  
  cat("Runtime:", round(elapsed[3], 2), "seconds\n")
  cat("95% CI: [", round(CI[[1]][1], 4), ", ", 
      round(CI[[1]][2], 4), "]\n")
}

# Summary
cat("\n=== SUMMARY ACROSS 4 RUNS ===\n")
cat("Mean runtime:", round(mean(runtime_results$runtime_seconds), 2), "seconds\n")
cat("SD runtime:", round(sd(runtime_results$runtime_seconds), 2), "seconds\n")
cat("Min runtime:", round(min(runtime_results$runtime_seconds), 2), "seconds\n")
cat("Max runtime:", round(max(runtime_results$runtime_seconds), 2), "seconds\n")

print(runtime_results)

