# ==============================================================================
# Simulation: Dose-Response Design with Treatment Effects
# Testing Coverage Rates for Different Instrument Strengths and Treatment Effects
# ==============================================================================

library(pbapply)

# Set directories and source functions
code_dir <- '/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run'
setwd(code_dir)

source(file.path(code_dir, 'B_gen_data.R'))
source(file.path(code_dir, 'B_AR_algo2.R'))

# ==============================================================================
# Data Generating Function
# ==============================================================================

gen_dose_response_data <- function(N, rho = 0.95, gamma = 0, tau = 0) {
  # Generate standard normal variables
  W1i <- rnorm(N, 0, 1)
  W2i <- rnorm(N, 0, 1)
  Xi <- rnorm(N, 0, 1)
  
  # Potential doses
  W_tilde_1i <- abs(Xi) * W1i  # |Xi|W1i
  d0 <- W_tilde_1i
  d1 <- gamma + d0
  
  # Potential outcomes
  # Control: Yi(di(0)) = ρW1i + sqrt(1-ρ²)W2i
  y0 <- rho * W1i + sqrt(1 - rho^2) * W2i
  
  # Treatment: Yi(di(1)) = Yi(di(0)) + tau * (d1 - d0)
  # FIXED: Effect is proportional to dose difference, not constant!
  y1 <- y0 + tau * (d1 - d0)
  
  return(data.frame(
    y1 = y1,
    y0 = y0,
    d1 = d1,
    d0 = d0
  ))
}

# ==============================================================================
# Single Simulation Run
# ==============================================================================

run_simulation <- function(N = 100, N1 = 50, rho = 0.95, 
                           gamma = 0, tau = 0, n_perms = 1000) {
  # Generate potential outcomes
  pot_outcomes <- gen_dose_response_data(N, rho, gamma, tau)
  
  # Random assignment
  assignment <- c(rep(1, N1), rep(0, N - N1))
  assignment <- sample(assignment)
  
  # Observed data
  y_obs <- ifelse(assignment == 1, pot_outcomes$y1, pot_outcomes$y0)
  d_obs <- ifelse(assignment == 1, pot_outcomes$d1, pot_outcomes$d0)
  
  data_table <- data.frame(
    Y_observed = y_obs,
    D_observed = d_obs,
    assignment = assignment
  )
  
  # Generate permutations
  zsim <- pbsapply(1:n_perms, gen_assignment_CR_index, N1 = N1, N0 = N0)
  
  # Run AR algorithm (time it)
  start_time <- proc.time()
  CI <- AR_algo2_custom(data_table, N1, N - N1, zsim, alpha = 0.95)
  runtime <- (proc.time() - start_time)[3]
  
  # Check if CI contains true value (tau)
  # Handle infinite CIs and ensure proper closed interval check [a,b]
  is_infinite_CI <- is.infinite(CI[[1]][1]) || is.infinite(CI[[1]][2])
  
  covers_true <- if (is_infinite_CI) {
    TRUE  # Infinite CI always covers any finite value
  } else {
    # Closed interval [lower, upper]: both endpoints included
    # Check: lower <= tau AND tau <= upper
    (CI[[1]][1] <= tau) && (tau <= CI[[1]][2])
  }
  
  return(list(
    covers = covers_true,
    CI_lower = CI[[1]][1],
    CI_upper = CI[[1]][2],
    CI_width = CI[[1]][2] - CI[[1]][1],
    runtime = runtime,
    is_infinite = is_infinite_CI
  ))
}

# ==============================================================================
# Run Full Simulation Grid
# ==============================================================================

cat("\n", strrep("=", 80), "\n")
cat("DOSE-RESPONSE SIMULATION: GAMMA x TAU GRID\n")
cat(strrep("=", 80), "\n\n")

# Parameters
gamma_values <- c(0, 0.5, 1)      # Instrument strength
tau_values <- c(0, 0.5, 1)        # Treatment effects
n_sims <- 100                      # Number of simulations per cell
n_perms <- 1000                    # Permutations per simulation

# Storage for results
results_matrix <- matrix(NA, nrow = length(gamma_values), ncol = length(tau_values))
rownames(results_matrix) <- paste("gamma =", gamma_values)
colnames(results_matrix) <- paste("tau =", tau_values)

coverage_matrix <- results_matrix
width_matrix <- results_matrix
runtime_matrix <- results_matrix
infinite_matrix <- results_matrix  # Track % of infinite CIs

# Set seed for reproducibility
set.seed(12345)

# Main simulation loop
total_cells <- length(gamma_values) * length(tau_values)
current_cell <- 0

for (i in 1:length(gamma_values)) {
  for (j in 1:length(tau_values)) {
    current_cell <- current_cell + 1
    gamma <- gamma_values[i]
    tau <- tau_values[j]
    
    cat("\n", strrep("-", 80), "\n")
    cat("Cell", current_cell, "of", total_cells, 
        "| gamma =", gamma, "| tau =", tau, "\n")
    cat(strrep("-", 80), "\n")
    
    # Run simulations for this cell
    sim_results <- replicate(n_sims, {
      run_simulation(N = 100, N1 = 50, gamma = gamma, tau = tau, n_perms = n_perms)
    }, simplify = FALSE)
    
    # Extract results
    coverage <- mean(sapply(sim_results, function(x) x$covers))
    mean_width <- mean(sapply(sim_results, function(x) x$CI_width))
    mean_runtime <- mean(sapply(sim_results, function(x) x$runtime))
    pct_infinite <- mean(sapply(sim_results, function(x) x$is_infinite))
    
    # Store in matrices
    coverage_matrix[i, j] <- coverage
    width_matrix[i, j] <- mean_width
    runtime_matrix[i, j] <- mean_runtime
    infinite_matrix[i, j] <- pct_infinite
    
    cat("Coverage rate:", round(coverage * 100, 1), "%\n")
    cat("Mean CI width:", round(mean_width, 3), "\n")
    cat("Mean runtime:", round(mean_runtime, 2), "seconds\n")
    cat("% Infinite CIs:", round(pct_infinite * 100, 1), "%\n")
  }
}

# ==============================================================================
# Create Results Tables
# ==============================================================================

cat("\n\n", strrep("=", 80), "\n")
cat("FINAL RESULTS: COVERAGE RATES (%)\n")
cat(strrep("=", 80), "\n\n")

coverage_df <- as.data.frame(round(coverage_matrix * 100, 1))
print(coverage_df)

cat("\n\n", strrep("=", 80), "\n")
cat("FINAL RESULTS: MEAN CI WIDTH\n")
cat(strrep("=", 80), "\n\n")

width_df <- as.data.frame(round(width_matrix, 3))
print(width_df)

cat("\n\n", strrep("=", 80), "\n")
cat("FINAL RESULTS: MEAN RUNTIME (SECONDS)\n")
cat(strrep("=", 80), "\n\n")

runtime_df <- as.data.frame(round(runtime_matrix, 2))
print(runtime_df)

cat("\n\n", strrep("=", 80), "\n")
cat("FINAL RESULTS: % INFINITE CIs\n")
cat(strrep("=", 80), "\n\n")

infinite_df <- as.data.frame(round(infinite_matrix * 100, 1))
print(infinite_df)

# ==============================================================================
# Create Combined Table for Vignette
# ==============================================================================

cat("\n\n", strrep("=", 80), "\n")
cat("COMBINED RESULTS TABLE (FOR VIGNETTE)\n")
cat(strrep("=", 80), "\n\n")

# Create long-format table
combined_results <- data.frame(
  Gamma = rep(gamma_values, each = length(tau_values)),
  Tau = rep(tau_values, times = length(gamma_values)),
  Coverage_Pct = as.vector(t(coverage_matrix)) * 100,
  Mean_Width = as.vector(t(width_matrix)),
  Mean_Runtime_Sec = as.vector(t(runtime_matrix)),
  Pct_Infinite = as.vector(t(infinite_matrix)) * 100
)

combined_results$Coverage_Pct <- round(combined_results$Coverage_Pct, 1)
combined_results$Mean_Width <- round(combined_results$Mean_Width, 3)
combined_results$Mean_Runtime_Sec <- round(combined_results$Mean_Runtime_Sec, 2)
combined_results$Pct_Infinite <- round(combined_results$Pct_Infinite, 1)

print(combined_results)

# ==============================================================================
# Save Results
# ==============================================================================

cat("\n\nSaving results...\n")

simulation_results <- list(
  coverage_matrix = coverage_matrix,
  width_matrix = width_matrix,
  runtime_matrix = runtime_matrix,
  infinite_matrix = infinite_matrix,
  combined_table = combined_results,
  parameters = list(
    gamma_values = gamma_values,
    tau_values = tau_values,
    n_sims = n_sims,
    n_perms = n_perms,
    N = 100,
    N1 = 50,
    rho = 0.95
  )
)

save(simulation_results, file = file.path(code_dir, 'dose_response_simulation_results.Rdata'))
cat("Results saved to 'dose_response_simulation_results.Rdata'\n")

# ==============================================================================
# Summary Statistics
# ==============================================================================

cat("\n\n", strrep("=", 80), "\n")
cat("SUMMARY STATISTICS\n")
cat(strrep("=", 80), "\n\n")

cat("Coverage rates:\n")
cat("  Min:", round(min(coverage_matrix) * 100, 1), "%\n")
cat("  Max:", round(max(coverage_matrix) * 100, 1), "%\n")
cat("  Mean:", round(mean(coverage_matrix) * 100, 1), "%\n")
cat("  Target: 95.0%\n\n")

cat("CI widths:\n")
cat("  Min:", round(min(width_matrix), 3), "\n")
cat("  Max:", round(max(width_matrix), 3), "\n")
cat("  Mean:", round(mean(width_matrix), 3), "\n\n")

cat("Runtimes:\n")
cat("  Min:", round(min(runtime_matrix), 2), "seconds\n")
cat("  Max:", round(max(runtime_matrix), 2), "seconds\n")
cat("  Mean:", round(mean(runtime_matrix), 2), "seconds\n\n")

cat("% Infinite CIs:\n")
cat("  Min:", round(min(infinite_matrix) * 100, 1), "%\n")
cat("  Max:", round(max(infinite_matrix) * 100, 1), "%\n")
cat("  Mean:", round(mean(infinite_matrix) * 100, 1), "%\n\n")

# ==============================================================================
# Interpretation Guide
# ==============================================================================

cat(strrep("=", 80), "\n")
cat("INTERPRETATION GUIDE\n")
cat(strrep("=", 80), "\n\n")

cat("Expected patterns:\n")
cat("  • gamma=0 (no IV): ~100% infinite CIs → coverage should be ~100%\n")
cat("  • gamma=0.5 (weak IV): High % infinite CIs → coverage ~95-100%\n")
cat("  • gamma=1 (strong IV): ~0% infinite CIs → coverage should be ~95%\n\n")

cat("Coverage interpretation:\n")
cat("  • When most CIs are infinite: coverage ≈ 100% (correct!)\n")
cat("  • When most CIs are finite: coverage ≈ 95% (correct!)\n")
cat("  • Low coverage with finite CIs: potential problem\n\n")

cat(strrep("=", 80), "\n")
cat("SIMULATION COMPLETE!\n")
cat(strrep("=", 80), "\n")
