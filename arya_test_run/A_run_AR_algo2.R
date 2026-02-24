# ==============================================================================
# A_run_AR_algo2.R
# Runner Script for Algorithm 2 (Fast Version) on ALO Dataset
# ==============================================================================

rm(list=ls())
set.seed(1)  # For replication

# Set directories
code_dir <- '/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run'
data_dir <- "/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run/data/ALO_data.csv"
setwd(data_dir)

# Load packages
library(pbapply)
library(ivreg)
library(lmtest)
library(sandwich)

# Load supporting functions
source('2_gen_data.R')  # For gen_assignment_CR_index

# Load Algorithm 2
source(file.path(code_dir, 'B_AR_algo2.R'))

cat("\n", strrep("=", 70), "\n")
cat("ALGORITHM 2: FAST AR PERMUTATION TEST ON ALO DATA\n")
cat(strrep("=", 70), "\n\n")

# ==============================================================================
# STEP 1: Load and prepare data
# ==============================================================================
cat("Loading ALO dataset...\n")
data_dir <- "/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run/data/ALO_data.csv"
data <- read.csv(data_dir, header=TRUE, sep=",")

# Filter to males
data_m <- data[which(data$sex=="M"), ]

# Create SSP treatment/control groups
data_ssp_m <- data_m[which(data_m$control==1 | data_m$ssp==1), ]

cat("  Total observations:", nrow(data), "\n")
cat("  Male observations:", nrow(data_m), "\n")
cat("  SSP sample:", nrow(data_ssp_m), "\n\n")

# ==============================================================================
# STEP 2: Create dummy variables for categorical covariate
# ==============================================================================

cat("Creating dummy variables for 'lastmin' covariate...\n")
data_ssp_m$lastminnever <- as.numeric(data_ssp_m$lastmin=='never')
data_ssp_m$lastminusual <- as.numeric(data_ssp_m$lastmin=='usually')
data_ssp_m$lastminoccasional <- as.numeric(data_ssp_m$lastmin=='occasionally')
data_ssp_m$lastminoften <- as.numeric(data_ssp_m$lastmin=='often')
data_ssp_m$lastminrarely <- as.numeric(data_ssp_m$lastmin=='rarely')

# ==============================================================================
# STEP 3: Prepare data table for Algorithm 2
# ==============================================================================

cat("Preparing data table...\n")

data_ssp_m1 <- data_ssp_m[, c('GPA_year1', 'ssp_p', 'ssp', 'gpa0', 
                              'lastminnever', 'lastminusual', 'lastminoccasional', 
                              'lastminoften', 'lastminrarely')]

data_ssp_m1 <- data_ssp_m1[!is.na(data_ssp_m1$GPA_year1), ]

colnames(data_ssp_m1) <- c('Y_observed', 'D_observed', 'assignment', 
                           'x1', 'x2', 'x3', 'x4', 'x5', 'x6')

# Center covariates
data_ssp_m1[,'x1'] <- data_ssp_m1[,'x1'] - mean(data_ssp_m1[,'x1'])
data_ssp_m1[,'x2'] <- data_ssp_m1[,'x2'] - mean(data_ssp_m1[,'x2'])
data_ssp_m1[,'x3'] <- data_ssp_m1[,'x3'] - mean(data_ssp_m1[,'x3'])
data_ssp_m1[,'x4'] <- data_ssp_m1[,'x4'] - mean(data_ssp_m1[,'x4'])
data_ssp_m1[,'x5'] <- data_ssp_m1[,'x5'] - mean(data_ssp_m1[,'x5'])
data_ssp_m1[,'x6'] <- data_ssp_m1[,'x6'] - mean(data_ssp_m1[,'x6'])

N1 <- sum(data_ssp_m1$assignment)
N0 <- nrow(data_ssp_m1) - N1

cat("  Final sample size:", nrow(data_ssp_m1), "\n")
cat("  Treated (N1):", N1, "\n")
cat("  Control (N0):", N0, "\n")
cat("  Compliance rate:", 
    sum(data_ssp_m1$D_observed * data_ssp_m1$assignment) / sum(data_ssp_m1$assignment), "\n\n")

# ==============================================================================
# STEP 4: Generate 1000 random permutations
# ==============================================================================

cat("Generating 1000 random permutations...\n")
zsim <- pbsapply(1:1000, gen_assignment_CR_index, N1=N1, N0=N0)
cat("  ✓ Permutations generated\n\n")

# ==============================================================================
# STEP 5: Run Algorithm 2 (Fast Version)
# ==============================================================================

cat("Running Algorithm 2 (this should be much faster!)...\n\n")

start_time <- proc.time()

AR_CI_algo2 <- AR_algo2_custom(
  data_table = data_ssp_m1,
  N1 = N1,
  N0 = N0,
  zsim = zsim,
  tol = 1e-8,
  alpha = 0.95
)

end_time <- proc.time()
elapsed <- end_time - start_time

cat("\nAlgorithm 2 completed in", round(elapsed[3]/60, 2), "minutes\n\n")

# ==============================================================================
# STEP 6: Compare with ivreg
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("COMPARISON WITH IVREG\n")
cat(strrep("=", 70), "\n\n")

fit_ivreg <- ivreg(GPA_year1 ~ ssp_p + gpa0 + lastmin | ssp + gpa0 + lastmin, 
                   data = data_ssp_m)

ivreg_coef <- coef(fit_ivreg)[2]
ivreg_ci <- confint(fit_ivreg, vcov = vcovHC(fit_ivreg, type = "HC1"))[2, ]

cat("IVREG Results (2SLS with robust SE):\n")
cat("  Point estimate:", round(ivreg_coef, 4), "\n")
cat("  95% CI: [", round(ivreg_ci[1], 4), ", ", round(ivreg_ci[2], 4), "]\n\n", sep="")

cat("Algorithm 2 Results (Randomization-based):\n")
if (length(AR_CI_algo2) == 0) {
  cat("  95% CI: EMPTY SET\n\n")
} else {
  cat("  95% Confidence Set:\n")
  for (i in 1:length(AR_CI_algo2)) {
    cat("    [", round(AR_CI_algo2[[i]][1], 4), ", ", 
        round(AR_CI_algo2[[i]][2], 4), "]\n", sep="")
  }
  cat("\n")
}

# ==============================================================================
# STEP 7: Save results
# ==============================================================================

cat("Saving results...\n")

results_algo2 <- list(
  AR_confidence_set = AR_CI_algo2,
  ivreg_estimate = ivreg_coef,
  ivreg_ci = ivreg_ci,
  sample_size = nrow(data_ssp_m1),
  N1 = N1,
  N0 = N0,
  n_permutations = 1000,
  elapsed_time_minutes = elapsed[3]/60,
  algorithm = "Algorithm 2 (Fast Jumping)"
)

save(results_algo2, file=file.path(code_dir, 'algo2_ALO_results.Rdata'))
cat("  ✓ Results saved to 'algo2_ALO_results.Rdata'\n\n")

cat(strrep("=", 70), "\n")
cat("DONE! Algorithm 2 is much faster than Algorithm 1!\n")
cat(strrep("=", 70), "\n")