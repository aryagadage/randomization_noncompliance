# ==============================================================================
# run_algo1_on_ALO_data.R
# Runner Script for Algorithm 1 on ALO Dataset
# ==============================================================================
#
# This script:
# 1. Loads the ALO dataset
# 2. Prepares the data (SSP first year GPA analysis for males)
# 3. Generates 1000 random permutations
# 4. Runs Algorithm 1 using the custom implementation
# 5. Compares with ivreg results
#
# ==============================================================================

rm(list=ls())
set.seed(1)  # For replication

# Set working directory
dir <- '/Users/ag5276/Documents/Github/Rand_IV/grace'
setwd(dir)

# Load packages
library(pbapply)
library(ivreg)
library(lmtest)
library(sandwich)

# Load supporting functions
source('2_gen_data.R')  # For gen_assignment_CR_index

# Load Algorithm 1
source('AR_algo1_custom.R')

cat("="*70, "\n")
cat("ALGORITHM 1: AR PERMUTATION TEST ON ALO DATA\n")
cat("="*70, "\n\n")

# ==============================================================================
# STEP 1: Load and prepare data
# ==============================================================================

cat("Loading ALO dataset...\n")
data_dir <- "/Users/ag5276/Documents/ALO_data.csv"
data <- read.csv(data_dir, header=TRUE, sep=",")

# Filter to males
data_m <- data[which(data$sex=="M"), ]

# Create SSP treatment/control groups
data_ssp_m <- data_m[which(data_m$control==1 | data_m$ssp==1), ]

cat("  Total observations:", nrow(data), "\n")
cat("  Male observations:", nrow(data_m), "\n")
cat("  SSP sample:", nrow(data_ssp_m), "\n\n")

# ==============================================================================
# STEP 2: Create dummy variables for categorical covariate (lastmin)
# ==============================================================================

cat("Creating dummy variables for 'lastmin' covariate...\n")
data_ssp_m$lastminnever <- as.numeric(data_ssp_m$lastmin=='never')
data_ssp_m$lastminusual <- as.numeric(data_ssp_m$lastmin=='usually')
data_ssp_m$lastminoccasional <- as.numeric(data_ssp_m$lastmin=='occasionally')
data_ssp_m$lastminoften <- as.numeric(data_ssp_m$lastmin=='often')
data_ssp_m$lastminrarely <- as.numeric(data_ssp_m$lastmin=='rarely')

# ==============================================================================
# STEP 3: Prepare data table for Algorithm 1
# ==============================================================================

cat("Preparing data table...\n")

# Select columns: outcome, treatment, assignment, covariates
data_ssp_m1 <- data_ssp_m[, c('GPA_year1', 'ssp_p', 'ssp', 'gpa0', 
                              'lastminnever', 'lastminusual', 'lastminoccasional', 
                              'lastminoften', 'lastminrarely')]

# Remove missing values
data_ssp_m1 <- data_ssp_m1[!is.na(data_ssp_m1$GPA_year1), ]

# Rename columns to standard format
colnames(data_ssp_m1) <- c('Y_observed', 'D_observed', 'assignment', 
                           'x1', 'x2', 'x3', 'x4', 'x5', 'x6')

# Center covariates (important for algorithm)
data_ssp_m1[,'x1'] <- data_ssp_m1[,'x1'] - mean(data_ssp_m1[,'x1'])
data_ssp_m1[,'x2'] <- data_ssp_m1[,'x2'] - mean(data_ssp_m1[,'x2'])
data_ssp_m1[,'x3'] <- data_ssp_m1[,'x3'] - mean(data_ssp_m1[,'x3'])
data_ssp_m1[,'x4'] <- data_ssp_m1[,'x4'] - mean(data_ssp_m1[,'x4'])
data_ssp_m1[,'x5'] <- data_ssp_m1[,'x5'] - mean(data_ssp_m1[,'x5'])
data_ssp_m1[,'x6'] <- data_ssp_m1[,'x6'] - mean(data_ssp_m1[,'x6'])

# Calculate sample sizes
N1 <- sum(data_ssp_m1$assignment)  # Treated
N0 <- nrow(data_ssp_m1) - N1       # Control

cat("  Final sample size:", nrow(data_ssp_m1), "\n")
cat("  Treated (N1):", N1, "\n")
cat("  Control (N0):", N0, "\n")
cat("  Compliance rate:", sum(data_ssp_m1$D_observed * data_ssp_m1$assignment) / sum(data_ssp_m1$assignment), "\n\n")

# ==============================================================================
# STEP 4: Generate 1000 random permutations
# ==============================================================================

cat("Generating 1000 random permutations...\n")
zsim <- pbsapply(1:1000, gen_assignment_CR_index, N1=N1, N0=N0)
cat("  ✓ Permutations generated\n\n")

# ==============================================================================
# STEP 5: Run Algorithm 1
# ==============================================================================

cat("Running Algorithm 1 (this may take several minutes)...\n\n")

start_time <- proc.time()

# Run Algorithm 1 with 95% confidence level
AR_CI <- AR_algo1_custom(
  data_table = data_ssp_m1,
  N1 = N1,
  N0 = N0,
  zsim = zsim,
  tol = 1e-8,
  alpha = 0.95
)

end_time <- proc.time()
elapsed <- end_time - start_time

cat("\nAlgorithm 1 completed in", elapsed[3], "seconds\n\n")

# ==============================================================================
# STEP 6: Compare with ivreg (for reference)
# ==============================================================================

cat("="*70, "\n")
cat("COMPARISON WITH IVREG\n")
cat("="*70, "\n\n")

# Run ivreg on same data
fit_ivreg <- ivreg(GPA_year1 ~ ssp_p + gpa0 + lastmin | ssp + gpa0 + lastmin, 
                   data = data_ssp_m)

# Get point estimate and 95% CI
ivreg_coef <- coef(fit_ivreg)[2]
ivreg_ci <- confint(fit_ivreg, vcov = vcovHC(fit_ivreg, type = "HC1"))[2, ]

cat("IVREG Results (2SLS with robust SE):\n")
cat("  Point estimate:", round(ivreg_coef, 4), "\n")
cat("  95% CI: [", round(ivreg_ci[1], 4), ", ", round(ivreg_ci[2], 4), "]\n\n", sep="")

cat("Algorithm 1 Results (Randomization-based):\n")
if (length(AR_CI) == 0) {
  cat("  95% CI: EMPTY SET\n\n")
} else {
  cat("  95% Confidence Set:\n")
  for (i in 1:length(AR_CI)) {
    cat("    [", round(AR_CI[[i]][1], 4), ", ", round(AR_CI[[i]][2], 4), "]\n", sep="")
  }
  cat("\n")
}

# ==============================================================================
# STEP 7: Save results
# ==============================================================================

cat("Saving results...\n")

results <- list(
  AR_confidence_set = AR_CI,
  ivreg_estimate = ivreg_coef,
  ivreg_ci = ivreg_ci,
  sample_size = nrow(data_ssp_m1),
  N1 = N1,
  N0 = N0,
  n_permutations = 1000,
  elapsed_time = elapsed[3]
)

save(results, file='algo1_ALO_results.Rdata')
cat("  ✓ Results saved to 'algo1_ALO_results.Rdata'\n\n")

cat("="*70, "\n")
cat("DONE!\n")
cat("="*70, "\n")
