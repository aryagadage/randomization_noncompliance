# ==============================================================================
# compare_algo1_vs_algo2.R
# Quick Comparison: Algorithm 1 vs Algorithm 2
# ==============================================================================
#
# This script runs BOTH algorithms with only 100 permutations
# to quickly compare their speed and verify they give similar results
#
# ==============================================================================

rm(list=ls())
set.seed(321)

# Set directories
code_dir <- '/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run'
data_dir <- "/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run/data/ALO_data.csv"

setwd(code_dir)

# Load packages
library(pbapply)
library(ivreg)

# Load supporting functions
source('B_gen_data.R')  # For gen_assignment_CR_index

cat("\n", strrep("=", 70), "\n")
cat("QUICK COMPARISON: ALGORITHM 1 vs ALGORITHM 2\n")
cat(strrep("=", 70), "\n\n")

# ==============================================================================
# Prepare Data (same for both)
# ==============================================================================

# From this:
data_file <- file.path(data_dir, "..", "ALO_data.csv")

# To this (full path to the CSV file):
data_file <- "/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run/data/ALO_data.csv"
data <- read.csv(data_file, header=TRUE, sep=",")

data_m <- data[which(data$sex=="M"), ]
data_ssp_m <- data_m[which(data_m$control==1 | data_m$ssp==1), ]

# Create dummy variables
data_ssp_m$lastminnever <- as.numeric(data_ssp_m$lastmin=='never')
data_ssp_m$lastminusual <- as.numeric(data_ssp_m$lastmin=='usually')
data_ssp_m$lastminoccasional <- as.numeric(data_ssp_m$lastmin=='occasionally')
data_ssp_m$lastminoften <- as.numeric(data_ssp_m$lastmin=='often')
data_ssp_m$lastminrarely <- as.numeric(data_ssp_m$lastmin=='rarely')

# Prepare data table
data_ssp_m1 <- data_ssp_m[, c('GPA_year1', 'ssp_p', 'ssp', 'gpa0', 
                              'lastminnever', 'lastminusual', 'lastminoccasional', 
                              'lastminoften', 'lastminrarely')]
data_ssp_m1 <- data_ssp_m1[!is.na(data_ssp_m1$GPA_year1), ]

colnames(data_ssp_m1) <- c('Y_observed', 'D_observed', 'assignment', 
                           'x1', 'x2', 'x3', 'x4', 'x5', 'x6')

# Center covariates
for (i in 4:9) {
  data_ssp_m1[, i] <- data_ssp_m1[, i] - mean(data_ssp_m1[, i])
}

N1 <- sum(data_ssp_m1$assignment)
N0 <- nrow(data_ssp_m1) - N1

cat("  Sample size:", nrow(data_ssp_m1), "\n")
cat("  Treated (N1):", N1, "\n")
cat("  Control (N0):", N0, "\n\n")

# ==============================================================================
# Generate 1000 permutations (for real comparison)
# ==============================================================================

cat("Generating 1000 random permutations...\n")
zsim <- pbsapply(1:1000, gen_assignment_CR_index, N1=N1, N0=N0)
cat("  Permutations generated\n\n")

# ==============================================================================
# RUN ALGORITHM 1
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("RUNNING ALGORITHM 1 (Exhaustive - checks all intervals)\n")
cat(strrep("=", 70), "\n\n")

source(file.path(code_dir, 'B_AR_algo1.R'))

start_algo1 <- proc.time()

CI_algo1 <- AR_algo1_custom(
  data_table = data_ssp_m1,
  N1 = N1,
  N0 = N0,
  zsim = zsim,
  tol = 1e-8,
  alpha = 0.95
)

time_algo1 <- proc.time() - start_algo1

cat("\nAlgorithm 1 completed in", round(time_algo1[3]/60, 2), "minutes\n\n")

# ==============================================================================
# RUN ALGORITHM 2
# ==============================================================================

cat(strrep("=", 70), "\n")
cat("RUNNING ALGORITHM 2 (Fast - jumps between crossings)\n")
cat(strrep("=", 70), "\n\n")

source(file.path(code_dir, 'B_AR_algo2.R'))

start_algo2 <- proc.time()

CI_algo2 <- AR_algo2_custom(
  data_table = data_ssp_m1,
  N1 = N1,
  N0 = N0,
  zsim = zsim,
  tol = 1e-8,
  alpha = 0.95
)

time_algo2 <- proc.time() - start_algo2

cat("\nAlgorithm 2 completed in", round(time_algo2[3]/60, 2), "minutes\n\n")

# ==============================================================================
# COMPARISON
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("COMPARISON RESULTS\n")
cat(strrep("=", 70), "\n\n")

cat("Time Comparison:\n")
cat("  Algorithm 1: ", round(time_algo1[3]/60, 2), " minutes\n", sep="")
cat("  Algorithm 2: ", round(time_algo2[3]/60, 2), " minutes\n", sep="")
cat("  Speedup: ", round(time_algo1[3] / time_algo2[3], 1), "x faster\n\n", sep="")

cat("Confidence Intervals:\n\n")

cat("Algorithm 1 (95% CI):\n")
if (length(CI_algo1) == 0) {
  cat("  EMPTY SET\n")
} else {
  for (i in 1:length(CI_algo1)) {
    cat("  [", round(CI_algo1[[i]][1], 4), ", ", round(CI_algo1[[i]][2], 4), "]\n", sep="")
  }
}

cat("\nAlgorithm 2 (95% CI):\n")
if (length(CI_algo2) == 0) {
  cat("  EMPTY SET\n")
} else {
  for (i in 1:length(CI_algo2)) {
    cat("  [", round(CI_algo2[[i]][1], 4), ", ", round(CI_algo2[[i]][2], 4), "]\n", sep="")
  }
}

cat("\n")

# Check if results are similar
if (length(CI_algo1) == length(CI_algo2)) {
  if (length(CI_algo1) == 0) {
    cat("✓ Both algorithms returned EMPTY SET (match!)\n")
  } else {
    # Compare intervals
    all_close <- TRUE
    for (i in 1:length(CI_algo1)) {
      diff_left <- abs(CI_algo1[[i]][1] - CI_algo2[[i]][1])
      diff_right <- abs(CI_algo1[[i]][2] - CI_algo2[[i]][2])
      
      if (diff_left > 0.01 || diff_right > 0.01) {
        all_close <- FALSE
        break
      }
    }
    
    if (all_close) {
      cat("Both algorithms produced SIMILAR confidence intervals!\n")
    } else {
      cat(" Intervals differ slightly \n")
    }
  }
} else {
  cat(" Different number of intervals \n")
}

cat("\n", strrep("=", 70), "\n")
cat("CONCLUSION\n")
cat(strrep("=", 70), "\n\n")

cat("With", ncol(zsim), "permutations:\n")
cat("  Algorithm 2 is ", round(time_algo1[3] / time_algo2[3], 1), "x faster\n\n", sep="")

# ==============================================================================
# Save comparison results
# ==============================================================================

comparison_results <- list(
  n_permutations = ncol(zsim),
  algo1_time_seconds = time_algo1[3],
  algo2_time_seconds = time_algo2[3],
  speedup = time_algo1[3] / time_algo2[3],
  algo1_CI = CI_algo1,
  algo2_CI = CI_algo2
)

save(comparison_results, file=file.path(code_dir, 'algo_comparison_results.Rdata'))

# ==============================================================================
# Display Results Table
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("SUMMARY TABLE\n")
cat(strrep("=", 70), "\n\n")

# Create comparison table
comparison_table <- data.frame(
  Metric = c(
    "Permutations",
    "Algo 1 Time",
    "Algo 2 Time",
    "Speedup",
    "Algo 1 CI",
    "Algo 2 CI"
  ),
  Value = c(
    ncol(zsim),
    paste(round(time_algo1[3] / 60, 2), "min"),
    paste(round(time_algo2[3] / 60, 2), "min"),
    paste(round(time_algo1[3] / time_algo2[3], 1), "x faster"),
    ifelse(length(CI_algo1) > 0, 
           paste0("[", round(CI_algo1[[1]][1], 3), ", ", round(CI_algo1[[1]][2], 3), "]"),
           "Empty"),
    ifelse(length(CI_algo2) > 0, 
           paste0("[", round(CI_algo2[[1]][1], 3), ", ", round(CI_algo2[[1]][2], 3), "]"),
           "Empty")
  )
)

print(comparison_table, row.names = FALSE)

cat("\n✓ Results saved to 'algo_comparison_results.Rdata'\n\n")