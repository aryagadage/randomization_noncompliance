# ==============================================================================
# compare_1000_vs_10000_permutations.R
# Compare Algorithm 2 with 1000 vs 10000 Permutations
# ==============================================================================

rm(list=ls())
set.seed(662)

# Set directories
code_dir <- '/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run'

setwd(code_dir)

# Load packages
library(pbapply)
library(ivreg)
library(lmtest)
library(sandwich)

# Load supporting functions
source(file.path(code_dir, 'B_gen_data.R'))
source(file.path(code_dir, 'B_AR_algo2.R'))

cat("\n", strrep("=", 80), "\n")
cat("ALGORITHM 2: COMPARING 1000 vs 10000 PERMUTATIONS\n")
cat(strrep("=", 80), "\n\n")

# ==============================================================================
# Load and prepare data (same for both)
# ==============================================================================

cat("Loading and preparing ALO dataset...\n")
data_file <- file.path(code_dir, "data", "ALO_data.csv")
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
cat("  Control (N0):", N0, "\n")
cat("  Compliance rate:", 
    sum(data_ssp_m1$D_observed * data_ssp_m1$assignment) / sum(data_ssp_m1$assignment), "\n\n")

# ==============================================================================
# RUN 1: Algorithm 2 with 1000 permutations
# ==============================================================================

cat(strrep("=", 80), "\n")
cat("RUN 1: ALGORITHM 2 WITH 1000 PERMUTATIONS\n")
cat(strrep("=", 80), "\n\n")

cat("Generating 1000 random permutations...\n")
zsim_1000 <- pbsapply(1:1000, gen_assignment_CR_index, N1=N1, N0=N0)
cat("  ✓ Permutations generated\n\n")

cat("Running Algorithm 2 with 1000 permutations...\n\n")
start_1000 <- proc.time()

CI_1000 <- AR_algo2_custom(
  data_table = data_ssp_m1,
  N1 = N1,
  N0 = N0,
  zsim = zsim_1000,
  tol = 1e-8,
  alpha = 0.95
)

time_1000 <- proc.time() - start_1000

cat("\nCompleted in", round(time_1000[3]/60, 2), "minutes\n\n")

# ==============================================================================
# RUN 2: Algorithm 2 with 10000 permutations
# ==============================================================================

cat(strrep("=", 80), "\n")
cat("RUN 2: ALGORITHM 2 WITH 10000 PERMUTATIONS\n")
cat(strrep("=", 80), "\n\n")

cat("Generating 10000 random permutations...\n")
zsim_10000 <- pbsapply(1:10000, gen_assignment_CR_index, N1=N1, N0=N0)
cat("  ✓ Permutations generated\n\n")

cat("Running Algorithm 2 with 10000 permutations...\n\n")
start_10000 <- proc.time()

CI_10000 <- AR_algo2_custom(
  data_table = data_ssp_m1,
  N1 = N1,
  N0 = N0,
  zsim = zsim_10000,
  tol = 1e-8,
  alpha = 0.95
)

time_10000 <- proc.time() - start_10000

cat("\nCompleted in", round(time_10000[3]/60, 2), "minutes\n\n")

# ==============================================================================
# COMPARISON
# ==============================================================================

cat("\n", strrep("=", 80), "\n")
cat("COMPARISON RESULTS\n")
cat(strrep("=", 80), "\n\n")

cat("Runtime Comparison:\n")
cat("  1000 permutations:  ", round(time_1000[3]/60, 2), " minutes\n", sep="")
cat("  10000 permutations: ", round(time_10000[3]/60, 2), " minutes\n", sep="")
cat("  Ratio: ", round(time_10000[3] / time_1000[3], 2), "x longer\n\n", sep="")

cat("Confidence Interval Comparison:\n\n")

cat("With 1000 permutations:\n")
if (length(CI_1000) == 0) {
  cat("  EMPTY SET\n")
} else {
  for (i in 1:length(CI_1000)) {
    cat("  [", round(CI_1000[[i]][1], 4), ", ", round(CI_1000[[i]][2], 4), "]\n", sep="")
  }
}

cat("\nWith 10000 permutations:\n")
if (length(CI_10000) == 0) {
  cat("  EMPTY SET\n")
} else {
  for (i in 1:length(CI_10000)) {
    cat("  [", round(CI_10000[[i]][1], 4), ", ", round(CI_10000[[i]][2], 4), "]\n", sep="")
  }
}

cat("\n")

# Compare interval widths
if (length(CI_1000) > 0 && length(CI_10000) > 0) {
  width_1000 <- CI_1000[[1]][2] - CI_1000[[1]][1]
  width_10000 <- CI_10000[[1]][2] - CI_10000[[1]][1]
  
  cat("Interval Width Comparison:\n")
  cat("  1000 permutations:  ", round(width_1000, 4), "\n", sep="")
  cat("  10000 permutations: ", round(width_10000, 4), "\n", sep="")
  cat("  Difference: ", round(abs(width_1000 - width_10000), 4), "\n\n", sep="")
  
  # Check if intervals are similar
  diff_left <- abs(CI_1000[[1]][1] - CI_10000[[1]][1])
  diff_right <- abs(CI_1000[[1]][2] - CI_10000[[1]][2])
  
  if (diff_left < 0.01 && diff_right < 0.01) {
    cat("✓ Both give very similar confidence intervals!\n")
  } else {
    cat("⚠ Intervals differ by more than 0.01\n")
  }
}

cat("\n", strrep("=", 80), "\n")
cat("SUMMARY TABLE\n")
cat(strrep("=", 80), "\n\n")

# Create comparison table
comparison_table <- data.frame(
  Metric = c(
    "Permutations",
    "Runtime",
    "CI Lower",
    "CI Upper",
    "CI Width"
  ),
  Run_1000 = c(
    "1000",
    paste(round(time_1000[3]/60, 2), "min"),
    ifelse(length(CI_1000) > 0, round(CI_1000[[1]][1], 4), "Empty"),
    ifelse(length(CI_1000) > 0, round(CI_1000[[1]][2], 4), "Empty"),
    ifelse(length(CI_1000) > 0, round(CI_1000[[1]][2] - CI_1000[[1]][1], 4), "N/A")
  ),
  Run_10000 = c(
    "10000",
    paste(round(time_10000[3]/60, 2), "min"),
    ifelse(length(CI_10000) > 0, round(CI_10000[[1]][1], 4), "Empty"),
    ifelse(length(CI_10000) > 0, round(CI_10000[[1]][2], 4), "Empty"),
    ifelse(length(CI_10000) > 0, round(CI_10000[[1]][2] - CI_10000[[1]][1], 4), "N/A")
  )
)

print(comparison_table, row.names = FALSE)

cat("\n")

# ==============================================================================
# Save results
# ==============================================================================

cat("Saving comparison results...\n")

comparison_results <- list(
  n_1000 = 1000,
  n_10000 = 10000,
  time_1000_seconds = time_1000[3],
  time_10000_seconds = time_10000[3],
  CI_1000 = CI_1000,
  CI_10000 = CI_10000,
  sample_size = nrow(data_ssp_m1),
  N1 = N1,
  N0 = N0
)

save(comparison_results, file=file.path(code_dir, 'comparison_1000_vs_10000.Rdata'))
cat("  Results saved to 'comparison_1000_vs_10000.Rdata'\n\n")

cat(strrep("=", 80), "\n")
cat("CONCLUSION\n")
cat(strrep("=", 80), "\n\n")


cat("DONE!\n")
cat(strrep("=", 80), "\n")