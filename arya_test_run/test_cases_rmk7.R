################################################################################
## Test Cases for AR_intersection
## Purpose: Verify AR_intersection correctly implements Remark 7
################################################################################

# Load the AR_intersection function
source("/Users/ag5276/Documents/Github/Rand_IV/arya test run/solve_coef_01.R")
source("/Users/ag5276/Documents/Github/Rand_IV/arya test run/calculate_intersections_02.R")
source("/Users/ag5276/Documents/Github/Rand_IV/arya test run/find_intervals_03.R")

################################################################################
## TEST 1: Linear functions with known intersection
################################################################################

test_1_linear <- function(){
  cat("\n===============================================\n")
  cat("TEST 1: Two linear functions\n")
  cat("===============================================\n")
  cat("AR1(β) = (2β + 1) / 1\n")
  cat("AR2(β) = (β + 3) / 1\n")
  cat("Expected: β = 2 (since 2β + 1 = β + 3 => β = 2)\n\n")
  
  coef1 <- c(0, 2, 1, 0, 0, 1)  # (2β + 1) / 1
  coef2 <- c(0, 1, 3, 0, 0, 1)  # (β + 3) / 1
  
  intersections <- AR_intersection(coef1, coef2)
  cat("Found", length(intersections), "intersection(s):", intersections, "\n\n")
  
  # Verify each intersection
  if(length(intersections) > 0){
    cat("Verification:\n")
    for(beta in intersections){
      result <- check_intersection(beta, coef1, coef2)
      cat(sprintf("  β = %.6f: AR1 = %.6f, AR2 = %.6f, diff = %.2e", 
                  result$beta, result$AR1, result$AR2, result$difference))
      if(result$is_valid){
        cat(" ✓\n")
      } else {
        cat(" ✗\n")
      }
    }
  }
  
  # Check if we found the expected intersection at β = 2
  if(length(intersections) == 1 && abs(intersections[1] - 2) < 1e-6){
    cat("\n✓ TEST 1 PASSED: Found expected intersection at β = 2\n")
    return(TRUE)
  } else {
    cat("\n✗ TEST 1 FAILED\n")
    return(FALSE)
  }
}

################################################################################
## TEST 2: Parallel curves (no intersection)
################################################################################

test_2_no_intersection <- function(){
  cat("\n===============================================\n")
  cat("TEST 2: Parallel quadratics (no intersection)\n")
  cat("===============================================\n")
  cat("AR1(β) = β² / 1\n")
  cat("AR2(β) = (β² + 1) / 1\n")
  cat("Expected: No intersections (curves never meet)\n\n")
  
  coef1 <- c(1, 0, 0, 0, 0, 1)  # β² / 1
  coef2 <- c(1, 0, 1, 0, 0, 1)  # (β² + 1) / 1
  
  intersections <- AR_intersection(coef1, coef2)
  cat("Found", length(intersections), "intersection(s)\n")
  
  if(length(intersections) == 0){
    cat("\n✓ TEST 2 PASSED: Correctly found no intersections\n")
    return(TRUE)
  } else {
    cat("\n✗ TEST 2 FAILED: Should have found no intersections\n")
    return(FALSE)
  }
}

################################################################################
## TEST 3: Identical functions
################################################################################

test_3_identical <- function(){
  cat("\n===============================================\n")
  cat("TEST 3: Identical functions\n")
  cat("===============================================\n")
  cat("AR1(β) = AR2(β) for all β\n")
  cat("Expected: Empty set (treated as no discrete intersections)\n\n")
  
  coef1 <- c(1, 2, 3, 4, 5, 6)
  coef2 <- c(1, 2, 3, 4, 5, 6)
  
  intersections <- AR_intersection(coef1, coef2)
  cat("Found", length(intersections), "intersection(s)\n")
  
  if(length(intersections) == 0){
    cat("\n✓ TEST 3 PASSED: Correctly returned empty set\n")
    return(TRUE)
  } else {
    cat("\n✗ TEST 3 FAILED: Should have returned empty set\n")
    return(FALSE)
  }
}

################################################################################
## TEST 4: Complex rational functions
################################################################################

test_4_complex <- function(){
  cat("\n===============================================\n")
  cat("TEST 4: Complex rational functions\n")
  cat("===============================================\n")
  cat("AR1(β) = (β² - 4) / (β² + 1)\n")
  cat("AR2(β) = (2β) / (β² + 1)\n")
  cat("Expected: Multiple intersections (up to 4)\n\n")
  
  coef1 <- c(1, 0, -4, 1, 0, 1)  # (β² - 4) / (β² + 1)
  coef2 <- c(0, 2, 0, 1, 0, 1)   # (2β) / (β² + 1)
  
  intersections <- AR_intersection(coef1, coef2)
  cat("Found", length(intersections), "intersection(s):", intersections, "\n\n")
  
  # Verify each intersection
  if(length(intersections) > 0){
    cat("Verification:\n")
    all_valid <- TRUE
    for(beta in intersections){
      result <- check_intersection(beta, coef1, coef2)
      cat(sprintf("  β = %.6f: AR1 = %.6f, AR2 = %.6f, diff = %.2e", 
                  result$beta, result$AR1, result$AR2, result$difference))
      if(result$is_valid){
        cat(" ✓\n")
      } else {
        cat(" ✗\n")
        all_valid <- FALSE
      }
    }
    
    if(all_valid){
      cat("\n✓ TEST 4 PASSED: All intersections are valid\n")
      return(TRUE)
    } else {
      cat("\n✗ TEST 4 FAILED: Some intersections are invalid\n")
      return(FALSE)
    }
  } else {
    cat("\n✗ TEST 4 FAILED: Should have found at least one intersection\n")
    return(FALSE)
  }
}

################################################################################
## TEST 5: Degree check (at most 4 intersections)
################################################################################

test_5_max_degree <- function(){
  cat("\n===============================================\n")
  cat("TEST 5: Verify at most 4 intersections (Remark 7)\n")
  cat("===============================================\n")
  cat("Remark 7 states: solution set has size at most 4\n")
  cat("Testing with random coefficients...\n\n")
  
  # Test with several random coefficient pairs
  max_found <- 0
  for(i in 1:10){
    coef1 <- runif(6, -10, 10)
    coef2 <- runif(6, -10, 10)
    
    intersections <- AR_intersection(coef1, coef2)
    n_intersections <- length(intersections)
    
    if(n_intersections > max_found){
      max_found <- n_intersections
    }
    
    if(n_intersections > 4){
      cat(sprintf("✗ Found %d intersections (> 4) in trial %d\n", n_intersections, i))
      cat("✗ TEST 5 FAILED: Violated degree constraint\n")
      return(FALSE)
    }
  }
  
  cat(sprintf("Maximum intersections found: %d (≤ 4)\n", max_found))
  cat("\n✓ TEST 5 PASSED: Never exceeded 4 intersections\n")
  return(TRUE)
}
################################################################################
## TEST 6: Check with real data 
################################################################################

################################################################################
## Test Case 6: Real Data from ALO Study
## Purpose: Verify AR_intersection works with actual coefficients from solve_coef_01
################################################################################

# Load required functions
source("calculate_intersections.R")
source("solve_coef_01.R")

################################################################################
## TEST 6: Real ALO Data
################################################################################

test_6_real_data <- function(){
  
  cat("\n===============================================\n")
  cat("TEST 6: Real ALO Data with solve_coef_01\n")
  cat("===============================================\n")
  
  # Load data
  data_dir <- "/Users/ag5276/Documents/ALO_data.csv"
  data <- read.csv(data_dir, header=TRUE, sep=",")

  
  # Filter: male students in SSP treatment or control
  data_ssp_m <- data[data$sex == "M" & (data$ssp == 1 | data$control == 1), ]
  
  # Remove missing GPA values
  data_ssp_m <- data_ssp_m[!is.na(data_ssp_m$GPA_year1), ]
  
  # Select relevant columns and rename
  data_subset <- data_ssp_m[, c('GPA_year1', 'ssp_p', 'ssp', 'gpa0')]
  colnames(data_subset) <- c('Y_observed', 'D_observed', 'assignment', 'x1')
  
  # Center covariate
  data_subset$x1 <- data_subset$x1 - mean(data_subset$x1)
  
  cat(sprintf("Sample size after filtering: %d\n", nrow(data_subset)))
  
  # Calculate N1 and N0
  N1 <- sum(data_subset$assignment)
  N0 <- nrow(data_subset) - N1
  
  cat(sprintf("Treated (N1): %d\n", N1))
  cat(sprintf("Control (N0): %d\n\n", N0))
  
  # Generate first coefficient vector using observed assignment
  cat("Generating coefficients from observed data...\n")
  coef1 <- solve_coef_01(data_subset, N1, N0)
  cat("coef1:", round(coef1, 6), "\n\n")
  
  # Generate second coefficient vector using permuted assignment
  cat("Generating coefficients from permuted assignment...\n")
  set.seed(123)  # For reproducibility
  data_subset2 <- data_subset
  data_subset2$assignment <- sample(data_subset$assignment)
  
  coef2 <- solve_coef_01(data_subset2, N1, N0)
  cat("coef2:", round(coef2, 6), "\n\n")
  
  # Find intersections
  cat("Finding intersections between AR functions...\n")
  intersections <- AR_intersection(coef1, coef2)
  
  cat(sprintf("Found %d intersection(s)\n", length(intersections)))
  
  if(length(intersections) > 0){
    cat("\nIntersection values:\n")
    print(intersections)
    
    # Verify each intersection (show first 3)
    cat("\nVerification by substitution:\n")
    n_to_check <- min(3, length(intersections))
    
    all_valid <- TRUE
    for(i in 1:n_to_check){
      beta <- intersections[i]
      result <- check_intersection(beta, coef1, coef2)
      
      cat(sprintf("  β = %.6f: AR1 = %.6f, AR2 = %.6f, diff = %.2e", 
                  result$beta, result$AR1, result$AR2, result$difference))
      
      if(result$is_valid){
        cat(" ✓\n")
      } else {
        cat(" ✗\n")
        all_valid <- FALSE
      }
    }
    
    if(length(intersections) > 3){
      cat(sprintf("  ... and %d more intersection(s)\n", length(intersections) - 3))
    }
    
    # Check degree constraint
    if(length(intersections) > 4){
      cat("\n✗ TEST 6 FAILED: More than 4 intersections found!\n")
      return(FALSE)
    }
    
    if(all_valid){
      cat("\n✓ TEST 6 PASSED: All checked intersections are valid\n")
      return(TRUE)
    } else {
      cat("\n✗ TEST 6 FAILED: Some intersections are invalid\n")
      return(FALSE)
    }
    
  } else {
    cat("\n✓ TEST 6 PASSED: No intersections found (valid result)\n")
    return(TRUE)
  }
}

################################################################################
## Run the test
################################################################################

if(interactive()){
  cat("\nReal data test loaded. Run: test_6_real_data()\n\n")
} else {
  # Run automatically if sourced
  test_6_real_data()
}

################################################################################
## TEST 7: find_intervals with simple coefficients
################################################################################

test_7_find_intervals_simple <- function(){
  cat("\n===============================================\n")
  cat("TEST 7: find_intervals with simple coefficients\n")
  cat("===============================================\n")
  cat("AR1(β) = (β² - 2β + 3) / (4β² - β + 2)\n")
  cat("AR2(β) = (2β² - β + 1) / (3β² - 2β + 1)\n")
  cat("Expected: Intervals where AR1(β) < AR2(β)\n\n")
  
  # Define coefficients
  coef1 <- c(1, -2, 3, 4, -1, 2)
  coef2 <- c(2, -1, 1, 3, -2, 1)
  
  # Find intervals
  result <- find_intervals(coef1, coef2)
  
  cat(sprintf("Found %d interval(s)\n", ncol(result)))
  
  if(ncol(result) > 0){
    cat("\nInterval boundaries:\n")
    for(i in 1:ncol(result)){
      cat(sprintf("  [%.6f, %.6f]\n", result[1, i], result[2, i]))
    }
    
    # Verify each interval
    cat("\nVerification (testing midpoint of each interval):\n")
    all_valid <- TRUE
    
    for(i in 1:ncol(result)){
      left <- result[1, i]
      right <- result[2, i]
      
      # Choose test point
      if(!is.infinite(left) && !is.infinite(right)){
        beta_test <- (left + right) / 2
      } else if(is.infinite(left)){
        beta_test <- right - 1
      } else {
        beta_test <- left + 1
      }
      
      # Check if AR1 < AR2 at test point
      is_below <- compare_AR(beta_test, coef1, coef2)
      
      cat(sprintf("  β = %.6f: AR1 < AR2? %s", beta_test, is_below))
      
      if(is_below){
        cat(" ✓\n")
      } else {
        cat(" ✗\n")
        all_valid <- FALSE
      }
    }
    
    if(all_valid){
      cat("\n✓ TEST 7 PASSED\n")
      return(TRUE)
    } else {
      cat("\n✗ TEST 7 FAILED\n")
      return(FALSE)
    }
    
  } else {
    cat("\nNo intervals found (AR1 never below AR2)\n")
    cat("\n✓ TEST 7 PASSED\n")
    return(TRUE)
  }
}

################################################################################
## TEST 8: find_intervals with no intersections (always below)
################################################################################

test_8_no_intersection_always_below <- function(){
  cat("\n===============================================\n")
  cat("TEST 8: No intersections (AR1 always below AR2)\n")
  cat("===============================================\n")
  cat("AR1(β) = (β² + 5) / (β² + 1)\n")
  cat("AR2(β) = (β² + 10) / (β² + 1)\n")
  cat("Expected: Single interval (-∞, ∞)\n\n")
  
  # Define coefficients
  coef1 <- c(1, 0, 5, 1, 0, 1)   # AR1 ≈ 5 for large |β|
  coef2 <- c(1, 0, 10, 1, 0, 1)  # AR2 ≈ 10 for large |β|
  
  # Check for intersections
  intersections <- AR_intersection(coef1, coef2)
  cat(sprintf("Intersections found: %d\n", length(intersections)))
  
  if(length(intersections) != 0){
    cat("\n✗ TEST 8 FAILED: Expected no intersections\n")
    return(FALSE)
  }
  
  # Find intervals
  result <- find_intervals(coef1, coef2)
  
  cat(sprintf("Found %d interval(s)\n", ncol(result)))
  
  # Should return single interval (-Inf, Inf)
  if(ncol(result) != 1){
    cat("\n✗ TEST 8 FAILED: Expected exactly 1 interval\n")
    return(FALSE)
  }
  
  if(!is.infinite(result[1, 1]) || result[1, 1] >= 0){
    cat("\n✗ TEST 8 FAILED: Left endpoint should be -Inf\n")
    return(FALSE)
  }
  
  if(!is.infinite(result[2, 1]) || result[2, 1] <= 0){
    cat("\n✗ TEST 8 FAILED: Right endpoint should be +Inf\n")
    return(FALSE)
  }
  
  # Verify at multiple test points
  cat("\nVerification at test points:\n")
  test_points <- c(-100, -10, 0, 10, 100)
  all_valid <- TRUE
  
  for(beta in test_points){
    is_below <- compare_AR(beta, coef1, coef2)
    cat(sprintf("  β = %.1f: AR1 < AR2? %s", beta, is_below))
    
    if(is_below){
      cat(" ✓\n")
    } else {
      cat(" ✗\n")
      all_valid <- FALSE
    }
  }
  
  if(all_valid){
    cat("\n✓ TEST 8 PASSED\n")
    return(TRUE)
  } else {
    cat("\n✗ TEST 8 FAILED\n")
    return(FALSE)
  }
}

