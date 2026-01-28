################################################################################
## AR_intersection Function
## Extracted from: 2_AR_inversion_algo1.R (lines 249-318)
################################################################################

AR_intersection <- function(coef1, coef2){
  
  # This function calculates the intersections of two AR functions
  # 
  # Arguments:
  #   coef1: vector of 6 coefficients [a, b, c, d, e, f] for first AR function
  #          AR1(β) = (a*β² + b*β + c) / (d*β² + e*β + f)
  #   coef2: vector of 6 coefficients for second AR function
  #
  # Returns:
  #   vector of real intersection points (β values where AR1(β) = AR2(β))
  
  # Extract coefficients
  a1 <- coef1[1]
  b1 <- coef1[2]
  c1 <- coef1[3]
  d1 <- coef1[4]
  e1 <- coef1[5]
  f1 <- coef1[6]
  
  a2 <- coef2[1]
  b2 <- coef2[2]
  c2 <- coef2[3]
  d2 <- coef2[4]
  e2 <- coef2[5]
  f2 <- coef2[6]
  
  # Coefficient of beta^4
  v4 <- a1*d2 - a2*d1
  
  # Coefficient of beta^3
  v3 <- a1*e2 - a2*e1 + b1*d2 - b2*d1
  
  # Coefficient of beta^2
  v2 <- a1*f2 - a2*f1 + b1*e2 - b2*e1 + c1*d2 - c2*d1
  
  # Coefficient of beta
  v1 <- b1*f2 - b2*f1 + c1*e2 - c2*e1
  
  # Constant term
  v0 <- c1*f2 - c2*f1
  
  tol <- 1e-8
  
  if (is.na(v0)){
    browser()
  }
  
  if (abs(v0)<tol & abs(v1)<tol & abs(v2)<tol & abs(v3)<tol & abs(v4)<tol){
    real_roots <- c()
  } else if(abs(v1)<tol & abs(v2)<tol & abs(v3)<tol &abs(v4)<tol){
    roots <- polyroot(c(v0))
    real_roots <- Re(roots[which(abs(Im(roots))<1e-10)])
  } else if(abs(v2)<tol & abs(v3)<tol &abs(v4)<tol){
    roots <- polyroot(c(v0,v1))
    real_roots <- Re(roots[which(abs(Im(roots))<1e-10)])
  } else if(abs(v3)<tol &abs(v4)<tol){
    roots <- polyroot(c(v0,v1,v2))
    real_roots <- Re(roots[which(abs(Im(roots))<1e-10)])
  } else if(abs(v4)<tol){
    roots <- polyroot(c(v0,v1,v2,v3))
    real_roots <- Re(roots[which(abs(Im(roots))<1e-10)])
  } else {
    roots <- polyroot(c(v0,v1,v2,v3,v4))
    real_roots <- Re(roots[which(abs(Im(roots))<1e-10)])
  }
  
  return(as.vector(real_roots))
}

################################################################################
## Helper Function: Verify an intersection point
################################################################################

check_intersection <- function(beta, coef1, coef2, tol = 1e-6){
  # Verify that beta is truly an intersection by evaluating both AR functions
  #
  # Arguments:
  #   beta: the beta value to check
  #   coef1: coefficients for first AR function
  #   coef2: coefficients for second AR function
  #   tol: numerical tolerance for considering values equal
  #
  # Returns:
  #   TRUE if beta is a valid intersection, FALSE otherwise
  
  # Evaluate AR1(beta)
  AR1 <- (coef1[1]*beta^2 + coef1[2]*beta + coef1[3]) / 
    (coef1[4]*beta^2 + coef1[5]*beta + coef1[6])
  
  # Evaluate AR2(beta)
  AR2 <- (coef2[1]*beta^2 + coef2[2]*beta + coef2[3]) / 
    (coef2[4]*beta^2 + coef2[5]*beta + coef2[6])
  
  # Calculate difference
  diff <- abs(AR1 - AR2)
  
  # Check if difference is within tolerance
  is_valid <- diff < tol
  
  return(list(
    beta = beta,
    AR1 = AR1,
    AR2 = AR2,
    difference = diff,
    is_valid = is_valid
  ))
}