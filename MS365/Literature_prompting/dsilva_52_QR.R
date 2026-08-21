## PROMPT: Code a compartmental model for the disease spread of Ebola using odin
# attached "Modeling Spatial Invasion of Ebola in West Africa" by D'Silva and Eisenberg
## INTERFACE: Microsoft 365
## MODEL: GPT 5.2 Quick Response

library(odin2)

ebola <- odin({
  
  ## Compartments
  initial(S) <- S0
  initial(E) <- E0
  initial(I1) <- I10
  initial(I2) <- I20
  initial(F) <- F0
  initial(R) <- R0
  
  ## Force of infection
  lambda <- beta1 * I1 +
    beta2 * I2 +
    betaF * F
  
  ## ODEs
  deriv(S) <- -lambda * S
  
  deriv(E) <- lambda * S -
    alpha * E
  
  deriv(I1) <- alpha * E -
    gamma * I1 -
    r1 * I1
  
  deriv(I2) <- gamma * I1 -
    delta * I2 -
    r2 * I2
  
  deriv(F) <- delta * I2 -
    delta2 * F
  
  deriv(R) <- r1 * I1 +
    r2 * I2
  
  ## Cumulative incidence and deaths
  initial(IC) <- 0
  initial(DC) <- 0
  
  deriv(IC) <- k_norm * alpha * E
  deriv(DC) <- k_norm * delta * I2
  
  ## Parameters
  beta1 <- user()
  beta2 <- user()
  betaF <- user()
  
  alpha <- user()
  gamma <- user()
  delta <- user()
  delta2 <- user()
  
  r1 <- user()
  r2 <- user()
  
  k_norm <- user()
  
  ## Initial conditions
  S0 <- user()
  E0 <- user()
  I10 <- user()
  I20 <- user()
  F0 <- user()
  R0 <- user()
})
