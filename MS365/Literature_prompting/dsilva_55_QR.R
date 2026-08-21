## PROMPT: Code a compartmental model for the disease spread of Ebola using odin
# attached "Modeling Spatial Invasion of Ebola in West Africa" by D'Silva and Eisenberg
## INTERFACE: Microsoft 365
## MODEL: GPT 5.5 Quick Response

ebola <- odin::odin({
  
  ## Initial conditions
  initial(S) <- S0
  initial(E) <- E0
  initial(I1) <- I10
  initial(I2) <- I20
  initial(F) <- F0
  initial(R) <- R0
  
  ## Population size
  N <- S + E + I1 + I2 + F + R
  
  ## Force of infection
  lambda <- (beta1 * I1 +
               beta2 * I2 +
               betaF * F) / N
  
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
  
  ## Incidence tracking
  deriv(cases) <- alpha * E
  initial(cases) <- 0
  
  deriv(deaths) <- delta * I2
  initial(deaths) <- 0
  
  ## Parameters
  beta1 <- user()
  beta2 <- user()
  betaF <- user()
  
  alpha  <- user()
  gamma  <- user()
  delta  <- user()
  delta2 <- user()
  
  r1 <- user()
  r2 <- user()
  
  S0  <- user()
  E0  <- user()
  I10 <- user()
  I20 <- user()
  F0  <- user()
  R0  <- user()
})

library(odin)
library(deSolve)

mod <- ebola$new(
  beta1 = 0.15,
  beta2 = 0.35,
  betaF = 0.40,
  alpha = 1/9,
  gamma = 1/5,
  delta = 1/4,
  delta2 = 1/2,
  r1 = 0.05,
  r2 = 0.03,
  S0 = 1e6,
  E0 = 10,
  I10 = 5,
  I20 = 0,
  F0 = 0,
  R0 = 0
)

tt <- seq(0, 365, by = 1)
out <- mod$run(tt)

matplot(
  out[, "time"],
  out[, c("S","E","I1","I2","F","R")],
  type = "l", lty = 1,
  xlab = "Time (days)",
  ylab = "Population"
)
legend(
  "right",
  legend = c("S","E","I1","I2","F","R"),
  col = 1:6,
  lty = 1
)