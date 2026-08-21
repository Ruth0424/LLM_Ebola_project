## PROMPT: Code up a compartmental model for the disease spread of Ebola using R in this file 
## INTERFACE: VS Code

# Ebola compartmental model in R
# This script implements a simple SEIRD model:
# S: susceptible
# E: exposed
# I: infectious
# R: recovered
# D: deceased

if (!requireNamespace("deSolve", quietly = TRUE)) {
  install.packages("deSolve")
}

library(deSolve)

ebola_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + E + I + R + D
    beta <- parameters["beta"]
    sigma <- parameters["sigma"]
    gamma <- parameters["gamma"]
    mu <- parameters["mu"]

    force_of_infection <- beta * I / N

    dS <- -force_of_infection * S
    dE <- force_of_infection * S - sigma * E
    dI <- sigma * E - (gamma + mu) * I
    dR <- gamma * I
    dD <- mu * I

    list(c(dS, dE, dI, dR, dD))
  })
}

# Initial conditions
initial_state <- c(
  S = 999999,
  E = 0,
  I = 1,
  R = 0,
  D = 0
)

# Time grid
times <- seq(0, 180, by = 1)

# Parameters
params <- c(
  beta = 0.55,
  sigma = 0.35,
  gamma = 0.15,
  mu = 0.08
)

# Solve the model
solution <- ode(
  y = initial_state,
  times = times,
  func = ebola_model,
  parms = params
)

solution <- as.data.frame(solution)

# Plot the results
matplot(
  solution$time,
  solution[, c("S", "E", "I", "R", "D")],
  type = "l",
  lty = 1,
  lwd = 2,
  col = c("blue", "orange", "red", "green", "black"),
  xlab = "Time (days)",
  ylab = "Population",
  main = "Ebola SEIRD Compartmental Model"
)
legend(
  "right",
  legend = c("Susceptible", "Exposed", "Infectious", "Recovered", "Deceased"),
  col = c("blue", "orange", "red", "green", "black"),
  lty = 1,
  lwd = 2
)

# Print the final state
cat("\nFinal state after 180 days:\n")
print(tail(solution, 1))

## PROMPT: Add code for a model that incorporates very specific aspects of disease spread of Ebola 

# A more detailed Ebola-specific model with:
# - community transmission
# - hospital transmission
# - funeral/burial transmission
# - recovery and death outcomes

ebola_specific_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S + E + I + H + F + R + D

    lambda_c <- beta_c * I / N
    lambda_h <- beta_h * H / N
    lambda_f <- beta_f * F / N
    lambda <- lambda_c + lambda_h + lambda_f

    dS <- -lambda * S
    dE <- lambda * S - sigma * E
    dI <- sigma * E - (alpha_h + alpha_d + gamma_i) * I
    dH <- alpha_h * I - (alpha_f + gamma_h + delta_h) * H
    dF <- alpha_f * H - (gamma_f + delta_f) * F
    dR <- gamma_i * I + gamma_h * H + gamma_f * F
    dD <- alpha_d * I + delta_h * H + delta_f * F

    list(c(dS, dE, dI, dH, dF, dR, dD))
  })
}

initial_state_specific <- c(
  S = 999999,
  E = 0,
  I = 1,
  H = 0,
  F = 0,
  R = 0,
  D = 0
)

times_specific <- seq(0, 180, by = 1)

params_specific <- c(
  beta_c = 0.45,
  beta_h = 0.25,
  beta_f = 0.60,
  sigma = 0.25,
  alpha_h = 0.20,
  alpha_d = 0.05,
  alpha_f = 0.10,
  gamma_i = 0.08,
  gamma_h = 0.12,
  gamma_f = 0.03,
  delta_h = 0.06,
  delta_f = 0.15
)

solution_specific <- ode(
  y = initial_state_specific,
  times = times_specific,
  func = ebola_specific_model,
  parms = params_specific
)

solution_specific <- as.data.frame(solution_specific)

matplot(
  solution_specific$time,
  solution_specific[, c("S", "E", "I", "H", "F", "R", "D")],
  type = "l",
  lty = 1,
  lwd = 2,
  col = c("blue", "orange", "red", "purple", "brown", "green", "black"),
  xlab = "Time (days)",
  ylab = "Population",
  main = "Detailed Ebola Model"
)
legend(
  "right",
  legend = c("Susceptible", "Exposed", "Infectious", "Hospitalized", "Funeral", "Recovered", "Deceased"),
  col = c("blue", "orange", "red", "purple", "brown", "green", "black"),
  lty = 1,
  lwd = 2
)

cat("\nFinal state for detailed Ebola model after 180 days:\n")
print(tail(solution_specific, 1))
