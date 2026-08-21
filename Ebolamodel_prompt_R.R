
## PROMPT: Code up a compartmental model for the disease spread of Ebola using R 
## INTERFACE: Microsoft 365
## MODEL: GPT 5.5 Think Deeper


# Ebola compartmental model in R
# Compartments:
# S = susceptible
# E = exposed, infected but not yet infectious
# I = infectious in the community
# H = hospitalised / isolated infectious cases
# F = deceased but not yet safely buried
# R = recovered
# B = safely buried deceased individuals

# Install if needed:
# install.packages("deSolve")

library(deSolve)

# Define the Ebola model
ebola_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Total living population
    N <- S + E + I + H + R
    
    # Force of infection
    # Transmission can occur from:
    # - infectious individuals in the community
    # - hospitalised individuals
    # - deceased individuals before safe burial
    lambda <- (beta_I * I + beta_H * H + beta_F * F) / N
    
    # Differential equations
    dS <- -lambda * S
    
    dE <- lambda * S - sigma * E
    
    dI <- sigma * E -
      gamma_H * theta * I -
      gamma_R * (1 - theta) * (1 - delta_I) * I -
      gamma_D * (1 - theta) * delta_I * I
    
    dH <- gamma_H * theta * I -
      gamma_HR * (1 - delta_H) * H -
      gamma_HD * delta_H * H
    
    dF <- gamma_D * (1 - theta) * delta_I * I +
      gamma_HD * delta_H * H -
      gamma_F * F
    
    dR <- gamma_R * (1 - theta) * (1 - delta_I) * I +
      gamma_HR * (1 - delta_H) * H
    
    dB <- gamma_F * F
    
    list(c(dS, dE, dI, dH, dF, dR, dB))
  })
}

# -----------------------------
# Initial conditions
# -----------------------------

initial_state <- c(
  S = 9990,  # susceptible
  E = 5,     # exposed
  I = 5,     # infectious
  H = 0,     # hospitalised
  F = 0,     # deceased but unburied
  R = 0,     # recovered
  B = 0      # safely buried
)

# -----------------------------
# Parameter values
# -----------------------------
# These are illustrative values.
# For serious inference or forecasting, parameters should be estimated
# from outbreak-specific data.

parameters <- c(
  beta_I = 0.25,   # transmission rate from community infectious cases
  beta_H = 0.10,   # transmission rate from hospitalised cases
  beta_F = 0.50,   # transmission rate from deceased/unburied cases
  
  sigma = 1 / 9.4, # incubation rate, inverse of average incubation period
  
  theta = 0.6,     # proportion of infectious cases hospitalised
  
  delta_I = 0.7,   # case fatality ratio outside hospital
  delta_H = 0.5,   # case fatality ratio in hospital
  
  gamma_H = 1 / 5,    # rate of hospitalisation from infectious class
  gamma_R = 1 / 10,   # recovery rate outside hospital
  gamma_D = 1 / 8,    # death rate outside hospital
  
  gamma_HR = 1 / 12,  # recovery rate in hospital
  gamma_HD = 1 / 7,   # death rate in hospital
  
  gamma_F = 1 / 2     # burial rate
)

# -----------------------------
# Time period
# -----------------------------

times <- seq(0, 200, by = 1)

# -----------------------------
# Solve ODE system
# -----------------------------

output <- ode(
  y = initial_state,
  times = times,
  func = ebola_model,
  parms = parameters
)

output <- as.data.frame(output)

# -----------------------------
# Plot results
# -----------------------------

matplot(
  output$time,
  output[, c("S", "E", "I", "H", "F", "R", "B")],
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "Time, days",
  ylab = "Number of individuals",
  main = "Compartmental Model of Ebola Transmission"
)

legend(
  "right",
  legend = c("S", "E", "I", "H", "F", "R", "B"),
  col = 1:7,
  lty = 1,
  lwd = 2,
  bty = "n"
)

# -----------------------------
# Plot infectious-related classes only
# -----------------------------

matplot(
  output$time,
  output[, c("E", "I", "H", "F")],
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "Time, days",
  ylab = "Number of individuals",
  main = "Ebola Epidemic Dynamics"
)

legend(
  "topright",
  legend = c("Exposed", "Infectious", "Hospitalised", "Funeral transmission risk"),
  col = 1:4,
  lty = 1,
  lwd = 2,
  bty = "n"
)

# Approximate daily new infections
output$new_infections <- c(NA, -diff(output$S))

plot(
  output$time,
  output$new_infections,
  type = "l",
  lwd = 2,
  xlab = "Time, days",
  ylab = "New infections per day",
  main = "Daily Ebola Incidence"
)

ebola_model_intervention <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Intervention reduces transmission after day 60
    if (time >= 60) {
      beta_I_eff <- beta_I * 0.6
      beta_H_eff <- beta_H * 0.5
      beta_F_eff <- beta_F * 0.3
    } else {
      beta_I_eff <- beta_I
      beta_H_eff <- beta_H
      beta_F_eff <- beta_F
    }
    
    N <- S + E + I + H + R
    
    lambda <- (beta_I_eff * I + beta_H_eff * H + beta_F_eff * F) / N
    
    dS <- -lambda * S
    dE <- lambda * S - sigma * E
    
    dI <- sigma * E -
      gamma_H * theta * I -
      gamma_R * (1 - theta) * (1 - delta_I) * I -
      gamma_D * (1 - theta) * delta_I * I
    
    dH <- gamma_H * theta * I -
      gamma_HR * (1 - delta_H) * H -
      gamma_HD * delta_H * H
    
    dF <- gamma_D * (1 - theta) * delta_I * I +
      gamma_HD * delta_H * H -
      gamma_F * F
    
    dR <- gamma_R * (1 - theta) * (1 - delta_I) * I +
      gamma_HR * (1 - delta_H) * H
    
    dB <- gamma_F * F
    
    list(c(dS, dE, dI, dH, dF, dR, dB))
  })
}

output_intervention <- ode(
  y = initial_state,
  times = times,
  func = ebola_model_intervention,
  parms = parameters
)

output_intervention <- as.data.frame(output_intervention)

plot(
  output$time,
  output$I,
  type = "l",
  lwd = 2,
  col = "red",
  xlab = "Time, days",
  ylab = "Infectious individuals",
  main = "Effect of Intervention on Ebola Transmission"
)

lines(
  output_intervention$time,
  output_intervention$I,
  lwd = 2,
  col = "blue"
)

legend(
  "topright",
  legend = c("No intervention", "Intervention after day 60"),
  col = c("red", "blue"),
  lty = 1,
  lwd = 2,
  bty = "n"
)