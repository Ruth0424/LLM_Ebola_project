## PROMPT: Code an SIR model using Odin
## INTERFACE: Microsoft 365
## MODEL: GPT 5.5 Think Deeper

# install.packages("odin")  # run once if needed
library(odin)

sir_generator <- odin::odin({
  # SIR equations, frequency-dependent transmission
  deriv(S) <- -beta * S * I / N
  deriv(I) <-  beta * S * I / N - gamma * I
  deriv(R) <-  gamma * I
  
  # Initial conditions
  initial(S) <- S0
  initial(I) <- I0
  initial(R) <- R0_init
  
  # Total population
  N <- S + I + R
  
  # Optional outputs
  output(incidence) <- beta * S * I / N
  output(R_eff) <- beta / gamma * S / N
  
  # User-supplied parameters, with defaults
  beta <- user(0.3)       # transmission rate
  gamma <- user(0.1)      # recovery rate
  S0 <- user(999)
  I0 <- user(1)
  R0_init <- user(0)
})

# Create model instance
mod <- sir_generator$new(
  beta = 0.3,
  gamma = 0.1,
  S0 = 999,
  I0 = 1,
  R0_init = 0
)

# Run simulation
times <- seq(0, 160, by = 1)
out <- mod$run(times)

# Inspect first rows
head(out)

# Plot S, I, R
matplot(
  out[, "t"],
  out[, c("S", "I", "R")],
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "Time",
  ylab = "Number of individuals",
  col = c("blue", "red", "darkgreen")
)

legend(
  "right",
  legend = c("Susceptible", "Infectious", "Recovered"),
  col = c("blue", "red", "darkgreen"),
  lty = 1,
  lwd = 2,
  bty = "n"
)

## PROMPT: Modify the code to fit a disease with flu-like symptoms which onset 
# 2-21 days after infection, is spread through body fluid contact with infected 
# bodies (alive and dead before buried), and treatment involves isolation in
# hospital

# install.packages("odin")  # run once if needed
library(odin)

ebola_like_generator <- odin::odin({
  
  # ------------------------------------------------------------
  # Force of infection
  # ------------------------------------------------------------
  # Transmission occurs from:
  #   I = symptomatic infectious people in the community
  #   H = hospitalised/isolated cases, with reduced infectiousness
  #   F = unburied dead bodies
  #
  # Exposed/incubating people are not infectious.
  lambda <- (beta_I * I + beta_H * H + beta_F * F) / N_alive
  
  # ------------------------------------------------------------
  # Susceptible population
  # ------------------------------------------------------------
  deriv(S) <- -lambda * S
  
  # ------------------------------------------------------------
  # Exposed/incubating classes
  # ------------------------------------------------------------
  # Three exposed stages give a more realistic incubation-period
  # distribution than one exponential E compartment.
  deriv(E1) <- lambda * S - sigma * E1
  deriv(E2) <- sigma * E1 - sigma * E2
  deriv(E3) <- sigma * E2 - sigma * E3
  
  # ------------------------------------------------------------
  # Symptomatic infectious cases in the community
  # ------------------------------------------------------------
  deriv(I) <- sigma * E3 -
    hosp_rate * I -
    gamma_I * I -
    mu_I * I
  
  # ------------------------------------------------------------
  # Hospitalised/isolated cases
  # ------------------------------------------------------------
  # Hospital isolation does not have to mean zero transmission:
  # beta_H can be set low to represent infection prevention/control.
  deriv(H) <- hosp_rate * I -
    gamma_H * H -
    mu_H * H
  
  # ------------------------------------------------------------
  # Funeral/unburied dead body compartment
  # ------------------------------------------------------------
  # People who die from I or H enter F until safe burial.
  deriv(F) <- mu_I * I + mu_H * H - burial_rate * F
  
  # ------------------------------------------------------------
  # Recovered and safely buried/dead compartments
  # ------------------------------------------------------------
  deriv(R) <- gamma_I * I + gamma_H * H
  deriv(B) <- burial_rate * F
  
  # ------------------------------------------------------------
  # Cumulative outputs
  # ------------------------------------------------------------
  deriv(cum_infections) <- lambda * S
  deriv(cum_hospitalised) <- hosp_rate * I
  deriv(cum_deaths) <- mu_I * I + mu_H * H
  
  # ------------------------------------------------------------
  # Initial conditions
  # ------------------------------------------------------------
  initial(S) <- S0
  initial(E1) <- E10
  initial(E2) <- E20
  initial(E3) <- E30
  initial(I) <- I0
  initial(H) <- H0
  initial(F) <- F0
  initial(R) <- R0_init
  initial(B) <- B0
  
  initial(cum_infections) <- E10 + E20 + E30 + I0 + H0 + R0_init
  initial(cum_hospitalised) <- H0
  initial(cum_deaths) <- F0 + B0
  
  # ------------------------------------------------------------
  # Population size
  # ------------------------------------------------------------
  N_alive <- S + E1 + E2 + E3 + I + H + R
  
  # ------------------------------------------------------------
  # Useful outputs
  # ------------------------------------------------------------
  output(E) <- E1 + E2 + E3
  output(infectious_total) <- I + H + F
  output(new_infections) <- lambda * S
  output(new_hospitalisations) <- hosp_rate * I
  output(new_deaths) <- mu_I * I + mu_H * H
  
  # ------------------------------------------------------------
  # User-supplied parameters
  # ------------------------------------------------------------
  
  # Transmission rates
  beta_I <- user(0.35)    # community symptomatic transmission
  beta_H <- user(0.05)    # hospital transmission, reduced by isolation/PPE
  beta_F <- user(0.70)    # transmission from unburied dead bodies
  
  # Incubation/onset period
  # Ebola-like onset is 2-21 days; mean often around 8-10 days.
  # With 3 exposed stages, use sigma = 3 / mean_incubation.
  mean_incubation <- user(9.0)
  sigma <- 3 / mean_incubation
  
  # Hospital isolation rate
  hosp_rate <- user(0.20) # rate symptomatic community cases enter hospital
  
  # Recovery rates
  gamma_I <- user(0.05)   # community recovery rate
  gamma_H <- user(0.08)   # hospital recovery rate
  
  # Mortality rates
  mu_I <- user(0.04)      # community death rate
  mu_H <- user(0.03)      # hospital death rate
  
  # Burial rate
  burial_rate <- user(0.50)
  
  # Initial population
  S0 <- user(999)
  E10 <- user(0)
  E20 <- user(0)
  E30 <- user(0)
  I0 <- user(1)
  H0 <- user(0)
  F0 <- user(0)
  R0_init <- user(0)
  B0 <- user(0)
})

mod <- ebola_like_generator$new(
  S0 = 999,
  I0 = 1,
  
  beta_I = 0.35,
  beta_H = 0.05,
  beta_F = 0.70,
  
  mean_incubation = 9.0,
  hosp_rate = 0.20,
  
  gamma_I = 0.05,
  gamma_H = 0.08,
  
  mu_I = 0.04,
  mu_H = 0.03,
  
  burial_rate = 0.50
)

times <- seq(0, 200, by = 0.25)
out <- mod$run(times)

head(out)

matplot(
  out[, "t"],
  out[, c("S", "E", "I", "H", "F", "R", "B")],
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "Time / days",
  ylab = "Number of individuals",
  col = c("blue", "orange", "red", "purple", "black", "darkgreen", "grey40")
)

legend(
  "right",
  legend = c(
    "Susceptible",
    "Exposed/incubating",
    "Infectious community",
    "Hospital isolated",
    "Dead unburied",
    "Recovered",
    "Buried/dead"
  ),
  col = c("blue", "orange", "red", "purple", "black", "darkgreen", "grey40"),
  lty = 1,
  lwd = 2,
  bty = "n"
)

