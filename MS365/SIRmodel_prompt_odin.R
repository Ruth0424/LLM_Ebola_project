
## PROMPT: Code an SIR model using odin
## INTERFACE: Microsoft 365
## MODEL: GPT 5.5 Think Deeper

# install.packages("odin")
# install.packages("ggplot2")

library(odin)
library(ggplot2)

sir_model <- odin({
  # Differential equations
  deriv(S) <- -beta * S * I / N
  deriv(I) <-  beta * S * I / N - gamma * I
  deriv(R) <-  gamma * I
  
  # Initial conditions
  initial(S) <- S_ini
  initial(I) <- I_ini
  initial(R) <- R_ini
  
  # Total population
  N <- S + I + R
  
  # User-defined parameters
  beta <- user()
  gamma <- user()
  
  S_ini <- user()
  I_ini <- user()
  R_ini <- user()
})

# Create model instance
sir <- sir_model$new(
  beta = 0.3,      # transmission rate
  gamma = 0.1,     # recovery rate
  S_ini = 999,
  I_ini = 1,
  R_ini = 0
)

# Run simulation
times <- seq(0, 160, by = 0.1)
out <- as.data.frame(sir$run(times))
out$time <- times

# Reshape for plotting
out_long <- reshape(
  out,
  varying = c("S", "I", "R"),
  v.names = "value",
  timevar = "compartment",
  times = c("S", "I", "R"),
  direction = "long"
)

# Plot
ggplot(out_long, aes(x = time, y = value, colour = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    title = "SIR model using odin",
    x = "Time",
    y = "Number of individuals",
    colour = "Compartment"
  ) +
  theme_minimal()

## PROMPT: Extend this SIR model to model key aspects of an Ebola outbreak 

library(odin)
library(ggplot2)

ebola_model <- odin({
  # ------------------------------------------------------------
  # Force of infection
  # ------------------------------------------------------------
  # Transmission can occur from:
  #   I = infectious community cases
  #   H = hospitalised / isolated cases
  #   F = infectious funerals or unsafe burials
  #
  # beta_I: community transmission rate
  # beta_H: hospital transmission rate
  # beta_F: funeral transmission rate
  
  lambda <- (beta_I * I + beta_H * H + beta_F * F) / N
  
  # ------------------------------------------------------------
  # Differential equations
  # ------------------------------------------------------------
  
  deriv(S) <- -lambda * S
  
  deriv(E) <- lambda * S - sigma * E
  
  deriv(I) <- sigma * E -
    gamma_I * I -
    hosp_rate * I -
    death_I * I
  
  deriv(H) <- hosp_rate * I -
    gamma_H * H -
    death_H * H
  
  deriv(F) <- death_I * I +
    death_H * H -
    burial_rate * F
  
  deriv(R) <- gamma_I * I +
    gamma_H * H
  
  deriv(B) <- burial_rate * F
  
  # ------------------------------------------------------------
  # Initial conditions
  # ------------------------------------------------------------
  
  initial(S) <- S_ini
  initial(E) <- E_ini
  initial(I) <- I_ini
  initial(H) <- H_ini
  initial(F) <- F_ini
  initial(R) <- R_ini
  initial(B) <- B_ini
  
  # ------------------------------------------------------------
  # Population size
  # ------------------------------------------------------------
  # B is excluded from the living population.
  # F is deceased but still potentially infectious, so exclude it
  # from N for density scaling.
  N <- S + E + I + H + R
  
  # ------------------------------------------------------------
  # User-defined parameters
  # ------------------------------------------------------------
  
  beta_I <- user()
  beta_H <- user()
  beta_F <- user()
  
  sigma <- user()       # rate from exposed to infectious
  gamma_I <- user()     # community recovery rate
  gamma_H <- user()     # hospital recovery rate
  
  hosp_rate <- user()   # rate of hospitalisation / isolation
  death_I <- user()     # death rate in community
  death_H <- user()     # death rate in hospital
  burial_rate <- user() # rate of safe burial
  
  S_ini <- user()
  E_ini <- user()
  I_ini <- user()
  H_ini <- user()
  F_ini <- user()
  R_ini <- user()
  B_ini <- user()
})

# example simulation

mod <- ebola_model$new(
  # Transmission rates
  beta_I = 0.5,   # community transmission
  beta_H = 0.2,   # reduced transmission in hospital/isolation
  beta_F = 0.2,   # funeral transmission can be high
  
  # Natural history
  sigma = 1 / 7.65,   # incubation period around 9 days
  
  # Recovery rates
  gamma_I = 163 / 2880,
  gamma_H = 1 / 13,
  
  # Hospitalisation / isolation
  hosp_rate = 1 / 6,
  
  # Mortality rates
  death_I = 25 / 576,
  death_H = 6 / 91,
  
  # Safe burial
  burial_rate = 1 / 2,
  
  # Initial conditions
  S_ini = 1000,
  E_ini = 0,
  I_ini = 1,
  H_ini = 0,
  F_ini = 0,
  R_ini = 0,
  B_ini = 0
)

times <- seq(0, 150, by = 0.25)

out <- as.data.frame(mod$run(times))
out$time <- times

out_long <- reshape(
  out,
  varying = c("S", "E", "I", "H", "F", "R", "B"),
  v.names = "value",
  timevar = "compartment",
  times = c("S", "E", "I", "H", "F", "R", "B"),
  direction = "long"
)

ggplot(out_long, aes(x = time, y = value, colour = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Extended Ebola outbreak model",
    x = "Time, days",
    y = "Number of individuals",
    colour = "Compartment"
  ) +
  theme_minimal()
tail(out)
