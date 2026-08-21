## PROMPT: Generate code using Odin from R in this file for a compartmental model to model an outbreak of Ebola 
## INTERFACE: VS Code

# Ebola outbreak model implemented with odin in R
# This script defines an SEIRD compartmental model for Ebola transmission.

if (!requireNamespace("odin", quietly = TRUE)) {
  install.packages("odin")
}

library(odin)

ebola_seird <- odin::odin({
  initial(S) <- N0 - E0 - I0 - R0 - D0
  initial(E) <- E0
  initial(I) <- I0
  initial(R) <- R0
  initial(D) <- D0

  # Transmission rate: probability of infection per infectious contact
  beta <- user(0.55)

  # Progression rate from exposed to infectious
  sigma <- user(0.35)

  # Recovery rate from infectious to recovered
  gamma <- user(0.15)

  # Death rate from infectious compartment
  mu <- user(0.08)

  # Initial population size
  N0 <- user(1000000)

  # Initial numbers in each compartment
  E0 <- user(0)
  I0 <- user(1)
  R0 <- user(0)
  D0 <- user(0)

  # Total population size
  N <- S + E + I + R + D

  # Force of infection, proportional to infectious individuals
  force_of_infection <- beta * I / N

  deriv(S) <- -force_of_infection * S
  deriv(E) <- force_of_infection * S - sigma * E
  deriv(I) <- sigma * E - (gamma + mu) * I
  deriv(R) <- gamma * I
  deriv(D) <- mu * I

  output(incidence) <- force_of_infection * S
})

model <- ebola_seird$new(
  # Total initial population
  N0 = 1000000,
  # Initial compartment values
  E0 = 0,
  I0 = 1,
  R0 = 0,
  D0 = 0,
  # Model parameters
  beta = 0.55,
  sigma = 0.35,
  gamma = 0.15,
  mu = 0.08
)


## PROMPT: Add code for a new model that also includes the transmission of infection by dead bodies and a hospitalised compartment 
# New model: SEIHRD with hospitalisation and transmission from dead bodies
ebola_seihrd <- odin::odin({
  initial(S) <- N0 - E0 - I0 - H0 - R0 - D0
  initial(E) <- E0
  initial(I) <- I0
  initial(H) <- H0
  initial(R) <- R0
  initial(D) <- D0

  # Transmission from infectious individuals
  beta_i <- user(0.55)

  # Transmission from infectious contact with dead bodies
  beta_d <- user(0.20)

  # Progression rate from exposed to infectious
  sigma <- user(0.35)

  # Hospitalisation rate of infectious individuals
  rho <- user(0.20)

  # Recovery rates
  gamma_i <- user(0.15)
  gamma_h <- user(0.10)

  # Mortality rates
  mu_i <- user(0.08)
  mu_h <- user(0.10)

  # Initial population size
  N0 <- user(1000000)

  # Initial compartment values
  E0 <- user(0)
  I0 <- user(1)
  H0 <- user(0)
  R0 <- user(0)
  D0 <- user(0)

  # Total population size
  N <- S + E + I + H + R + D

  # Force of infection from infectious individuals and dead bodies
  force_of_infection <- (beta_i * I + beta_d * D) / N

  deriv(S) <- -force_of_infection * S
  deriv(E) <- force_of_infection * S - sigma * E
  deriv(I) <- sigma * E - (rho + gamma_i + mu_i) * I
  deriv(H) <- rho * I - (gamma_h + mu_h) * H
  deriv(R) <- gamma_i * I + gamma_h * H
  deriv(D) <- mu_i * I + mu_h * H

  output(incidence) <- force_of_infection * S
})

model_hosp <- ebola_seihrd$new(
  # Total initial population
  N0 = 1000000,
  # Initial compartment values
  E0 = 0,
  I0 = 1,
  H0 = 0,
  R0 = 0,
  D0 = 0,
  # Model parameters
  beta_i = 0.55,
  beta_d = 0.20,
  sigma = 0.35,
  rho = 0.20,
  gamma_i = 0.15,
  gamma_h = 0.10,
  mu_i = 0.08,
  mu_h = 0.10
)

## PROMPT: Add code for a new model that includes infection from hospitalised infectious individuals and models the dead bodies no longer being infectious once safely buried 
## MESSAGE: “Response cleared due to possible match to public code, retrying with modified prompt.” 

# New model: hospital transmission with safe burial of dead bodies
# Deceased individuals are infectious only before burial, then move to a buried compartment
# where they no longer contribute to transmission.
ebola_hosp_burial <- odin::odin({
  initial(S) <- N0 - E0 - I0 - H0 - R0 - D0 - B0
  initial(E) <- E0
  initial(I) <- I0
  initial(H) <- H0
  initial(R) <- R0
  initial(D) <- D0
  initial(B) <- B0

  # Transmission from community infectious cases
  beta_i <- user(0.45)

  # Transmission from hospitalised infectious cases
  beta_h <- user(0.25)

  # Transmission from dead bodies before burial
  beta_d <- user(0.10)

  # Progression rate from exposed to infectious
  sigma <- user(0.35)

  # Rate of hospitalisation of infectious cases
  rho <- user(0.20)

  # Recovery rates
  gamma_i <- user(0.15)
  gamma_h <- user(0.10)

  # Death rates from community and hospitalised infectious states
  mu_i <- user(0.08)
  mu_h <- user(0.10)

  # Burial rate for deceased individuals
  psi <- user(0.05)

  # Initial population size
  N0 <- user(1000000)

  # Initial compartment values
  E0 <- user(0)
  I0 <- user(1)
  H0 <- user(0)
  R0 <- user(0)
  D0 <- user(0)
  B0 <- user(0)

  # Total population size
  N <- S + E + I + H + R + D + B

  # Force of infection from infectious, hospitalised and pre-burial dead individuals
  force_of_infection <- (beta_i * I + beta_h * H + beta_d * D) / N

  deriv(S) <- -force_of_infection * S
  deriv(E) <- force_of_infection * S - sigma * E
  deriv(I) <- sigma * E - (rho + gamma_i + mu_i) * I
  deriv(H) <- rho * I - (gamma_h + mu_h) * H
  deriv(R) <- gamma_i * I + gamma_h * H
  deriv(D) <- mu_i * I + mu_h * H - psi * D
  deriv(B) <- psi * D

  output(incidence) <- force_of_infection * S
})

model_hosp_burial <- ebola_hosp_burial$new(
  N0 = 1000000,
  E0 = 0,
  I0 = 1,
  H0 = 0,
  R0 = 0,
  D0 = 0,
  B0 = 0,
  beta_i = 0.45,
  beta_h = 0.25,
  beta_d = 0.10,
  sigma = 0.35,
  rho = 0.20,
  gamma_i = 0.15,
  gamma_h = 0.10,
  mu_i = 0.08,
  mu_h = 0.10,
  psi = 0.05
)

times <- seq(0, 180, by = 1)
solution <- as.data.frame(model$run(times))
solution_hosp <- as.data.frame(model_hosp$run(times))
solution_hosp_burial <- as.data.frame(model_hosp_burial$run(times))

matplot(
  solution[, "t"],
  solution[, c("S", "E", "I", "R", "D")],
  type = "l",
  lty = 1,
  lwd = 2,
  col = c("blue", "orange", "red", "green", "black"),
  xlab = "Time (days)",
  ylab = "Population",
  main = "Ebola SEIRD outbreak model (ODIN)"
)
legend(
  "right",
  legend = c("Susceptible", "Exposed", "Infectious", "Recovered", "Deceased"),
  col = c("blue", "orange", "red", "green", "black"),
  lty = 1,
  lwd = 2
)

matplot(
  solution_hosp[, "t"],
  solution_hosp[, c("S", "E", "I", "H", "R", "D")],
  type = "l",
  lty = 1,
  lwd = 2,
  col = c("blue", "orange", "red", "purple", "green", "black"),
  xlab = "Time (days)",
  ylab = "Population",
  main = "Ebola SEIHRD outbreak model with dead-body transmission (ODIN)"
)
legend(
  "right",
  legend = c("Susceptible", "Exposed", "Infectious", "Hospitalised", "Recovered", "Deceased"),
  col = c("blue", "orange", "red", "purple", "green", "black"),
  lty = 1,
  lwd = 2
)

matplot(
  solution_hosp_burial[, "t"],
  solution_hosp_burial[, c("S", "E", "I", "H", "R", "D", "B")],
  type = "l",
  lty = 1,
  lwd = 2,
  col = c("blue", "orange", "red", "purple", "green", "black", "gray"),
  xlab = "Time (days)",
  ylab = "Population",
  main = "Ebola hospital-burial model with non-infectious buried bodies (ODIN)"
)
legend(
  "right",
  legend = c("Susceptible", "Exposed", "Infectious", "Hospitalised", "Recovered", "Deceased", "Buried"),
  col = c("blue", "orange", "red", "purple", "green", "black", "gray"),
  lty = 1,
  lwd = 2
)

cat("\nFinal state after 180 days for the original SEIRD model:\n")
print(tail(solution, 1))

cat("\nFinal state after 180 days for the new SEIHRD model:\n")
print(tail(solution_hosp, 1))

cat("\nFinal state after 180 days for the hospital-burial model:\n")
print(tail(solution_hosp_burial, 1))
