## PROMPT: Generate code using odin for a compartmental model to
# model an outbreak of Ebola. In the model account for: 
# an incubation period, infection in community and hospital, 
# transmission of infection by dead bodies, removal of dead 
# bodies by safe burial 
## INTERFACE: VS Code

# Ebola outbreak model implemented with odin in R
# Compartments:
# S  - susceptible
# E  - exposed (incubation period)
# Ic - infectious in the community
# Ih - infectious in hospital
# R  - recovered
# D  - dead bodies awaiting safe burial
# B  - buried bodies

if (!requireNamespace("odin", quietly = TRUE)) {
  install.packages("odin")
}

library(odin)

ebola_odin <- odin::odin({
  initial(S) <- N0 - E0 - Ic0 - Ih0 - R0 - D0 - B0
  initial(E) <- E0
  initial(Ic) <- Ic0
  initial(Ih) <- Ih0
  initial(R) <- R0
  initial(D) <- D0
  initial(B) <- B0

  # Transmission rates
  beta_c <- user(0.35)  # community transmission
  beta_h <- user(0.20)  # hospital transmission
  beta_d <- user(0.10)  # transmission from dead bodies

  # Progression and removal rates
  sigma <- user(0.20)   # incubation / exposed -> infectious
  rho <- user(0.15)     # hospitalisation of infectious individuals
  gamma_c <- user(0.10) # recovery in the community
  gamma_h <- user(0.08) # recovery in hospital
  mu_c <- user(0.05)    # death in the community
  mu_h <- user(0.08)    # death in hospital
  psi <- user(0.10)     # safe burial rate for dead bodies

  # Initial conditions and population size
  N0 <- user(10000)
  E0 <- user(0)
  Ic0 <- user(1)
  Ih0 <- user(0)
  R0 <- user(0)
  D0 <- user(0)
  B0 <- user(0)

  N <- S + E + Ic + Ih + R + D + B

  # Force of infection from community, hospital and dead bodies
  lambda <- (beta_c * Ic + beta_h * Ih + beta_d * D) / N

  deriv(S) <- -lambda * S
  deriv(E) <- lambda * S - sigma * E
  deriv(Ic) <- sigma * E - (rho + gamma_c + mu_c) * Ic
  deriv(Ih) <- rho * Ic - (gamma_h + mu_h) * Ih
  deriv(R) <- gamma_c * Ic + gamma_h * Ih
  deriv(D) <- mu_c * Ic + mu_h * Ih - psi * D
  deriv(B) <- psi * D

  output(incidence) <- lambda * S
})

model <- ebola_odin$new(
  N0 = 1001,
  E0 = 0,
  Ic0 = 1,
  Ih0 = 0,
  R0 = 0,
  D0 = 0,
  B0 = 0,
  beta_c = 0.5,
  beta_h = 0.2,
  beta_d = 0.2,
  sigma = 1/7.65,
  rho = 1/6,
  gamma_c = 163/2880,
  gamma_h = 1/13,
  mu_c = 25/576,
  mu_h = 6/91,
  psi = 0.5
)

times <- seq(0, 180, by = 1)
out <- as.data.frame(model$run(times))

print(head(out))
print(tail(out))

if (interactive()) {
  matplot(
    out$t,
    out[, c("S", "E", "Ic", "Ih", "R", "D", "B")],
    type = "l",
    lty = 1,
    lwd = 2,
    col = c("blue", "orange", "red", "purple", "green", "black", "gray"),
    xlab = "Time (days)",
    ylab = "Population",
    main = "Ebola outbreak model with odin"
  )
  legend(
    "right",
    legend = c("Susceptible", "Exposed", "Community infectious", "Hospital infectious", "Recovered", "Dead bodies", "Buried"),
    col = c("blue", "orange", "red", "purple", "green", "black", "gray"),
    lty = 1,
    lwd = 2
  )
}
