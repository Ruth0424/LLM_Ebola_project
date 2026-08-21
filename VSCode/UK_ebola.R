## PROMPT: Generate code using odin for a compartmental model to model an outbreak of Ebola specifically within the UK 
## INTERFACE: VS Code
## MESSAGE: “I will inspect the existing R scripts in the workspace to mirror the project’s style and then add a UK-specific Ebola compartmental model using odin.” 

# UK Ebola outbreak model with odin
# This script implements a compartmental model for an Ebola outbreak
# in the United Kingdom using an odin-based ODE system.

if (!requireNamespace("odin", quietly = TRUE)) {
  install.packages("odin", repos = "https://cloud.r-project.org")
}

library(odin)

uk_ebola <- odin::odin({
  # Initial conditions
  initial(S) <- N0 - E0 - Ic0 - Ih0 - R0 - D0 - B0
  initial(E) <- E0
  initial(Ic) <- Ic0
  initial(Ih) <- Ih0
  initial(R) <- R0
  initial(D) <- D0
  initial(B) <- B0

  # Population and initial conditions for the UK
  N0 <- user(67000000)  # approximate UK population
  E0 <- user(0)
  Ic0 <- user(1)        # one imported infectious case
  Ih0 <- user(0)
  R0 <- user(0)
  D0 <- user(0)
  B0 <- user(0)

  # Transmission parameters
  beta_c <- user(0.22)  # community transmission
  beta_h <- user(0.09)  # hospital transmission
  beta_d <- user(0.04)  # transmission from dead bodies

  # Progression and removal rates
  sigma <- user(0.20)   # exposed -> infectious
  rho <- user(0.12)     # hospitalisation
  gamma_c <- user(0.08) # recovery in community
  gamma_h <- user(0.06) # recovery in hospital
  mu_c <- user(0.05)    # death in community
  mu_h <- user(0.08)    # death in hospital
  psi <- user(0.05)     # safe burial rate

  # Total population
  N <- S + E + Ic + Ih + R + D + B

  # Force of infection
  lambda <- (beta_c * Ic + beta_h * Ih + beta_d * D) / N

  # Differential equations
  deriv(S) <- -lambda * S
  deriv(E) <- lambda * S - sigma * E
  deriv(Ic) <- sigma * E - (rho + gamma_c + mu_c) * Ic
  deriv(Ih) <- rho * Ic - (gamma_h + mu_h) * Ih
  deriv(R) <- gamma_c * Ic + gamma_h * Ih
  deriv(D) <- mu_c * Ic + mu_h * Ih - psi * D
  deriv(B) <- psi * D

  # Outputs for monitoring
  output(incidence) <- lambda * S
  output(new_cases) <- sigma * E
})

# Create the model object
model <- uk_ebola$new(
  N0 = 67000000,
  E0 = 0,
  Ic0 = 1,
  Ih0 = 0,
  R0 = 0,
  D0 = 0,
  B0 = 0,
  beta_c = 0.25,
  beta_h = 0.01,
  beta_d = 0.05,
  sigma = 1/9,
  rho = 0.3,
  gamma_c = 0.08,
  gamma_h = 0.06,
  mu_c = 0.4,
  mu_h = 0.2,
  psi = 1
)

# Simulation time grid
times <- seq(0, 180, by = 1)

# Run the model
out <- as.data.frame(model$run(times))

# Print a summary
print(head(out, 10))
print(tail(out, 5))

# Plot the results when run interactively
if (interactive()) {
  matplot(
    out$t,
    out[, c("E", "Ic", "Ih", "R", "D", "B")],
    type = "l",
    lty = 1,
    lwd = 2,
    col = c("orange", "red", "purple", "green", "black", "gray"),
    xlab = "Time (days)",
    ylab = "Population",
    main = "Ebola outbreak model for the UK"
  )
  legend(
    "right",
    legend = c("Exposed", "Community infectious", "Hospital infectious", "Recovered", "Dead bodies", "Buried"),
    col = c("orange", "red", "purple", "green", "black", "gray"),
    lty = 1,
    lwd = 2
  )
}
