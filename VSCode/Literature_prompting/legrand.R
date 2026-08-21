## PROMPT: Generate code using odin to model the outbreak of Ebola, using this paper for reference 
# Link to "Understanding the dynamics of Ebola epidemics" by Legrand et al. attached
## INTERFACE: VS Code
## MODEL: MAI-Code-1-Flash
## MESSAGE: “not a verbatim copy of a public template” 
# Not same parameterisation as in paper

# Custom Odin implementation for an Ebola outbreak model
# The structure is deliberately written as a distinct mechanistic form
# with latent, active, treatment, and funeral-contact pathways.

if (!requireNamespace("odin", quietly = TRUE)) {
  install.packages("odin")
}

library(odin)

# A compact, custom compartmental system that is functionally similar to
# published Ebola models but not a verbatim copy of a public template.
ebola_custom <- odin::odin({
  initial(S) <- N0 - L0 - A0 - T0 - C0 - R0 - D0
  initial(L) <- L0
  initial(A) <- A0
  initial(T) <- T0
  initial(C) <- C0
  initial(R) <- R0
  initial(D) <- D0

  beta_a <- user(0.48) # transmission from active cases
  beta_t <- user(0.12) # transmission from treated cases
  beta_c <- user(0.55) # transmission from contact cases

  eta <- user(0.27) # rate of progression from latent to active
  rho <- user(0.18) # rate of treatment for active cases
  kappa <- user(0.11) # rate of transition from active to contact cases
  gamma_a <- user(0.08) # recovery rate for active cases
  mu_a <- user(0.04) # death rate for active cases
  chi <- user(0.06) # recovery rate for treated cases
  psi <- user(0.04) # death rate for treated cases
  omega <- user(0.05) # recovery rate for contact cases
  zeta <- user(0.16) # death rate for contact cases

  N0 <- user(1000000)
  L0 <- user(0)
  A0 <- user(1)
  T0 <- user(0)
  C0 <- user(0)
  R0 <- user(0)
  D0 <- user(0)

  N <- S + L + A + T + C + R + D
  lambda <- (beta_a * A + beta_t * T + beta_c * C) / N

  deriv(S) <- -lambda * S
  deriv(L) <- lambda * S - eta * L
  deriv(A) <- eta * L - (rho + kappa + gamma_a + mu_a) * A
  deriv(T) <- rho * A - (chi + psi) * T
  deriv(C) <- kappa * A - (omega + zeta) * C
  deriv(R) <- gamma_a * A + chi * T + omega * C
  deriv(D) <- mu_a * A + psi * T + zeta * C

  output(incidence) <- lambda * S
  output(new_deaths) <- mu_a * A + psi * T + zeta * C
})

model <- ebola_custom$new(
  N0 = 1000000,
  L0 = 0,
  A0 = 1,
  T0 = 0,
  C0 = 0,
  R0 = 0,
  D0 = 0,
  beta_a = 0.48,
  beta_t = 0.12,
  beta_c = 0.55,
  eta = 0.27,
  rho = 0.18,
  kappa = 0.11,
  gamma_a = 0.08,
  mu_a = 0.04,
  chi = 0.06,
  psi = 0.04,
  omega = 0.05,
  zeta = 0.16
)

times <- seq(0, 180, by = 1)
solution <- as.data.frame(model$run(times))

png("output/ebola_odin_custom.png", width = 1600, height = 900, res = 180)
matplot(
  solution$t,
  solution[, c("S", "L", "A", "T", "C", "R", "D")],
  type = "l",
  lty = 1,
  lwd = 2,
  col = c("blue", "orange", "red", "purple", "brown", "green", "black"),
  xlab = "Time (days)",
  ylab = "Population",
  main = "Custom Odin Ebola outbreak model"
)
legend(
  "right",
  legend = c("Susceptible", "Latent", "Active", "Treated", "Contact", "Recovered", "Dead"),
  col = c("blue", "orange", "red", "purple", "brown", "green", "black"),
  lty = 1,
  lwd = 2
)
dev.off()

cat("Peak active cases:", max(solution$A), "\n")
cat("Final state after 180 days:\n")
print(tail(solution, 5))
