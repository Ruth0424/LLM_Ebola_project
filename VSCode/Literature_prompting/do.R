## PROMPT: Generate code using odin to model the outbreak of Ebola, using this paper for reference 
# Link to "Modeling the Spread of Ebola" by Do and Lee attached
## INTERFACE: VS Code
## MODEL: MAI-Code-1-Flash
# Not same system as in paper (should be SLIRD)
## PROMPT: Write code for a compartmental model which follows the same structure as given in this paper 
## PROMPT: Modify the code to use the same compartments as in this paper (prompted with url rather than downloaded pdf) 
# Format still differs to paper

# Ebola outbreak model implemented with odin in R
# This version follows the same compartment structure as the Legrand-style
# Ebola paper: S, E, I, H, F, R, D.
#
# Compartments:
# S: susceptible
# E: exposed
# I: infectious in the community
# H: hospitalised
# F: funeral / infectious contacts
# R: recovered
# D: deceased

if (!requireNamespace("odin", quietly = TRUE)) {
  install.packages("odin", repos = "https://cloud.r-project.org")
}

library(odin)

ebola_paper_structure <- odin::odin({
  initial(S) <- N0 - E0 - I0 - H0 - F0 - R0 - D0
  initial(E) <- E0
  initial(I) <- I0
  initial(H) <- H0
  initial(F) <- F0
  initial(R) <- R0
  initial(D) <- D0

  # Transmission parameters
  beta_i <- user(0.50)   # transmission from infectious community cases
  beta_h <- user(0.25)   # transmission from hospitalised cases
  beta_f <- user(0.70)   # transmission from funeral contacts

  # Transition and removal rates
  sigma <- user(0.25)    # exposed -> infectious
  alpha_h <- user(0.15)  # hospitalisation of infectious cases
  alpha_d <- user(0.05)  # death from infectious community state
  alpha_f <- user(0.10)  # movement from hospital to funeral state
  gamma_i <- user(0.07)  # recovery from infectious state
  gamma_h <- user(0.10)  # recovery from hospitalised state
  gamma_f <- user(0.03)  # recovery from funeral state
  delta_h <- user(0.06)  # death during hospitalisation
  delta_f <- user(0.12)  # death during funeral state

  # Population and initial conditions
  N0 <- user(1000000)
  E0 <- user(0)
  I0 <- user(1)
  H0 <- user(0)
  F0 <- user(0)
  R0 <- user(0)
  D0 <- user(0)

  N <- S + E + I + H + F + R + D

  # Force of infection from the three infectious pathways
  lambda <- (beta_i * I + beta_h * H + beta_f * F) / N

  deriv(S) <- -lambda * S
  deriv(E) <- lambda * S - sigma * E
  deriv(I) <- sigma * E - (alpha_h + alpha_d + gamma_i) * I
  deriv(H) <- alpha_h * I - (alpha_f + gamma_h + delta_h) * H
  deriv(F) <- alpha_f * H - (gamma_f + delta_f) * F
  deriv(R) <- gamma_i * I + gamma_h * H + gamma_f * F
  deriv(D) <- alpha_d * I + delta_h * H + delta_f * F

  output(incidence) <- lambda * S
  output(total_infectious) <- I + H + F
})

model <- ebola_paper_structure$new(
  N0 = 1000000,
  E0 = 0,
  I0 = 1,
  H0 = 0,
  F0 = 0,
  R0 = 0,
  D0 = 0,
  beta_i = 0.50,
  beta_h = 0.25,
  beta_f = 0.70,
  sigma = 0.25,
  alpha_h = 0.15,
  alpha_d = 0.05,
  alpha_f = 0.10,
  gamma_i = 0.07,
  gamma_h = 0.10,
  gamma_f = 0.03,
  delta_h = 0.06,
  delta_f = 0.12
)

times <- seq(0, 180, by = 1)
out <- as.data.frame(model$run(times))

peak_idx <- which.max(out$total_infectious)
peak_day <- out$t[peak_idx]
peak_total_infectious <- out$total_infectious[peak_idx]
final_state <- tail(out, 1)

cat("Peak total infectious burden:", round(peak_total_infectious, 2), "on day", round(peak_day, 1), "\n")
cat("Final state after 180 days:\n")
print(final_state)

if (!dir.exists("output")) {
  dir.create("output")
}

png("output/ebola_paper_structure.png", width = 1400, height = 900, res = 150)
matplot(
  out$t,
  out[, c("S", "E", "I", "H", "F", "R", "D")],
  type = "l",
  lty = 1,
  lwd = 2,
  col = c("blue", "orange", "red", "purple", "brown", "green", "black"),
  xlab = "Time (days)",
  ylab = "Population",
  main = "Ebola outbreak model with paper-style compartments"
)
legend(
  "right",
  legend = c("Susceptible", "Exposed", "Community infectious", "Hospitalised", "Funeral", "Recovered", "Deceased"),
  col = c("blue", "orange", "red", "purple", "brown", "green", "black"),
  lty = 1,
  lwd = 2
)
dev.off()

cat("\nPlot saved to output/ebola_paper_structure.png\n")
