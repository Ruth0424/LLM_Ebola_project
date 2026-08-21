## PROMPT: Generate code for an SIR compartmental model using odin from R 
## INTERFACE: VS Code

# (PROMPT: Adjust the plot code to accommodate S and R)
# SIR compartmental model implemented with odin in R
# This script simulates a simple susceptible-infectious-recovered epidemic.

if (!requireNamespace("odin", quietly = TRUE)) {
  install.packages("odin")
}

library(odin)

sir_model <- odin::odin({
  initial(S) <- N - I0 - R0
  initial(I) <- I0
  initial(R) <- R0

  beta <- user(0.3)
  gamma <- user(0.1)

  N <- user(1000)
  I0 <- user(1)
  R0 <- user(0)

  lambda <- beta * I / N

  deriv(S) <- -lambda * S
  deriv(I) <- lambda * S - gamma * I
  deriv(R) <- gamma * I

  output(incidence) <- lambda * S
})

model <- sir_model$new(
  N = 1000,
  I0 = 1,
  R0 = 0,
  beta = 0.3,
  gamma = 0.1
)

times <- seq(0, 200, by = 0.1)
solution <- as.data.frame(model$run(times))

matplot(
  solution$t,
  solution[, c("S", "I", "R")],
  type = "l",
  lty = 1,
  lwd = 2,
  col = c("blue", "red", "green"),
  xlab = "Time",
  ylab = "Population",
  main = "SIR model with odin"
)
legend(
  "right",
  legend = c("Susceptible", "Infectious", "Recovered"),
  col = c("blue", "red", "green"),
  lty = 1,
  lwd = 2
)

cat("\nFinal state after 200 time units:\n")
print(tail(solution, 1))

## PROMPT: Add code for a new model that extends the SIR model to account for the key aspects of an Ebola outbreak 

# Ebola-style extension of the SIR model
# Adds exposed, recovered, and deceased compartments to reflect
# latent infection, recovery, and mortality observed in Ebola outbreaks.
ebola_model <- odin::odin({
  initial(S) <- N - E0 - I0 - R0 - D0
  initial(E) <- E0
  initial(I) <- I0
  initial(R) <- R0
  initial(D) <- D0

  beta <- user(0.55)
  sigma <- user(0.35)
  gamma <- user(0.15)
  mu <- user(0.08)

  N <- user(1000000)
  E0 <- user(0)
  I0 <- user(1)
  R0 <- user(0)
  D0 <- user(0)

  lambda <- beta * I / N

  deriv(S) <- -lambda * S
  deriv(E) <- lambda * S - sigma * E
  deriv(I) <- sigma * E - (gamma + mu) * I
  deriv(R) <- gamma * I
  deriv(D) <- mu * I

  output(incidence) <- lambda * S
})

model_ebola <- ebola_model$new(
  N = 1000000,
  E0 = 0,
  I0 = 1,
  R0 = 0,
  D0 = 0,
  beta = 0.55,
  sigma = 0.35,
  gamma = 0.15,
  mu = 0.08
)

times_ebola <- seq(0, 180, by = 1)
solution_ebola <- as.data.frame(model_ebola$run(times_ebola))

matplot(
  solution_ebola$t,
  solution_ebola[, c("S", "E", "I", "R", "D")],
  type = "l",
  lty = 1,
  lwd = 2,
  col = c("blue", "orange", "red", "green", "black"),
  xlab = "Time (days)",
  ylab = "Population",
  main = "Ebola-style extension of the SIR model"
)
legend(
  "right",
  legend = c("Susceptible", "Exposed", "Infectious", "Recovered", "Deceased"),
  col = c("blue", "orange", "red", "green", "black"),
  lty = 1,
  lwd = 2
)

cat("\nFinal state after 180 days:\n")
print(tail(solution_ebola, 1))
