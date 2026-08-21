## INTERFACE: VS Code
## MODEL: MAI-Code-1-Flash

## PROMPT: Generate code using odin to model the outbreak of Ebola, using this paper for reference 
## Link to "A simple mathematical model for Ebola in Africa" by Berge et al. attached
# Generates code for SEIHRDB model
## PROMPT: Modify the code to use the same compartments as in this paper (prompted with url rather than downloaded pdf) 
# Generates code instead for SEIHRD model
## PROMPT: Match the exact notation used in the article
# Only now requested access to link for paper

# Ebola outbreak model implemented with odin in R
# This version follows the exact notation used in Berge et al. (2016):
# dS/dt = pi - (beta1 I + beta2 D + lambda P) S - mu S
# dI/dt = (beta1 I + beta2 D + lambda P) S - (mu + delta + gamma) I
# dR/dt = gamma I - mu R
# dD/dt = (mu + delta) I - b D
# dP/dt = sigma + xi I + alpha D - eta P

if (!requireNamespace("odin", quietly = TRUE)) {
  install.packages("odin")
}

library(odin)

# State variables from the paper:
# S - susceptible
# I - infectious
# R - recovered
# D - deceased / removed
# P - pathogen in the environment

ebola_odin <- odin::odin({
  initial(S) <- S0
  initial(I) <- I0
  initial(R) <- R0
  initial(D) <- D0
  initial(P) <- P0

  # Parameters from the paper
  pi <- user(1e-4)
  beta1 <- user(0.02)
  beta2 <- user(0.002)
  lambda <- user(0.02)
  mu <- user(0.002)
  delta <- user(0.05)
  gamma <- user(0.05)
  b <- user(0.02)
  sigma <- user(0.001)
  xi <- user(0.001)
  alpha <- user(0.001)
  eta <- user(0.01)

  # Initial conditions
  S0 <- user(1000)
  I0 <- user(1)
  R0 <- user(0)
  D0 <- user(0)
  P0 <- user(0)

  force_of_infection <- beta1 * I + beta2 * D + lambda * P

  deriv(S) <- pi - force_of_infection * S - mu * S
  deriv(I) <- force_of_infection * S - (mu + delta + gamma) * I
  deriv(R) <- gamma * I - mu * R
  deriv(D) <- (mu + delta) * I - b * D
  deriv(P) <- sigma + xi * I + alpha * D - eta * P

  output(incidence) <- force_of_infection * S
})

simulate_ebola <- function(params, times = seq(0, 365, by = 1)) {
  model <- ebola_odin$new(
    S0 = params[["S0"]],
    I0 = params[["I0"]],
    R0 = params[["R0"]],
    D0 = params[["D0"]],
    P0 = params[["P0"]],
    pi = params[["pi"]],
    beta1 = params[["beta1"]],
    beta2 = params[["beta2"]],
    lambda = params[["lambda"]],
    mu = params[["mu"]],
    delta = params[["delta"]],
    gamma = params[["gamma"]],
    b = params[["b"]],
    sigma = params[["sigma"]],
    xi = params[["xi"]],
    alpha = params[["alpha"]],
    eta = params[["eta"]]
  )

  as.data.frame(model$run(times))
}

params <- c(
  S0 = 1000,
  I0 = 1,
  R0 = 0,
  D0 = 0,
  P0 = 0,
  pi = 1e-4,
  beta1 = 0.02,
  beta2 = 0.002,
  lambda = 0.02,
  mu = 0.002,
  delta = 0.05,
  gamma = 0.05,
  b = 0.02,
  sigma = 0.001,
  xi = 0.001,
  alpha = 0.001,
  eta = 0.01
)

out <- simulate_ebola(params)

print(head(out))
print(tail(out))

if (interactive()) {
  matplot(
    out$t,
    out[, c("S", "I", "R", "D", "P")],
    type = "l",
    lty = 1,
    lwd = 2,
    col = c("blue", "red", "green", "black", "orange"),
    xlab = "Time (days)",
    ylab = "State value",
    main = "Ebola outbreak model following Berge et al. (2016)"
  )

  legend(
    "right",
    legend = c("Susceptible", "Infectious", "Recovered", "Deceased", "Pathogen in environment"),
    col = c("blue", "red", "green", "black", "orange"),
    lty = 1,
    lwd = 2
  )
}
