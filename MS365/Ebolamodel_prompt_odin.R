
## PROMPT: Code up a compartmental model for the disease spread of Ebola using Odin
## INTERFACE: Microsoft 365
## MODEL: GPT 5.5 Think Deeper

# install.packages("odin")
# install.packages("deSolve")

library(odin)
library(deSolve)

ebola_model <- odin::odin({
  ## Force of infection
  lambda <- beta_I * I / N + beta_D * D / N
  
  ## ODE system
  deriv(S) <- -lambda * S
  
  deriv(E) <-  lambda * S - sigma * E
  
  deriv(I) <-  sigma * E - gamma * I
  
  deriv(R) <- (1 - cfr) * gamma * I
  
  deriv(D) <-  cfr * gamma * I - delta * D
  
  deriv(B) <-  delta * D
  
  ## Total living population used in transmission denominator
  N <- S + E + I + R
  
  ## Initial conditions
  initial(S) <- S0
  initial(E) <- E0
  initial(I) <- I0
  initial(R) <- 0
  initial(D) <- 0
  initial(B) <- 0
  
  ## User parameters
  S0 <- user()
  E0 <- user()
  I0 <- user()
  
  beta_I <- user()   # transmission rate from infectious living cases
  beta_D <- user()   # transmission rate from dead but unburied cases
  
  sigma <- user()    # rate exposed -> infectious; approx 1 / incubation period
  gamma <- user()    # rate infectious -> recovered/dead; approx 1 / infectious period
  delta <- user()    # rate unsafe dead -> safely buried
  
  cfr <- user()      # case fatality ratio
})

## Parameter values below are illustrative only.
## You should calibrate them to outbreak-specific data.

pars <- list(
  S0 = 1e5,
  E0 = 5,
  I0 = 1,
  
  beta_I = 0.25,
  beta_D = 0.50,
  
  sigma = 1 / 9,     # mean latent period around 9 days
  gamma = 1 / 7,     # mean infectious period around 7 days
  delta = 1 / 2,     # mean time to safe burial around 2 days
  
  cfr = 0.50
)

mod <- ebola_model$new(user = pars)

times <- seq(0, 200, by = 0.2)

out <- as.data.frame(mod$run(times))

head(out)

matplot(
  out$t,
  out[, c("S", "E", "I", "R", "D", "B")],
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "Time / days",
  ylab = "Number of individuals"
)

legend(
  "right",
  legend = c("S", "E", "I", "R", "D", "B"),
  col = 1:6,
  lty = 1,
  lwd = 2,
  bty = "n"
)

library(ggplot2)
library(tidyr)

out_long <- pivot_longer(
  out,
  cols = c(S, E, I, R, D, B),
  names_to = "compartment",
  values_to = "number"
)

ggplot(out_long, aes(x = t, y = number, colour = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time / days",
    y = "Number of individuals",
    colour = "Compartment"
  ) +
  theme_minimal()

