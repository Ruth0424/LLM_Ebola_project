## PROMPT: Code an SIR model using odin
## INTERFACE: Microsoft 365
## MODEL: auto function

## ISSUES:
# errors to do with dimensions
# labelling error - needed out$t but gave out$time

install.packages("odin")     # run once
install.packages("dde")      # solver backend

library(odin)
library(dde)

sir_model <- odin::odin({
  
  ## Differential equations
  deriv(S) <- -beta * S * I / N
  deriv(I) <- beta * S * I / N - gamma * I
  deriv(R) <- gamma * I
  
  ## Total population (constant)
  N <- S + I + R
  
  ## Parameters
  beta <- user()    # transmission rate
  gamma <- user()   # recovery rate
  
  ## Initial conditions
  initial(S) <- S0
  initial(I) <- I0
  initial(R) <- R0
  
  S0 <- user()
  I0 <- user()
  R0 <- user()
})

model <- sir_model$new(
  beta = 0.3,
  gamma = 0.1,
  S0 = 999,
  I0 = 1,
  R0 = 0
)

times <- seq(0, 160, by = 1)

out <- model$run(times)
out <- as.data.frame(out)

matplot(
  out$time,
  out[, c("S", "I", "R")],
  type = "l",
  lty = 1,
  col = c("blue", "red", "green"),
  xlab = "Time",
  ylab = "Population",
  main = "SIR Model (odin)"
)

legend(
  "right",
  legend = c("S", "I", "R"),
  col = c("blue", "red", "green"),
  lty = 1
)

# Error in dim(x) <- c(n, 1L) : attempt to set an attribute on NULL
## PROMPT: The code produced an error: Error in dim(x) <- c(n, 1L) : attempt to set an attribute on NULL

library(odin)
library(dde)

sir_model <- odin::odin({
  
  ## Time variable (important for odin)
  t <- time
  
  ## Total population
  N <- S + I + R
  
  ## Differential equations
  deriv(S) <- -beta * S * I / N
  deriv(I) <- beta * S * I / N - gamma * I
  deriv(R) <- gamma * I
  
  ## Parameters
  beta <- user()
  gamma <- user()
  
  ## Initial conditions
  initial(S) <- S0
  initial(I) <- I0
  initial(R) <- R0
  
  S0 <- user()
  I0 <- user()
  R0 <- user()
})

mod <- sir_model$new(
  beta = 0.3,
  gamma = 0.1,
  S0 = 999,
  I0 = 1,
  R0 = 0
)

times <- seq(0, 160, by = 1)

out <- mod$run(times)

# Ensure it's not NULL
stopifnot(!is.null(out))

out <- as.data.frame(out)

matplot(
  out$time,
  out[, c("S", "I", "R")],
  type = "l",
  lty = 1,
  col = c("blue", "red", "green"),
  xlab = "Time",
  ylab = "Population"
)

legend(
  "right",
  legend = c("S", "I", "R"),
  col = c("blue", "red", "green"),
  lty = 1
)

# Error: Reserved name 't' for lhs : t <- time # (line 1)
## PROMPT: The code produced an error: Error: Reserved name 't' for lhs : t <- time # (line 1)

library(odin)
library(dde)

sir_model <- odin::odin({
  
  ## Total population
  N <- S + I + R
  
  ## Differential equations
  deriv(S) <- -beta * S * I / N
  deriv(I) <- beta * S * I / N - gamma * I
  deriv(R) <- gamma * I
  
  ## Parameters
  beta <- user()
  gamma <- user()
  
  ## Initial conditions
  initial(S) <- S0
  initial(I) <- I0
  initial(R) <- R0
  
  S0 <- user()
  I0 <- user()
  R0 <- user()
})

mod <- sir_model$new(
  beta = 0.3,
  gamma = 0.1,
  S0 = 999,
  I0 = 1,
  R0 = 0
)

times <- seq(0, 160, by = 1)

out <- mod$run(times)
out <- as.data.frame(out)

matplot(
  out$time,
  out[, c("S", "I", "R")],
  type = "l",
  lty = 1,
  col = c("blue", "red", "green"),
  xlab = "Time",
  ylab = "Population"
)

legend(
  "right",
  legend = c("S", "I", "R"),
  col = c("blue", "red", "green"),
  lty = 1
)

# Error in dim(x) <- c(n, 1L) : attempt to set an attribute on NULL
## PROMPT: The code produced an error: Reserved name 't' for lhs : t <- time # (line 1)

install.packages("odin")     # run once
install.packages("dde")      # solver backend

library(odin)
library(dde)

sir_model <- odin::odin({
  
  ## Differential equations
  deriv(S) <- -beta * S * I / N
  deriv(I) <- beta * S * I / N - gamma * I
  deriv(R) <- gamma * I
  
  ## Total population (constant)
  N <- S + I + R
  
  ## Parameters
  beta <- user()    # transmission rate
  gamma <- user()   # recovery rate
  
  ## Initial conditions
  initial(S) <- S0
  initial(I) <- I0
  initial(R) <- R0
  
  S0 <- user()
  I0 <- user()
  R0 <- user()
})

model <- sir_model$new(
  beta = 0.3,
  gamma = 0.1,
  S0 = 999,
  I0 = 1,
  R0 = 0
)

times <- seq(0, 160, by = 1)

out <- model$run(times)
out <- as.data.frame(out)

matplot(
  out$t,
  out[, c("S", "I", "R")],
  type = "l",
  lty = 1,
  col = c("blue", "red", "green"),
  xlab = "Time",
  ylab = "Population",
  main = "SIR Model (odin)"
)

legend(
  "right",
  legend = c("S", "I", "R"),
  col = c("blue", "red", "green"),
  lty = 1
)

