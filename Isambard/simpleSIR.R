
## PROMPT: Provide R code for an SIR compartmental model 
## INTERFACE: Isambard vllm

#' SIR Compartmental Model
#'
#' Solves the classic SIR ODE system:
#'   dS/dt = -beta * S * I / N
#'   dI/dt = beta * S * I / N - gamma * I
#'   dR/dt = gamma * I
#'
#' @param times   Named vector of time points
#' @param S0      Initial susceptible count
#' @param I0      Initial infected count
#' @param R0      Initial recovered count
#' @param beta    Transmission rate
#' @param gamma   Recovery rate
#' @return        Data frame of S, I, R over time
#' @export
sir_model <- function(times = seq(0, 100, by = 1),
                      S0 = 999000,
                      I0 = 1000,
                      R0 = 0,
                      beta = 0.3,
                      gamma = 0.1) {
  
  library(deSolve)
  
  # Total population (conserved)
  N <- S0 + I0 + R0
  
  # Initial states
  y <- c(S = S0, I = I0, R = R0)
  
  # ODE system
  sir_ode <- function(t, state, parameters) {
    with(as.list(state), {
      dS <- -beta * S * I / N
      dI <-  beta * S * I / N - gamma * I
      dR <-  gamma * I
      list(c(dS, dI, dR))
    })
  }
  
  # Solve
  out <- ode(y = y, times = times, func = sir_ode, parms = list(beta = beta, gamma = gamma))
  
  as.data.frame(out)
}


#' Basic reproduction number R0
#'
#' @param beta Transmission rate
#' @param gamma Recovery rate
#' @return R0 value
reproduction_number <- function(beta, gamma) {
  beta / gamma
}


#' Quick demonstration with plotting
#'
#' @param S0,I0,R0,beta,gamma Passed to sir_model
#' @export
demo_sir <- function(S0 = 999000, I0 = 1000, R0 = 0,
                     beta = 0.3, gamma = 0.1) {
  
  library(ggplot2)
  library(tidyr)
  
  res <- sir_model(S0 = S0, I0 = I0, R0 = R0, beta = beta, gamma = gamma)
  
  res_long <- res %>%
    pivot_longer(-time, names_to = "compartment", values_to = "count")
  
  ggplot(res_long, aes(x = time, y = count, colour = compartment)) +
    geom_line(size = 1.2) +
    labs(title = "SIR Compartmental Model",
         subtitle = paste0("R0 = ", round(reproduction_number(beta, gamma), 2)),
         x = "Time",
         y = "Population",
         colour = "Compartment") +
    scale_colour_brewer(palette = "Set1") +
    theme_minimal()
}

# Example usage:
res  <- sir_model(beta = 0.3, gamma = 0.1)
demo_sir()
cat("R0 =", reproduction_number(0.3, 0.1), "\n")


## PROMPT: deSolve::ode() passes the argument named parms, not params as has been given in this model code 

#' SIR Compartmental Model
#'
#' Solves the classic SIR ODE system:
#'   dS/dt = -beta * S * I / N
#'   dI/dt = beta * S * I / N - gamma * I
#'   dR/dt = gamma * I
#'
#' @param times   Named vector of time points
#' @param S0      Initial susceptible count
#' @param I0      Initial infected count
#' @param R0      Initial recovered count
#' @param beta    Transmission rate
#' @param gamma   Recovery rate
#' @return        Data frame of S, I, R over time
sir_model <- function(times = seq(0, 100, by = 1),
                      S0 = 999000,
                      I0 = 1000,
                      R0 = 0,
                      beta = 0.3,
                      gamma = 0.1) {
  
  library(deSolve)
  
  N <- S0 + I0 + R0
  
  y <- c(S = S0, I = I0, R = R0)
  
  sir_ode <- function(t, state, parameters) {
    with(as.list(state), {
      dS <- -parameters$beta * S * I / N
      dI <-  parameters$beta * S * I / N - parameters$gamma * I
      dR <-  parameters$gamma * I
      list(c(dS, dI, dR))
    })
  }
  
  out <- ode(y = y, times = times, func = sir_ode,
             parms = list(beta = beta, gamma = gamma))   # <-- parms, not params
  
  as.data.frame(out)
}

reproduction_number <- function(beta, gamma) beta / gamma

demo_sir <- function(S0 = 999000, I0 = 1000, R0 = 0,
                     beta = 0.3, gamma = 0.1) {
  library(ggplot2)
  library(tidyr)
  
  res <- sir_model(S0 = S0, I0 = I0, R0 = R0, beta = beta, gamma = gamma)
  
  res %>%
    pivot_longer(-time, names_to = "compartment", values_to = "count") %>%
    ggplot(aes(x = time, y = count, colour = compartment)) +
    geom_line(size = 1.2) +
    labs(title = "SIR Compartmental Model",
         subtitle = paste0("R0 = ", round(reproduction_number(beta, gamma), 2)),
         x = "Time", y = "Population") +
    scale_colour_brewer(palette = "Set1") +
    theme_minimal()
}

res  <- sir_model(beta = 0.3, gamma = 0.1)
demo_sir()
cat("R0 =", reproduction_number(0.3, 0.1), "\n")
