
## PROMPT: Provide R code using Odin for this SIR compartmental model 
## Isambard vllm

## odin_sir.odin (the Odin model file): ############################
# SIR compartmental model in Odin

parameters {
  beta  :  # transmission rate
    gamma :  # recovery rate
    N     :  # total population
}

state {
  S :  # susceptible
    I :  # infected
    R :  # recovered
}

init {
  S <- N * 0.999
  I <- N * 0.001
  R <- 0
}

derivatives {
  dS <- -beta * S * I / N
  dI <-  beta * S * I / N - gamma * I
  dR <-  gamma * I
}

include "odin_base.odin"
####################################################################


library(odin)

#' Compile and run the SIR Odin model
#'
#' @param beta  Transmission rate
#' @param gamma Recovery rate
#' @param N     Total population
#' @param t     Time span (e.g. 0:200)
#' @param init  Named list of initial conditions (overwrites defaults)
#' @return      Data frame with columns time, S, I, R
#' @export
sir_odin <- function(beta = 0.3,
                     gamma = 0.1,
                     N = 1e6,
                     t = 0:200,
                     init = NULL) {
  
  # Compile the Odin model from the .odin file
  model <- odin("(location of .odin file)")
  
  # Run the simulation
  out <- model$run(t = t, beta = beta, gamma = gamma, N = N)
  
  # Merge time back in
  result <- data.frame(time = t, out)
  
  if (!is.null(init)) {
    out <- model$run(t = t, beta = beta, gamma = gamma, N = N,
                     S0 = init$S0, I0 = init$I0, R0 = init$R0)
    result <- data.frame(time = t, out)
  }
  
  result
}


#' Basic reproduction number
reproduction_number <- function(beta, gamma) beta / gamma


#' Demo with plotting
#'
#' @param beta,gamma,N,t Passed to sir_odin
#' @export
demo_sir_odin <- function(beta = 0.3, gamma = 0.1, N = 1e6, t = 0:200) {
  
  library(ggplot2)
  library(tidyr)
  
  res <- sir_odin(beta = beta, gamma = gamma, N = N, t = t)
  
  res %>%
    pivot_longer(-time, names_to = "compartment", values_to = "count") %>%
    ggplot(aes(x = time, y = count, colour = compartment)) +
    geom_line(size = 1.2) +
    labs(title = "SIR Compartmental Model (Odin)",
         subtitle = paste0("R0 = ", round(reproduction_number(beta, gamma), 2)),
         x = "Time", y = "Population") +
    scale_colour_brewer(palette = "Set1") +
    theme_minimal()
}

# Example usage:
res  <- sir_odin(beta = 0.3, gamma = 0.1, N = 1e6)
demo_sir_odin()


## PROMPT: Odin syntax is more R-like than this block-based organisation. 
## Odin uses deriv() to define the system of ODEs, initial() for the initial 
## conditions, and parameter assignments like beta <- user(). 

library(odin)
library(ggplot2)
library(tidyr)

#' Compile and run the SIR Odin model
#'
#' @param beta  Transmission rate (default 0.3)
#' @param gamma Recovery rate (default 0.1)
#' @param S0    Initial susceptible count (default 999000)
#' @param I0    Initial infected count (default 1000)
#' @param t     Time span (e.g. 0:200)
#' @return      Data frame with columns t, S, I, R
#' @export
sir_odin <- function(beta = 0.3,
                     gamma = 0.1,
                     S0 = 999000,
                     I0 = 1000,
                     t = 0:200) {
  
  model <- odin(x = "(location of Odin file"))
  
  instance <- model$new(beta = beta, gamma = gamma, S0 = S0, I0 = I0)
  res <- instance$run(t = t)
  as.data.frame(res)
}


reproduction_number <- function(beta, gamma) beta / gamma


demo_sir_odin <- function(beta = 0.3, gamma = 0.1, S0 = 999000,
                          I0 = 1000, t = 0:200) {
  
  res <- sir_odin(beta = beta, gamma = gamma, S0 = S0, I0 = I0, t = t)
  
  res %>%
    pivot_longer(-t, names_to = "compartment", values_to = "count") %>%
    ggplot(aes(x = t, y = count, colour = compartment)) +
    geom_line(size = 1.2) +
    labs(title = "SIR Compartmental Model (Odin)",
         subtitle = paste0("R0 = ", round(reproduction_number(beta, gamma), 2)),
         x = "Time", y = "Population") +
    scale_colour_brewer(palette = "Set1") +
    theme_minimal()
}

demo_sir_odin()


## PROMPT: Embed the Odin model directly within the odin() call instead of 
## defining it in a .odin file

library(odin)
library(ggplot2)
library(tidyr)

sir_odin <- function(beta = 0.3,
                     gamma = 0.1,
                     S0 = 999000,
                     I0 = 1000,
                     t = 0:200) {
  
  model <- odin(x = c(
    "deriv(S) <- -beta * S * I / N",
    "deriv(I) <-  beta * S * I / N - gamma * I",
    "deriv(R) <-  gamma * I",
    "",
    "initial(S) <- S0",
    "initial(I) <- I0",
    "initial(R) <- 0",
    "",
    "N <- S0 + I0",
    "",
    "beta  <- user(0.3)",
    "gamma <- user(0.1)",
    "S0    <- user(999000)",
    "I0    <- user(1000)"
  ))
  
  instance <- model$new(beta = beta, gamma = gamma, S0 = S0, I0 = I0)
  as.data.frame(instance$run(t = t))
}


reproduction_number <- function(beta, gamma) beta / gamma


demo_sir_odin <- function(beta = 0.3, gamma = 0.1, S0 = 999000,
                          I0 = 1000, t = 0:200) {
  
  res <- sir_odin(beta = beta, gamma = gamma, S0 = S0, I0 = I0, t = t)
  
  res %>%
    pivot_longer(-t, names_to = "compartment", values_to = "count") %>%
    ggplot(aes(x = t, y = count, colour = compartment)) +
    geom_line(size = 1.2) +
    labs(title = "SIR Compartmental Model (Odin)",
         subtitle = paste0("R0 = ", round(reproduction_number(beta, gamma), 2)),
         x = "Time", y = "Population") +
    scale_colour_brewer(palette = "Set1") +
    theme_minimal()
}
demo_sir_odin()



## PROMPT: Instead, use this code for defining the model: model <- odin({ deriv(S)... I0 <- user()}) 

library(odin)
library(ggplot2)
library(tidyr)

sir_odin <- function(beta = 0.3,
                     gamma = 0.1,
                     S0 = 999000,
                     I0 = 1000,
                     t = 0:200) {
  
  model <- odin({
    deriv(S) <- -beta * S * I / N
    deriv(I) <-  beta * S * I / N - gamma * I
    deriv(R) <-  gamma * I
    
    initial(S) <- S0
    initial(I) <- I0
    initial(R) <- 0
    
    N <- S0 + I0
    
    beta  <- user(0.3)
    gamma <- user(0.1)
    S0    <- user(999000)
    I0    <- user(1000)
  })
  
  instance <- model$new(beta = beta, gamma = gamma, S0 = S0, I0 = I0)
  as.data.frame(instance$run(t = t))
}


reproduction_number <- function(beta, gamma) beta / gamma


demo_sir_odin <- function(beta = 0.3, gamma = 0.1, S0 = 999000,
                          I0 = 1000, t = 0:200) {
  
  res <- sir_odin(beta = beta, gamma = gamma, S0 = S0, I0 = I0, t = t)
  
  res %>%
    pivot_longer(-t, names_to = "compartment", values_to = "count") %>%
    ggplot(aes(x = t, y = count, colour = compartment)) +
    geom_line(size = 1.2) +
    labs(title = "SIR Compartmental Model (Odin)",
         subtitle = paste0("R0 = ", round(reproduction_number(beta, gamma), 2)),
         x = "Time", y = "Population") +
    scale_colour_brewer(palette = "Set1") +
    theme_minimal()
}
demo_sir_odin()
