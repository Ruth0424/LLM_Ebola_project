
## PROMPT: Describe a compartmental model that would be suitable for modelling an outbreak of Ebola 
## INTERFACE: Isambard vllm

## Include a rate of burial to remove individuals from the D compartment, 
## and use fatality probabilities and overall removal rates for the 
## parameters leaving I and H, so that units are consistent (day-1) 

library(odin)

#' Ebola SEIR-HPD compartmental model
#'
#' Compartments:
#'   S - Susceptible
#'   E - Exposed (incubating, not infectious)
#'   I - Infectious (symptomatic, community-transmitting)
#'   H - Hospitalised (isolated, reduced transmission)
#'   R - Recovered (immune)
#'   D - Deceased (infectious through unsafe burial)
#'
#' The ODE system:
#'   dS/dt = -beta * S * (I + theta_H * H + theta_D * D) / N + omega * R
#'   dE/dt =  beta * S * (I + theta_H * H + theta_D * D) / N - sigma * E
#'   dI/dt =  sigma * E - lambda_I * I
#'   dH/dt =  p_H * lambda_I * I - lambda_H * H
#'   dR/dt =  p_R * lambda_I * I + p_RH * lambda_H * H - omega * R
#'   dD/dt =  p_D * lambda_I * I + p_DH * lambda_H * H - lambda_B * D
#'
#' where p_R = 1 - p_H - p_D and p_RH = 1 - p_DH are derived internally.
#'
#' @param beta   Transmission rate (day^-1)
#' @param sigma  Incubation rate (day^-1), mean incubation = 1/sigma
#' @param lambda_I  Overall removal rate from I (day^-1)
#' @param p_H       Probability of hospitalisation given removal from I
#' @param p_D       Probability of death given removal from I
#' @param lambda_H  Overall removal rate from H (day^-1)
#' @param p_DH      Probability of death given removal from H
#' @param theta_H   Relative transmission while hospitalised
#' @param theta_D   Relative transmission from infectious corpses
#' @param lambda_B  Burial rate (day^-1)
#' @param omega     Immunity waning rate (day^-1)
#' @param N0        Total population
#' @param S0        Initial susceptible
#' @param E0        Initial exposed
#' @param I0        Initial infectious
#' @param H0        Initial hospitalised
#' @param R0        Initial recovered
#' @param D0        Initial deceased
#' @param t         Time span
#' @return          Data frame with columns t, S, E, I, H, R, D
#' @export
ebola_odin <- function(beta = 0.3,
                       sigma = 1/8,
                       lambda_I = 1/10,
                       p_H = 0.3,
                       p_D = 0.4,
                       lambda_H = 1/7,
                       p_DH = 0.15,
                       theta_H = 0.05,
                       theta_D = 0.02,
                       lambda_B = 1/5,
                       omega = 0,
                       N0 = 1e6,
                       S0 = NA,
                       E0 = 0,
                       I0 = 100,
                       H0 = 0,
                       R0 = 0,
                       D0 = 0,
                       t = 0:200) {
  
  if (is.na(S0)) {
    S0 <- N0 - I0 - E0 - H0 - R0 - D0
  }
  
  # Derived probabilities
  p_R  <- 1 - p_H - p_D       # recovery probability from I
  p_RH <- 1 - p_DH             # recovery probability from H
  
  # Sanity check: probabilities must sum to <= 1
  if (p_H + p_D > 1) stop("p_H + p_D must be <= 1")
  if (p_DH > 1) stop("p_DH must be <= 1")
  
  model <- odin({
    deriv(S) <- -beta * S * (I + theta_H * H + theta_D * D) / N + omega * R
    
    deriv(E) <-  beta * S * (I + theta_H * H + theta_D * D) / N - sigma * E
    
    deriv(I) <-  sigma * E - lambda_I * I
    
    deriv(H) <-  p_H * lambda_I * I - lambda_H * H
    
    deriv(R) <-  p_R * lambda_I * I + p_RH * lambda_H * H - omega * R
    
    deriv(D) <-  p_D * lambda_I * I + p_DH * lambda_H * H - lambda_B * D
    
    initial(S) <- S0
    initial(E) <- E0
    initial(I) <- I0
    initial(H) <- H0
    initial(R) <- R0
    initial(D) <- D0
    
    N <- S0 + I0 + E0 + H0 + R0 + D0
    
    beta       <- user(0.3)
    sigma      <- user(1/8)
    lambda_I   <- user(1/10)
    p_H        <- user(0.3)
    p_D        <- user(0.4)
    lambda_H   <- user(1/7)
    p_DH       <- user(0.15)
    theta_H    <- user(0.05)
    theta_D    <- user(0.02)
    lambda_B   <- user(1/5)
    omega      <- user(0)
    S0         <- user(999900)
    E0         <- user(0)
    I0         <- user(100)
    H0         <- user(0)
    R0         <- user(0)
    D0         <- user(0)
  })
  
  instance <- model$new(beta = beta, sigma = sigma,
                        lambda_I = lambda_I, p_H = p_H, p_D = p_D,
                        lambda_H = lambda_H, p_DH = p_DH,
                        theta_H = theta_H, theta_D = theta_D,
                        lambda_B = lambda_B, omega = omega,
                        S0 = S0, E0 = E0, I0 = I0, H0 = H0,
                        R0 = R0, D0 = D0)
  as.data.frame(instance$run(t = t))
}


#' Basic reproduction number
#'
#' Derivation: a single infectious individual in a fully susceptible
#' population generates:
#'   - beta / sigma * 1 direct infections during incubation delay
#'   - Of those, (1 - p_H - p_D) recover, (p_H) get hospitalised,
#'     (p_D) die.
#'   - Each hospitalised individual transmits further: p_H * beta / lambda_H
#'   - Each deceased transmits further: p_D * theta_D * lambda_D / lambda_B
#'
#' R0 = beta / sigma * (1 + p_H * theta_H / lambda_H + p_D * theta_D / lambda_B)
#'
#' @param beta  Transmission rate
#' @param sigma Incubation rate
#' @param p_H   Probability of hospitalisation
#' @param p_D   Probability of death
#' @param lambda_H Overall removal rate from hospitalised
#' @param lambda_B Burial rate
#' @param theta_H Relative transmission while hospitalised
#' @param theta_D Relative transmission from corpses
#' @return R0
reproduction_number <- function(beta, sigma,
                                p_H, p_D,
                                lambda_H, lambda_B,
                                theta_H = 0.05, theta_D = 0.02) {
  
  exposure_duration <- 1 / sigma
  community_transmission <- 1 + p_H * theta_H / lambda_H + p_D * theta_D / lambda_B
  beta / sigma * community_transmission
}


#' Ebola demonstration with plotting
#'
#' @param ... Passed to [ebola_odin()]
#' @export
demo_ebola <- function(beta = 0.3, sigma = 1/8, lambda_I = 1/10,
                       p_H = 0.3, p_D = 0.4, lambda_H = 1/7, p_DH = 0.15,
                       theta_H = 0.05, theta_D = 0.02, lambda_B = 1/5,
                       omega = 0, S0 = 999900, E0 = 0, I0 = 100,
                       H0 = 0, R0 = 0, D0 = 0, t = 0:200) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    stop("Requires ggplot2 and tidyr packages")
  }
  
  library(ggplot2)
  library(tidyr)
  
  res <- ebola_odin(beta = beta, sigma = sigma, lambda_I = lambda_I,
                    p_H = p_H, p_D = p_D, lambda_H = lambda_H, p_DH = p_DH,
                    theta_H = theta_H, theta_D = theta_D, lambda_B = lambda_B,
                    omega = omega, S0 = S0, E0 = E0, I0 = I0, H0 = H0,
                    R0 = R0, D0 = D0, t = t)
  
  res <- res %>%
    tidyr::pivot_longer(-t, names_to = "compartment", values_to = "count")
  
  ggplot(res, aes(x = t, y = count, colour = compartment)) +
    geom_line(size = 1.2) +
    labs(title = "Ebola SEIR-HPD Model",
         subtitle = paste0("R0 = ", round(reproduction_number(beta, sigma,
                                                              p_H, p_D,
                                                              lambda_H,
                                                              lambda_B,
                                                              theta_H,
                                                              theta_D), 2),
                           "\nMean incubation = ", round(1/sigma, 1),
                           " days, Mean time-to-removal (I) = ", round(1/lambda_I, 1),
                           " days"),
         x = "Time (days)", y = "Population") +
    scale_colour_brewer(palette = "Set2") +
    theme_minimal()
}
demo_ebola()


## PROMPT: Define p_R and p_RH differently, so that they are recognised as known variables within odin({})

library(odin)

#' Ebola SEIR-HPD compartmental model
#'
#' Compartments:
#'   S - Susceptible
#'   E - Exposed (incubating, not infectious)
#'   I - Infectious (symptomatic, community-transmitting)
#'   H - Hospitalised (isolated, reduced transmission)
#'   R - Recovered (immune)
#'   D - Deceased (infectious through unsafe burial)
#'
#' The ODE system:
#'   dS/dt = -beta * S * (I + theta_H * H + theta_D * D) / N + omega * R
#'   dE/dt =  beta * S * (I + theta_H * H + theta_D * D) / N - sigma * E
#'   dI/dt =  sigma * E - lambda_I * I
#'   dH/dt =  p_H * lambda_I * I - lambda_H * H
#'   dR/dt =  p_R * lambda_I * I + p_RH * lambda_H * H - omega * R
#'   dD/dt =  p_D * lambda_I * I + p_DH * lambda_H * H - lambda_B * D
#'
#' where p_R = 1 - p_H - p_D and p_RH = 1 - p_DH are derived internally.
#'
#' @param beta   Transmission rate (day^-1)
#' @param sigma  Incubation rate (day^-1), mean incubation = 1/sigma
#' @param lambda_I  Overall removal rate from I (day^-1)
#' @param p_H       Probability of hospitalisation given removal from I
#' @param p_D       Probability of death given removal from I
#' @param lambda_H  Overall removal rate from H (day^-1)
#' @param p_DH      Probability of death given removal from H
#' @param theta_H   Relative transmission while hospitalised
#' @param theta_D   Relative transmission from infectious corpses
#' @param lambda_B  Burial rate (day^-1)
#' @param omega     Immunity waning rate (day^-1)
#' @param N0        Total population
#' @param S0        Initial susceptible
#' @param E0        Initial exposed
#' @param I0        Initial infectious
#' @param H0        Initial hospitalised
#' @param R0        Initial recovered
#' @param D0        Initial deceased
#' @param t         Time span
#' @return          Data frame with columns t, S, E, I, H, R, D
#' @export
ebola_odin <- function(beta = 0.3,
                       sigma = 1/8,
                       lambda_I = 1/10,
                       p_H = 0.3,
                       p_D = 0.4,
                       lambda_H = 1/7,
                       p_DH = 0.15,
                       theta_H = 0.05,
                       theta_D = 0.02,
                       lambda_B = 1/5,
                       omega = 0,
                       N0 = 1e6,
                       S0 = NA,
                       E0 = 0,
                       I0 = 100,
                       H0 = 0,
                       R0 = 0,
                       D0 = 0,
                       t = 0:200) {
  
  if (is.na(S0)) {
    S0 <- N0 - I0 - E0 - H0 - R0 - D0
  }
  
  model <- odin({ ## had to change spacing slightly due to new lines being considered new equations
    deriv(S) <- -beta * S * (I + theta_H * H + theta_D * D) / N + omega * R
    
    deriv(E) <-  beta * S * (I + theta_H * H + theta_D * D) / N - sigma * E
    
    deriv(I) <-  sigma * E - lambda_I * I
    
    deriv(H) <-  p_H * lambda_I * I - lambda_H * H
    
    deriv(R) <-  p_R * lambda_I * I + p_RH * lambda_H * H - omega * R
    
    deriv(D) <-  p_D * lambda_I * I + p_DH * lambda_H * H - lambda_B * D
    
    initial(S) <- S0
    initial(E) <- E0
    initial(I) <- I0
    initial(H) <- H0
    initial(R) <- R0
    initial(D) <- D0
    
    N <- S0 + I0 + E0 + H0 + R0 + D0
    
    # Derived probabilities (known to odin)
    p_R  <- 1 - p_H - p_D
    p_RH <- 1 - p_DH
    
    beta       <- user(0.3)
    sigma      <- user(1/8)
    lambda_I   <- user(1/10)
    p_H        <- user(0.3)
    p_D        <- user(0.4)
    lambda_H   <- user(1/7)
    p_DH       <- user(0.15)
    theta_H    <- user(0.05)
    theta_D    <- user(0.02)
    lambda_B   <- user(1/5)
    omega      <- user(0)
    S0         <- user(999900)
    E0         <- user(0)
    I0         <- user(100)
    H0         <- user(0)
    R0         <- user(0)
    D0         <- user(0)
  })
  
  instance <- model$new(beta = beta, sigma = sigma,
                        lambda_I = lambda_I, p_H = p_H, p_D = p_D,
                        lambda_H = lambda_H, p_DH = p_DH,
                        theta_H = theta_H, theta_D = theta_D,
                        lambda_B = lambda_B, omega = omega,
                        S0 = S0, E0 = E0, I0 = I0, H0 = H0,
                        R0 = R0, D0 = D0)
  as.data.frame(instance$run(t = t))
}


#' Basic reproduction number
#'
#' Derivation: a single infectious individual in a fully susceptible
#' population generates direct infections during incubation delay.
#' Of those, some get hospitalised (transmit further at reduced rate),
#' some die (transmit through unsafe burial).
#'
#' R0 = beta / sigma * (1 + p_H * theta_H / lambda_H
#'                      + p_D * theta_D / lambda_B)
#'
#' @param beta  Transmission rate
#' @param sigma Incubation rate
#' @param p_H   Probability of hospitalisation
#' @param p_D   Probability of death
#' @param lambda_H Overall removal rate from hospitalised
#' @param lambda_B Burial rate
#' @param theta_H Relative transmission while hospitalised
#' @param theta_D Relative transmission from corpses
#' @return R0
reproduction_number <- function(beta, sigma,
                                p_H, p_D,
                                lambda_H, lambda_B,
                                theta_H = 0.05, theta_D = 0.02) {
  
  beta / sigma * (1 + p_H * theta_H / lambda_H + p_D * theta_D / lambda_B)
}


#' Ebola demonstration with plotting
#'
#' @param ... Passed to [ebola_odin()]
#' @export
demo_ebola <- function(beta = 0.3, sigma = 1/8, lambda_I = 1/10,
                       p_H = 0.3, p_D = 0.4, lambda_H = 1/7, p_DH = 0.15,
                       theta_H = 0.05, theta_D = 0.02, lambda_B = 1/5,
                       omega = 0, S0 = 999900, E0 = 0, I0 = 100,
                       H0 = 0, R0 = 0, D0 = 0, t = 0:200) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    stop("Requires ggplot2 and tidyr packages")
  }
  
  library(ggplot2)
  library(tidyr)
  
  res <- ebola_odin(beta = beta, sigma = sigma, lambda_I = lambda_I,
                    p_H = p_H, p_D = p_D, lambda_H = lambda_H, p_DH = p_DH,
                    theta_H = theta_H, theta_D = theta_D, lambda_B = lambda_B,
                    omega = omega, S0 = S0, E0 = E0, I0 = I0, H0 = H0,
                    R0 = R0, D0 = D0, t = t)
  
  res <- res %>%
    tidyr::pivot_longer(-t, names_to = "compartment", values_to = "count")
  
  ggplot(res, aes(x = t, y = count, colour = compartment)) +
    geom_line(size = 1.2) +
    labs(title = "Ebola SEIR-HPD Model",
         subtitle = paste0("R0 = ", round(reproduction_number(beta, sigma,
                                                              p_H, p_D,
                                                              lambda_H,
                                                              lambda_B,
                                                              theta_H,
                                                              theta_D), 2),
                           "\nMean incubation = ", round(1/sigma, 1),
                           " days, Mean time-to-removal (I) = ",
                           round(1/lambda_I, 1), " days"),
         x = "Time (days)", y = "Population") +
    scale_colour_brewer(palette = "Set2") +
    theme_minimal()
}
demo_ebola()



