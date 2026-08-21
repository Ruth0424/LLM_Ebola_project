
## PROMPT: Extend this code, still using Odin’s inline R syntax, to produce a 
## compartmental model better suited to a disease which has an incubation period. 
## INTERFACE: Isambard vllm

library(odin)

seir_odin <- function(beta = 0.3,
                      gamma = 0.1,
                      sigma = 1/5,
                      N0 = 1e6,
                      S0 = NA,
                      E0 = 0,
                      I0 = 1000,
                      R0 = 0,
                      t = 0:200) {
  
  if (is.na(S0)) {
    S0 <- N0 - I0 - E0 - R0
  }
  
  model <- odin({
    deriv(S) <- -beta * S * I / N
    deriv(E) <-  beta * S * I / N - sigma * E
    deriv(I) <-  sigma * E - gamma * I
    deriv(R) <-  gamma * I
    
    initial(S) <- S0
    initial(E) <- E0
    initial(I) <- I0
    initial(R) <- R0
    
    N <- S0 + I0 + E0 + R0
    
    beta  <- user(0.3)
    gamma <- user(0.1)
    sigma <- user(1/5)
    S0    <- user(999000)
    E0    <- user(0)
    I0    <- user(1000)
    R0    <- user(0)
  })
  
  instance <- model$new(beta = beta, gamma = gamma, sigma = sigma,
                        S0 = S0, E0 = E0, I0 = I0, R0 = R0)
  as.data.frame(instance$run(t = t))
}


reproduction_number <- function(beta, gamma) beta / gamma


demo_seir <- function(beta = 0.3, gamma = 0.1, sigma = 1/5,
                      S0 = 999000, E0 = 1000, I0 = 0, R0 = 0,
                      t = 0:200) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("tidyr", quietly = TRUE)) {
    stop("Requires ggplot2 and tidyr packages")
  }
  
  library(ggplot2)
  library(tidyr)
  
  res <- seir_odin(beta = beta, gamma = gamma, sigma = sigma,
                   S0 = S0, E0 = E0, I0 = I0, R0 = R0, t = t)
  
  res %>%
    tidyr::pivot_longer(-t, names_to = "compartment", values_to = "count") %>%
    ggplot(aes(x = t, y = count, colour = compartment)) +
    geom_line(size = 1.2) +
    labs(title = "SEIR Compartmental Model",
         subtitle = paste0("R0 = ", round(reproduction_number(beta, gamma), 2),
                           "  |  Mean incubation = ", round(1/sigma, 1), " days"),
         x = "Time", y = "Population") +
    scale_colour_brewer(palette = "Set1") +
    theme_minimal()
}

demo_seir()
