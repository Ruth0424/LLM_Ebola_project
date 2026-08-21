
## PROMPT: Code a compartmental model for the disease spread of Ebola using odin
# attached "Modeling Spatial Invasion of Ebola in West Africa" by D'Silva and Eisenberg
## INTERFACE: Microsoft 365
## MODEL: GPT 5.5 Think Deeper

library(odin)
library(deSolve)
library(ggplot2)
library(reshape2)

ebola_model <- odin({
  
  ## Core transmission model
  deriv(S)  <- -lambda * S
  deriv(E)  <-  lambda * S - alpha * E
  deriv(I1) <-  alpha * E - gamma * I1 - r1 * I1
  deriv(I2) <-  gamma * I1 - delta * I2 - r2 * I2
  deriv(F)  <-  delta * I2 - delta2 * F
  deriv(R)  <-  r1 * I1 + r2 * I2
  
  ## Cumulative reported outcomes
  deriv(IC) <- k_norm * alpha * E
  deriv(DC) <- k_norm * delta * I2
  
  ## Force of infection
  lambda <- beta1 * I1 + beta2 * I2 + betaF * F
  
  ## Initial conditions
  initial(S)  <- S0
  initial(E)  <- E0
  initial(I1) <- I10
  initial(I2) <- I20
  initial(F)  <- F0
  initial(R)  <- R0
  initial(IC) <- IC0
  initial(DC) <- DC0
  
  ## User parameters
  S0     <- user()
  E0     <- user()
  I10    <- user()
  I20    <- user()
  F0     <- user()
  R0     <- user()
  IC0    <- user()
  DC0    <- user()
  
  beta1  <- user()
  beta2  <- user()
  betaF  <- user()
  
  alpha  <- user()
  gamma  <- user()
  delta  <- user()
  delta2 <- user()
  r1     <- user()
  r2     <- user()
  
  k_norm <- user()
})

mod <- ebola_model$new(
  S0 = 1e6 - 20,
  E0 = 10,
  I10 = 5,
  I20 = 5,
  F0 = 0,
  R0 = 0,
  IC0 = 10,
  DC0 = 0,
  
  ## Transmission coefficients
  beta1 = 2.0e-7,
  beta2 = 4.0e-7,
  betaF = 4.0e-7,
  
  ## Disease progression rates
  alpha  = 1 / 9,   # exposed to infectious
  gamma  = 1 / 5,   # I1 to I2
  delta  = 1 / 4,   # I2 to funeral/death
  delta2 = 1 / 2,   # burial rate
  r1     = 1 / 10,  # recovery from I1
  r2     = 1 / 5,   # recovery from I2
  
  ## Reporting / population-at-risk scaling
  k_norm = 1
)

times <- seq(0, 200, by = 1)
out <- as.data.frame(mod$run(times))

head(out)

out_long <- reshape2::melt(
  out,
  id.vars = "t",
  measure.vars = c("S", "E", "I1", "I2", "F", "R", "IC", "DC")
)

ggplot(out_long, aes(x = t, y = value, colour = variable)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time / days",
    y = "Number of individuals",
    colour = "Compartment",
    title = "Compartmental Ebola model"
  ) +
  theme_minimal() # + facet_wrap(~ variable, scales = "free_y")

ebola_gravity_model <- odin({
  ## State equations for each patch
  deriv(S[])  <- -lambda[i] * S[i]
  deriv(E[])  <-  lambda[i] * S[i] - alpha * E[i]
  deriv(I1[]) <-  alpha * E[i] - gamma[i] * I1[i] - r1[i] * I1[i]
  deriv(I2[]) <-  gamma[i] * I1[i] - delta * I2[i] - r2[i] * I2[i]
  deriv(F[])  <-  delta * I2[i] - delta2 * F[i]
  deriv(R[])  <-  r1[i] * I1[i] + r2[i] * I2[i]
  
  ## Cumulative incidence and deaths
  deriv(IC[]) <- k_norm[i] * alpha * E[i]
  deriv(DC[]) <- k_norm[i] * delta * I2[i]
  
  ## Total force of infection into patch i
  lambda[i] <- local_lambda[i] + spatial_lambda[i]
  
  ## Local transmission within patch i
  local_lambda[i] <- beta1[i] * I1[i] +
    beta2[i] * I2[i] +
    betaF[i] * F[i]
  
  ## Spatial transmission from all other patches j into patch i
  spatial_lambda[i] <- sum(theta[i, ] *
                             (beta1[i] * I1[] +
                                beta2[i] * I2[] +
                                betaF[i] * F[]))
  
  ## Gravity term
  theta[i, j] <- if (i == j) 0 else
    kappa[i] * rho[i] * rho[j] / distance[i, j]^iota
  
  ## Initial conditions
  initial(S[])  <- S0[i]
  initial(E[])  <- E0[i]
  initial(I1[]) <- I10[i]
  initial(I2[]) <- I20[i]
  initial(F[])  <- F0[i]
  initial(R[])  <- R0[i]
  initial(IC[]) <- IC0[i]
  initial(DC[]) <- DC0[i]
  
  ## Dimensions
  dim(S)  <- n_patch
  dim(E)  <- n_patch
  dim(I1) <- n_patch
  dim(I2) <- n_patch
  dim(F)  <- n_patch
  dim(R)  <- n_patch
  dim(IC) <- n_patch
  dim(DC) <- n_patch
  
  dim(lambda)         <- n_patch
  dim(local_lambda)   <- n_patch
  dim(spatial_lambda) <- n_patch
  
  dim(theta)    <- c(n_patch, n_patch)
  dim(distance) <- c(n_patch, n_patch)
  
  dim(S0)  <- n_patch
  dim(E0)  <- n_patch
  dim(I10) <- n_patch
  dim(I20) <- n_patch
  dim(F0)  <- n_patch
  dim(R0)  <- n_patch
  dim(IC0) <- n_patch
  dim(DC0) <- n_patch
  
  dim(beta1)  <- n_patch
  dim(beta2)  <- n_patch
  dim(betaF)  <- n_patch
  dim(gamma)  <- n_patch
  dim(r1)     <- n_patch
  dim(r2)     <- n_patch
  dim(kappa)  <- n_patch
  dim(rho)    <- n_patch
  dim(k_norm) <- n_patch
  
  ## User inputs
  n_patch <- user()
  
  S0  <- user()
  E0  <- user()
  I10 <- user()
  I20 <- user()
  F0  <- user()
  R0  <- user()
  IC0 <- user()
  DC0 <- user()
  
  beta1  <- user()
  beta2  <- user()
  betaF  <- user()
  
  alpha  <- user()
  gamma  <- user()
  delta  <- user()
  delta2 <- user()
  r1     <- user()
  r2     <- user()
  
  kappa    <- user()
  rho      <- user()
  distance <- user()
  iota     <- user()
  k_norm   <- user()
})

n_patch <- 3

rho <- c(
  Guinea = 12e6,
  SierraLeone = 7e6,
  Liberia = 4.5e6
)

distance <- matrix(
  c(
    1,    340,  590,
    340,  1,    520,
    590,  520,  1
  ),
  nrow = 3,
  byrow = TRUE
)

mod_g <- ebola_gravity_model$new(
  n_patch = n_patch,
  
  S0  = rho - c(20, 5, 5),
  E0  = c(10, 2, 2),
  I10 = c(5, 1, 1),
  I20 = c(5, 1, 1),
  F0  = c(0, 0, 0),
  R0  = c(0, 0, 0),
  IC0 = c(10, 2, 2),
  DC0 = c(0, 0, 0),
  
  beta1 = c(2.0e-8, 2.5e-8, 3.0e-8),
  beta2 = c(4.0e-8, 5.0e-8, 6.0e-8),
  betaF = c(4.0e-8, 5.0e-8, 6.0e-8),
  
  alpha  = 1 / 9,
  gamma  = c(1 / 5, 1 / 5, 1 / 5),
  delta  = 1 / 4,
  delta2 = 1 / 2,
  r1     = c(1 / 10, 1 / 10, 1 / 10),
  r2     = c(1 / 5, 1 / 5, 1 / 5),
  
  kappa = c(1e-17, 1e-17, 1e-17),
  rho = rho,
  distance = distance,
  iota = 2,
  k_norm = c(1, 1, 1)
)

times <- seq(0, 200, by = 1)
out_g <- as.data.frame(mod_g$run(times))

head(out_g)