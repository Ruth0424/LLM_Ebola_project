## PROMPT: Code a compartmental model for the disease spread of Ebola using odin
# attached "Modeling Spatial Invasion of Ebola in West Africa" by D'Silva and Eisenberg
## INTERFACE: Microsoft 365
## MODEL: GPT 5.2 Think Deeper


library(odin)
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

ebola_model <- odin::odin({
  ## ------------------------------------------------------------
  ## Ebola compartmental model
  ##
  ## Compartments:
  ## S  : susceptible
  ## E  : exposed / latent
  ## I1 : infectious, stage 1
  ## I2 : infectious, stage 2
  ## F  : dead but not yet buried; funeral transmission
  ## R  : recovered
  ## IC : cumulative reported infections
  ## DC : cumulative reported deaths
  ##
  ## Structure based on D'Silva & Eisenberg's EVD model:
  ## S, E, I1, I2, F, R with reported cumulative cases/deaths.
  ## ------------------------------------------------------------
  
  ## Force of infection
  lambda <- beta1 * I1 + beta2 * I2 + betaF * F
  
  ## Differential equations
  deriv(S)  <- -lambda * S
  deriv(E)  <-  lambda * S - alpha * E
  deriv(I1) <-  alpha * E - gamma * I1 - r1 * I1
  deriv(I2) <-  gamma * I1 - delta * I2 - r2 * I2
  deriv(F)  <-  delta * I2 - delta2 * F
  deriv(R)  <-  r1 * I1 + r2 * I2
  
  ## Cumulative reported cases and deaths
  deriv(IC) <- k_norm * alpha * E
  deriv(DC) <- k_norm * delta * I2
  
  ## Initial conditions
  initial(S)  <- S0
  initial(E)  <- E0
  initial(I1) <- I10
  initial(I2) <- I20
  initial(F)  <- F0
  initial(R)  <- R0
  initial(IC) <- IC0
  initial(DC) <- DC0
  
  ## User-defined parameters
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

## Create model object
mod <- ebola_model$new(
  user = list(
    ## Initial conditions
    S0  = 1e5 - 20,
    E0  = 10,
    I10 = 5,
    I20 = 3,
    F0  = 2,
    R0  = 0,
    IC0 = 0,
    DC0 = 0,
    
    ## Transmission parameters
    beta1 = 1.0e-6,
    beta2 = 2.0e-6,
    betaF = 2.0e-6,
    
    ## Progression rates
    ## alpha: exposed -> infectious
    ## gamma: I1 -> I2
    ## delta: I2 -> F
    ## delta2: burial/removal from F
    ## r1, r2: recovery rates from I1 and I2
    alpha  = 1 / 9,
    gamma  = 1 / 5,
    delta  = 1 / 4,
    delta2 = 1 / 2,
    r1     = 1 / 10,
    r2     = 1 / 10,
    
    ## Reporting / normalisation factor
    k_norm = 1
  )
)

## Run simulation
times <- seq(0, 180, by = 1)
out <- as.data.frame(mod$run(times))

head(out)

out_long <- out %>%
  pivot_longer(
    cols = c(S, E, I1, I2, F, R),
    names_to = "compartment",
    values_to = "value"
  )

ggplot(out_long, aes(x = t, y = value, colour = compartment)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time, days",
    y = "Number of individuals",
    colour = "Compartment",
    title = "Ebola compartmental model"
  ) +
  theme_minimal()

out_cum <- out %>%
  select(t, IC, DC) %>%
  pivot_longer(
    cols = c(IC, DC),
    names_to = "quantity",
    values_to = "value"
  )

ggplot(out_cum, aes(x = t, y = value, colour = quantity)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Time, days",
    y = "Cumulative reported count",
    colour = "Quantity",
    title = "Reported cumulative Ebola cases and deaths"
  ) +
  theme_minimal()

# frequency-dependent version and gravity coupling options mentioned and code provided