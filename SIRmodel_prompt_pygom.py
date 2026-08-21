## PROMPT: Code an SIR model using PyGOM
## INTERFACE: Microsoft 365
## MODEL: GPT 5.5 Think Deeper

from pygom import common_models
import numpy as np
import matplotlib.pyplot as plt
import math

# -----------------------------
# Parameters
# -----------------------------

N = 10_000          # total population
I0 = 1              # initial infected
R0_initial = 0      # initial recovered
S0 = N - I0 - R0_initial

gamma = 1 / 4       # recovery rate: average infectious period = 4 days
R_0 = 1.3           # basic reproduction number
beta = R_0 * gamma  # transmission rate

params = {
    "beta": beta,
    "gamma": gamma,
    "N": N
}

# -----------------------------
# Set up PyGOM SIR model
# -----------------------------

ode = common_models.SIR(params)

# Initial conditions: [S, I, R]
x0 = [S0, I0, R0_initial]

# Time grid
tmax = 365
dt = 0.1
n_timestep = math.ceil(tmax / dt)
t = np.linspace(0, tmax, n_timestep)

# PyGOM requires initial values and initial time
ode.initial_values = (x0, t[0])

# -----------------------------
# Solve deterministic model
# -----------------------------

solution = ode.solve_determ(t[1:])

# Extract compartments
S = solution[:, 0]
I = solution[:, 1]
R = solution[:, 2]

# -----------------------------
# Plot results
# -----------------------------

plt.figure(figsize=(10, 6))

plt.plot(t, S, label="Susceptible")
plt.plot(t, I, label="Infected")
plt.plot(t, R, label="Recovered")

plt.xlabel("Time")
plt.ylabel("Number of people")
plt.title("Deterministic SIR Model using PyGOM")
plt.legend()
plt.tight_layout()
plt.show()

## PROMPT: Extend this SIR model to model key aspects of an Ebola outbreak

from pygom import SimulateOde, Transition, TransitionType
import numpy as np
import matplotlib.pyplot as plt
import math

# ------------------------------------------------------------
# States and parameters
# ------------------------------------------------------------

state_list = ["S", "E", "I", "H", "F", "R", "D"]

param_list = [
    "N",
    "beta_I",      # transmission from community infectious cases
    "beta_H",      # transmission from hospitalised cases
    "beta_F",      # transmission from funeral/body contact
    "sigma",       # rate from exposed to infectious
    "theta",       # hospitalisation/isolation rate
    "gamma_I",     # recovery rate from community infection
    "gamma_H",     # recovery rate from hospital
    "delta_I",     # death rate from community infection
    "delta_H",     # death rate from hospital
    "phi"          # safe burial rate
]

# ------------------------------------------------------------
# Force of infection
# ------------------------------------------------------------
#
# lambda = (beta_I I + beta_H H + beta_F F) / N
#
# New infections move S -> E.
#
# Ebola transmission here is split into:
#   - community transmission from I
#   - healthcare/hospital transmission from H
#   - funeral/body-contact transmission from F
#

transitions = [

    # Infection: susceptible people become exposed
    Transition(
        origin="S",
        destination="E",
        equation="S*(beta_I*I + beta_H*H + beta_F*F)/N",
        transition_type=TransitionType.T
    ),

    # End of incubation period: exposed become infectious
    Transition(
        origin="E",
        destination="I",
        equation="sigma*E",
        transition_type=TransitionType.T
    ),

    # Infectious community cases are hospitalised/isolated
    Transition(
        origin="I",
        destination="H",
        equation="theta*I",
        transition_type=TransitionType.T
    ),

    # Community infectious cases recover
    Transition(
        origin="I",
        destination="R",
        equation="gamma_I*I",
        transition_type=TransitionType.T
    ),

    # Community infectious cases die
    Transition(
        origin="I",
        destination="F",
        equation="delta_I*I",
        transition_type=TransitionType.T
    ),

    # Hospitalised cases recover
    Transition(
        origin="H",
        destination="R",
        equation="gamma_H*H",
        transition_type=TransitionType.T
    ),

    # Hospitalised cases die
    Transition(
        origin="H",
        destination="F",
        equation="delta_H*H",
        transition_type=TransitionType.T
    ),

    # Safe burial removes funeral transmission risk
    Transition(
        origin="F",
        destination="D",
        equation="phi*F",
        transition_type=TransitionType.T
    )
]

# ------------------------------------------------------------
# Build PyGOM model
# ------------------------------------------------------------

model = SimulateOde(
    state=state_list,
    param=param_list,
    transition=transitions
)

# ------------------------------------------------------------
# Example parameter values chosen to show clear epidemic dynamics
# within the first 250 days
# ------------------------------------------------------------

N = 1001

params = [
    ("N", N),

    # Transmission rates
    ("beta_I", 0.5),     # community transmission
    ("beta_H", 0.2),     # hospital transmission, lower than community
    ("beta_F", 0.2),     # funeral/body-contact transmission

    # Disease progression
    ("sigma", 1 / 7.65),   # mean incubation period around 7 days

    # Movement from community infection to hospital/isolation
    ("theta", 1 / 6.0),   # mean time to hospitalisation around 6 days

    # Community outcomes
    ("gamma_I", 163 / 2880),   # community recovery rate
    ("delta_I", 25 / 576),   # community death rate

    # Hospital outcomes
    ("gamma_H", 1 / 13),   # hospital recovery rate
    ("delta_H", 6 / 91),   # hospital death rate

    # Burial process
    ("phi", 1 / 2)      # mean time to safe burial around 3 days
]

model.parameters = params

# ------------------------------------------------------------
# Initial conditions
# ------------------------------------------------------------

E0 = 0
I0 = 1
H0 = 0
F0 = 0
R0 = 0
D0 = 0
S0 = N - E0 - I0 - H0 - F0 - R0 - D0

x0 = [S0, E0, I0, H0, F0, R0, D0]

# ------------------------------------------------------------
# Solve model
# ------------------------------------------------------------
# ------------------------------------------------------------
# Time grid
# ------------------------------------------------------------

tmax = 150
dt = 0.1

t = np.arange(0, tmax + dt, dt)

model.initial_values = (x0, t[0])

# Ask PyGOM to solve over the full time grid
solution = model.solve_determ(t)

# ------------------------------------------------------------
# Make sure time and solution lengths match
# ------------------------------------------------------------

print("t shape:", t.shape)
print("solution shape:", solution.shape)

# If PyGOM returns one extra row, trim it
if solution.shape[0] > len(t):
    solution = solution[:len(t), :]

# If PyGOM returns one fewer row, prepend initial condition
elif solution.shape[0] == len(t) - 1:
    solution = np.vstack([x0, solution])

# Final safety check
assert solution.shape[0] == len(t), (
    f"Mismatch: len(t)={len(t)}, solution rows={solution.shape[0]}"
)

# ------------------------------------------------------------
# Extract compartments
# ------------------------------------------------------------

S = solution[:, 0]
E = solution[:, 1]
I = solution[:, 2]
H = solution[:, 3]
F = solution[:, 4]
R = solution[:, 5]
D = solution[:, 6]
print(solution[-1])
# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

plt.figure(figsize=(11, 7))

plt.plot(t, S, label="Susceptible")
plt.plot(t, E, label="Exposed")
plt.plot(t, I, label="Infectious community")
plt.plot(t, H, label="Hospitalised / isolated")
plt.plot(t, F, label="Funeral infectious")
plt.plot(t, R, label="Recovered")
plt.plot(t, D, label="Safely buried deaths")

plt.xlabel("Time / days")
plt.ylabel("Number of people")
plt.title("Ebola-style SEIHFR Model using PyGOM")
plt.legend()
plt.tight_layout()
plt.show()