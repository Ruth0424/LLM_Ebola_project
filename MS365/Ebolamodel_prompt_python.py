## PROMPT: Code up a compartmental model for the disease spread of Ebola using python
## INTERFACE: Microsoft 365
## MODEL: GPT 5.5 Think Deeper


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def ebola_model(t, y, params):
    """
    Ebola compartmental model.

    Compartments:
        S : susceptible
        E : exposed, infected but not yet infectious
        I : infectious in the community
        H : hospitalized / isolated infectious individuals
        F : deceased but not yet safely buried
        R : recovered
        D : safely buried / removed deaths

    Parameters:
        beta_I : transmission rate from community infectious individuals
        beta_H : transmission rate from hospitalized individuals
        beta_F : transmission rate from unsafe funeral/burial contact
        sigma  : rate from exposed to infectious, 1/incubation period
        gamma_I: recovery rate from community infection
        gamma_H: recovery rate from hospitalization
        theta  : hospitalization rate
        delta_I: death rate from community infection
        delta_H: death rate from hospitalization
        gamma_F: safe burial rate
    """

    S, E, I, H, F, R, D = y

    beta_I = params["beta_I"]
    beta_H = params["beta_H"]
    beta_F = params["beta_F"]
    sigma = params["sigma"]
    gamma_I = params["gamma_I"]
    gamma_H = params["gamma_H"]
    theta = params["theta"]
    delta_I = params["delta_I"]
    delta_H = params["delta_H"]
    gamma_F = params["gamma_F"]

    N = S + E + I + H + R  # living population

    # Force of infection
    lambda_t = (beta_I * I + beta_H * H + beta_F * F) / N

    # ODE system
    dSdt = -lambda_t * S

    dEdt = lambda_t * S - sigma * E

    dIdt = sigma * E - gamma_I * I - theta * I - delta_I * I

    dHdt = theta * I - gamma_H * H - delta_H * H

    dFdt = delta_I * I + delta_H * H - gamma_F * F

    dRdt = gamma_I * I + gamma_H * H

    dDdt = gamma_F * F

    return [dSdt, dEdt, dIdt, dHdt, dFdt, dRdt, dDdt]


# -----------------------------
# Baseline parameters
# -----------------------------

params = {
    # Transmission rates
    "beta_I": 0.35,    # community transmission
    "beta_H": 0.10,    # hospital transmission
    "beta_F": 0.50,    # unsafe funeral transmission

    # Progression rate from exposed to infectious
    # Ebola incubation is commonly 2-21 days, with average onset around 8-10 days.
    # Here we use 1/9 per day.
    "sigma": 1 / 9,

    # Recovery rates
    "gamma_I": 1 / 10,  # recovery from community infection
    "gamma_H": 1 / 12,  # recovery from hospital

    # Hospitalization rate
    "theta": 1 / 5,

    # Death rates
    "delta_I": 1 / 12,
    "delta_H": 1 / 14,

    # Safe burial rate
    "gamma_F": 1 / 2
}


# -----------------------------
# Initial conditions
# -----------------------------

N0 = 100_000

E0 = 10
I0 = 5
H0 = 0
F0 = 0
R0 = 0
D0 = 0
S0 = N0 - E0 - I0 - H0 - F0 - R0 - D0

y0 = [S0, E0, I0, H0, F0, R0, D0]


# -----------------------------
# Simulation setup
# -----------------------------

t_start = 0
t_end = 400
t_eval = np.linspace(t_start, t_end, 1000)


solution = solve_ivp(
    fun=lambda t, y: ebola_model(t, y, params),
    t_span=(t_start, t_end),
    y0=y0,
    t_eval=t_eval,
    method="RK45"
)


# -----------------------------
# Extract results
# -----------------------------

S, E, I, H, F, R, D = solution.y
t = solution.t


# -----------------------------
# Plot results
# -----------------------------

plt.figure(figsize=(12, 7))

plt.plot(t, S, label="Susceptible")
plt.plot(t, E, label="Exposed")
plt.plot(t, I, label="Infectious community")
plt.plot(t, H, label="Hospitalized")
plt.plot(t, F, label="Unsafe burial")
plt.plot(t, R, label="Recovered")
plt.plot(t, D, label="Safely buried deaths")

plt.xlabel("Time / days")
plt.ylabel("Number of people")
plt.title("Compartmental Model of Ebola Spread")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# -----------------------------
# Print useful summary statistics
# -----------------------------

peak_I = np.max(I)
peak_I_day = t[np.argmax(I)]

peak_H = np.max(H)
peak_H_day = t[np.argmax(H)]

final_size = R[-1] + D[-1]

print(f"Peak community infections: {peak_I:.0f} on day {peak_I_day:.1f}")
print(f"Peak hospitalizations: {peak_H:.0f} on day {peak_H_day:.1f}")
print(f"Final outbreak size, recovered + deaths: {final_size:.0f}")
print(f"Final recovered: {R[-1]:.0f}")
print(f"Final safely buried deaths: {D[-1]:.0f}")
