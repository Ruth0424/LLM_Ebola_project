
## PROMPT: Provide Python code for an SIR compartmental model
## INTERFACE: Isambard vllm

"""SIR Compartmental Model — Susceptible → Infected → Recovered."""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def sir_model(t: float, y: list[float], beta: float, gamma: float) -> list[float]:
    """Right-hand side of the SIR ODEs.

       dS/dt = -β · S · I / N
       dI/dt =  β · S · I / N - γ · I
       dR/dt =  γ · I
       """
    S, I, R = y
    N = S + I + R
    dSdt = -beta * S * I / N
    dIdt =  beta * S * I / N - gamma * I
    dRdt =  gamma * I
    return [dSdt, dIdt, dRdt]


def simulate_sir(
    S0: float = 990_000,
    I0: float = 1_000,
    R0: float = 0.0,
    beta: float = 0.3,
    gamma: float = 0.1,
    t_max: float = 100.0,
    steps: int = 1000) -> tuple:
    """Simulate the SIR model and return time + compartments.

    Returns (t, S, I, R).
    """
    N = S0 + I0 + R0
    y0 = [S0, I0, R0]
    t_span = (0, t_max)
    t_eval = np.linspace(t_span[0], t_span[1], steps)
    sol = solve_ivp(sir_model, t_span, y0, args=(beta, gamma), t_eval=t_eval)
    return sol.t, sol.y[0], sol.y[1], sol.y[2]


def plot_sir(t, S, I, R, beta: float, gamma: float, title: str = "SIR Model"):
    """Plot the SIR trajectories and compute key statistics."""
    R0 = beta / gamma  # Basic reproduction number

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t, S, label="Susceptible", linewidth=2)
    ax.plot(t, I, label="Infected", linewidth=2)
    ax.plot(t, R, label="Recovered", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel("Population")
    ax.set_title(f"{title}\nβ={beta}, γ={gamma}, R₀={R0:.2f}")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    t, S, I, R = simulate_sir()
    plot_sir(t, S, I, R, beta=0.3, gamma=0.1)
