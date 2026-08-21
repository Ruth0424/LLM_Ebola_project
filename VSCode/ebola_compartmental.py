## PROMPT: Code up a compartmental model for the disease spread of Ebola using Python in this file 
## INTERFACE: VS Code

"""Ebola compartmental model in Python.

This script implements an Ebola-inspired SEIHRD model with:
- S: susceptible
- E: exposed
- I: infectious in the community
- H: hospitalized
- R: recovered
- D: dead bodies awaiting burial
- B: buried bodies

The model uses a simple ordinary differential equation system and a basic
Euler solver so it can run without external scientific packages.
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

try:
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover - optional dependency
    plt = None


def ebola_model(t: float, y: List[float], params: Dict[str, float]) -> List[float]:
    """Return derivatives for the Ebola compartmental model."""
    S, E, I, H, R, D, B = y

    N = S + E + I + H + R + D + B
    if N <= 0:
        raise ValueError("Population size became non-positive.")

    beta_i = params["beta_i"]
    beta_h = params["beta_h"]
    beta_d = params["beta_d"]
    sigma = params["sigma"]
    rho = params["rho"]
    gamma_i = params["gamma_i"]
    gamma_h = params["gamma_h"]
    mu_i = params["mu_i"]
    mu_h = params["mu_h"]
    psi = params["psi"]

    force_of_infection = (beta_i * I + beta_h * H + beta_d * D) / N

    dS = -force_of_infection * S
    dE = force_of_infection * S - sigma * E
    dI = sigma * E - (rho + gamma_i + mu_i) * I
    dH = rho * I - (gamma_h + mu_h) * H
    dR = gamma_i * I + gamma_h * H
    dD = mu_i * I + mu_h * H - psi * D
    dB = psi * D

    return [dS, dE, dI, dH, dR, dD, dB]


def solve_model(
    params: Dict[str, float],
    initial_state: List[float],
    days: int = 180,
    dt: float = 0.1,
) -> Tuple[List[float], List[List[float]]]:
    """Solve the ODE system using a simple explicit Euler method."""
    steps = int(math.ceil(days / dt))
    times = [i * dt for i in range(steps + 1)]
    states: List[List[float]] = [[0.0] * len(initial_state) for _ in range(steps + 1)]
    states[0] = list(initial_state)

    for i in range(steps):
        derivative = ebola_model(times[i], states[i], params)
        states[i + 1] = [
            states[i][j] + dt * derivative[j]
            for j in range(len(initial_state))
        ]

    return times, states


def plot_results(times: List[float], states: List[List[float]]) -> None:
    """Plot the model compartments over time."""
    if plt is None:
        print("matplotlib is not installed; skipping plot output.")
        return

    labels = ["S", "E", "I", "H", "R", "D", "B"]
    colors = ["tab:blue", "tab:orange", "tab:red", "tab:purple", "tab:green", "black", "gray"]

    plt.figure(figsize=(10, 6))
    for idx, label in enumerate(labels):
        values = [state[idx] for state in states]
        plt.plot(times, values, label=label, color=colors[idx], linewidth=1.8)

    plt.xlabel("Time (days)")
    plt.ylabel("Population")
    plt.title("Ebola compartmental model")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.show()


def main() -> None:
    params = {
        "beta_i": 0.45,
        "beta_h": 0.15,
        "beta_d": 0.10,
        "sigma": 0.20,
        "rho": 0.15,
        "gamma_i": 0.06,
        "gamma_h": 0.08,
        "mu_i": 0.03,
        "mu_h": 0.06,
        "psi": 0.12,
    }

    initial_state = [999_999, 0, 1, 0, 0, 0, 0]

    times, states = solve_model(params, initial_state, days=180, dt=0.1)

    labels = ["S", "E", "I", "H", "R", "D", "B"]
    final_state = dict(zip(labels, states[-1]))

    print("Final state after 180 days:")
    for name, value in final_state.items():
        print(f"  {name}: {value:,.0f}")

    print("\nPeak infectious count:")
    infectious_values = [state[2] for state in states]
    peak_infectious = max(infectious_values)
    peak_day = times[infectious_values.index(peak_infectious)]
    print(f"  {peak_infectious:,.0f} on day {peak_day:.1f}")

    plot_results(times, states)


if __name__ == "__main__":
    main()
