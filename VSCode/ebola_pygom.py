## PROMPT: Generate code for a compartmental model to model an outbreak of Ebola using the classes SimulateOde, Transition, TransitionType from PyGOM 
## INTERFACE: VS Code

"""Ebola outbreak compartmental model implemented with PyGOM.

The model uses a seven-compartment structure:
- S: susceptible
- E: exposed
- I: infectious in the community
- H: hospitalized
- R: recovered
- D: dead bodies awaiting burial
- B: buried bodies

This version uses the PyGOM classes SimulateOde, Transition, and TransitionType.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from pygom import SimulateOde, Transition, TransitionType


def build_ebola_model() -> SimulateOde:
    """Build an Ebola outbreak model using PyGOM transitions."""
    state_names = ["S", "E", "I", "H", "R", "D", "B"]
    param_names = [
        "beta_i",
        "beta_h",
        "beta_d",
        "sigma",
        "rho",
        "gamma_i",
        "gamma_h",
        "mu_i",
        "mu_h",
        "psi",
        "N",
    ]

    ode_transitions = [
        Transition(
            origin="S",
            equation="-beta_i*S*I/N - beta_h*S*H/N - beta_d*S*D/N",
            transition_type=TransitionType.ODE,
        ),
        Transition(
            origin="E",
            equation="beta_i*S*I/N + beta_h*S*H/N + beta_d*S*D/N - sigma*E",
            transition_type=TransitionType.ODE,
        ),
        Transition(
            origin="I",
            equation="sigma*E - rho*I - gamma_i*I - mu_i*I",
            transition_type=TransitionType.ODE,
        ),
        Transition(
            origin="H",
            equation="rho*I - gamma_h*H - mu_h*H",
            transition_type=TransitionType.ODE,
        ),
        Transition(
            origin="R",
            equation="gamma_i*I + gamma_h*H",
            transition_type=TransitionType.ODE,
        ),
        Transition(
            origin="D",
            equation="mu_i*I + mu_h*H - psi*D",
            transition_type=TransitionType.ODE,
        ),
        Transition(
            origin="B",
            equation="psi*D",
            transition_type=TransitionType.ODE,
        ),
    ]

    return SimulateOde(state_names, param_names, ode=ode_transitions)


def simulate_model(
    params: Dict[str, float],
    initial_state: List[float],
    days: int = 180,
    dt: float = 0.1,
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve the Ebola model over time with PyGOM."""
    model = build_ebola_model()
    model.parameters = params

    times = np.arange(0.0, days + dt, dt)
    model.initial_values = (initial_state, times[0])

    solution = np.asarray(model.integrate(times[1:]))
    if solution.ndim == 1:
        solution = solution.reshape(1, -1)

    if solution.shape[0] == len(times) + 1:
        solution = solution[:-1]
    elif solution.shape[0] == len(times) - 1:
        solution = np.vstack([np.asarray(initial_state, dtype=float)[None, :], solution])
    elif solution.shape[0] != len(times):
        if solution.shape[0] > len(times):
            solution = solution[: len(times)]
        elif solution.shape[0] < len(times):
            pad = np.repeat(solution[-1:][None, :], len(times) - solution.shape[0], axis=0)
            solution = np.vstack([solution, pad])

    states = np.asarray(solution, dtype=float)
    return times, states


def plot_results(times: np.ndarray, states: np.ndarray) -> None:
    """Plot the model compartments over time."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover - optional dependency
        print("matplotlib is not installed; skipping plot output.")
        return

    labels = ["S", "E", "I", "H", "R", "D", "B"]
    colors = ["tab:blue", "tab:orange", "tab:red", "tab:purple", "tab:green", "black", "gray"]

    plt.figure(figsize=(10, 6))
    for idx, label in enumerate(labels):
        plt.plot(times, states[:, idx], label=label, color=colors[idx], linewidth=1.8)

    plt.xlabel("Time (days)")
    plt.ylabel("Population")
    plt.title("Ebola compartmental model")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.show()


def main() -> None:
    population = 1000.0
    params = {
        "beta_i": 0.5,
        "beta_h": 0.2,
        "beta_d": 0.2,
        "sigma": 1/7.65,
        "rho": 1/6,
        "gamma_i": 163/2880,
        "gamma_h": 1/13,
        "mu_i": 25/576,
        "mu_h": 6/91,
        "psi": 0.5,
        "N": population,
    }

    initial_state = [population - 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]

    times, states = simulate_model(params, initial_state, days=180, dt=0.1)

    labels = ["S", "E", "I", "H", "R", "D", "B"]
    final_state = dict(zip(labels, states[-1]))

    print("PyGOM Ebola model simulation complete.")
    print("\nFinal state after 180 days:")
    for name, value in final_state.items():
        print(f"  {name}: {value:,}") # removed .0f to show decimal values

    infectious_values = states[:, 2]
    peak_infectious = np.max(infectious_values)
    peak_day = times[np.argmax(infectious_values)]
    print(f"\nPeak infectious count: {peak_infectious:,.0f} on day {peak_day:.1f}")

    plot_results(times, states)


if __name__ == "__main__":
    main()
