## PROMPT: Generate code for an SIR compartmental model using the classes SimulateOde, Transition, TransitionType from PyGOM 
## INTERFACE: VS Code

## PROMPT: Add code which creates a plot to view the outbreak dynamics 

"""SIR compartmental model implemented with PyGOM.

This script builds a simple susceptible-infected-recovered model using
PyGOM's SimulateOde, Transition, and TransitionType classes.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

from pygom import SimulateOde, Transition, TransitionType


def build_sir_model() -> SimulateOde:
    """Create a PyGOM SIR model object."""
    state = ["S", "I", "R"]
    param_list = ["beta", "gamma", "N"]
    transition = [
        Transition(
            origin="S",
            destination="I",
            equation="beta*S*I/N",
            transition_type=TransitionType.T,
        ),
        Transition(
            origin="I",
            destination="R",
            equation="gamma*I",
            transition_type=TransitionType.T,
        ),
    ]
    return SimulateOde(state, param_list, transition=transition)


def simulate_sir(
    days: float = 200.0,
    dt: float = 0.1,
    beta: float = 0.3,
    gamma: float = 0.1,
    population: float = 1000.0,
    initial_infected: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Simulate the SIR model over time and return the time and state arrays."""
    model = build_sir_model()
    params = {"beta": beta, "gamma": gamma, "N": population}
    model.parameters = params

    times = np.linspace(0.0, float(days), int(float(days) / float(dt)) + 1)
    initial_state = [population - initial_infected, initial_infected, 0.0]
    model.initial_values = (initial_state, times[0])

    solution = np.asarray(model.integrate(times[1:]))
    if solution.ndim == 1:
        solution = solution.reshape(1, -1)

    if solution.shape[0] == len(times):
        solution = solution[1:]
    elif solution.shape[0] != len(times) - 1:
        if solution.shape[0] > len(times) - 1:
            solution = solution[: len(times) - 1]
        else:
            pad = np.repeat(solution[-1:][None, :], len(times) - 1 - solution.shape[0], axis=0)
            solution = np.vstack([solution, pad])

    states = np.vstack([np.array(initial_state, dtype=float)[None, :], solution])

    return times, states

def plot_results(times: np.ndarray, states: np.ndarray) -> None:
    """Plot the SIR outbreak dynamics for the susceptible, infected, and recovered compartments."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping plot output.")
        return

    labels = ["S", "I", "R"]
    colors = ["tab:blue", "tab:red", "tab:green"]

    plt.figure(figsize=(8, 4.5))
    for idx, label in enumerate(labels):
        plt.plot(times, states[:, idx], label=label, color=colors[idx], linewidth=1.8)

    plt.xlabel("Time (days)")
    plt.ylabel("Population")
    plt.title("SIR outbreak dynamics")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig("sir_outbreak_dynamics.png", dpi=200, bbox_inches="tight")
    plt.show()


def main() -> None:
    times, states = simulate_sir()
    labels = ["S", "I", "R"]
    final_state = dict(zip(labels, states[-1]))

    print("SIR model simulation complete.")
    print("Final state:")
    for name, value in final_state.items():
        print(f"  {name}: {value:.2f}")

    infected = states[:, 1]
    peak_day = times[np.argmax(infected)]
    peak_infected = infected.max()
    print(f"Peak infected: {peak_infected:.2f} on day {peak_day:.1f}")

    plot_results(times, states)


if __name__ == "__main__":
    main()
