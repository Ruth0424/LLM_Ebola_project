## PROMPT: Generate code for an SIR compartmental model using PyGOM 
## INTERFACE: VS Code

## PROMPT: Amend the code to use class DeterministicOde from the pygom package rather than Model, without the SciPy solver 

"""SIR compartmental model implemented with PyGOM when available.

This script builds a simple susceptible-infected-recovered model with:
- S: susceptible
- I: infectious
- R: recovered

It creates a PyGOM model object if the library is installed and otherwise
falls back to a direct SciPy solver so the example remains runnable.
"""

from __future__ import annotations

import warnings
from typing import Dict, List, Tuple

import numpy as np

try:
    from pygom import DeterministicOde, Transition
except Exception as exc:  # pragma: no cover - depends on environment
    DeterministicOde = None
    Transition = None
    PYGOM_IMPORT_ERROR = exc
else:
    PYGOM_IMPORT_ERROR = None


def sir_rhs(t: float, y: np.ndarray, params: Dict[str, float]) -> np.ndarray:
    """Return the SIR ODE system derivatives."""
    S, I, R = y
    N = S + I + R

    if N <= 0:
        raise ValueError("Population size became non-positive.")

    beta = params["beta"]
    gamma = params["gamma"]

    force_of_infection = beta * S * I / N

    dS = -force_of_infection
    dI = force_of_infection - gamma * I
    dR = gamma * I

    return np.array([dS, dI, dR], dtype=float)


def build_pygom_model():
    """Create a PyGOM deterministic ODE model for the SIR system."""
    if Transition is None or DeterministicOde is None:
        return None

    transitions = [
        Transition(
            origin="S",
            destination="I",
            equation="beta * S * I / (S + I + R)",
        ),
        Transition(origin="I", destination="R", equation="gamma * I"),
    ]

    candidate_kwargs = [
        {"state": ["S", "I", "R"], "param": ["beta", "gamma"], "transition": transitions},
        {"state_names": ["S", "I", "R"], "param_names": ["beta", "gamma"], "transitions": transitions},
    ]

    for kwargs in candidate_kwargs:
        try:
            return DeterministicOde(**kwargs)
        except TypeError:
            continue
        except Exception as exc:  # pragma: no cover - depends on underlying PyGOM version
            warnings.warn(f"PyGOM model construction failed: {exc}")
            return None

    warnings.warn("PyGOM deterministic model construction did not match the available API.")
    return None


def solve_model(
    model,
    params: Dict[str, float],
    initial_state: List[float],
    days: int = 200,
    dt: float = 0.1,
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve the SIR model over time using PyGOM's deterministic ODE solver."""
    if model is None:
        raise RuntimeError("A PyGOM deterministic model is required for solving.")

    t_eval = np.arange(0, days + dt, dt)
    initial_state_array = np.array(initial_state, dtype=float)

    solution = None
    for kwargs in (
        {"t": t_eval, "y0": initial_state_array, "params": params},
        {"t": t_eval, "y0": initial_state_array, "theta": params},
        {"t": t_eval, "y0": initial_state_array, "param": params},
    ):
        try:
            solution = model.solve(**kwargs)
            break
        except TypeError:
            continue

    if solution is None:
        raise RuntimeError("PyGOM deterministic solver did not accept the available argument names.")

    solution_array = np.asarray(solution)
    if solution_array.ndim == 1:
        states = solution_array.reshape(1, -1)
        return t_eval[:1], states

    if solution_array.shape[0] == len(t_eval):
        states = solution_array
    elif solution_array.shape[1] == len(t_eval):
        states = solution_array.T
    else:
        raise RuntimeError(f"Unexpected PyGOM solution shape: {solution_array.shape}")

    return t_eval, states


def plot_results(times: np.ndarray, states: np.ndarray) -> None:
    """Plot the model compartments over time."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover - optional dependency
        print("matplotlib is not installed; skipping plot output.")
        return

    labels = ["S", "I", "R"]
    colors = ["tab:blue", "tab:red", "tab:green"]

    plt.figure(figsize=(8, 4.5))
    for idx, label in enumerate(labels):
        plt.plot(times, states[:, idx], label=label, color=colors[idx], linewidth=1.8)

    plt.xlabel("Time (days)")
    plt.ylabel("Population")
    plt.title("SIR compartmental model")
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.show()


def main() -> None:
    params = {"beta": 0.30, "gamma": 0.10}
    initial_state = [999, 1, 0]

    if PYGOM_IMPORT_ERROR is not None:
        raise RuntimeError(f"PyGOM is not available in this environment: {PYGOM_IMPORT_ERROR}")

    model = build_pygom_model()
    if model is None:
        raise RuntimeError("PyGOM deterministic model construction was not successful.")

    print("PyGOM deterministic model object created successfully.")

    times, states = solve_model(model, params, initial_state, days=200, dt=0.1)

    labels = ["S", "I", "R"]
    final_state = dict(zip(labels, states[-1]))

    print("\nFinal state after 200 days:")
    for name, value in final_state.items():
        print(f"  {name}: {value:,.0f}")

    infectious_values = states[:, 1]
    peak_infectious = np.max(infectious_values)
    peak_day = times[np.argmax(infectious_values)]
    print(f"\nPeak infectious count: {peak_infectious:,.0f} on day {peak_day:.1f}")

    plot_results(times, states)


if __name__ == "__main__":
    main()
