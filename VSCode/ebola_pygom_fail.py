## PROMPT: Generate code PyGOM in this file for a compartmental model to model an outbreak of Ebola 
## INTERFACE: VS Code

## PROMPT: Amend the code to use the classes SimulateOde, Transition, TransitionType from the pygom package rather than Model, and remove the SciPy solver 

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

try:
    from pygom import SimulateOde, Transition, TransitionType
except Exception as exc:  # pragma: no cover - depends on environment
    SimulateOde = None
    Transition = None
    TransitionType = None
    PYGOM_IMPORT_ERROR = exc
else:
    PYGOM_IMPORT_ERROR = None


def build_pygom_model():
    """Build the Ebola epidemic model with PyGOM transitions."""
    if Transition is None or TransitionType is None or SimulateOde is None:
        raise ImportError(f"PyGOM is not available: {PYGOM_IMPORT_ERROR}")

    transitions = [
        Transition(
            origin="S",
            destination="E",
            equation="beta_i * S * I / N + beta_h * S * H / N + beta_d * S * D / N",
            transition_type=TransitionType.ODE,
        ),
        Transition(origin="E", destination="I", equation="sigma * E", transition_type=TransitionType.ODE),
        Transition(origin="I", destination="H", equation="rho * I", transition_type=TransitionType.ODE),
        Transition(origin="I", destination="R", equation="gamma_i * I", transition_type=TransitionType.ODE),
        Transition(origin="H", destination="R", equation="gamma_h * H", transition_type=TransitionType.ODE),
        Transition(origin="I", destination="D", equation="mu_i * I", transition_type=TransitionType.ODE),
        Transition(origin="H", destination="D", equation="mu_h * H", transition_type=TransitionType.ODE),
        Transition(origin="D", destination="B", equation="psi * D", transition_type=TransitionType.ODE),
    ]

    state_names = ["S", "E", "I", "H", "R", "D", "B"]
    param_names = ["beta_i", "beta_h", "beta_d", "sigma", "rho", "gamma_i", "gamma_h", "mu_i", "mu_h", "psi"]

    model = SimulateOde(state_names=state_names, param_names=param_names, transitions=transitions)
    return model


def solve_model(
    params: Dict[str, float],
    initial_state: List[float],
    days: int = 180,
    dt: float = 0.1,
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve the Ebola model over time with PyGOM."""
    model = build_pygom_model()
    times = np.arange(0, days + dt, dt)
    solution = model.solve(init_state=initial_state, theta=params, time=times)
    return times, np.array(solution)


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

    if PYGOM_IMPORT_ERROR is not None:
        raise RuntimeError(f"PyGOM is not available in this environment: {PYGOM_IMPORT_ERROR}")

    print("PyGOM model object created successfully.")

    times, states = solve_model(params, initial_state, days=180, dt=0.1)

    labels = ["S", "E", "I", "H", "R", "D", "B"]
    final_state = dict(zip(labels, states[-1]))

    print("\nFinal state after 180 days:")
    for name, value in final_state.items():
        print(f"  {name}: {value:,.0f}")

    infectious_values = states[:, 2]
    peak_infectious = np.max(infectious_values)
    peak_day = times[np.argmax(infectious_values)]
    print(f"\nPeak infectious count: {peak_infectious:,.0f} on day {peak_day:.1f}")

    plot_results(times, states)


if __name__ == "__main__":
    main()
