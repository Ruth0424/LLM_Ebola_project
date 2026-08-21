
## PROMPT: Provide Python code using the PyGOM package for this SIR compartmental model
## INTERFACE: Isambard vllm
# Makes use of a class ‘Event’, which is not a PyGOM class, instead of including an ‘equation =’ for defining the movement between compartments.

"""SIR compartmental model built with PyGOM."""

import numpy as np
import matplotlib.pyplot as plt
from pygom.model import SimulateOde, Transition, TransitionType, Event


def create_sir_model():
    """Build a standard SIR model using PyGOM's transition-based API.

    dS/dt = -β·S·I / N
    dI/dt =  β·S·I / N - γ·I
    dR/dt =  γ·I
    """
    # --- Define compartments (state, (lower_bound, upper_bound)) ---
    state_list = [('S', (0, None)), ('I', (0, None)), ('R', (0, None))]

    # --- Define parameters ---
    param_list = ['beta', 'gamma', 'N']

    # --- Define transitions (events) ---
    # 1) Infection: S → I at rate β·S·I/N
    trans_infection = Transition(
        origin='S',
        destination='I',
        transition_type=TransitionType.T
    )
    event_infection = Event(
        rate='beta * S * I / N',
        transition_list=[trans_infection]
    )

    # 2) Recovery: I → R at rate γ·I
    trans_recovery = Transition(
        origin='I',
        destination='R',
        transition_type=TransitionType.T
    )
    event_recovery = Event(
        rate='gamma * I',
        transition_list=[trans_recovery]
    )

    # --- Assemble the model ---
    model = SimulateOde(
        state=state_list,
        param=param_list,
        event=[event_infection, event_recovery]
    )

    return model


# ------------------------------------------------------------------ #
#  Convenience: use PyGOM's built-in SIR shortcut
# ------------------------------------------------------------------ #
def create_sir_model_builtin():
    """Build SIR using pygom.common_models.SIR()."""
    from pygom.model.common_models import SIR
    return SIR()


# ------------------------------------------------------------------ #
#  Simulation & plotting
# ------------------------------------------------------------------ #
def simulate(model, t_max=100, steps=1000):
    """Run deterministic simulation and return (t, S, I, R)."""
    t = np.linspace(0, t_max, steps)
    model.initial_values = ([990_000.0, 1_000.0, 0.0], t[0])
    solution = model.integrate(t[1:])          # skip t=0 (duplicate)
    return t, solution[:, 0], solution[:, 1], solution[:, 2]


def plot_sir(t, S, I, R, beta, gamma):
    """Plot trajectories and display R₀."""
    R0 = beta / gamma

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t, S, label='Susceptible', linewidth=2)
    ax.plot(t, I, label='Infected', linewidth=2)
    ax.plot(t, R, label='Recovered', linewidth=2)
    ax.set_xlabel('Time')
    ax.set_ylabel('Population')
    ax.set_title(f'SIR Model (PyGOM)\nβ={beta}, γ={gamma}, R₀={R0:.2f}')
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


# ------------------------------------------------------------------ #
#  Demo
# ------------------------------------------------------------------ #
if __name__ == '__main__':
    # Option A: build from transitions (full PyGOM style)
    model = create_sir_model()
    # Option B: use built-in shortcut
    # model = create_sir_model_builtin()

    # Set parameter values
    model.parameters = {
        'beta': 0.3,
        'gamma': 0.1,
        'N': 1_000_000,
    }

    t, S, I, R = simulate(model)
    plot_sir(t, S, I, R, beta=0.3, gamma=0.1)