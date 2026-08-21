## INTERFACE: Isambard vllm
## PROMPT: Provide Python code using the PyGOM package for this SIR compartmental model.
## Use these classes from PyGOM:
## DeterministicOde, Transition, TransitionType

import numpy as np
import matplotlib.pyplot as plt
from pygom import DeterministicOde, Transition, TransitionType


def create_sir_model():
    # Define transitions between compartments
    transitions = [
        # Infection: S -> I  (rate = beta * S * I)
        Transition(
            origin='S',
            equation='beta * S * I',
            transition_type=TransitionType.T,
            destination='I',
            name='infection'
        ),
        # Recovery: I -> R  (rate = gamma * I)
        Transition(
            origin='I',
            equation='gamma * I',
            transition_type=TransitionType.T,
            destination='R',
            name='recovery'
        ),
    ]

    # Build the model
    model = DeterministicOde(
        state=['S', 'I', 'R'],
        param=['beta', 'gamma'],
        transition=transitions
    )
    return model



def simulate_sir(
    N = 1000,
    beta=0.3,
    gamma=0.1,
    S0=999.0,
    I0=1.0,
    R0=0.0,
    t_max=50,
    t_points=500,
):
    """
    Simulate the SIR model and return the time series.

    Parameters
    ----------
    beta : float
        Transmission rate (default 0.3).
    gamma : float
        Recovery rate (default 0.1).
    S0, I0, R0 : float
        Initial compartment sizes (default: 999, 1, 0).
    t_max : float
        Final time for simulation (default 50).
    t_points : int
        Number of output time points (default 500).

    Returns
    -------
    time : ndarray
        Time points.
    solution : ndarray
        Array of shape (t_points, 3) with columns [S, I, R].
    model : DeterministicOde
        The PyGOM model object.
    """
    # Create model
    model = create_sir_model()

    # Set initial conditions and parameters
    model.initial_state = [S0, I0, R0]
    model.initial_time = 0.0
    model.parameters = [beta, gamma]

    # Print model summary
    print("=" * 50)
    print("SIR Model Summary")
    print("=" * 50)
    print(f"States   : {model.state_list}")
    print(f"Parameters: {model.param_list}")
    print(f"R0 estimate: {beta / gamma:.2f}")
    print()
    print("Transitions:")
    for t in model.transition_list:
        print(f"  {t}")
    print()
    print("ODE equations:")
    model.print_ode()
    print()

    # Integrate (model prepends t0 to the time array)
    time = np.linspace(0, t_max, t_points)
    solution = model.integrate(time)
    # solution includes the prepended t0, so slice it
    solution = solution[1:]
    # time is also extended by prepending t0, so slice
    time = np.linspace(0, t_max, t_points + 1)[1:]

    # Print final state
    print(f"Final values (t={t_max:.0f}):")
    print(f"  S = {solution[-1, 0]:.1f}")
    print(f"  I = {solution[-1, 1]:.1f}")
    print(f"  R = {solution[-1, 2]:.1f}")
    print("=" * 50)

    return time, solution, model


def plot_sir(time, solution):
    """Plot the SIR simulation results."""
    fig, ax = plt.subplots(figsize=(10, 6))

    S = solution[:, 0]
    I = solution[:, 1]
    R = solution[:, 2]

    ax.plot(time, S, label='Susceptible (S)', linewidth=2)
    ax.plot(time, I, label='Infected (I)', linewidth=2)
    ax.plot(time, R, label='Recovered (R)', linewidth=2)

    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Population', fontsize=12)
    ax.set_title('SIR Compartmental Model (PyGOM)', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('sir_results.png', dpi=150)
    try:
        plt.show()
    except Exception:
        pass  # Non-interactive backend (e.g., headless server)
    print("Plot saved to sir_results.png")

time, solution, model = simulate_sir()

simulate_sir()

plot_sir(time, solution)