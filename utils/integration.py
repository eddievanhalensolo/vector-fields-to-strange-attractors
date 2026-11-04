"""
Numerical integration and analysis tools for dynamical systems.
"""

import numpy as np
from scipy.integrate import odeint, solve_ivp
from typing import Callable, Tuple, List

__all__ = [
    'integrate_system', 'find_fixed_points_1d', 'jacobian_numerical',
    'classify_fixed_point_2d', 'lyapunov_exponent_1d_map',
    'lyapunov_exponents_flow', 'poincare_section', 'continuation_1d',
    'box_counting_dimension',
]


# ═══════════════════════════════════════════════════════════════════════════════
# NUMERICAL INTEGRATION
# ═══════════════════════════════════════════════════════════════════════════════

def integrate_system(f: Callable, initial_state: np.ndarray,
                     t_span: Tuple[float, float], num_points: int = 1000,
                     method: str = 'odeint') -> Tuple[np.ndarray, np.ndarray]:
    """
    Integrate a system of ODEs.

    Parameters:
    -----------
    f : callable
        Function that returns derivatives: f(state) -> d(state)/dt
    initial_state : array-like
        Initial conditions
    t_span : tuple
        (t_start, t_end)
    num_points : int
        Number of time points
    method : str
        'odeint' (LSODA) or 'RK45'

    Returns:
    --------
    t : ndarray
        Time points
    sol : ndarray
        Solution array (num_points x num_variables)
    """
    t = np.linspace(t_span[0], t_span[1], num_points)

    if method == 'odeint':
        sol = odeint(lambda state, t: f(state), initial_state, t)
    else:
        result = solve_ivp(lambda t, state: f(state), t_span, initial_state,
                          method=method, t_eval=t)
        sol = result.y.T

    return t, sol


# ═══════════════════════════════════════════════════════════════════════════════
# FIXED POINT ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def find_fixed_points_1d(f: Callable, x_range: Tuple[float, float],
                         num_initial: int = 100, tol: float = 1e-10) -> List[float]:
    """
    Find fixed points of a 1D system dx/dt = f(x) using Newton's method.
    """
    from scipy.optimize import fsolve

    fixed_points = []
    x_initial = np.linspace(x_range[0], x_range[1], num_initial)

    for x0 in x_initial:
        try:
            x_star = fsolve(f, x0, full_output=True)
            if x_star[2] == 1:  # Solution found
                x_fp = x_star[0][0]
                if x_range[0] <= x_fp <= x_range[1]:
                    # Check if it's actually a fixed point
                    if abs(f(x_fp)) < tol:
                        # Check if we already found this one
                        is_new = True
                        for existing in fixed_points:
                            if abs(x_fp - existing) < 1e-6:
                                is_new = False
                                break
                        if is_new:
                            fixed_points.append(x_fp)
        except (ValueError, RuntimeError):
            pass

    return sorted(fixed_points)


def jacobian_numerical(f: Callable, state: np.ndarray, h: float = 1e-8) -> np.ndarray:
    """
    Compute the Jacobian matrix numerically using central differences.

    The Jacobian tells you everything about local stability:
    - Its eigenvalues determine whether perturbations grow or decay
    - Its eigenvectors show the directions of growth/decay
    """
    state = np.array(state, dtype=float)
    n = len(state)
    J = np.zeros((n, n))

    for i in range(n):
        state_plus = state.copy()
        state_minus = state.copy()
        state_plus[i] += h
        state_minus[i] -= h

        f_plus = np.array(f(state_plus))
        f_minus = np.array(f(state_minus))

        J[:, i] = (f_plus - f_minus) / (2 * h)

    return J


def classify_fixed_point_2d(eigenvalues: np.ndarray) -> str:
    """
    Classify a 2D fixed point based on eigenvalues.

    This is the taxonomy of equilibria:
    - Nodes: real eigenvalues, same sign
    - Saddles: real eigenvalues, opposite signs
    - Spirals/Foci: complex eigenvalues
    - Centers: purely imaginary eigenvalues
    """
    ev = eigenvalues

    # Check if complex
    if np.any(np.iscomplex(ev)):
        real_part = np.real(ev[0])
        if abs(real_part) < 1e-10:
            return "Center (marginally stable)"
        elif real_part < 0:
            return "Stable Focus (spiral sink)"
        else:
            return "Unstable Focus (spiral source)"

    # Real eigenvalues
    ev = np.real(ev)

    if ev[0] * ev[1] < 0:
        return "Saddle Point"
    elif ev[0] < 0 and ev[1] < 0:
        if abs(ev[0] - ev[1]) < 1e-10:
            return "Stable Star Node"
        else:
            return "Stable Node"
    elif ev[0] > 0 and ev[1] > 0:
        if abs(ev[0] - ev[1]) < 1e-10:
            return "Unstable Star Node"
        else:
            return "Unstable Node"
    else:
        return "Degenerate (zero eigenvalue)"


# ═══════════════════════════════════════════════════════════════════════════════
# LYAPUNOV EXPONENTS - THE MEASURE OF CHAOS
# ═══════════════════════════════════════════════════════════════════════════════

def lyapunov_exponent_1d_map(f: Callable, df: Callable, x0: float,
                              num_iterations: int = 10000,
                              transient: int = 1000) -> float:
    """
    Compute the Lyapunov exponent for a 1D map x_{n+1} = f(x_n).

    λ = lim_{n→∞} (1/n) Σ log|f'(x_i)|

    If λ > 0: Chaos (nearby trajectories diverge exponentially)
    If λ < 0: Stable (nearby trajectories converge)
    If λ = 0: Marginal (periodic orbit or bifurcation point)
    """
    x = x0

    # Transient
    for _ in range(transient):
        x = f(x)

    # Compute Lyapunov exponent
    lyap_sum = 0.0
    for _ in range(num_iterations):
        derivative = abs(df(x))
        if derivative > 0:
            lyap_sum += np.log(derivative)
        x = f(x)

    return lyap_sum / num_iterations


def lyapunov_exponents_flow(f: Callable, jacobian: Callable,
                            initial_state: np.ndarray,
                            t_total: float = 100,
                            dt: float = 0.01) -> np.ndarray:
    """
    Compute all Lyapunov exponents for a continuous flow.

    Uses the standard algorithm of integrating tangent vectors
    and periodically orthonormalizing via Gram-Schmidt.

    This reveals the full spectrum of chaos:
    - Positive exponents: chaotic directions
    - Zero exponent: along the flow (always present in flows)
    - Negative exponents: contracting directions

    Sum of exponents = average divergence (phase space volume change)
    """
    n = len(initial_state)
    num_steps = int(t_total / dt)
    reorthonormalize_every = 10

    # Initialize state and tangent vectors
    state = np.array(initial_state, dtype=float)
    Q = np.eye(n)  # Orthonormal tangent vectors

    lyap_sums = np.zeros(n)

    def rhs(s, q):
        return np.array(f(s)), jacobian(s) @ q

    for step in range(1, num_steps + 1):
        # Integrate state and tangent vectors together with RK4
        k1_s, k1_q = rhs(state, Q)
        k2_s, k2_q = rhs(state + 0.5 * dt * k1_s, Q + 0.5 * dt * k1_q)
        k3_s, k3_q = rhs(state + 0.5 * dt * k2_s, Q + 0.5 * dt * k2_q)
        k4_s, k4_q = rhs(state + dt * k3_s, Q + dt * k3_q)
        state = state + (dt / 6) * (k1_s + 2*k2_s + 2*k3_s + k4_s)
        Q = Q + (dt / 6) * (k1_q + 2*k2_q + 2*k3_q + k4_q)

        # Periodically reorthonormalize
        if step % reorthonormalize_every == 0 or step == num_steps:
            Q, R = np.linalg.qr(Q)
            lyap_sums += np.log(np.abs(np.diag(R)))

    t = num_steps * dt
    return lyap_sums / t


# ═══════════════════════════════════════════════════════════════════════════════
# POINCARÉ SECTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def poincare_section(f: Callable, initial_state: np.ndarray,
                     section_condition: Callable,
                     section_direction: int = 1,
                     t_total: float = 1000,
                     dt: float = 0.01,
                     transient: float = 100) -> np.ndarray:
    """
    Compute a Poincaré section.

    The Poincaré section reduces a continuous flow to a discrete map
    by recording where trajectories cross a chosen surface.

    This is how you see periodic orbits as fixed points
    and chaos as scattered points.

    Parameters:
    -----------
    f : callable
        System dynamics
    initial_state : array
        Starting point
    section_condition : callable
        Function that is zero on the section surface
    section_direction : int
        +1 for crossing from negative to positive, -1 for opposite
    """
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    crossings = []
    state = np.array(initial_state)
    t = 0

    # Skip transient
    result = solve_ivp(lambda t, s: f(s), (0, transient), state,
                      method='RK45', dense_output=True)
    state = result.y[:, -1]

    prev_condition = section_condition(state)

    # Find crossings
    while t < t_total:
        result = solve_ivp(lambda t, s: f(s), (0, dt), state,
                          method='RK45', dense_output=True)
        new_state = result.y[:, -1]
        new_condition = section_condition(new_state)

        # Check for crossing
        if prev_condition * new_condition < 0:
            if (new_condition - prev_condition) * section_direction > 0:
                # Interpolate to find exact crossing
                def crossing_func(tau):
                    s = result.sol(tau)
                    return section_condition(s)

                try:
                    tau_cross = brentq(crossing_func, 0, dt)
                    crossing_state = result.sol(tau_cross)
                    crossings.append(crossing_state)
                except ValueError:
                    pass

        state = new_state
        prev_condition = new_condition
        t += dt

    return np.array(crossings)


# ═══════════════════════════════════════════════════════════════════════════════
# BIFURCATION ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def continuation_1d(f: Callable, x_start: float, r_start: float,
                    r_end: float, num_steps: int = 1000) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simple parameter continuation for 1D systems.

    Follows a branch of fixed points as parameter changes.
    """
    from scipy.optimize import fsolve

    r_values = np.linspace(r_start, r_end, num_steps)
    x_values = np.zeros(num_steps)

    x = x_start
    for i, r in enumerate(r_values):
        try:
            x = fsolve(lambda x: f(x, r), x)[0]
            x_values[i] = x
        except (ValueError, RuntimeError):
            x_values[i] = np.nan

    return r_values, x_values


# ═══════════════════════════════════════════════════════════════════════════════
# FRACTAL DIMENSION
# ═══════════════════════════════════════════════════════════════════════════════

def box_counting_dimension(points: np.ndarray,
                           epsilon_range: Tuple[float, float] = None,
                           num_scales: int = 20) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Estimate the box-counting (Minkowski) dimension of a point set.

    D = lim_{ε→0} log(N(ε)) / log(1/ε)

    where N(ε) is the number of boxes of size ε needed to cover the set.

    The fit is only meaningful inside the scaling range, so choose epsilon_range
    for your data. With too few points N(ε) saturates at small ε and the estimate
    comes out too low: a Lorenz trajectory gives about 1.9 even with 2 million
    points, against a literature value of about 2.06.
    """
    # Normalize points to unit hypercube
    points = np.array(points)
    min_vals = points.min(axis=0)
    max_vals = points.max(axis=0)
    scale = max_vals - min_vals
    scale[scale == 0] = 1
    points_normalized = (points - min_vals) / scale

    if epsilon_range is None:
        epsilon_range = (0.001, 0.5)

    epsilons = np.logspace(np.log10(epsilon_range[0]),
                           np.log10(epsilon_range[1]), num_scales)
    counts = []

    for eps in epsilons:
        # Discretize points to grid
        grid_indices = np.floor(points_normalized / eps).astype(int)
        # Count unique boxes
        unique_boxes = len(set(map(tuple, grid_indices)))
        counts.append(unique_boxes)

    counts = np.array(counts)

    # Fit line to log-log plot
    log_eps = np.log(1 / epsilons)
    log_N = np.log(counts)

    # Linear regression
    coeffs = np.polyfit(log_eps, log_N, 1)
    dimension = coeffs[0]

    return dimension, epsilons, counts
