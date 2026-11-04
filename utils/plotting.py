"""
Visualization utilities for dynamical systems.
Crafted for deep understanding, not just pretty pictures.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from typing import Callable, Tuple, Optional, List
from scipy.integrate import odeint


# ═══════════════════════════════════════════════════════════════════════════════
# STYLE CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

def set_style():
    """Set a clean, publication-quality style."""
    plt.rcParams.update({
        'figure.figsize': (10, 8),
        'figure.dpi': 100,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'lines.linewidth': 1.5,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'font.family': 'serif',
    })


# ═══════════════════════════════════════════════════════════════════════════════
# 1D SYSTEMS
# ═══════════════════════════════════════════════════════════════════════════════

def plot_1d_phase_line(f: Callable, x_range: Tuple[float, float],
                       num_points: int = 1000, ax=None,
                       show_flow: bool = True, show_fixed_points: bool = True):
    """
    Plot the phase line for a 1D system dx/dt = f(x).

    This is the fundamental tool for understanding 1D dynamics:
    - Where f(x) > 0: flow goes right
    - Where f(x) < 0: flow goes left
    - Where f(x) = 0: fixed points
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 5))

    x = np.linspace(x_range[0], x_range[1], num_points)
    y = np.array([f(xi) for xi in x])

    # Plot f(x)
    ax.plot(x, y, 'b-', linewidth=2, label=r"$f(x) = \dot{x}$")
    ax.axhline(y=0, color='k', linewidth=0.5)
    ax.axvline(x=0, color='k', linewidth=0.5)

    # Find and classify fixed points
    if show_fixed_points:
        # Find sign changes (zero crossings)
        sign_changes = np.where(np.diff(np.sign(y)))[0]
        for idx in sign_changes:
            # Linear interpolation to find zero
            x_star = x[idx] - y[idx] * (x[idx+1] - x[idx]) / (y[idx+1] - y[idx])

            # Compute derivative at fixed point (for stability)
            h = 1e-6
            df = (f(x_star + h) - f(x_star - h)) / (2 * h)

            if df < 0:
                ax.plot(x_star, 0, 'ko', markersize=12, label='Stabil' if idx == sign_changes[0] else '')
            else:
                ax.plot(x_star, 0, 'ko', markersize=12, fillstyle='none',
                       markeredgewidth=2, label='Instabil' if idx == sign_changes[0] else '')

    # Show flow direction on x-axis
    if show_flow:
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks([])

        # Draw arrows on a line below
        arrow_y = ax.get_ylim()[0] * 0.9
        for i in range(0, len(x)-1, len(x)//20):
            if abs(y[i]) > 0.01:
                direction = 1 if y[i] > 0 else -1
                ax.annotate('', xy=(x[i] + direction*0.1, arrow_y),
                           xytext=(x[i], arrow_y),
                           arrowprops=dict(arrowstyle='->', color='red', lw=1.5))

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$\dot{x} = f(x)$')
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


def plot_1d_potential(V: Callable, x_range: Tuple[float, float],
                      num_points: int = 1000, ax=None):
    """
    Plot the potential V(x) for a gradient system dx/dt = -dV/dx.

    The ball-rolling-downhill analogy:
    - System always moves downhill in V
    - Stable fixed points = local minima
    - Unstable fixed points = local maxima
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    x = np.linspace(x_range[0], x_range[1], num_points)
    v = np.array([V(xi) for xi in x])

    ax.plot(x, v, 'b-', linewidth=2)
    ax.fill_between(x, v, ax.get_ylim()[0], alpha=0.1)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$V(x)$')
    ax.set_title('Potential Landscape')
    ax.grid(True, alpha=0.3)

    return ax


def plot_1d_time_series(f: Callable, x0_list: List[float],
                        t_span: Tuple[float, float], ax=None):
    """Plot solutions x(t) for different initial conditions."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    t = np.linspace(t_span[0], t_span[1], 1000)

    colors = plt.cm.viridis(np.linspace(0, 1, len(x0_list)))

    for x0, color in zip(x0_list, colors):
        sol = odeint(lambda x, t: f(x), x0, t)
        ax.plot(t, sol[:, 0], color=color, label=f'$x_0 = {x0}$')

    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$x(t)$')
    ax.legend()
    ax.grid(True, alpha=0.3)

    return ax


# ═══════════════════════════════════════════════════════════════════════════════
# 2D SYSTEMS - VECTOR FIELDS AND PHASE PORTRAITS
# ═══════════════════════════════════════════════════════════════════════════════

def plot_vector_field(f: Callable, x_range: Tuple[float, float],
                      y_range: Tuple[float, float], density: int = 20,
                      ax=None, normalize: bool = True, cmap='coolwarm'):
    """
    Plot the vector field for a 2D system.

    f: function(state) -> [dx/dt, dy/dt]

    This is the foundation: at every point (x,y), there's a vector
    telling you where the system wants to go.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    x = np.linspace(x_range[0], x_range[1], density)
    y = np.linspace(y_range[0], y_range[1], density)
    X, Y = np.meshgrid(x, y)

    U = np.zeros_like(X)
    V = np.zeros_like(Y)

    for i in range(density):
        for j in range(density):
            derivatives = f([X[i,j], Y[i,j]])
            U[i,j] = derivatives[0]
            V[i,j] = derivatives[1]

    # Magnitude for coloring
    M = np.sqrt(U**2 + V**2)

    if normalize:
        # Normalize arrows to same length (direction only)
        U_norm = U / (M + 1e-10)
        V_norm = V / (M + 1e-10)
        ax.quiver(X, Y, U_norm, V_norm, M, cmap=cmap, alpha=0.8)
    else:
        ax.quiver(X, Y, U, V, M, cmap=cmap, alpha=0.8)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    return ax


def plot_streamlines(f: Callable, x_range: Tuple[float, float],
                     y_range: Tuple[float, float], density: int = 50,
                     ax=None, color_by_speed: bool = True):
    """
    Plot streamlines (integral curves of the vector field).

    Streamlines are tangent to the vector field everywhere -
    they show the actual paths that solutions follow.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    x = np.linspace(x_range[0], x_range[1], density)
    y = np.linspace(y_range[0], y_range[1], density)
    X, Y = np.meshgrid(x, y)

    U = np.zeros_like(X)
    V = np.zeros_like(Y)

    for i in range(density):
        for j in range(density):
            derivatives = f([X[i,j], Y[i,j]])
            U[i,j] = derivatives[0]
            V[i,j] = derivatives[1]

    speed = np.sqrt(U**2 + V**2)

    if color_by_speed:
        strm = ax.streamplot(X, Y, U, V, color=speed, cmap='coolwarm',
                            density=1.5, linewidth=1, arrowsize=1.2)
        plt.colorbar(strm.lines, ax=ax, label='Speed')
    else:
        ax.streamplot(X, Y, U, V, density=1.5, linewidth=1, arrowsize=1.2)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.set_aspect('equal')

    return ax


def plot_nullclines(f: Callable, x_range: Tuple[float, float],
                    y_range: Tuple[float, float], density: int = 500, ax=None):
    """
    Plot the nullclines of a 2D system.

    x-nullcline: where dx/dt = 0 (vertical arrows only)
    y-nullcline: where dy/dt = 0 (horizontal arrows only)
    Intersections: fixed points!

    This is THE tool for qualitative analysis.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    x = np.linspace(x_range[0], x_range[1], density)
    y = np.linspace(y_range[0], y_range[1], density)
    X, Y = np.meshgrid(x, y)

    dX = np.zeros_like(X)
    dY = np.zeros_like(Y)

    for i in range(density):
        for j in range(density):
            derivatives = f([X[i,j], Y[i,j]])
            dX[i,j] = derivatives[0]
            dY[i,j] = derivatives[1]

    # Plot nullclines as contours at level 0
    ax.contour(X, Y, dX, levels=[0], colors='blue', linewidths=2, linestyles='-')
    ax.contour(X, Y, dY, levels=[0], colors='red', linewidths=2, linestyles='-')

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)

    # Add legend manually
    ax.plot([], [], 'b-', linewidth=2, label=r'$\dot{x} = 0$ (x-nullcline)')
    ax.plot([], [], 'r-', linewidth=2, label=r'$\dot{y} = 0$ (y-nullcline)')
    ax.legend()

    return ax


def plot_trajectory(f: Callable, initial_conditions: List[Tuple[float, float]],
                    t_span: Tuple[float, float], ax=None, num_points: int = 2000,
                    show_start: bool = True, show_end: bool = True):
    """
    Plot trajectories from given initial conditions.

    This is where you see the actual motion of the system.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))

    t = np.linspace(t_span[0], t_span[1], num_points)
    colors = plt.cm.viridis(np.linspace(0, 1, len(initial_conditions)))

    for (x0, y0), color in zip(initial_conditions, colors):
        sol = odeint(lambda state, t: f(state), [x0, y0], t)
        ax.plot(sol[:, 0], sol[:, 1], color=color, linewidth=1.5)

        if show_start:
            ax.plot(x0, y0, 'o', color=color, markersize=8)
        if show_end:
            ax.plot(sol[-1, 0], sol[-1, 1], 's', color=color, markersize=6)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    return ax


def plot_phase_portrait(f: Callable, x_range: Tuple[float, float],
                        y_range: Tuple[float, float],
                        initial_conditions: Optional[List[Tuple[float, float]]] = None,
                        t_span: Tuple[float, float] = (0, 20),
                        show_nullclines: bool = True,
                        show_vector_field: bool = True,
                        figsize: Tuple[int, int] = (12, 10)):
    """
    Complete phase portrait: vector field + nullclines + trajectories.

    The complete picture of the system's dynamics.
    """
    fig, ax = plt.subplots(figsize=figsize)

    if show_vector_field:
        plot_vector_field(f, x_range, y_range, density=20, ax=ax)

    if show_nullclines:
        plot_nullclines(f, x_range, y_range, ax=ax)

    if initial_conditions is not None:
        plot_trajectory(f, initial_conditions, t_span, ax=ax)

    ax.set_title('Phase Portrait')

    return fig, ax


# ═══════════════════════════════════════════════════════════════════════════════
# 3D SYSTEMS AND CHAOS
# ═══════════════════════════════════════════════════════════════════════════════

def plot_3d_trajectory(f: Callable, initial_conditions: List[Tuple[float, float, float]],
                       t_span: Tuple[float, float], ax=None, num_points: int = 10000,
                       alpha: float = 0.8):
    """
    Plot 3D trajectories - essential for visualizing chaotic attractors.
    """
    if ax is None:
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

    t = np.linspace(t_span[0], t_span[1], num_points)
    colors = plt.cm.plasma(np.linspace(0, 1, len(initial_conditions)))

    for ic, color in zip(initial_conditions, colors):
        sol = odeint(lambda state, t: f(state), ic, t)
        ax.plot(sol[:, 0], sol[:, 1], sol[:, 2],
               color=color, linewidth=0.5, alpha=alpha)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$z$')

    return ax


def plot_lorenz_attractor(sigma: float = 10, rho: float = 28, beta: float = 8/3,
                          initial: Tuple[float, float, float] = (1, 1, 1),
                          t_max: float = 50, num_points: int = 10000):
    """
    The Lorenz attractor - the icon of chaos theory.

    Named after Edward Lorenz, who discovered it while studying
    atmospheric convection. The butterfly shape is no coincidence -
    it's a strange attractor with fractal dimension ~2.06.
    """
    def lorenz(state, t):
        x, y, z = state
        return [sigma * (y - x),
                x * (rho - z) - y,
                x * y - beta * z]

    t = np.linspace(0, t_max, num_points)
    sol = odeint(lorenz, initial, t)

    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Color by time for visual effect
    points = sol.reshape(-1, 1, 3)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    ax.plot(sol[:, 0], sol[:, 1], sol[:, 2], linewidth=0.5, alpha=0.8)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$z$')
    ax.set_title(f'Lorenz Attractor ($\\sigma={sigma}, \\rho={rho}, \\beta={beta:.2f}$)')

    return fig, ax, sol


# ═══════════════════════════════════════════════════════════════════════════════
# BIFURCATION DIAGRAMS
# ═══════════════════════════════════════════════════════════════════════════════

def plot_bifurcation_diagram_1d(f: Callable, r_range: Tuple[float, float],
                                x_range: Tuple[float, float],
                                num_r: int = 500, num_x: int = 200,
                                figsize: Tuple[int, int] = (12, 8)):
    """
    Bifurcation diagram for a 1D map x_{n+1} = f(x_n, r).

    Shows how the long-term behavior changes with parameter r.
    The route to chaos becomes visible.
    """
    fig, ax = plt.subplots(figsize=figsize)

    r_values = np.linspace(r_range[0], r_range[1], num_r)

    for r in r_values:
        x = 0.5  # Initial condition

        # Transient
        for _ in range(500):
            x = f(x, r)

        # Collect attractors
        x_values = []
        for _ in range(num_x):
            x = f(x, r)
            x_values.append(x)

        ax.plot([r] * len(x_values), x_values, 'k,', markersize=0.5)

    ax.set_xlabel(r'Parameter $r$')
    ax.set_ylabel(r'$x^*$ (Attractor)')
    ax.set_title('Bifurcation Diagram')
    ax.set_xlim(r_range)
    ax.set_ylim(x_range)

    return fig, ax


def plot_logistic_bifurcation(r_min: float = 2.5, r_max: float = 4.0,
                              figsize: Tuple[int, int] = (14, 8)):
    """
    The logistic map bifurcation diagram - perhaps the most famous
    image in chaos theory.

    x_{n+1} = r * x_n * (1 - x_n)

    Watch order turn to chaos through period-doubling.
    """
    def logistic(x, r):
        return r * x * (1 - x)

    return plot_bifurcation_diagram_1d(logistic, (r_min, r_max), (0, 1),
                                        num_r=1000, num_x=300, figsize=figsize)
