"""Check the numerical tools against results that are known analytically."""

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pytest

from utils import (
    box_counting_dimension,
    classify_fixed_point_2d,
    continuation_1d,
    find_fixed_points_1d,
    integrate_system,
    jacobian_numerical,
    lyapunov_exponent_1d_map,
    lyapunov_exponents_flow,
    plot_1d_phase_line,
    plot_1d_potential,
    plot_1d_time_series,
    plot_3d_trajectory,
    plot_logistic_bifurcation,
    plot_lorenz_attractor,
    plot_phase_portrait,
    plot_streamlines,
    poincare_section,
)


def lorenz(u, sigma=10, rho=28, beta=8/3):
    return np.array([sigma * (u[1] - u[0]), u[0] * (rho - u[2]) - u[1], u[0] * u[1] - beta * u[2]])


def lorenz_jacobian(u, sigma=10, rho=28, beta=8/3):
    return np.array([[-sigma, sigma, 0], [rho - u[2], -1, -u[0]], [u[1], u[0], -beta]])


def test_integrate_system_matches_exponential_decay():
    for method in ('odeint', 'RK45'):
        _, sol = integrate_system(lambda s: -s, [1.0], (0, 1), 11, method=method)
        assert sol[-1, 0] == pytest.approx(np.exp(-1), abs=1e-3)


def test_find_fixed_points_1d():
    assert find_fixed_points_1d(np.sin, (-4, 4)) == pytest.approx([-np.pi, 0, np.pi], abs=1e-8)


@pytest.mark.parametrize('state', [[1.0, 2.0], [1, 2]], ids=['float', 'int'])
def test_jacobian_numerical(state):
    def f(s):
        return np.array([s[0] * s[1], s[0] + s[1]])

    assert jacobian_numerical(f, np.array(state)) == pytest.approx(np.array([[2, 1], [1, 1]]), abs=1e-5)


@pytest.mark.parametrize('eigenvalues, expected', [
    ([1j, -1j], 'Center'),
    ([-1 + 1j, -1 - 1j], 'Stable Focus'),
    ([1 + 1j, 1 - 1j], 'Unstable Focus'),
    ([1, -1], 'Saddle'),
    ([-1, -2], 'Stable Node'),
    ([1, 2], 'Unstable Node'),
    ([0, -1], 'Degenerate'),
])
def test_classify_fixed_point_2d(eigenvalues, expected):
    assert expected in classify_fixed_point_2d(np.array(eigenvalues))


def test_lyapunov_exponent_logistic_map():
    chaotic = lyapunov_exponent_1d_map(lambda x: 4 * x * (1 - x), lambda x: 4 - 8 * x, 0.3)
    assert chaotic == pytest.approx(np.log(2), abs=1e-3)

    # Stable fixed point x* = 1 - 1/r with multiplier 2 - r
    stable = lyapunov_exponent_1d_map(lambda x: 2.5 * x * (1 - x), lambda x: 2.5 - 5 * x, 0.3)
    assert stable == pytest.approx(np.log(0.5), abs=1e-3)


def test_lyapunov_spectrum_of_lorenz():
    spectrum = lyapunov_exponents_flow(lorenz, lorenz_jacobian, [1.0, 1.0, 1.0], t_total=300)

    assert spectrum[0] == pytest.approx(0.906, abs=0.1)
    assert spectrum[1] == pytest.approx(0.0, abs=0.1)
    assert spectrum[2] == pytest.approx(-14.57, abs=0.1)
    # The sum must equal the constant divergence -(sigma + 1 + beta)
    assert spectrum.sum() == pytest.approx(-(10 + 1 + 8/3), abs=0.02)


def test_poincare_section_lies_on_the_section():
    points = poincare_section(lorenz, np.array([1.0, 1.0, 1.0]), lambda u: u[2] - 27,
                              t_total=20, transient=5)
    assert len(points) > 0
    assert points[:, 2] == pytest.approx(27, abs=1e-3)


def test_continuation_follows_saddle_node_branch():
    r, x = continuation_1d(lambda x, r: r - x**2, 1.0, 1.0, 4.0, 4)
    assert x == pytest.approx(np.sqrt(r))


def test_box_counting_dimension_of_a_line():
    points = np.column_stack([np.linspace(0, 1, 100000)] * 2)
    dimension, _, _ = box_counting_dimension(points, epsilon_range=(0.005, 0.2))
    assert dimension == pytest.approx(1, abs=0.05)


def test_box_counting_dimension_of_a_filled_square():
    points = np.random.default_rng(0).random((200000, 2))
    dimension, _, _ = box_counting_dimension(points, epsilon_range=(0.02, 0.2))
    assert dimension == pytest.approx(2, abs=0.1)


@pytest.mark.parametrize('draw', [
    lambda: plot_1d_phase_line(lambda x: x**2 - 1, (-2, 2)),
    lambda: plot_1d_potential(lambda x: x**4 / 4 - x**2 / 2, (-2, 2)),
    lambda: plot_1d_time_series(lambda x: x - x**3, [0.1, -0.1], (0, 5)),
    lambda: plot_phase_portrait(lambda s: [s[1], -s[0]], (-2, 2), (-2, 2), [(1, 0)]),
    lambda: plot_streamlines(lambda s: [s[1], -s[0]], (-2, 2), (-2, 2)),
    lambda: plot_3d_trajectory(lorenz, [(1, 1, 1)], (0, 5), num_points=500),
    lambda: plot_lorenz_attractor(t_max=5, num_points=500),
    lambda: plot_logistic_bifurcation(),
], ids=['phase_line', 'potential', 'time_series', 'phase_portrait', 'streamlines',
        '3d_trajectory', 'lorenz', 'logistic'])
def test_plots_render(draw):
    draw()
    plt.close('all')


def test_star_import_does_not_leak_dependencies():
    namespace = {}
    exec('from utils import *', namespace)
    assert 'lyapunov_exponents_flow' in namespace
    assert 'np' not in namespace and 'plt' not in namespace
