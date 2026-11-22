"""Analytical open-channel inflow profiles (rigid lid, smooth bed).

Used for the inlet ``boundaryData`` and for the initial fields of the
channel tests. Heights are measured from the local bed elevation z_b;
values below the bed are zero (velocity) or small floors (k, omega).

Laminar: half-Poiseuille profile of an open channel with a slip lid,
    u(eta) = 1.5 U_b (2 eta - eta^2),  eta = (z - z_b)/h,
    bed shear stress tau_b/rho = 3 nu U_b / h.
Turbulent: 1/7 power law for u with depth-mean U_b, a linear k profile
    k = u_*^2/sqrt(C_mu) (1 - eta) and omega = sqrt(k)/(C_mu^0.25 l),
    with a Blackadar-type mixing length l = kappa y/(1 + kappa y/(0.09 h)).
    u_* follows from the depth-integrated smooth-wall log law.
"""

import math
from typing import Tuple

import numpy as np

KAPPA = 0.41      # von Karman constant [-]
B_SMOOTH = 5.2    # smooth-wall log-law intercept [-]
C_MU = 0.09       # k-omega model constant [-]


def friction_velocity(u_bulk: float, depth: float, nu: float) -> float:
    """Return u_* [m/s] of a smooth open channel from the depth-mean log law.

    Solves U_b/u_* = ln(h u_*/nu)/kappa + B - 1/kappa by fixed-point
    iteration.

    Parameters
    ----------
    u_bulk : float
        Depth-averaged velocity U_b [m/s].
    depth : float
        Flow depth h [m].
    nu : float
        Kinematic viscosity [m^2/s].
    """
    u_star = 0.05 * u_bulk
    for _ in range(100):
        ratio = math.log(depth * u_star / nu) / KAPPA + B_SMOOTH - 1.0 / KAPPA
        u_new = u_bulk / ratio
        if abs(u_new - u_star) < 1e-12:
            break
        u_star = u_new
    return u_star


def laminar_profile(y: np.ndarray, depth: np.ndarray,
                    u_bulk: float) -> np.ndarray:
    """Return u [m/s] of the laminar open-channel profile.

    Parameters
    ----------
    y : array
        Height above the bed [m] (negative below the bed).
    depth : array
        Local flow depth h [m].
    u_bulk : float
        Depth-averaged velocity [m/s].
    """
    eta = np.clip(y / depth, 0.0, 1.0)
    u = 1.5 * u_bulk * (2.0 * eta - eta ** 2)
    return np.where(y > 0.0, u, 0.0)


def turbulent_profile(
    y: np.ndarray, depth: np.ndarray, u_bulk: float, nu: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (u [m/s], k [m^2/s^2], omega [1/s]) of the turbulent profile.

    Parameters
    ----------
    y : array
        Height above the bed [m] (negative below the bed).
    depth : array
        Local flow depth h [m].
    u_bulk : float
        Depth-averaged velocity [m/s].
    nu : float
        Kinematic viscosity [m^2/s].
    """
    depth = np.asarray(depth, dtype=float) * np.ones_like(y)
    u_star = np.array([friction_velocity(u_bulk, h, nu) for h in
                       np.atleast_1d(depth).ravel()]).reshape(depth.shape)
    eta = np.clip(y / depth, 1e-6, 1.0)
    u = u_bulk * 8.0 / 7.0 * eta ** (1.0 / 7.0)
    k_floor = 1e-3 * u_star ** 2
    k = u_star ** 2 / math.sqrt(C_MU) * (1.0 - eta) + k_floor
    y_pos = np.maximum(y, 1e-4 * depth)
    mix = KAPPA * y_pos / (1.0 + KAPPA * y_pos / (0.09 * depth))
    omega = np.sqrt(k) / (C_MU ** 0.25 * mix)
    below = y <= 0.0
    u = np.where(below, 0.0, u)
    k = np.where(below, k_floor, k)
    omega = np.where(below, np.sqrt(k_floor) / (C_MU ** 0.25 * KAPPA
                                                * 1e-4 * depth), omega)
    return u, k, omega


def log_law_u_plus(y_plus: np.ndarray) -> np.ndarray:
    """Return u+ of the smooth-wall log law (linear for y+ < 11.25)."""
    y_plus = np.asarray(y_plus, dtype=float)
    log_part = np.log(np.maximum(y_plus, 1e-12)) / KAPPA + B_SMOOTH
    return np.where(y_plus < 11.25, y_plus, log_part)
