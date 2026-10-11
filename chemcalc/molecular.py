"""Molecular-scale models: Lennard-Jones potential, Maxwell-Boltzmann speeds."""
import numpy as np

R_GAS = 8.314462618  # J/(mol*K)
KB_KJ = 0.008314462618  # kJ/(mol*K)


def lennard_jones(r, epsilon: float, sigma: float):
    """V(r) = 4*eps*[(sigma/r)^12 - (sigma/r)^6]"""
    sr6 = (sigma / np.asarray(r, dtype=float)) ** 6
    return 4 * epsilon * (sr6 ** 2 - sr6)


def lj_force(r, epsilon: float, sigma: float):
    """F(r) = -dV/dr"""
    r = np.asarray(r, dtype=float)
    return 24 * epsilon / r * (2 * (sigma / r) ** 12 - (sigma / r) ** 6)


def lj_minimum(sigma: float) -> float:
    """Equilibrium separation r_min = 2^(1/6) * sigma (where V = -epsilon)."""
    return 2 ** (1 / 6) * sigma


def mb_speeds(molar_mass_g: float, temp_k: float) -> dict:
    """Most probable, mean and RMS speeds (m/s)."""
    m = molar_mass_g / 1000.0
    return {
        "most_probable": float(np.sqrt(2 * R_GAS * temp_k / m)),
        "mean": float(np.sqrt(8 * R_GAS * temp_k / (np.pi * m))),
        "rms": float(np.sqrt(3 * R_GAS * temp_k / m)),
    }


def mb_distribution(v, molar_mass_g: float, temp_k: float):
    """Maxwell-Boltzmann probability density f(v) in s/m."""
    m = molar_mass_g / 1000.0
    v = np.asarray(v, dtype=float)
    a = m / (2 * np.pi * R_GAS * temp_k)
    return 4 * np.pi * a ** 1.5 * v ** 2 * np.exp(-m * v ** 2 / (2 * R_GAS * temp_k))
