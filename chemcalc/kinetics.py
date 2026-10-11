"""Chemical kinetics: rate laws, half-life, Arrhenius, consecutive reactions."""
import numpy as np

R = 8.314462618  # J/(mol*K)


def concentration(order: int, a0: float, k: float, t):
    """[A](t) for a 0th, 1st or 2nd order reaction."""
    t = np.asarray(t, dtype=float)
    if order == 0:
        return np.maximum(a0 - k * t, 0.0)
    if order == 1:
        return a0 * np.exp(-k * t)
    if order == 2:
        return 1.0 / (1.0 / a0 + k * t)
    raise ValueError("order must be 0, 1 or 2")


def half_life(order: int, a0: float, k: float) -> float:
    if order == 0:
        return a0 / (2 * k)
    if order == 1:
        return float(np.log(2) / k)
    if order == 2:
        return 1.0 / (k * a0)
    raise ValueError("order must be 0, 1 or 2")


def arrhenius(a_factor: float, ea_j_mol: float, temp_k: float) -> float:
    """k = A * exp(-Ea / RT)"""
    return float(a_factor * np.exp(-ea_j_mol / (R * temp_k)))


def activation_energy(k1: float, t1: float, k2: float, t2: float) -> float:
    """Ea (J/mol) from rate constants at two temperatures (K)."""
    return R * np.log(k2 / k1) / (1.0 / t1 - 1.0 / t2)


def consecutive(a0: float, k1: float, k2: float, t):
    """A -> B -> C (both first order). Returns (A, B, C) arrays."""
    t = np.asarray(t, dtype=float)
    a = a0 * np.exp(-k1 * t)
    if np.isclose(k1, k2):
        b = a0 * k1 * t * np.exp(-k1 * t)
    else:
        b = a0 * k1 / (k2 - k1) * (np.exp(-k1 * t) - np.exp(-k2 * t))
    return a, b, a0 - a - b
