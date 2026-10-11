"""Gas laws: ideal gas (PV = nRT) and van der Waals real gases."""
import numpy as np

R = 0.082057  # L*atm/(mol*K)


def ideal_gas(p=None, v=None, n=None, t=None):
    """Solve PV = nRT for the one argument left as None.
    Units: P (atm), V (L), n (mol), T (K). Works with scalars or NumPy arrays."""
    vals = {"p": p, "v": v, "n": n, "t": t}
    missing = [k for k, x in vals.items() if x is None]
    if len(missing) != 1:
        raise ValueError("Leave exactly one value empty")
    for k, x in vals.items():
        if x is not None and np.any(np.asarray(x) <= 0):
            raise ValueError(f"{k.upper()} must be > 0")
    k = missing[0]
    if k == "p": return n * R * t / v
    if k == "v": return n * R * t / p
    if k == "n": return p * v / (R * t)
    return p * v / (n * R)


# van der Waals constants: a (L^2 atm / mol^2), b (L / mol)
VDW_CONSTANTS = {
    "He": (0.0346, 0.0238), "H2": (0.2476, 0.02661), "N2": (1.370, 0.0387),
    "O2": (1.382, 0.03186), "Ar": (1.355, 0.03201), "CH4": (2.253, 0.04278),
    "CO2": (3.640, 0.04267), "NH3": (4.225, 0.0371), "H2O": (5.536, 0.03049),
}


def van_der_waals_pressure(v, n: float, t: float, a: float, b: float):
    """P = nRT/(V - nb) - a n^2 / V^2   (atm)."""
    v = np.asarray(v, dtype=float)
    return n * R * t / (v - n * b) - a * n ** 2 / v ** 2


def compressibility(p, v, n: float, t: float):
    """Z = PV / nRT (equals 1 for an ideal gas)."""
    return np.asarray(p) * np.asarray(v) / (n * R * t)
