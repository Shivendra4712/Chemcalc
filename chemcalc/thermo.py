"""Thermodynamics (Gibbs energy, K, van 't Hoff) and electrochemistry (Nernst)."""
import math

R = 8.314462618       # J/(mol K)
F = 96485.33212       # C/mol


def gibbs(dh_kj: float, ds_j: float, temp_k: float) -> float:
    """dG = dH - T dS  (kJ/mol)."""
    return dh_kj - temp_k * ds_j / 1000.0


def crossover_temperature(dh_kj: float, ds_j: float):
    """Temperature where dG = 0 (None if it does not exist)."""
    if ds_j == 0:
        return None
    t = dh_kj * 1000.0 / ds_j
    return t if t > 0 else None


def k_from_dg(dg_kj: float, temp_k: float) -> float:
    x = -dg_kj * 1000.0 / (R * temp_k)
    return math.exp(x) if x < 700 else math.inf


def dg_from_k(k: float, temp_k: float) -> float:
    return -R * temp_k * math.log(k) / 1000.0


def van_t_hoff(k1: float, t1: float, dh_kj: float, t2: float) -> float:
    """K at T2 from K at T1: ln(K2/K1) = -dH/R (1/T2 - 1/T1)."""
    return k1 * math.exp(-dh_kj * 1000.0 / R * (1.0 / t2 - 1.0 / t1))


def nernst(e0: float, n: int, q: float, temp_k: float = 298.15) -> float:
    """E = E0 - (RT/nF) ln Q"""
    if q <= 0:
        raise ValueError("Q must be > 0")
    return e0 - (R * temp_k / (n * F)) * math.log(q)


# Standard reduction potentials (V, 25 C): name -> (E0, electrons)
POTENTIALS = {
    "Li+/Li": (-3.04, 1), "Ca2+/Ca": (-2.87, 2), "Na+/Na": (-2.71, 1),
    "Mg2+/Mg": (-2.37, 2), "Al3+/Al": (-1.66, 3), "Zn2+/Zn": (-0.76, 2),
    "Fe2+/Fe": (-0.44, 2), "Ni2+/Ni": (-0.25, 2), "Pb2+/Pb": (-0.13, 2),
    "H+/H2": (0.00, 2), "Cu2+/Cu": (0.34, 2), "I2/I-": (0.54, 2),
    "Fe3+/Fe2+": (0.77, 1), "Ag+/Ag": (0.80, 1), "Br2/Br-": (1.07, 2),
    "Cl2/Cl-": (1.36, 2), "Au3+/Au": (1.50, 3),
}


def galvanic_cell(cathode: str, anode: str) -> dict:
    """Standard cell potential, dG and K for a galvanic cell."""
    ec, nc = POTENTIALS[cathode]
    ea, na = POTENTIALS[anode]
    n = nc * na // math.gcd(nc, na)
    e0 = ec - ea
    dg = -n * F * e0 / 1000.0
    return {"E0": e0, "n": n, "dG_kJ": dg, "K": k_from_dg(dg, 298.15),
            "spontaneous": e0 > 0}
