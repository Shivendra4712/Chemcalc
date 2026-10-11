"""Solution chemistry: molarity, dilution, pH."""
import math


def molarity(moles: float, volume_l: float) -> float:
    if volume_l <= 0:
        raise ValueError("Volume must be > 0")
    return moles / volume_l


def dilution(m1=None, v1=None, m2=None, v2=None) -> float:
    """M1*V1 = M2*V2. Leave exactly one argument as None to solve for it."""
    vals = {"m1": m1, "v1": v1, "m2": m2, "v2": v2}
    missing = [k for k, v in vals.items() if v is None]
    if len(missing) != 1:
        raise ValueError("Leave exactly one value empty")
    k = missing[0]
    if k == "m1": return m2 * v2 / v1
    if k == "v1": return m2 * v2 / m1
    if k == "m2": return m1 * v1 / v2
    return m1 * v1 / m2


def ph_from_h(h_conc: float) -> float:
    if h_conc <= 0:
        raise ValueError("[H+] must be > 0")
    return -math.log10(h_conc)


def ph_strong_acid(conc: float, protons: int = 1) -> float:
    return ph_from_h(conc * protons)


def ph_weak_acid(conc: float, ka: float) -> float:
    """Exact quadratic solution of Ka = x^2 / (C - x)."""
    x = (-ka + math.sqrt(ka ** 2 + 4 * ka * conc)) / 2
    return ph_from_h(x)


def ph_buffer(pka: float, base_conc: float, acid_conc: float) -> float:
    """Henderson-Hasselbalch: pH = pKa + log10([A-]/[HA])."""
    return pka + math.log10(base_conc / acid_conc)


def ph_strong_base(conc: float, hydroxides: int = 1) -> float:
    return 14.0 - (-math.log10(conc * hydroxides))


def ph_weak_base(conc: float, kb: float) -> float:
    """Exact quadratic solution of Kb = x^2 / (C - x), x = [OH-]."""
    x = (-kb + math.sqrt(kb ** 2 + 4 * kb * conc)) / 2
    return 14.0 + math.log10(x)
