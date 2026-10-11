"""Acid-base titration curves via exact charge-balance solving."""
import numpy as np

KW = 1e-14


def titration_ph(ca: float, va_ml: float, cb: float, vb_ml: float, ka=None) -> float:
    """pH after adding vb_ml of strong base (cb, M) to va_ml of acid (ca, M).
    ka=None -> strong acid; otherwise monoprotic weak acid with that Ka.
    Solves [Na+] + [H+] = [A-] + [OH-] by bisection (exact, no approximations)."""
    vt = va_ml + vb_ml
    c_acid = ca * va_ml / vt
    c_na = cb * vb_ml / vt

    def f(h):
        a_minus = c_acid if ka is None else c_acid * ka / (ka + h)
        return c_na + h - a_minus - KW / h

    lo, hi = -15.0, 1.0  # log10([H+])
    for _ in range(100):
        mid = (lo + hi) / 2
        if f(10 ** mid) > 0:
            hi = mid
        else:
            lo = mid
    return -(lo + hi) / 2


def titration_curve(ca, va_ml, cb, vmax_ml, ka=None, points=300):
    vb = np.linspace(0, vmax_ml, points)
    ph = np.array([titration_ph(ca, va_ml, cb, v, ka) for v in vb])
    return vb, ph


def equivalence_volume(ca: float, va_ml: float, cb: float) -> float:
    return ca * va_ml / cb
