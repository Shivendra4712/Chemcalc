"""Ideal gas law: PV = nRT (P in atm, V in L, n in mol, T in K)."""
R = 0.082057  # L*atm/(mol*K)


def ideal_gas(p=None, v=None, n=None, t=None) -> float:
    vals = {"p": p, "v": v, "n": n, "t": t}
    missing = [k for k, x in vals.items() if x is None]
    if len(missing) != 1:
        raise ValueError("Leave exactly one value empty")
    for k, x in vals.items():
        if x is not None and x <= 0:
            raise ValueError(f"{k.upper()} must be > 0")
    k = missing[0]
    if k == "p": return n * R * t / v
    if k == "v": return n * R * t / p
    if k == "n": return p * v / (R * t)
    return p * v / (n * R)
