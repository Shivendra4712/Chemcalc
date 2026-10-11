"""Chemical equilibrium: reaction parsing, Q vs K, ICE-table solver."""
import math
import re


def parse_reaction(text: str):
    """'H2 + I2 <=> 2HI' -> ([(1,'H2'),(1,'I2')], [(2,'HI')])"""
    for tok in ("<=>", "<->", "\u21cc", "->", "=>", "="):
        text = text.replace(tok, "|")
    sides = text.split("|")
    if len(sides) != 2:
        raise ValueError("Write the reaction like: N2 + 3H2 <=> 2NH3")

    def terms(side):
        out = []
        for t in side.split("+"):
            t = t.strip()
            if not t:
                continue
            m = re.match(r"^(\d+)?\s*(.+)$", t)
            out.append((int(m.group(1) or 1), m.group(2).strip()))
        if not out:
            raise ValueError("Empty side in reaction")
        return out

    return terms(sides[0]), terms(sides[1])


def _ln_q(reactants, products, x):
    s = 0.0
    for coef, c0 in products:
        c = c0 + coef * x
        if c <= 0:
            return -math.inf
        s += coef * math.log(c)
    for coef, c0 in reactants:
        c = c0 - coef * x
        if c <= 0:
            return math.inf
        s -= coef * math.log(c)
    return s


def reaction_quotient(reactants, products, x: float = 0.0) -> float:
    """Q after extent x. reactants/products: lists of (coefficient, concentration)."""
    v = _ln_q(reactants, products, x)
    return 0.0 if v == -math.inf else math.inf if v == math.inf else math.exp(v)


def direction(q: float, k: float) -> str:
    if math.isclose(q, k, rel_tol=1e-6):
        return "At equilibrium"
    return "Shifts right (towards products)" if q < k else "Shifts left (towards reactants)"


def solve_equilibrium(reactants, products, k: float):
    """Find extent x so that Q(x) = K. Returns (x, reactant_conc, product_conc)."""
    if k <= 0:
        raise ValueError("K must be > 0")
    lo = -min((c0 / coef for coef, c0 in products), default=0.0)
    hi = min(c0 / coef for coef, c0 in reactants)
    if hi <= lo:
        raise ValueError("Need some initial reactant or product present")
    ln_k = math.log(k)
    for _ in range(300):
        mid = (lo + hi) / 2
        if _ln_q(reactants, products, mid) > ln_k:
            hi = mid
        else:
            lo = mid
    x = (lo + hi) / 2
    return (x, [c0 - coef * x for coef, c0 in reactants],
            [c0 + coef * x for coef, c0 in products])
