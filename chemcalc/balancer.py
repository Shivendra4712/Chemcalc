"""Chemical equation balancer using exact fraction linear algebra."""
from fractions import Fraction
from math import gcd
from functools import reduce
from .formula import parse_formula


def _nullspace(M):
    rows, cols = len(M), len(M[0])
    M = [r[:] for r in M]
    pivots, r = [], 0
    for c in range(cols):
        if r >= rows:
            break
        p = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
    basis = []
    for f in (c for c in range(cols) if c not in pivots):
        v = [Fraction(0)] * cols
        v[f] = Fraction(1)
        for row, pc in enumerate(pivots):
            v[pc] = -M[row][f]
        basis.append(v)
    return basis


def balance(equation: str) -> str:
    """'Fe + O2 -> Fe2O3'  ->  '4 Fe + 3 O2 -> 2 Fe2O3'"""
    sep = "->" if "->" in equation else "=" if "=" in equation else None
    if sep is None:
        raise ValueError("Use '->' between reactants and products")
    lhs, rhs = equation.split(sep)
    reac = [s.strip() for s in lhs.split("+") if s.strip()]
    prod = [s.strip() for s in rhs.split("+") if s.strip()]
    if not reac or not prod:
        raise ValueError("Need reactants and products")
    comps = [parse_formula(s) for s in reac + prod]
    elements = sorted({e for c in comps for e in c})
    M = [[Fraction(c.get(e, 0) * (1 if j < len(reac) else -1))
          for j, c in enumerate(comps)] for e in elements]
    basis = _nullspace(M)
    if len(basis) != 1:
        raise ValueError("Equation cannot be balanced uniquely")
    vec = basis[0]
    lcm = reduce(lambda a, b: a * b // gcd(a, b), (x.denominator for x in vec))
    ints = [int(x * lcm) for x in vec]
    if all(x < 0 for x in ints):
        ints = [-x for x in ints]
    if any(x <= 0 for x in ints):
        raise ValueError("No valid positive balance (check the species)")
    g = reduce(gcd, ints)
    ints = [x // g for x in ints]
    fmt = lambda names, co: " + ".join(f"{n} {s}" if n != 1 else s
                                       for n, s in zip(co, names))
    return f"{fmt(reac, ints[:len(reac)])} -> {fmt(prod, ints[len(reac):])}"
