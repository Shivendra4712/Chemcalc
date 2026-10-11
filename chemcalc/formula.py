"""Chemical formula parsing, molar mass and percent composition."""
import re
from collections import Counter
from .elements import ATOMIC_MASSES


def _parse_part(s: str) -> Counter:
    stack = [Counter()]
    i = 0
    while i < len(s):
        c = s[i]
        if c in "([":
            stack.append(Counter())
            i += 1
        elif c in ")]":
            if len(stack) < 2:
                raise ValueError("Unbalanced brackets in formula")
            i += 1
            m = re.match(r"\d+", s[i:])
            mult = int(m.group()) if m else 1
            i += len(m.group()) if m else 0
            group = stack.pop()
            for el, n in group.items():
                stack[-1][el] += n * mult
        elif c.isupper():
            sym = re.match(r"[A-Z][a-z]?", s[i:]).group()
            if sym not in ATOMIC_MASSES:
                raise ValueError(f"Unknown element: {sym}")
            i += len(sym)
            m = re.match(r"\d+", s[i:])
            n = int(m.group()) if m else 1
            i += len(m.group()) if m else 0
            stack[-1][sym] += n
        else:
            raise ValueError(f"Invalid character in formula: {c!r}")
    if len(stack) != 1:
        raise ValueError("Unbalanced brackets in formula")
    return stack[0]


def parse_formula(formula: str) -> dict:
    """'Ca(OH)2' -> {'Ca': 1, 'O': 2, 'H': 2}. Supports hydrates: CuSO4.5H2O"""
    formula = formula.replace(" ", "")
    if not formula:
        raise ValueError("Empty formula")
    total = Counter()
    for part in re.split(r"[\u00b7.*]", formula):
        m = re.match(r"^(\d*)(.*)$", part)
        coeff = int(m.group(1)) if m.group(1) else 1
        if not m.group(2):
            raise ValueError("Invalid formula")
        for el, n in _parse_part(m.group(2)).items():
            total[el] += n * coeff
    return dict(total)


def molar_mass(formula: str) -> float:
    return sum(ATOMIC_MASSES[el] * n for el, n in parse_formula(formula).items())


def percent_composition(formula: str) -> dict:
    mass = molar_mass(formula)
    return {el: ATOMIC_MASSES[el] * n / mass * 100
            for el, n in parse_formula(formula).items()}
