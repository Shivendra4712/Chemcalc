"""Limiting reagent, theoretical yield and percent yield."""
from .balancer import balance_species, format_equation
from .formula import molar_mass


def limiting_reagent(equation: str, masses_g: dict) -> dict:
    """Given reactant masses (g), find the limiting reagent and product yields."""
    reac, prod, co = balance_species(equation)
    coeff = dict(zip(reac + prod, co))
    missing = [s for s in reac if s not in masses_g]
    if missing:
        raise ValueError("Give a mass for: " + ", ".join(missing))
    moles = {s: masses_g[s] / molar_mass(s) for s in reac}
    extent = {s: moles[s] / coeff[s] for s in reac}
    limiting = min(extent, key=extent.get)
    xi = extent[limiting]
    return {
        "balanced": format_equation(reac, prod, co),
        "coefficients": coeff,
        "limiting": limiting,
        "extent": xi,
        "moles": moles,
        "products": {p: {"moles": coeff[p] * xi, "grams": coeff[p] * xi * molar_mass(p)}
                     for p in prod},
        "excess": {s: {"left_moles": moles[s] - coeff[s] * xi,
                       "left_grams": (moles[s] - coeff[s] * xi) * molar_mass(s)}
                   for s in reac if s != limiting},
    }


def percent_yield(actual: float, theoretical: float) -> float:
    if theoretical <= 0:
        raise ValueError("Theoretical yield must be > 0")
    return actual / theoretical * 100
