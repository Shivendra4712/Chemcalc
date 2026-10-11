"""Unit conversion and physical constants."""

FACTORS = {
    "Pressure": {"Pa": 1, "kPa": 1e3, "MPa": 1e6, "bar": 1e5, "atm": 101325.0,
                 "mmHg": 133.322368, "torr": 101325.0 / 760, "psi": 6894.757},
    "Energy": {"J": 1, "kJ": 1e3, "cal": 4.184, "kcal": 4184.0,
               "eV": 1.602176634e-19, "L\u00b7atm": 101.325, "kWh": 3.6e6},
    "Volume": {"L": 1, "mL": 1e-3, "m\u00b3": 1e3, "cm\u00b3": 1e-3, "dm\u00b3": 1},
    "Mass": {"g": 1, "kg": 1e3, "mg": 1e-3, "\u00b5g": 1e-6,
             "lb": 453.59237, "oz": 28.349523125},
    "Length": {"m": 1, "cm": 1e-2, "mm": 1e-3, "\u00b5m": 1e-6, "nm": 1e-9,
               "pm": 1e-12, "\u00c5": 1e-10, "in": 0.0254},
    "Temperature": {"K": None, "\u00b0C": None, "\u00b0F": None},
}

CONSTANTS = {
    "Gas constant R": (8.314462618, "J/(mol\u00b7K)"),
    "R (L\u00b7atm)": (0.082057366, "L\u00b7atm/(mol\u00b7K)"),
    "Avogadro N_A": (6.02214076e23, "1/mol"),
    "Boltzmann k_B": (1.380649e-23, "J/K"),
    "Faraday F": (96485.33212, "C/mol"),
    "Planck h": (6.62607015e-34, "J\u00b7s"),
    "Speed of light c": (299792458.0, "m/s"),
    "Atomic mass unit": (1.66053907e-27, "kg"),
}


def _to_kelvin(v, unit):
    return {"K": v, "\u00b0C": v + 273.15, "\u00b0F": (v - 32) * 5 / 9 + 273.15}[unit]


def _from_kelvin(k, unit):
    return {"K": k, "\u00b0C": k - 273.15, "\u00b0F": (k - 273.15) * 9 / 5 + 32}[unit]


def convert(value: float, frm: str, to: str, category: str) -> float:
    if category not in FACTORS:
        raise ValueError(f"Unknown category: {category}")
    if category == "Temperature":
        return _from_kelvin(_to_kelvin(value, frm), to)
    table = FACTORS[category]
    return value * table[frm] / table[to]
