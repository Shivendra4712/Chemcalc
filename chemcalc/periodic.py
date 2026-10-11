"""Periodic table data (Z = 1-54) and electron configurations."""
from typing import NamedTuple, Optional

_CATEGORIES = {"nm": "Nonmetal", "ng": "Noble gas", "am": "Alkali metal",
               "ae": "Alkaline earth metal", "tm": "Transition metal",
               "pt": "Post-transition metal", "md": "Metalloid", "ha": "Halogen"}

_RAW = """
1 H Hydrogen 2.20 1 1 nm
2 He Helium - 18 1 ng
3 Li Lithium 0.98 1 2 am
4 Be Beryllium 1.57 2 2 ae
5 B Boron 2.04 13 2 md
6 C Carbon 2.55 14 2 nm
7 N Nitrogen 3.04 15 2 nm
8 O Oxygen 3.44 16 2 nm
9 F Fluorine 3.98 17 2 ha
10 Ne Neon - 18 2 ng
11 Na Sodium 0.93 1 3 am
12 Mg Magnesium 1.31 2 3 ae
13 Al Aluminium 1.61 13 3 pt
14 Si Silicon 1.90 14 3 md
15 P Phosphorus 2.19 15 3 nm
16 S Sulfur 2.58 16 3 nm
17 Cl Chlorine 3.16 17 3 ha
18 Ar Argon - 18 3 ng
19 K Potassium 0.82 1 4 am
20 Ca Calcium 1.00 2 4 ae
21 Sc Scandium 1.36 3 4 tm
22 Ti Titanium 1.54 4 4 tm
23 V Vanadium 1.63 5 4 tm
24 Cr Chromium 1.66 6 4 tm
25 Mn Manganese 1.55 7 4 tm
26 Fe Iron 1.83 8 4 tm
27 Co Cobalt 1.88 9 4 tm
28 Ni Nickel 1.91 10 4 tm
29 Cu Copper 1.90 11 4 tm
30 Zn Zinc 1.65 12 4 tm
31 Ga Gallium 1.81 13 4 pt
32 Ge Germanium 2.01 14 4 md
33 As Arsenic 2.18 15 4 md
34 Se Selenium 2.55 16 4 nm
35 Br Bromine 2.96 17 4 ha
36 Kr Krypton 3.00 18 4 ng
37 Rb Rubidium 0.82 1 5 am
38 Sr Strontium 0.95 2 5 ae
39 Y Yttrium 1.22 3 5 tm
40 Zr Zirconium 1.33 4 5 tm
41 Nb Niobium 1.6 5 5 tm
42 Mo Molybdenum 2.16 6 5 tm
43 Tc Technetium 1.9 7 5 tm
44 Ru Ruthenium 2.2 8 5 tm
45 Rh Rhodium 2.28 9 5 tm
46 Pd Palladium 2.20 10 5 tm
47 Ag Silver 1.93 11 5 tm
48 Cd Cadmium 1.69 12 5 tm
49 In Indium 1.78 13 5 pt
50 Sn Tin 1.96 14 5 pt
51 Sb Antimony 2.05 15 5 md
52 Te Tellurium 2.1 16 5 md
53 I Iodine 2.66 17 5 ha
54 Xe Xenon 2.6 18 5 ng
"""


class Element(NamedTuple):
    z: int
    symbol: str
    name: str
    electronegativity: Optional[float]   # Pauling scale
    group: int
    period: int
    category: str


def _load():
    out = []
    for line in _RAW.strip().splitlines():
        z, sym, name, en, grp, per, cat = line.split()
        out.append(Element(int(z), sym, name, None if en == "-" else float(en),
                           int(grp), int(per), _CATEGORIES[cat]))
    return out


ELEMENTS = _load()
CATEGORY_NAMES = list(_CATEGORIES.values())


def get_element(key) -> Element:
    """Look up by symbol ('Fe'), name ('iron') or atomic number (26)."""
    k = str(key).strip().lower()
    for e in ELEMENTS:
        if k in (e.symbol.lower(), e.name.lower(), str(e.z)):
            return e
    raise ValueError(f"Element not found (supported: Z = 1-54): {key}")


_AUFBAU = [("1s", 2), ("2s", 2), ("2p", 6), ("3s", 2), ("3p", 6), ("4s", 2),
           ("3d", 10), ("4p", 6), ("5s", 2), ("4d", 10), ("5p", 6)]
_EXCEPTIONS = {24: {"4s": 1, "3d": 5}, 29: {"4s": 1, "3d": 10},
               41: {"5s": 1, "4d": 4}, 42: {"5s": 1, "4d": 5},
               44: {"5s": 1, "4d": 7}, 45: {"5s": 1, "4d": 8},
               46: {"5s": 0, "4d": 10}, 47: {"5s": 1, "4d": 10}}
_NOBLE = [(36, "Kr"), (18, "Ar"), (10, "Ne"), (2, "He")]


def _occupancy(z: int) -> dict:
    occ, left = {}, z
    for orb, cap in _AUFBAU:
        n = min(cap, left)
        if n:
            occ[orb] = n
        left -= n
        if left <= 0:
            break
    occ.update(_EXCEPTIONS.get(z, {}))
    return {k: v for k, v in occ.items() if v > 0}


def _fmt(occ: dict) -> str:
    key = lambda o: (int(o[0]), "spdf".index(o[1]))
    return " ".join(f"{o}{occ[o]}" for o in sorted(occ, key=key))


def electron_configuration(z: int, short: bool = False) -> str:
    """Ground-state configuration, e.g. 26 -> '1s2 2s2 ... 3d6 4s2'."""
    if not 1 <= z <= 54:
        raise ValueError("Supported range: Z = 1-54")
    occ = _occupancy(z)
    if short:
        for zc, sym in _NOBLE:
            if zc < z:
                core = _occupancy(zc)
                rest = {o: occ[o] - core.get(o, 0) for o in occ
                        if occ[o] - core.get(o, 0) > 0}
                return f"[{sym}] {_fmt(rest)}"
    return _fmt(occ)
