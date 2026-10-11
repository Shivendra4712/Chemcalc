# 🧪 ChemCalc — Computational Chemistry Toolkit

A Python toolkit with a **Streamlit web GUI** and a CLI for chemistry calculations and simulations.

## Features
| Module | What it does |
|---|---|
| `formula.py` | Formula parser (nested brackets, hydrates), molar mass, percent composition |
| `balancer.py` | Equation balancer using exact fraction linear algebra (nullspace) |
| `solutions.py` | Molarity, dilution (M1V1=M2V2), pH of strong/weak acids, buffers |
| `gas.py` | Ideal gas law solver (PV = nRT) |
| `titration.py` | Exact titration curves (strong & weak acid) via charge-balance root solving |
| `kinetics.py` | 0/1/2-order rate laws, half-life, Arrhenius, A→B→C consecutive reactions |
| `molecular.py` | Lennard-Jones potential & force, Maxwell–Boltzmann speed distribution |

## Setup
```bash
git clone https://github.com/<your-username>/chemcalc.git
cd chemcalc
python -m venv .venv
.venv\Scripts\activate        # Windows  (Mac/Linux: source .venv/bin/activate)
pip install -r requirements.txt
```

## Run
```bash
streamlit run app.py                     # GUI (7 tabs with live plots)
python -m chemcalc                       # CLI interactive menu
python -m chemcalc mass "Ca(OH)2"
python -m chemcalc balance "C3H8 + O2 -> CO2 + H2O"
```

## Use as a library
```python
from chemcalc import molar_mass, balance
from chemcalc import titration, kinetics, molecular

molar_mass("H2O")                                  # 18.015
balance("Fe + O2 -> Fe2O3")                        # 4 Fe + 3 O2 -> 2 Fe2O3
titration.titration_ph(0.1, 25, 0.1, 12.5, 1.8e-5) # 4.745 (= pKa at half-equivalence)
molecular.mb_speeds(28.014, 300)["rms"]            # ~517 m/s for N2
```

## Methods / theory
- **Balancer**: builds the element-by-species matrix and finds its nullspace with `fractions.Fraction` (no floating-point error).
- **Titration**: solves the charge balance `[Na+] + [H+] = [A-] + [OH-]` by bisection — no Henderson–Hasselbalch approximation.
- **Lennard-Jones**: `V(r) = 4ε[(σ/r)¹² − (σ/r)⁶]`, minimum at `r = 2^(1/6)σ` with depth `−ε`.
- **Maxwell–Boltzmann**: `f(v) = 4π (M/2πRT)^(3/2) v² exp(−Mv²/2RT)`.
- **Arrhenius**: `k = A·exp(−Ea/RT)`.

## Tests
```bash
pytest        # 16 tests: validated against known values (pKa at half-equivalence, N2 v_rms, LJ minimum, etc.)
```

## Structure
```
app.py                # Streamlit GUI
chemcalc/             # core package (formula, balancer, solutions, gas, titration, kinetics, molecular, cli)
tests/
requirements.txt
```

## Author
Shivendra Pratap Singh — B.Tech CSE (AI & ML), VIT Bhopal
