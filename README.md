# 🧪 ChemCalc — Computational Chemistry Toolkit

![tests](https://github.com/Shivendra4712/Chemcalc/actions/workflows/tests.yml/badge.svg)
![python](https://img.shields.io/badge/python-3.9%2B-blue)
![license](https://img.shields.io/badge/license-MIT-green)

A Python toolkit **and** a dark-themed Streamlit web app with interactive Plotly graphs, covering
stoichiometry, solutions, gases, thermodynamics, equilibrium, electrochemistry, titration, kinetics,
molecular simulation, the periodic table and drug-likeness analysis.

> **Live demo:** (https://chemcalcgitcdchemcalc-5vxdveckvardaqtsrbg8fm.streamlit.app/)

## ✨ Features

| Page | What you get |
|---|---|
| ⚗️ **Formula & Stoichiometry** | Molar mass (brackets, hydrates), composition pie chart, mole/particle conversion, equation balancer with atom check, **limiting reagent & % yield** |
| 🧪 **Solutions & pH** | Dilution solver, molarity from mass, strong/weak acids & bases, **live pH gauge**, buffers (Henderson–Hasselbalch) |
| 💨 **Gas Laws** | Ideal gas solver + isotherm, **van der Waals real gases** (Z factor, deviation), **3D P–V–T surface** |
| 🔥 **Thermo & Equilibrium** | ΔG = ΔH − TΔS, K, van 't Hoff, **ICE-table solver**, Q vs K direction, **Nernst equation & galvanic cells** |
| 📉 **Titration Lab** | Exact strong/weak acid curves, derivative plot, indicator ranges, CSV export |
| ⏱️ **Kinetics** | 0/1/2-order rate laws, **automatic order fitting from data**, A→B→C, Arrhenius & Ea |
| ⚛️ **Molecular Simulation** | Lennard-Jones potential/force, Maxwell–Boltzmann, **live 2D molecular dynamics** with RDF and energy conservation |
| 💠 **Periodic Table** | Interactive heatmap (Z = 1–54), electron configurations, electronegativity trends |
| 🧬 **Molecule Lab** | SMILES → 2D structure, descriptors, **Lipinski Rule of Five** (RDKit) |
| 🛠️ **Tools** | Unit converter, physical constants |

## 🚀 Quick start

```bash
git clone https://github.com/Shivendra4712/Chemcalc.git
cd Chemcalc
python -m pip install -r requirements.txt
python -m streamlit run app.py
```
Windows: just double-click `run_app.bat`.

CLI (no GUI needed):
```bash
python -m chemcalc                         # interactive menu
python -m chemcalc mass "Ca(OH)2"
python -m chemcalc balance "C3H8 + O2 -> CO2 + H2O"
python -m chemcalc element Fe
```

## 📚 Use as a library
```python
from chemcalc import molar_mass, balance
from chemcalc import thermo, equilibrium, titration, md, periodic
from chemcalc.stoich import limiting_reagent

balance("Fe + O2 -> Fe2O3")                          # 4 Fe + 3 O2 -> 2 Fe2O3
limiting_reagent("H2 + O2 -> H2O", {"H2": 4, "O2": 32})["limiting"]   # 'H2'
equilibrium.solve_equilibrium([(1, 1.0), (1, 1.0)], [(2, 0.0)], 50)   # ICE solver
periodic.electron_configuration(26, short=True)      # [Ar] 3d6 4s2
res = md.run_md(n_side=6, steps=2000)                # Lennard-Jones MD
```

## 🔬 Methods
- **Balancer** – nullspace of the element matrix using exact `fractions.Fraction` arithmetic.
- **Titration** – bisection on the full charge balance `[Na⁺] + [H⁺] = [A⁻] + [OH⁻]` (no approximations).
- **ICE solver** – bisection on ln Q(x) = ln K within the physically allowed extent range.
- **Molecular dynamics** – 2D Lennard-Jones fluid, velocity-Verlet, periodic boundaries, minimum-image
  convention, shifted potential, velocity-rescaling thermostat then NVE; energy drift reported.
- **Order fitting** – linear regression of [A], ln[A] and 1/[A] against time; best R² wins.
- **Equations** – PV = nRT · (P + a n²/V²)(V − nb) = nRT · ΔG = −RT ln K · E = E° − (RT/nF) ln Q ·
  k = A e^(−Ea/RT) · V(r) = 4ε[(σ/r)¹² − (σ/r)⁶] · f(v) = 4π(M/2πRT)^{3/2} v² e^(−Mv²/2RT)

## ✅ Tests
```bash
pytest        # 28 tests validated against known values
```
Examples: pH = pKa at half-equivalence, N₂ v_rms ≈ 517 m/s, Cu/Zn cell E° = 1.10 V,
H₂ + I₂ ⇌ 2HI equilibrium, MD energy conservation, aspirin Lipinski profile.
CI runs automatically on every push (GitHub Actions).

## ☁️ Deploy free on Streamlit Community Cloud
1. Push this repo to GitHub.
2. Go to [share.streamlit.io](https://share.streamlit.io) → **New app** → pick the repo, branch `main`, file `app.py`.
3. Click **Deploy** and paste the link in the *Live demo* line above.

## 🗂️ Structure
```
app.py                  Streamlit entry point (sidebar navigation)
ui/                     one module per page + shared theme helpers
chemcalc/               core package (pure Python + NumPy)
  formula · balancer · stoich · solutions · gas · thermo · equilibrium
  titration · kinetics · molecular · md · periodic · units · molecule · cli
tests/                  pytest suite
.streamlit/config.toml  dark theme
.github/workflows/      CI
```

## 🗺️ Roadmap
Real-gas critical constants · 3D molecular dynamics · quantum chemistry (PySCF) · more elements (Z > 54)

## 👤 Author
**Shivendra Pratap Singh** — B.Tech CSE (AI & ML), VIT Bhopal · MIT License
