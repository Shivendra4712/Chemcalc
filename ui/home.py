import streamlit as st
from chemcalc import molar_mass, __version__
from chemcalc.periodic import ELEMENTS
from chemcalc.thermo import POTENTIALS
from chemcalc.gas import VDW_CONSTANTS
from . import common as c

FEATURES = [
    ("\u2697\ufe0f", "Stoichiometry", "Molar mass, equation balancer (exact nullspace), limiting reagent & % yield"),
    ("\U0001f9ea", "Solutions & pH", "Dilution, molarity from mass, strong/weak acids & bases, buffers with live pH gauge"),
    ("\U0001f4a8", "Gas Laws", "Ideal gas solver, van der Waals real gases, 3D P\u2013V\u2013T surface"),
    ("\U0001f525", "Thermo & Equilibrium", "Gibbs energy, K, van 't Hoff, ICE-table solver, Nernst & galvanic cells"),
    ("\U0001f4c9", "Titration", "Exact curves (strong/weak acid), derivative plot, indicator ranges, CSV export"),
    ("\u23f1\ufe0f", "Kinetics", "Rate laws, automatic order fitting, A\u2192B\u2192C, Arrhenius & Ea"),
    ("\u269b\ufe0f", "Molecular Simulation", "Lennard-Jones, Maxwell\u2013Boltzmann, live 2D molecular dynamics with RDF"),
    ("\U0001f4a0", "Periodic Table", "Interactive heatmap, electron configurations, periodic trends"),
    ("\U0001f9ec", "Molecule Lab", "SMILES \u2192 2D structure, descriptors & Lipinski drug-likeness (RDKit)"),
    ("\U0001f6e0\ufe0f", "Tools", "Unit converter and physical constants"),
]


def render():
    c.header("\U0001f9ea", "ChemCalc", f"Computational chemistry toolkit \u00b7 v{__version__}")
    m = st.columns(3)
    m[0].metric("Elements in database", len(ELEMENTS))
    m[1].metric("Standard half-cells", len(POTENTIALS))
    m[2].metric("Real-gas (vdW) species", len(VDW_CONSTANTS))

    st.markdown("#### What's inside")
    cols = st.columns(3)
    for i, (icon, title, text) in enumerate(FEATURES):
        cols[i % 3].markdown(c.card(icon, title, text), unsafe_allow_html=True)

    st.markdown("#### Quick molar mass")
    f = st.text_input("Formula", "C6H12O6", key="home_formula")
    try:
        st.metric(f, f"{molar_mass(f):.3f} g/mol")
    except ValueError as e:
        st.error(str(e))
    st.caption("Pick a module from the sidebar. All core calculations are pure Python "
               "(NumPy) and covered by unit tests.")
