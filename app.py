"""ChemCalc - premium Streamlit app.   Run:  python -m streamlit run app.py"""
import streamlit as st

st.set_page_config(page_title="ChemCalc", page_icon="\U0001f9ea", layout="wide")

from chemcalc import __version__
from ui import common, home, formula, solutions, gas, thermo, titration, kinetics, molecular, periodic, molecule, tools

common.setup()

PAGES = {
    "\U0001f3e0 Home": home.render,
    "\u2697\ufe0f Formula & Stoichiometry": formula.render,
    "\U0001f9ea Solutions & pH": solutions.render,
    "\U0001f4a8 Gas Laws": gas.render,
    "\U0001f525 Thermo & Equilibrium": thermo.render,
    "\U0001f4c9 Titration Lab": titration.render,
    "\u23f1\ufe0f Kinetics": kinetics.render,
    "\u269b\ufe0f Molecular Simulation": molecular.render,
    "\U0001f4a0 Periodic Table": periodic.render,
    "\U0001f9ec Molecule Lab": molecule.render,
    "\U0001f6e0\ufe0f Tools": tools.render,
}

with st.sidebar:
    st.markdown("## \U0001f9ea ChemCalc")
    st.caption(f"Computational Chemistry Toolkit \u00b7 v{__version__}")
    choice = st.radio("Navigate", list(PAGES), label_visibility="collapsed")
    st.divider()
    st.caption("Built by Shivendra Pratap Singh")

PAGES[choice]()
