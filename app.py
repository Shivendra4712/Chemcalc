"""ChemCalc - Streamlit GUI.  Run:  streamlit run app.py"""
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt

from chemcalc import (molar_mass, percent_composition, balance, dilution,
                      molarity, ph_strong_acid, ph_weak_acid, ph_buffer, ideal_gas)
from chemcalc import kinetics, molecular, titration

st.set_page_config(page_title="ChemCalc", page_icon="🧪", layout="wide")
st.title("🧪 ChemCalc — Computational Chemistry Toolkit")
st.caption("Calculators + simulations: stoichiometry, solutions, gases, titration, kinetics, molecular models")

tabs = st.tabs(["Formula", "Balancer", "Solutions & pH", "Gas Law",
                "Titration", "Kinetics", "Molecular"])

# ---------------------------------------------------------------- Formula
with tabs[0]:
    f = st.text_input("Chemical formula", "CuSO4.5H2O", key="formula")
    try:
        c1, c2 = st.columns(2)
        c1.metric("Molar mass", f"{molar_mass(f):.3f} g/mol")
        comp = percent_composition(f)
        c2.dataframe({"Element": list(comp), "Mass %": [round(v, 2) for v in comp.values()]})
        st.bar_chart(comp)
    except ValueError as e:
        st.error(str(e))

# --------------------------------------------------------------- Balancer
with tabs[1]:
    eq = st.text_input("Equation", "KMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2", key="eq")
    try:
        st.success(balance(eq))
    except ValueError as e:
        st.error(str(e))

# ------------------------------------------------------- Solutions & pH
with tabs[2]:
    mode = st.radio("Calculator", ["Dilution (M1V1=M2V2)", "Molarity",
                                   "pH strong acid", "pH weak acid", "Buffer pH"],
                    horizontal=True)
    try:
        if mode.startswith("Dilution"):
            solve = st.selectbox("Solve for", ["M1", "V1", "M2", "V2"], index=3)
            vals = {k: (None if k == solve else st.number_input(k, value=1.0, key="d" + k))
                    for k in ["M1", "V1", "M2", "V2"]}
            res = dilution(vals["M1"], vals["V1"], vals["M2"], vals["V2"])
            st.metric(solve, f"{res:.4g}")
        elif mode == "Molarity":
            n = st.number_input("Moles", value=0.5)
            v = st.number_input("Volume (L)", value=0.25)
            st.metric("Molarity", f"{molarity(n, v):.4g} M")
        elif mode == "pH strong acid":
            c = st.number_input("Concentration (M)", value=0.01, format="%.6f")
            p = st.selectbox("Protons released", [1, 2, 3])
            st.metric("pH", f"{ph_strong_acid(c, p):.3f}")
        elif mode == "pH weak acid":
            c = st.number_input("Concentration (M)", value=0.1, format="%.6f")
            ka = st.number_input("Ka", value=1.8e-5, format="%.3e")
            st.metric("pH", f"{ph_weak_acid(c, ka):.3f}")
        else:
            pka = st.number_input("pKa", value=4.76)
            a = st.number_input("[A-] (M)", value=0.1)
            ha = st.number_input("[HA] (M)", value=0.1)
            st.metric("pH", f"{ph_buffer(pka, a, ha):.3f}")
    except (ValueError, ZeroDivisionError) as e:
        st.error(str(e))

# ---------------------------------------------------------------- Gas law
with tabs[3]:
    solve = st.selectbox("Solve for", ["P (atm)", "V (L)", "n (mol)", "T (K)"], index=1)
    defaults = {"P (atm)": 1.0, "V (L)": 22.4, "n (mol)": 1.0, "T (K)": 273.15}
    vals = {k: (None if k == solve else st.number_input(k, value=defaults[k], key="g" + k))
            for k in defaults}
    try:
        res = ideal_gas(*vals.values())
        st.metric(solve, f"{res:.4f}")
        n = vals["n (mol)"] or (res if solve.startswith("n") else 1.0)
        t = vals["T (K)"] or (res if solve.startswith("T") else 273.15)
        v = np.linspace(5, 60, 200)
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(v, ideal_gas(v=v, n=n, t=t), lw=2)
        ax.set(xlabel="V (L)", ylabel="P (atm)", title=f"Isotherm (Boyle's law) n={n:.3g} mol, T={t:.3g} K")
        ax.grid(alpha=.3)
        st.pyplot(fig)
    except (ValueError, ZeroDivisionError) as e:
        st.error(str(e))

# -------------------------------------------------------------- Titration
with tabs[4]:
    c1, c2, c3 = st.columns(3)
    ca = c1.number_input("Acid conc (M)", value=0.1, format="%.4f")
    va = c2.number_input("Acid volume (mL)", value=25.0)
    cb = c3.number_input("NaOH conc (M)", value=0.1, format="%.4f")
    kind = st.radio("Acid type", ["Strong (HCl)", "Weak (HA)"], horizontal=True)
    ka = None
    if kind.startswith("Weak"):
        ka = st.number_input("Ka", value=1.8e-5, format="%.3e")
    veq = titration.equivalence_volume(ca, va, cb)
    vb, ph = titration.titration_curve(ca, va, cb, 2 * veq, ka)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(vb, ph, lw=2)
    ax.axvline(veq, ls="--", c="r", label=f"Equivalence = {veq:.2f} mL")
    ax.set(xlabel="Volume NaOH (mL)", ylabel="pH", title="Titration curve")
    ax.legend(); ax.grid(alpha=.3)
    st.pyplot(fig)
    st.metric("pH at equivalence", f"{titration.titration_ph(ca, va, cb, veq, ka):.2f}")

# --------------------------------------------------------------- Kinetics
with tabs[5]:
    st.subheader("Integrated rate laws")
    c1, c2, c3, c4 = st.columns(4)
    order = c1.selectbox("Order", [0, 1, 2], index=1)
    a0 = c2.number_input("[A]0 (M)", value=1.0)
    k = c3.number_input("k", value=0.05, format="%.4f")
    tmax = c4.number_input("Time range", value=100.0)
    t = np.linspace(0, tmax, 300)
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(t, kinetics.concentration(order, a0, k, t), lw=2)
    ax.set(xlabel="time", ylabel="[A] (M)"); ax.grid(alpha=.3)
    st.pyplot(fig)
    st.metric("Half-life", f"{kinetics.half_life(order, a0, k):.4g}")

    st.subheader("Consecutive reaction A → B → C")
    d1, d2 = st.columns(2)
    k1 = d1.number_input("k1", value=0.2, format="%.4f")
    k2 = d2.number_input("k2", value=0.05, format="%.4f")
    a, b, c = kinetics.consecutive(a0, k1, k2, t)
    fig, ax = plt.subplots(figsize=(7, 3.5))
    for y, lab in zip((a, b, c), "ABC"):
        ax.plot(t, y, lw=2, label=lab)
    ax.set(xlabel="time", ylabel="Conc (M)"); ax.legend(); ax.grid(alpha=.3)
    st.pyplot(fig)

    st.subheader("Arrhenius equation")
    e1, e2, e3 = st.columns(3)
    A = e1.number_input("A (pre-exponential)", value=1e10, format="%.3e")
    ea = e2.number_input("Ea (kJ/mol)", value=50.0)
    temp = e3.number_input("T (K)", value=300.0)
    st.metric("k", f"{kinetics.arrhenius(A, ea * 1000, temp):.4e}")
    T = np.linspace(250, 500, 200)
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(1000 / T, np.log(A * np.exp(-ea * 1000 / (kinetics.R * T))), lw=2)
    ax.set(xlabel="1000/T (1/K)", ylabel="ln k", title="Arrhenius plot (slope = -Ea/R)")
    ax.grid(alpha=.3)
    st.pyplot(fig)

# -------------------------------------------------------------- Molecular
with tabs[6]:
    st.subheader("Lennard-Jones potential")
    presets = {"Ne": (35.6, 2.75), "Ar": (119.8, 3.405), "Kr": (164.0, 3.65), "Xe": (221.0, 4.10)}
    gas = st.selectbox("Noble gas", list(presets))
    eps_k, sigma = presets[gas]
    eps = eps_k * molecular.KB_KJ  # kJ/mol
    r = np.linspace(0.9 * sigma, 3 * sigma, 400)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(r, molecular.lennard_jones(r, eps, sigma), lw=2)
    ax.axhline(0, c="k", lw=.5)
    ax.axvline(molecular.lj_minimum(sigma), ls="--", c="r",
               label=f"r_min = {molecular.lj_minimum(sigma):.2f} Å, depth = {eps:.3f} kJ/mol")
    ax.set(xlabel="r (Å)", ylabel="V(r) (kJ/mol)", ylim=(-1.5 * eps, 2 * eps),
           title=f"{gas}–{gas} interaction")
    ax.legend(); ax.grid(alpha=.3)
    st.pyplot(fig)

    st.subheader("Maxwell–Boltzmann speed distribution")
    g = st.selectbox("Gas", ["He", "H2", "N2", "O2", "CO2", "Ar"], index=2)
    temps = st.multiselect("Temperatures (K)", [100, 200, 300, 500, 800, 1200], default=[300, 800])
    M = molar_mass(g)
    v = np.linspace(0, 3000, 600)
    fig, ax = plt.subplots(figsize=(7, 4))
    for T in temps:
        ax.plot(v, molecular.mb_distribution(v, M, T), lw=2, label=f"{T} K")
    ax.set(xlabel="speed (m/s)", ylabel="f(v)", title=f"{g} (M = {M:.2f} g/mol)")
    ax.legend(); ax.grid(alpha=.3)
    st.pyplot(fig)
    if temps:
        st.dataframe({"T (K)": temps,
                      **{k: [round(molecular.mb_speeds(M, T)[k], 1) for T in temps]
                         for k in ("most_probable", "mean", "rms")}})
