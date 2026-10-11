import numpy as np
import plotly.graph_objects as go
import streamlit as st
from chemcalc import molar_mass
from chemcalc import molecular as mol
from chemcalc.md import run_md
from . import common as c

PRESETS = {"Ne": (35.6, 2.75), "Ar": (119.8, 3.405), "Kr": (164.0, 3.65), "Xe": (221.0, 4.10)}


@st.cache_data(show_spinner="Running molecular dynamics\u2026")
def _run(n_side, density, temperature, steps, thermostat):
    return run_md(n_side=n_side, density=density, temperature=temperature,
                  steps=steps, thermostat=thermostat)


def render():
    c.header("\u269b\ufe0f", "Molecular Simulation",
             "Intermolecular potentials, gas speed distributions and live molecular dynamics")
    t1, t2, t3 = st.tabs(["Lennard-Jones", "Maxwell\u2013Boltzmann", "Molecular dynamics (2D)"])

    with t1:
        gas = st.selectbox("Noble gas", list(PRESETS), index=1, key="m_gas")
        eps_k, sigma = PRESETS[gas]
        eps = eps_k * mol.KB_KJ
        r = np.linspace(0.9 * sigma, 3 * sigma, 400)
        rmin = mol.lj_minimum(sigma)
        m = st.columns(3)
        m[0].metric("\u03c3", f"{sigma} \u00c5"); m[1].metric("r_min", f"{rmin:.3f} \u00c5")
        m[2].metric("Well depth \u03b5", f"{eps:.3f} kJ/mol")
        a, b = st.columns(2)
        f1 = go.Figure(go.Scatter(x=r, y=mol.lennard_jones(r, eps, sigma), mode="lines", line=dict(width=3)))
        f1.add_vline(x=rmin, line_dash="dash", line_color="#f43f5e")
        f1.update_yaxes(range=[-1.5 * eps, 2 * eps])
        with a:
            c.show(c.style(f1, "Potential V(r)", "r (\u00c5)", "kJ/mol", 360), "m_lj")
        f2 = go.Figure(go.Scatter(x=r, y=mol.lj_force(r, eps, sigma), mode="lines", line=dict(width=3, color="#f59e0b")))
        f2.update_yaxes(range=[-3 * eps / sigma, 3 * eps / sigma])
        with b:
            c.show(c.style(f2, "Force F(r) = \u2212dV/dr", "r (\u00c5)", "kJ/mol/\u00c5", 360), "m_ljf")

    with t2:
        a, b = st.columns(2)
        g = a.selectbox("Gas", ["He", "H2", "N2", "O2", "CO2", "Ar"], index=2, key="m_mbgas")
        temps = b.multiselect("Temperatures (K)", [100, 200, 300, 500, 800, 1200], default=[300, 800], key="m_temps")
        M = molar_mass(g)
        v = np.linspace(0, 3000, 600)
        fig = go.Figure()
        for T in temps:
            fig.add_trace(go.Scatter(x=v, y=mol.mb_distribution(v, M, T), name=f"{T} K", line=dict(width=3)))
        c.show(c.style(fig, f"{g} speed distribution (M = {M:.2f} g/mol)", "speed (m/s)", "f(v)"), "m_mb")
        if temps:
            import pandas as pd
            st.dataframe(pd.DataFrame({"T (K)": temps,
                                       **{k: [round(mol.mb_speeds(M, T)[k], 1) for T in temps]
                                          for k in ("most_probable", "mean", "rms")}}), hide_index=True)

    with t3:
        st.caption("Argon-like Lennard-Jones fluid in reduced units (\u03b5 = \u03c3 = m = k_B = 1), velocity-Verlet, "
                   "periodic boundaries. Thermostat for the first half, then energy-conserving NVE.")
        with st.form("md_form"):
            a, b, cc, d = st.columns(4)
            n_side = a.slider("Particles per side", 4, 10, 6)
            dens = b.slider("Density \u03c1*", 0.1, 0.9, 0.7)
            temp = cc.slider("Temperature T*", 0.3, 2.5, 1.0)
            steps = d.slider("Steps", 500, 4000, 1500, step=250)
            st.form_submit_button("\u25b6 Run simulation")
        res = _run(n_side, dens, temp, steps, True)
        m = st.columns(4)
        m[0].metric("Particles", res.n)
        m[1].metric("Mean T* (2nd half)", f"{res.temperature[res.thermostat_steps:].mean():.3f}")
        m[2].metric("Mean PE / particle", f"{res.potential[res.thermostat_steps:].mean() / res.n:.3f}")
        m[3].metric("NVE energy drift / particle", f"{res.energy_drift_per_particle:.2e}")
        cols = st.columns(len(res.snapshots))
        for col, (step, pos) in zip(cols, res.snapshots):
            f = go.Figure(go.Scatter(x=pos[:, 0], y=pos[:, 1], mode="markers",
                                     marker=dict(size=14, color="#14b8a6", line=dict(width=1, color="#042f2e"))))
            f.update_xaxes(range=[0, res.box], showgrid=False)
            f.update_yaxes(range=[0, res.box], showgrid=False, scaleanchor="x")
            with col:
                c.show(c.style(f, f"step {step}", height=300), f"m_snap{step}")
        f1 = go.Figure()
        for y, lab in ((res.kinetic, "Kinetic"), (res.potential, "Potential"), (res.total, "Total")):
            f1.add_trace(go.Scatter(x=res.time, y=y, name=lab, line=dict(width=2.5)))
        f1.add_vrect(x0=0, x1=res.thermostat_steps * res.dt, fillcolor="#f59e0b", opacity=0.1,
                     annotation_text="thermostat on")
        a, b = st.columns(2)
        with a:
            c.show(c.style(f1, "Energy", "time (reduced)", "E", 340), "m_energy")
        f2 = go.Figure(go.Scatter(x=res.rdf_r, y=res.rdf_g, mode="lines", line=dict(width=3, color="#818cf8")))
        f2.add_hline(y=1, line_dash="dash", line_color="#94a3b8")
        with b:
            c.show(c.style(f2, "Radial distribution function g(r)", "r (\u03c3)", "g(r)", 340), "m_rdf")
