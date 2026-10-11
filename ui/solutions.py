import numpy as np
import plotly.graph_objects as go
import streamlit as st
from chemcalc import (molar_mass, molarity, dilution, ph_strong_acid, ph_weak_acid, ph_buffer)
from chemcalc.solutions import ph_strong_base, ph_weak_base
from . import common as c


def gauge(ph: float, key: str):
    fig = go.Figure(go.Indicator(
        mode="gauge+number", value=ph, number={"valueformat": ".2f"},
        gauge={"axis": {"range": [0, 14]}, "bar": {"color": "#f8fafc", "thickness": 0.25},
               "steps": [{"range": [0, 3], "color": "#ef4444"}, {"range": [3, 6], "color": "#f59e0b"},
                         {"range": [6, 8], "color": "#22c55e"}, {"range": [8, 11], "color": "#38bdf8"},
                         {"range": [11, 14], "color": "#6366f1"}]}))
    fig.update_layout(height=250, paper_bgcolor="rgba(0,0,0,0)", font={"color": "#e2e8f0"},
                      margin=dict(l=30, r=30, t=30, b=0))
    c.show(fig, key)


def render():
    c.header("\U0001f9ea", "Solutions & pH", "Dilution, molarity, acids/bases and buffers")
    t1, t2, t3, t4 = st.tabs(["Dilution", "Molarity", "pH & pOH", "Buffer"])

    with t1:
        solve = st.radio("Solve for", ["M1", "V1", "M2", "V2"], index=3, horizontal=True, key="s_solve")
        cols = st.columns(4)
        vals = {}
        for k, col in zip(["M1", "V1", "M2", "V2"], cols):
            vals[k] = None if k == solve else col.number_input(k, value=1.0, min_value=1e-9, key="s_d" + k)
        try:
            st.metric(solve, f"{dilution(vals['M1'], vals['V1'], vals['M2'], vals['V2']):.5g}")
        except (ValueError, ZeroDivisionError) as e:
            st.error(str(e))
        st.caption("M1\u00b7V1 = M2\u00b7V2 (any consistent volume unit)")

    with t2:
        a, b = st.columns(2)
        with a:
            st.markdown("**Molarity from mass**")
            f = st.text_input("Solute formula", "NaCl", key="s_f")
            g = st.number_input("Mass (g)", value=5.85, min_value=0.0, key="s_g")
            v = st.number_input("Solution volume (mL)", value=100.0, min_value=1e-6, key="s_v")
            try:
                mm = molar_mass(f)
                st.metric("Molarity", f"{molarity(g / mm, v / 1000):.4f} M")
                st.caption(f"{g / mm:.4f} mol of {f} ({mm:.3f} g/mol)")
            except ValueError as e:
                st.error(str(e))
        with b:
            st.markdown("**Mass needed for a target molarity**")
            m = st.number_input("Target (M)", value=0.5, min_value=0.0, key="s_tm")
            vv = st.number_input("Volume (mL)", value=250.0, min_value=1e-6, key="s_tv")
            try:
                st.metric("Mass to weigh", f"{m * vv / 1000 * molar_mass(f):.3f} g")
            except ValueError:
                pass

    with t3:
        kind = st.radio("Type", ["Strong acid", "Strong base", "Weak acid", "Weak base"],
                        horizontal=True, key="s_kind")
        conc = st.number_input("Concentration (M)", value=0.1, min_value=1e-12, format="%.6f", key="s_conc")
        try:
            if kind == "Strong acid":
                n = st.selectbox("Protons released", [1, 2, 3], key="s_np")
                ph = ph_strong_acid(conc, n)
            elif kind == "Strong base":
                n = st.selectbox("OH\u207b released", [1, 2, 3], key="s_nb")
                ph = ph_strong_base(conc, n)
            elif kind == "Weak acid":
                ka = st.number_input("Ka", value=1.8e-5, format="%.3e", min_value=1e-20, key="s_ka")
                ph = ph_weak_acid(conc, ka)
                st.metric("% ionisation", f"{10 ** -ph / conc * 100:.3f} %")
            else:
                kb = st.number_input("Kb", value=1.8e-5, format="%.3e", min_value=1e-20, key="s_kb")
                ph = ph_weak_base(conc, kb)
            m1, m2, m3 = st.columns(3)
            m1.metric("pH", f"{ph:.3f}"); m2.metric("pOH", f"{14 - ph:.3f}")
            m3.metric("[H\u207a]", f"{10 ** -ph:.3e} M")
            gauge(ph, "s_gauge")
        except (ValueError, ZeroDivisionError) as e:
            st.error(str(e))

    with t4:
        presets = {"Acetic acid / acetate": 4.76, "Formic acid / formate": 3.75,
                   "Carbonic / bicarbonate": 6.35, "H2PO4\u207b / HPO4\u00b2\u207b": 7.21,
                   "Ammonium / ammonia": 9.25, "Custom": 4.76}
        p = st.selectbox("Buffer system", list(presets), key="s_buf")
        pka = st.number_input("pKa", value=presets[p], key="s_pka")
        x, y = st.columns(2)
        a_ = x.number_input("[A\u207b] conjugate base (M)", value=0.1, min_value=1e-9, key="s_a")
        ha = y.number_input("[HA] weak acid (M)", value=0.1, min_value=1e-9, key="s_ha")
        st.metric("Buffer pH", f"{ph_buffer(pka, a_, ha):.3f}")
        ratio = np.logspace(-2, 2, 200)
        fig = go.Figure(go.Scatter(x=ratio, y=pka + np.log10(ratio), mode="lines", line=dict(width=3)))
        fig.add_vrect(x0=0.1, x1=10, fillcolor="#14b8a6", opacity=0.12, annotation_text="effective range")
        fig.add_trace(go.Scatter(x=[a_ / ha], y=[ph_buffer(pka, a_, ha)], mode="markers",
                                 marker=dict(size=13, color="#f59e0b"), name="your buffer"))
        fig.update_xaxes(type="log")
        c.show(c.style(fig, "Henderson\u2013Hasselbalch", "[A\u207b]/[HA]", "pH", 360), "s_hh")
