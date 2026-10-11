import numpy as np
import plotly.graph_objects as go
import streamlit as st
from chemcalc import ideal_gas
from chemcalc.gas import R, VDW_CONSTANTS, van_der_waals_pressure, compressibility
from . import common as c


def render():
    c.header("\U0001f4a8", "Gas Laws", "Ideal gas, van der Waals real gases and P\u2013V\u2013T surfaces")
    t1, t2, t3 = st.tabs(["Ideal gas", "Real gas (van der Waals)", "P\u2013V\u2013T surface"])

    with t1:
        solve = st.selectbox("Solve for", ["P (atm)", "V (L)", "n (mol)", "T (K)"], index=1, key="g_solve")
        defaults = {"P (atm)": 1.0, "V (L)": 22.4, "n (mol)": 1.0, "T (K)": 273.15}
        cols = st.columns(4)
        vals = {k: (None if k == solve else col.number_input(k, value=d, min_value=1e-6, key="g_" + k))
                for (k, d), col in zip(defaults.items(), cols)}
        try:
            res = ideal_gas(*vals.values())
        except (ValueError, ZeroDivisionError) as e:
            st.error(str(e))
        else:
            st.metric(solve, f"{res:.4f}")
            n = vals["n (mol)"] if vals["n (mol)"] is not None else res
            t = vals["T (K)"] if vals["T (K)"] is not None else res
            v = np.linspace(1, 80, 300)
            fig = go.Figure(go.Scatter(x=v, y=ideal_gas(v=v, n=n, t=t), mode="lines", line=dict(width=3)))
            c.show(c.style(fig, f"Isotherm (Boyle's law): n={n:.3g} mol, T={t:.4g} K", "V (L)", "P (atm)"), "g_iso")

    with t2:
        a, b, cc, d = st.columns(4)
        gas = a.selectbox("Gas", list(VDW_CONSTANTS), index=6, key="g_vdw")
        n = b.number_input("n (mol)", value=1.0, min_value=0.01, key="g_vn")
        t = cc.number_input("T (K)", value=280.0, min_value=10.0, key="g_vt")
        vpt = d.number_input("V for readout (L)", value=1.0, min_value=0.05, key="g_vv")
        aa, bb = VDW_CONSTANTS[gas]
        vmin = max(n * bb * 1.05, 0.05)
        v = np.linspace(vmin, 6.0, 500)
        p_id, p_vdw = n * R * t / v, van_der_waals_pressure(v, n, t, aa, bb)
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=v, y=p_id, name="Ideal", line=dict(dash="dash", width=2)))
        fig.add_trace(go.Scatter(x=v, y=p_vdw, name="van der Waals", line=dict(width=3)))
        fig.update_yaxes(range=[min(0, float(np.nanmin(p_vdw))) * 1.2 if np.nanmin(p_vdw) < 0 else 0,
                                float(np.percentile(p_id, 90))])
        c.show(c.style(fig, f"{gas}: ideal vs real", "V (L)", "P (atm)"), "g_vdwplot")
        if vpt > n * bb:
            pi, pv = n * R * t / vpt, float(van_der_waals_pressure(vpt, n, t, aa, bb))
            m = st.columns(4)
            m[0].metric("P ideal", f"{pi:.3f} atm"); m[1].metric("P van der Waals", f"{pv:.3f} atm")
            m[2].metric("Deviation", f"{(pv - pi) / pi * 100:+.1f} %")
            m[3].metric("Z", f"{float(compressibility(pv, vpt, n, t)):.3f}")
        else:
            st.warning("V must be larger than n\u00b7b for the van der Waals equation.")

    with t3:
        n = st.slider("n (mol)", 0.1, 3.0, 1.0, key="g_sn")
        v = np.linspace(5, 50, 40)
        tt = np.linspace(200, 600, 40)
        V, T = np.meshgrid(v, tt)
        fig = go.Figure(go.Surface(x=v, y=tt, z=n * R * T / V, colorscale="Tealgrn"))
        fig.update_layout(scene=dict(xaxis_title="V (L)", yaxis_title="T (K)", zaxis_title="P (atm)"))
        c.show(c.style(fig, "P(V, T) = nRT/V", height=520), "g_surf")
