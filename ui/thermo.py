import math
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from chemcalc import thermo as th
from chemcalc import equilibrium as eqm
from . import common as c


def render():
    c.header("\U0001f525", "Thermodynamics & Equilibrium",
             "Gibbs energy, ICE tables, Nernst equation and galvanic cells")
    t1, t2, t3 = st.tabs(["Gibbs energy", "Equilibrium (ICE)", "Electrochemistry"])

    with t1:
        a, b, cc = st.columns(3)
        dh = a.number_input("\u0394H (kJ/mol)", value=-92.2, key="t_dh")
        ds = b.number_input("\u0394S (J/mol\u00b7K)", value=-198.7, key="t_ds")
        t = cc.number_input("T (K)", value=298.15, min_value=1.0, key="t_t")
        dg = th.gibbs(dh, ds, t)
        k = th.k_from_dg(dg, t)
        m = st.columns(3)
        m[0].metric("\u0394G", f"{dg:.2f} kJ/mol")
        m[1].metric("K", f"{k:.3e}")
        m[2].metric("Spontaneous?", "Yes" if dg < 0 else "No")
        tc = th.crossover_temperature(dh, ds)
        if tc:
            st.info(f"\u0394G changes sign at T \u2248 {tc:.0f} K.")
        T = np.linspace(100, 1000, 200)
        fig = go.Figure(go.Scatter(x=T, y=[th.gibbs(dh, ds, x) for x in T], mode="lines", line=dict(width=3)))
        fig.add_hline(y=0, line_dash="dash", line_color="#94a3b8")
        fig.add_trace(go.Scatter(x=[t], y=[dg], mode="markers", marker=dict(size=13, color="#f59e0b"), name="T"))
        c.show(c.style(fig, "\u0394G vs temperature", "T (K)", "\u0394G (kJ/mol)", 360), "t_dg")
        t2_ = st.number_input("K at a second temperature T2 (K)", value=400.0, min_value=1.0, key="t_t2")
        st.metric(f"K at {t2_:.0f} K (van 't Hoff)", f"{th.van_t_hoff(k, t, dh, t2_):.3e}")

    with t2:
        rx_text = st.text_input("Reaction", "H2 + I2 <=> 2HI", key="t_rx")
        kc = st.number_input("Equilibrium constant K", value=50.0, min_value=1e-12, format="%.4g", key="t_k")
        try:
            rx, pr = eqm.parse_reaction(rx_text)
        except ValueError as e:
            st.error(str(e)); rx = None
        if rx:
            st.markdown("**Initial concentrations (M)**")
            cols = st.columns(len(rx) + len(pr))
            r0 = [(co, cols[i].number_input(sp, value=1.0, min_value=0.0, key=f"t_r{i}_{sp}"))
                  for i, (co, sp) in enumerate(rx)]
            p0 = [(co, cols[len(rx) + i].number_input(sp, value=0.0, min_value=0.0, key=f"t_p{i}_{sp}"))
                  for i, (co, sp) in enumerate(pr)]
            try:
                q = eqm.reaction_quotient(r0, p0)
                x, r_eq, p_eq = eqm.solve_equilibrium(r0, p0, kc)
            except (ValueError, ZeroDivisionError) as e:
                st.error(str(e))
            else:
                m = st.columns(3)
                m[0].metric("Q (initial)", f"{q:.4g}"); m[1].metric("K", f"{kc:.4g}")
                m[2].metric("Direction", eqm.direction(q, kc).split(" (")[0])
                names = [s for _, s in rx] + [s for _, s in pr]
                init = [v for _, v in r0] + [v for _, v in p0]
                chg = [f"\u2212{co}x" for co, _ in rx] + [f"+{co}x" for co, _ in pr]
                eq_ = r_eq + p_eq
                st.markdown(f"**ICE table** (x = {x:.5g})")
                st.dataframe(pd.DataFrame({"Species": names, "Initial": init, "Change": chg,
                                           "Equilibrium": [round(v, 6) for v in eq_]}), hide_index=True)
                fig = go.Figure([go.Bar(name="Initial", x=names, y=init),
                                 go.Bar(name="Equilibrium", x=names, y=eq_)])
                fig.update_layout(barmode="group")
                c.show(c.style(fig, "Concentrations", "", "M", 340), "t_ice")

    with t3:
        names = list(th.POTENTIALS)
        a, b = st.columns(2)
        cat = a.selectbox("Cathode (reduction)", names, index=names.index("Cu2+/Cu"), key="t_cat")
        an = b.selectbox("Anode (oxidation)", names, index=names.index("Zn2+/Zn"), key="t_an")
        cell = th.galvanic_cell(cat, an)
        m = st.columns(4)
        m[0].metric("E\u00b0 cell", f"{cell['E0']:.2f} V"); m[1].metric("n", cell["n"])
        m[2].metric("\u0394G\u00b0", f"{cell['dG_kJ']:.1f} kJ"); m[3].metric("K", f"{cell['K']:.2e}")
        st.success("Spontaneous galvanic cell") if cell["spontaneous"] else st.warning(
            "Not spontaneous in this direction (swap cathode/anode)")
        st.markdown("**Nernst equation**")
        lq = st.slider("log\u2081\u2080 Q", -6.0, 6.0, 0.0, key="t_lq")
        temp = st.number_input("T (K)", value=298.15, min_value=1.0, key="t_nt")
        st.metric("E cell", f"{th.nernst(cell['E0'], cell['n'], 10 ** lq, temp):.4f} V")
        x = np.linspace(-6, 6, 100)
        fig = go.Figure(go.Scatter(x=x, y=[th.nernst(cell["E0"], cell["n"], 10 ** v, temp) for v in x],
                                   mode="lines", line=dict(width=3)))
        c.show(c.style(fig, "E vs log Q", "log\u2081\u2080 Q", "E (V)", 340), "t_nernst")
