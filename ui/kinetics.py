import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots
from chemcalc import kinetics as kn
from . import common as c

ORDER_LABEL = {0: "Zero order", 1: "First order", 2: "Second order"}


def render():
    c.header("\u23f1\ufe0f", "Chemical Kinetics", "Rate laws, order determination, consecutive reactions, Arrhenius")
    t1, t2, t3, t4 = st.tabs(["Rate laws", "Find the order", "A \u2192 B \u2192 C", "Arrhenius"])

    with t1:
        a, b, cc, d = st.columns(4)
        order = a.selectbox("Order", [0, 1, 2], index=1, key="k_order")
        a0 = b.number_input("[A]\u2080 (M)", value=1.0, min_value=1e-6, key="k_a0")
        k = cc.number_input("k", value=0.05, min_value=1e-9, format="%.4f", key="k_k")
        tmax = d.number_input("Time range", value=100.0, min_value=1.0, key="k_tmax")
        t = np.linspace(0, tmax, 300)
        fig = go.Figure()
        for o in (0, 1, 2):
            fig.add_trace(go.Scatter(x=t, y=kn.concentration(o, a0, k, t), name=ORDER_LABEL[o],
                                     line=dict(width=4 if o == order else 1.5)))
        th = kn.half_life(order, a0, k)
        fig.add_vline(x=th, line_dash="dot", annotation_text="t\u00bd")
        c.show(c.style(fig, "[A] vs time (selected order in bold)", "time", "[A] (M)"), "k_rate")
        st.metric(f"Half-life ({ORDER_LABEL[order].lower()})", f"{th:.4g}")

    with t2:
        st.caption("Edit the table (or paste your own data). The app fits all three linearised rate laws.")
        t = np.arange(0, 60, 10.0)
        df0 = pd.DataFrame({"Time": t, "Concentration": np.round(np.exp(-0.03 * t), 4)})
        df = st.data_editor(df0, num_rows="dynamic", key="k_table").dropna()
        try:
            best, out = kn.fit_order(df["Time"], df["Concentration"])
        except ValueError as e:
            st.error(str(e))
        else:
            m = st.columns(2)
            m[0].metric("Best fit", ORDER_LABEL[best])
            m[1].metric("Rate constant k", f"{out[best]['k']:.4g}")
            st.dataframe(pd.DataFrame([{"Order": ORDER_LABEL[o], "k": round(v["k"], 5), "R\u00b2": round(v["r2"], 5)}
                                       for o, v in out.items()]), hide_index=True)
            x = df["Time"].to_numpy(float); cn = df["Concentration"].to_numpy(float)
            ys = {0: ("[A]", cn), 1: ("ln [A]", np.log(cn)), 2: ("1/[A]", 1 / cn)}
            fig = make_subplots(rows=1, cols=3, subplot_titles=[f"{ORDER_LABEL[o]}  (R\u00b2={out[o]['r2']:.4f})" for o in ys])
            for i, (o, (lab, y)) in enumerate(ys.items(), 1):
                sl, ic = np.polyfit(x, y, 1)
                fig.add_trace(go.Scatter(x=x, y=y, mode="markers", marker=dict(size=9), showlegend=False), 1, i)
                fig.add_trace(go.Scatter(x=x, y=sl * x + ic, mode="lines", showlegend=False), 1, i)
                fig.update_yaxes(title_text=lab, row=1, col=i)
            c.show(c.style(fig, "", "", "", 340), "k_fit")

    with t3:
        a, b, cc = st.columns(3)
        k1 = a.number_input("k\u2081 (A\u2192B)", value=0.2, min_value=1e-6, format="%.4f", key="k_k1")
        k2 = b.number_input("k\u2082 (B\u2192C)", value=0.05, min_value=1e-6, format="%.4f", key="k_k2")
        a0 = cc.number_input("[A]\u2080 (M)", value=1.0, min_value=1e-6, key="k_ca0")
        t = np.linspace(0, 8 / min(k1, k2), 400)
        A, B, C = kn.consecutive(a0, k1, k2, t)
        fig = go.Figure()
        for y, lab in zip((A, B, C), "ABC"):
            fig.add_trace(go.Scatter(x=t, y=y, name=lab, line=dict(width=3)))
        c.show(c.style(fig, "Concentration profiles", "time", "M"), "k_abc")
        if not np.isclose(k1, k2):
            tm = np.log(k1 / k2) / (k1 - k2)
            st.metric("B peaks at t =", f"{tm:.3g}  ([B]max = {kn.consecutive(a0, k1, k2, tm)[1]:.3f} M)")

    with t4:
        a, b, cc = st.columns(3)
        A_ = a.number_input("A (pre-exponential)", value=1e10, format="%.3e", min_value=1.0, key="k_A")
        ea = b.number_input("Ea (kJ/mol)", value=50.0, min_value=0.1, key="k_ea")
        T = cc.number_input("T (K)", value=300.0, min_value=50.0, key="k_T")
        st.metric("k = A\u00b7exp(\u2212Ea/RT)", f"{kn.arrhenius(A_, ea * 1000, T):.4e}")
        Ts = np.linspace(250, 600, 200)
        lnk = np.log(A_) - ea * 1000 / (kn.R * Ts)
        fig = go.Figure(go.Scatter(x=1000 / Ts, y=lnk, mode="lines", line=dict(width=3)))
        c.show(c.style(fig, "Arrhenius plot (slope = \u2212Ea/R)", "1000/T (1/K)", "ln k", 340), "k_arr")
        st.markdown("**Find Ea from two rate constants**")
        x, y, z, w = st.columns(4)
        k1_ = x.number_input("k\u2081", value=1.5e-3, format="%.3e", min_value=1e-12, key="k_ea_k1")
        t1_ = y.number_input("T\u2081 (K)", value=300.0, min_value=1.0, key="k_ea_t1")
        k2_ = z.number_input("k\u2082", value=1.2e-2, format="%.3e", min_value=1e-12, key="k_ea_k2")
        t2_ = w.number_input("T\u2082 (K)", value=320.0, min_value=1.0, key="k_ea_t2")
        if t1_ != t2_:
            st.metric("Ea", f"{kn.activation_energy(k1_, t1_, k2_, t2_) / 1000:.2f} kJ/mol")
