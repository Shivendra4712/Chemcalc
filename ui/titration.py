import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from chemcalc import titration as tt
from . import common as c

INDICATORS = {"None": None, "Phenolphthalein": (8.2, 10.0), "Methyl orange": (3.1, 4.4),
              "Methyl red": (4.4, 6.2), "Bromothymol blue": (6.0, 7.6)}


def render():
    c.header("\U0001f4c9", "Titration Lab", "Exact titration curves of an acid with NaOH (charge-balance solver)")
    a, b, cc = st.columns(3)
    ca = a.number_input("Acid conc (M)", value=0.1, min_value=1e-4, format="%.4f", key="ti_ca")
    va = b.number_input("Acid volume (mL)", value=25.0, min_value=1.0, key="ti_va")
    cb = cc.number_input("NaOH conc (M)", value=0.1, min_value=1e-4, format="%.4f", key="ti_cb")
    x, y = st.columns(2)
    kind = x.radio("Acid", ["Strong (HCl)", "Weak (HA)"], horizontal=True, key="ti_kind")
    ka = y.number_input("Ka", value=1.8e-5, format="%.3e", min_value=1e-14, key="ti_ka") \
        if kind.startswith("Weak") else None
    ind = st.selectbox("Indicator", list(INDICATORS), index=1, key="ti_ind")

    veq = tt.equivalence_volume(ca, va, cb)
    vb, ph = tt.titration_curve(ca, va, cb, 2 * veq, ka)
    m = st.columns(4)
    m[0].metric("Equivalence volume", f"{veq:.2f} mL")
    m[1].metric("Initial pH", f"{ph[0]:.2f}")
    m[2].metric("pH at half-equivalence", f"{tt.titration_ph(ca, va, cb, veq / 2, ka):.2f}")
    m[3].metric("pH at equivalence", f"{tt.titration_ph(ca, va, cb, veq, ka):.2f}")

    fig = go.Figure(go.Scatter(x=vb, y=ph, mode="lines", line=dict(width=3), name="pH"))
    fig.add_vline(x=veq, line_dash="dash", line_color="#f43f5e", annotation_text="equivalence")
    if INDICATORS[ind]:
        lo, hi = INDICATORS[ind]
        fig.add_hrect(y0=lo, y1=hi, fillcolor="#f59e0b", opacity=0.18, annotation_text=ind)
    c.show(c.style(fig, "Titration curve", "Volume NaOH (mL)", "pH"), "ti_curve")

    d = np.gradient(ph, vb)
    f2 = go.Figure(go.Scatter(x=vb, y=d, mode="lines", line=dict(width=3, color="#818cf8")))
    c.show(c.style(f2, "First derivative dpH/dV (peak = equivalence point)", "Volume NaOH (mL)", "dpH/dV", 300),
           "ti_deriv")
    c.csv_button(pd.DataFrame({"V_NaOH_mL": vb, "pH": ph}), "titration_curve.csv", "ti_csv")
