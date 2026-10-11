import pandas as pd
import streamlit as st
from chemcalc.units import FACTORS, CONSTANTS, convert
from . import common as c


def render():
    c.header("\U0001f6e0\ufe0f", "Tools", "Unit converter and physical constants")
    t1, t2 = st.tabs(["Unit converter", "Constants"])
    with t1:
        cat = st.selectbox("Quantity", list(FACTORS), key="u_cat")
        units = list(FACTORS[cat])
        a, b, cc = st.columns(3)
        v = a.number_input("Value", value=1.0, key="u_v")
        frm = b.selectbox("From", units, key=f"u_from_{cat}")
        to = cc.selectbox("To", units, index=min(1, len(units) - 1), key=f"u_to_{cat}")
        st.metric(f"{v:g} {frm} =", f"{convert(v, frm, to, cat):.6g} {to}")
    with t2:
        st.dataframe(pd.DataFrame([{"Constant": k, "Value": f"{v:.6e}", "Unit": u}
                                   for k, (v, u) in CONSTANTS.items()]), hide_index=True)
