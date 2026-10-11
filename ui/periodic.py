import numpy as np
import plotly.graph_objects as go
import streamlit as st
from chemcalc.elements import ATOMIC_MASSES
from chemcalc.periodic import ELEMENTS, CATEGORY_NAMES, get_element, electron_configuration
from . import common as c

CAT_COLORS = dict(zip(CATEGORY_NAMES, ["#22c55e", "#818cf8", "#f43f5e", "#f59e0b", "#38bdf8", "#14b8a6", "#a3e635", "#e879f9"]))


def render():
    c.header("\U0001f4a0", "Periodic Table", "Interactive heatmap, electron configurations and trends (Z = 1\u201354)")
    prop = st.radio("Colour by", ["Electronegativity", "Atomic mass", "Atomic number"], horizontal=True, key="p_prop")
    getv = {"Electronegativity": lambda e: e.electronegativity,
            "Atomic mass": lambda e: ATOMIC_MASSES[e.symbol],
            "Atomic number": lambda e: e.z}[prop]
    z = np.full((5, 18), np.nan); text = [[""] * 18 for _ in range(5)]; hover = [[""] * 18 for _ in range(5)]
    for e in ELEMENTS:
        r, col = e.period - 1, e.group - 1
        v = getv(e)
        z[r, col] = np.nan if v is None else v
        text[r][col] = e.symbol
        hover[r][col] = f"{e.name} (Z={e.z})<br>{prop}: {v if v is not None else 'n/a'}<br>{e.category}"
    fig = go.Figure(go.Heatmap(z=z, x=list(range(1, 19)), y=[1, 2, 3, 4, 5], text=text, texttemplate="%{text}",
                               hovertext=hover, hoverinfo="text", colorscale="Tealgrn",
                               xgap=3, ygap=3, colorbar=dict(title=prop)))
    fig.update_yaxes(autorange="reversed", title="Period", dtick=1)
    fig.update_xaxes(title="Group", dtick=1)
    c.show(c.style(fig, "", height=420), "p_heat")

    labels = [f"{e.symbol} \u2013 {e.name}" for e in ELEMENTS]
    sel = st.selectbox("Element details", labels, index=25, key="p_sel")
    e = get_element(sel.split(" ")[0])
    m = st.columns(4)
    m[0].metric("Atomic number", e.z); m[1].metric("Atomic mass", f"{ATOMIC_MASSES[e.symbol]} g/mol")
    m[2].metric("Electronegativity", e.electronegativity if e.electronegativity else "n/a")
    m[3].metric("Group / Period", f"{e.group} / {e.period}")
    st.markdown(f'<span class="pill">{e.category}</span>', unsafe_allow_html=True)
    st.code(f"Short: {electron_configuration(e.z, short=True)}\nFull : {electron_configuration(e.z)}")

    st.markdown("#### Periodic trend: electronegativity")
    fig = go.Figure()
    for cat in CATEGORY_NAMES:
        els = [x for x in ELEMENTS if x.category == cat and x.electronegativity]
        fig.add_trace(go.Scatter(x=[x.z for x in els], y=[x.electronegativity for x in els], mode="markers+text",
                                 text=[x.symbol for x in els], textposition="top center", name=cat,
                                 marker=dict(size=10, color=CAT_COLORS[cat])))
    c.show(c.style(fig, "", "Atomic number Z", "Pauling electronegativity", 420), "p_trend")
