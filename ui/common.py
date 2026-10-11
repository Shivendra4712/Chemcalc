"""Shared UI helpers: theme CSS, headers, Plotly styling."""
import streamlit as st

PALETTE = ["#14b8a6", "#f59e0b", "#818cf8", "#f43f5e", "#22c55e", "#38bdf8"]

CSS = """
<style>
.block-container {padding-top: 1.4rem; max-width: 1250px;}
.hero {background: linear-gradient(135deg,#0f766e 0%,#1e3a8a 100%);
       padding: 1.5rem 1.8rem; border-radius: 16px; margin-bottom: 1.1rem;}
.hero h1 {margin: 0; font-size: 2rem; color: #fff; padding: 0;}
.hero p {margin: .3rem 0 0; color: #cbd5e1; font-size: 1rem;}
.card {background:#111a2e; border:1px solid #1f2b47; border-radius:14px;
       padding:1rem 1.1rem; min-height: 120px; margin-bottom: .8rem;}
.card h4 {margin:0 0 .35rem 0; color:#5eead4; font-size:1.05rem;}
.card p {margin:0; color:#94a3b8; font-size:.88rem; line-height:1.4;}
[data-testid="stMetric"] {background:#111a2e; border:1px solid #1f2b47;
       padding:.7rem 1rem; border-radius:12px;}
.pill {display:inline-block; padding:.15rem .7rem; border-radius:999px;
       background:#0f766e; color:#ecfeff; font-size:.8rem; margin-right:.4rem;}
</style>
"""


def setup():
    st.markdown(CSS, unsafe_allow_html=True)


def header(icon: str, title: str, subtitle: str):
    st.markdown(f'<div class="hero"><h1>{icon} {title}</h1><p>{subtitle}</p></div>',
                unsafe_allow_html=True)


def card(icon: str, title: str, text: str) -> str:
    return f'<div class="card"><h4>{icon} {title}</h4><p>{text}</p></div>'


def style(fig, title: str = "", x: str = "", y: str = "", height: int = 420):
    fig.update_layout(
        template="plotly_dark", title=title, height=height, colorway=PALETTE,
        paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="#0f172a",
        margin=dict(l=10, r=10, t=55 if title else 25, b=10),
        xaxis_title=x, yaxis_title=y,
        legend=dict(orientation="h", y=1.1, x=0))
    fig.update_xaxes(gridcolor="#1e293b", zeroline=False)
    fig.update_yaxes(gridcolor="#1e293b", zeroline=False)
    return fig


def show(fig, key: str):
    st.plotly_chart(fig, key=key)


def csv_button(df, filename: str, key: str):
    st.download_button("\u2b07 Download CSV", df.to_csv(index=False).encode(),
                       filename, "text/csv", key=key)
