import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from chemcalc import parse_formula, molar_mass, percent_composition
from chemcalc.balancer import balance_species, format_equation
from chemcalc.elements import ATOMIC_MASSES
from chemcalc.stoich import limiting_reagent, percent_yield
from . import common as c

NA = 6.02214076e23


def render():
    c.header("\u2697\ufe0f", "Formula & Stoichiometry",
             "Molar mass, composition, equation balancing, limiting reagent")
    t1, t2, t3 = st.tabs(["Molar mass", "Equation balancer", "Limiting reagent & yield"])

    with t1:
        f = st.text_input("Chemical formula", "CuSO4.5H2O", key="f_formula")
        st.caption("Supports brackets and hydrates: Al2(SO4)3, CuSO4.5H2O")
        try:
            counts, mm, comp = parse_formula(f), molar_mass(f), percent_composition(f)
        except ValueError as e:
            st.error(str(e))
        else:
            a, b = st.columns([1, 2])
            a.metric("Molar mass", f"{mm:.3f} g/mol")
            grams = a.number_input("Mass (g)", value=10.0, min_value=0.0, key="f_grams")
            mol = grams / mm
            a.metric("Moles", f"{mol:.4f} mol")
            a.metric("Particles", f"{mol * NA:.3e}")
            df = pd.DataFrame({"Element": list(counts), "Atoms": list(counts.values()),
                               "Mass (g/mol)": [round(ATOMIC_MASSES[e] * n, 3) for e, n in counts.items()],
                               "Mass %": [round(comp[e], 2) for e in counts]})
            b.dataframe(df, hide_index=True)
            fig = go.Figure(go.Pie(labels=df["Element"], values=df["Mass %"], hole=0.55))
            c.show(c.style(fig, "Mass composition", height=330), "f_pie")

    with t2:
        eq = st.text_input("Unbalanced equation", "KMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2", key="f_eq")
        try:
            reac, prod, co = balance_species(eq)
        except ValueError as e:
            st.error(str(e))
        else:
            st.success(format_equation(reac, prod, co))
            n = len(reac)
            rows = {}
            for sp, k in zip(reac, co[:n]):
                for el, a in parse_formula(sp).items():
                    rows.setdefault(el, [0, 0])[0] += a * k
            for sp, k in zip(prod, co[n:]):
                for el, a in parse_formula(sp).items():
                    rows.setdefault(el, [0, 0])[1] += a * k
            st.markdown("**Atom check**")
            st.dataframe(pd.DataFrame([{"Element": e, "Left": l, "Right": r, "Balanced": "\u2705" if l == r else "\u274c"}
                                       for e, (l, r) in rows.items()]), hide_index=True)

    with t3:
        eq = st.text_input("Reaction", "N2 + H2 -> NH3", key="f_lim_eq")
        try:
            reac, prod, co = balance_species(eq)
        except ValueError as e:
            st.error(str(e))
            return
        st.info("Balanced: " + format_equation(reac, prod, co))
        cols = st.columns(len(reac))
        defaults = {"N2": 28.0, "H2": 10.0}
        masses = {sp: col.number_input(f"Mass of {sp} (g)", value=defaults.get(sp, 10.0),
                                       min_value=0.0, key=f"f_m_{sp}")
                  for sp, col in zip(reac, cols)}
        try:
            r = limiting_reagent(eq, masses)
        except ValueError as e:
            st.error(str(e))
            return
        st.metric("Limiting reagent", r["limiting"])
        rows = [{"Species": p, "Theoretical moles": round(v["moles"], 4),
                 "Theoretical mass (g)": round(v["grams"], 3)} for p, v in r["products"].items()]
        rows += [{"Species": f"{s} (excess left)", "Theoretical moles": round(v["left_moles"], 4),
                  "Theoretical mass (g)": round(v["left_grams"], 3)} for s, v in r["excess"].items()]
        st.dataframe(pd.DataFrame(rows), hide_index=True)
        first = next(iter(r["products"]))
        actual = st.number_input(f"Actual mass of {first} obtained (g)", value=0.0, min_value=0.0, key="f_actual")
        if actual > 0:
            st.metric("Percent yield", f"{percent_yield(actual, r['products'][first]['grams']):.1f} %")
