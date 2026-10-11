import streamlit as st
from chemcalc import molecule as mo
from . import common as c


def render():
    c.header("\U0001f9ec", "Molecule Lab", "SMILES \u2192 structure, descriptors and Lipinski drug-likeness")
    if not mo.HAS_RDKIT:
        st.warning("RDKit is not installed. Run `pip install rdkit` and restart the app.")
        return
    ex = st.selectbox("Example molecule", list(mo.EXAMPLES), key="mo_ex")
    smi = st.text_input("SMILES", mo.EXAMPLES[ex], key=f"mo_smi_{ex}")
    try:
        d = mo.descriptors(smi)
        img = mo.draw(smi)
    except ValueError as e:
        st.error(str(e)); return
    a, b = st.columns([1, 1])
    a.image(img, caption=d["formula"])
    m = b.columns(2)
    m[0].metric("Molecular weight", f"{d['mol_weight']:.2f}"); m[1].metric("logP", f"{d['logp']:.2f}")
    m = b.columns(2)
    m[0].metric("H-bond donors", d["hbd"]); m[1].metric("H-bond acceptors", d["hba"])
    m = b.columns(2)
    m[0].metric("TPSA (\u00c5\u00b2)", f"{d['tpsa']:.1f}"); m[1].metric("Rotatable bonds", d["rotatable_bonds"])
    st.markdown("#### Lipinski's Rule of Five")
    if d["lipinski_pass"]:
        st.success("Drug-like: at most one violation." if d["violations"] else "Drug-like: no violations.")
    else:
        st.error("Likely poor oral bioavailability.")
    for v in d["violations"]:
        st.write("\u2022 " + v)
    st.caption("Rule of Five: MW \u2264 500, logP \u2264 5, H-bond donors \u2264 5, acceptors \u2264 10.")
