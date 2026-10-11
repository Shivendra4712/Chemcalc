import math
import numpy as np
import pytest
from chemcalc import (balance, solutions, gas, kinetics, thermo, equilibrium,
                      periodic, units, stoich, md)


def test_limiting_reagent():
    r = stoich.limiting_reagent("H2 + O2 -> H2O", {"H2": 4.0, "O2": 32.0})
    assert r["balanced"] == "2 H2 + O2 -> 2 H2O"
    assert r["limiting"] == "H2"
    assert r["products"]["H2O"]["grams"] == pytest.approx(4.0 / 2.016 * 18.015, rel=1e-3)
    assert stoich.percent_yield(30, 40) == pytest.approx(75)


def test_equilibrium_hi():
    rx, pr = equilibrium.parse_reaction("H2 + I2 <=> 2HI")
    assert rx == [(1, "H2"), (1, "I2")] and pr == [(2, "HI")]
    x, r, p = equilibrium.solve_equilibrium([(1, 1.0), (1, 1.0)], [(2, 0.0)], 50.0)
    assert x == pytest.approx(0.7795, abs=1e-3)
    assert p[0] == pytest.approx(1.559, abs=2e-3)
    q = equilibrium.reaction_quotient([(1, 1.0), (1, 1.0)], [(2, 0.5)])
    assert "right" in equilibrium.direction(q, 50.0)


def test_thermo():
    assert thermo.gibbs(-92.2, -198.7, 298.15) == pytest.approx(-32.96, abs=0.05)
    assert thermo.crossover_temperature(-92.2, -198.7) == pytest.approx(464.0, abs=1)
    k = thermo.k_from_dg(-10.0, 300)
    assert thermo.dg_from_k(k, 300) == pytest.approx(-10.0)
    assert thermo.van_t_hoff(1.0, 300, 0.0, 400) == pytest.approx(1.0)


def test_electrochem():
    cell = thermo.galvanic_cell("Cu2+/Cu", "Zn2+/Zn")
    assert cell["E0"] == pytest.approx(1.10) and cell["n"] == 2 and cell["spontaneous"]
    assert cell["dG_kJ"] == pytest.approx(-212.3, abs=0.5)
    assert thermo.nernst(1.10, 2, 1) == pytest.approx(1.10)
    assert thermo.nernst(1.10, 2, 10) == pytest.approx(1.0704, abs=1e-3)


def test_periodic():
    assert periodic.electron_configuration(26) == "1s2 2s2 2p6 3s2 3p6 3d6 4s2"
    assert periodic.electron_configuration(26, short=True) == "[Ar] 3d6 4s2"
    assert periodic.electron_configuration(24, short=True) == "[Ar] 3d5 4s1"
    assert periodic.electron_configuration(29, short=True) == "[Ar] 3d10 4s1"
    assert periodic.electron_configuration(46, short=True) == "[Kr] 4d10"
    assert periodic.get_element("fe").z == 26 and periodic.get_element(8).symbol == "O"
    assert len(periodic.ELEMENTS) == 54
    with pytest.raises(ValueError):
        periodic.get_element("Unobtainium")


def test_units():
    assert units.convert(1, "atm", "kPa", "Pressure") == pytest.approx(101.325)
    assert units.convert(100, "\u00b0C", "K", "Temperature") == pytest.approx(373.15)
    assert units.convert(32, "\u00b0F", "\u00b0C", "Temperature") == pytest.approx(0)
    assert units.convert(1, "\u00c5", "nm", "Length") == pytest.approx(0.1)


def test_more_ph():
    assert solutions.ph_strong_base(1e-3) == pytest.approx(11.0)
    assert solutions.ph_weak_base(0.1, 1.8e-5) == pytest.approx(11.12, abs=0.01)


def test_vdw():
    a, b = gas.VDW_CONSTANTS["CO2"]
    p = gas.van_der_waals_pressure(22.414, 1, 273.15, a, b)
    assert float(p) == pytest.approx(0.9946, abs=1e-3)
    assert float(gas.compressibility(p, 22.414, 1, 273.15)) < 1.0


def test_fit_order():
    t = np.arange(0, 60, 10.0)
    best, out = kinetics.fit_order(t, np.exp(-0.03 * t))
    assert best == 1 and out[1]["k"] == pytest.approx(0.03, rel=1e-6)
    best, out = kinetics.fit_order(t, 1 / (1 + 0.2 * t))
    assert best == 2 and out[2]["k"] == pytest.approx(0.2, rel=1e-6)


def test_md_energy_conservation():
    res = md.run_md(n_side=5, steps=400, thermostat=False, density=0.6)
    assert res.temperature[0] == pytest.approx(1.0, abs=1e-6)
    assert res.energy_drift_per_particle < 0.02
    assert len(res.snapshots) == 3 and res.rdf_g.shape == res.rdf_r.shape


def test_molecule():
    pytest.importorskip("rdkit")
    from chemcalc import molecule
    d = molecule.descriptors("CC(=O)Oc1ccccc1C(=O)O")
    assert d["formula"] == "C9H8O4" and d["mol_weight"] == pytest.approx(180.16, abs=0.05)
    assert d["hbd"] == 1 and d["lipinski_pass"]
    with pytest.raises(ValueError):
        molecule.descriptors("not-a-smiles(((")


def test_ideal_gas_arrays():
    v = np.array([10.0, 20.0])
    p = gas.ideal_gas(v=v, n=1, t=300)
    assert p.shape == (2,) and p[0] == pytest.approx(2 * p[1])
    with pytest.raises(ValueError):
        gas.ideal_gas(p=1, v=-1, n=1)
