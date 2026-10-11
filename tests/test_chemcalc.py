import pytest
from chemcalc import *


def test_molar_mass():
    assert molar_mass("H2O") == pytest.approx(18.015, abs=0.01)
    assert molar_mass("Ca(OH)2") == pytest.approx(74.09, abs=0.01)
    assert molar_mass("CuSO4.5H2O") == pytest.approx(249.68, abs=0.05)


def test_parse_nested():
    assert parse_formula("Al2(SO4)3") == {"Al": 2, "S": 3, "O": 12}


def test_bad_formula():
    with pytest.raises(ValueError):
        parse_formula("Xx2")
    with pytest.raises(ValueError):
        parse_formula("Ca(OH2")


def test_balance():
    assert balance("Fe + O2 -> Fe2O3") == "4 Fe + 3 O2 -> 2 Fe2O3"
    assert balance("C3H8 + O2 -> CO2 + H2O") == "C3H8 + 5 O2 -> 3 CO2 + 4 H2O"


def test_gas():
    assert ideal_gas(p=1, n=1, t=273.15) == pytest.approx(22.414, abs=0.01)


def test_ph():
    assert ph_strong_acid(1e-3) == pytest.approx(3.0)
    assert ph_weak_acid(0.1, 1.8e-5) == pytest.approx(2.88, abs=0.01)
    assert ph_buffer(4.76, 0.1, 0.1) == pytest.approx(4.76)


def test_dilution():
    assert dilution(m1=2, v1=0.1, m2=0.5) == pytest.approx(0.4)
