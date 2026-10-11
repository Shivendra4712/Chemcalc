import numpy as np
import pytest
from chemcalc import kinetics, molecular, titration


def test_first_order_half_life():
    k, a0 = 0.05, 1.0
    t_half = kinetics.half_life(1, a0, k)
    assert kinetics.concentration(1, a0, k, t_half) == pytest.approx(a0 / 2)


def test_second_order_half_life():
    t_half = kinetics.half_life(2, 0.5, 0.2)
    assert kinetics.concentration(2, 0.5, 0.2, t_half) == pytest.approx(0.25)


def test_arrhenius_roundtrip():
    k1 = kinetics.arrhenius(1e10, 50000, 300)
    k2 = kinetics.arrhenius(1e10, 50000, 320)
    assert kinetics.activation_energy(k1, 300, k2, 320) == pytest.approx(50000)


def test_consecutive_mass_balance():
    t = np.linspace(0, 50, 100)
    a, b, c = kinetics.consecutive(1.0, 0.2, 0.1, t)
    assert np.allclose(a + b + c, 1.0)
    assert np.allclose(kinetics.consecutive(1, 0.1, 0.1, t)[1],
                       1 * 0.1 * t * np.exp(-0.1 * t))


def test_lennard_jones():
    eps, sig = 0.996, 3.405
    rmin = molecular.lj_minimum(sig)
    assert molecular.lennard_jones(rmin, eps, sig) == pytest.approx(-eps)
    assert molecular.lennard_jones(sig, eps, sig) == pytest.approx(0, abs=1e-12)
    assert molecular.lj_force(rmin, eps, sig) == pytest.approx(0, abs=1e-9)


def test_mb_speeds_n2():
    s = molecular.mb_speeds(28.014, 300)
    assert s["rms"] == pytest.approx(517, abs=2)
    assert s["most_probable"] < s["mean"] < s["rms"]


def test_mb_normalised():
    v = np.linspace(0, 5000, 20001)
    area = np.trapz(molecular.mb_distribution(v, 28.014, 300), v) if hasattr(np, "trapz") \
        else np.trapezoid(molecular.mb_distribution(v, 28.014, 300), v)
    assert area == pytest.approx(1.0, abs=1e-3)


def test_titration_strong():
    assert titration.titration_ph(0.1, 25, 0.1, 25) == pytest.approx(7.0, abs=0.01)
    assert titration.titration_ph(0.1, 25, 0.1, 0) == pytest.approx(1.0, abs=0.01)


def test_titration_weak():
    ka = 1.8e-5
    assert titration.titration_ph(0.1, 25, 0.1, 0, ka) == pytest.approx(2.88, abs=0.02)
    assert titration.titration_ph(0.1, 25, 0.1, 12.5, ka) == pytest.approx(4.745, abs=0.02)
    assert titration.titration_ph(0.1, 25, 0.1, 25, ka) == pytest.approx(8.72, abs=0.03)
