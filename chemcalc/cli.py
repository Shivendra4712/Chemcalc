"""Command line + interactive menu for ChemCalc."""
import argparse
import sys
from .elements import ATOMIC_MASSES
from . import (molar_mass, percent_composition, balance, dilution,
               ph_strong_acid, ph_weak_acid, ph_buffer, ideal_gas, molarity)


def _num(prompt):
    s = input(prompt).strip()
    return float(s) if s else None


def menu():
    options = """
=== ChemCalc ===
1. Molar mass            5. pH (strong acid)
2. Percent composition   6. pH (weak acid)
3. Balance equation      7. Buffer pH
4. Dilution (M1V1=M2V2) 8. Ideal gas law
9. Molarity              0. Exit"""
    while True:
        print(options)
        ch = input("Choice: ").strip()
        try:
            if ch == "1":
                f = input("Formula: ")
                print(f"Molar mass = {molar_mass(f):.3f} g/mol")
            elif ch == "2":
                for el, p in percent_composition(input("Formula: ")).items():
                    print(f"  {el}: {p:.2f}%")
            elif ch == "3":
                print(balance(input("Equation (e.g. H2 + O2 -> H2O): ")))
            elif ch == "4":
                print("Leave ONE value blank to solve for it.")
                print("Result =", dilution(_num("M1: "), _num("V1: "),
                                           _num("M2: "), _num("V2: ")))
            elif ch == "5":
                print(f"pH = {ph_strong_acid(_num('Conc (M): '), int(_num('Protons (1/2): ') or 1)):.3f}")
            elif ch == "6":
                print(f"pH = {ph_weak_acid(_num('Conc (M): '), _num('Ka: ')):.3f}")
            elif ch == "7":
                print(f"pH = {ph_buffer(_num('pKa: '), _num('[A-]: '), _num('[HA]: ')):.3f}")
            elif ch == "8":
                print("Units: atm, L, mol, K. Leave ONE blank.")
                print("Result =", ideal_gas(_num("P: "), _num("V: "),
                                            _num("n: "), _num("T: ")))
            elif ch == "9":
                print(f"M = {molarity(_num('Moles: '), _num('Volume (L): ')):.4f} mol/L")
            elif ch == "0":
                return
            else:
                print("Invalid choice.")
        except (ValueError, TypeError, ZeroDivisionError) as e:
            print("Error:", e)


def main(argv=None):
    p = argparse.ArgumentParser(prog="chemcalc", description="Computational chemistry calculator")
    sub = p.add_subparsers(dest="cmd")
    sub.add_parser("menu", help="interactive menu")
    sub.add_parser("mass").add_argument("formula")
    sub.add_parser("composition").add_argument("formula")
    sub.add_parser("balance").add_argument("equation")
    sub.add_parser("element").add_argument("symbol")
    g = sub.add_parser("gas")
    for k in "pvnt":
        g.add_argument(f"-{k}", type=float)
    d = sub.add_parser("dilute")
    for k in ("m1", "v1", "m2", "v2"):
        d.add_argument(f"--{k}", type=float)
    a = p.parse_args(argv)
    try:
        if a.cmd in (None, "menu"):
            menu()
        elif a.cmd == "mass":
            print(f"{molar_mass(a.formula):.3f} g/mol")
        elif a.cmd == "composition":
            for el, pc in percent_composition(a.formula).items():
                print(f"{el}: {pc:.2f}%")
        elif a.cmd == "balance":
            print(balance(a.equation))
        elif a.cmd == "element":
            from .periodic import get_element, electron_configuration
            e = get_element(a.symbol)
            print(f"{e.name} ({e.symbol})  Z={e.z}  mass={ATOMIC_MASSES[e.symbol]}  "
                  f"EN={e.electronegativity}  group {e.group}, period {e.period}  [{e.category}]")
            print("Configuration:", electron_configuration(e.z, short=True))
        elif a.cmd == "gas":
            print(ideal_gas(a.p, a.v, a.n, a.t))
        elif a.cmd == "dilute":
            print(dilution(a.m1, a.v1, a.m2, a.v2))
    except ValueError as e:
        sys.exit(f"Error: {e}")


if __name__ == "__main__":
    main()
