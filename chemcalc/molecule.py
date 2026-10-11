"""Molecule properties & Lipinski drug-likeness from SMILES (needs RDKit)."""
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Draw
    HAS_RDKIT = True
except ImportError:  # pragma: no cover
    HAS_RDKIT = False

EXAMPLES = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "Ethanol": "CCO",
    "Benzene": "c1ccccc1",
    "Glucose": "OCC1OC(O)C(O)C(O)C1O",
}


def _mol(smiles: str):
    if not HAS_RDKIT:
        raise ImportError("RDKit is not installed. Run: pip install rdkit")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ValueError("Invalid SMILES string")
    return mol


def descriptors(smiles: str) -> dict:
    mol = _mol(smiles)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    violations = []
    if mw > 500: violations.append("Molecular weight > 500")
    if logp > 5: violations.append("logP > 5")
    if hbd > 5: violations.append("H-bond donors > 5")
    if hba > 10: violations.append("H-bond acceptors > 10")
    return {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mol_weight": mw, "logp": logp, "hbd": hbd, "hba": hba,
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "violations": violations,
        "lipinski_pass": len(violations) <= 1,
    }


def draw(smiles: str, size=(420, 320)):
    return Draw.MolToImage(_mol(smiles), size=size)
