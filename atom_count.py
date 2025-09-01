from rdkit import Chem

mol = Chem.MolFromSmiles("NC(=O)[C@@H]1CC[C@@H](N1)c1ccc(cc1)OCc1ccccc1F.Cl")

num_all_atoms = mol.GetNumAtoms(onlyExplicit=False)
print(f"Number of all atoms: {num_all_atoms}")
