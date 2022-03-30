import rdkit.Chem as Chem

### 1: Creating Mol Objects ###
print("\n----- Creating Mol Objects -----", end="\n\n")

inchi_adenine = "InChI=1S/C5H5N5/c6-4-3-5(9-1-7-3)10-2-8-4/h1-2H,(H3,6,7,8,9,10)"
canonical_smiles_adenine = "C1=NC2=NC=NC(=C2N1)N"

mol = Chem.MolFromInchi(inchi_adenine)
mol2 = Chem.MolFromSmiles(canonical_smiles_adenine)

#Can output molecules as inchi/smiles strings...
print("Adenine inchi:", Chem.MolToInchi(mol2))
print("Adenine smiles:", Chem.MolToSmiles(mol), end="\n\n")

#Can also output as mol blocks...
print("Adenine mol:", Chem.MolToMolBlock(mol), end="\n\n")

#Hs arn't implicitly added - must be implictly added
print("Number of atoms:", mol.GetNumAtoms())
mol3 = Chem.AddHs(mol)
print("Number of atoms with Hs:", mol3.GetNumAtoms(), end="\n\n")

print("Adenine mol with Hs:", Chem.MolToMolBlock(mol3), end="\n\n")

# ### 2: Describing Mol Objects ###
# import rdkit.Chem.rdMolDescriptors as Desc

# print("\n----- Describing Mol Objects -----", end="\n\n")

# #Can get a lot of (standardized) information about a molecule, such as...

# #...Molecular weight, with or without Hs (The Desc library includes Hs by default)
# print("Molecular Weight:", Desc.CalcExactMolWt(mol))
# print("Molecular Weight (no Hs):", Desc.CalcExactMolWt(mol, onlyHeavy=True))

# #...chemical formula
# print("Formula:", Desc.CalcMolFormula(mol))

# #...Sterocenters (more on this later, I don't trust this method all the time)
# print("Chiral Centers:", Desc.CalcNumAtomStereoCenters(mol))

# #We also can get structural information about the molecule, such as individual atoms, bonds, or rings
# atoms = mol3.GetAtoms()
# for atom in atoms:  #Atomic information
#     print(atom.GetAtomicNum(), end=" ")
# print()
# for bond in mol3.GetBonds():  #Bond information
#     print(bond.GetBondType(), end=" ")
# print()

# #Detailed information about a specific atomc
# atom = mol3.GetAtomWithIdx(0)
# print("Selected atom:", atom.GetAtomicNum())
# print("Neighbors:", end=" ")
# for neighbor in atom.GetNeighbors():
#     print(neighbor.GetAtomicNum(), end=" ")
# print()
# print("Bond between atom & neighbor 0:",
#       atom.GetNeighbors()[0].GetBonds()[0].GetBondType())

# #Ring information
# print("Number of rings:", mol3.GetRingInfo().NumRings())

# ### 3. Substructure Search
# print("\n----- Substructures -----", end="\n\n")

# adenine_mol = mol3
# guanine_mol = Chem.MolFromInchi(
#     "InChI=1S/C5H5N5O/c6-5-9-3-2(4(11)10-5)7-1-8-3/h1H,(H4,6,7,8,9,10,11)")
# guanine_mol = Chem.AddHs(guanine_mol)

# #Can find specific patterns within a molecule, using SMARTS
# patt = Chem.MolFromSmarts("c-N")
# print("Is there a C-N bond in Adenine?", adenine_mol.HasSubstructMatch(patt))
# print("Which atoms are involved?",
#       adenine_mol.GetSubstructMatch(patt),
#       end="\n\n")

# #Note: there is a `useChirality` which can be useful for chiral molecules

# # Can find substructures (Maximal Common Substructure, to be precise) between two mol objects
# from rdkit.Chem import rdFMCS

# print("MCS between adenine & guanine:",
#       rdFMCS.FindMCS([adenine_mol, guanine_mol]).smartsString)

# ### 4. Fingerprints & Similarity
# print("\n----- Fingerprints -----", end="\n\n")

# #Fingerprints are unique identifers of a compound, usually in vector form
# #RDKit has a native fingerprint structure, there are many ways to do this
# adenine_fp = Chem.RDKFingerprint(adenine_mol)
# guanine_fp = Chem.RDKFingerprint(guanine_mol)

# from rdkit.DataStructs.cDataStructs import ExplicitBitVect

# print("Adenine fingerprint:",
#       ExplicitBitVect.ToBitString(adenine_fp),
#       end="\n\n")
# print("Guanine fingerprint:",
#       ExplicitBitVect.ToBitString(guanine_fp),
#       end="\n\n")

# #These fingerprints can be compared to each other (there's many ways to do this, Tanimoto is Euclidian distance)
# from rdkit import DataStructs

# print("Tanimoto Similarity between Adenine & Guanine:",
#       DataStructs.FingerprintSimilarity(adenine_fp, guanine_fp))
# """ 4.5: Coding challenge time!! 

# Find two molecules of your choice, and:
# a) Find the molecular weight & formula
# b) Find the Maximal Common Substructure between them
# c) Find the Tanimoto Similarity between them

# """
# ######
# #Code Goes here
# #####

# ### 5. Chiral Embedding (over LUCA compounds)

# print("\n----- Chiral Analysis -----", end="\n\n")

# #Read in LUCA compounds, find smiles strings associated with each compound
# import json
# import pandas as pd

# fpath = "Luca_ec_compound_unique_list_clean.json"
# with open(fpath) as f:
#     data = json.load(f)
# luca_cpds = data["Luca"]

# kegg_df = pd.read_csv("chiral_molweight_formula_labels.csv")
# luca_smiles = kegg_df[kegg_df["C"].isin(luca_cpds)]["S"].tolist()

# print("Subset of LUCA:", luca_smiles[0:10])

# from rdkit.Chem import AllChem
# import re
# from tqdm import tqdm

# luca_cpds_chiral = []
# print("\n--- Analyzing LUCA compounds ---")
# for cpd in tqdm(luca_smiles):
#     try:
#         #Find canonical smiles
#         m = Chem.MolFromSmiles(cpd)
#         Chem.SanitizeMol(m)  #Cleans 3D structures
#         m2 = Chem.AddHs(m)  #Should be added in order to form proper geometry
#         AllChem.EmbedMolecule(m2)  #Embeds molecule in 3D space
#         Chem.AssignAtomChiralTagsFromStructure(m2)  #Assigns chiral tags
#         chiral_smiles = Chem.MolToSmiles(
#             m2, isomericSmiles=True
#         )  #Ensures chiral tags are included in smiles output

#         #Replace @@ (L) for @
#         chiral_smiles = re.sub("@@", '@', chiral_smiles)

#         luca_cpds_chiral.append(chiral_smiles)
#     except:
#         pass

# print("\nChiral subset of LUCA:", luca_cpds_chiral[0:10], end="\n\n")

# #Get chiral percentage - this is the percentage of smiles strings with a '@'
# luca_chiral_count = 0
# for cpd in luca_cpds_chiral:
#     if "@" in cpd:
#         luca_chiral_count += 1

# print("LUCA chiral percentage:",
#       luca_chiral_count / float(len(luca_cpds_chiral)))

# ### 6. Drawing
