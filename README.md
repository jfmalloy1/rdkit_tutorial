# RDKit Tutorial

Basic RDKit tutorial, focusing on creating mol objects, converting from mol to strings, and visualizing chemicals.

[Official RDKit Tutorial](https://www.rdkit.org/docs/GettingStartedInPython.html)
[

## RDKit Setup

1. *Via pip (check system first!):* https://pypi.org/project/rdkit-pypi/
2. *Via conda:* https://www.rdkit.org/docs/Install.html

## 1: Creating Mol Objects

- Create a mol object
- Inchi to smiles & back
- Add/remove Hs

## 2: Describe Mol Objects
- Use [Descriptors](https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html) to give various properties 
- Identify specific atoms & bondds

## 3: Substructure Searche
- Finding specific substructures using [SMARTS](https://www.ics.uci.edu/~dock/manuals/DaylightTheoryManual/theory.smarts.html)
- Maximal Common Substructure

## 4: Fingerprints
- Building Fingerprints 
- Finding Tanimoto similarity

## 5: Chirality
- Chiral embedding & output over KEGG

## 6: Visualization
- Basic drawing
- Highlighting patterns
- Mol to picture and back "metadata in molecule images"
