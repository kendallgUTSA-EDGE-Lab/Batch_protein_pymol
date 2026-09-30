# Batch_protein_pymol
Pipeline to run batch proteins to visualize catalytic residues in pyMOL
# DUF PyMOL Visualization Pipeline

This repository contains a PyMOL batch visualization pipeline used to generate
structural figures for Drosophila proteins containing domains of unknown
function (DUFs).

## Visualization scheme

For each protein:

- Entire structure is displayed using spheres
- Protein outside the DUF is colored `gray50`
- DUF region is colored `tv_blue`
- Predicted catalytic residues are colored `tv_orange`
- Images are rendered at 2400 × 2400 pixels
- PNG images are saved with transparent backgrounds

## Input

The pipeline expects AlphaFold PDB files in:

```text
pdb_files/
