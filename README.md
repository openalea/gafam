# GAFAM: Growing AgroForestry with Apple in Mediterranean climate

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14253434.svg)](https://doi.org/10.5281/zenodo.14253434) [![NBViewer](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.org/github/openalea/gafam/blob/paper/notebook/Paper.ipynb)

## Installation

conda install -c conda-forge -c openalea3 openalea.mtg pandas notebook ipython

## Analysis

Here is the variables that can be extracted
![Multiscale appeltree extraction](./extraction.jpg)

### Tree (per year)
- Sum( leaf area )
- Sum( GUs Lengths )
- Count( number of GUs )
- Count( FL, VG, LB, EX )
- Mean( angle of +B)
- Mean( nb_leaves )
- Mean( nb_fruits )
- Mean( nb_flowers )
- Mean( fruit_drops )

### Trunk/ Main axis (per year) (symbol A)
- Sum( A Length )
- Count( number of branches “+B”)
- Mean( A diameter )

### Branch (per year)
- Sum( leaf area )
- Sum( GU Lengths )
- Count ( GUs )
- Count( FL, VG, LB, EX )
- Angle of the branch
- Height of the branching internode (or its identity/node number)

### GU (per year and branching order: A, B, C, D)
- Sum( Leaf area )
- Sum( Length )
- Count( number of GUs )
- Count( FL, VG, LB, EX )


## Notebook for Reproducibility

[![NBViewer](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.org/github/openalea/gafam/blob/paper/notebook/Paper.ipynb)
