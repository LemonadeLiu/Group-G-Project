# Exploring materials with PySCF: Electronic Structures and Optical Properties 

## Final Project, Computational Physics, Fall 2024

_Yichen Guo_, _Xiaohui Liu_, _Xinyue Peng_


## Introduction 

Scientists have long sought to solve the quantum mechanical many-electron Schrödinger equation to elucidate electronic structures and predict the properties of molecules and materials.

Solving this equation for real systems is an extremely challenging task due to the complex interactions between electrons. The exact solution is computationally infeasible for systems with more than a few particles, making approximate methods like **Hartree-Fock Self-Consistent Field** (HF-SCF) and **Kohn-Sham Density Functional Theory** (KS-DFT) indispensable.

HF-SCF and KS-DFT are two foundational methods in quantum chemistry and condensed matter physics with fundamentally different approximations and frameworks.

### Hartree-Fock Self-Consistent Field (HF-SCF)
The HF-SCF method provides an approximate solution to the many-electron problem by representing the many-body wavefunction as a single Slater determinant constructed from molecular orbitals. This method simplifies the complex many-body interactions by treating the repulsion between electrons in an average sense, leading to a mean-field description of electron-electron interactions.

### Kohn-Sham Density Functional Theory (KS-DFT)
KS-DFT is a more contemporary and widely used approach, based on the Hohenberg-Kohn theorems and the Kohn-Sham formalism. Unlike HF-SCF, which explicitly treats wavefunctions, KS-DFT reformulates the many-electron problem in terms of the electron density, a simpler quantity.

## Installation

## Example: Raman for H2O...

## Example: DFT calculation of Graphene

## Example: Strain-induced effect in band structures and DOS

## Resources
