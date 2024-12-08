# Exploring materials with PySCF: Electronic Structures and Optical Properties 

## Final Project, Computational Physics, Fall 2024

_Yichen Guo_, _Xiaohui Liu_, _Xinyue Peng_


## Introduction 

Scientists have long sought to solve the quantum mechanical many-electron Schrödinger equation to elucidate electronic structures and predict the properties of molecules and materials.

Solving this equation for real systems is an extremely challenging task due to the complex interactions between electrons. The exact solution is computationally infeasible for systems with more than a few particles, making approximate methods like **Hartree-Fock Self-Consistent Field** (HF-SCF) and **Kohn-Sham Density Functional Theory** (KS-DFT) indispensable.

HF-SCF and KS-DFT are two foundational methods in quantum chemistry and condensed matter physics with fundamentally different approximations and frameworks.

### Hartree-Fock Self-Consistent Field (HF-SCF)
The HF-SCF method provides an approximate solution to the many-electron problem by representing the ground state wavefunction as a single Slater determinant $$\Psi_0$$ constructed from molecular orbitals $$\psi_0$$
```math
\Psi_0 = A|\psi_1(1)\psi_2(2)...\psi_N(N)|
```
Then the ground state energy could be represented as
```math
E = \braket{\Psi_0| \hat{H} | \Psi_0}
```
The minimization of the total energy leads to a self-consistent equation

$$\mathbf{FC} = \mathbf{SCE}$$

Here C is the matrix of molecular orbital coefficients, E is a diagonal matrix of the corresponding eigenenergies, and S is the atomic orbital overlap matrix. The Fock matrix F could be defined as 

$$\mathbf{F} = \mathbf{T}+\mathbf{V}+\mathbf{J}+\mathbf{K}$$

where $\mathbf{T}$ is the kinetic energy matrix, $\mathbf{V}$ is the external potential, $\mathbf{J}$ is the Coulomb matrix, and $\mathbf{K}$ is the exchange matrix.

HF-SCF simplifies the complex many-body interactions by treating the repulsion between electrons in an average sense, leading to a mean-field description of electron-electron interactions.

### Kohn-Sham Density Functional Theory (KS-DFT)
KS-DFT is a more contemporary and widely used approach based on the Hohenberg-Kohn theorems and the Kohn-Sham formalism. It uses the same formulation as HF-SCF and calculates the same self-consistent equation, but replaces the Fock matrix with a different Fock potential

```math
\mathbf{F} = \mathbf{T}+\mathbf{V}+\mathbf{J}+\mathbf{E_{xc}}
```

Where $\mathbf{E_{xc}}$ is the exchange-correlation(xc) energy and is approximated by a density functional approximation. Based on different structures and interactions in different systems, there're lots of different density functional approximation. For example, local density approximation (LDA) depends only on the electron density $\rho$ and generalized gradient approximation (GGA) depends on both electron density $\rho$ and the density gradient $|\nabla\rho|$.

## Installation

## Example: Raman for H2O...

## Example: DFT calculation of Graphene



## Example: Strain-induced effect in band structures and DOS ...

## Resources

[PySCF](https://pyscf.org/user/dft.html)
[ASE](https://wiki.fysik.dtu.dk/ase/index.html)
