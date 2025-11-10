# Spectra of m-Banded Toeplitz operators

**Authors:** E. O. HILTUNEN and Y. DE BRUIJN  
**Date:** 9.11.2025

---

In this computational notebook, we provide the MATLAB code for the computations in [1].

## II.1 The complex band structure
We plot the phase and the magnitude of the roots to the polynomial $f_m(z) - \lambda = 0$ as a function of real valued $\lambda$. The code `CBScontinuous.m` is implemented via a root tracking algorithm and plots the roots directly as $z(\lambda)$. Another method is to compute the roots sparately for each $\lambda$, this is also implented in `CBSscatter.m`.

- `CBScontinuous.m`
<p align="center"> <img src="CBSandDecay.png" alt="BandedCBS" width="400"/> </p>

## II.2 Open spectrum

This animation illustrates how the open spectrum comprises the interfection of the spectra of Topelitz operators evalueted on the $r$-scaled torus, i.e.

<p align="center"> <img src="OpenLimitDefinition.png" alt="OpenLimit" width="400"/> </p>
$\lim_{n \to\infty}\sigma\bigl(\mathbf{T}_n(f)\bigr) = \bigcap_{r \in (0, \infty)} \sigma\Bigl(\mathbf{T}\bigl(f(r\mathbb{T})\bigr)\Bigr)$

- `CollapseSymbolMovie.m`

<p align="center"> <img src="animation_r_variation.gif" width="500"/>  <p align="center"> <img src="animation_r_variation_complex.gif" width="500"/>  


