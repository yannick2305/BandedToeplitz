<h1 align="center">Spectra of <i>m</i>-Banded Toeplitz Operators</h1>
<p align="center">
  <b>Y. DE BRUIJN</b> and <b>E. O. HILTUNEN</b><br>
  <sub><i>University of Oslo</i></sub><br>
  <sub>November 9, 2025</sub>
</p>
<p align="center">
  <img src="https://img.shields.io/badge/MATLAB-R2023b-orange" alt="MATLAB">
</p>


**Abstract:** We provide the complete computational framework supporting the theoretical results in [1].


## II.1 The Complex Band Structure

We plot the phase and magnitude of the roots of f<sub>m</sub>(z) − λ = 0 using two approaches:

**Root Tracking Algorithm** (`CBScontinuous.m`)  
Continuously tracks roots as λ varies, plotting z(λ) directly. Blue crosses overlay the complex band structure, denoting numerically computed exponential decay rates of Toeplitz matrix eigenvectors.

<p align="center"> 
  <img src="Figures/CBSandDecay.png" alt="Complex Band Structure and Decay" width="350"/> 
</p>

**Discrete Root Finding** (`CBSscatter.m`)  
Computes roots independently for each λ value.

 <p align="center"> 
  <img src="Figures/CBSscatter.png" alt="Complex Band Structure and Decay" width="350"/> 
</p>


## II.2 Convergence of Floquet Parameters

The exponential decay length of eigenvectors for non-Hermitian *m*-banded Toeplitz operators depends on the bandwidth. We numerically demonstrate that decay length approaches zero as bandwidth increases.

**Convergence of Floquet Parameters** (`ConvergenceFloquetParameter.m`)

<p align="center"> 
  <img src="Figures/Screenshot 2025-11-11 at 14.17.17.png" alt="Convergence of Floquet Parameter" width="700"/> 
</p>


## II.3 Convergence of Pseudospectra

We demonstrate the convergence behavior of pseudospectra for finite truncations of *m*-banded Toeplitz operators.

**Pseudospectrum** (`PseudospectrumConvergence.m`)

<p align="center"> 
  <img src="Figures/PseudospectrumConvergence.png" alt="Pseudospectrum Convergence" width="700"/> 
</p>

## II.4 Spectrum of Open Limit

This animation illustrates how the open spectrum comprises the intersection of spectra of Toeplitz operators evaluated on the *r*-scaled torus:

$$\lim_{n\to\infty} \sigma\left(\mathbf{T}_n(f_m)\right) = \bigcap_{r>0} \sigma\left(\mathbf{T}\left(f_m(r\mathbb{T})\right)\right) = \left\lbrace \lambda \in \mathbb{C} ~:~ |z_{m}(\lambda)| = |z_{m+1}(\lambda)| \right\rbrace$$

**Collapsed symbol Movie** (`CollapseSymbolMovie.m`)

<p align="center"> <img src="Figures/animation_r_variation.gif" alt="Open spectrum collapse (real)" width="480"/> <br> <em>Figure 1: Open spectrum collapse (real)</em> </p> <p align="center"> <img src="Figures/animation_r_variation_complex.gif" alt="Open spectrum collapse (complex)" width="480"/> <br> <em>Figure 2: Open spectrum collapse (complex)</em> </p>


<p align="center"> 
  <img src="Figures/animation_r_variation.gif" alt="Open spectrum collapse (real)" width="450"/>
  <img src="Figures/animation_r_variation_complex.gif" alt="Open spectrum collapse (complex)" width="450"/>
</p>

## II.5 Reality of the Open Limit

We verify numerically that the open limit produces real-valued spectra providet that $\Lambda(f_m)$ is traced out by a polar curve.

**Set $\Lambda(f)$** (`OpenLimit.m`)

<p align="center"> 
  <img src="Figures/LambdaOfF.png" alt="Real-valued open limit" width="700"/> 
</p>


## II.6 Defect Modes

We numerically illustrate composite decay bounds acting on defect modes, demonstrating localization phenomena in perturbed systems.

**Eigenvector Decay Estimates** (`JaffardCBSEstimate.m`)

<p align="center"> 
  <img src="Figures/JaffardCBS.png" alt="Defect mode decay bounds" width="700"/> 
</p>

## II.7 Complex-Valued Frequencies

For the open spectrum of pristine Toeplitz operators—and more generally for defect modes—eigenvalues are no longer restricted to the real line. We extend the analysis from `CBScontinuous.m` to complex-valued frequencies.

**Method:**  
For each λ, we compute roots of *f*<sub>*m*</sub>(*z*) − λ = 0 and sort them in ascending order by magnitude: |*z*₁| ≤ ⋯ ≤ |*z*<sub>*m*</sub>| ≤ |*z*<sub>*m*+1</sub>| ≤ ⋯ ≤ |*z*₂<sub>*m*</sub>|. We then plot the decay parameter β where *e*<sup>−β</sup> = |*z*<sub>*m*+1</sub>|. For frequencies λ ∈ σ<sub>wind</sub>, we have |*z*<sub>*m*+1</sub>| < 1, corresponding to the region where β > 0.

**Complex Band Structure** (`CBScomplexLambda.m`)

<p align="center"> 
  <img src="Figures/CBSgeneralLambda.png" alt="Complex band structure for general frequencies" width="700"/> 
</p>


## II.8 Hermitian Matrices

For symbol functions with symmetric coefficients, the symbol can be expressed using trigonometric functions:

$$f_m(e^{i(\alpha + i \beta)}) = a_0 + 2 \sum_{k = 1}^m a_k\bigl(\cos(\alpha k)\cosh(\beta k) - 2i\sin(\alpha k)\sinh(\beta k) \bigr)$$

The imaginary part vanishes along specific paths in the complex plane. Consequently, allowed quasimomenta must be restricted to these contours.

**$\alpha$ and $\beta$ path** (`RealSymbolContour.m`)

<p align="center"> 
  <img src="Figures/RealSymbolFunctionHermitian.png" alt="Real symbol function contours for Hermitian case" width="700"/> 
</p>


## III. Non-Hermitian Skin Effect

Code for simulating the non-Hermitian skin effect in 3-dimensional systems is available at:  
**Repository:** https://github.com/jinghaocao/skin_effect

See also the accompanying paper [2] for theoretical background and additional numerical methods.


## IV. References

> [1] Davies, B., De Bruijn, Y., Dupuy, S. and Hiltunen, E.O. (2025), *TODO*

> [2] Habib Ammari, Silvio Barandun, Jinghao Cao, Bryn Davies, Erik Orvehed Hiltunen, Ping Liu, *The non-Hermitian skin effect with three-dimensional long-range coupling.* J. Eur. Math. Soc. (2025), https://ems.press/journals/jems/articles/14299016

## Citation

If you use this code in your research, please cite:

> Davies, B., De Bruijn, Y., Dupuy S. and Hiltunen, E.O. (2025), *TODO*

