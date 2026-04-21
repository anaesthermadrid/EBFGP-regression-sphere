# Time-adaptive Functional Gaussian Process Regression

MATLAB research code accompanying the arXiv preprint:

**Ruiz-Medina, M.D., Madrid, A.E., Torres-Signes, A., Angulo, J.M. (2026)**  
*Time-adaptive functional Gaussian Process regression*  
arXiv:2603.21144

---

## Description

This repository contains MATLAB implementations of the **time-adaptive Empirical Bayes
Functional Gaussian Process (EBFGP) regression in the special case of compact manifolds 
defined by the sphere n the three dimensional Real line.

The proposed methodological approach is developed in the arXiv preprint.

These  two codes provide the implementation of EBFGP in the simulations
for two subfamilies (Cauchy and Matérn) in the spatiotemporal Gneiting class restricted to the sphere.

---

## Important note

⚠️ **This is research code**, provided primarily for transparency and reproducibility purposes.

- The current version corresponds to the **original simulation and inference scripts**
  used to generate the numerical results and figures presented in the manuscript.
- The code structure closely follows the implementation strategy described in the paper.
- At this stage, the code is **not modularized** and **not intended as a general-purpose
  software package**.
- A cleaned, modular (covering a wider range of truncation and spatial sampling 
  frequency parameter values, as well as functional sample sizes),
  and reusable version of the implementation may be released in the near future.

---

## Repository contents

- `EBFGP_Cauchy_Sphere.m` – Cauchy (LRD) covariance model implementation.
- `EBFGP_Matern_Sphere.m` – Matérn (SRD) covariance model implementation.

Each script is self-contained and can be executed independently.

---

## Requirements

- MATLAB R2020b or later
- Statistics and Machine Learning Toolbox

---

## Running the code

```matlab
EBFGP_Cauchy_Sphere
```

or

```matlab
EBFGP_Matern_Sphere
```

---

## Reproducibility remarks

Numerical results depend on Monte Carlo sampling. For exact reproducibility,
set the MATLAB random seed explicitly.

---

## Citation

If you use this code, please cite:

```bibtex
@article{RuizMedina2026EBFGP,
  title   = {Time-adaptive functional Gaussian Process regression},
  author  = {Ruiz-Medina, M.D. and Madrid, A.E. and Torres-Signes, A. and Angulo, J.M.},
  journal = {arXiv preprint arXiv:2603.21144},
  year    = {2026}
}
```

---

## License

This code is released under the MIT License.
