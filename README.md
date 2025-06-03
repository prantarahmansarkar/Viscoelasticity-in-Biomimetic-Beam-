# 📘 Viscoelasticity in Biomimetic Beam

This repository contains the MATLAB source codes and raw datasets used for generating the results presented in the manuscript:

**"Viscoelasticity of Biomimetic Scale Beams from Trapped Complex Fluids"**

---

## 📂 MATLAB Scripts

- **`Contour_plots.m`**  
  Generates the contour plots shown in **Figure 7**. Associated raw data:
  - `Figs. 7(a)-(b).xlsx`
  - `Fig. 7(c).xlsx`
  - `Fig. 7(d).xlsx`

- **`Couette_flow_beam.m`**  
  Implements Eq. (5) for constant shear rate `ψ̇`.

- **`Couette_flow_beam_sinusoidal_curvature.m`**  
  Implements Eq. (5) for sinusoidal curvature input `ψ(t)`.

- **`Couette_flow_beam_variable_Omega_Omega_n.m`**  
  Simulates RED behavior for varying `Ω/Ωₙ`.

- **`Couette_flow_beam_variable_alphaL.m`**  
  Computes RED factor for different values of `αL`.

- **`Couette_flow_beam_variable_alphaL_delta_L.m`**  
  Used to generate **Figure 6**, examining combined variation of `αL` and `δL`.

- **`Couette_flow_beam_variable_deltaL.m`**  
  Calculates RED for varying `δL`.

- **`Couette_flow_beam_variable_eta_new.m`**  
  Computes RED results for varying scale overlap `η`.

---

## Excel Data Files

- `Fig. 2.xlsx`  
- `Fig. 3(a).xlsx`  
- `Fig. 3(b).xlsx`  
- `Fig. 3(c)-(d).xlsx`  
- `Figs. 4-5.xlsx`  

These files contain the raw numerical data used to plot Figures 2–5.

---

## Notes

- All figures and trends in the manuscript can be reproduced using the provided scripts and datasets.
- Scripts are commented for clarity and use standard MATLAB functions.
- No special toolboxes are required beyond core MATLAB functionality.

---

## Citation

If you use this code, please cite the corresponding manuscript once published.

