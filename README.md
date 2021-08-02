# Pattern Formation
- Explorations of pattern formation in partial differential equations
- May. 2, 2021 ~ Present

----------
## Diffusive prey-predator systems with Holling II functional response
<sup>∂<sup>&eta;</sup>x</sup>/<sub>∂t</sub> = &Delta; u  + u(1-u) - <sup>uv</sup>/<sub>u+&alpha;</sub>

### 1. Temporal **first**-derivative system
Forward Euler method in time and Central differences in space (FTCS) | [Code](https://github.com/OH-Seoyoung/Pattern_Formation/blob/master/Diffusive_prey-predator_systems/FTCS_first-derivative_Holling_II_functional_response.m)  

### 2. Temporal **fractional**-derivative system
- It can form steadily spatial patterns even though its first-derivative counterpart can't exhibit any steady pattern.
- Letnikov method in time (Caputo's derivative by Grunwald) and Central differences in space
- Zero-flux boundary condition
  
## References
```
[1] Yin, Hongwei, and Xiaoqing Wen. "Pattern formation through temporal fractional derivatives." Scientific reports 8.1 (2018): 1-9.
[2] Ciesielski, Mariusz, and Jacek Leszczynski. "Numerical simulations of anomalous diffusion." arXiv preprint math-ph/0309007 (2003).
```
