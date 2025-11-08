# Modeling Mucociliary Mixing and Transport at Tissue Scale

**Ling Xu¹, Pejman Sanaei², Yi Jiang²**  
¹Department of Mathematics and Statistics, North Carolina A&T State University, Greensboro, NC  
²Department of Mathematics and Statistics, Georgia State University, Atlanta, GA  

*To appear in PLoS Comp Bio, 2025*

## Overview

This repository contains MATLAB codes for 3D simulations of mucociliary clearance at tissue scale. The model uses the regularized stokeslet method to simulate cilia beating as rigid rods with prescribed motion, driving Stokes flow in mucus and tracking massless particle transport.

## Requirements

- MATLAB R2020b or later
- Parallel Computing Toolbox (optional, for faster computation)

## Installation

```bash
git clone https://github.com/Jiang-Lab/Mucociliary_Clearance-Tissue_Scale_model
cd Mucociliary_Clearance-Tissue_Scale_model
```

Add to MATLAB path:
```matlab
addpath(genpath('src'));
```

## Repository Structure

```
├── src/                    # Core simulation code
│   
├── scripts/                # Scripts and data for generating each of the figures 
│   
└── example/               # Example usage
```

## Quick Start

```matlab
% Single cilia cluster velocity field (Fig. 11)
scripts/codes_1cluster/plot_1cluster_1rod_grid_vel.m

% Three clusters with spacing sweep (Fig. 16)
scripts/figure16/plot_fig16_3cluster_varyD_nrods_mp_1_shift_center.m
```

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Rod length | 7 μm | Maximum cilium length |
| Beating frequency | 18 Hz | Cilium oscillation |
| Cilia spacing | 4 μm | Within cluster (D) |
| Cluster spacing | 3D-12D | Between clusters (Dc) |
| Ciliary density | 0.1, 0.2, 0.4 | Fraction of area (ν) |
| Domain size | 100×100 μm | Tissue patch |

## Main Findings

1. **Swirl size scales with ciliary density** - Higher density enhances transport
2. **Single cluster generates directional transport** - Horizontal and upward movement
3. **Optimal cluster spacing exists** - Maximum transport at Dc = 5D
4. **Metachronal waves hinder horizontal transport** - But enhance vertical mixing
5. **Spatially inhomogeneous diffusion** - Particles aggregate into discrete groups


## Citation

```bibtex
@article{xu2025mucociliary,
  title={Modeling mucociliary mixing and transport at tissue scale},
  author={Xu, Ling and Sanaei, Pejman and Jiang, Yi},
  year={2025},
  note={TBA}
}
```

## Contact

- Ling Xu: lxu@ncat.edu
- Yi Jiang: yjiang12@gsu.edu

## License

MIT LicenseMIT License

Copyright (c) 2025 Jiang Lab, Georgia State University


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)


