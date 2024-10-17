# Nonlinear Mechanical System Analysis: From Identification to Simulation

## Description
This project implements the analysis of a nonlinear mechanical system with two degrees of freedom (DOFs), characterized by nonlinearities in the stiffness and damping. The analysis covers the detection, characterization, and simulation of the system’s nonlinear behavior using various algorithms. Numerical methods like the **Shooting Algorithm** and **Sequential Continuation** are used to compute the **Frequency Response Curves (FRCs)**, and the **Nonlinear Normal Modes (NNMs)** are computed to explore the system's dynamic responses.

The advisor of the project was professor Gaëtan Kerschen from University of Liège while the co-advisor was Dr. Ghislain Raze.

## Requirements
- **MATLAB**
- **NI2D Software** for validation of numerical methods
- **Wavelet Toolbox** (for harmonic analysis)

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/alessandropase/NonLin_Mech_analysis.git
   ```
2. Ensure you have MATLAB installed along with the required toolboxes.


## Project Structure
- `Identification.m` — Identification part: Testing campaign, Characterization and Estimation of the nonlinearities
- `Simulation.m` — Simulation part: FRCs, NNMs and Basins of attraction computation
- `Functions` — MATLAB Functions: Shooting algorith, Continuation algorithm, Basins of attraction computation …
- `Datas` — Contains all the data from tests and simulations useful for plots
- `README.md` — Project documentation

## External Resources
- [NI2D Software](https://nolisys.com/?page_id=97) — Software used for validation.
