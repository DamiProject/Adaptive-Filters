# Adaptive-Filters

## Objective

**Adaptive Filtering for Stochastic System Identification:** This project implements and compares adaptive filtering algorithms for identifying unknown systems under white-noise and coloured-noise input conditions. The simulations evaluate convergence behavior, tracking behavior, ensemble-averaged mean-squared error, and adaptive weight trajectories for both time-invariant and time-varying systems. It also shows the importance of selecting appropriate step-size and forgetting-factor parameters for LMS-based and RLS-based adaptive algorithms. These parameters affect convergence speed, tracking ability, steady-state error, and computational complexity tradeoffs.

## Project Status
This repository is under active development. The core adaptive filter algorithms are implemented, while future improvements may include code refactoring, unit testing, additional stochastic test scenarios, fixed-point implementation, and hardware-oriented optimization.

## Proof of Concept 

**Simulation Structure:** The adaptive filters used in this project are Normalized Least Mean Square (NLMS), Recursive Least Square (RLS), and Discrete Fourier Transform Least Mean Square (DFT-LMS). These adaptive filters follow the same basic system identification workflow, as shown in Figure 1. However, DFT-LMS includes a preprocessing block that transforms the input signal before adaptive filtering. This helps reduce input correlation and improves convergence behavior, especially for coloured noise inputs, as shown in Figure 2.

<img width="400" height="300" alt="System identification drawio" src="https://github.com/user-attachments/assets/a52f8cbe-6ccb-461f-9bac-846532129c92" /> 

Figure 1: System identification block diagram using adaptive filters.

<img width="400" height="300" alt="DFT system identification drawio" src="https://github.com/user-attachments/assets/196c80e8-ccd6-40d2-a489-e7463e33c45f" />

Figure 2: System identification block diagram using DFT-LMS.

**Simulation:** 

1. **Step-size Impact on Tradeoffs:** The DFT-LMS adaptive filter was used to estimate the impulse response of a time-varying unknown system modeled as a first-order Markov process. Three different step sizes were evaluated to demonstrate how the step size affects convergence speed, tracking ability, and steady-state error, as shown in the table and figures below.

<img width="1037" height="611" alt="image" src="https://github.com/user-attachments/assets/fad5c564-ed29-4e3b-9fb8-427c68f519f1" />

Figure 3: 

<img width="715" height="393" alt="image" src="https://github.com/user-attachments/assets/66245002-c287-4258-8418-f85800c81986" />










