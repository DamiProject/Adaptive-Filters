# Adaptive-Filters

## Objective

**Adaptive Filtering for Stochastic System Identification:** This project implements and compares adaptive filtering algorithms for identifying unknown systems under white-noise and coloured-noise input conditions. The simulations evaluate convergence behavior, tracking behavior, ensemble-averaged mean-squared error, and adaptive weight trajectories for both time-invariant and time-varying systems. It also shows the importance of selecting appropriate step-size and forgetting-factor parameters for LMS-based and RLS-based adaptive algorithms. These parameters affect convergence speed, tracking ability, steady-state error, and computational complexity tradeoffs.

## Project Status
This repository presents a proof of concept MATLAB implementation of adaptive filtering algorithms for stochastic system identification. Current simulations demonstrate convergence, tracking, MSE behavior, and algorithm tradeoffs. Future work includes additional unit tests and stochastic parameter tests, fixed-point implementation, and hardware-oriented optimization.

## Proof of Concept 

**Simulation Structure:** The adaptive filters used in this project are Normalized Least Mean Square (NLMS), Recursive Least Square (RLS), and Discrete Fourier Transform Least Mean Square (DFT-LMS). These adaptive filters follow the same basic system identification workflow, as shown in Figure 1. However, DFT-LMS includes a preprocessing block that transforms the input signal before adaptive filtering. This helps reduce input correlation and improves convergence behavior, especially for coloured noise inputs, as shown in Figure 2.

<img width="400" height="300" alt="System identification drawio" src="https://github.com/user-attachments/assets/a52f8cbe-6ccb-461f-9bac-846532129c92" /> 

Figure 1: System identification block diagram using adaptive filters.

<img width="400" height="300" alt="DFT system identification drawio" src="https://github.com/user-attachments/assets/196c80e8-ccd6-40d2-a489-e7463e33c45f" />

Figure 2: System identification block diagram using DFT-LMS.

**Simulation:** 

1. **Step-size Impact on Tradeoffs:** A 5-tap DFT-LMS adaptive filter was used to estimate the impulse response of a time-varying unknown system modeled as a first-order Markov process. Three different step sizes were evaluated to demonstrate how the step size affects convergence speed, tracking ability, and steady-state error, as shown in the table and figures below.

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/fad5c564-ed29-4e3b-9fb8-427c68f519f1" />

Figure 3: Time-varying estimate & tracking of the unknown system excited by white noise

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/0557b917-a083-4579-bc96-e4005afc7906" />

Figure 4: Time-varying estimate & tracking of the unknown system excited by coloured noise

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/f9bc7cd4-ac9e-4596-bf41-558e4c60d1ea" />

Figure 5: mu = 0.01, ensembled average of the time-varying estimate & tracking of the unknown system 

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/bb122401-d177-4994-9c77-fc82c4aba2ac" />

Figure 6: Time-varying estimate & tracking of the unknown system excited by white noise

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/397cabb4-698e-48b4-bd48-8a05821fed97" />

Figure 7: Time-varying estimate & tracking of the unknown system excited by coloured noise

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/a9e10352-88bc-42c7-8b66-2373e491a5a0" />

Figure 8: mu = 0.05, ensembled average error of the time-varying estimate & tracking of the unknown system 

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/3defe9ce-dd3c-4099-915f-a9e01b784c68" />

Figure 9: Time-varying estimate & tracking of the unknown system excited by white noise

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/23312372-3431-484d-bac6-c4ef6d1ef237" />

Figure 10: Time-varying estimate & tracking of the unknown system excited by coloured noise

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/6ade5cd7-827b-4133-8f8a-8d7a78bdf4c1" />

Figure 11: mu = 1, ensembled average error of the time-varying estimate & tracking of the unknown system 

<img width="450" height="357" alt="image" src="https://github.com/user-attachments/assets/74361326-e0fb-4639-8239-7ff3b5f7c95b" />

Table 1:  Effect of step size on convergence speed, tracking ability, and steady-state error. 

**Conclusion** 
Based on the MSE results, a step size of 0.05 provided the best balance between convergence speed, tracking accuracy, and steady-state error. A convergence factor of 0.01 produced slower convergence and higher tracking lag, while a step size of 1 was too large and caused divergence. Also, with the exception of the divergent case where μ = 1, the TD-LMS adaptive filter performed better in white-noise environments than in coloured-noise environments. This is because white noise has lower input correlation, allowing faster and more reliable adaptation, while coloured noise increases correlation and can slow convergence and reduce tracking accuracy.

2. **Input Signal Impact Tradeoffs:** White-noise and coloured-noise inputs were used with the 5-tap variable-RLS and 5-tap NLMS adaptive filters to demonstrate the impact of input-signal correlation on convergence speed and steady-state error during the identification of the impulse response of a time-invariant unknown system.

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/97b9e389-b421-4a2c-8d97-0eac846611a4" />

Figure 12: Ensembled average error of the time-invariant estimate of the unknown system impulse response excited by coloured noise.

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/aaa9571e-1d8d-454a-ad2b-6fd114edef2f" />

Figure 13: Ensembled average error of the time-invariant estimate of the unknown system impulse response excited by white noise.

**Conclusion**
For highly correlated input signals such as speech, careful adaptive-filter selection is important when identifying a time-invariant unknown system. Figures 12 and 13 show that both NLMS and variable-RLS can identify the system coefficients, but they differ in convergence speed, steady-state MSE, and computational cost.Where variable-RLS provides faster convergence and reaches its steady-state MSE in fewer samples compared to NLMS but with higher computational cost. 
Although these tradeoff considerations were demonstrated using a time-invariant system, the same principles apply to tracking time-varying systems. A variable-RLS algorithm will generally track time-varying system changes much faster than NLMS due to its superior convergence properties, though at the expense of significantly higher computational complexity.

<img width="420" height="320" alt="image" src="https://github.com/user-attachments/assets/9de231c3-6af0-4860-9308-d4b84731f16d" />

Table 2: Tradeoff Considerations


## Repository Structure
```
Adaptive-Filters/
|-- Time-invariant System/       # Time-invariant system identification 
|-- Time-varying System/         # Time-varrying system identification
```
---

## Requirements
MATLAB

- ---

## How to Run

Clone the repository:

```bash
git clone https://github.com/DamiProject/Adaptive-Filters.git
```

Open MATLAB and navigate to the project root directory.


Individual modules can be explored and run from the `Time-invariant System/` and `Time-varying System/`  folder.

---

## Author

**Damilola Awotunde**

MEng, Communications & Signal Processing - Western University | [LinkedIn](https://www.linkedin.com/in/damilola-awotunde) 












