# Ensemble-Kalman-Filter
EnKF in Matlab for a system of two PDEs modeling burglar behavior (model of Short et al.) 

## Note:
I have experienced several bugs in Matlab's solvepde() function. 
For instance, the stated restart from the last time point does not work for systems of PDEs.
My workaround consists in defining an anonymous function that takes the value of the last (point of time of the) solution. This function is taken as the new initial value in the solution process.

# Files
```
Short.m: Solves Short's et al. PDE system with the following coefficients and initial value:
acoeffunction.m: Coefficient a function for our PDE system
ccoeffunction.m: Coefficient c function for our PDE system
fcoeffunction.m: Coefficient f function for our PDE system
u0fun.m: Initial value u0 function for our PDE system

EnKF.m: Ensemble Kalman-Filter performing the aggregation between last aggregated solution simulated 
        to the next time point with the observational data for that point.
EnKF_ens_create.m: creates the ensemble
EnKF_main.m: MAIN

save_fig.m: saves all open plot windows as portable network graphics
subplot4.m: Arranges the solutions as four subplots: artificial solution, 
            noisy observation, solution prior to EnKF and aggregated solution (EnKF).
short_bsp.m: Example script for the unfiltered solution of Short's model.
```

