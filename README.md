# Ensemble-Kalman-Filter
EnKF in Matlab for a system of two PDEs modeling burglar behavior (model of Short et al.) 

## Note:
I have experienced several bugs in Matlab's solvepde() function. 
For instance, the stated restart from the last time point does not work for systems of PDEs.
My workaround consists in defining an anonymous function that takes the value of the last (point of time of the) solution. This function is taken as the new initial value in the solution process.
