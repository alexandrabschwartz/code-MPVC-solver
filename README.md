# code-MPVC-solver
Some solvers to tackle optimization problems with vanishing constraints. The included functions are described below.

I did test this code, but it may still contain bugs. If you find errors or have suggestions for improvements, please let me know.

## solveMPVC

This function takes a nonlinear optimization problem with vanishing constraints as well as optional options as input. Within the options, you can specify, which solver for vanishing constraints you want to use. At the moment, DIRECT, RELAXATION and RELAXATION_posLB are possible. The function then passes the problem on to the specified solver. If you want to add more solvers, you have to register them here.

## solveMPVC_direct

This function takes a nonlinear optimization problem with vanishing constraints as well as optional options as input. It reformulates the problem into an equivalent nonlinear optimization problem. The reformulation is then solved using an NLP solver of your choice.

## solveMPVC_relaxation

This function takes a nonlinear optimization problem with vanishing constraints as well as an optional options as input. It reformulates the problem into an equivalent nonlinear optimization problem and relaxes the vanishing constraint using a relaxation function of your choice. A sequence of relaxed problems is then solved using an NLP solver of your choice.

Some of these algorithms -- partially in a version for cardinality-constrained problems -- are described here:
* T. Hoheisel, C. Kanzow and A. Schwartz: *Convergence of a local regularization approach for mathematical programmes with complementarity or vanishing constraints* , Optimization Methods and Software  27, 483-512, 2012
* T. Hoheisel, C. Kanzow and A. Schwartz: *Mathematical programs with vanishing constraints: a new regularization approach with strong convergence properties*, Optimization  61, 619-636, 2012
* T. Hoheisel, B. Pablos, A. Pooladian, A. Schwartz and L. Steverango: *A study of one-parameter regularization methods for mathematical programs with vanishing constraints*, Optimization Methods and Software, 2020

## solveMPVC_relaxation_posLB

This function takes a nonlinear optimization problem with vanishing constraints as well as an optional options as input. It reformulates the problem into an equivalent nonlinear optimization problem and relaxes the vanishing constraint using a relaxation function of your choice. A sequence of relaxed problems is then solved using an NLP solver of your choice. Incontrast to solveMPVC_relaxation, this algorithm introduces a small positive lower bound on the function H, which is often used in truss design to ensure regularity of the stiffness matrix.

## setupMPVC_missingData

This function takes a nonlinear optimization problem with vanishing constraints as input. It checks the problem data for completeness and inserts missing data -- if possible -- using default values. E.g. if you did not specify box constraints on the variable, it inserts -inf/inf as lower/upper bounds.

## setupMPVC_defaultOptions

This function takes an options struct as input and sets up missing options using default values. E.g. if you did not specify the solver for the MPVC, it chooses DIRECT.

## relaxationMPVC

This function takes two vectors, a relaxation parameter and a type of relaxation functions as input and passes this data on to the specified relaxation function. At the moment, SCHOLTES, STEFFENSEN, SCHWARTZ and KADRANI are possible. If you want to add more types, you have to register them here.

## relaxationMPVC_scholtes, relaxationMPVC_steffensen, relaxationMPVC_schwartz, relaxationMPVC_kadrani

These functions take two vectors and a relaxation parameter as input and evaulate the relaxation function. The relaxation functions are called after one of the authors of the respective initial papers:
* S. Scholtes: *Convergence properties of a regularization scheme for mathematical programs with complementarity constraints*, SIAM J. Optim. 11, 918–936, 2001
* S. Steffensen and M. Ulbrich: *A new regularization scheme for mathematical programs with equilibrium constraints*, SIAM J. Optim. 20, 2504–2539, 2010
* C. Kanzow and A. Schwartz: *A new regularization method for mathematical programs with complementarity constraints with strong convergence properties*, SIAM J. Optim. 23, 770–798, 2013
* A. Kadrani et al: *A new regularization scheme for mathematical programs with complementarity constraints*, SIAM J. Optim. 20, 78–103, 2009

Note that relaxationMPVC_kadrani and relaxationMPVC_schwartz are modified in comparison to the "standard" complementarity version of these relaxation functions by shifting the kink onto the H-axis.

## testMPVC

This script contains a toy problem to illustrate how the functions are called.
