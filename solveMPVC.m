function [x_opt, f_opt, information] = solveMPVC(problem, options)

% This function is given an optimization problem with vanishing constraints
% of the form
%    min f(x)  s.t. xl <=   x  <= xu
%                   bl <=  A*x <= bu
%                   cl <= c(x) <= cu
%                   H(x) >= 0, G(x) .* H(x) <= 0
% and solves it using different algorithms.

% The problem should be provided as a struct with the following fields: 
    % problem.objective = @objective
    % problem.xl = xl
    % problem.xu = xu
    % problem.A = A
    % problem.bl = bl
    % problem.bu = bu
    % problem.nlcons = @nlcons
    % problem.cl = cl
    % problem.cu = cu
    % problem.vancons = @vancons
    % problem.x_start = x_start
    % problem.dimension = n_x 
% For the objective function and the nonlinear/vanishing constraints the 
% respective functions can either only return the function value or 
% additionally the gradients (oriented row-wise). The default assumption is
% that no gradients are provided. 

% If you want to use gradient information or slack variables for the 
% vanishing constraints, specify the MPVC algorithm or the  NLP solver, 
% additionally provide an options struct with the fields
    % options.objectiveGradient = true or false
    % options.constraintsJacobian = true or false
    % options.NLPsolver = 'fmincon' or 'snopt'
    % options.algorithm = 'direct' or 'relaxation' or 'relaxation_posLB'
    % options.relaxation = 'scholtes' or 'steffensen' or 'schwartz' or 'kadrani'
    % options.slacks =  = true or false

% The function returns
    % x_opt                    computed solution
    % f_opt                    objective function value in x_opt
    % information.message      exit message of the solver
    % information.maxVio_box   maximum violation of box constraints
    % information.maxVio_lin   maximum violation of linear constraints
    % information.maxVio_nln   maximum violation of nonlinear constraints
    % information.maxVio_van   maximum violation of vanishing constraints
    % information.iterations   number of NLPs solved


%% set up missing options using default values

if nargin == 1
    options = [];
end
options = setupMPVC_defaultOptions(options);



%% call the specified solution algorithm for MPVC problems

switch options.algorithm
    case 'direct'
        [x_opt, f_opt, information] = solveMPVC_direct(problem, options);
    case 'relaxation'
        [x_opt, f_opt, information] = solveMPVC_relaxation(problem, options);
    case 'relaxation_posLB'
        [x_opt, f_opt, information] = solveMPVC_relaxation_posLB(problem, options);
    otherwise
        disp('Unknown MPVC algorithm, will use direct NLP reformulation instead')
        [x_opt, f_opt, information] = solveMPVC_direct(problem, options);
end
        
        

