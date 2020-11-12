function [x_opt, f_opt, information] = solveMPVC_relaxation(problem, options)

% This function is given an optimization problem with vanishing constraints
% of the form
%    min f(x)  s.t. xl <=   x  <= xu
%                   bl <=  A*x <= bu
%                   cl <= c(x) <= cu
%                   H(x) >= 0, G(x) .* H(x) <= 0
% and solves it by replacing the vanishing constraints with
%                   H(x) >= 0,  phi(G(x),H(x);t) <= 0
% and solving the resulting NLP for decreasing values of t > 0.

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
% respective functions   can either only return the function value or 
% additionally the gradients (oriented row-wise). The default assumption is
% that no gradients are provided. 

% If you want to use gradient information or slack variables for the 
% vanishing constraints, specify the relaxation function phi or the NLP 
% solver, additionally provide an options struct with the fields
    % options.objectiveGradient = true or false
    % options.constraintsJacobian = true or false
    % options.NLPsolver = 'fmincon' or 'snopt'
    % options.slacks =  = true or false
    % options.relaxation = 'scholtes' or 'steffensen' or 'schwartz' or 'kadrani'

% The function returns
    % x_opt                    computed solution
    % f_opt                    objective function value in x_opt
    % information.message      final exit message of the NLP solver
    % information.maxVio_box   maximum violation of box constraints
    % information.maxVio_lin   maximum violation of linear constraints
    % information.maxVio_nln   maximum violation of nonlinear constraints
    % information.maxVio_van   maximum violation of vanishing constraints
    % information.iterations   number of NLPs solved

    
%% parameters
 
t_start = 1; % initial regularization parameter
sigma = 0.1; % reducing factor for the regularization parameter in each iteration
t_min = 10^-8; % lower bound for the regularization parameter

van_tol = 10^-6; % tolerance for G(x) .* H(y) <= tol

iterations = 0; % counter for the number of iterations
iterations_max = 100;  % maximum number of iterations


%% set up missing options using default values

if nargin == 1
    options = [];
end
options = setupMPVC_defaultOptions(options);


%% gather problem data

[problem, n_x, n_lin, n_nln, n_van] = setupMPVC_missingData(problem);


%% define nonlinear reformulation

if options.slacks ~= true
    % define nonlinear reformulation without slack variables
    NLPproblem.objective = problem.objective;
    NLPproblem.xl = problem.xl;
    NLPproblem.xu = problem.xu;
    NLPproblem.A = problem.A;
    NLPproblem.bl = problem.bl;
    NLPproblem.bu = problem.bu;
    NLPproblem.nlcons = @nlcons_with_vanishing; % [c(x); H(x); phi(G(x),H(x);t)]
    NLPproblem.cl = [problem.cl; zeros(n_van,1); -inf(n_van,1)];
    NLPproblem.cu = [problem.cu; inf(n_van,1); zeros(n_van,1)];
    NLPproblem.x_start = problem.x_start;
    NLPproblem.dimension = problem.dimension;   
    
else
    [G_start, H_start] = problem.vancons(problem.x_start);
    
    % define nonlinear reformulation with slack variables y = G(x), z = H(x)
    NLPproblem.objective = @objective_with_slacks;
    NLPproblem.xl = [problem.xl; -inf(n_van,1); zeros(n_van,1)];
    NLPproblem.xu = [problem.xu; inf(n_van,1); inf(n_van,1)];
    NLPproblem.A = [problem.A zeros(n_lin, 2*n_van)];
    NLPproblem.bl = problem.bl;
    NLPproblem.bu = problem.bu;
    NLPproblem.nlcons = @nlcons_with_slacks; % [c(x); G(x)-y; H(x)-z; phi(y,z;t)]
    NLPproblem.cl = [problem.cl; zeros(n_van,1); zeros(n_van,1); -inf(n_van,1)];
    NLPproblem.cu = [problem.cu; zeros(n_van,1); zeros(n_van,1); zeros(n_van,1)];
    NLPproblem.x_start = [problem.x_start; G_start; H_start];
    NLPproblem.dimension = n_x + 2*n_van;  
end


%% solve the nonlinear reformulation for decreasing relaxation parameters

t = t_start;

while (iterations == 0) || ((iterations < iterations_max) && (t > t_min) ...
        && any(G_opt .* H_opt > van_tol))
    
    % solve the relaxed problem
    [X_opt, f_opt, NLPinformation] = solveNLP(NLPproblem, options);
    
    % compute return values
    X_opt = X_opt(:);
    x_opt = X_opt(1:n_x);
    [G_opt, H_opt] = problem.vancons(x_opt);
    
    % update the initial point
    if options.slacks ~= true
        NLPproblem.x_start = x_opt;
    else
        NLPproblem.x_start = [x_opt; G_opt; H_opt];
    end
    
    % decrease the relaxation parameter
    t = t*sigma;
    
    % update the iteration counter
    iterations = iterations + 1;   
end


%% compute return values 

information.iterations = iterations;
information.message = NLPinformation.message;
information.maxVio_box = max([max(x_opt-problem.xu, 0);...
                               max(problem.xl-x_opt, 0)]);
information.maxVio_lin = max([max(problem.A*x_opt-problem.bu, 0);...
                               max(problem.bl-problem.A*x_opt, 0)]);
information.maxVio_nln = max([max(problem.nlcons(x_opt)-problem.cu, 0);...
                               max(problem.cl-problem.nlcons(x_opt), 0)]);
information.maxVio_van = max([max(-H_opt, 0);...
                               max(G_opt .* H_opt, 0)]);


%% auxiliary functions

function [C, DC] = nlcons_with_vanishing(x)
    % rewritten nonlinear/vanishing constraints [c(x); H(x); phi(G(x), H(x); t)]

    if nargout == 1
        c = problem.nlcons(x);
        [G,H] = problem.vancons(x);
        phi = relaxationMPVC(G, H, t, options.relaxation);
        C = [c; H; phi];
        
    elseif nargout > 1
        [c, Dc] = problem.nlcons(x);
        [G, H, DG, DH] = problem.vancons(x);
        [phi, Dphi] = relaxationMPVC(G, H, t, options.relaxation);
        C = [c; H; phi];
        % gradients of the constraints are row vectors
        DC = [Dc;...
              DH;...
              diag(Dphi(:,1))*DG + diag(Dphi(:,2))*DH];
    end
end


function [f, Df] = objective_with_slacks(X)
    % rewritten objective function with slack variables
    % X = (x,y,z)
    X = X(:);
    x = X(1:n_x);
    
    if nargout == 1
        f = problem.objective(x);
        
    elseif nargout > 1
        [f, Df] = problem.objective(x);
        % gradient of the objective is a row vector
        Df = [Df zeros(1,2*n_van)];
    end
end


function [C,DC] = nlcons_with_slacks(X)
    % X = (x,y,z)
    % rewritten nonlinear/vanishing constraints [c(x); G(x)-y; H(x)-z; phi(y,z;t)]
    
    X = X(:);
    x = X(1:n_x);
    y = X(n_x+1:n_x+n_van);
    z = X(n_x+n_van+1:end);

    if nargout == 1
        c = problem.nlcons(x);
        [G, H] = problem.vancons(x);
        phi = relaxationMPVC(y, z, t, options.relaxation);
        C = [c; G-y; H-z; phi];
        
    elseif nargout > 1
        [c, Dc] = problem.nlcons(x);
        [G, H, DG, DH] = problem.vancons(x);
        [phi, Dphi] = relaxationMPVC(y, z, t, options.relaxation);
        C = [c; G-y; H-z; phi];
        % gradients of the constraints are row vectors
        DC = [Dc               zeros(n_nln,n_van) zeros(n_nln,n_van);...
              DG               -eye(n_van)        zeros(n_van,n_van);...
              DH               zeros(n_van,n_van) -eye(n_van);...
              zeros(n_van,n_x) diag(Dphi(:,1))    diag(Dphi(:,2))];
    end
end

end