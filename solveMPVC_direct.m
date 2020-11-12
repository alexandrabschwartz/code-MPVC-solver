function [x_opt, f_opt, information] = solveMPVC_direct(problem, options)

% This function is given an optimization problem with vanishing constraints
% of the form
%    min f(x)  s.t. xl <=   x  <= xu
%                   bl <=  A*x <= bu
%                   cl <= c(x) <= cu
%                   H(x) >= 0, G(x) .* H(x) <= 0
% and solves it by intepretting as a standard NLP.

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

% If you want to use gradient information of slack variables for the 
% vanishing constraints, specify the  NLP solver, additionally provide an 
% options struct with the fields
    % options.objectiveGradient = true or false
    % options.constraintsJacobian = true or false
    % options.NLPsolver = 'fmincon' or 'snopt'
    % options.slacks =  = true or false

% The function returns
    % x_opt                    computed solution
    % f_opt                    objective function value in x_opt
    % information.message      exit message of the NLP solver
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
    NLPproblem.nlcons = @nlcons_with_vanishing; % [c(x); H(x); G(x).*H(x)]
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
    NLPproblem.nlcons = @nlcons_with_slacks; % [c(x); G(x)-y; H(x)-z; y.*z]
    NLPproblem.cl = [problem.cl; zeros(n_van,1); zeros(n_van,1); -inf(n_van,1)];
    NLPproblem.cu = [problem.cu; zeros(n_van,1); zeros(n_van,1); zeros(n_van,1)];
    NLPproblem.x_start = [problem.x_start; G_start; H_start];
    NLPproblem.dimension = n_x + 2*n_van;  
end



%% solve the nonlinear reformulation
[X_opt, f_opt, NLPinformation] = solveNLP(NLPproblem, options);


%% compute return values
X_opt = X_opt(:);
x_opt = X_opt(1:n_x);

information.iterations = 1;
information.message = NLPinformation.message;
information.maxVio_box = max([max(x_opt-problem.xu, 0);...
                               max(problem.xl-x_opt, 0)]);
information.maxVio_lin = max([max(problem.A*x_opt-problem.bu, 0);...
                               max(problem.bl-problem.A*x_opt, 0)]);
information.maxVio_nln = max([max(problem.nlcons(x_opt)-problem.cu, 0);...
                               max(problem.cl-problem.nlcons(x_opt), 0)]);
[G_opt, H_opt] = problem.vancons(x_opt);                          
information.maxVio_van = max([max(-H_opt, 0);...
                               max(G_opt .* H_opt, 0)]);                           

%% auxiliary functions

function [C, DC] = nlcons_with_vanishing(x)
    % rewritten nonlinear/vanishing constraints [c(x); H(x); G(x).*H(x)]

    if nargout == 1
        c = problem.nlcons(x);
        [G, H] = problem.vancons(x);
        C = [c; H; G.*H];
        
    elseif nargout > 1
        [c, Dc] = problem.nlcons(x);
        [G, H, DG, DH] = problem.vancons(x);
        C = [c; H; G.*H];
        % gradients of the constraints are row vectors
        DC = [Dc;...
              DH;...
              diag(H)*DG + diag(G)*DH];
    end
end

function [f, Df] = objective_with_slacks(X)
    % rewritten objective function with slack variables
    % X = [x; y; z]
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
    % X = [x; y; z]
    % rewritten nonlinear/vanishing constraints [c(x); G(x)-y; H(x)-z; y.*z]
    
    X = X(:);
    x = X(1:n_x);
    y = X(n_x+1:n_x+n_van);
    z = X(n_x+n_van+1:end);

    if nargout == 1
        c = problem.nlcons(x);
        [G, H] = problem.vancons(x);
        C = [c; G-y; H-z; y.*z];
        
    elseif nargout > 1
        [c, Dc] = problem.nlcons(x);
        [G, H, DG, DH] = problem.vancons(x);
        C = [c; G-y; H-z; y.*z];
        % gradients of the constraints are row vectors
        DC = [Dc               zeros(n_nln,n_van) zeros(n_nln,n_van);...
              DG               -eye(n_van)        zeros(n_van,n_van);...
              DH               zeros(n_van,n_van) -eye(n_van);...
              zeros(n_van,n_x) diag(z)            diag(y)];
    end
end

            
end