%%  toy example from truss optimization to test MPVC solvers

%   min 4x_1 + 2x_2 s.t. x_1 >= 0,
%                        (5*sqrt(2) - x_1 - x_2) * x_1 <= 0, 
%                        x_2 >= 0,
%                        (5 - x_1 - x_2) * x_2 <= 0

% [0; 0] is the global minimum and an isolated feasible point
% [0; 5] is a local minimum
% [0; 5*sqrt(2)] is not a local solution

problem.objective = @objective_truss;
problem.xl = [0; 0];
problem.xu = [inf; inf]; 
problem.A =[];
problem.bl = [];
problem.bu = [];
problem.nlcons = [];
problem.cl = [];
problem.cu = [];
problem.vancons = @vancons_truss;
% problem.x_start = [3; 3];
problem.dimension = 2 ;

options.objectiveGradient = true;
options.constraintsJacobian = true;
options.solver = 'fmincon';
options.slacks = false;


%% call of the MPVC solvers

tol = 10^-2;
figure

% direct solver

subplot(2,5,1)
title('direct')
xlim([-5,20]);
ylim([-5,20]);
hold on

counter_global = 0;
counter_local = 0;
counter_nonopt = 0;
counter_other = 0;
counter_iterations = 0;

options.algorithm = 'direct';
for x = -5:20
    for y = -5:20
        problem.x_start = [x; y];
        [x_opt, f_opt, information] = solveMPVC(problem,options);
        counter_iterations = counter_iterations + information.iterations;

        if norm(x_opt - [0; 0]) <= tol
            % x^0 = [0; 0] is the global minimum
            plot(x,y, 'bo')
            counter_global = counter_global + 1;
        elseif norm(x_opt - [0; 5]) <= tol
            % x^* = [0; 5] is the local minimum
            plot(x,y, 'r.')
            counter_local = counter_local + 1;
        elseif norm(x_opt - [0; 5*sqrt(2)]) <= tol
            % x^+ = [0; 5*sqrt(2)] is not a local solution
            plot(x,y, 'k+')
            counter_nonopt = counter_nonopt + 1;
        else
            counter_other = counter_other + 1;
        end
    end
end

disp([options.algorithm  ' ==========='])
[counter_global  counter_local counter_nonopt counter_other counter_iterations]

% relaxation algorithms

options.algorithm = 'relaxation';
relaxations = {'scholtes', 'steffensen', 'schwartz', 'kadrani'};
for relax = 1:4
    options.relaxation = relaxations{relax};
    
    subplot(2,5,relax+1)
    title(['relaxation: ' relaxations{relax}])
    xlim([-5,20]);
    ylim([-5,20]);
    hold on
    
    counter_global = 0;
    counter_local = 0;
    counter_nonopt = 0;
    counter_other = 0;
    counter_iterations = 0;
    
    for x = -5:20
        for y = -5:20
            problem.x_start = [x; y];
            [x_opt, f_opt, information] = solveMPVC(problem,options);
            counter_iterations = counter_iterations + information.iterations;

            if norm(x_opt - [0; 0]) <= tol
                % x^0 = [0; 0] is the global minimum
                plot(x,y, 'bo')
                counter_global = counter_global + 1;
            elseif norm(x_opt - [0; 5]) <= tol
                % x^* = [0; 5] is the local minimum
                plot(x,y, 'r.')
                counter_local = counter_local + 1;
            elseif norm(x_opt - [0; 5*sqrt(2)]) <= tol
                % x^+ = [0; 5sqrt(2)] is not a local solution
                plot(x,y, 'k+')
                counter_nonopt = counter_nonopt + 1;
            else
                counter_other = counter_other + 1;
            end
        end
    end
    
    disp([options.algorithm ' ' options.relaxation ' ==========='])
    [counter_global  counter_local counter_nonopt counter_other counter_iterations]
end


% relaxation algorithms with positive lower bound on H

options.algorithm = 'relaxation_posLB';
relaxations = {'scholtes', 'steffensen', 'schwartz', 'kadrani'};
for relax = 1:4
    options.relaxation = relaxations{relax};
    
    subplot(2,5,relax+6)
    title(['relaxation with posLB: ' relaxations{relax}])
    xlim([-5,20]);
    ylim([-5,20]);
    hold on
    
    counter_global = 0;
    counter_local = 0;
    counter_nonopt = 0;
    counter_other = 0;
    counter_iterations = 0;
    
    for x = -5:20
        for y = -5:20
            problem.x_start = [x; y];
            [x_opt, f_opt, information] = solveMPVC(problem,options);
            counter_iterations = counter_iterations + information.iterations;

            if norm(x_opt - [0; 0]) <= tol
                % x^0 = [0; 0] is the global minimum
                plot(x,y, 'bo')
                counter_global = counter_global + 1;
            elseif norm(x_opt - [0; 5]) <= tol
                % x^* = [0; 5] is the local minimum
                plot(x,y, 'r.')
                counter_local = counter_local + 1;
            elseif norm(x_opt - [0; 5*sqrt(2)]) <= tol
                % x^+ = [0; 5sqrt(2)] is not a local solution
                plot(x,y, 'k+')
                counter_nonopt = counter_nonopt + 1;
            else
                counter_other = counter_other + 1;
            end
        end
    end
    
    disp([options.algorithm ' ' options.relaxation ' ==========='])
    [counter_global  counter_local counter_nonopt counter_other counter_iterations]
end





%% objective functions and nonlinear constraints

function [f, Df] = objective_truss(x)
    f = [4 2] * x;
    if nargout > 1
        Df = [4 2];
    end
end

function [G, H, DG, DH] = vancons_truss(x)
    G = [5*sqrt(2) - x(1) - x(2);...
         5 - x(1) - x(2)];
     H = x;
    if nargout > 1
        DG = [-1 -1; -1 -1];
        DH = [1 0; 0 1];
    end
end


        