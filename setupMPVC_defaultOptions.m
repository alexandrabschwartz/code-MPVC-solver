function options = setupMPVC_defaultOptions(options)

% This function is given an empty argument or an options struct with some
% or all of the following fields:
    % options.objectiveGradient = true or false
    % options.constraintsJacobian = true or false
    % options.NLPsolver = 'fmincon' or 'snopt'
    % options.algorithm = 'direct' or 'relaxation' of 'relaxation_posLB'
    % options.relaxation = 'scholtes' or 'steffensen' or 'schwartz' or 'kadrani'
    % options.slacks = true or false
    
% It sets up missing or empty fields using the following default values:
    % options.objectiveGradient = false
    % options.constraintsJacobian = false
    % options.NLPsolver = 'fmincon'
    % options.algorithm = 'direct'
    % options.relaxation = 'scholtes'
    % options.slacks = false

%%

if isempty(options)
    options.objectiveGradient = [];
    options.constraintsJacobian = [];
    options.NLPsolver = [];
    options.algorithm = [];
    options.relaxation = [];
    options.slacks = [];
end

if ~isfield(options, 'objectiveGradient') || isempty(options.objectiveGradient)
    % default value is no gradient information
    options.objectiveGradient = false;
end

if ~isfield(options, 'constraintsJacobian') || isempty(options.constraintsJacobian)
    % default value is no gradient information
    options.constraintsJacobian = false;
end

if ~isfield(options, 'NLPsolver') || isempty(options.NLPsolver)
    % default solver is fmincon
    options.NLPsolver = 'fmincon';
end

if ~isfield(options, 'algorithm') || isempty(options.algorithm)
    % default algorithm is an NLP reformulation
    options.algorithm = 'direct';
end

if ~isfield(options, 'relaxation') || isempty(options.relaxation)
    % default relaxation is scholtes
    options.relaxation = 'scholtes';
end

if ~isfield(options, 'slacks') || isempty(options.slacks)
    % default is no slack variables for the vanishing constraints
    options.slacks = false;
end
