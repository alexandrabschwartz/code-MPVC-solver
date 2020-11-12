function [phi, Dphi] = relaxationMPVC_scholtes(a,b,t)

% For two vectors a, b of the same length and a positive scalar t, this 
% function computes
    % phi = a.*b -t
% and the gradient Dphi = [b a]

phi = a.*b - t;
      
if nargout > 1
    % (n_ab x 2) gradient of the relaxation function
    Dphi =  [b a];
end