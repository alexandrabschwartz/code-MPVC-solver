function [phi, Dphi] = relaxationMPVC_kadrani(a,b,t)

% For two vectors a, b of the same length n and a positive scalar t, this
% function computes
    % phi = (a).*(b-t)
% and the gradient Dphi = [b-t a]

phi = a.*(b-t);
      
if nargout > 1
    % (n_ab x 2) gradient of the relaxation function
    Dphi =  [b-t a];
end