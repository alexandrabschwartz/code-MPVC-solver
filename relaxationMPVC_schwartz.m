function [phi, Dphi] = relaxationMPVC_schwartz(a,b,t)

% For two vectors a, b of the same length and a positive scalar t, this
% function computes
    % phi(a,b) = a*(b-t)               if    a+b >= t
    %            -0.5*(a^2 - (b-t)^2)  else
% and the gradient Dphi = [Dphi_a Dphi_b]

phi =   a.*(b-t).*(a+b >= t) ...
      - 0.5*(a.^2 + (b-t).^2).*(a+b < t) ;
      
if nargout > 1
    % (n_ab x 2) gradient of the relaxation function
    Dphi =   [b-t a].*repmat((a+b >= t),1,2) ...
           - [a   b-t].*repmat((a+b < t),1,2);
end