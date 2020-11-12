function [phi, Dphi] = relaxationMPVC_steffensen(a,b,t)

% For two vectors a, b of the same length and a positive scalar t, this function computes
    % phi = a+b - |a-b|                                    if |a-b| >= t
    %       a+b - t*(2/pi*sin(pi/2 * (a-b)/t + 3pi/2) + 1) else
% and the gradient Dphi = [Dphi_a Dphi_b]

phi =   (a + b - abs(a-b)).*(abs(a-b) >= t) ...
      + (a + b - t*(2/pi*sin(pi/2*(a-b)/t + 3*pi/2) + 1) ).*(abs(a-b) < t) ;

if nargout > 1
    % (n_ab x 2) gradient of the relaxation function
    Dphi =   [zeros(length(a),1)                            2*ones(length(b),1)].*repmat((a-b >= t),1,2)...
           + [2*ones(length(a),1)                           zeros(length(b),1)].*repmat((a-b <= -t),1,2)...
           + [ones(length(a),1)-cos(pi/2*(a-b)/t + 3*pi/2)  ones(length(b),1)+cos(pi/2*(a-b)/t + 3*pi/2)].*repmat((abs(a-b) < t),1,2);
end