function [phi, Dphi] = relaxationMPVC(a, b, t, relaxation)
% For two vectors a and b of the same length and a positive scalar t this 
% function evaluates the specified relaxation function.
    
a = a(:);
b = b(:);

switch lower(relaxation)
    case 'scholtes'
        [phi, Dphi] = relaxationMPVC_scholtes(a,b,t);
    case 'steffensen'
        [phi, Dphi] = relaxationMPVC_steffensen(a,b,t);
    case 'kadrani'
        [phi, Dphi] = relaxationMPVC_kadrani(a,b,t);
    case 'schwartz'
        [phi, Dphi] = relaxationMPVC_schwartz(a,b,t);
end