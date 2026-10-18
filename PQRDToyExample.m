function [A,Q,R] = PQRDToyExample(seed)
% Reconstruct Example 6 from equation (52) of Khattak et al., TSP 2026.
% The original helper was not supplied; the seed controls the PEVD random Q.
if nargin<1, seed=0; end
Q = PUPolyMatRand(2,2,seed,'complex');
R = zeros(2,2,2);
R(:,:,1) = [1 .5; 0 1];
R(:,:,2) = [1i -.5; 0 -1];
A = PolyMatConvFast(Q,R);
end
