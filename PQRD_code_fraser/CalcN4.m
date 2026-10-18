function N = CalcN4(A)
%  N = FNorm3D(A);
%
%  Returns the N4 of the 3D matrix A, simply the sum of the squared F-norm 
%  for all lags.
%
%  Input parameters
%      F      M x K x L MIMO system matrix 
%             M   output dimension
%             K   input dimension
%             L   lenght of FIR filters
%
%  Output parameters
%      N      The F-Norm of the input matrix A
%  
%  J. Corr, University of Strathclyde, 25/11/2013

% [~,~,L] = size(A);
% 
% Norms = zeros(L,1);
% 
% for i= 1:1:L
%     Norms(i) = norm(A(:,:,i),'fro')^2;
% end
% N = sum(Norms);

Renergy = sum(sum(sum(A.*conj(A))));
N = sqrt(Renergy);