function N = FNorm3D(A)
%  N = FNorm3D(A);
%
%  Returns the F-Norm of the 3D matrix A, simply the sum of the F-norm for
%  all lags.
%
%  Input parameters
%      F      MIMO system matrix 
%
%  Output parameters
%      N      The F-Norm of the input matrix A
%  
%  J. Corr, University of Strathclyde, 25/11/2013
%  Made correct by F. Coutts, University of Strathclyde, 13/01/2016

[~,~,L] = size(A);

Norms = zeros(L,1);

for i= 1:1:L
    Norms(i) = norm(A(:,:,i),'fro')^2;
end
N = sum(Norms);