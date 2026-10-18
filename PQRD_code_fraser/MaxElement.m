function [val,m,n,l] = MaxElement(X);
% [val,m,n,l] = MaxElement(X);
%
% Return the value and coordinates of the maximum value
% for a 3-d array X.
%
% Input parameters:
%     X           3-d array
%
% Output parameters:
%     val         maximum value
%     m           column index
%     n           row index
%     l           time slices index
%
% S Weiss, UoS, 2/11/2004


% Find the matrix contianing te maximum rows and columns using the 3rd dim
[X3,Index3] = max(X,[],3);
% find the maximum vector from the previous matrix
[X2,Index2] = max(X3,[],2);
% find the ultimate maximum from the previous vector
[X1,Index1] = max(X2,[],1);

val = X1;
m = Index1;         % index1 is a singular vaule
n = Index2(m);      % n is the mth element in Index2
l = Index3(m,n);    % l is the m,nth element in Index3
