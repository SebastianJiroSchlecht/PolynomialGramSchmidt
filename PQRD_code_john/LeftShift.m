function G = LeftShift(F, K, Shift)
% apply delay matrix from left hand side

% Inputs:
%       F       the matrix that the shift is being applied to.
%       K       the col index that the shift is to be applied to. 
%       Shift   the delay that we want shift F by. 

[M,N,L] = size(F);
if Shift == 0
    G = F;
elseif Shift > 0
    G = zeros(M,N,L+Shift);
    G(:,:,1:L) = F;
    G(K,:,1+Shift:L+Shift) = F(K,:,:);
    G(K,:,1:Shift) = zeros(1,N,Shift);
elseif Shift < 0
    G = zeros(M,N,L-Shift);
    G(:,:,1-Shift:L-Shift) = F;
    G(K,:,1:L) = F(K,:,:);
    G(K,:,L+1:L-Shift) = zeros(1,N,-Shift);
end;