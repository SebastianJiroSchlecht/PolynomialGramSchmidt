function [EnergyBD,EnergyAD,EnergyZerolag,EnergyDiag] = CalcEnergies(A,zerolag);

% Input parameters:
%
%   A             polynomial matrix  
%   zerolag       position of zero lag coefficient matrix in A
%
% Output parameters:
%
%   EnergyBD        squared F-norm of all elements below the diag of A
%   EnergyAD        squared F-norm of all elements above the diag of A
%   EnergyZerolag   squared F-norm of the zerolag coefficient matrix of A
%   EnergyDiag      squared F-norm of each diagonal elements of A in a vector

[a1,a2,a3] = size(A);

MaskOffD = toeplitz([0;ones(a1-1,1)],[0;ones(a2-1,1)]');   % Ones on all off-diag elements, zeros on diagonal
MaskOnD = ones(a1,a2) - MaskOffD;
MaskLT = toeplitz([0;ones(a1-1,1)],zeros(1,a2));           % Ones below the diagonal, zeros everywhere else
MaskUT = toeplitz(zeros(a1,1),[0;ones(a2-1,1)]');          % Ones above the diagonal, zeros everywhere else

for t = 1:size(A,3)
    EnA(:,:,t) = A(:,:,t).*MaskUT;
    EnB(:,:,t) = A(:,:,t).*MaskLT;
    EnD(:,:,t) = A(:,:,t).*MaskOnD;
end;

EnergyBD = sum(sum(sum(EnB.*conj(EnB))));     % Squared F-norm of all element below the diagonal
EnergyAD = sum(sum(sum(EnA.*conj(EnA))));     % Squared F-norm of all element above the diagonal
EnergyZerolag = sum(sum(A(:,:,zerolag).*conj(A(:,:,zerolag))));     % Sq F-Norm of zerolag coefficient matrix
EnergyDiag = diag(sum(EnD.*conj(EnD),3))';    % vector containin all squared F-norm of all diagonal elements of A

% Check:
% Total = sum(sum(sum(A.*conj(A))))
% Below = EnergyBD
% Above = EnergyAD
% Diag = EnergyDiag
% Total2 = EnergyBD+EnergyAD+sum(EnergyDiag)
