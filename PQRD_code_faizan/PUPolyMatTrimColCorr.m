function [H2,rho] = PUPolyMatTrimColCorr(H1,Mu)
% [H2,rho] = PUPolyMatTrimRowCorr(H1,Mu)
%
%  H2=PUPolyMatTrimRpwCorr(H1) trims the time dimension of a paraunitary
%  matrix by removing outer matrix coefficients from each row individually.
%  The trimming is performed from both ends, and will remove the matrix 
%  coefficient with the smallest Frobenius norm, such that the total energy
%  of the trimmed parts does not exceed one permille of the total energy in 
%  H1.
% 
%  H2=PUPolyMatTrimRowCorr(H1,gamma) trims such that the ratio between the 
%  removed and total energy of H1 is less than gamma, with 0<=gamma<1.
%
%  [H2,rho]=PUPolyMatTrimRowCorr(H1) or [H2,rho]=PUPolyMatTrimRowCorr(H1,r)
%  additionally returns the error in paraunitarity, i.e. the norm of  
%  (H2*~H2-I) , after trimming in the variable rho.
%  
%  This function is the Row-Shift 'trim' operation used in [1] implemented
%  by applying the PUPolyMatTrim to each row of the paraunitary matrix H1.
%
%  Input parameters:
%     H1      paraunitary matrix
%     gamma   maximum ratio between the maximum energy removed at outer lags 
%             and the total energy 
%             default: 1/1000
%             Compared to PUPolyMatTrim gamma typically be set 5 times
%             higher to obtain similar error in paraunitary (rho)
%  
%  Output parameters:
%     H2      trimmed parahermitian matrix
%     rho     error in paraunitarity
%
%  References:
%  [1] J. Corr, K. Thompson, S. Weiss, I.K. Proudler, J.G. McWhirter,
%      "Row-Shift Corrected Truncation of Paraunitary Matrices for PEVD 
%      Algorithms", 23rd European Signal Processing Conference, Nice,
%      France, Sept. 2015.
%
% See also: PUPolyMatTrim

% J. Corr, University of Strathclyde, 11/12/2014
%
% Last modified 26/8/2015

[M,N,L] = size(H1);
RowMu = Mu/M; % quite simplistic there may be a better approach...
H2 = zeros(M,N,L);
for i=1:N
    Row = H1(:,i,:);             % Extract Col
    Row = PUPolyMatTrim(Row,RowMu); % Truncate Col
    [~,~,RowLen] = size(Row);
    H2(:,i,1:RowLen) = Row;      % Place new row in old row (at the beginning of the matrix)
    %Deltas(i) = RDelta;
end
% Remove Zeros
[H2,~] = PUPolyMatTrim(H2,0);

% prepare additional output parameter if error in paraunitarity is requested
if nargout>1,
   HH2 = PolyMatConv(H2,ParaHerm(H2));
   L2 = (size(HH2,3)+1)/2;
   HH2(:,:,L2) = HH2(:,:,L2)-eye(N);
   rho = PolyMatNorm(HH2);
end;

end
