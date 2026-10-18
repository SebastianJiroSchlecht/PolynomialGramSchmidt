function [Qhat] = PolyMatGramSchmidt(A,threshold,Smoothing,DiffOrd);
% [Q] = PolyMatGramSchmidt(A,threshold,Smoothing,DiffOrder);
%
% Orthonalisation of an analytic matrix A using a modified Gram-Schmidt method.
% For now, A is assumed to be square and full rank.
%
% Input parameters:
%     A          polynomial matrix, MxMxN 
%     threshold  threshold for normalisation and trimming
%     Smoothing  phase smoothing: 'on' or 'off'
%     DiffOrd    differential order for smoothing of normalised vectors
%                (optional; default 3)
%
% Output parameters:
%     Q          MxMxN2 paraunitary matrix
%     R          MxMxN3 upper right triangular polynomial matrix

% S. Weiss, 17/7/25

M = size(A,1);
Aw = A;     
if nargin<4, DiffOrd=3; end;
%disp(sprintf('initial energy in A: %f',PolyMatNorm(A)));
%----------------------------------------------------
%   Iterative modified Gram-Schmidt 
%----------------------------------------------------
for m = 1:M, 
%  disp(sprintf('normalising column %d of %d',[m M]));
  q = PolyVecNormalisation(Aw(:,1,:),threshold,1024,Smoothing,DiffOrd);
  % PolyVecNormalisation needed an internally sharper trimming
  if m==1,         % first column: only perform normalisation
     Qhat = q;
  else             % further columns: accumulate results
     Qold = Qhat;
     LQold = size(Qold,3); Lq = size(q,3);
     Qhat = zeros(M,m,max([LQold, Lq]));
     Qhat(:,1:m-1,1:LQold) = Qold; Qhat(:,m,1:Lq) = q;
  end;   

  % projecting the last columns of A, unless we are in the final iteration
  if m<M,
%     Qperp = - PolyMatTrim( PolyMatConvFast(q,ParaHerm(q)) , threshold );
     Qperp = - PolyMatTrim( PolyMatConvPH2Fast(q,q) , threshold );
     L = size(Qperp,3);
     Qperp(:,:,(L+1)/2) = Qperp(:,:,(L+1)/2)+eye(M);
     Aw = PolyMatTrim( PolyMatConvFast(Qperp,Aw(:,2:end,:)) , threshold );
%     disp(sprintf('energy in A_bar: %f',PolyMatNorm(Aw)));
%     save A_bar.mat Aw

  end;
end;

