function [q,error_P] = PolyVecNormalisation(a,eps_P,MaxDFT,Smoothing,DiffOrd);
% q = PolyVecNormalisation(a,eps_P,MaxDFT)
% 
% q = PolyVecNormalisation(a,MaxDFT) performs the normalisation of a polynomial 
% vector a. The result is guaranteed to be analytic [1], and the order of the
% approximation is refined until either a maximum DFT size MaxDFT is reached 
% or the projection error [1] falls below the threshold eps_P [1]. The return
% parameter q is the normalised vector.
% 
% [q,error_P] = PolyVecNormalisation(a,MaxDFT) additionally returns the reached
% projection error in the variable error_P.
% 
% The implementation reflects Algorithm 2 in [1].
%
% Input parameters
%   a          Mx1xN1 polynomial vector
%   eps_P      threshold for projection error
%   MaxDFT     maximum DFT size
%   Smoothing  phase smoothing: 'on' or 'off'
%   DiffOrd    differential order for smoothing (opt., default 3)
%
% Output parameter
%   q        normalised vector of dimension Mx1xN2
%   error_P  final projection error
%
% Reference:
% [1] F.A. Khattak, I.K. Proudler, S.J. Schlecht, and S. Weiss: "..."

% S. Weiss, UoS, 7/12/2024
%  trimming added, 4/7/2025  
%  DFT implementation of projection criterion 27/1/26: 

%------------------------------------------------------------------------------
% parameters and initialisation
%------------------------------------------------------------------------------
% ensure that we can perform at least one iteration at the maximum DFT size
if size(a,3)>MaxDFT/2, 
   a = PolyMatCentreTrim(a,MaxDFT/2); 
   disp('previous solution truncated to continue iteration with next column');
end;   
[M,L,N1] = size(a);
if L>1,
   disp('input dimension incorrect');
end;
a = squeeze(a);
K = 2^(ceil(log2(N1))-1);
error_P = 10*eps_P;

if nargin<4, DiffOrd = 3; end;

K = 2*K;                             % starting DFT size
af = fft(a,K,2);             
%------------------------------------------------------------------------------
% iteration for sufficient DFT size
%------------------------------------------------------------------------------
while and( (error_P>eps_P),(K<=MaxDFT/2) ),
   DiffOrd = max(11-ceil(log2(K)),2);    % high differential order for shorter DFTs
   qf = PolyVecDFTNormalisation(af,Smoothing,DiffOrd);  
   % check projection error
   q2f = fft(ifft(qf,K,2),2*K,2);   % interpolate qf for DFT-domain convolution
   K = 2*K;
   af = fft(a,K,2);                 % for next iteration, but also for check
   % check of projection criterion
   qaf = sum(conj(q2f).*af);
   error_P = norm(af - q2f.*(qaf(ones(M,1),:)),'fro')^2;
   % disp(sprintf('DFT size K = %d :  projecton error %0.5g',[K error_P]));                  
end;       
DiffOrd = max(11-ceil(log2(K/2)),2);    % high differential order for shorter DFTs
% trim
q = PolyMatTrim(permute(ifft(qf,K/2,2),[1 3 2]), 1e-10);

%disp('done');

%------------------------------------------------------------------------------
% trimming of leading or trailing zeros
%------------------------------------------------------------------------------
theta = eps_P;          % threshold values;  everything below threshold is set to zero
NonZeroIndices = find(squeeze(sum(abs(q),1))/M > theta);
q = q(:,:,NonZeroIndices(1):NonZeroIndices(end));


