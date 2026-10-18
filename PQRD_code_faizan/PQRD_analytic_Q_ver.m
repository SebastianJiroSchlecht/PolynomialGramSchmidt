function [Q,chi,R] = PQRD_analytic_Q_ver(R,Nmax,P,Thresh1,Thresh2);
% function provided by Faizan
%
%   Input parameters
%       R        MxMxN polynomial matrix (time-domain)
%       Nmax     maximum FFT length for iteration (default 256)
%       P        Derivative order for smoothness metric
%       Th1      threshold for orthonormality (default 5e-5)
%       Th2      threshold for trimming (default 1e-4)
%
%   Output parameters
%       Q        MxMx? matrix of eigenvectors (time domain)
%       chi      error in paraunitarity

%[Q,chi] = AnalyticEigVecExtract(R,Lambda,Nmax,Th1,Th2)
%
%   Extract approximate analytic eigenvectors from a parahermitian matrix R, 
%   given analytic eigenvalues in Lambda. These are used to order the eigen-
%   vectors in each bin, extract 1D eigenspaces across algebraic multipli-
%   citities, and obtain an analytic eigenvector in each such 1d subspace
%   through finding the phase in each bin that creates the smoothes possible
%   function, or shortest time-domain support.

  


%------------------------------------------------------------------------------
%   extract parameters
%------------------------------------------------------------------------------


[M,~,N] = size(R);
L=ceil(size(R,3));
Nfft = 2^ceil(log2(max(N,L))+1);          % next higher power of 2 as DFT length
Nfft_mem=Nfft;
% P = 3;                                  % derivative order
max_iter=150;
if nargin>1,
  NfftMax = Nmax;
else
  NfftMax = 256;                          % threshold for error in orthonormality
end;

if nargin >2
else
    P=3;
end



if nargin>3,
  Thres1 = Thresh1;
else
  % changed from 2e-5   by SW 27/8/25
  Thres1 = 2e-5;                          % threshold for error in orthonormality
end;


if nargin>4,
  Thres2 = Thresh2;
else
  % changed from 1e-6 by SW 27/8/25
  Thres2 = 1e-6;                          % threshold for error in orthonormality
end;

quiet = 0;

%------------------------------------------------------------------------------
%   iteration for DFT length
%------------------------------------------------------------------------------
Crit = 1; LastCrit = 1; iter = 1;
while Crit == 1,

     % performs QR decomposition in bin-wise
     Qf_new=PolyMat_QR_DFT(R,Nfft);
  
   %  perform phase smoothing within each eigenspace
   xi = zeros(M,1);
   for m = 1:M,
      q_in = squeeze(Qf_new(:,m,:));
%      disp(sprintf('before phase smoothing for element %d',m));
      [q,dummy,chi_phi,iter] = PhaseSmoothing_QR(q_in,P,max_iter);
%      disp('after phase smoothing');
      fcn_val(m)=dummy(end);
      xi(m) = chi_phi;
      for k = 1:Nfft,
         Qf_new(:,m,k) = q(:,k);
      end;
   end;
   xi=mean(xi);
   % check on error
   Qlong=ifft(Qf_new,[],3);
   %chi=reconst_error_NPH(Qlong,ifft(Qf_new(:,:,1:2:end),[],3));
   %chi=paraunitarity_error(Qlong);
   QQ = PolyMatConv(ParaHerm(Qlong),Qlong);
   QQ(:,:,(size(QQ,3)+1)/2) = QQ(:,:,(size(QQ,3)+1)/2) - eye(M);
   chi = PolyMatNorm(QQ);

%if (2*Nfft<=NfftMax) & (chi>Thres1 | xi>1e-5),
if (2*Nfft<=NfftMax) & (chi>Thres1 | xi>1e-5),    % changed by SW, 27/8/25
      
      Nfft = 2*Nfft;
    
else
   
       Crit = 0;
    
    
   end;
end;




%------------------------------------------------------------------------------
%   trim result in the time domain
%------------------------------------------------------------------------------
Q=Qlong;
R=PolyMatConv(ParaHerm(Q),R);




















end

%% Bin-wise QR decomposition performed
function [Q,R]=PolyMat_QR_DFT(A,K)
% assuming given polynomial matrix is not laurent
Af=fft(A,K,3);
for k=1:K
[Q(:,:,k),R(:,:,k)]=qr(Af(:,:,k));
end


end


function R=QR_Constrained(A,Q,N)
% assuming given polynomial matrix is not laurent
K=size(A,3);
if K<N
    K=N;
end
Af=fft(A,K,3);
Qf=Q;
for i=1:K
R(:,:,i)=Qf(:,:,i)'*Af(:,:,i);
end


end

function R=QR_DFT_Constrained(A,Q,N)
% assuming given polynomial matrix is not laurent
K=size(A,3);
if K<N
    K=N;
end
Af=fft(A,K,3);
Qf=fft(Q,K,3);
for i=1:K
R(:,:,i)=Qf(:,:,i)'*Af(:,:,i);
R(:,:,i)=triu(R(:,:,i));
end


end

