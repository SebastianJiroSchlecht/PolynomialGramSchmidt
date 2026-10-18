function [Q,chi,R] = PQRD_analytic_R_ver(A,Nmax,P,Thresh1,Thresh2);
%[Q,chi] = AnalyticEigVecExtract(R,Lambda,Nmax,Th1,Th2)
%
%   Extract approximate analytic eigenvectors from a parahermitian matrix R, 
%   given analytic eigenvalues in Lambda. These are used to order the eigen-
%   vectors in each bin, extract 1D eigenspaces across algebraic multipli-
%   citities, and obtain an analytic eigenvector in each such 1d subspace
%   through finding the phase in each bin that creates the smoothes possible
%   function, or shortest time-domain support.
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
  
% S. Weiss, UoS, 06/11/2021

%------------------------------------------------------------------------------
%   extract parameters
%------------------------------------------------------------------------------


[M,~,N] = size(A);
L=ceil(size(A,3))+5;
Nfft = 2^ceil(log2(max(N,L)));          % next higher power of 2 as DFT length
Nfft_mem=Nfft;
% P = 3;                                  % derivative order

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
  Thres1 = 1e-5;                          % threshold for error in orthonormality
end;


if nargin>4,
  Thres2 = Thresh2;
else
  Thres2 = 1e-5;                          % threshold for error in orthonormality
end;

quiet = 0;

%------------------------------------------------------------------------------
%   iteration for DFT length
%------------------------------------------------------------------------------
Crit = 1; LastCrit = 1; iter = 1;
while Crit == 1,
   %  extract smooth eigenspaces
     % performs QR decomposition in bin-wise
     [~,Rf_new]=PolyMat_QR_DFT(A,Nfft);

  
   %  perform phase smoothing within each eigenspace
   xi = zeros(M,1);
   for m = 1:M,
      for m2=1:M
      q_in(m2,:) = squeeze(Rf_new(m,m2,:));
      end
      [r,dummy] = PhaseSmoothing_QR(q_in,P);
      xi(m) = dummy(end);
      for k = 1:Nfft,
         Rf_new(m,:,k) = r(:,k);
      end;
   end;
   % check on error
   R2 = ifft(Rf_new,Nfft,3);
   R3 = ifft(Rf_new(:,:,1:2:end),Nfft/2,3);
   
   chi = reconst_error_NPH(R2,R3);
   if quiet == 0,   
      disp(sprintf('iteration %d: smoothness %f;    orthonormality %f',[iter sum(xi) chi]));
   end;
   if (2*Nfft<=NfftMax)&(chi>Thres1),
      Rf_old = Rf_new;
      Nfft = 2*Nfft;
      iter = iter+1;
   else
       Crit = 0;
   end;
end;



%------------------------------------------------------------------------------
%   trim result in the time domain
%------------------------------------------------------------------------------
Rlong = ifft(Rf_new,[],3);
Qf_new=QR_DFT_Constrained2(A,Rlong,Nfft_mem);
R=Rf_new;
Nlong = size(Rlong,3);
NormVec = zeros(N,1);
for n = 1:Nlong,
   NormVec(n) = norm(Rlong(:,:,Nlong-n+1),'fro')^2;
end;
NormVec = cumsum(NormVec);
NormVec = NormVec/NormVec(end);
N1 = sum(NormVec<Thres2);
R = Rlong(:,:,1:Nlong-N1);


%------------------------------------------------------------------------------
%   trim result in the time domain
%------------------------------------------------------------------------------
Qlong = ifft(Qf_new,[],3);
Q=Qf_new;
Nlong = size(Qlong,3);
NormVec = zeros(N,1);
for n = 1:Nlong,
   NormVec(n) = norm(Qlong(:,:,Nlong-n+1),'fro')^2;
end;
NormVec = cumsum(NormVec);
NormVec = NormVec/NormVec(end);
N1 = sum(NormVec<Thres2);
Q = Qlong(:,:,1:Nlong-N1);










end

%% Bin-wise QR decomposition performed
function [Q,R]=PolyMat_QR_DFT(A,K)
% assuming given polynomial matrix is not laurent
Af=fft(A,K,3);
for k=1:K
[Q(:,:,k),R(:,:,k)]=qr(Af(:,:,k));
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
end


end

function Q=QR_DFT_Constrained2(A,R,N)
% assuming given polynomial matrix is not laurent
K=size(A,3);
if K<N
    K=N;
end
Af=fft(A,K,3);
Rf=fft(R,K,3);
for i=1:K
Q(:,:,i)=Af(:,:,i)*inv(Rf(:,:,i));
end


end
