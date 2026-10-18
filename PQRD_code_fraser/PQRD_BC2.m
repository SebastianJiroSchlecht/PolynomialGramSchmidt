function [Q,Gamma] = PQRD_BC2(R,maxiter,epsilon,PHMu,PUMu)

% [Q, R] = PQRD_BC2(H,maxiter,thresh)
%
% Implementation of the Polynomial QR Decomposition by Columns algorithm
% (PQRD_BC). Computes H(z) = Q(z)R(z) where Q(z) is a paraunitary
% polynomial matrix and R(z) is an upper right triangular polynomial
% matrix.
%
% The implementation is based on the paper "A Novel Algorithm for 
% calculating the QR Decomposition of a Polynomial Matrix" by McWhirter et 
% al. but with a few differences. The Givens Rotation is  modified to deal 
% with complex channels. This implementation also does not include the 
% paraunitary inverse delay matrix.
%  
% Input parameters:
%    H            channel polynomial matrix, H(z)
%    maxiter      maximum number of iterations for each Elementary 
%                 Polynomial Givens Rotation (EPGR)
%                 (default 50)  
%    thresh       after each EPGR the Q(z) and R(z) matrices are shortened.
%                 The polynomials are trimmed by discarding the first and
%                 last elements whose absolute value is below the variable
%                 thresh. Setting thresh to 0 means no elements are
%                 discarded
%                 (default 10^-5)  
%
% Output parameters:
%    Q            paraunitary polynomial matrix, Q(z)
%    R            upper right triangular polynomial matrix, R(z)
%
%-----------------------  FUNCTIONS REQUIRED  ----------------------------%
%                                                                         %
%                         ParaHerm , MIMOConv                             %
%                                                                         %
%-------------------------------------------------------------------------%
% Fraser K Coutts, 
% Centre for Signal and Image Processing
% University of Strathclyde, 11/05/2016

% DM = zeros(maxiter+1,1);
% T = zeros(maxiter+1,1);
% O = zeros(maxiter+1,1);
N = size(R,1);

if ~exist('maxiter'),
	maxiter = 400;
end;
if ~exist('epsilon'),
	epsilon = 0.0001;
end;

[M,N,L] = size(R); % get channel dimensions
Q(:,:,1) = eye(M);

% O(1,1) = 1;
% O(1,2) = L;

% Calculate N4 (the F-Norm of the entire system squared)
% PHref = (FNorm3D(R))^2;
% PUref = (FNorm3D(Q))^2;
PHref = FNorm3D(R);
PUref = FNorm3D(Q);

% Calculate constants for masking process
MPlus1 = M+1;
NSlice = M*N;

% Initilaise iteration variables
stopcrit = 0;
iter = 0;
% totaliter = 0;
% DM(1) = DiagMeasUpperTri(R);

% append zeros to the channel matrix so the zero lag is in the centre %
% Gamma = zeros(M,N,2*L-1);
% Gamma(:,:,L:end) = R;
Gamma = R;

totalmeas = 0;
tt = tic;
   
for n = 1:min(M-1,N), % loop over columns - starting from first through to last

    while stopcrit == 0,
       iter = iter+1;
%        totaliter = totaliter+1;
       ttmeas = tic;
%        DM(totaliter+1) = DiagMeasUpperTri(Gamma);
%        totalmeas = totalmeas + toc(ttmeas);
%        T(totaliter+1) = toc(tt)-totalmeas; 
       
       %Gamma = R;
       [~,~,L] = size(Gamma);

       % Zeolag needs re-aligned each iteration as R grows
       zerolag = (L+1)/2; 

       R1 = Gamma;
       NT = NSlice*L;

       for start = 1:M:NSlice
           for l = start:(MPlus1):NSlice
              R1(l:NSlice:NT)=0;
           end
       end
       R1(:,[1:n-1 n+1:N],:) = 0;

       [maxEle,m,~,l] = MaxElement(abs(R1)); % Find location and value of max

       if maxEle>epsilon,
           
          delay = zerolag-l;
          
          a10 = Gamma(n,n,zerolag);
          a2t = Gamma(m,n,l);
          
          if abs(a10) > 0
              theta = atan(abs(a2t)/abs(a10));
          else
              theta = pi/2;
          end
          alpha = -angle(a10);
          phi = -angle(a2t);

          Gamma = AdvanceDelayRow(Gamma,m,delay);
          Q1 = AdvanceDelayRow(Q,m,delay);
          
          c = cos(theta); s = sin(theta);
          Qmm = c*exp(sqrt(-1)*alpha);        Qmn = s*exp(sqrt(-1)*phi); 
          Qnm = -s*exp(-sqrt(-1)*phi);        Qnn = c*exp(-sqrt(-1)*alpha);
          
          R = Gamma;
          R(n,:,:) = Qmm*Gamma(n,:,:) + Qmn*Gamma(m,:,:);
          R(m,:,:) = Qnm*Gamma(n,:,:) + Qnn*Gamma(m,:,:);
          
          Q = Q1;
          Q(n,:,:) = Qmm*Q1(n,:,:) + Qmn*Q1(m,:,:);
          Q(m,:,:) = Qnm*Q1(n,:,:) + Qnn*Q1(m,:,:);
          
          R = AdvanceDelayRow(R,m,-delay);
          Q = AdvanceDelayRow(Q,m,-delay);
          
       else
          stopcrit = 1;
          R = Gamma;
       end;

       if (iter >= maxiter),
          stopcrit = 1;
       end;

       % Truncation step
       [R,PHref] = TruncatePH_QR_2(R,PHMu,PHref);

       % Remove any unnecessary zeros from PU matrix
       [Q,PUref] = TruncatePU_2(Q,PUMu,PUref);

%        O(totaliter+1,1) = size(Q,3);
%        O(totaliter+1,2) = size(R,3);

       Gamma = R;
    end
    
    iter = 0;
    stopcrit = 0;
    
end





