function [Q,Gamma,DM,T,O] = PQRD_BC2_MS13_MultShift_Search(R,maxiter,epsilon,PHMu,PUMu)

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
%
% Concatenates Q and applies as matrix multiplication. Concatenated row
% shifts.

[M,N,L] = size(R); % get channel dimensions

DM = zeros(N*maxiter+1,1);
T = zeros(N*maxiter+1,1);
O = zeros(N*maxiter+1,1);

if ~exist('maxiter'),
	maxiter = 400;
end;
if ~exist('epsilon'),
	epsilon = 0.0001;
end;

Q(:,:,1) = eye(M);

O(1,1) = 1;
O(1,2) = L;

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
totaliter = 0;
DM(1) = DiagMeasUpperTri(R);

% append zeros to the channel matrix so the zero lag is in the centre %
% Gamma = zeros(M,N,2*L-1);
% Gamma(:,:,L:end) = R;
Gamma = R;

totalmeas = 0;
tt = tic;
   
for n = 1:min(M-1,N), % loop over columns - starting from first through to last

    while stopcrit == 0,
       iter = iter+1;
       totaliter = totaliter+1;
       ttmeas = tic;
       DM(totaliter+1) = DiagMeasUpperTri(Gamma);
       totalmeas = totalmeas + toc(ttmeas);
       T(totaliter+1) = toc(tt)-totalmeas; 
       
       %Gamma = R;
       [~,~,L] = size(Gamma);

       % Zeolag needs re-aligned each iteration as R grows
       zerolag = ceil(L+1)/2; 

       R1 = Gamma;
       NT = NSlice*L;

       for start = 1:M:NSlice
           for l = start:(MPlus1):NSlice
              R1(l:NSlice:NT)=0;
           end
       end
       R1(:,[1:n-1 n+1:N],:) = 0;
       
       m = n+1:M;
       [maxEle,l] = max(abs(R1(m,n,:)),[],3);

       if n < M-1
           zerolag = round(max(zerolag,1));
           delays = [zeros(n,1); zerolag-l];
           Gamma = AdvanceDelayRow_MultShift(Gamma,delays);
           Q1 = AdvanceDelayRow_MultShift(Q,delays);   
           if min(delays) < 0
               newZL = zerolag - min(delays);
           else
               newZL = zerolag;
           end

           Q2 = eye(M,M);
           GZL = Gamma(:,:,newZL);
           for shift = M-n:-1:1
              a10 = GZL(n,n);
              a2t = GZL(m(shift),n);

              if abs(a10) > 0
                  theta = atan(abs(a2t)/abs(a10));
              else
                  theta = pi/2;
              end
              alpha = -angle(a10);
              phi = -angle(a2t);

              c = cos(theta); s = sin(theta);
              Qmm = c*exp(sqrt(-1)*alpha);        Qmn = s*exp(sqrt(-1)*phi); 
              Qnm = -s*exp(-sqrt(-1)*phi);        Qnn = c*exp(-sqrt(-1)*alpha);

              A = [GZL(n,:); GZL(m(shift),:)];
              GZL(n,:) = Qmm*A(1,:) + Qmn*A(2,:);
              GZL(m(shift),:) = Qnm*A(1,:) + Qnn*A(2,:);

              A = [Q2(n,:); Q2(m(shift),:)];
              Q2(n,:) = Qmm*A(1,:) + Qmn*A(2,:);
              Q2(m(shift),:) = Qnm*A(1,:) + Qnn*A(2,:);
           end

           L = size(Gamma,3);
           Lq = size(Q1,3);
           R = reshape(Q2*reshape(Gamma,M,L*N,1),M,N,L);
           Q = reshape(Q2*reshape(Q1,M,Lq*M,1),M,M,Lq);

           R = AdvanceDelayRow_MultShift(R,-delays);
           Q = AdvanceDelayRow_MultShift(Q,-delays);
       
       else
       
%        [maxEle,m,~,l] = MaxElement(abs(R1)); % Find location and value of max

       if maxEle>epsilon,
          zerolag = round(max(zerolag,1));
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
    
       end
%        else
%           stopcrit = 1;
%           R = Gamma;
%        end;

       if (iter >= maxiter/(M-n)),
%        if (iter >= maxiter),
          stopcrit = 1;
       end;

       % Truncation step
       [R,PHref] = TruncatePH_QR_2(R,PHMu,PHref);

       % Remove any unnecessary zeros from PU matrix
       [Q,PUref] = TruncatePU_2(Q,PUMu,PUref);

%        [~,~,Len] = size(Q);
%        O(totaliter+1) = Len;
%        [~,~,Len] = size(R);
%        O(totaliter+1) = Len;
       O(totaliter+1,1) = size(Q,3);
       O(totaliter+1,2) = size(R,3);

       Gamma = R;
%         figure(100); clf;
%         PolyMatDisplay(Gamma);
%         xx = 1;
    end
    
    iter = 0;
    stopcrit = 0;
    
end





