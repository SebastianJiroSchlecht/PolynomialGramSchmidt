function [Q,Gamma,iter] = PQRD_SM_1(R,maxiter,epsilon,PHMu,PUMu)
% PQRD method for Sequential PQRD proposed by Hassan Diyari

% Faizan Ahmad Khattak, 
% Centre for Signal and Image Processing
% University of Strathclyde, 11/05/2016


[M,N,L] = size(R);

if ~exist('maxiter'),
	maxiter = 400;
end;
if ~exist('epsilon'),
	epsilon = 0.0001;
end;
Q=eye(M);

[M,N,L] = size(R); % get channel dimensions
Q(:,:,1) = eye(M);


% Calculate N4 (the F-Norm of the entire system squared)
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
zerolag=1;
[Q,~]=qr(R(:,:,zerolag));

for i=1:size(R,3)
Gamma(:,:,i)=Q'*Gamma(:,:,i);
end
Q=Q';





    while stopcrit == 0,
       iter = iter+1;
       L=size(Gamma,3);
       %zerolag = (L+1)/2; 
       R1 = Gamma;
       NT = NSlice*L;

       for start = 1:M:NSlice
           for l = start:(MPlus1):NSlice
              R1(l:NSlice:NT)=0;
           end
       end

      % L2 norm
      [l,m,val]=max_row(abs(R1));
       if val>epsilon,
           
          delay = zerolag-l;
          if delay<0
          zerolag=zerolag+abs(delay);
          end
          
          
          Gamma = AdvanceDelayRow(Gamma,m,delay);
          Q = AdvanceDelayRow(Q,m,delay);
            
          [P,~]=qr(Gamma(:,:,zerolag));
          for i=1:size(Gamma,3)
          Gamma(:,:,i)=P'*Gamma(:,:,i);
          end
        
          for i=1:size(Q,3)
          Q(:,:,i)=P'*Q(:,:,i);
          end
          
          

    % Truncation step

         [Gamma,zerolag] = TruncBoth(Gamma,zerolag,PHref,PHMu);
%        [Gamma,PHref] = TruncatePH_QR_2(Gamma,1e-6,PHref);

       % Remove any unnecessary zeros from PU matrix
       [Q,PUref] = TruncatePU_2(Q,PUMu,PUref);

         
       else
          stopcrit = 1;
         
       end;

       if (iter >= maxiter),
          stopcrit = 1;
       end;

       % Truncation step
%        [Gamma,PHref] = TruncatePH_QR_2(Gamma,PHMu,PHref);

       % Remove any unnecessary zeros from PU matrix
%        [Q,PUref] = TruncatePU_2(Q,PUMu,PUref);

      
       
    end
    
    
    stopcrit = 0;
    


end

%%
function [tau,row,val]=max_row(R)
for m=2:size(R,1)
[v(m-1),l(m-1)]=max(vecnorm(R(m,:,:),2));
end
[val,row]=max(v);
tau=l(row);
row=row+1;


end


function G = AdvanceColumn(F,Mcol,Advance)
% G = AdvanceColumn(F,Mcol,Advance);
%
% Advances the M-th column of the polynomial matrix F by
% Advance samples.
%
%

% Obtain the dimensions of F
[M,N,L] = size(F);

TL=L+abs(Advance);
vec=1:N;
vec(vec==Mcol)=[];
G=zeros(M,N,TL);

if Advance>0
G(:,Mcol,Advance+1:TL)=F(:,Mcol,:);
G(:,vec,1:L)=F(:,vec,1:L);
else

G(:,Mcol,1:L)=F(:,Mcol,:);
G(:,vec,abs(Advance)+1:TL)=F(:,vec,:);


end
end


