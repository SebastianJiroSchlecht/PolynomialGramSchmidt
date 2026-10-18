% QR Algorithm for Polynomial Matrices
%
% Algorithm consists of:-
%                           (1) inital phase adjustment
% Then at each iteration:
%                           (1) Locates max lower triangular element in kth column
%                           (2) Left delay
%                           (3) Left rotate 
%                               Uses the rotation matrix:
%                                              G = [cos(theta)e^(i*alpha)       sin(theta)e^i*phi]
%                                                  [-sin(theta)e^-i*phi     cos(theta)e^(-i*alpha)]
%                           (4) Optional Left Undelay (set the parameter undelay equal to 1, otherwise
%                                                      set it equal to 0.)
% The maximum lower triangular element in A(z) in chosen to be the element that
% is rotated. 

function [A,zerolag,Q,zl_Q,NumIts2,EBDAll,NumSweeps] = PQRDColsMultSweeps(A,maxiter,undelay,crop_param,SCond,zerolag,MaxSweeps);
    
% Input parameters:
%    A            polynomial matrix 
%    maxiter      maximum number of iterations
%                 (default 1000)  
%    undelay      if set equal to 1, then an additional undelay is applied,
%                 if set equal to zero (as default) no undelay is applied
%                 at each iteration
%    crop_param   maximum fraction of initial squared F-norm which may be lost due
%                 due to trimming the order of Dz after each iteration
%                 (default 0.0001)
%    SCond        stop if smallest coefficient associated with an element 
%                 beneath the diagonal is
%                 smaller than epsilon*reference
%                 (default 0.0001)  
%    zerolag      position of the zerolag coefficient matrix in A
%    MaxSweeps    maximum number of sweeps
%                 (default 1000)  

% Output parameters:
%    A            approximately upper-triangular polynomial matrix 
%    zerolag      position of the zerolag coefficient matrix in A
%    Q            paraunitary polynomial matrix
%    zl_Q         position of the zerolag coefficient matrix in Q
%    NumIts2      Number of Iterations for the algorithm to converge
%    EBDAll       Squared F-norm of elements beneath the diagonal of A over
%                 the series of iterations of the algorithm
%    NumSweeps    Number of sweeps required by the PQRD-BC algorithm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Default inputs if not specified:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
if ~exist('MaxSweeps'),
   MaxSweeps = 1000;
end;   
if ~exist('maxiter'),
   maxiter = 1000;
end;   
if ~exist('undelay'),
   undelay = 0;
end;   
if ~exist('epsilon'),
   epsilon = 0.0001;
end;
if ~exist('mu'),
   mu = 0.0001;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Initialise the paraunitary polynomial matrix:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
NumIts2 = 0; NumSweeps = 0; ParametersAll = []; 
MaskLT = toeplitz([0; ones(size(A,1)-1,1)],zeros(1,size(A,2)));
for k = 1:size(A,3)
    En(:,:,k) = A(:,:,k).*MaskLT;
end;
EBDAll = sum(sum(sum(En.*conj(En))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

for Sweep = 1:MaxSweeps
    
    [row1,col1,t1,maxmag1,DomCoef] = MaxElementBelowDiag(A);
    
    if abs(maxmag1) > SCond         % If magnitude of largest coefficient beneath the diagonal is greater
                                    % than stopping condition require another sweep
                                    % of the algorithm
        % Implement a single sweep of the PQRD-BC algorithm:
       [A,zerolag,Parameters,NumIts,LTE] = PQR_cols(A,maxiter,crop_param,SCond,zerolag,undelay);
       NumSweeps = NumSweeps + 1;
       NumIts2 = NumIts + NumIts2;
       ParametersAll = [ParametersAll,Parameters];
       EBDAll = [EBDAll,LTE];
    elseif maxmag1 <= SCond
        break
    end;
end;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the paraunitary matrix Q:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = eye(size(A,1)); zl_Q = 1;
[Q,zl_Q] = CalcQ(ParametersAll,undelay,crop_param,Q,zl_Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%