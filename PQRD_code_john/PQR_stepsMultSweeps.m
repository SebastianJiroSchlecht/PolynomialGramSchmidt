%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs PQRD by Steps algorithm, but allows for multiple sweeps if required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,zerolag,Q,zl_Q,NumIts,g1,LTE1,NumSweeps] = PQR_stepsMultSweeps(A,maxiter,undelay,crop_param,SCond,zerolag,MaxSweeps);
    
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
%    LTE1         Squared F-norm of elements beneath the diagonal of A over
%                 the series of iterations of the algorithm
%    NumSweeps    Number of sweeps required by the PQRD-BS algorithm

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

g1 = []; ParametersAll = []; NumIts = 0; NumSweeps = 0;
MaskLT = toeplitz([0; ones(size(A,1)-1,1)],zeros(1,size(A,2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Calculate energy beneath the diagonal of input matrix:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
for k = 1:size(A,3)
    En(:,:,k) = A(:,:,k).*MaskLT;
end;
LTE1 = sum(sum(sum(En.*conj(En))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

for Sweep = 1:MaxSweeps
    
    [row1,col1,t1,maxmag1,DomCoef] = MaxElementBelowDiag(A);
    
    if maxmag1 > SCond      % If magnitude of largest coefficient beneath the diagonal is greater
                            % than stopping condition require another sweep
                            % of the algorithm
       % Implement a single sweep of the PQRD-BS algorithm:                     
       [A,zerolag,Parameters,NumIts2,g,LTE] = PQR_steps(A,maxiter,crop_param,SCond,zerolag,undelay);
       NumSweeps = NumSweeps + 1;
       ParametersAll = [ParametersAll,Parameters];
       g1 = [g1,g]; LTE1 = [LTE1,LTE];
       NumIts = NumIts + NumIts2;
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