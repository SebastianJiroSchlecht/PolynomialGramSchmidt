%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE PQRD BY SEQUENTIAL BEST ROTATION ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code calculates a QR factorisation of a polynomial matrix.
% The algorithm consists of a left delay, followed by a left rotate.
% An optional left undelay after the rotation step can be applied - note that if an undelay
% is not used then convergence of the algorithm is more erratic. 
% To run with the undelay step set the parameter undelay equal to 1, otherwise
% set it equal to 0.
% The coefficient that is rotated at each iteration of the algorithm is chosen to be the
% largest off-diagonal coefficient from any of the polynomial elements beneath the
% diagonal of the matrix. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,zerolag,Q,zl_Q,NumIts,g,LTE] = PQRD_SBR(A,maxiter,undelay,crop_param,SCond,zerolag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    A            polynomial matrix 
%    maxiter      maximum number of iterations
%                 (default 1000)  
%    undelay      if set equal to 1 (as default), then an additional 
%                 undelay is applied, if set equal to zero no undelay 
%                 is applied. This is different to the PQRD-BS and -BC
%                 algorithms where an undelay is not required. 
%                 at each iteration
%    crop_param   maximum fraction of initial squared F-norm which may be lost due
%                 due to trimming the order of Dz after each iteration
%                 (default 0.0001)
%    SCond        stop if smallest coefficient associated with an element 
%                 beneath the diagonal is
%                 smaller than epsilon*reference
%                 (default 0.0001)  
%    zerolag      position of the zerolag coefficient matrix in Ainput
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    A            approximately upper-triangular polynomial matrix 
%    zerolag      position of the zerolag coefficient matrix in A
%    Q            paraunitary polynomial matrix
%    zl_Q         position of the zerolag coefficient matrix in Q
%    NumIts       Number of Iterations for the algorithm to converge
%    g            A vector containing the series of coefficients of the
%                 dominant coefficient (i.e. the one zeroed) over the
%                 series of iterations.
%    LTE          The sqaured F-norm of all elements beneath the diagonal of A
%                 at each iteration of the algorithm.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Default inputs if not specified:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if ~exist('maxiter'),
   maxiter = 1000;
end;   
if ~exist('undelay'),
   undelay = 1;        % This is different to the PQRD-BS and -BC algorithms where 
                       % the default setting is undelay = 0;
end;   
if ~exist('epsilon'),
   epsilon = 0.0001;
end;
if ~exist('mu'),
   mu = 0.0001;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

stopcrit = 0;
Parameters = [];
Mask = toeplitz([0; ones(size(A,1)-1,1)],zeros(1,size(A,2)));
maxiter = maxiter-1;
iter = 0;

while iter <= maxiter & stopcrit == 0
   
    [row,col,t,maxmag,maxval] = MaxElementBelowDiag(A);
    shift = t-zerolag;
    g(iter+1) = abs(A(row,col,t));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the energy below the diagonal at each iteration:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EnBD = zeros(size(A,1),size(A,2),size(A,3));
    for h = 1:(size(A,3))
        EnBD(:,:,h) = A(:,:,h).*Mask;   % R1 has zeros on and above the diagonal 
    end;
    LTE(iter+1) =  sum(sum(sum(EnBD.*conj(EnBD))));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if  maxmag > SCond    % If dominant coefficient beneath diagonal of A > stopping 
                          % condition apply EPGR to drive it to zero
                          
        iter = iter + 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DELAY:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A = LeftShift(A,row,-shift);               
        if shift > 0 
           zerolag = zerolag + shift;
        end;        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ROTATE:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Akk = A(col,col,zerolag);
        Ajk = A(row,col,zerolag);
        rot = Ajk/abs(Ajk);
        if Akk ~= 0 
           alpha = (Akk/abs(Akk));
           theta = atan2(abs(Ajk),abs(Akk));
        elseif Akk == 0 
           theta = pi/2;
           alpha = 0;
        end;
        A = Rotate(A,row,col,theta,alpha,rot);        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Parameters(:,iter) = [row,col,shift,theta,alpha,rot].';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % UNDELAY:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
        if undelay == 1 
             A = LeftShift(A,row,shift);
             if shift < 0
               zerolag = zerolag - shift;
             end;
        end;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TRUNCATE A:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [A,zerolag] = TruncBoth(A,zerolag,(sum(sum(sum(A.*(conj(A)))))),crop_param/2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    else   
        stopcrit = 1;
    end;
end;                                            % End of while loop.
NumIts = iter;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the paraunitary matrix Q:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = eye(size(A,1)); zl_Q = 1;
[Q,zl_Q] = CalcQ(Parameters,undelay,crop_param,Q,zl_Q);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = Rotate(A,row_index,col_index,theta,alpha,rot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applies rotation matrix to A from the left hand side
c=cos(theta);
s=sin(theta);
Gkk = c*conj(alpha);                 % G(k,k) (1,1)   c*e(i*alpha)
Gjj = c*(alpha);                     % G(j,j) (2,2)   c*e(-i*alpha)
Gkj = s*conj(rot);                   % G(k,j) (1,2)   se^i*phi
Gjk = -s*(rot);                      % G(j,k) (2,1)   -se^-i*phi
F1 = A(row_index,:,:);
F2 = A(col_index,:,:);
A(row_index,:,:) = Gjj*F1 + Gjk*F2;
A(col_index,:,:) = Gkj*F1 + Gkk*F2;   
clear F1 F2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
