%
% Ordered PQR (PQRD by Steps Algorithm)
% This code runs the PQR ordered as with the Givens technique for computing
% the QRD of a scalar matrix
%
% Decomposition performed is Q*A = R
%
% INPUT PARAMETERS:
% Ainput        Matrix to be decomposed
% maxiter       maximum number of iterations to allow each step to run for
%               (This is not the total amount of iterations for the algorithm!)
% crop_param    Truncation parameter (proportion of energy to be permitted to be
%               lost each time it the truncation method is used). This
%               stops the order of the matrices becoming unneccessarily
%               large. 
% SCond         Stopping condition for each step. Stops eahc step when all
%               coefficients from the polynomial element are less than
%               SCond.
% zerolag       Position of the zero-lag coefficient matrix in Ainput.
%
% OUTPUT PARAMETERS:
% A             Upper triangular matrix
% zerolag       Position of the zero-lag coefficient matrix
% Parameters    Matrix containing all the shift and rotation parameters
%               used to claculate A. This is needed for calculating Q
% iter          Number of iterations in this sweep
% g             Series of dominant coefficients at each iteration
% LTE           Squared Fnorm of elements beneath the diagonal of the
%               matrix A.

function [A,zerolag,Parameters,iter,g,LTE] = PQR_steps(A,maxiter,crop_param,SCond,zerolag,undelay)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialising any variables and vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter = 0; Parameters = [];
CM = 0;
Mask = toeplitz([0; ones(size(A,1)-1,1)],zeros(1,size(A,2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for row_index = 2:size(A,1)
        
    qq = min((row_index-1),size(A,2)); 
    for col_index = 1:qq
        stopcrit = 0;
        step_iter = 0;
        
        while step_iter <= maxiter & stopcrit == 0 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Finds optimal coefficient in A(row_index,col_index,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear R
            R1 = abs(A(row_index,col_index,:));
            [max_value,lag_index] = max(R1);
            shift = lag_index - zerolag;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if max_value > SCond  % If dominant coefficient in element > stopping condition 
                                  % apply EPGR to drive it to zero
                
                iter = iter + 1;
                step_iter = step_iter + 1;
                g(iter) = max_value;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Apply left shift:
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                A = LeftShift(A,row_index,-shift);            
                if shift > 0 
                   zerolag = zerolag + shift;
                end;
			    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Apply the Left Rotate Function:
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Ajk = A(row_index,col_index,zerolag);
                Akk = A(col_index,col_index,zerolag);
                rot = exp(i*angle(Ajk));
                alpha = exp(i*angle(Akk));
                theta = atan2(abs(Ajk),abs(Akk));
                A = Rotate(A,row_index,col_index,theta,alpha,rot);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                Parameters(:,iter) = [row_index,col_index,shift,theta,alpha,rot].';   % should i do conj of everything inside to counteract the '????
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Undelay step: (optional)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				if undelay == 1 
                   A = LeftShift(A,row_index,shift);
                   if shift < 0
                      zerolag = zerolag - shift;
                   end;
                end;
		        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Trim A matrices:
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [A,zerolag] = TruncBoth(A,zerolag,(sum(sum(sum(A.*(conj(A)))))),crop_param/2);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
                % CALCULATE ENERGY BELOW THE DIAGONAL:
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                clear Tt
                for tt = 1:size(A,3)
                    Tt(:,:,tt) = Mask.*A(:,:,tt);
                end;
                LTE(iter) = sum(sum(sum(Tt.*conj(Tt))));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            else 
                stopcrit = 1;
            end;                        % end of if max coefficient < SCond loop
        end;                            % End of while loop
    

    end;                % End of col index increment for loop
end;                    % End of row index increment for loop



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