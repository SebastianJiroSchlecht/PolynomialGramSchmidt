
function [A,zerolag,Parameters,NumIts,LTE] = PQR_cols(A,maxiter,crop_param,SCond,zerolag,undelay);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialising any variables and vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters = []; LTE = []; [a1,a2,a3] = size(A);
Mask = toeplitz([0; ones(a1-1,1)],zeros(1,a2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate how many columms with elements beneath the diagonal in the 
% polynomial matrix:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if a2 == a1
  col_indices = [1:(a2-1)];
elseif a1>a2
   col_indices = [1:a2];
else
    col_indices = [1:(a1-1)];
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform a series of EPGRs on each column in turn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for hh = 1:length(col_indices)
    
    col_index = col_indices(hh);
    iter = 0;
    stopcrit = 0;
    
    while iter <= maxiter & stopcrit == 0 
        
        iter = iter + 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find optimal coefficient below diagonal in the kth col:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        R1 = zeros(size(A,1),size(A,2),size(A,3));
        for h = 1:(size(A,3))
            R1(:,:,h) = A(:,:,h).*Mask;         % Elements below the diagonal only
        end;
        [max_values_each_lag,row_indices] = max(abs(R1(:,col_index,:)));
        [max_mag,lag_index]=(max(max(abs(R1(:,col_index,:)))));
        row_index = row_indices(lag_index);
        maxval = A(row_index,col_index,lag_index);
        shift = lag_index - zerolag;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        
        if abs(maxval) > SCond      % If dominant coefficient in col > stopping 
                                    % condition apply EPGR to drive it to zero

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % APPLY DELAYS:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            A = LeftShift(A,row_index,-shift);               
            if shift > 0
               zerolag = zerolag + shift;
            end;                                         
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % APPLY ROTATIONS:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            Ajk = A(row_index,col_index,zerolag);
            Akk = A(col_index,col_index,zerolag);                
            rot = exp(i*angle(Ajk));
            alpha = exp(i*angle(Akk));
            theta = atan2(abs(Ajk),abs(Akk));
            A = Rotate(A,row_index,col_index,theta,alpha,rot);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Parameters(:,size(Parameters,2)+1) = [row_index,col_index,shift,theta,alpha,rot].';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Undelay step: (optional)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if undelay == 1 
               A = LeftShift(A,row_index,shift);
               if shift < 0
                  zerolag = zerolag - shift;
               end;
            end;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TRUNCATE ORDER OF THE MATRICES (NB: Require crop_param/2 due to the way the functions have been written)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [A,zerolag] = TruncBoth(A,zerolag,(sum(sum(sum(A.*(conj(A)))))),crop_param/2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % CALCULATE ENERGY BELOW THE DIAGONAL
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            R1 = zeros(size(A,1),size(A,2),size(A,3));
            for h = 1:(size(A,3))
                R1(:,:,h) = A(:,:,h).*Mask;         % Elements below the diagonal only
            end;
            LTE(length(LTE)+1) = sum(sum(sum(R1.*conj(R1))));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%          [EnergyBD(lp+1),EnergyAD(lp+1),EnergyZerolag,EnergyDiag] = CalcEnergies(A,zerolag);  % Function to calculate different energy mesaures.
            
      else 
            stopcrit = 1;
        end;            % end of if loop (if maxmag > SCond)
    end;                % end of while loop.
    
    Iterations_Taken(hh) = (iter);

end;     % end of for loop (goes through the col_indices)

NumIts = sum(Iterations_Taken-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions used within this algorithm:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = LeftShift(F, K, Shift)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply delay matrix from left hand side

% Inputs:
%       F       the matrix that the shift is being applied to.
%       K       the col index that the shift is to be applied to. 
%       Shift   the delay that we want shift F by. 

[M,N,L] = size(F);
if Shift == 0
    G = F;
elseif Shift > 0
    G = zeros(M,N,L+Shift);
    G(:,:,1:L) = F;
    G(K,:,1+Shift:L+Shift) = F(K,:,:);
    G(K,:,1:Shift) = zeros(1,N,Shift);
elseif Shift < 0
    G = zeros(M,N,L-Shift);
    G(:,:,1-Shift:L-Shift) = F;
    G(K,:,1:L) = F(K,:,:);
    G(K,:,L+1:L-Shift) = zeros(1,N,-Shift);
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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