% Calculates the paraunitary matrix Q from the PQRD by Steps algorithm

function [Q,zl_Q] = CalcQ(Parameters,undelay,crop_param,Q,zl_Q)

% Q = eye(M);
% zl_Q = 1;

for t = 1:size(Parameters,2)
    row_index = Parameters(1,t);
    col_index = Parameters(2,t);
    shift = Parameters(3,t);
    theta = Parameters(4,t);
    alpha = Parameters(5,t);
    rot = Parameters(6,t);
    
    Q = LeftShift(Q,row_index,-shift); 
    
    if shift > 0 
       zl_Q = zl_Q + shift;
    end;

    c=cos(theta);
    s=sin(theta);

    Gkk = c*conj(alpha);                 % G(k,k) (1,1)   c*e(i*alpha)
    Gjj = c*(alpha);                     % G(j,j) (2,2)   c*e(-i*alpha)
    Gkj = s*conj(rot);                   % G(k,j) (1,2)   se^i*phi
    Gjk = -s*(rot);                      % G(j,k) (2,1)   -se^-i*phi

    F1 = Q(row_index,:,:);
    F2 = Q(col_index,:,:);
    Q(row_index,:,:) = Gjj*F1 + Gjk*F2;
    Q(col_index,:,:) = Gkj*F1 + Gkk*F2;   
    clear F1 F2
   
    if undelay == 1 
       Q = LeftShift(Q,row_index,shift);
       if shift < 0
          zl_Q = zl_Q - shift;
       end;
    end;
    
    % TRUNCATE:
    en = sum(sum(sum(Q.*conj(Q))));
    q3 = size(Q,3);
    Q = TruncLeft(Q,en,crop_param/2);
    zl_Q = zl_Q - q3 + size(Q,3);
    Q = TruncRight(Q,en,crop_param/2);
    
end;
