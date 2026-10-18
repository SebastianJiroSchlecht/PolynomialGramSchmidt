% TRUNCATION METHOD - TRUNCATE POLYNOMIAL MATRIX R(z) FROM RIGHT 
% USE WITH THE FUNCTION 'TruncLeft'

% This function truncates the polynomial matrix R(z) allowing at most c percent of N4 to be discarded
% from the matrix. N4 can be set as either the energy of the initial input matrix (N4(1)) or 
% alternatively the total energy within the matrix at that iteration, i.e.
% N4(iter).

%  Input parameters:
%
% R - polynomial matrix to be truncated.
% N4 - comparative measure to decide on how to truncate the polynomial
% matrix. Typically, set this measure equal to the squared Frobenius norm
% of the polynomial matrix R, i.e. sum(sum(sum(R.*conj(R)))).
% c - the proportion of N4 that is permitted to be truncated through a
% single application of this function.

%  Output parameters:
%
% Tr_R - the truncated version of the polynomial matrix R

function Tr_R = TruncRight(R,N4,c)

[a1,a2,a3] = size(R);
En = R.*(conj(R));
Q = 0;

for h = 1:a3
    p = [1:(a3-h)];
    Q = Q + sum(sum(En(:,:,(a3-(h-1)))));      
    if (Q/N4) > c   
      t1 = a3-(h-1);
      Tr_R = R(:,:,[1:t1]);
      return
    end;
end;