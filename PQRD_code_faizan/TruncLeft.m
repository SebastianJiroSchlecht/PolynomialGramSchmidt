% TRUNCATION METHOD - TRUNCATE POLYNOMIAL MATRIX R(z) FROM LEFT
% USE WITH THE FUNCTION 'TruncRight'

% This function truncates the polynomial matrix R(z) allowing at most c percent of N4 to be discarded
% from the matrix. N4 can be set as either the energy of the initial input matrix (N4(1)) or 
% alternatively the total energy within the matrix at that iteration, i.e. N4(iter).

function Tr_R = TruncLeft(R,N4,c)

[a1,a2,a3] = size(R);
En = R.*(conj(R));
Q = 0;

for h = 1:a3       
    Q = Q + (sum(sum(sum(En(:,:,h)))));
  if (Q/N4) > c                                
      Tr_R = R(:,:,[h:a3]);
      return
    end;
end;