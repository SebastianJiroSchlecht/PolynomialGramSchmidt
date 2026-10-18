function G = AdvanceDelayRow(F,Mrow,shift)
% G = AdvanceDelayRow(F,MRow,shift)
% 
% Advances or delays the M-th row of the polynomial matrix F by shift
% samples. A Positive shift gives an delay of the row.
%
% Last Update: J. Corr, University of Strathclyde, 17/1/14

% Obtain the dimensions of F
  [M,N,L] = size(F);
  % Preallocate the expanded sixe of F
  G = zeros(M,N,L+abs(shift));
  if(shift>0) % Delay
      G(:,:,1:L) = F;
      LT = L+shift;
      % Possibility to reduce to only shift non zero elements...
  
      % Copy the row plane that the max element is on and shift it back by
      % 'delay' to get the max elemant on the zerolag
      G(Mrow,:,shift+1:LT) = G(Mrow,:,1:LT-shift);
      % Complete the shift by zeroing the rows that were copied but not
      % overwriten - for preallocation this could be made more efficint so that
      % it does not zero sections that we know are zero...
      G(Mrow,:,1:shift) = zeros(1,N,shift);
  else % Advance
      shift = -shift;
      G(:,:,shift+1:L+shift) = F;
      LT=L+shift; 
      
      % Can be optimised to only shift non zero elements...
  
      % The column to be advanced is copied forward
      G(Mrow,:,1:LT-shift) = G(Mrow,:,shift+1:LT);
      % The elements not overwriten are now zerod
      G(Mrow,:,LT-(shift-1):LT) = zeros(1,N,shift);
  end