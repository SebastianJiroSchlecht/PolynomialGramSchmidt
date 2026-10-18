function G = AdvanceDelayRow_MultShift(F,shift)
% G = AdvanceDelayRow(F,MRow,shift)
% 
% Advances or delays the M-th row of the polynomial matrix F by shift
% samples. A Positive shift gives an delay of the row.
%
% Last Update: F. Coutts, University of Strathclyde, 16/05/16

% Obtain the dimensions of F
  [M,N,L] = size(F);
  
  maxDelay = max(shift);
  minDelay = min(shift);
  if maxDelay >= 0
      posMax = 1;
  else
      posMax = 0;
  end
  
  if minDelay >= 0
      posMin = 1;
  else
      posMin = 0;
  end
  
  %{
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
  %}
if posMax && posMin
    G = zeros(M,N,(L+maxDelay));
    G(:,:,1:L) = F;
    LT = (L+maxDelay);
    for ii = 1:length(shift)
      G(ii,:,shift(ii)+1:LT) = G(ii,:,1:LT-shift(ii));
      G(ii,:,1:shift(ii)) = zeros(1,N,shift(ii));
    end
elseif posMax && ~posMin
    G = zeros(M,N,(L+(maxDelay-minDelay)));
    G(:,:,(-minDelay)+1:L+(-minDelay)) = F;
    LT = (L+(maxDelay-minDelay));
    absShift = abs(shift);
    for ii = 1:length(shift)
      if shift(ii) >= 0
          G(ii,:,absShift(ii)+1:LT) = G(ii,:,1:LT-absShift(ii));
          G(ii,:,1:absShift(ii)) = zeros(1,N,absShift(ii));
      else
          G(ii,:,1:LT-absShift(ii)) = G(ii,:,absShift(ii)+1:LT);
          G(ii,:,LT-(absShift(ii)-1):LT) = zeros(1,N,absShift(ii));
      end
    end          
else
    G = zeros(M,N,(L+(-minDelay)));
    G(:,:,(-minDelay)+1:L+(-minDelay)) = F;
    LT = (L-minDelay);
    absShift = abs(shift);
    for ii = 1:length(shift)
      G(ii,:,1:LT-absShift(ii)) = G(ii,:,absShift(ii)+1:LT);
      G(ii,:,LT-(absShift(ii)-1):LT) = zeros(1,N,absShift(ii));
    end    
end  