function B = PolyMatTrim(A,p);
% B = PolyMatTrim(A,p);

E = squeeze(sum(sum(abs(A).^2,1),2));
[Emax,EmaxIndex] = max(E);

% assume maximum is centered, decay towards beginning and end 
E1_zeroed = find(E>(Emax*p));
E1_length = E1_zeroed(end)-E1_zeroed(1);

% assume maximum is centered, decay towards beginning and end 
Shift = ceil(length(E)/2);
E2 = circshift(E,-Shift);
E2_zeroed = find(E2>(Emax*p));
E2_length = E2_zeroed(end)-E2_zeroed(1);

if E1_length<=E2_length,
   B = A(:,:,E1_zeroed(1):E1_zeroed(end));
else
   disp('consider non-causal');
   E1_zeroed = find(circshift(E,Shift)>(Emax*p));
   E1_length = E1_zeroed(end)-E1_zeroed(1);
   A = circshift(A,[0 0 -Shift]);  
   B = A(:,:,E1_zeroed(1):E1_zeroed(end));
end;
