% Function to calculate the coefficient situated beneath the diagonal within
% a polynomial matrix with maximum magnitude. This coefficient is referred 
% to as the dominant coefficient. 

function [row,col,t,maxmag,maxval] = MaxElementBelowDiag(A);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A - polynomial matrix whose largest coefficient beneath the diagonal (in
% magnitude) is to be located
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT PARAMETERS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% row - row index of dominant coefficient (beneath diagonal)
% col - col index of dominant coefficient
% t - lag index of dominant coefficient
% Dominant coefficent is A(row,col,t)
% maxmag - magnitude of dominant coefficient
% maxval - Dominant coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a1,a2,a3] = size(A);
Mask = toeplitz([0; ones(a1-1,1)],zeros(1,a2)); % Matrix with zeros on and 
                                                % above the diagonal, zeros
                                                % everywhere else. 
for h = 1:(size(A,3))                  
    R(:,:,h) = A(:,:,h).*Mask;
end;

[a,b] = max(abs(R));    % Maximum element in each column of each slice (i.e. coefficient matrix) 
                        % of the polynomial matrix OD
                        % a contains the value of each maximum element, b
                        % contains the row position.
[c,d] = max(a);         % Maximum element in each vector a(:,:,k) c contains 
                        % the positions  
[maxmag,t] = max(c);
col = d(:,:,t);
row = b(:,col,t);
maxval = A(row,col,t);




