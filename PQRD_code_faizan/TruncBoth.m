
function [Atr,zerolag] = TruncBoth(A,zerolag,Aen,mu);

a = size(A,3);
A = TruncLeft(A,Aen,mu);  %NB: this function crops by mu/2
b = size(A,3);
zerolag = zerolag - (a-b);
Atr = TruncRight(A,Aen,mu);