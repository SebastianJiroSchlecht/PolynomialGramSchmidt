% This code calculates the paraconjugate of a polynomial matrix.

function Hpc = ParaConj(H)

[a,b,c] = size(H);

for t = 1:c
    Hpc(:,:,t) = H(:,:,c-(t-1))';
end;