function lambda = DiagMeasUpperTri(R)
% lambda = DiagMeasUpperTri(R)
% 
% This function returns the upper-triangularisation measure of the
% matrix R. Ideally for a fully upper-triangular matrix lambda would be 0.
%
% F. Coutts, University of Strathclyde, 11/05/16

[M,N,L] = size(R);

R1 = R;
MPlus1 = M+1;
NSlice = M*N;
NT = NSlice*L;

for start = 1:M:NSlice
    for l = start:(MPlus1):NSlice
       R1(l:NSlice:NT)=0;
    end
end

R1N4 = CalcN4(R1);
N4 = CalcN4(R);
lambda = R1N4/N4;


end

