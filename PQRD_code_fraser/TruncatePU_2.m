function [A,ref] = TruncatePU_2(A,Mu,ref)
% A = TruncatePU_2d2(A,MuRatio)
% 
% Truncates a Para-Unitary Matrix removing the end lags whose sum of their
% F-Norm^2 is less than or equal to the vaulue given for MuRatio. Where A
% is the PU Matrix 
%
% F. Coutts, University of Strathclyde, 15/12/2015
%
% Last modified 15/12/2015

[M,N,L] = size(A);

N4Trunc = 0;
N4Trunc_temp = 0;
Tlim = 0;

% ref = (CalcN4(A));
% if abs(1-ref/2) > 0.0001
%     test = 1;
% end
MuRatio = Mu*ref/2;
% Trim the higher order coefficients
while N4Trunc_temp <= MuRatio
    N4Trunc = N4Trunc_temp;
%     N4Trunc_temp = N4Trunc + norm(A(:,end-M*(Tlim+1)+1:end-Tlim*M),'fro')^2;
    N4Trunc_temp = N4Trunc + norm(A(:,:,end-Tlim),'fro')^2;
    Tlim = Tlim + 1;
end
if Tlim>1
%     A = A(:,1:end-(Tlim-1)*M);
    A = A(:,:,1:end-(Tlim-1));
    ref = ref - N4Trunc;
end


% N4Trunc = norm(A(:,1:M),'fro')^2;
N4Trunc = norm(A(:,:,1),'fro')^2;
N4Trunc_temp = N4Trunc;
NTlim = 1;
% ALen = size(A,2)/size(A,1);
ALen = size(A,3);
% Trim the lower order coefficients
while (N4Trunc_temp <= MuRatio)&&(NTlim<ALen)
    N4Trunc = N4Trunc_temp;
    NTlim = NTlim + 1;
%     N4Trunc_temp = N4Trunc + norm(A(:,(NTlim-1)*M+1:NTlim*M),'fro')^2;
    N4Trunc_temp = N4Trunc + norm(A(:,:,NTlim),'fro')^2;
end
if NTlim>1
%     A = A(:,(NTlim-1)*M+1:end);
    A = A(:,:,NTlim:end);
    ref = ref - N4Trunc;
end

