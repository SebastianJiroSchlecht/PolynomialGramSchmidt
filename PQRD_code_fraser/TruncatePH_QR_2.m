function [R,ref] = TruncatePH_QR_2(R,Mu,ref)
% R = TruncatePH_2(R,MuRatio)
% 
% Truncates a Para-Hermitian Matrix removing the end lags whose sum of their
% F-Norm^2 is less than or equal to the vaulue given for MuRatio. Where R
% is the PH Matrix and MuRatio is Mu*N4 for the whole PH matrix. Further
% details of the 'trim function' can be found in J. G. McWhirter, 
% P. D. Baxter, T. Cooper, S. Redif, and J. Foster, ?An EVD Algorithm for 
% Para-Hermitian Polynomial Matrices,? IEEE Transactions on Signal 
% Processing, vol. 55, no. 5, pp. 2158?2169, May 2007.
%
% F. Coutts, University of Strathclyde, 15/12/2015
%
% Last modified 15/12/2015

% [M,N] = size(R);
% L = N/M;

MuRatio = Mu*ref/2;

% N4Trunc = norm(R(:,1:M),'fro')^2;
N4Trunc1 = norm(R(:,:,1),'fro')^2;
N4Trunc2 = norm(R(:,:,end),'fro')^2;
N4Trunc_temp1 = N4Trunc1;
N4Trunc_temp2 = N4Trunc2;
Tlim = 1;
while N4Trunc_temp1 <= MuRatio && N4Trunc_temp2 <= MuRatio
    N4Trunc1 = N4Trunc_temp1;
    Tlim = Tlim + 1;
%     N4Trunc_temp = N4Trunc + norm(R(:,(Tlim-1)*M+1:Tlim*M),'fro')^2;
    N4Trunc_temp1 = N4Trunc1 + norm(R(:,:,Tlim),'fro')^2;
    N4Trunc_temp2 = N4Trunc2 + norm(R(:,:,end-Tlim+1),'fro')^2;
end
if Tlim>1
%     R = R(:,(Tlim-1)*M+1:end);
    R = R(:,:,Tlim:end-(Tlim-1));
    ref = ref - N4Trunc1 - N4Trunc2;
end