function C = PolyMatConvPHFast(A,B)
% C = PolyMatConvPHFast(A,B)
%
% PolyMatConvFast() is a fast convolution of two polynomial matrices A^P and B,
% noting the parahermitian transpose of the first argument.
% It replaces the command PolyMatConv(ParaHerm(A),B), with a DFT domain 
% implementation of the convolution.
%
% Input parameters:
%    A           NxMxL1 input matrix
%    B           NxKxL2 input matrix
% Output parameter:
%    C           MxKx(L1+L2-1) product

% S. Weiss, UoS, 26/1/2026

%------------------------------------------------------------------------------
% parameters
%------------------------------------------------------------------------------
[LA1,LA2,LA3] = size(A); 
[LB1,LB2,LB3] = size(B);
Nfft = LA3+LB3-1;

%------------------------------------------------------------------------------
% convolution executed in the DFT domain
%------------------------------------------------------------------------------
Af = fft(A,Nfft,3);  Bf = fft(B,Nfft,3);
Cf = zeros(LA2,LB2,Nfft);
for k=1:LA3+LB3-1, 
   Cf(:,:,k) = Af(:,:,k)'*Bf(:,:,k);
end;
C = circshift(ifft(Cf,Nfft,3),[0 0 LA3-1]);
