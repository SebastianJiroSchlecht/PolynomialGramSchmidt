function C = PolyMatConvPH2Fast(A,B)
% C = PolyMatConvPH2Fast(A,B)
%
% PolyMatConvFast() is a fast convolution of two polynomial matrices A and B^P,
% noting the parahermitian transpose of the second argument.
% It replaces the command PolyMatConv(A,ParaHerm(B)), with a DFT domain 
% implementation of the convolution.
%
% Input parameters:
%    A           MxNxL1 input matrix
%    B           KxNxL2 input matrix
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
Cf = zeros(LA1,LB1,Nfft);
for k=1:LA3+LB3-1, 
   Cf(:,:,k) = Af(:,:,k)*Bf(:,:,k)';
end;
C = circshift(ifft(Cf,Nfft,3),[0 0 LB3-1]);
