function H = MIMOconv(F,G);
%  Y = MIMOconv(F,G);
%
%  Returns the convolution of two MIMO
%  system matrices F * G. This convolution operator
%  is NOT commutative, and if representig polynomial
%  matrices, F is left multiplied to G.
%
%  Input parameters
%      F      M x K x L1 MIMO system matrix 
%             M   output dimension
%             K   input dimension
%             L1   lenght of FIR filters
%      G      K x N x L1 MIMO system matrix 
%             K   output dimension
%             N   input dimension
%             L2   lenght of FIR filters
%
%  Output parameters
%      H      M x N x (L1+L2-1) MIMO system matrix
%  
%  S Weiss, Univ of Southampton, 15/7/2004
  
[M1,N1,L1] = size(F);
[M2,N2,L2] = size(G);
if N1 ~= M2,
  error('MIMOconv: dimensions do not agree');
end;
H = zeros(M1,N2,L1+L2-1);
for m = 1:M1,
  for n = 1:N2,
    for k = 1:N1,
      a = shiftdim(F(m,k,:),1);
      b = shiftdim(G(k,n,:),1);
      c(1,:,1) = conv(a,b);
      H(m,n,:) = H(m,n,:) + shiftdim(c,-1);
    end;   
  end;  
end;

