function Y = PolyMatCentreTrim(X,L);
% Y = PolyMatCentreTrim(X,N);
%
% For an MxNxK matrix X, PolyMatCentreTrim(X,N) extracts the central L 
% coefficients, returning an MxNxL matrix Y, such that energy in the 
% outer, trimmer coefficients is minimum.
% If N > L, Y is zero-padded.
%
% Input parameters:
%     X           MxNxK input matrix
%     L           output length
%
% Output parameter:
%     Y           MxNxL output matrix

% S. Weiss, 24/7/25

[M,N,K] = size(X);

if L==K,
   Y = X;
elseif L>K,
   Y = zeros(M,N,L); 
   Y(:,:,1:K) = X;
else
   Energy = squeeze(sum(sum(abs(X).^2,1),2));
   IndexStart = 1; IndexEnd = K;
   for k = 1:(K-L),
      if Energy(IndexStart)>Energy(IndexEnd),
         IndexEnd = IndexEnd-1;
      else
         IndexStart = IndexStart+1;
      end;
   end;
   Y = X(:,:,IndexStart:IndexEnd);      
end;      

