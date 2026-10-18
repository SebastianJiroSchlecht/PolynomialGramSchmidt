function error=reconst_error_NPH(H,Hest)
N=max(size(H,3),size(Hest,3));
N=nextpow2(N);
Rf=fft(H,N,3);
Rfhat=fft(Hest,N,3);

Errord = 0;
ErrorNorm = 0;
for n = 1:N,
   Errord = Errord + norm(Rfhat(:,:,n)-Rf(:,:,n),'fro')^2;
   ErrorNorm = ErrorNorm + norm(Rf(:,:,n),'fro')^2; 
end;
%disp(sprintf('reconstruction error: %f',Error0));
error = Errord/ErrorNorm;







end