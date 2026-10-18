function error=paraunitarity_error(X)
[M,N,L]=size(X);

K=3*L;
Xf=fft(X,K,3);

if N==M
        Xf=mtimesx(Xf,permute(conj(Xf),[2,1,3]))-eye(N);
        error=sum(abs(Xf).^2,'all')/K;
else

    Xf=(vecnorm(Xf,2,1).^2-1).^2;
    error=sum(Xf,3)/K;
end

end
