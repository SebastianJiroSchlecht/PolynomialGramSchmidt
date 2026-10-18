% TSP26_PolyGramSchmidt_Figs12.m
%
% generate Figures 3 and 4 for Gram-Schmidt paper [1].
% 
% [1] Faizan A. Khattak, Ian K. Proudler, Stephan Weiss, and Sebastian J. Schlecht, 
%     "Analytic Gram-Schmidt Orthogonalisation and QR Decomposition for Polynomial 
%     Matrices," IEEE Transactions on Signal Processing, vol. 74, pp. 3597-3611, 
%     July 2026.

% S. Weiss, 22/7/2025

clear all; close all;

[A,Q,R] = PQRDToyExample(0);

%------------------------------------------------------------------------------
%  Figure 3
%------------------------------------------------------------------------------
figure(3);
subplot(221);
stem(0:3,squeeze(real(A(1,1,:))),'bo'); hold on; stem(0:3,squeeze(imag(A(1,1,:))),'r*'); 
axis([-.5 3.5 -1 1]); grid on; 
title('$m=1$','interpreter','latex');
ylabel('$A_{1,m}[n]$','interpreter','latex');
%
subplot(222);
stem(0:3,squeeze(real(A(1,2,:))),'bo'); hold on; stem(0:3,squeeze(imag(A(1,2,:))),'r*'); 
axis([-.5 3.5 -1 1]); grid on; 
title('$m=2$','interpreter','latex');
%
subplot(223);
stem(0:3,squeeze(real(A(2,1,:))),'bo'); hold on; stem(0:3,squeeze(imag(A(2,1,:))),'r*'); 
axis([-.5 3.5 -1 1]); grid on;
ylabel('$A_{2,m}[n]$','interpreter','latex');
xlabel('time index $n$','interpreter','latex');
%
subplot(224);
stem(0:3,squeeze(real(A(2,2,:))),'bo'); hold on; stem(0:3,squeeze(imag(A(2,2,:))),'r*'); 
axis([-.5 3.5 -1 1]); grid on;
xlabel('time index $n$','interpreter','latex');
%
set(gcf,'OuterPosition',[230 250 570 320]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Figure3.eps

%------------------------------------------------------------------------------
%  Figure 4  norm of columns and determinant of A
%------------------------------------------------------------------------------
figure(4);
Af = fft(A,1024,3); Rf = fft(R,1024,3);
for k = 1:1024,
   detA(k) = det(Af(:,:,k));
   detR(k) = det(Rf(:,:,k));
   na1(k) = norm(Af(:,1,k)); na2(k) = norm(Af(:,2,k));
   nr1(k) = norm(Rf(:,1,k)); nr2(k) = norm(Rf(:,2,k));
end;   
f = (0:1023)/1024;
plot(f,na1,'b-'); hold on;
plot(f,na2,'r--')'
plot(f,abs(detA),'-.','color',[0 0.5 0]);
axis([0 1 0 3.5]); grid on;
xlabel('norm.~angular freq.~$\Omega$','interpreter','latex'); 
ylabel('norm $\|\mathbf{a}_m(\mathrm{e}^{\mathrm{j}\Omega})\|_2$','interpreter','latex');
legend({'$m=1$','$m=2$','$|\mathrm{det}\{\mathbf{A}(\mathrm{e}^{\mathrm{j}\Omega})|$'},...
'interpreter','latex','location','NorthEast');
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$',...
      '$3\pi/2$','$7\pi/4$','$2\pi$'},...
    'YTick',[0 1 2 3],'YTickLabel',...
     {'$0$','$1$','$2$','$3$'});
%set(gcf,'OuterPosition',[230 250 570 250]);
set(gcf,'OuterPosition',[230 250 570 280]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Figure4.eps

%------------------------------------------------------------------------------
%  PQRD Results --- proposed method
%------------------------------------------------------------------------------
disp('proposed method');
Qhat = PolyMatGramSchmidt(A,1e-12,'off',7);
Rhat = PolyMatTrim(PolyMatConv(ParaHerm(Qhat),A),1e-10);
% 
% error in paraunitarity
M=2;
QE = PolyMatConv(ParaHerm(Qhat),Qhat);
MidIndex = (size(QE,3)+1)/2;
QE(:,:,MidIndex) = QE(:,:,MidIndex) - eye(M);
disp(sprintf('error in paraunitarity:   %f dB',10*log10(PolyMatNorm(QE))));
% triangularisation --- enery in upper and lower triangle
upper = 0; lower = 0;
for m = 1:M,
   lower = lower + PolyMatNorm(Rhat(m,1:m-1,:));
   upper = upper + PolyMatNorm(Rhat(m,m:M,:));
end;
TriEnergy = 10*log10(lower/PolyMatNorm(A));
disp(sprintf('error in triangularisation:   %f dB',TriEnergy));
%
QhatQ = PolyMatConv(ParaHerm(Qhat),Q);
% in principle, this could be a diagonal matrix of allpass filters. Because of
% the simplicity of the setup, we have an ordinary matrix of phase shifts, 
% which we need to find
EnergyQ3 = squeeze(sum(sum(abs(QhatQ).^2,1),2));   % energy in 3rd dim.
[~,Q3Index] = max(EnergyQ3);
Phi = diag(diag(QhatQ(:,:,Q3Index)));
% display angles:
disp(sprintf('angle of phi (rad/pi): [%f, %f]',angle(diag(Phi))/pi));
% metric of QhatQ
QEabs = abs(QhatQ);
QEabs(:,:,Q3Index) = QEabs(:,:,Q3Index) - eye(2);
disp(sprintf('norm of QhatQ-Phi:   %f dB',10*log10(PolyMatNorm(QEabs))));
RE = Rhat - PolyMatConv(Phi,R);
disp(sprintf('error between R and Rhat:   %f dB',10*log10(PolyMatNorm(RE))));
% reconstruction error
Rhat2 = Rhat;
for k = 1:M,    Rhat2(m+1:M,m,:) = 0; end;
Ahat2 = PolyMatConv(Qhat,Rhat2);
[~,etaRC] = PolyMatAlign(A,Ahat2);
disp(sprintf('reconstruction error:         %f dB',10*log10(etaRC)));

%------------------------------------------------------------------------------
%  PQRD-BC Results [Foster 2010, Foster 2012]
%------------------------------------------------------------------------------
disp('-----');
disp('benchmark: [Foster 2010, Foster 2012]:');
path(path,'PQRD_code_john');
Mask = toeplitz([1; zeros(size(A,1)-1,1)],ones(1,size(A,2)));
SCond = 1e-3;
crop_param = 1e-6;
undelay = 0;
A_inputlag=1;
maxiter = 2000;
[A1,zerolag1,Q1,zl_Q1,NumIts1,g1,LTE1,Sweeps1] = PQR_stepsMultSweeps(A,maxiter,undelay,crop_param,SCond,A_inputlag);
% error in paraunitarity
M=2;
QE = PolyMatConv(ParaHerm(Q1),Q1);
MidIndex = (size(QE,3)+1)/2;
QE(:,:,MidIndex) = QE(:,:,MidIndex) - eye(M);
disp(sprintf('error in paraunitarity:   %f dB',10*log10(PolyMatNorm(QE))));
% triangularisation --- enery in upper and lower triangle
upper = 0; lower = 0;
for m = 1:M,
   lower = lower + PolyMatNorm(A1(m,1:m-1,:));
   upper = upper + PolyMatNorm(A1(m,m:M,:));
end;
TriEnergy = 10*log10(lower/PolyMatNorm(A));
disp(sprintf('error in triangularisation:   %f dB',TriEnergy));
% reconstruction error
for k = 1:size(A1,3),    A1(:,:,k) = A1(:,:,k).*Mask; end;
Ahat2 = PolyMatConv(ParaHerm(Q1),A1);
[~,etaRC] = PolyMatAlign(A,Ahat2);
disp(sprintf('reconstruction error:         %f dB',10*log10(etaRC)));
