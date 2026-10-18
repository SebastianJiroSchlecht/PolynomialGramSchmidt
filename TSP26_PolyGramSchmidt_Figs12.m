% TSP26_PolyGramSchmidt_Figs12.m
%
% generate Figures 1 and 2 for Gram-Schmidt paper [1].
% 
% [1] Faizan A. Khattak, Ian K. Proudler, Stephan Weiss, and Sebastian J. Schlecht, 
%     "Analytic Gram-Schmidt Orthogonalisation and QR Decomposition for Polynomial 
%     Matrices," IEEE Transactions on Signal Processing, vol. 74, pp. 3597-3611, 
%     July 2026.

clear all; close all;

% vector from Example 2
a = zeros(2,1,3);
a(:,1,1) = [1; 0]; a(:,1,2) = [1; 1]; a(:,1,3) = [0; 1];
a = squeeze(a);

NDFT = 64;
A = fft(a,NDFT,2);
A(:,NDFT+1) = A(:,1);

% norm in bin at pi will be zero --- shift
OmegaPi = pi+2*pi/NDFT/8;
A(:,NDFT/2+1) = a*exp(-1i*OmegaPi*(0:2)');
Omegas = (0:NDFT-1)*2*pi/NDFT;
Omegas(NDFT/2+1) = OmegaPi;

NormA = zeros(NDFT+1,1);
U = zeros(2,NDFT+1); 
for k = 1:NDFT+1,
  NormA(k) = norm(A(:,k));
  U(:,k) = A(:,k)/NormA(k);
end; 

%-----------------------------------------------------------------------
% norm and non-uniform bins
%-----------------------------------------------------------------------
figure(1);
NormA_highres = sqrt(sum(abs(fft(a,1024,2)).^2));
NormA_highres(1025) = NormA_highres(1);
plot((0:NDFT)/NDFT,NormA_highres(1:16:end),'bo');
Omegas2 = Omegas; Omegas2(NDFT+1) = 2*pi;
hold on;
plot(Omegas2/2/pi,NormA,'r*');
plot((0:1024)/1024,NormA_highres,'b-');
plot((0:NDFT)/NDFT,NormA_highres(1:16:end),'bo');
Omegas2 = Omegas; Omegas2(NDFT+1) = 2*pi;
plot(Omegas2/2/pi,NormA,'r*');
grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/8:1),'XTickLabel',{'$0$','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$',...
      '$5\pi/4$','$3\pi/2$','$7\pi/4$','$2\pi$'});
xlabel('normalised angular frequency $\Omega$, $\Omega_k$','interpreter','latex');   
ylabel('$\|$\emph{\textbf{a}}$(\mathrm{e}^{\mathrm{j}\Omega})\|_2$, $\|$\emph{\textbf{a}}${}_k\|_2$','interpreter','latex');  
legend({'uniform','non-uniform'},'interpreter','latex','location','SouthWest');
Box1 = [15.375/32 16.625/32 0 0.2];
plot(Box1([1 1 2 2 1]),Box1([3 4 4 3 3]),'k-','linewidth',1);
axes('position',[.6 .5 .3 .35]);
  plot((0:NDFT)/NDFT,NormA_highres(1:16:end),'bo');
  hold on;
  plot(Omegas2/2/pi,NormA,'r*');
  plot((0:1024)/1024,NormA_highres,'b-');
  plot((0:NDFT)/NDFT,NormA_highres(1:16:end),'bo');
  plot(Omegas2/2/pi,NormA,'r*');
  grid on;
  axis(Box1);
  set(gca,'TickLabelInterpreter','latex',...
    'XTick',(31:33)/64,'XTickLabel',{'$\frac{31\pi}{32}$','$\pi$','$\frac{33\pi}{32}$'});
set(gcf,'OuterPosition',[230 250 570 250]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Figure1.eps

%-----------------------------------------------------------------------
% normalisation steps of Algorithm 1
%-----------------------------------------------------------------------
figure(2);
subplot(121);
plot3((0:NDFT)/NDFT,-real(U(1,:)),imag(U(1,:)),'bo-');
hold on; 
plot3((0:NDFT)/NDFT,-real(U(2,:)),imag(U(2,:)),'r*--');
axis([0 1 -1/sqrt(2) 1/sqrt(2) -1/sqrt(2) 1/sqrt(2)]); 
view(-70,15);
grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/4:1),'XTickLabel',{'$0$','$\pi/2$','$\pi$',...
      '$3\pi/2$','$2\pi$'},... 
    'YTick',[-1/sqrt(2) 0 1/sqrt(2)],'YTickLabel',...
     {'$\frac{1}{\sqrt{2}}$,','$0$','$-\frac{1}{\sqrt{2}}$'},...
    'ZTick',[-1/sqrt(2) 0 1/sqrt(2)],'ZTickLabel',...
     {'$-\frac{1}{\sqrt{2}}$,','$0$','$\frac{1}{\sqrt{2}}$'});
xlabel('$\Omega_k$','interpreter','latex');
ylabel('$\Re\{u_{m,k}\}$','interpreter','latex');
zlabel('$\Im\{u_{m,k}\}$','interpreter','latex');



subplot(122);
d = ones(1,NDFT+1); d(NDFT/2+1:NDFT) = -d(NDFT/2+1:NDFT);
d1 = d; d1(end) = -d1(end);
plot3((0:NDFT)/NDFT,-real(d1.*U(1,:)),imag(d1.*U(1,:)),'bo-');
hold on; 
plot3((0:NDFT)/NDFT,-real(d1.*U(2,:)),imag(d1.*U(2,:)),'r*--');
axis([0 1 -1/sqrt(2) 1/sqrt(2) -1/sqrt(2) 1/sqrt(2)]); 
view(-70,15);
grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/4:1),'XTickLabel',{'$0$','$\pi/2$','$\pi$',...
      '$3\pi/2$','$2\pi$'},... 
    'YTick',[-1/sqrt(2) 0 1/sqrt(2)],'YTickLabel',...
     {'$\frac{1}{\sqrt{2}}$,','$0$','-$\frac{1}{\sqrt{2}}$'},...
    'ZTick',[-1/sqrt(2) 0 1/sqrt(2)],'ZTickLabel',...
     {'$-\frac{1}{\sqrt{2}}$,','$0$','$\frac{1}{\sqrt{2}}$'});
xlabel('$\Omega_k$','interpreter','latex');
ylabel('$\Re\{d_k u_{m,k}\}$','interpreter','latex');
zlabel('$\Im\{d_k u_{m,k}\}$','interpreter','latex');

set(gcf,'OuterPosition',[230 250 570 310]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Figure2ab.eps

figure(3);
subplot(121);
d1=d;
d1(1:NDFT) = d1(1:NDFT).*exp(-1i*pi*(0:NDFT-1)/NDFT);
plot3((0:NDFT)/NDFT,-real(d1.*U(1,:)),imag(d1.*U(1,:)),'bo-');
hold on; 
plot3((0:NDFT)/NDFT,-real(d1.*U(2,:)),imag(d1.*U(2,:)),'r*--');
axis([0 1 -1/sqrt(2) 1/sqrt(2) -1/sqrt(2) 1/sqrt(2)]); 
view(-70,15);
grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/4:1),'XTickLabel',{'$0$','$\pi/2$','$\pi$',...
      '$3\pi/2$','$2\pi$'},... 
    'YTick',[-1/sqrt(2) 0 1/sqrt(2)],'YTickLabel',...
     {'$\frac{1}{\sqrt{2}}$,','$0$','$-\frac{1}{\sqrt{2}}$'},...
    'ZTick',[-1/sqrt(2) 0 1/sqrt(2)],'ZTickLabel',...
     {'$-\frac{1}{\sqrt{2}}$,','$0$','$\frac{1}{\sqrt{2}}$'});
xlabel('$\Omega_k$','interpreter','latex');
ylabel('$\Re\{q_{m,k}\}$','interpreter','latex');
zlabel('$\Im\{q_{m,k}\}$','interpreter','latex');

subplot(122);
d1=d;
d1(1:NDFT) = d1(1:NDFT).*exp(1i*pi*(0:NDFT-1)/NDFT);
plot3((0:NDFT)/NDFT,-real(d1.*U(1,:)),imag(d1.*U(1,:)),'bo-');
hold on; 
plot3((0:NDFT)/NDFT,-real(d1.*U(2,:)),imag(d1.*U(2,:)),'r*--');
grid on;
axis([0 1 -1/sqrt(2) 1/sqrt(2) -1/sqrt(2) 1/sqrt(2)]); 
view(-70,15);
grid on;
set(gca,'TickLabelInterpreter','latex',...
    'XTick',(0:1/4:1),'XTickLabel',{'$0$','$\pi/2$','$\pi$',...
      '$3\pi/2$','$2\pi$'},... 
    'YTick',[-1/sqrt(2) 0 1/sqrt(2)],'YTickLabel',...
     {'$\frac{1}{\sqrt{2}}$,','$0$','$-\frac{1}{\sqrt{2}}$'},...
    'ZTick',[-1/sqrt(2) 0 1/sqrt(2)],'ZTickLabel',...
     {'$-\frac{1}{\sqrt{2}}$,','$0$','$\frac{1}{\sqrt{2}}$'});
xlabel('$\Omega_k$','interpreter','latex');
ylabel('$\Re\{q_{m,k}\}$','interpreter','latex');
zlabel('$\Im\{q_{m,k}\}$','interpreter','latex');

set(gcf,'OuterPosition',[230 250 570 310]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Figure2cd.eps

