% TSP26_PolyGramSchmidt_Figs56.m
%
% generate Figures 5 and 6 for Gram-Schmidt paper [1].
% 
% [1] Faizan A. Khattak, Ian K. Proudler, Stephan Weiss, and Sebastian J. Schlecht, 
%     "Analytic Gram-Schmidt Orthogonalisation and QR Decomposition for Polynomial 
%     Matrices," IEEE Transactions on Signal Processing, vol. 74, pp. 3597-3611, 
%     July 2026.

% S. Weiss, 23/7/2025

clear all; close all;

randn('seed',0);
M = 5;
A = randn(M,M,4) + 1i*randn(M,M,4);
tA = (0:3);
Aaxis = [-.5 3.5 -3 3];
%------------------------------------------------------------------------------
%  Figure 5 --- matrix A(z)
%------------------------------------------------------------------------------
figure(5);
for Nrow=1:M,
   for Ncol = 1:M,
      subplot(M,M,Ncol+M*(Nrow-1));
      stem(tA,squeeze(real(A(Nrow,Ncol,:))),'bo'); hold on; stem(tA,squeeze(imag(A(Nrow,Ncol,:))),'r*'); 
      axis(Aaxis); grid on; 
      if Nrow==1,
         title(sprintf('$m=%d$',Ncol),'interpreter','latex');
      end;
      if Ncol==1,      
         ylabel(sprintf('$A_{%d,m}[n]$',Nrow),'interpreter','latex');
      end;   
      if Nrow == M,
         xlabel('$n$','interpreter','latex');
      end;
      set(gca,'TickLabelInterpreter','latex');
   end;
end;     
set(gcf,'OuterPosition',[230 250 570 500]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Figure5.eps


%------------------------------------------------------------------------------
%  Analytic QRD and metrics
%------------------------------------------------------------------------------
if exist('Qhat.mat')==2,
   load Qhat.mat
else   
   Qhat = PolyMatGramSchmidt(A,1e-12,'on',4);
   save Qhat.mat Qhat
end;       
% error in paraunitarity
QE = PolyMatConv(ParaHerm(Qhat),Qhat);
MidIndex = (size(QE,3)+1)/2;
QE(:,:,MidIndex) = QE(:,:,MidIndex) - eye(M);
disp(sprintf('error in paraunitarity:   %f dB',10*log10(PolyMatNorm(QE))));
% triangularisation --- enery in upper and lower triangle
Rhat = PolyMatConv(ParaHerm(Qhat),A);
upper = 0; lower = 0;
for m = 1:M,
   lower = lower + PolyMatNorm(Rhat(m,1:m-1,:));
   upper = upper + PolyMatNorm(Rhat(m,m:M,:));
end;
TriEnergy = 10*log10(lower/PolyMatNorm(A));
disp(sprintf('error in triangularisation:   %f dB',TriEnergy));

% determine the support of Rhat down to -50~dB coefficient size
RR = zeros(M,M,50);
for m = 1:M,
    RR(m,:,:) = PolyMatCentreTrim(Rhat(m,:,:),50);
end;    
REn = squeeze(sum(sum(abs(RR).^2,1),2));

%------------------------------------------------------------------------------
%  display Q and Rhat
%------------------------------------------------------------------------------
Qh = zeros(M,M,21); Rh = zeros(M,M,15);
for m = 1:M,
   Qh(:,m,:) = PolyMatCentreTrim(Qhat(:,m,:),21);
   Rh(m,:,:) = PolyMatCentreTrim(Rhat(m,:,:),15);
end;   

%------------------------------------------------------------------------------
%  old Figure 6 --- matrix Q(z)
%------------------------------------------------------------------------------
figure(6);
tA = (0:20);
Aaxis = [0 20 -.5 .5];
for Nrow=1:M,
   for Ncol = 1:M,
      subplot(M,M,Ncol+M*(Nrow-1));
      stem(tA,squeeze(real(Qh(Nrow,Ncol,:))),'bo'); hold on; stem(tA,squeeze(imag(Qh(Nrow,Ncol,:))),'r*'); 
      axis(Aaxis); grid on; 
      if Nrow==1,
         title(sprintf('$m=%d$',Ncol),'interpreter','latex');
      end;
      if Ncol==1,      
         ylabel(sprintf('$Q_{%d,m}[n]$',Nrow),'interpreter','latex');
      end;   
      if Nrow == M,
         xlabel('$n$','interpreter','latex');
      end;
   end;
end;     
set(gcf,'OuterPosition',[230 250 570 500]);
set(gca,'LooseInset',get(gca,'TightInset'));

%------------------------------------------------------------------------------
%  Figure 6 --- matrix R(z)
%------------------------------------------------------------------------------
figure(6); clf;
tA = (0:14);
Rh(3,:,:) = -Rh(3,:,:);       % for display only ... effectively sign inversion of q2(z)
Aaxis = [0 14 -2 7];
for Nrow=1:M,
   for Ncol = 1:M,
      subplot(M,M,Ncol+M*(Nrow-1));
      stem(tA,squeeze(real(Rh(Nrow,Ncol,:))),'bo'); hold on; stem(tA,squeeze(imag(Rh(Nrow,Ncol,:))),'r*'); 
      axis(Aaxis); grid on; 
      if Nrow==1,
         title(sprintf('$m=%d$',Ncol),'interpreter','latex');
      end;
      if Ncol==1,      
         ylabel(sprintf('$\\hat{R}_{%d,m}[n]$',Nrow),'interpreter','latex');
      end;   
      if Nrow == M,
         xlabel('$n$','interpreter','latex');
      end;
      set(gca,'TickLabelInterpreter','latex','YTick',[-2 0 3 6]);
   end;
end;     
set(gcf,'OuterPosition',[230 250 570 500]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Figure6.eps
      

