% TSP26_PolyGramSchmidt_Fig7.m
%
% generate Figure 7 for Gram-Schmidt paper [1].
% 
% [1] Faizan A. Khattak, Ian K. Proudler, Stephan Weiss, and Sebastian J. Schlecht, 
%     "Analytic Gram-Schmidt Orthogonalisation and QR Decomposition for Polynomial 
%     Matrices," IEEE Transactions on Signal Processing, vol. 74, pp. 3597-3611, 
%     July 2026.

clear all; close all;

A = dlmread('PQRDEnsembleResults.txt');
% extract data for different methods
Results = zeros(5,5,4);
Upper = Results; Lower = Results;
M = 2.^(1:5);
for method = 1:4
   dummy = A(find(A(:,1)==method),3:end);
   % remaining data
   %   (1) M
   %   (2) L
   %   (3) paraunitary error
   %   (4) triangularisation error
   %   (5) reconstruction erroor
   %   (6) support of Q up to 60dB drop
   %   (7) support of R up to 60dB drop
   %   (8) time of execution
   
   for m = 1:length(M),
      if (m==5)&(method==2),
         % no result
      else        
         dummy2 = dummy(find(dummy(:,1)==M(m)),3:8);
         % add orders of Q and R together
         dummy2(:,4) = dummy2(:,4) + dummy2(:,5)-2;
         if or(method==1,method==3)  
            dummy2 = sort(dummy2(1:10000,[1 2 3 4 6]));               
            Results(m,:,method) = dummy2(5000,:);
            Lower(m,:,method) = dummy2(501,:);
            Upper(m,:,method) = dummy2(9500,:);
         else   
            dummy2 = sort(dummy2(1:1000,[1 2 3 4 6]));               
            Results(m,:,method) = dummy2(500,:);
            Lower(m,:,method) = dummy2(51,:);
            Upper(m,:,method) = dummy2(950,:);
         end;   
      end;
   end;
end;      

LW = 1;
for r = 1:3,
  subplot(1,5,r);
     if r==1,   % plot curves for legend only
         plot(-[1 1],[1 1],'d-.','color',[0 0.5 0]); hold on;
         plot(-[1 1],[1 1],'k:+');
         plot(-[1 1],[1 1],'r--*');
         plot(-[1 1],[1 1],'b-o');
     end;    
        
%  errorbar(0.98*M,Stats(8,:),(Stats(8,:)-Stats(7,:)),(Stats(9,:)-Stats(8,:)),'b.');
  errorbar(log2(M)-0.05,10*log10(Results(:,r,1)),...
           10*log10(Results(:,r,1))-10*log10(Lower(:,r,1)),...
           10*log10(Upper(:,r,1))-10*log10(Results(:,r,1)),'b.')
  hold on;         
  plot(log2(M),10*log10(Results(:,r,1)),'b-o','linewidth',LW); hold on;
  errorbar(log2(M)-0.05,10*log10(Results(:,r,2)),...
           10*log10(Results(:,r,2))-10*log10(Lower(:,r,2)),...
           10*log10(Upper(:,r,2))-10*log10(Results(:,r,2)),'.','color',[0 0.5 0])
  plot(log2(M),10*log10(Results(:,r,2)),'d-.','color',[0 0.5 0],'linewidth',LW);
  errorbar(log2(M)+0.05,10*log10(Results(:,r,3)),...
           10*log10(Results(:,r,3))-10*log10(Lower(:,r,3)),...
           10*log10(Upper(:,r,3))-10*log10(Results(:,r,3)),'r.')
  plot(log2(M),10*log10(Results(:,r,3)),'r--*','linewidth',LW);
  errorbar(log2(M)+0.05,10*log10(Results(:,r,4)),...
           10*log10(Results(:,r,4))-10*log10(Lower(:,r,4)),...
           10*log10(Upper(:,r,4))-10*log10(Results(:,r,4)),'k.')
  plot(log2(M),10*log10(Results(:,r,4)),'k:+','linewidth',LW);

  set(gca,'TickLabelInterpreter','latex','XTick',(1:5),'XTickLabel',...
   {'$2$','$4$','$8$','$16$','$32$'});
  axis([0.5 5.5 -350 0]); grid on;
  xlabel('$M$','interpreter','latex');
  if r==1, 
     text(1.2,-15,'(a)','interpreter','latex');
     ylabel('paraunitarity error, $10\log_{10}\{\eta_{\mathrm{PU}} \}$','interpreter','latex');
     legend({'PQRD-BC [17,20]','PQRD-MS [21]','PQRD-DFT [23]','proposed'},'interpreter','latex',...
             'location','West');
  elseif r==2,   
     text(1.2,-15,'(b)','interpreter','latex');
     ylabel('triangularisation error, $10\log_{10}\{\eta_{\mathrm{TR}} \}$','interpreter','latex');
  else,   
     text(1.2,-15,'(c)','interpreter','latex');
     ylabel('reconstruction error, $10\log_{10}\{\eta_{\mathrm{REC}} \}$','interpreter','latex');
  end;   
end;

subplot(1,5,5);
r = 5;
semilogy(log2(M),Results(:,r,1),'b-o','linewidth',LW); hold on;
errorbar(log2(M)-0.05,(Results(:,r,1)),...
           (Results(:,r,1))-(Lower(:,r,1)),...
           (Upper(:,r,1))-(Results(:,r,1)),'b.');
semilogy(log2(M),Results(:,r,2),'d-.','color',[0 0.5 0],'linewidth',LW); hold on;
errorbar(log2(M)-0.05,(Results(:,r,2)),...
           (Results(:,r,2))-(Lower(:,r,2)),...
           (Upper(:,r,2))-(Results(:,r,2)),'.','color',[0 0.5 0]);
semilogy(log2(M),Results(:,r,3),'r--*','linewidth',LW);
errorbar(log2(M)+0.05,(Results(:,r,3)),...
           (Results(:,r,3))-(Lower(:,r,3)),...
           (Upper(:,r,3))-(Results(:,r,3)),'r.');
semilogy(log2(M),Results(:,r,4),'k:+','linewidth',LW);
errorbar(log2(M)+0.05,(Results(:,r,4)),...
           (Upper(:,r,4))-(Results(:,r,4)),'k.');
semilogy(log2(M),Results(:,r,4),'k:');
  set(gca,'TickLabelInterpreter','latex','XTick',(1:5),'XTickLabel',...
   {'$2$','$4$','$8$','$16$','$32$'});
  xlabel('$M$','interpreter','latex'); 
  ylabel('execution time /[s]','interpreter','latex');
  axis([0.5 5.5 0.0005 300]); grid on;
text(1.2,160,'(e)','interpreter','latex');

subplot(1,5,4);
r = 4;
semilogy(log2(M),Results(:,r,1),'b-o','linewidth',LW); hold on;
errorbar(log2(M)-0.05,(Results(:,r,1)),...
           (Results(:,r,1))-(Lower(:,r,1)),...
           (Upper(:,r,1))-(Results(:,r,1)),'b.');
semilogy(log2(M),Results(:,r,2),'s-.','color',[0 0.5 0],'linewidth',LW); 
errorbar(log2(M)-0.05,(Results(:,r,2)),...
           (Results(:,r,2))-(Lower(:,r,2)),...
           (Upper(:,r,2))-(Results(:,r,2)),'.','color',[0 0.5 0]);
semilogy(log2(M),Results(:,r,3),'r--*','linewidth',LW);
errorbar(log2(M)+0.05,(Results(:,r,3)),...
           (Results(:,r,3))-(Lower(:,r,3)),...
           (Upper(:,r,3))-(Results(:,r,3)),'r.');
semilogy(log2(M),Results(:,r,4),'k:+','linewidth',LW);
errorbar(log2(M)+0.05,(Results(:,r,4)),...
           (Results(:,r,4))-(Lower(:,r,4)),...
           (Upper(:,r,4))-(Results(:,r,4)),'k.');
semilogy(log2(M),Results(:,r,4),'k:');
set(gca,'TickLabelInterpreter','latex','XTick',(1:5),'XTickLabel',...
   {'$2$','$4$','$8$','$16$','$32$'},...
   'YTick',[2 6 20 60 200 600 2000],'YTickLabel',...
   {'$2$','$6$','$20$','$60$','$200$','$600$','$2000$'})
xlabel('$M$','interpreter','latex'); 
ylabel('cumulative order, $\mathcal{O}\{\hat{\textbf{\textit{Q}}}(z)\} + \mathcal{O}\hat{\{\textbf{\textit{R}}}(z)\} $','interpreter','latex');
axis([0.5 5.5 1.5 2000]); grid on;
text(1.2,1450,'(d)','interpreter','latex');

set(gcf,'OuterPosition',[230 250 1500 600]);
set(gca,'LooseInset',get(gca,'TightInset'));
print -depsc Figure7.eps

