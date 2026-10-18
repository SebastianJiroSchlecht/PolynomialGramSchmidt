close all

K = 14;        %  order of paraunitary mixing matrix

maxiter = 50;
epsilon = 1e-6;
PHMu = 1e-6;
PUMu = 1e-6;
ensemble = 10^3;

% allM = [5 7 9 11 13];
allM = [5];
% allL = [5 7 9 11 13];
allL = [5];
numM = length(allM);

N = max(allL);
D = zeros(N*maxiter+1,ensemble,numM);
T = zeros(N*maxiter+1,ensemble,numM);
O = zeros(N*maxiter+1,ensemble,numM,2);
D2 = zeros(N*maxiter+1,ensemble,numM);
T2 = zeros(N*maxiter+1,ensemble,numM);
O2 = zeros(N*maxiter+1,ensemble,numM,2);

for iii = 1:ensemble
    Miter = 1;
    for M = allM
        L = allL(Miter);
        
        H = (randn(M,L,K+1)+sqrt(-1).*randn(M,L,K+1))./sqrt(2);

        [U2,R2,DM2_0,T2_0,O2_0] = PQRD_BC2_MS13_MultShift_Search(H,maxiter,epsilon,PHMu,PUMu);
        [U,R,DM_0,T_0,O_0] = PQRD_BC2(H,maxiter,epsilon,PHMu,PUMu);
    
        T(1:size(T_0,1),iii,Miter) = T_0; T2(1:size(T2_0,1),iii,Miter) = T2_0;
        D(1:size(DM_0,1),iii,Miter) = DM_0; D2(1:size(DM2_0,1),iii,Miter) = DM2_0;
        O(1:size(O_0,1),iii,Miter,:) = O_0; O2(1:size(O2_0,1),iii,Miter,:) = O2_0;
        Miter = Miter + 1;
    
    end
    iii
end



%%

meanEnsembleT = squeeze(mean(T,2));
meanEnsembleT2 = squeeze(mean(T2,2));
meanEnsembleD = squeeze(mean(D,2));
meanEnsembleD2 = squeeze(mean(D2,2));
meanEnsembleO = squeeze(mean(O,2));
meanEnsembleO2 = squeeze(mean(O2,2));

for iter = 1:numM
figure(200+iter);
plot(meanEnsembleT(meanEnsembleT(:,iter)>0,iter),10*log10(meanEnsembleD(meanEnsembleT(:,iter)>0,iter)))
hold all;
plot(meanEnsembleT2(meanEnsembleT2(:,iter)>0,iter),10*log10(meanEnsembleD2(meanEnsembleT2(:,iter)>0,iter)),'--r')
xlabel('Time / s','fontsize',14,...
  'Interpreter','Latex');
% ylabel('Ratio of below-diagonal to diagonal plus above-diagonal energy / [dB]','fontsize',14)
ylabel('$5\mathrm{log}_{10}\mathcal{E}\{E_{\mathrm{tri}}^{(i)} \}$ / [dB]','fontsize',14,...
  'Interpreter','Latex');
legend('standard','proposed');
set(legend(),'Interpreter','Latex')
t = sprintf('M = %d',allM(iter));
title(t);
grid on
hold off;
end

for iter = 1:numM
figure(400+iter);
plot(meanEnsembleT(meanEnsembleT(:,iter)>0,iter),(meanEnsembleO(meanEnsembleT(:,iter)>0,iter,1)))
hold all;
plot(meanEnsembleT2(meanEnsembleT2(:,iter)>0,iter),(meanEnsembleO2(meanEnsembleT2(:,iter)>0,iter,1)),'--r')
xlabel('Time / s','fontsize',14,...
  'Interpreter','Latex');
% ylabel('Ratio of below-diagonal to diagonal plus above-diagonal energy / [dB]','fontsize',14)
ylabel('Paraunitary order','fontsize',14,...
  'Interpreter','Latex');
legend('standard','proposed');
set(legend(),'Interpreter','Latex')
t = sprintf('M = %d',allM(iter));
title(t);
grid on
hold off;
end

