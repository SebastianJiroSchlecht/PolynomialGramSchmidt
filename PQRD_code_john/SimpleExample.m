% Simple example to demonstrate the three different algorithms for 
% calculating the PQRD. 

clear 
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load polynomial matrix:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ainput = randn(5,3,3);
% load A1; Ainput = A1;                     % This is example 2 from Chapter 6 in my thesis
% load A2; Ainput = A2;                     % Example 3 from Chapter 6
% load ChanMatComplex; Ainput = C;          % Chapter 8, example 1
Ainput_zerolag = 1;
PlotPolyMat(abs(Ainput),Ainput_zerolag); [ax,h1]=suplabel('Lag');  [ax,h2]=suplabel('Magnitude','y'); [ax,h3]=suplabel('Input'  ,'t');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise variables:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mask = toeplitz([1; zeros(size(Ainput,1)-1,1)],ones(1,size(Ainput,2)));
SCond = 1e-3;
crop_param = 1e-6;
undelay = 0;
maxiter = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PQRD-BS (By Steps) algorithm:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[A1,zerolag1,Q1,zl_Q1,NumIts1,g1,LTE1,Sweeps1] = PQR_stepsMultSweeps(Ainput,maxiter,undelay,crop_param,SCond,Ainput_zerolag);
TIME_TAKEN_BS = toc;
PlotPolyMat(abs(Q1),zl_Q1);    [ax,h1]=suplabel('Lag');  [ax,h2]=suplabel('Magnitude','y'); [ax,h3]=suplabel('PQRD-BS: Q'  ,'t');
PlotPolyMat(abs(A1),zerolag1); [ax,h1]=suplabel('Lag');  [ax,h2]=suplabel('Magnitude','y'); [ax,h3]=suplabel('PQRD-BS: R'  ,'t');
% Find largest coefficient beneath the diagonal:
[row1,col1,t1,maxmag1,maxval1] = MaxElementBelowDiag(A1);
% Calculate the relative error:
[Erel1,Ahat1,zerolag_Ahat1] = CheckInvDecPQR(A1,Q1,Ainput,zerolag1,zl_Q1,Ainput_zerolag);
A1_new = zeros(size(A1,1),size(A1,2),size(A1,3));
for kk = 1:size(A1,3)
    A1_new(:,:,kk) = A1(:,:,kk).*Mask;
end;
[Erel1_all,Ahat1_all,zerolag_Ahat1_all] = CheckInvDecPQR(A1_new,Q1,Ainput,zerolag1,zl_Q1,Ainput_zerolag);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PQRD-BC (By Columns) algorithm:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
[A2,zerolag2,Q2,zl_Q2,NumIts2,LTE2,Sweeps2] = PQRDColsMultSweeps(Ainput,maxiter,undelay,crop_param,SCond,Ainput_zerolag);
TIME_TAKEN_BC = toc;
PlotPolyMat(abs(Q2),zl_Q2);    [ax,h1]=suplabel('Lag');  [ax,h2]=suplabel('Magnitude','y'); [ax,h3]=suplabel('PQRD-BC: Q'  ,'t');
PlotPolyMat(abs(A2),zerolag2); [ax,h1]=suplabel('Lag');  [ax,h2]=suplabel('Magnitude','y'); [ax,h3]=suplabel('PQRD-BC: R'  ,'t');
% Find largest coefficient beneath the diagonal:
[row2,col2,t2,maxmag2,maxval2] = MaxElementBelowDiag(A2);
% Calculate the relative error:
[Erel2,Ahat2,zerolag_Ahat2] = CheckInvDecPQR(A2,Q2,Ainput,zerolag2,zl_Q2,Ainput_zerolag);
A2_new = zeros(size(A2,1),size(A2,2),size(A2,3));
for kk = 1:size(A2,3)
    A2_new(:,:,kk) = A2(:,:,kk).*Mask;
end;
[Erel2_all,Ahat2_all,zerolag_Ahat2_all] = CheckInvDecPQR(A2_new,Q2,Ainput,zerolag2,zl_Q2,Ainput_zerolag);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PQRD-SBR (by Sequential Best Rotation) algorithm:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
undelay = 1;
tic
[A3,zerolag3,Q3,zl_Q3,NumIts3,g3,LTE3] = PQRD_SBR(Ainput,maxiter,undelay,crop_param,SCond,Ainput_zerolag);
TIME_TAKEN_SBR = toc;
PlotPolyMat(abs(Q3),zl_Q3);    [ax,h1]=suplabel('Lag');  [ax,h2]=suplabel('Magnitude','y'); [ax,h3]=suplabel('PQRD-SBR: Q'  ,'t'); 
PlotPolyMat(abs(A3),zerolag3); [ax,h1]=suplabel('Lag');  [ax,h2]=suplabel('Magnitude','y'); [ax,h3]=suplabel('PQRD-SBR: R'  ,'t');
% Find largest coefficient beneath the diagonal:
[row3,col3,t3,maxmag3,maxval3] = MaxElementBelowDiag(A3);
% Calculate the relative error:
[Erel3,Ahat3,zerolag_Ahat3] = CheckInvDecPQR(A3,Q3,Ainput,zerolag3,zl_Q3,Ainput_zerolag);
A3_new = zeros(size(A3,1),size(A3,2),size(A3,3));
for kk = 1:size(A3,3)
    A3_new(:,:,kk) = A3(:,:,kk).*Mask;
end;
[Erel3_all,Ahat3_all,zerolag_Ahat3_all] = CheckInvDecPQR(A3_new,Q3,Ainput,zerolag3,zl_Q3,Ainput_zerolag);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot squared F-norm of all elements beneath the diagonal of the
% polynomial matrix A over the series of iterations over each of the PQRD
% algorithms:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot([0:length(LTE1)-1],sqrt(LTE1),'k'); hold on
plot([0:length(LTE2)-1],sqrt(LTE2),'r:'); hold on
plot([0:length(LTE3)-1],sqrt(LTE3),'b:')
legend('PQRD-BS','PQRD-BC','PQRD-SBR'); set(gcf,'Color','white')
xlabel('Iteration'); ylabel('F-norm of elements below the diagonal of A')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparing the different methods:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'PQRD by....                           :      Steps            Cols         SBR \n')
fprintf(1,'Number of Iterations                  :   %10.0f   %10.0f   %10.0f\n', [NumIts1,NumIts2,NumIts3])
fprintf(1,'Number of Sweeps                      :   %10.0f   %10.0f         -\n', [Sweeps1,Sweeps2])
fprintf(1,'Fnorm of elements beneath the diagonal:   %10.6f   %10.6f   %10.6f\n',[sqrt(LTE1(end)),sqrt(LTE2(end)),sqrt(LTE3(end))])
fprintf(1,'Order of R                            :   %10.0f   %10.0f   %10.0f\n',[(size(A1,3)-1),(size(A2,3)-1),(size(A3,3)-1)])
fprintf(1,'Order of Q                            :   %10.0f   %10.0f   %10.0f\n',[(size(Q1,3)-1),(size(Q2,3)-1),(size(Q3,3)-1)])
fprintf(1,'Magnitude of Dominant Coefficient     :   %10.6f   %10.6f   %10.6f\n',[maxmag1,maxmag2,maxmag3])
fprintf(1,'Rel Error (Trun Only)                 :   %10.6f   %10.6f   %10.6f\n',[Erel1,Erel2,Erel3])
fprintf(1,'Rel Error (Trun and Approx)           :   %10.6f   %10.6f   %10.6f\n',[Erel1_all,Erel2_all,Erel3_all])
fprintf(1,'Computational Time (Secs)             :   %10.2f   %10.2f   %10.2f\n',[TIME_TAKEN_BS,TIME_TAKEN_BC,TIME_TAKEN_SBR])
fprintf(['\n'])
fprintf(['\n'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
