% test file for analytic PQRD algorithm
clear all
mult=3;
for L=[1:5]

for seed=1:100
randn('seed',seed)
R=randn(4,4,L*mult+1)+j*randn(4,4,L*mult+1);
R=R/sqrt(2);
R(2:4,1,:)=0;
R(3:4,2,:)=0;
R(4,3,:)=0;

P=4;
[~,Q,~,~,~]=random_STCM_SUM(seed,mult*L,4,0);
A=PolyMatConv(Q,R);
tic
[Q_out,chi,R_out] = PQRD_analytic_Q_ver(A,128,P);
t_ana(seed,L)=toc;
% order computation
ord_Q(seed,L)=size(Q_out,3)-1;
ord_R(seed,L)=size(R_out,3)-1;
% triangularization computation
t_ratio_prop(seed,L)=triangularization(A,R_out);

%% By columns
tic
[Q_BC,R_BC] = PQRD_BC2(A,500,1e-8,1e-6,1e-6);
t_BS(seed,L)=toc;
ord_Q_BC(seed,L)=size(Q_BC,3)-1;
ord_R_BC(seed,L)=size(R_BC,3)-1;
t_ratio_BC(seed,L)=triangularization(A,R_BC);

%% By sequential
tic
[Q_sec,R_sec,iter_sec(seed)] = PQRD_SM_1(A,500,1e-8,1e-6,1e-6);
t_seq(seed,L)=toc;
ord_Q_seq(seed,L)=size(Q_sec,3)-1;
ord_R_seq(seed,L)=size(R_sec,3)-1;
t_ratio_seq(seed,L)=triangularization(A,R_sec);

%%
% tic
% [Q_MS,R_MS] = PQRD_BC2_MS13_MultShift_Search(A,200,1e-8,1e-6,1e-6);
% t_MS(seed,L)=toc;
% ord_Q_MS(seed,L)=size(Q_MS,3)-1;
% ord_R_MS(seed,L)=size(R_MS,3)-1;
% t_ratio_MS(seed,L)=triangularization(A,R_MS);

seed
end
L
end


%% 


function error=reconst_tricky(A,Ar)
z=Ar-A(:,:,1);
for i=1:size(z,3)
z2(i)=norm(z(:,:,i));
end
[~,zlag]=min(z2);

L=max([size(Ar,3),size(A,3)])*2;
Af=fft(A,L,3);
Arf=ParaHermDFT(Ar,L,zlag);

Errord = 0;
ErrorNorm = 0;
for n = 1:L,
   Errord = Errord + norm(Arf(:,:,n)-Af(:,:,n),'fro')^2;
   ErrorNorm = ErrorNorm + norm(Af(:,:,n),'fro')^2; 
end;
%disp(sprintf('reconstruction error: %f',Error0));
error = Errord/ErrorNorm;

end


%% 

function tria=triangularization(A,Ar)
L=size(Ar,3);
for i=1:L
    Ar(:,:,i)=triu(Ar(:,:,i));
end
tria=1-FNorm3D(Ar)/FNorm3D(A);

end


