function Results = PQRD_Faizan_CAMSAP2(A);
% call Faizan's function for implementation of the CAMSAP'23 paper

%-------------------------------------------------------------------
%   generate matrix
%-------------------------------------------------------------------
[M,~,L] = size(A);

%-------------------------------------------------------------------
%   Faizan's
%-------------------------------------------------------------------
tic;
[Qhat,~,Rhat] = PQRD_analytic_Q_ver(A,1024,3);     % max DFT size 1024, differ. order for smoothing 4
%[Qhat,~,Rhat] = AnalyticQRnew(A,1024,3);     % max DFT size 1024, differ. order for smoothing 4
ExecTime = toc;

%Qhat = ParaHerm(Qhat);

% create a purely upper triangular version
Rphat = zeros(size(Rhat));            % R prime in paper
for n = 1:size(Rhat,3), Rphat(:,:,n) = triu(Rhat(:,:,n)); end;

%-------------------------------------------------------------------
%   metrics 
%-------------------------------------------------------------------
% loss of paraunitarity
QQ = PolyMatConv(ParaHerm(Qhat),Qhat);
QQ(:,:,(size(QQ,3)+1)/2) = QQ(:,:,(size(QQ,3)+1)/2) - eye(M);
etaPU = PolyMatNorm(QQ);
% upper trangularisation
etaTR = PolyMatNorm(Rhat-Rphat)/PolyMatNorm(Rphat);
% reconstruction error
Ahat = PolyMatConv(Qhat,Rphat);
[~,etaRC] = PolyMatAlign(A,Ahat);

% 60dB decay for Qhat and Rphat
Qord60 = size(PolyMatTrim(Qhat,1e-12),3);
Rord60 = size(PolyMatTrim(Rphat,1e-12),3);

% return 
Results = [etaPU etaTR etaRC Qord60 Rord60 ExecTime];

