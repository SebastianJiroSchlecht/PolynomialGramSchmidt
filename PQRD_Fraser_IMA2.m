function Results = PQRD_Fraser_IMA2(A);
% calculate PQRD using Fraser Coutts' algorithm from IMA'16.

%-------------------------------------------------------------------
%   generate matrix
%-------------------------------------------------------------------
[M,~,L] = size(A);

%-------------------------------------------------------------------
%   Faizan's
%-------------------------------------------------------------------
tic;
[Qhat, Rhat] = PQRD_BC2_MS13_MultShift_Search(A,500,1e-8,1e-8,1e-8);
Qhat = ParaHerm(Qhat);
ExecTime = toc;
 
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
etaTR = PolyMatNorm(Rhat-Rphat)/PolyMatNorm(A);

% reconstruction error
Ahat = PolyMatConv(Qhat,Rphat);
[~,etaRC] = PolyMatAlign(A,Ahat);

% 60dB decay for Qhat and Rphat
Qord60 = size(PolyMatTrim(Qhat,1e-12),3);
Rord60 = size(PolyMatTrim(Rphat,1e-12),3);

% return 
Results = [etaPU etaTR etaRC Qord60 Rord60 ExecTime];

