function Results = PQRDBCRandomMatrixSingleTest2(A);
% Results = PQRDBCRandomMatrixSingleTest(SeedVal,M,L)
% 
% Perform one simulation run within an ensemble of tests using randomised
% matrices for a polynomial QRD using the iterative by-column approach in
% [2,3] as a benchmark for the results in [1].
%
%    Results(1)       eta_PU (loss of paraunitarity)
%    Results(2)       eta_TR (error in triangularisation)
%    Results(3)       eta_RC (reconstruction error)
%    Results(4)       Qord60     60dB decay order of Qhat
%    Results(5)       Rord60     60dB decay order of Rhat
%    Results(6)       execution time
% 
% Input parameter:
%       SeedVal       seed value for random number generator
%       M             spatial dimensions
%       L             temporal length (order plus one)
%
% Output parameter:
%       Results              vector containing various metrics
%
% [1] F.A. Khattak, I.K. Proudler, S.J. Schlecht. and S. Weiss: "Analytic 
%     Gram-Schmidt Orthogonalisation and QR Decomposition for Polynomial 
%     Matrices, to be submitted to IEEETrans. Signal Process., March 2025.
% [2] J. Foster, TSP, 2010
% [3] J. Foster, IET, 2012

%-------------------------------------------------------------------
%   generate matrix
%-------------------------------------------------------------------
[M,~,L] = size(A);

%-------------------------------------------------------------------
%   calculate PQRD-BC 
%-------------------------------------------------------------------
% parameters
Mask = toeplitz([1; zeros(size(A,1)-1,1)],ones(1,size(A,2)));
SCond = 1e-3;
crop_param = 1e-10;
undelay = 0;
maxiter = 500;
Ainput_zerolag = 1;
% PQRD-BC routine
tic;
[Rhat,zerolag2,Qhat,zl_Q2,NumIts2,LTE2,Sweeps2] = ...
     PQRDColsMultSweeps(A,maxiter,undelay,crop_param,SCond,Ainput_zerolag);
Qhat = ParaHerm(Qhat);
ExecTime = toc;

% create a purely upper triangular version
Rphat = zeros(size(Rhat));            % R prime in paper
for n = 1:size(Rhat,3), Rphat(:,:,n) = triu(Rhat(:,:,n)); end;

%-------------------------------------------------------------------
%   metrics 
%-------------------------------------------------------------------
% loss of paraunitarity
QQ = PolyMatConvPHFast(Qhat,Qhat);
QQ(:,:,(size(QQ,3)+1)/2) = QQ(:,:,(size(QQ,3)+1)/2) - eye(M);
etaPU = PolyMatNorm(QQ);

% upper trangularisation
etaTR = PolyMatNorm(Rhat-Rphat)/PolyMatNorm(A);

% reconstruction error
Ahat = PolyMatConvFast(Qhat,Rphat);
[~,etaRC] = PolyMatAlign(A,Ahat);

% 60dB decay for Qhat and Rphat
Qord60 = size(PolyMatTrim(Qhat,1e-12),3);
Rord60 = size(PolyMatTrim(Rphat,1e-12),3);

% return 
Results = [etaPU etaTR etaRC Qord60 Rord60 ExecTime];


