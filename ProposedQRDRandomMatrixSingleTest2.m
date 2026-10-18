function Results = ProposedQRDRandomMatrixSingleTest2(A,threshold,P);
% Results = ProposedQRDRandomMatrixSingleTest(SeedVal,M,L,threshold,P)
% 
% Perform one simulation run within an ensemble of tests using randomised
% matrices with ground truth Procrustes solution. The randomisation is 
% initialised by the seed value SeedVal. The function returns a number of
% metrics for the simulations in [1] in the variable Results:
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
%       threshold     threshold for GS procedure
%       P             derivative order to phase smoothing operation inside GS
%
% Output parameter:
%       Results              vector containing various metrics
%
% [1] F.A. Khattak, I.K. Proudler, S.J. Schlecht. and S. Weiss: "Analytic 
%     Gram-Schmidt Orthogonalisation and QR Decomposition for Polynomial 
%     Matrices, to be submitted to IEEETrans. Signal Process., March 2025.

%-------------------------------------------------------------------
%   generate matrix
%-------------------------------------------------------------------
[M,~,L] = size(A);

%-------------------------------------------------------------------
%   calculate PQRD 
%-------------------------------------------------------------------
tic;
Qhat = PolyMatGramSchmidt(A,threshold,'on',P);
Rhat = PolyMatConvPHFast(Qhat,A);
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


