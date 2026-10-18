% August 2007

% Checks the decomposition in the application results, checks the inverse
% decompsotion with the input matrix

% To do: Check also inverse decomposition when all elements beneath the
% diagonal are set equal to zero.

function [Erel,Ahat,zerolag_Ahat] = CheckInvDecPQR(A,Q,Ainput,zerolag,zl_Q,Ainput_zerolag);

Ahat = MIMOConv(ParaConj(Q),A);
min_lag_A = zerolag - 1;
max_lag_A = size(A,3) - zerolag;
min_lag_Q = size(Q,3) - zl_Q;           % This is the opposite way round the the A expressions as we are doing the paraconjgate(Q)*A
max_lag_Q = zl_Q - 1;
min_lag_Ahat = min_lag_Q + min_lag_A;
zerolag_Ahat = min_lag_Ahat + 1;
min_lag_Ainput = Ainput_zerolag - 1;
max_lag_Ainput = size(Ainput,3) - Ainput_zerolag;

e_lower = zerolag_Ahat- min_lag_Ainput;
e_upper = zerolag_Ahat + max_lag_Ainput;
Extract = Ahat(:,:,[e_lower:e_upper]);

% figure
% NpH_DispBlack(abs(Extract),1)

for t = 1:size(Ainput,3)
    Diff(:,:,t) = Ainput(:,:,t) - Extract(:,:,t);
end;

T1 = Ahat(:,:,[1:(e_lower-1)]);
T2 = Ahat(:,:,[(e_upper+1):size(Ahat,3)]);
D = sum(sum(sum(T1.*conj(T1)))) + sum(sum(sum(T2.*conj(T2))));
Erel = sqrt(sum(sum(sum(Diff.*conj(Diff)))) + D);
Erel = Erel/(sqrt(sum(sum(sum(Ainput.*conj(Ainput))))));

% fprintf(['Relative Error for the QRD:     ', num2str(Erel)])
% sprintf('\n')
