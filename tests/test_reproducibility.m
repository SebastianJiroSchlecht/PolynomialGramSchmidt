function test_reproducibility
SetPaths;
set(groot,'defaultFigureVisible','off');
cleanup = onCleanup(@() close('all'));
rng(1);
A = randn(2,3,4)+1i*randn(2,3,4);
B = randn(3,2,3)+1i*randn(3,2,3);
assert(PolyMatNorm(PolyMatConvFast(A,B)-PolyMatConv(A,B))<1e-24);
assert(PolyMatNorm(PolyMatConvPHFast(A,A)-PolyMatConv(ParaHerm(A),A))<1e-24);
assert(PolyMatNorm(PolyMatConvPH2Fast(A,A)-PolyMatConv(A,ParaHerm(A)))<1e-24);

% Odd-length circular trimming must keep both boundary coefficients.
X = reshape([10 0 0 0 9],1,1,5);
Y = PolyMatTrim(X,1e-12);
assert(size(Y,3)==2 && abs(PolyMatNorm(Y)-181)<1e-12);
assert(~any(PolyMatTrim(zeros(2,2,3),1e-12),'all'));

% Both positive and negative phase jumps must be corrected (Algorithm 1).
af = [1 exp(-3i) exp(-3i) 1; zeros(1,4)];
qf = PolyVecDFTNormalisation(af,'off');
assert(all(real(sum(conj(qf(:,1:end-1)).*qf(:,2:end),1))>=0));

% Constant single-channel normalization and omitted smoothing argument.
q = PolyVecNormalisation(reshape([2 0],1,1,2),1e-8,16);
assert(abs(PolyMatNorm(q)-1)<1e-10);
try
    PolyVecDFTNormalisation(zeros(2,8));
    error('Zero input should have been rejected.');
catch exception
    assert(strcmp(exception.identifier,'PolyVecDFTNormalisation:InvalidVector'));
end

[A,Q,R] = PQRDToyExample(0);
assert(PolyMatNorm(A-PolyMatConv(Q,R))<1e-24);
Qhat = PolyMatGramSchmidt(A,1e-12,'off',7);
E = PolyMatConvPHFast(Qhat,Qhat);
E(:,:,(size(E,3)+1)/2) = E(:,:,(size(E,3)+1)/2)-eye(2);
assert(PolyMatNorm(E)<1e-20,'Toy example must yield a paraunitary factor.');

% Exercise every method through the paper's ensemble generator at two dimensions.
output = [tempname '.txt'];
outputCleanup = onCleanup(@() delete(output));
PQRDRandomMatrixEnsembleTest(output,[2 4],1,1:4);
r = readmatrix(output);
assert(size(r,1)==8 && size(unique(r(:,1:4),'rows'),1)==8);
assert(all(isfinite(r),'all') && all(r(:,5:end)>=0,'all'));
assert(isequal(sort(unique(r(:,3))),[2;4]));
assert(isequal(sort(unique(r(:,1))),[1;2;3;4]));
fprintf('Numerical regressions passed.\n');
end
