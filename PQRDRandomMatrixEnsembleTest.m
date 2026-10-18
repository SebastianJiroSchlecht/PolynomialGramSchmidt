function PQRDRandomMatrixEnsembleTest(FileName,SpaceDimM,SeedVals,Methods)
% Regenerate Section VI-C; keep the shipped archive unchanged.
% Optional arguments allow small validation runs with the same code path.
if nargin<1, FileName='PQRDEnsembleResults-recomputed.txt'; end
if nargin<2, SpaceDimM=[2 4 8 16 32]; end
if nargin<3, SeedVals=1:10000; end
if nargin<4, Methods=1:4; end
assert(~isfile(FileName),'Move the existing output before starting a fresh ensemble.');
SetPaths;
L = 2;
for dimensionIndex = 1:numel(SpaceDimM)
    M = SpaceDimM(dimensionIndex);
    for methodIndex = 1:numel(Methods)
        method = Methods(methodIndex);
        assert(ismember(method,1:4),'Unknown method.');
        if method==2 && M==32, continue; end
        seeds = SeedVals;
        if ismember(method,[2 4]), seeds=seeds(1:min(1000,numel(seeds))); end
        RESULTS = zeros(numel(seeds),10);
        valid = true(numel(seeds),1);
        for seedIndex = 1:numel(seeds)
            seed = seeds(seedIndex);
            Q = PUPolyMatRand(M,1,seed,'complex');
            randn('seed',seed);
            R = (randn(M,M,2)+1i*randn(M,M,2))/sqrt(2);
            for coefficientIndex = 1:2
                R(:,:,coefficientIndex) = triu(R(:,:,coefficientIndex));
            end
            % |det(Q)|=1: test the paper's determinant exclusion using R.
            Rf = fft(R,4096,3);
            determinant = ones(4096,1);
            for columnIndex = 1:M
                determinant = determinant.*squeeze(Rf(columnIndex,columnIndex,:));
            end
            if any(abs(determinant)<1e-12)
                valid(seedIndex)=false;
                warning('Excluded seed %d at M=%d: determinant below 1e-12.',seed,M);
                continue;
            end
            A = PolyMatConvFast(Q,R);
            if method==1
                result = ProposedQRDRandomMatrixSingleTest2(A,1e-8,2);
            elseif method==2
                result = PQRDBCRandomMatrixSingleTest2(A);
            elseif method==3
                result = PQRD_Faizan_CAMSAP2(A);
            else
                result = PQRD_Fraser_IMA2(A);
            end
            RESULTS(seedIndex,:) = [method seed M L result];
        end
        if isfile(FileName)
            dlmwrite(FileName,RESULTS(valid,:),'-append','precision','%.17g');
        else
            dlmwrite(FileName,RESULTS(valid,:),'precision','%.17g');
        end
    end
end
end
