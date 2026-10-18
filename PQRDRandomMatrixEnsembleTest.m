% PQRDRandomMatrixEnsembleTest.m
%
% Matlab script file to generate the ensemble simulation.
%
% [1] Faizan / Ian / Sebastian / Stephan, TSP


%-------------------------------------------------------------------------------
%  Parameters
%-------------------------------------------------------------------------------
FileName = 'PQRDEnsembleResults.txt';
% code for benchmark methods
%SetPaths;                                % set paths for PEVD toolbox files
path(path,'PQRD_code_john');              % benchmark algorithms
path(path,'PQRD_code_faizan');
path(path,'PQRD_code_fraser');

SpaceDimM = [2 4 8 16 32];                % spatial dimensions
Lengths = [2];                            % it says 'Lengths' but is the order of the polynomials
SeedVals = (1:10000);                     % seed values for random number generator
Methods = [1 2 3 4];                      % different methods
%   1 = proposed;  2 = Joanne Foster's; 3 = Faizan's CAMSAP; 4 = Fraser Coutts' 

%-------------------------------------------------------------------------------
%  Ensemble simulation
%-------------------------------------------------------------------------------
for i = 1:length(SpaceDimM),
   M = SpaceDimM(i);
   disp(sprintf('****** spatial dimension %d ******',M));
   for ell = 1:length(Lengths),
      L = Lengths(ell);
      disp(sprintf('***** polynomial length %d *****',L));
      for m = 1:length(Methods),
         method = Methods(m);
         disp(sprintf('*** Method %d *****',method));    
         RESULTS = zeros(1000,10);
         parfor i = 1:length(SeedVals),
            disp(sprintf('Seedvalue: %d',SeedVals(i)));
            %------------------------------------------------------------------- 
	    %  Generate random matrix
            %-------------------------------------------------------------------
            % paraunitary matrix
            Q = PUPolyMatRand(M,max(floor(L/2),1),SeedVals(i),'complex');
            randn('seed',SeedVals(i));
            % upper right triangular matrix
            R = (randn(M,M,floor(L/2)+1) + 1i*randn(M,M,floor(L/2)+1) )/sqrt(2);
            for m = 1:M-1, 
               for n = m+1:M,
                   R(n,m,:) = zeros(1,1,size(R,3));
               end;
            end;
            A = PolyMatConvFast(Q,R);       
            %------------------------------------------------------------------- 
	    %  call different PQRD algorithms
            %-------------------------------------------------------------------
            if method==1,              % proposed method
               Results = ProposedQRDRandomMatrixSingleTest2(A,1e-8,2);
            elseif method==2,   
               Results = PQRDBCRandomMatrixSingleTest2(A);
            elseif method==3,   
               Results = PQRD_Faizan_CAMSAP2(A);
            elseif method==4,   
               Results = PQRD_Fraser_IMA2(A);
            else
               disp('method not implemented');   
            end;
            %--------------------------------------------------------
            %  write results to a file
            %--------------------------------------------------------
%            if exist(FileName,'file') ~= 2,
%               disp('new results file created');
%               dlmwrite(FileName,[method SeedVals(i) M L Results]);
%            else  
%               disp('results appended');
%               dlmwrite(FileName,[method SeedVals(i) M L Results],'-append');
%            end;
             RESULTS(i,:) = [method SeedVals(i) M L Results];             
         end;      % loop for ensemble / seed values
         if exist(FileName,'file') ~= 2,
            disp('new results file created');
            dlmwrite(FileName,RESULTS);
         else   
            disp('results appended');
            dlmwrite(FileName,RESULTS,'-append');
         end;
      end;     % loop for methods      
   end;    % loop for temporal dimension L 
end;   % loop for spatial dimension M
