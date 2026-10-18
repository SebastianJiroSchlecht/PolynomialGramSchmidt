function PlotPolyMat(H,zerolag,Hmin,Hmax);
  
% Plot the MIMO system described by H with zerolag coefficient matrix
% positioned at zerolag. 

% Input Parameters
%
% H - the polynomial matrix (dimensions p x q x T)
% zerolag - the position of the zerolag coefficient matrix (in the third
% dimension). I.e. The zerolag coefficient matrix is H(:,:,zerolag)
% Hmin - Sets the minimum y-axis value for each of the stem plots (overrides the default settings 
% of stem and can be useful if comparing two polynomial matrices). This is
% an optional input argument and does not need to be used if not required. 
% Hmax - Sets the maximum value of y-axis for each of the stem plots.
% Again, this is an optional input argument.

% Output Parameters
%
% A figure containing p x q stem plots represetning the series of
% coefficients associated with each of the polynomial elements of the
% matrix. The position of the stem plot relates to the position of the
% polynomial element that it represents. Each stem plot consist of T points.

% clf;
figure
[K,L,M] = size(H);

if ~exist('Hmin'),
   Hmin = min(min(min(H)));
end;  
if ~exist('Hmax'),
   Hmax = max(max(max(H)));
end;  

if Hmin>0, 
  Hmin=0; 
end;
Lx = size(H,3);
cc_u = (Lx - zerolag);
cc_l = zerolag - 1;
i = 0;
for k = 1:K,
  for l = 1:L,
    i = i + 1;
    subplot(K,L,i);
    h = stem(((-cc_l):cc_u),shiftdim(H(k,l,:),1),'k');
    hold on
    Q = [(-cc_l):cc_u];
    plot(Q,zeros(1,length(Q)),'k')
    axis([-cc_l cc_u Hmin*1.1 Hmax*1.1]);
    box on
  end;
end;
set(gcf,'Color','white')


