function [q,chi,chi_phi,n] = PhaseSmoothing_QR(q_in,P,max_iter);
%PhaseSmoothing(q_in,P);
%
%   With q_in an MxNfft matrix, this function aims to change the phase of every
%   column of q_in such that its inverse DFT yields a support that is as 
%   compact as possible. This is equivalent to smoothing the phase across the 
%   column vectors of q_in. 
%
%   This smoothing is performed by minimising the power in the P-the derivative 
%   of the Dirichlet interpolations through all M functions.   
%
%   Input parameters:
%      q_in       MxNfft matrix, representing a frequency-dependent M-element 
%                    vector
%      P          derivative order for power calculation
%      q_opt      optimal phase adjustment (optional & for comparison only)
%
%   Output parameters:
%      q          frequency-dependent vector with phase adjustment, such that 
%                 ifft(q,[],2) has maximally short support.

%   S. Weiss, UoS, 22/7/2020
%         updated with modulation correction, 31/7/2020
%         updated to account for correct Hessian, and working with only (Nfft-1) angles, 5/1/22

%-----------------------------------------
%  derivative and internal parameters
%-----------------------------------------
[M,Nfft] = size(q_in);

%-----------------------------------------
%  set up smoothness metric
%-----------------------------------------
c = (0:(Nfft-1)).^(P);
W = dftmtx(Nfft)/Nfft;
C = W*diag(c.^2)*W'; 
DISPLAY=0;
%-----------------------------------------
%  display interpolations of starting point
%-----------------------------------------


%-----------------------------------------
%  starting metrics
%-----------------------------------------

  chi_opt = 1.0;


%-----------------------------------------
%  initialisation
%-----------------------------------------
q_in = q_in*exp(-j*angle(q_in(1,1)));   % first reference element is real
q = zeros(M,Nfft);

%-----------------------------------------
%  iteration
%-----------------------------------------
phi2=0;
q = q_in;
prof = zeros(Nfft,1);
chi(1) = norm(ifft(q,Nfft,2)*diag(c),'fro').^2/chi_opt;
% if DISPLAY == 1, disp(sprintf('iteration 0, relative smoothness  %f',chi(1))); end;
D = zeros(Nfft,Nfft);
for m = 1:M, D = D + diag(q(m,:))'*C*diag(q(m,:)); end;
phi = zeros(Nfft,1)+pi/3;
mu = 2; crit1 = 1; crit2 = 1; n = 1; chi_old = 1/eps; Mod = 0;
% introduced a maximum for n into Faizan's modification
while (crit1 == 1)&(n<max_iter), 
  n = n + 1;
  % update as proposed in SSPD'20 paper
  Phi = diag(exp(j*phi));
  a = exp(j*phi);  
  Ar = [zeros(1,Nfft-1); diag(a(2:Nfft))];
  dummy = Ar'*D*a;
  H = 2*real(Ar'*D*Ar); % Hessian Matrix %- 2*diag(real(dummy));
  Hinv = pinv(H);
%  phi = phi - mu*Hinv*imag(conj(Phi)*D*Phi*ones(Nfft,1));
  phi(2:Nfft) = phi(2:Nfft) - mu*Hinv*imag(dummy);
%    phi = phi - mu*imag(conj(Phi)*D*Phi*ones(Nfft,1));
  phi = phi - phi(1);    % first reference element is real
  chi_phi=norm(phi-phi2);
  phi2=phi;



  for k = 1:Nfft,
    q(:,k) = q_in(:,k)*exp(j*phi(k));
  end;
  % the two following evaluations of the cost are equivalent
  %  chi(n) = ones(1,Nfft)*conj(Phi)*D*Phi*ones(Nfft,1);   
  chi(n) = norm(ifft(q,Nfft,2)*diag(c),'fro').^2/chi_opt;
  ratio = chi(n)/chi(n-1);
%   semilogy(n,chi(n),'k--o')
%   hold on
  % now need to check whether we have ended up in a critical point, and require modulation
  if (crit2 == 1) & (n<max_iter),
    if (ratio<1) & (ratio>0.99)   % learning curve has levelled out
      if (chi(n) < chi_old),   % smaller than at previous critical point -> modulation
        phi_old = phi;
        chi_old = chi(n);
        phi = phi+(2*pi/Nfft*(0:(Nfft-1))');
        Mod = 1;
      end;
    elseif (ratio>1000)&(Mod==1),
      % cost function is now worse than at the previous critical point
      %   [NEED TO FIND A MORE RELIABLE TEST FOR THE RATIO!]
        phi = phi_old;  % undo last modulation
        crit2 = 0;      % suppress further modulation steps
       
    else
      % do nothing
    end;
  else  % now iterate a bit longer until desired accuracy is reached
    if (ratio>0.99),    % changed from 0.9999
        if chi_phi<1e-4 crit1 = 0; end        % changed from 1e-5
    end;     % max number of iterations
  end;
end;
for k = 1:Nfft,
  q(:,k) = q_in(:,k)*exp(j*phi(k));
  
end;
%plot(phi)
