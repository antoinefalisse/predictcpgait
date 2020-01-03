% Model of dF/dt feedback
%
% INPUTS
%   edFf: muscle excitation from dF/dt feedback (state)
%   dFtilde: dF/dt (control) 
%   tau: time delay
%   g: feedback gain
%   b: parameter determining transition smoothness (tanh)
%   threshold: feedback threshold
%
% OUTPUT
%   dedFfdt: derivative of dF/dt

function dedFfdt = dFdtFDynamics(edFf,dFtilde,tau,g,b,threshold)

c1 = -edFf./tau;
c2 = (g.*ones(size(c1,2),1))./tau;

f = 0.5*tanh(b.*(dFtilde-threshold));
dedFfdt = c1+ c2.*dFtilde.*(f+0.5);

end




