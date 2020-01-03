% Explicit model of dF/dt feedback
%
% INPUTS
%   t: time
%   dFtilde: dF/dt
%   tau: time delay
%   g: feedback gain
%   b: parameter determining transition smoothness (tanh)
%   threshold: feedback threshold
%
% OUTPUT
%   edFf: muscle excitation from dF/dt feedback

function edFf = dFdtFDynamicsExplicit(t,dFtilde,tau,g,b,threshold)

dFtilde_ppe = spline(t,dFtilde);

e0 = 0;

AbsTol = 1e-6;
RelTol = 1e-3;
options = odeset('AbsTol',AbsTol,'RelTol',RelTol);
uvf0 = e0;

[~,edFf] = ode45(@SpindleDynamicsdFtOde,t,uvf0,options,dFtilde_ppe,tau,g,b,threshold);

return

function dedFfdt = SpindleDynamicsdFtOde(t,edFf,dFtilde_ppe,tau,g,b,threshold)

dFtilde = ppval(dFtilde_ppe,t);

c1 = -edFf./tau;
c2 = (g.*ones(size(c1,2),1))./tau;

f = 0.5*tanh(b.*(dFtilde-threshold));
dedFfdt = c1+ c2.*dFtilde.*(f+0.5);

return