% Explicit model of muscle-tendon force feedback
%
% INPUTS
%   t: time
%   Ftilde: muscle-tendon force
%   tau: time delay
%   g: feedback gain
%   b: parameter determining transition smoothness (tanh)
%   threshold: feedback threshold
%
% OUTPUT
%   eFf: muscle excitation from muscle-tendon force feedback

function eFf = forceFDynamicsExplicit(t,Ftilde,tau,g,b,threshold)

Ftilde_ppe = spline(t,Ftilde);

e0 = 0;

AbsTol = 1e-6;
RelTol = 1e-3;
options = odeset('AbsTol',AbsTol,'RelTol',RelTol);
ulf0 = e0;

[~,eFf] = ode45(@SpindleDynamicsFtOde,t,ulf0,options,Ftilde_ppe,tau,g,b,threshold);

return

function deFfdt = SpindleDynamicsFtOde(t,eFf,Ftilde_ppe,tau,g,b,threshold)

Ftilde = ppval(Ftilde_ppe,t);  

c1 = -eFf./tau;
c2 = (g.*ones(size(c1,2),1))./tau;

f = 0.5*tanh(b.*(Ftilde-threshold));
deFfdt = c1+ c2.*Ftilde.*(f+0.5);

return