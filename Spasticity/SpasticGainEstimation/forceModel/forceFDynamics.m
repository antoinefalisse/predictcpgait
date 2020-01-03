% Model of muscle-tendon force feedback
%
% INPUTS
%   eFf: muscle excitation from muscle-tendon force feedback (state)
%   Ftilde: muscle-tendon force (control) 
%   tau: time delay
%   g: feedback gain
%   b: parameter determining transition smoothness (tanh)
%   threshold: feedback threshold
%
% OUTPUT
%   deFfdt: derivative of muscle-tendon force

function deFfdt = forceFDynamics(eFf,Ftilde,tau,g,b,threshold)

c1 = -eFf./tau;
c2 = (g.*ones(size(c1,2),1))./tau;

f = 0.5*tanh(b.*(Ftilde-threshold));
deFfdt = c1+ c2.*Ftilde.*(f+0.5);

end

