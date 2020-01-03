% Model of sensory feedback
%
% INPUTS
%   a_sf: muscle activation from sensory feedback (state)
%   s: sensory information (control) 
%   tau: time delay
%   g: feedback gain
%   b: parameter determining transition smoothness (tanh)
%   threshold: feedback threshold
%
% OUTPUT
%   da_sfdt: derivative of muscle activation from sensory feedback

function da_sfdt = spindleDynamics(a_sf,s,tau,g,b,threshold)

c1 = -a_sf./tau;
c2 = (g.*ones(size(c1,2),1))./tau;

f = 0.5*tanh(b.*(s-threshold));
da_sfdt = c1+ c2.*s.*(f+0.5);

end




