% This function computes muscle fiber lengths from muscle-tendon forces.
% More details in De Groote et al. (2016): DOI: 10.1007/s10439-016-1591-9
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function [lM,lMtilde ] = FiberLength_TendonForce(FTtilde,params,lMT)

% Ftilde = FT./params(1,:);
lMo = ones(size(FTtilde,1),1)*params(2,:);
lTs = ones(size(FTtilde,1),1)*params(3,:);
alphao = ones(size(FTtilde,1),1)*params(4,:);

% Non-linear tendon
lTtilde = (log(5*(FTtilde + 0.25))/35 + 0.995);

% Hill-model relationship
lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
lMtilde = lM./lMo;
end

