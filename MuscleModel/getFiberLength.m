% This function computes the muscle fiber length from the tendon force

function [lM,lMtilde ] = getFiberLength(FTtilde,params,lMT)

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

