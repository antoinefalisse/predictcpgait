function [vM,vMtilde] = getFiberVelocity(FTtile,dFTtilde,params,lMT,vMT)


% Ftilde = FT./params(1,:);
lMo = ones(size(FTtile,1),1)*params(2,:);
lTs = ones(size(FTtile,1),1)*params(3,:);
alphao = ones(size(FTtile,1),1)*params(4,:);
vMmax = ones(size(FTtile,1),1)*params(5,:);
Atendon = 35;

% Inverse tendon force-length characteristic
lTtilde = log(5*(FTtile + 0.25))./Atendon + 0.995;

% Hill-type muscle model: geometric relationships
lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);

vT = lTs.*dFTtilde./(7*exp(35*(lTtilde-0.995)));
cos_alpha = (lMT-lTs.*lTtilde)./lM;
vM = (vMT-vT).*cos_alpha;
vMtilde = vM./vMmax;