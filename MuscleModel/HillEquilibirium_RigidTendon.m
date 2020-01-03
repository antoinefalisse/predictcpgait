function [err, FT, FM, lMtilde, FMltilde, FMvtilde, FMpasltilde, cos_alpha] = HillEquilibirium_RigidTendon(a,lMT,vMT,Parameters,Fvparam,Faparam,Kpe)

Fmax = Parameters(1,:);
lMopt = Parameters(2,:);
lTs = Parameters(3,:);
alphaopt = Parameters(4,:);
vMmax = Parameters(5,:);

% Geometric relationships
% Rigid Tendon: lT = lTs
w = lMopt.*sin(alphaopt);
lM = sqrt((lMT-lTs).^2+w^2);

lMtilde = lM./lMopt;

lT = lTs*ones(size(lMtilde));
cos_alpha = (lMT-lT)./lM;
% TENDON 
% FTFL = TendonForceLengthRelationship(lTtilde);
% FT = Fmax*FTFL;

vMtilde = (vMT./vMmax).*cos_alpha;

% MUSCLE
% Active muscle force-length characteristic
b11 = Faparam(1);
b21 = Faparam(2);
b31 = Faparam(3);
b41 = Faparam(4);
b12 = Faparam(5);
b22 = Faparam(6);
b32 = Faparam(7);
b42 = Faparam(8);

b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
num3 = lMtilde-b23;
den3 = b33+b43*lMtilde;
FMtilde3 = b13*exp(-0.5*num3.^2./den3.^2);

num1 = lMtilde-b21;
den1 = b31+b41*lMtilde;
FMtilde1 = b11*exp(-0.5*num1.^2./den1.^2);

num2 = lMtilde-b22;
den2 = b32+b42*lMtilde;
FMtilde2 = b12*exp(-0.5*num2.^2./den2.^2);

FMltilde = FMtilde1+FMtilde2+FMtilde3;

% Active muscle force-velocity characteristic
FMvtilde = Fvparam(1)*log((Fvparam(2)*vMtilde+Fvparam(3))+sqrt((Fvparam(2)*vMtilde+Fvparam(3)).^2+1))+Fvparam(4);

FMacttilde = a.*FMltilde.*FMvtilde;

% Passive muscle force-length characteristic
e0 = 0.6; t50 = exp(Kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(Kpe); pp2 = (t7 - 0.10e1);
Fpparam = [pp1;pp2];

e0 = 0.6;
t5 = exp(Kpe * (lMtilde - 0.10e1) / e0);
FMpasltilde = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2);

FMtilde = (FMacttilde+FMpasltilde);
FM = Fmax*FMtilde;

FT = FM.*cos_alpha; % TBC
err =  FM.*cos_alpha-FT; % TBC

return