%--------------------------------------------------------------------------
function dFTdt = Tendon_ForceLengthRelationshipODE(t,FT,ppa,pplMT,ppvMT,Parameters,Fvparam,Faparam,Kpe)
% Solve tendon force ODE for tendon force time history using numerical
% integration

% Calculate time-varying inputs at current time
a = ppval(ppa,t);
lMT = ppval(pplMT,t);
vMT = ppval(ppvMT,t);

FMax = Parameters(1);
lMopt = Parameters(2);
lTs = Parameters(3);
alphaopt = Parameters(4);
vMmax = Parameters(5);

% TENDON
FTtilde = FT/FMax;
lTtilde = real(log(5*(FTtilde + 0.25))/35 + 0.995);

% MUSCLE
% Geometric relationships
% Non Rigid Tendon: lT != lTs
w = lMopt*sin(alphaopt);
lMcosalpha = lMT-lTs*lTtilde;
lM = sqrt(w^2+lMcosalpha.^2);
alpha = asin(w/lM)';
lMtilde = lM/lMopt;

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

% Passive muscle force-length characteristic
e0 = 0.6; t50 = exp(Kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(Kpe); pp2 = (t7 - 0.10e1);
Fpparam = [pp1;pp2];

e0 = 0.6;
t5 = exp(Kpe * (lMtilde - 0.10e1) / e0);
FMpasltilde = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2);

% Active muscle force-velocity characteristic
FMce = FT./cos(alpha)-FMax*FMpasltilde;
FMvtilde = FMce./(a.*FMax.*FMltilde);

vMtilde = 1/Fvparam(2)*(sinh((FMvtilde-Fvparam(4))/Fvparam(1))-Fvparam(3));

vMtilde(vMtilde>1)=1;
vMtilde(vMtilde<-1)=-1;

vM = vMtilde*vMmax;
vT = vMT-vM./cos(alpha); 

dFTdt = FMax*7*exp(35*(lTtilde - 0.995)).*vT/lTs;

end