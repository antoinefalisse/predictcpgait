function [FM,Fce,Fpe,lMtilde]=HillModel_RigidTendon(a,lMT,vMT,alpha0,lTs,lMo,Fvparam,Fpparam,Faparam)

Nfr=size(a,1);
lMo = ones(Nfr,1)*lMo;
lTs = ones(Nfr,1)*lTs;
alphao = ones(Nfr,1)*alpha0;
vMmax = 10*lMo;

% Hill-type muscle model: geometric relationships
w = lMo.*sin(alphao);
lM = sqrt((lMT-lTs).^2+w.^2); % Rigid Tendon: lT = lTs
lMtilde=lM./lMo;
lT = lTs.*ones(size(lMtilde));
cos_alpha = (lMT-lT)./lM;

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
vMtilde = (vMT./vMmax).*cos_alpha;
e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);

FMvtilde = e1*log((e2*vMtilde+e3)+sqrt((e2*vMtilde+e3).^2+1))+e4;

% Active muscle force
% d is damping coefficient
d = 0.01;
Fce = a.*FMltilde.*FMvtilde +d*vMtilde;

% Passive muscle force-length characteristic
e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtilde - 0.10e1) / e0);
Fpe = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2);

FM=(Fce+Fpe).*cos_alpha;
Fce=Fce.*cos_alpha;
Fpe=Fpe.*cos_alpha;

end
