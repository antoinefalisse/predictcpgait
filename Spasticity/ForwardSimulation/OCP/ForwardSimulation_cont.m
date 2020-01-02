function phaseout = ForwardSimulation_cont(input)

% Get constant input data
NMuscles        = input.auxdata.NMuscles;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
splinestruct    = input.auxdata.splinestruct;
Faparam         = input.auxdata.Faparam;
Fpparam         = input.auxdata.Fpparam;
Fvparam         = input.auxdata.Fvparam;
b               = input.auxdata.b;
EMG = splinestruct.EMG;
lMT = splinestruct.lMT;
vMT = splinestruct.vMT;

% Get optimization variables
numColPoints    = size(input.phase.state,1);
% Get controls
e   = input.phase.control(:,1:NMuscles);
dFtilde  = 10*input.phase.control(:,NMuscles+1:2*NMuscles);
% Get states
a       = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:2*NMuscles);

% DYNAMIC CONSTRAINTS
% Activation dynamics
dadt = ones(numColPoints,NMuscles);
for m = 1:NMuscles
    dadt(:,m) = ActivationDynamics(e(:,m),a(:,m),tauAct,tauDeact,b);
end
% Contraction dynamics is implicit
phaseout.dynamics = [dadt dFtilde];

% PATH CONSTRAINTS
% Hill-equilibrium
[Hilldiff, ~] = HillEquilibrium_FTtildeState(a,Ftilde,dFtilde,lMT,vMT,...
    params,Fvparam,Fpparam,Faparam);
phaseout.path = Hilldiff;

% OBJECTIVE FUNCTION
Ediff = EMG - e;
ediffNorm = sum((Ediff).^2,2);
phaseout.integrand = ediffNorm;
