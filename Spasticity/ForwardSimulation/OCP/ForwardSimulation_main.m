% This function performs EMG-driven forward simulations using GPOPS

function [outputGPOPS,HilldiffGPOPS,M_musclesGPOPS,M_kneeGPOPS,FTGPOPS,...
    FTtildeGPOPS,dFtildeGPOPS,MAGPOS,lMtildeGPOPS,vMtildeGPOPS] = ...
    ForwardSimulation_main(time,EMG,parameters,lMT,MA,pathMuscleModel)

auxdata.NMuscles = size(EMG,2);
auxdata.tauAct = 0.015;
auxdata.tauDeact = 0.06;
auxdata.params = parameters;
auxdata.b = 0.1;

load([pathMuscleModel,'Fvparam.mat'],'Fvparam');
load([pathMuscleModel,'Faparam.mat'],'Faparam');

% Parameters of active muscle force-velocity characteristic
auxdata.Fvparam = Fvparam;
% Parameters of active muscle force-length characteristic
auxdata.Faparam = Faparam;
% Parameters of passive muscle force-length characteristic
e0 = 0.6; t50 = exp(4 * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(4); pp2 = (t7 - 0.10e1);
Fpparam = [pp1;pp2];
auxdata.Fpparam = Fpparam;

% Problem bounds
% Controls
e_min = 0; e_max = 1;       % bounds on muscle excitation
dF_min = -50;               % bounds on normalized derivative of force
dF_max = 50;      
% States
a_min = 0; a_max = 1;       % bounds on muscle activation
F_min = 0;                  % bounds on normalized force
F_max = 5;  
% Time bounds
t0 = time(1); tf = time(end);
bounds.phase.initialtime.lower = t0; 
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf; 
bounds.phase.finaltime.upper = tf;
% Controls bounds
umin = e_min*ones(1,auxdata.NMuscles); 
umax = e_max*ones(1,auxdata.NMuscles);
dFtildeMin = dF_min*ones(1,auxdata.NMuscles); 
dFtildeMax = dF_max*ones(1,auxdata.NMuscles);
bounds.phase.control.lower = [umin dFtildeMin]; 
bounds.phase.control.upper = [umax dFtildeMax];
% States bounds
actMin = a_min*ones(1,auxdata.NMuscles); 
actMax = a_max*ones(1,auxdata.NMuscles);
FtildeMin = F_min*ones(1,auxdata.NMuscles); 
FtildeMax = F_max*ones(1,auxdata.NMuscles);
bounds.phase.initialstate.lower = [actMin, FtildeMin,]; 
bounds.phase.initialstate.upper = [actMax, FtildeMax];
bounds.phase.state.lower = [actMin, FtildeMin]; 
bounds.phase.state.upper = [actMax, FtildeMax];
bounds.phase.finalstate.lower = [actMin, FtildeMin]; 
bounds.phase.finalstate.upper = [actMax, FtildeMax];
% Integral bounds
bounds.phase.integral.lower = 0; 
bounds.phase.integral.upper = 10000*(tf-t0);
% Path constraints
HillEquil = zeros(1,auxdata.NMuscles);
bounds.phase.path.lower = HillEquil; bounds.phase.path.upper = HillEquil;

% Initial guess (without damping)
N = size(time,1);
ind_muscles = 1:auxdata.NMuscles;
EMG_IG = EMG;
EMG_IG(EMG_IG<0.005)=0.005;
a_guess = zeros(size(time,1),auxdata.NMuscles);
for m = 1:auxdata.NMuscles
    a_guess(:,m) = integrateActivationDynamics(time,EMG_IG(:,m),...
        auxdata.tauAct,auxdata.tauDeact,auxdata.b);
end
[~,~,FT_IG] = generateMuscleMoments(time,auxdata.params,a_guess,lMT,MA,...
    auxdata.Fvparam,auxdata.Faparam,ind_muscles,[4 4 4]);
FTtilde_IG = FT_IG./repmat(auxdata.params(1,:),N,1);
LMT_IG = zeros(size(lMT,1),auxdata.NMuscles);
VMT_IG = zeros(size(lMT,1),auxdata.NMuscles);
for m = 1:auxdata.NMuscles
    ppy = spline(time,lMT(:,m));
    [LMT_IG(:,m),VMT_IG(:,m),~] = SplineEval_ppuval(ppy,time,1);
end
[dFtilde_IG] = getdFdt(FTtilde_IG,auxdata.params,a_guess,LMT_IG,VMT_IG,...
    auxdata.Fvparam,auxdata.Faparam,4);
dFtilde_IG(1,:) = dFtilde_IG(2,:);

guess.phase.time = time;
guess.phase.control = [EMG,dFtilde_IG];
guess.phase.state =  [a_guess,FTtilde_IG];
guess.phase.integral = 0;
for m = 1:auxdata.NMuscles
    auxdata.EMGSpline(m) = spline(time,EMG(:,m));
    auxdata.LMTSpline(m) = spline(time,lMT(:,m));     
    MASpline(m) = spline(time,MA(:,m));
end

% GPOPS setup        
setup.name = 'ForwardSimulation';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.derivatives.derivativelevel = 'second';
setup.nlp.ipoptoptions.tolerance = 1e-6;
setup.nlp.ipoptoptions.maxiterations = 10000;
setup.derivatives.supplier = 'adigator';
setup.scales.method = 'none';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-3;
setup.mesh.maxiterations = 0;
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.method = 'RPM-integration';
setup.displaylevel = 2;
Mesh_Frequency = 100;
NMeshIntervals = round((tf-t0)*Mesh_Frequency);
setup.mesh.phase.colpoints = 3*ones(1,NMeshIntervals);
setup.mesh.phase.fraction = (1/(NMeshIntervals))*ones(1,NMeshIntervals);
setup.functions.continuous = @ForwardSimulation_cont;
setup.functions.endpoint = @ForwardSimulation_end;

% ADiGator setup
global splinestruct
input.auxdata = auxdata;
tdummy = guess.phase.time;
splinestruct = ForwardSimulation_SplineInputData(tdummy,input);
splinenames = fieldnames(splinestruct);
for Scount = 1:length(splinenames)
  secdim = size(splinestruct.(splinenames{Scount}),2);
  splinestructad.(splinenames{Scount}) = ...
      adigatorCreateAuxInput([Inf,secdim]);
  splinestruct.(splinenames{Scount}) = zeros(0,secdim);
end
setup.auxdata.splinestruct = splinestructad;
adigatorGenFiles4gpops2(setup)

setup.functions.continuous = @Wrap4For_ForwardSimulation_cont;
setup.adigatorgrd.continuous = @ForwardSimulation_contGrdWrap;
setup.adigatorgrd.endpoint   = @ForwardSimulation_endADiGatorGrd;
setup.adigatorhes.continuous = @ForwardSimulation_contHesWrap;
setup.adigatorhes.endpoint   = @ForwardSimulation_endADiGatorHes;

outputGPOPS = gpops2(setup);
% Pre-allocation
lMTout = zeros(size(outputGPOPS.result.solution.phase.time,1),...
    auxdata.NMuscles);
vMTout = zeros(size(outputGPOPS.result.solution.phase.time,1),...
    auxdata.NMuscles);
MAout = zeros(size(outputGPOPS.result.solution.phase.time,1),...
    auxdata.NMuscles);
for m = 1:auxdata.NMuscles
     [lMTout(:,m),vMTout(:,m),~] = SplineEval_ppuval(...
         auxdata.LMTSpline(m),outputGPOPS.result.solution.phase.time,1);
     MAout(:,m) = ppval(MASpline(m),...
         outputGPOPS.result.solution.phase.time);
end
aGPOPS = outputGPOPS.result.solution.phase.state(:,1:auxdata.NMuscles);
FTtildeGPOPS = outputGPOPS.result.solution.phase.state(:,...
    auxdata.NMuscles+1:2*auxdata.NMuscles);
dFtildeGPOPS = 10*outputGPOPS.result.solution.phase.control(:,...
    auxdata.NMuscles+1:2*auxdata.NMuscles);
MAGPOS = MAout;
[HilldiffGPOPS, FTGPOPS] = HillEquilibrium_FTtildeState(aGPOPS,...
    FTtildeGPOPS,dFtildeGPOPS,lMTout,vMTout,auxdata.params,...
    auxdata.Fvparam,auxdata.Fpparam,auxdata.Faparam);
M_musclesGPOPS = MAout.*FTGPOPS;   
M_kneeGPOPS = sum(M_musclesGPOPS,2);
[~,lMtildeGPOPS] = getFiberLength(FTtildeGPOPS,auxdata.params,lMTout);
[~,vMtildeGPOPS] = getFiberVelocity(FTtildeGPOPS,dFtildeGPOPS,...
    auxdata.params,lMTout,vMTout);

delete ForwardSimulation_endADiGatorHes.mat
delete ForwardSimulation_endADiGatorHes.m
delete ForwardSimulation_endADiGatorGrd.mat
delete ForwardSimulation_endADiGatorGrd.m
delete ForwardSimulation_contADiGatorHes.mat
delete ForwardSimulation_contADiGatorHes.m
delete ForwardSimulation_contADiGatorGrd.mat
delete ForwardSimulation_contADiGatorGrd.m
delete ForwardSimulationIPOPTinfo.txt

end
