% This script estimates the spastic gains of the force-related spasticity
% model based on sensory information from fast passive stretch motions

function [output,DatStore] = forceModel_main(IK_path,ID_path,...
    MuscleAnalysis_path,EMG_path_agonist,EMG_path_antagonist,Misc)

%% Filter parameters 
% Default low-pass filter: Butterworth order: 6, Cutoff frequency: 6Hz
% Inverse Dynamics
if ~isfield(Misc,'f_cutoff_ID') || isempty(Misc.f_cutoff_ID)
    Misc.f_cutoff_ID=6;
end
if ~isfield(Misc,'f_order_ID') || isempty(Misc.f_order_ID)
    Misc.f_order_ID=6;
end
% Muscle-tendon lengths
if ~isfield(Misc,'f_cutoff_lMT') || isempty(Misc.f_cutoff_lMT)
    Misc.f_cutoff_lMT=6;
end
if ~isfield(Misc,'f_order_lMT') || isempty(Misc.f_order_lMT)
    Misc.f_order_lMT=6;
end
% Moment arms
if ~isfield(Misc,'f_cutoff_dM') || isempty(Misc.f_cutoff_dM)
    Misc.f_cutoff_dM=6;
end
if ~isfield(Misc,'f_order_dM') || isempty(Misc.f_order_dM)
    Misc.f_order_dM=6;
end
% Inverse Kinematics
if ~isfield(Misc,'f_cutoff_IK') || isempty(Misc.f_cutoff_IK)
    Misc.f_cutoff_IK=6;
end
if ~isfield(Misc,'f_order_IK') || isempty(Misc.f_order_IK)
    Misc.f_order_IK=6;
end
% Mesh Frequency
if ~isfield(Misc,'Mesh_Frequency') || isempty(Misc.Mesh_Frequency)
    Misc.Mesh_Frequency=300;
end

%%  Get muscle information 
Nsegment = Misc.Nsegment;
% Constant parameters
switch Misc.joint
    case 'knee'
        Muscles_spas = [1 8 9]; % hamstrings (BFLH, SM, ST)
        Muscles_antagonist = 6;
    case 'ankle'
        Muscles_spas = [1 2];   % gastrocnemii (GM,GL)
        Muscles_antagonist = 7;
end
% Path to muscle analysis results
Misc.MuscleAnalysisPath=MuscleAnalysis_path;
% Pre-allocation
DatStore = struct('MS',[]);
emg1processedstr = struct('MS',[]);
emg2processedstr = struct('MS',[]);
EMG_ordered = struct('MS',[]);
for ms = 1:Nsegment    
    % Get # dofs, lMT, MA
    [~,Misc.trialName(ms).MS,~]=fileparts(IK_path(ms).MS);
    if ~isfield(Misc,'MuscleNames_Input')||isempty(Misc.MuscleNames_Input)    
        Misc=getMusclesDofs_IPSA(Misc);
    end
    DatStore(ms).MS = getMuscleInfo_IPSA(IK_path(ms).MS,ID_path(ms).MS,...
        Misc,ms);
    % Get MT-parameters
    DatStore(ms).MS.params  = Misc.params_scaled;
    % Get EMG
    temp1 = load(EMG_path_agonist(ms).MS);
    temp2 = load(EMG_path_antagonist(ms).MS);
    emg1processedstr(ms).MS = ...
        temp1.emg1processednorm(Misc.range(ms).MS.range,:);
    emg2processedstr(ms).MS = ...
        temp2.emg2processednorm(Misc.range(ms).MS.range,:);  
    EMG_ordered(ms).MS = zeros(...
        size(emg1processedstr(ms).MS,1),DatStore(ms).MS.nMuscles);   
    EMG_ordered(ms).MS(:,Muscles_spas) = emg1processedstr(ms).MS;
    EMG_ordered(ms).MS(:,Muscles_antagonist) = emg2processedstr(ms).MS;
    % Remove beginning and end of motion to account for bad EMG 
    DatStore(ms).MS.time = DatStore(ms).MS.time(10:end-10,:);
    DatStore(ms).MS.dM = DatStore(ms).MS.dM(10:end-10,:,:);
    DatStore(ms).MS.LMT = DatStore(ms).MS.LMT(10:end-10,:);
    DatStore(ms).MS.T_exp = DatStore(ms).MS.T_exp(10:end-10,:);
    DatStore(ms).MS.EMG = emg1processedstr(ms).MS; % no need antagonists
    DatStore(ms).MS.EMG = DatStore(ms).MS.EMG(10:end-10,:);
    EMG_ordered(ms).MS = EMG_ordered(ms).MS(10:end-10,:);    
end

%% Load sensory information from EMG-driven forward simulations
load([Misc.savefolderForwardSimulation,'\outputGPOPSFTtilde_',...
    Misc.joint]);
load([Misc.savefolderForwardSimulation,'\M_musclesGPOPSFTtilde_',...
    Misc.joint]);
load([Misc.savefolderForwardSimulation,'\M_kneeGPOPSFTtilde_',...
    Misc.joint]);
load([Misc.savefolderForwardSimulation,'\FTGPOPSFTtilde_',...
    Misc.joint]);
load([Misc.savefolderForwardSimulation,'\FtildeGPOPSFTtilde_',...
    Misc.joint]);
load([Misc.savefolderForwardSimulation,'\dFtildeGPOPSFTtilde_',...
    Misc.joint]);
% Pre-allocation
M_muscles = struct('MS',[]);
M_knee = struct('MS',[]);
FT = struct('MS',[]);
aGPOPS = struct('MS',[]);
a = struct('MS',[]);
eGPOPS = struct('MS',[]);
e = struct('MS',[]);
Ftilde = struct('MS',[]);
dFtilde = struct('MS',[]);
for ms = 1:Nsegment
    trial = ['segment_' int2str(Misc.segment_sel(ms))];
    time_GPOPS = outputGPOPS.(trial).result.solution.phase.time;
    M_muscles(ms).MS = ...
        interp1(time_GPOPS,M_musclesGPOPS.(trial),DatStore(ms).MS.time);
    M_knee(ms).MS = ...
        interp1(time_GPOPS,M_kneeGPOPS.(trial),DatStore(ms).MS.time);
    FT(ms).MS = interp1(time_GPOPS,FTGPOPS.(trial),DatStore(ms).MS.time);  
    aGPOPS(ms).MS = outputGPOPS.(trial).result.solution.phase.state(:,...
        1:DatStore(1).MS.nMuscles);
    a(ms).MS = interp1(time_GPOPS,aGPOPS(ms).MS,DatStore(ms).MS.time);
    eGPOPS(ms).MS = outputGPOPS.(trial).result.solution.phase.control(:,...
        1:DatStore(1).MS.nMuscles);
    e(ms).MS = interp1(time_GPOPS,eGPOPS(ms).MS,DatStore(ms).MS.time);
    Ftilde(ms).MS = ...
        interp1(time_GPOPS,FtildeGPOPS.(trial),DatStore(ms).MS.time);
    dFtilde(ms).MS = ...
        interp1(time_GPOPS,dFtildeGPOPS.(trial),DatStore(ms).MS.time);
end

%% Select spastic thresholds 20ms (200Hz, 4 samples) before EMG onset
% We also adjust manually here with the '-9'.
% Pre-allocation
onoff1_range = struct('MS',[]);
onset1_range = struct('MS',[]);
onset1_range_adj = struct('MS',[]);
threshold_Ff = struct('MS',[]);
for ms = 1:Nsegment
    % Select thresholds based on EMG onset
    onoff1_range(ms).MS = Misc.onoff1man(ms).MS(Misc.range(ms).MS.range);
    onset1_range(ms).MS = find(diff(onoff1_range(ms).MS),1,'first')+1;
    switch Misc.joint
        case 'knee'
            onset1_range_adj(ms).MS = onset1_range(ms).MS - 4 - 9;
        case 'ankle'
            onset1_range_adj(ms).MS = onset1_range(ms).MS - 4 - 9;
    end
    threshold_Ff(ms).MS = Ftilde(ms).MS(onset1_range_adj(ms).MS,:);
end

%% Optimal control problem formulation
% Input arguments
auxdata.Nsegment        = Nsegment; % number of segments
auxdata.NMuscles        = DatStore(1).MS.nMuscles; % number of muscles
auxdata.NMuscles_spas   = length(Muscles_spas); % number of spastic muscles
auxdata.Muscles_spas    = Muscles_spas; % indices spastic muscles
auxdata.params          = DatStore(1).MS.params; % Muscle-tendon parameters
auxdata.Ndof            = DatStore(1).MS.nDOF; % number of dofs
auxdata.DOFNames        = DatStore(1).MS.DOFNames; % names of dofs
auxdata.tauAct          = 0.015; % activation time constant
auxdata.tauDeact        = 0.06; % deactivation time constant
auxdata.b               = 0.1; % parameter determining smoothness tanh
auxdata.joint           = Misc.joint; % joint being considered
auxdata.bspas           = 100;
% indices antagonist spastic muscles
auxdata.Muscles_antagonist = Muscles_antagonist;
% Thresholds of spastic models
% Muscle-tendon force threshold is 20ms before EMG onset
% dF/dt threshold is 0
for ms = 1:Nsegment
    for m = 1:auxdata.NMuscles_spas
        mm = auxdata.Muscles_spas(m);
        auxdata.threshold_dFf(ms,m) = 0;
        auxdata.threshold_Ff(ms,m) = threshold_Ff(ms).MS(mm);
    end
end
% Time delay spasticity models
auxdata.Ff_td = 0.03;
auxdata.dFf_td = auxdata.Ff_td;

%% Bounds
% Controls
FTtilde_min = 0;    % bounds on normalized tendon force 
FTtilde_max = 5;    % bounds on normalized tendon force 
dFTtilde_min = -50; % bounds on normalized tendon force derivative
dFTtilde_max = 50;  % bounds on normalized tendon force derivative 
FTtildeMin = FTtilde_min*ones(1,auxdata.NMuscles_spas);     
FTtildeMax = FTtilde_max*ones(1,auxdata.NMuscles_spas);
dFTtildeMin = dFTtilde_min*ones(1,auxdata.NMuscles_spas);   
dFTtildeMax = dFTtilde_max*ones(1,auxdata.NMuscles_spas);
% States 
eFf_min = 0; % bounds on muscle excitations (muscle-tendon force feedback)   
eFf_max = 1; % bounds on muscle excitations (muscle-tendon force feedback)
edFf_min = 0; % bounds on muscle excitations (dF/dt feedback)   
edFf_max = 1; % bounds on muscle excitations (dF/dt feedback)
eFfmin = eFf_min*ones(1,auxdata.NMuscles_spas);     
eFfmax = eFf_max*ones(1,auxdata.NMuscles_spas);
edFfmin = edFf_min*ones(1,auxdata.NMuscles_spas);   
edFfmax = edFf_max*ones(1,auxdata.NMuscles_spas);
% Parameters
g_dFf_min = 0; % bounds on gain dF/dt feedback
g_dFf_max = 15; % bounds on gain dF/dt feedback
g_Ff_min = 0; % bounds on gain muscle-tendon force feedback
g_Ff_max = 1; % bounds on gain muscle-tendon force feedback
gFfmin = g_Ff_min*ones(1,auxdata.NMuscles_spas);    
gFfmax = g_Ff_max*ones(1,auxdata.NMuscles_spas);
gdFfmin = g_dFf_min*ones(1,auxdata.NMuscles_spas);  
gdFfmax = g_dFf_max*ones(1,auxdata.NMuscles_spas); 
bounds.parameter.lower = [gFfmin,gdFfmin]; 
bounds.parameter.upper = [gFfmax,gdFfmax];
% Path constraints
FTtildediff  = zeros(1,auxdata.NMuscles_spas);
dFTtildediff = zeros(1,auxdata.NMuscles_spas);
% Additional constraint: simulated muscle excitation cannot exceed 
% EMG + 0.5*EMG at beginning / end of optimization interval
max_EMG_ordered_all = struct('MS',[]);
max_EMG_ordered_last = struct('MS',[]);
max_EMG_ordered_beg = struct('MS',[]);
for ms = 1:Nsegment    
    max_EMG_ordered_all(ms).MS = max(EMG_ordered(ms).MS(:,Muscles_spas));
    max_EMG_ordered_last(ms).MS = max(EMG_ordered(ms).MS(...
        Misc.range_optimization(ms).MS(end-1:end),Muscles_spas));
    max_EMG_ordered_beg(ms).MS = max(EMG_ordered(ms).MS(...
        Misc.range_optimization(ms).MS(1:10),Muscles_spas));
    % Baseline activation is min EMG from 100 first ms
    auxdata.min_EMG_ordered_beg_mat(ms,:) = ...
        min(EMG_ordered(ms).MS(1:100,Muscles_spas));
end
% Multiple phase optimal control problem
for ms = 1:Nsegment
    % Time
    DatStore(ms).MS.timeopt = ...
        DatStore(ms).MS.time(Misc.range_optimization(ms).MS);
    t0 = DatStore(ms).MS.timeopt(1); 
    tf = DatStore(ms).MS.timeopt(end);
    bounds.phase(ms).initialtime.lower = t0; 
    bounds.phase(ms).initialtime.upper = t0;
    bounds.phase(ms).finaltime.lower = tf; 
    bounds.phase(ms).finaltime.upper = tf;
    % Controls
    bounds.phase(ms).control.lower = [FTtildeMin dFTtildeMin]; 
    bounds.phase(ms).control.upper = [FTtildeMax dFTtildeMax];
    % States
    bounds.phase(ms).initialstate.lower = [eFfmin, edFfmin]; 
    bounds.phase(ms).initialstate.upper = [eFfmax, edFfmax];
    bounds.phase(ms).state.lower = [eFfmin, edFfmin]; 
    bounds.phase(ms).state.upper = [eFfmax, edFfmax];
    bounds.phase(ms).finalstate.lower = [eFfmin, edFfmin]; 
    bounds.phase(ms).finalstate.upper = [eFfmax, edFfmax];
    % Integral bounds
    bounds.phase(ms).integral.lower = 0; 
    bounds.phase(ms).integral.upper = 10000*(tf-t0);
    % Path constraints
    bounds.phase(ms).path.lower = ...
        [-FTtildediff,-dFTtildediff,-max_EMG_ordered_all(ms).MS]; 
    bounds.phase(ms).path.upper = ...
        [FTtildediff,dFTtildediff,max_EMG_ordered_all(ms).MS];
    % Event groups
    temp_eventgroupbeg(:,(ms-1)*auxdata.NMuscles_spas+1:...
        ms*auxdata.NMuscles_spas) = ...
        max_EMG_ordered_beg(ms).MS + 0.5*max_EMG_ordered_beg(ms).MS;
    temp_eventgroupend(:,(ms-1)*auxdata.NMuscles_spas+1:...
        ms*auxdata.NMuscles_spas) = ...
        max_EMG_ordered_last(ms).MS + 0.5*max_EMG_ordered_last(ms).MS;    
end
bounds.eventgroup.lower = [zeros(1,Nsegment*auxdata.NMuscles_spas),...
    zeros(1,Nsegment*auxdata.NMuscles_spas)];
bounds.eventgroup.upper = [temp_eventgroupbeg,temp_eventgroupend];

%% Initial guess
% Spastic gains
IG_Ff = 0.01*ones(1,auxdata.NMuscles_spas);
IG_dFf = 0.05*ones(1,auxdata.NMuscles_spas);
% Time delays
IG_td = auxdata.Ff_td;
% Pre-allocation
uFf_IG = struct('MS',[]);
udFf_IG = struct('MS',[]);
for ms = 1:Nsegment
    for m = 1:auxdata.NMuscles_spas
        mm = auxdata.Muscles_spas(m);
        uFf_IG(ms).MS(:,m) = forceFDynamicsExplicit(...
            DatStore(ms).MS.time,Ftilde(ms).MS(:,mm),IG_td,IG_Ff(m),...
            100,auxdata.threshold_Ff(ms,m));
        udFf_IG(ms).MS(:,m) = dFdtFDynamicsExplicit(...
            DatStore(ms).MS.time,dFtilde(ms).MS(:,mm),IG_td,IG_dFf(m),...
            100,auxdata.threshold_dFf(ms,m));        
    end
end
% Multiple phase optimal control problem
for ms = 1:Nsegment
    guess.phase(ms).time = DatStore(ms).MS.timeopt;
    guess.phase(ms).control = ...
        [Ftilde(ms).MS(Misc.range_optimization(ms).MS,...
        auxdata.Muscles_spas),...
        dFtilde(ms).MS(Misc.range_optimization(ms).MS,...
        auxdata.Muscles_spas)];
    guess.phase(ms).state = ...
        [uFf_IG(ms).MS(Misc.range_optimization(ms).MS,:),...
        udFf_IG(ms).MS(Misc.range_optimization(ms).MS,:)];
    guess.phase(ms).integral = 0;
end
guess.parameter = [IG_Ff,IG_dFf];

%% Spline structures
for ms = 1:Nsegment
    for mm = 1:auxdata.NMuscles_spas
        m = Muscles_spas(mm);
        auxdata.FtildeSpline(ms).MS(mm) = spline(...
            DatStore(ms).MS.timeopt,...
            Ftilde(ms).MS(Misc.range_optimization(ms).MS,m));
        auxdata.dFtildeSpline(ms).MS(mm) = spline(...
            DatStore(ms).MS.timeopt,...
            dFtilde(ms).MS(Misc.range_optimization(ms).MS,m));
        auxdata.EMGSpline(ms).MS(mm) = spline(...
            DatStore(ms).MS.timeopt,...
            DatStore(ms).MS.EMG(Misc.range_optimization(ms).MS,mm));
    end
end

%% GPOPS setup        
setup.name = 'forceModel';
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
NMeshIntervals = struct('MS',[]);
for ms = 1:Nsegment
    NMeshIntervals(ms).MS = round((DatStore(ms).MS.timeopt(end) - ...
        DatStore(ms).MS.timeopt(1))*Misc.Mesh_Frequency);
    setup.mesh.phase(ms).colpoints = 3*ones(1,NMeshIntervals(ms).MS);
    setup.mesh.phase(ms).fraction = (1/(NMeshIntervals(ms).MS))*ones(1,...
        NMeshIntervals(ms).MS);
end
setup.functions.continuous = @forceModel_cont;
setup.functions.endpoint = @forceModel_end;
    
% ADiGator setup
global splinestruct
input.auxdata = auxdata;

for ms = 1:Nsegment
    tdummy = guess.phase(ms).time;
    splinestruct(ms).MS = forceModel_SplineInputData(tdummy,input,ms);
    splinenames = fieldnames(splinestruct(ms).MS);
    for Scount = 1:length(splinenames)
      secdim = size(splinestruct(ms).MS.(splinenames{Scount}),2);
      splinestructad.(splinenames{Scount}) = ...
          adigatorCreateAuxInput([Inf,secdim]);
      splinestruct(ms).MS.(splinenames{Scount}) = zeros(0,secdim);
    end
    setup.auxdata.splinestruct(ms).MS = splinestructad;
end
adigatorGenFiles4gpops2(setup)
setup.functions.continuous      = @Wrap4forceModel_cont;
setup.adigatorgrd.continuous    = @forceModel_contGrdWrap;
setup.adigatorgrd.endpoint      = @forceModel_endADiGatorGrd;
setup.adigatorhes.continuous    = @forceModel_contHesWrap;
setup.adigatorhes.endpoint      = @forceModel_endADiGatorHes;

output = gpops2(setup);

delete forceModel_endADiGatorHes.mat
delete forceModel_endADiGatorHes.m
delete forceModel_endADiGatorGrd.mat
delete forceModel_endADiGatorGrd.m
delete forceModel_contADiGatorHes.mat
delete forceModel_contADiGatorHes.m
delete forceModel_contADiGatorGrd.mat
delete forceModel_contADiGatorGrd.m

end