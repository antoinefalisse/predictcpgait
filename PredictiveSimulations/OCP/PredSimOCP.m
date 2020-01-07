%%  Three-dimensional muscle-driven predictive simulations of walking in
% children with cerebral palsy
%
% Author: Antoine Falisse
% Date: 1/3/2020
%
% DEPENDENCY: please install CasADi (https://web.casadi.org/)

clear all;
clc
close all;

%% User inputs
% This script can be run to solve the optimal control problems but also to
% analyze and process the results. The user does not need to re-run the
% optimal control problems to analyze the results. Therefore, the user can
% select some settings beforehand through the variable num_set. For
% example, when num_set(1) is set to 0, the script will not run the
% optimal control problem. Here is a brief description of num_set:
% num_set(1): set to 1 to solve problem
% num_set(2): set to 1 to analyze results
% num_set(3): set to 1 to load results
% num_set(4): set to 1 to save results
% num_set(5): set to 1 to visualize guess vs bounds 
% num_set(6): set to 1 to write motion file starting at right heel strike
% num_set(7): set to 1 to write motion file starting at left heel strike

num_set = [1,0,0,0,0,0,0]; % This configuration solves the problem
% num_set = [0,1,1,1,0,1,0]; % This configuration analyzes the results

% The variable 'settings', loaded through PredSimOCP_settings in the following
% section, sets some parameters of the problem (e.g., weights in the cost
% function). Through the variable idx_settings, the user can select which row of
% parameters is used.
idx_settings = 1; % Index row in matrix settings

%% Settings
import casadi.*
subject = 'subject1';

solveProblem    = num_set(1); % set to 1 to solve problem
analyseResults  = num_set(2); % set to 1 to analyze results
loadResults     = num_set(3); % set to 1 to load results
saveResults     = num_set(4); % set to 1 to save results
checkBoundsIG   = num_set(5); % set to 1 to visualize guess vs bounds 
writeMotion_r   = num_set(6); % set to 1 to write motion file starting at right heel strike
writeMotion_l   = num_set(7); % set to 1 to write motion file starting at left heel strike

% Paths
pathmain = pwd;
[pathPredictiveSimulations,~,~] = fileparts(pathmain);
p = mfilename('fullpath');
[~,namescript,~] = fileparts(p);
[pathRepo,~,~] = fileparts(pathPredictiveSimulations);
pathSettings = [pathPredictiveSimulations,'/Settings'];
addpath(genpath(pathSettings));
PredSimOCP_settings

% Loop over set of parameters used in the problem
for www = 1:length(idx_settings)
ww = idx_settings(www);
% Variable parameters
W.TrackTD	= settings(ww,1);	% weight TD kinematics tracking
W.Act       = settings(ww,2);   % weight muscle activations (fatigue)
N           = settings(ww,3);   % number of mesh intervals
tol_ipopt   = settings(ww,4);   % NLP error tolerance: 1*10^(-tol_ipopt)
NSyn        = settings(ww,5);   % number of synergies per leg
paramsi     = settings(ww,6);   % MT-parameters (generic or personalized)
IGi         = settings(ww,7);   % initial guess identifier
W.MER       = settings(ww,8);   % weight metabolic energy rate
exp_E       = settings(ww,9);   % exponent metabolic energy rate
exp_A       = settings(ww,10);  % exponent muscle activations
spasi       = settings(ww,11);  % whether to take spasticity into account
W.Acc       = settings(ww,12);  % weight joint accelerations
gaitSpeed	= settings(ww,13);  % imposed gait speed (m s-1)
W.TrunkExc 	= settings(ww,14);  % weight trunk excitations
W.PassT     = settings(ww,15);  % weight passive joint torques
IGc         = settings(ww,16);  % initial guess case
% Fixed parameter
W.u = 0.001; % weight slack variables

% The problem formulation can be "easily" adjusted for use with multiple
% phases. For instance, the user may want to track the kinematics of
% multiple gait trials while optimizing a set of parameters. This is not
% the case here, so we set NPhases to 1. Yet we keep the multiple phase
% structure in case users are interested in that feature.
NPhases = 1;

% The filename is associated with the index of the row of parameters
% selected in 'settings'.
savename = ['_c',num2str(ww)];
savenamePhase = struct('mp',[]); 
phaseID = {'a','b','c','d'}; % in case multiple phases are used.
for mpi = 1:NPhases
    savenamePhase(mpi).mp = ['_c',num2str(ww),phaseID{mpi}];
end

%% External function
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.
pathExternalFunctions = [pathPredictiveSimulations,'/ExternalFunctions'];
% Loadi external function   
cd(pathExternalFunctions);
F = external('F','PredSim_CP.dll');
cd(pathmain);
% This is an example of how to call an external function with numerical values.
% inF = rand(63,1);
% outF = full(F(inF));

% Indices elements in external function F
% First: joint torques. 
jointi.pelvis.list  = 1; 
jointi.pelvis.rot   = 2; 
jointi.pelvis.tilt  = 3; 
jointi.pelvis.tx    = 4;
jointi.pelvis.ty    = 5;
jointi.pelvis.tz    = 6;
jointi.hip_flex.l   = 7;
jointi.hip_add.l    = 8;
jointi.hip_rot.l    = 9;
jointi.hip_flex.r   = 10;
jointi.hip_add.r    = 11;
jointi.hip_rot.r    = 12;
jointi.knee.l       = 13;
jointi.knee.r       = 14;
jointi.ankle.l      = 15;
jointi.ankle.r      = 16;
jointi.subt.l       = 17;
jointi.subt.r       = 18;
jointi.trunk.ext    = 19;
jointi.trunk.ben    = 20;
jointi.trunk.rot    = 21;
% Helper variables
% Vectors of indices
jointi.resi = jointi.pelvis.list:jointi.trunk.rot; % all 
jointi.res_bgpi = jointi.hip_flex.l:jointi.trunk.rot; % all but ground-pelvis
jointi.res_bptyi = [jointi.pelvis.list:jointi.pelvis.tx,...
    jointi.pelvis.tz:jointi.trunk.rot]; % all but pelvis_ty
jointi.res_trunki = jointi.trunk.ext:jointi.trunk.rot; % trunk
jointi.res_roti = [jointi.pelvis.list:jointi.pelvis.tilt,...
    jointi.hip_flex.l:jointi.trunk.rot]; % rotational DoFs
jointi.res_tri = jointi.pelvis.tx:jointi.pelvis.tz; % translational DoFs
jointi.res_gpi = jointi.pelvis.list:jointi.pelvis.tz; % ground-pelvis
% Numbers of DoFs
nq.res = length(jointi.resi); % all 
nq.res_gp = length(jointi.res_gpi); % ground-pelvis
nq.res_bgp = nq.res-nq.res_gp;% all but ground-pelvis
nq.leg = 6; % #joints needed for polynomials
nq.res_trunk = length(jointi.res_trunki); % trunk
% In the external function, joint positions and velocities are intertwined.
% Joint positions are at locations 1:2:end.
Qsi = 1:2:2*nq.res; % indices Qs only 
% Second: ground reaction forces (GRFs)
GRFi.r      = 22:24;
GRFi.l      = 25:27;
% Third: ground reaction moments (GRMs)
GRMi.r      = 28:30;
GRMi.l      = 31:33;    
GRFi.all 	= [GRFi.r,GRFi.l];
nGRF     	= length(GRFi.all);
GRMi.all   	= [GRMi.r,GRMi.l];
nGRM       	= length(GRMi.all);
% Fourth: joint center locations in the transversal plane
% Calcaneus
calcOr.r    = 34:35;
calcOr.l    = 36:37;  
calcOr.all  = [calcOr.r,calcOr.l];
NcalcOr     = length(calcOr.all);
% Tibias
tibiaOr.r   = 38:39;
tibiaOr.l   = 40:41;        
tibiaOr.all = [tibiaOr.r,tibiaOr.l];
NtibiaOr    = length(tibiaOr.all);

%% Subject mass
body_mass = 3.69727 + 4.64393 + 1.43323 + 0.01986 + 0.39058 + 0.04303 + ...
    4.64393 + 1.43323 + 0.01986 + 0.39058 + 0.04303 + 16.25541;

%% Collocation scheme
% We use a pseudospectral direct collocation method, i.e. we use Lagrange
% polynomials to approximate the state derivatives at the collocation
% points in each mesh interval. We use d=3 collocation points per mesh
% interval and Radau collocation points. 
pathCollocationScheme = [pathRepo,'/CollocationScheme'];
addpath(genpath(pathCollocationScheme));
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[tau_root,C,D,B] = CollocationScheme(d,method);

%% Muscle-tendon parameters 
% Muscle names from right leg
muscleNames = {'glut_max1_r','glut_max2_r','glut_max3_r','glut_med1_r',...
    'glut_med2_r','glut_med3_r','glut_min1_r','glut_min2_r',...
    'glut_min3_r','add_long_r','add_brev_r','add_mag1_r','add_mag2_r',...
    'add_mag3_r','pectineus_r','iliacus_r','psoas_r','quad_fem_r',...
    'gemellus_r','piri_r','TFL_r','gracilis_r','semimem_r','semiten_r',...
    'bi_fem_lh_r','bi_fem_sh_r','sartorius_r','rectus_fem_r',...
    'vas_med_r','vas_int_r','vas_lat_r','gas_med_r',...
    'gas_lat_r','soleus_r','tib_post_r','tib_ant_r','ext_dig_r',...
    'ext_hal_r','flex_dig_r','flex_hal_r','per_brev_r','per_long_r',...
    'per_tert_r'};
NMuscles = length(muscleNames)*2; % total number of muscles
% Muscle indices for later use
pathmusclemodel = [pathRepo,'/MuscleModel'];
addpath(genpath(pathmusclemodel));    
musi = muscleIndices(muscleNames);
% Load muscle-tendon parameters
% Muscle-tendon parameters. Row 1: maximal isometric forces; Row 2: optimal
% fiber lengths; Row 3: tendon slack lengths; Row 4: optimal pennation 
% angles; Row 5: maximal contraction velocities
pathOpenSimModel = [pathRepo,'/OpenSimModel/',subject,'/'];
pathParameterEstimation = [pathOpenSimModel,'ParameterEstimation/'];    
if paramsi == 1 % generic parameters
    load([pathParameterEstimation,'MTParameters_generic.mat']);     
elseif paramsi == 2 % personalized parameters
    load([pathParameterEstimation,'MTParameters_personalized.mat']);  
end
% Indices of the muscles actuating the different joints for later use
tl = load([pathOpenSimModel,'Polynomials/muscle_spanning_joint_INFO_',...
    subject,'_r.mat']);
[~,mai] = getMomentArmIndices(muscleNames,tl.muscle_spanning_joint_INFO_r);
% The original code was set up such that the force-length-velocity (FLV)
% relationships of some muscles could be ignored. This is not used in this code
% but we maintain this feature in case it is of interest for future users.
% Indices of muscles whose FLV relationships are not used
muscleNames_noFLV = {};
NMuscles_noFLV_r  = length(muscleNames_noFLV); % one side
musi_noFLV_r = zeros(1,NMuscles_noFLV_r);
for i = 1:NMuscles_noFLV_r
    musi_noFLV_r(i) = find(strcmp(muscleNames,muscleNames_noFLV{i}));
end
musi_noFLV = [musi_noFLV_r,musi_noFLV_r+NMuscles/2]; % both sides
NMuscles_noFLV = 2*NMuscles_noFLV_r; % both sides
NMuscles_FLV = NMuscles-NMuscles_noFLV;
musi_FLV_r = musi;
musi_FLV_r(musi_noFLV_r) = [];
musi_FLV = [musi,musi+NMuscles/2];
musi_FLV(musi_noFLV) = [];
MTParameters_FLV = MTParameters(:,musi_FLV);
MTParameters_noFLV = MTParameters(:,musi_noFLV);
% Parameters for activation dynamics
tact = 0.015;
tdeact = 0.06;

%% Metabolic energy model
% We extract the specific tensions and slow twitch rations.
pathMetabolicEnergy = [pathPredictiveSimulations,'/MetabolicEnergy'];
addpath(genpath(pathMetabolicEnergy));
tension = getSpecificTensions(muscleNames); 
tensions = [tension;tension];
pctst = getSlowTwitchRatios(muscleNames); 
pctsts = [pctst;pctst];

%% Spasticity
muscleNames_Spas = ...
    {'bi_fem_lh_r','semimem_r','semiten_r','gas_med_r','gas_lat_r'};
NMuscles_Spas_r  = length(muscleNames_Spas); % one side
musi_Spas_r = zeros(1,NMuscles_Spas_r);
musi_Spas_In_FLV_r = zeros(1,NMuscles_Spas_r);
for i = 1:NMuscles_Spas_r
    musi_Spas_r(i) = find(strcmp(muscleNames,muscleNames_Spas{i}));
    musi_Spas_In_FLV_r(i) = ...
        find(strcmp(muscleNames(musi_FLV_r),muscleNames_Spas{i}));    
end
musi_Spas = [musi_Spas_r,musi_Spas_r+NMuscles/2]; % both sides
musi_Spas_In_FLV = [musi_Spas_In_FLV_r,musi_Spas_In_FLV_r+NMuscles_FLV/2];
NMuscles_Spas = 2*NMuscles_Spas_r; % both sides
musi_noSpas = [musi,musi+NMuscles/2];
musi_noSpas(musi_Spas) = [];
% Load feedback / spasticity model parameters in case spasticity is taken
% into account
if spasi == 1 
    formulation = 'forceModel';
    pathSpasticity = ...
        [pathOpenSimModel,'Spasticity/SpasticGainEstimation/',formulation];     
    switch subject
        case 'subject1'
            segment_sel_hamstings = [6 7 8];
            segment_sel_gastrocs = [19 20 21];
    end          
    % Hamstrings
    joint = 'knee'; 
    number_save_fM = [formulation '_' int2str(segment_sel_hamstings)];
    load([pathSpasticity,'/output_',joint,'_',number_save_fM]);         
    NMuscles_spas = output.result.setup.auxdata.NMuscles_spas;
    gFf.ham = output.result.solution.parameter(:,1:NMuscles_spas);
    gdFf.ham = ...
        output.result.solution.parameter(:,NMuscles_spas+1:2*NMuscles_spas);
    tauFf.ham = output.result.setup.auxdata.Ff_td;
    taudFf.ham = output.result.setup.auxdata.dFf_td;         
    threshold_Ff.ham = output.result.setup.auxdata.threshold_Ff;
    threshold_dFf.ham = output.result.setup.auxdata.threshold_dFf;   
    % The thresholds during gait are the lowest thresholds during passive
    % motions (IPSA)
    threshold_dFf_gait.ham = min(threshold_dFf.ham);
    threshold_Ff_gait.ham = min(threshold_Ff.ham);
    % Gastrocs
    joint = 'ankle';
    number_save_fM = [formulation '_' int2str(segment_sel_gastrocs)];
    load([pathSpasticity,'/output_',joint,'_',number_save_fM]);  
    % Extract data        
    NMuscles_spas = output.result.setup.auxdata.NMuscles_spas;
    gFf.gas = output.result.solution.parameter(:,1:NMuscles_spas);
    gdFf.gas = ...
        output.result.solution.parameter(:,NMuscles_spas+1:2*NMuscles_spas);
    tauFf.gas = output.result.setup.auxdata.Ff_td;
    taudFf.gas = output.result.setup.auxdata.dFf_td;         
    threshold_Ff.gas = output.result.setup.auxdata.threshold_Ff;
    threshold_dFf.gas = output.result.setup.auxdata.threshold_dFf;   
    % The thresholds during gait are the lowest thresholds during passive
    % motions (IPSA)
    threshold_dFf_gait.gas = min(threshold_dFf.gas);
    threshold_Ff_gait.gas = min(threshold_Ff.gas);
    % Combine both muscle groups
    gFf.side = [gFf.ham,gFf.gas];
    gdFf.side = [gdFf.ham,gdFf.gas];
    tauFf.all = tauFf.ham;
    taudFf.all = taudFf.ham;
    threshold_dFf_gait.side = [threshold_dFf_gait.ham,threshold_dFf_gait.gas];
    threshold_Ff_gait.side = [threshold_Ff_gait.ham,threshold_Ff_gait.gas];
    % For now, we assume bilateral spasticity with same involvment
    gFf.all = [gFf.side,gFf.side];
    gdFf.all = [gdFf.side,gdFf.side];
    threshold_dFf_gait.all = [threshold_dFf_gait.side,threshold_dFf_gait.side];
    threshold_Ff_gait.all = [threshold_Ff_gait.side,threshold_Ff_gait.side];
    % Parameter determining the smoothness of the tanh approximation
    bspas = 100;        
end

%% Synergies
pathSynergies = [pathOpenSimModel,'Synergies/'];
synData = load([pathSynergies,'/Synergies.mat']);
if NSyn == 99 % no synergies imposed
    NMact = NMuscles;
    synW.l = synData.Left.W{1,4};
    synW.r = synData.Right.W{1,4};
else
    NMact = 2*NSyn;
    synW.l = synData.Left.W{1,NSyn};
    synW.r = synData.Right.W{1,NSyn};
end

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathPredictiveSimulations,'/CasADiFunctions'];
addpath(genpath(pathCasADiFunctions));
% We load some variables for the polynomial approximations
load([pathOpenSimModel,'Polynomials/muscle_spanning_joint_INFO_',...
    subject,'_r.mat']);
load([pathOpenSimModel,'Polynomials/muscle_spanning_joint_INFO_',...
    subject,'_l.mat']);
load([pathOpenSimModel,'Polynomials/MuscleInfo_',subject,'_r.mat']);
load([pathOpenSimModel,'Polynomials/MuscleInfo_',subject,'_l.mat']);
% For the polynomials, we want all independent muscles. So we do not need 
% the muscles from both legs, since we assume bilateral symmetry, but want
% all muscles from the trunk (indices 47:49).
musi_pol = musi;
NMuscles_pol = NMuscles/2;
PredSimOCP_CasADiFunctions

%% Passive joint torques
% Parameters for the passive torques of the lower limbs and trunk
pathPassiveTorques = [pathPredictiveSimulations,'/PassiveTorques'];
addpath(genpath(pathPassiveTorques));
PassiveTorquesData

%% Data from a tracking simulation of a TD walking trial.
% This data is used for settings bounds, initial guesses, and as part of the
% optimal control problem in some cases. This data represents what is
% achievable by the model in terms of tracking experimental data with null
% pelvis residuals.
joints = {'lower_torso_RX','lower_torso_RY','lower_torso_RZ',...
    'lower_torso_TX','lower_torso_TY','lower_torso_TZ','hip_flex_l',...
    'hip_add_l','hip_rot_l','hip_flex_r','hip_add_r','hip_rot_r',...
    'knee_flex_l','knee_flex_r','ankle_flex_l','ankle_flex_r',...
    'subt_angle_l','subt_angle_r','lumbar_pitch','lumbar_roll',...
    'lumbar_yaw'};
pathReferenceTracking = [pathOpenSimModel,'ReferenceTracking'];
load([pathReferenceTracking,'/ReferenceTracking.mat']);

Qs = struct('mp',[]); 
GRF = struct('mp',[]); 
ID = struct('mp',[]);
for mpi = 1:NPhases
    % Joint angles
    Qs_temp = ReferenceTracking.Results.Qs_opt(mpi).mp;
    Qs_temp(:,end+1:end+nq.res_trunk) = zeros(size(Qs_temp,1),nq.res_trunk);
    Qs_temp(:,jointi.res_roti) = Qs_temp(:,jointi.res_roti)*pi/180;
    % Ground reaction forces
    GRF_temp = ReferenceTracking.Results.GRFs_opt(mpi).mp;
    % Ground reaction moments
    GRM_temp = ReferenceTracking.Results.GRMs_opt(mpi).mp;
    % Joint torques
    ID_temp = ReferenceTracking.Results.Ts_opt(mpi).mp;
    ID_temp(:,end+1:end+nq.res_trunk) = zeros(size(ID_temp,1),nq.res_trunk);
    % Apply a low-pass filter to the trajectories
    order = 4;
    cutoff_low = 20;
    fs=1/mean(diff(ReferenceTracking.Results.tgrid(mpi).mp));
    [af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
    Qs_tempfilt = filtfilt(af,bf,Qs_temp);  
    Qs_temp_filt = [ReferenceTracking.Results.tgrid(mpi).mp,Qs_tempfilt];
    GRF_tempfilt = filtfilt(af,bf,GRF_temp);  
    GRF_temp_filt = [ReferenceTracking.Results.tgrid(mpi).mp,GRF_tempfilt];
    GRM_tempfilt = filtfilt(af,bf,GRM_temp);  
    GRM_temp_filt = [ReferenceTracking.Results.tgrid(mpi).mp,GRM_tempfilt];
    ID_tempfilt = filtfilt(af,bf,ID_temp);  
    ID_temp_filt = [ReferenceTracking.Results.tgrid(mpi).mp,ID_tempfilt];       
    % Adjust the time vector (shorter because N controls and N+1 states)
    time_opt(1) = ReferenceTracking.Results.tgrid(mpi).mp(1);
    time_opt(2) = ReferenceTracking.Results.tgrid(mpi).mp(end);                
    % Interpolation to number of mesh intervals
    step = (time_opt(2)-time_opt(1))/(N-1);
    interval = time_opt(1):step:time_opt(2);        
    Qs(mpi).mp.allinterpfilt = ...
        interp1(Qs_temp_filt(:,1),Qs_temp_filt,interval);
    Qs(mpi).mp.colheaders = ['Time',ReferenceTracking.colheaders.joints];
    ID(mpi).mp.allinterp = interp1(ID_temp_filt(:,1),ID_temp_filt,interval);        
    GRF(mpi).mp.val.allinterp = ...
        interp1(GRF_temp_filt(:,1),GRF_temp_filt,interval);
    GRF(mpi).mp.val.all = GRF(mpi).mp.val.allinterp;
    GRF(mpi).mp.MorGF.allinterp = ...
        interp1(GRM_temp_filt(:,1),GRM_temp_filt,interval);           
end  

%% Bounds
pathBounds = [pathPredictiveSimulations,'/Bounds'];
addpath(genpath(pathBounds));
pathVariousFunctions = [pathRepo,'/VariousFunctions'];
addpath(genpath(pathVariousFunctions));
bounds = struct('mp',[]); 
scaling = struct('mp',[]); 
Qs_CP = struct('mp',[]);
for mpi = 1:NPhases
    % Bounds accounting for both TD and CD walking patterns
    pathIK = [pathOpenSimModel,'IK/Gait/KS_Gait_average_r.mot'];
    Qs_CP(mpi).mp = getIK(pathIK,joints);           
    step = (Qs_CP(mpi).mp.time(end)-Qs_CP(mpi).mp.time(1))/(N-1);
    interval = Qs_CP(mpi).mp.time(1):step:Qs_CP(mpi).mp.time(end);        
    Qs_CP(mpi).mp.allinterpfilt = interp1(Qs_CP(mpi).mp.allfilt(:,1),...
        Qs_CP(mpi).mp.allfilt,interval);  
    [bounds(mpi).mp,scaling(mpi).mp] = getBounds(Qs(mpi).mp,Qs_CP(mpi).mp,...
        NMuscles,NMuscles_FLV,nq,jointi,GRF(mpi).mp,N,NSyn,NMuscles_Spas);
end

%% Initial guess
pathBounds = [pathPredictiveSimulations,'/Guess'];
addpath(genpath(pathBounds));
% Data-informed initial guess
guess = struct('mp',[]); 
Qs_CP = struct('mp',[]);
Qs_ig = struct('mp',[]);
for mpi = 1:NPhases
    if IGi == 1
        % IG is CP walking pattern
        pathIK = [pathOpenSimModel,'IK/Gait/KS_Gait_average_r.mot'];
        Qs_CP(mpi).mp = getIK(pathIK,joints);           
        step = (Qs_CP(mpi).mp.time(end)-Qs_CP(mpi).mp.time(1))/(N-1);
        interval = Qs_CP(mpi).mp.time(1):step:Qs_CP(mpi).mp.time(end);        
        Qs_CP(mpi).mp.allinterpfilt = interp1(Qs_CP(mpi).mp.allfilt(:,1),...
            Qs_CP(mpi).mp.allfilt,interval);        
        guess(mpi).mp = getGuess(Qs_CP(mpi).mp,nq,N,NMuscles,NMuscles_FLV,...
            jointi,scaling(mpi).mp,NSyn,NMuscles_Spas);   
    elseif IGi == 2
        % IG is TD walking pattern
        guess(mpi).mp = getGuess(Qs(mpi).mp,nq,N,NMuscles,NMuscles_FLV,...
            jointi,scaling(mpi).mp,NSyn,NMuscles_Spas);
    elseif IGi == 3
        % IG is existing simulated movement      
        leg = 'r';
        savename_ig = ['_c',num2str(IGc),'a_',leg];
        pathIK_ig = [pathPredictiveSimulations,'/Results/',namescript,...
            '/IK',savename_ig,'.mot'];
        Qs_ig(mpi).mp = getIK_MRI(pathIK_ig,joints);           
        step = (Qs_ig(mpi).mp.time(end)-Qs_ig(mpi).mp.time(1))/(N-1);        
        interval = Qs_ig(mpi).mp.time(1):step:Qs_ig(mpi).mp.time(end);        
        Qs_ig(mpi).mp.allinterpfilt = interp1(Qs_ig(mpi).mp.allfilt(:,1),...
            Qs_ig(mpi).mp.allfilt,interval);  
        guess(mpi).mp =  getGuess(Qs_ig(mpi).mp,nq,N,NMuscles,NMuscles_FLV,...
            jointi,scaling(mpi).mp,NSyn,NMuscles_Spas);  
    end
end
% This allows visualizing the initial guess and the bounds
if checkBoundsIG
    for mpi = 1:NPhases
        plot_BoundsVSInitialGuess
    end
end

%% Formulate the NLP
if solveProblem
    % Start with an empty NLP
    w   = {}; % design variables
    w0  = []; % initial guess for design variables
    lbw = []; % lower bounds for design variables
    ubw = []; % upper bounds for design variables
    J   = 0;  % initial value of cost function
    g   = {}; % constraints
    lbg = []; % lower bounds for constraints
    ubg = []; % upper bounds for constraints  
    % Final time
    tf              = MX.sym('tf',1);
    w               = [w {tf}];
    lbw             = [lbw; bounds(1).mp.tf.lower];
    ubw             = [ubw; bounds(1).mp.tf.upper];
    w0              = [w0;  guess(1).mp.tf];
    % When synergies are imposed, synergy weights are static parameters
    if NSyn ~= 99 % no synergies
        % Synergy weights: left side
        syn_wl          = MX.sym('syn_wl',NMuscles/2*NSyn);
        w               = [w {syn_wl}];
        lbw             = [lbw; bounds(1).mp.synw.lower.l];
        ubw             = [ubw; bounds(1).mp.synw.upper.l];
        w0              = [w0;  guess(1).mp.synw.l(:)];
        % Synergy weights: right side
        syn_wr          = MX.sym('syn_wr',NMuscles/2*NSyn);
        w               = [w {syn_wr}];
        lbw             = [lbw; bounds(1).mp.synw.lower.r];
        ubw             = [ubw; bounds(1).mp.synw.upper.r];
        w0              = [w0;  guess(1).mp.synw.r(:)];  
    end
    % Loop over phases
    for mpi = 1:NPhases  
        % Define states at first mesh point
        % Muscle activations
        a0              = MX.sym('a0',NMact);
        w               = [w {a0}];
        lbw             = [lbw; bounds(mpi).mp.a.lower'];
        ubw             = [ubw; bounds(mpi).mp.a.upper'];
        w0              = [w0;  guess(mpi).mp.a(1,:)'];
        % Muscle-tendon forces
        FTtilde0        = MX.sym('FTtilde0',NMuscles_FLV);
        w               = [w {FTtilde0}];
        lbw             = [lbw; bounds(mpi).mp.FTtilde.lower'];
        ubw             = [ubw; bounds(mpi).mp.FTtilde.upper'];
        w0              = [w0;  guess(mpi).mp.FTtilde(1,:)'];
        % Qs and Qdots
        X0              = MX.sym('X0',2*nq.res);
        w               = [w {X0}];    
        lbw             = [lbw; bounds(mpi).mp.QsQdots.lower(1,:)'];
        ubw             = [ubw; bounds(mpi).mp.QsQdots.upper(1,:)'];    
        w0              = [w0;  guess(mpi).mp.QsQdots(1,:)']; 
        % When spasticity is taken into account, feedback muscle
        % activations are states.
        if spasi == 1
            % Muscle activations from muscle-tendon force feedback
            a_Ff0           = MX.sym('a_Ff0',NMuscles_Spas);
            w               = [w {a_Ff0}];
            lbw             = [lbw; bounds(mpi).mp.a_Ff.lower'];
            ubw             = [ubw; bounds(mpi).mp.a_Ff.upper'];
            w0              = [w0;  guess(mpi).mp.a_Ff(1,:)'];    
            % Muscle activations from muscle-tendon force rate feedback
            a_dFf0           = MX.sym('a_dFf0',NMuscles_Spas);
            w               = [w {a_dFf0}];
            lbw             = [lbw; bounds(mpi).mp.a_dFf.lower'];
            ubw             = [ubw; bounds(mpi).mp.a_dFf.upper'];
            w0              = [w0;  guess(mpi).mp.a_dFf(1,:)'];   
        end
        % Trunk activations
        a_b0            = MX.sym('a_b0',nq.res_trunk);
        w               = [w {a_b0}];
        lbw             = [lbw; bounds(mpi).mp.a_b.lower'];
        ubw             = [ubw; bounds(mpi).mp.a_b.upper'];
        w0              = [w0;  guess(mpi).mp.a_b(1,:)'];    
        % We pre-allocate some of the states so that we can provide an
        % expression for the distance traveled
        for k=0:N
            Xk{k+1,1} = MX.sym(['X_' num2str(k+1)], 2*nq.res);
        end   
        % "Lift" initial conditions
        ak = a0;
        FTtildek = FTtilde0;
        Xk{1,1} = X0;
        if spasi == 1
            a_Ffk = a_Ff0; 
            a_dFfk = a_dFf0;
        end
        a_bk = a_b0; 
        % Provide expression for the distance traveled
        % initial position pelvis_tx
        pelvis_tx0 = Xk{1,1}(2*jointi.pelvis.tx-1,1).*...
            scaling(mpi).mp.QsQdots(2*jointi.pelvis.tx-1);  
         % final position pelvis_tx: N and not N+1 to match experimental data 
        pelvis_txf = Xk{N,1}(2*jointi.pelvis.tx-1,1).*...
            scaling(mpi).mp.QsQdots(2*jointi.pelvis.tx-1);   
        % distance traveled  
        dist_trav_tot = pelvis_txf-pelvis_tx0;
        % Time step
        h = tf/N;
        % Loop over mesh points
        for k=0:N-1
            % Define controls at mesh point (piecewise-constant in interval) 
            % Time derivative of muscle activations (states)
            vAk                 = MX.sym(['vA_' num2str(k)],NMact);
            w                   = [w {vAk}];
            lbw                 = [lbw; bounds(mpi).mp.vA.lower'];
            ubw                 = [ubw; bounds(mpi).mp.vA.upper'];
            w0                  = [w0; guess(mpi).mp.vA(k+1,:)'];
            % Time derivative of muscle-tendon forces (states)
            dFTtildek           = MX.sym(['dFTtilde_' num2str(k)],NMuscles_FLV);
            w                   = [w {dFTtildek}];
            lbw                 = [lbw; bounds(mpi).mp.dFTtilde.lower'];
            ubw                 = [ubw; bounds(mpi).mp.dFTtilde.upper'];
            w0                  = [w0; guess(mpi).mp.dFTtilde(k+1,:)'];  
            % Time derivative of Qdots (states) 
            Ak                  = MX.sym(['A_' num2str(k)],nq.res);
            w                   = [w {Ak}];
            lbw                 = [lbw; bounds(mpi).mp.Qdotdots.lower'];
            ubw                 = [ubw; bounds(mpi).mp.Qdotdots.upper'];
            w0                  = [w0; guess(mpi).mp.Qdotdots(k+1,:)'];
            % Trunk excitations
            e_bk                = MX.sym(['e_b_' num2str(k)],nq.res_trunk);
            w                   = [w {e_bk}];
            lbw                 = [lbw; bounds(mpi).mp.e_b.lower'];
            ubw                 = [ubw; bounds(mpi).mp.e_b.upper'];
            w0                  = [w0; guess(mpi).mp.e_b(k+1,:)'];                
            % Define states at collocation points    
            % Muscle activations
            akj = {};
            for j=1:d
                akj{j}  = MX.sym(['	a_' num2str(k) '_' num2str(j)],NMact);
                w       = {w{:}, akj{j}};
                lbw     = [lbw; bounds(mpi).mp.a.lower'];
                ubw     = [ubw; bounds(mpi).mp.a.upper'];
                w0      = [w0;  guess(mpi).mp.a(k+1,:)'];
            end   
            % Muscle-tendon forces
            FTtildekj = {};
            for j=1:d
                FTtildekj{j} = ...
                    MX.sym(['FTtilde_' num2str(k) '_' num2str(j)],NMuscles_FLV);
                w            = {w{:}, FTtildekj{j}};
                lbw          = [lbw; bounds(mpi).mp.FTtilde.lower'];
                ubw          = [ubw; bounds(mpi).mp.FTtilde.upper'];
                w0           = [w0;  guess(mpi).mp.FTtilde(k+1,:)'];
            end
            % Qs and Qdots        
            Xkj = {};
            for j=1:d
                Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)],2*nq.res);
                w      = {w{:}, Xkj{j}};
                lbw    = [lbw; bounds(mpi).mp.QsQdots.lower(k+1,:)'];
                ubw    = [ubw; bounds(mpi).mp.QsQdots.upper(k+1,:)'];
                w0     = [w0;  guess(mpi).mp.QsQdots(k+1,:)'];
            end   
            if spasi == 1
                % Muscle activations from muscle-tendon force feedback
                a_Ffkj = {};
                for j=1:d
                    a_Ffkj{j} = MX.sym(['a_Ff_' num2str(k) '_' num2str(j)], ...
                        NMuscles_Spas);
                    w       = {w{:}, a_Ffkj{j}};
                    lbw     = [lbw; bounds(mpi).mp.a_Ff.lower'];
                    ubw     = [ubw; bounds(mpi).mp.a_Ff.upper'];
                    w0      = [w0;  guess(mpi).mp.a_Ff(k+1,:)'];
                end
                % Muscle activations from muscle-tendon force rate feedback
                a_dFfkj = {};
                for j=1:d
                    a_dFfkj{j}= MX.sym(['a_dFf_' num2str(k) '_' num2str(j)], ...
                        NMuscles_Spas);
                    w       = {w{:}, a_dFfkj{j}};
                    lbw     = [lbw; bounds(mpi).mp.a_dFf.lower'];
                    ubw     = [ubw; bounds(mpi).mp.a_dFf.upper'];
                    w0      = [w0;  guess(mpi).mp.a_dFf(k+1,:)'];
                end
            end
            % Trunk activations
            a_bkj = {};
            for j=1:d
                a_bkj{j}= MX.sym(['	a_b_' num2str(k) '_' num2str(j)], ...
                    nq.res_trunk);
                w       = {w{:}, a_bkj{j}};
                lbw     = [lbw; bounds(mpi).mp.a_b.lower'];
                ubw     = [ubw; bounds(mpi).mp.a_b.upper'];
                w0      = [w0;  guess(mpi).mp.a_b(k+1,:)'];
            end 
            % Unscale variables for later use
            Xk_nsc          = Xk{k+1,1}.*scaling(mpi).mp.QsQdots';
            FTtildek_nsc    = FTtildek.*(scaling(mpi).mp.FTtilde');
            Ak_nsc          = Ak.*scaling(mpi).mp.Qdotdots';
            for j=1:d
                Xkj_nsc{j} = Xkj{j}.*scaling(mpi).mp.QsQdots';
                FTtildekj_nsc{j} = FTtildekj{j}.*scaling(mpi).mp.FTtilde';
            end                                  
            % Call external function
            [Tk] = F([Xk_nsc;Ak_nsc]);
            Tauk_abs = Tk(jointi.res_gpi,1);
            % Get muscle-tendon lengths, velocities, and moment arms
            % Left leg
            qin_l = [Xk_nsc(jointi.hip_flex.l*2-1,1),...
                Xk_nsc(jointi.hip_add.l*2-1,1), ...
                Xk_nsc(jointi.hip_rot.l*2-1,1), ...
                Xk_nsc(jointi.knee.l*2-1,1), ...
                Xk_nsc(jointi.ankle.l*2-1,1),...
                Xk_nsc(jointi.subt.l*2-1,1)];  
            qdotin_l = [Xk_nsc(jointi.hip_flex.l*2,1),...
                Xk_nsc(jointi.hip_add.l*2,1),...
                Xk_nsc(jointi.hip_rot.l*2,1),...
                Xk_nsc(jointi.knee.l*2,1),...
                Xk_nsc(jointi.ankle.l*2,1),...
                Xk_nsc(jointi.subt.l*2,1)];  
            [lMTk_l,vMTk_l,MA_l] = f_lMT_vMT_dM_l(qin_l,qdotin_l);    
            MA.hip_flex.l   =  MA_l(mai(1).mus.l',1);
            MA.hip_add.l    =  MA_l(mai(2).mus.l',2);
            MA.hip_rot.l    =  MA_l(mai(3).mus.l',3);
            MA.knee.l       =  MA_l(mai(4).mus.l',4);
            MA.ankle.l      =  MA_l(mai(5).mus.l',5);  
            MA.subt.l       =  MA_l(mai(6).mus.l',6); 
            % Right leg
            qin_r = [Xk_nsc(jointi.hip_flex.r*2-1,1),...
                Xk_nsc(jointi.hip_add.r*2-1,1),...
                Xk_nsc(jointi.hip_rot.r*2-1,1),...
                Xk_nsc(jointi.knee.r*2-1,1),...
                Xk_nsc(jointi.ankle.r*2-1,1),...
                Xk_nsc(jointi.subt.r*2-1,1)];  
            qdotin_r = [Xk_nsc(jointi.hip_flex.r*2,1),...
                Xk_nsc(jointi.hip_add.r*2,1),...
                Xk_nsc(jointi.hip_rot.r*2,1),...
                Xk_nsc(jointi.knee.r*2,1),...
                Xk_nsc(jointi.ankle.r*2,1),...
                Xk_nsc(jointi.subt.r*2,1)];      
            [lMTk_r,vMTk_r,MA_r] = f_lMT_vMT_dM_r(qin_r,qdotin_r);
            % Here we take the indices from left since the vector is 1:43
            MA.hip_flex.r   =  MA_r(mai(1).mus.l',1);
            MA.hip_add.r    =  MA_r(mai(2).mus.l',2);
            MA.hip_rot.r    =  MA_r(mai(3).mus.l',3);
            MA.knee.r       =  MA_r(mai(4).mus.l',4);
            MA.ankle.r      =  MA_r(mai(5).mus.l',5);
            MA.subt.r       =  MA_r(mai(6).mus.l',6);
            % Both legs
            lMTk_lr     = [lMTk_l;lMTk_r];
            vMTk_lr     = [vMTk_l;vMTk_r];   
            % Reconstruct muscle activations using synergies    
            if NSyn == 99
                asynk = ak;
            else
                alk = f_SynergyProduct(syn_wl,ak(1:NSyn,1));
                ark = f_SynergyProduct(syn_wr,ak(NSyn+1:2*NSyn,1));            
                asynk = [alk;ark];
            end                   
            % Muscle activations are combined feedback and feedforward
            a_totk = MX(NMuscles,1); 
            if spasi == 1
                a_totk(musi_noSpas) = asynk(musi_noSpas);
                a_totk(musi_Spas) = asynk(musi_Spas) + a_Ffk + a_dFfk;
            else
                a_totk = asynk;
            end
            % Get muscle-tendon forces and derive Hill-equilibrium       
            [Hilldiffk,FTk_FLV,Fcek,Fisok,~,massMk,Fpassk,~] = ...
                f_forceEquilibrium_FtildeState_all(a_totk(musi_FLV),...
                FTtildek.*scaling(mpi).mp.FTtilde',...
                dFTtildek.*scaling(mpi).mp.dFTtilde,lMTk_lr(musi_FLV),...
                vMTk_lr(musi_FLV),tensions(musi_FLV));
            FTk_noFLV = f_force_noFLV(a_totk(musi_noFLV));
            FTk = MX(NMuscles,1);
            FTk(musi_FLV) = FTk_FLV;
            if ~isempty(musi_noFLV)
                FTk(musi_noFLV) = FTk_noFLV;
            end
            % Get metabolic energy rate if in the cost function   
            if W.MER ~= 0    
                % Get muscle fiber lengths
                [~,lMtildek] = f_FiberLength_TendonForce(...
                    FTtildek.*scaling(mpi).mp.FTtilde',lMTk_lr(musi_FLV)); 
                % Get muscle fiber velocities
                [vMk,~] = f_FiberVelocity_TendonForce(...
                    FTtildek.*scaling(mpi).mp.FTtilde',...
                    dFTtildek.*scaling(mpi).mp.dFTtilde,...
                    lMTk_lr(musi_FLV),vMTk_lr(musi_FLV));
                % Get metabolic energy rate: Bhargava et al. (2004)
                [e_tot,~,~,~,~,~] = fgetMetabolicEnergySmooth2004all(...
                    a_totk(musi_FLV),a_totk(musi_FLV),lMtildek,vMk,Fcek,...
                    Fpassk,massMk,pctsts(musi_FLV),Fisok,...
                    MTParameters_FLV(1,:)',body_mass,10);
            end
            % Get passive torques
            Tau_passk.hip.flex.l    = f_PassiveTorques(k_pass.hip.flex,...
                theta.pass.hip.flex,Xk_nsc(jointi.hip_flex.l*2-1,1),...
                Xk_nsc(jointi.hip_flex.l*2,1));
            Tau_passk.hip.flex.r    = f_PassiveTorques(k_pass.hip.flex,...
                theta.pass.hip.flex,Xk_nsc(jointi.hip_flex.r*2-1,1),...
                Xk_nsc(jointi.hip_flex.r*2,1));
            Tau_passk.hip.add.l     = f_PassiveTorques(k_pass.hip.add,...
                theta.pass.hip.add,Xk_nsc(jointi.hip_add.l*2-1,1),...
                Xk_nsc(jointi.hip_add.l*2,1));
            Tau_passk.hip.add.r     = f_PassiveTorques(k_pass.hip.add,...
                theta.pass.hip.add,Xk_nsc(jointi.hip_add.r*2-1,1),...
                Xk_nsc(jointi.hip_add.r*2,1));
            Tau_passk.hip.rot.l     = f_PassiveTorques(k_pass.hip.rot,...
                theta.pass.hip.rot,Xk_nsc(jointi.hip_rot.l*2-1,1),...
                Xk_nsc(jointi.hip_rot.l*2,1));
            Tau_passk.hip.rot.r     = f_PassiveTorques(k_pass.hip.rot,...
                theta.pass.hip.rot,Xk_nsc(jointi.hip_rot.r*2-1,1),...
                Xk_nsc(jointi.hip_rot.r*2,1));
            Tau_passk.knee.l        = f_PassiveTorques(k_pass.knee,...
                theta.pass.knee,Xk_nsc(jointi.knee.l*2-1,1),...
                Xk_nsc(jointi.knee.l*2,1));
            Tau_passk.knee.r        = f_PassiveTorques(k_pass.knee,...
                theta.pass.knee,Xk_nsc(jointi.knee.r*2-1,1),...
                Xk_nsc(jointi.knee.r*2,1));
            Tau_passk.ankle.l       = f_PassiveTorques(k_pass.ankle,...
                theta.pass.ankle,Xk_nsc(jointi.ankle.l*2-1,1),...
                Xk_nsc(jointi.ankle.l*2,1));
            Tau_passk.ankle.r       = f_PassiveTorques(k_pass.ankle,...
                theta.pass.ankle,Xk_nsc(jointi.ankle.r*2-1,1),...
                Xk_nsc(jointi.ankle.r*2,1));        
            Tau_passk.subt.l       = f_PassiveTorques(k_pass.subt,...
                theta.pass.subt,Xk_nsc(jointi.subt.l*2-1,1),...
                Xk_nsc(jointi.subt.l*2,1));
            Tau_passk.subt.r       = f_PassiveTorques(k_pass.subt,...
                theta.pass.subt,Xk_nsc(jointi.subt.r*2-1,1),...
                Xk_nsc(jointi.subt.r*2,1));  
            Tau_passk.trunk.ext     = f_PassiveTorques(k_pass.trunk.ext,...
                theta.pass.trunk.ext,Xk_nsc(jointi.trunk.ext*2-1,1),...
                Xk_nsc(jointi.trunk.ext*2,1));
            Tau_passk.trunk.ben     = f_PassiveTorques(k_pass.trunk.ben,...
                theta.pass.trunk.ben,Xk_nsc(jointi.trunk.ben*2-1,1),...
                Xk_nsc(jointi.trunk.ben*2,1));
            Tau_passk.trunk.rot     = f_PassiveTorques(k_pass.trunk.rot,...
                theta.pass.trunk.rot,Xk_nsc(jointi.trunk.rot*2-1,1),...
                Xk_nsc(jointi.trunk.rot*2,1));     
            Tau_passk_all = [Tau_passk.hip.flex.l,Tau_passk.hip.flex.r,...
                Tau_passk.hip.add.l,Tau_passk.hip.add.r,...
                Tau_passk.hip.rot.l,Tau_passk.hip.rot.r,...
                Tau_passk.knee.l,Tau_passk.knee.r,Tau_passk.ankle.l,...
                Tau_passk.ankle.r,Tau_passk.subt.l,Tau_passk.subt.r,...
                Tau_passk.trunk.ext,Tau_passk.trunk.ben,...
                Tau_passk.trunk.rot]'; 
            
            % Loop over collocation points
            Xk_nsc_end          = D(1)*Xk_nsc;
            FTtildek_nsc_end    = D(1)*FTtildek_nsc;
            ak_end              = D(1)*ak;
            if spasi == 1
                a_Ffk_end = D(1)*a_Ffk;
                a_dFfk_end = D(1)*a_dFfk;
            end
            a_bk_end = D(1)*a_bk; 
            for j=1:d
                % Expression for the state derivatives at the collocation point
                xp_nsc          = C(1,j+1)*Xk_nsc;
                FTtildep_nsc    = C(1,j+1)*FTtildek_nsc;
                ap              = C(1,j+1)*ak;
                if spasi == 1
                    a_Ffp = C(1,j+1)*a_Ffk;
                    a_dFfp = C(1,j+1)*a_dFfk;
                end
                a_bp = C(1,j+1)*a_bk;
                for r=1:d
                    xp_nsc       = xp_nsc + C(r+1,j+1)*Xkj_nsc{r};
                    FTtildep_nsc = FTtildep_nsc + C(r+1,j+1)*FTtildekj_nsc{r};
                    ap           = ap + C(r+1,j+1)*akj{r};
                    if spasi == 1
                        a_Ffp = a_Ffp + C(r+1,j+1)*a_Ffkj{r};
                        a_dFfp = a_dFfp + C(r+1,j+1)*a_dFfkj{r};
                    end
                    a_bp = a_bp + C(r+1,j+1)*a_bkj{r};
                end 
                % Append collocation equations
                % Dynamic constraints are scaled using the same scale
                % factors as was used to scale the states
                % Activation dynamics (implicit formulation)  
                g       = {g{:}, (h*vAk.*scaling(mpi).mp.vA - ap)./ ...
                    scaling(mpi).mp.a};
                lbg     = [lbg; zeros(NMact,1)];
                ubg     = [ubg; zeros(NMact,1)]; 
                % Contraction dynamics (implicit formulation)               
                g       = {g{:}, (h*dFTtildek.*scaling(mpi).mp.dFTtilde - ...
                    FTtildep_nsc)./(scaling(mpi).mp.FTtilde')};
                lbg     = [lbg; zeros(NMuscles_FLV,1)];
                ubg     = [ubg; zeros(NMuscles_FLV,1)];
                if spasi == 1
                    % Spindle dynamics (explicit formulation)
                    % Force feedback            
                    da_Ffdt = f_spindleDynamics(a_Ffkj{j},...
                        FTtildek_nsc(musi_Spas_In_FLV),tauFf.all,gFf.all,...
                        bspas,threshold_Ff_gait.all);
                    g       = {g{:}, (h*da_Ffdt - a_Ffp)};
                    lbg     = [lbg; zeros(NMuscles_Spas,1)];
                    ubg     = [ubg; zeros(NMuscles_Spas,1)]; 
                    % Time derivative of force feedback feedback            
                    da_dFfdt = f_spindleDynamics(a_dFfkj{j},...
                        dFTtildek(musi_Spas_In_FLV)...
                        .*scaling(mpi).mp.dFTtilde,...
                        taudFf.all,gdFf.all,bspas,threshold_dFf_gait.all);
                    g       = {g{:}, (h*da_dFfdt - a_dFfp)};
                    lbg     = [lbg; zeros(NMuscles_Spas,1)];
                    ubg     = [ubg; zeros(NMuscles_Spas,1)];
                end
                % Skeleton dynamics (implicit formulation)        
                xj_nsc  = [...
                    Xkj_nsc{j}(2); Ak_nsc(1); Xkj_nsc{j}(4); Ak_nsc(2);...
                    Xkj_nsc{j}(6); Ak_nsc(3); Xkj_nsc{j}(8); Ak_nsc(4);...
                    Xkj_nsc{j}(10); Ak_nsc(5); Xkj_nsc{j}(12); Ak_nsc(6);...
                    Xkj_nsc{j}(14); Ak_nsc(7); Xkj_nsc{j}(16); Ak_nsc(8);...
                    Xkj_nsc{j}(18); Ak_nsc(9); Xkj_nsc{j}(20); Ak_nsc(10);...
                    Xkj_nsc{j}(22); Ak_nsc(11); Xkj_nsc{j}(24); Ak_nsc(12);...
                    Xkj_nsc{j}(26); Ak_nsc(13); Xkj_nsc{j}(28); Ak_nsc(14);...
                    Xkj_nsc{j}(30); Ak_nsc(15); Xkj_nsc{j}(32); Ak_nsc(16);...
                    Xkj_nsc{j}(34); Ak_nsc(17); Xkj_nsc{j}(36); Ak_nsc(18);...
                    Xkj_nsc{j}(38); Ak_nsc(19); Xkj_nsc{j}(40); Ak_nsc(20);...
                    Xkj_nsc{j}(42); Ak_nsc(21)];                   
                g	= {g{:}, (h*xj_nsc - xp_nsc)./(scaling(mpi).mp.QsQdots')};
                lbg	= [lbg; zeros(2*nq.res,1)];
                ubg	= [ubg; zeros(2*nq.res,1)];   
                % Trunk activation dynamics (explicit formulation)   
                dadt    = f_TrunkActivationDynamics(e_bk,a_bkj{j});
                g       = {g{:}, (h*dadt - a_bp)./scaling(mpi).mp.a_b};
                lbg     = [lbg; zeros(nq.res_trunk,1)];
                ubg     = [ubg; zeros(nq.res_trunk,1)];  
                % Add contribution to the end state
                Xk_nsc_end       = Xk_nsc_end + D(j+1)*Xkj_nsc{j};
                FTtildek_nsc_end = FTtildek_nsc_end + D(j+1)*FTtildekj_nsc{j};
                ak_end           = ak_end + D(j+1)*akj{j};
                if spasi == 1
                    a_Ffk_end = a_Ffk_end + D(j+1)*a_Ffkj{j};   
                    a_dFfk_end = a_dFfk_end + D(j+1)*a_dFfkj{j}; 
                end                 
                a_bk_end = a_bk_end + D(j+1)*a_bkj{j};  
                % Add contribution to quadrature function
                % Tracking term(s)
                % TD kinematics tracking
                Qs_costk = ...
                    B(j+1)*(f_JNq_bpty(Xk{k+1,1}(Qsi(jointi.res_bptyi))-...
                    Qs(mpi).mp.allinterpfilt(k+1,jointi.res_bptyi+1)'... 
                    ./scaling(mpi).mp.Qs(jointi.res_bptyi)'))*h;
                track_terms = W.TrackTD*Qs_costk;
                % Motor control terms
                % Muscle activations
                 a_costk = B(j+1)*(f_JNM_exp(asynk,exp_A))*h;
                % Metabolic energy rate
                if W.MER == 0
                    mE_costk = 0;
                else
                    mE_costk = B(j+1)*(f_JNM_FLVexp(e_tot,exp_E))/body_mass*h;
                end
                % Trunk excitations
                trunk_costk = B(j+1)*(f_Jnq_trunk(e_bk))*h;     
                % Joint accelerations
                Qdotdot_costk = B(j+1)*(f_JNq_all(Ak))*h;
                % Passive joint torques
                passT_costk = B(j+1)*(f_JNq_act(Tau_passk_all))*h;
                % Time derivative of muscle activations / muscle-tendon forces
                vA_costk = B(j+1)*(f_JNMact(vAk))*h;
                dFTtilde_costk = B(j+1)*(f_JNM_FLV(dFTtildek))*h;
                % All combined
                mc_terms = W.Act*a_costk + W.MER*mE_costk +...
                    W.Acc*Qdotdot_costk + W.PassT*passT_costk + ...
                    W.TrunkExc*trunk_costk + W.u*vA_costk + W.u*dFTtilde_costk;   
                % Quadrature function   
                dist_norm = dist_trav_tot;
                J = J + track_terms + 1/dist_norm*mc_terms;                  
            end                              
            % Add path constraints
            % Pelvis residuals
            g   = {g{:},(Tauk_abs)./scaling(mpi).mp.T(1)};
            % Null pelvis residuals
            lbg = [lbg; zeros(size(Tauk_abs,1),1)];        
            ubg = [ubg; zeros(size(Tauk_abs,1),1)];     
            % Muscle-driven joint torques for the lower limbs
            % Hip flexion, left
            Ft_hip_flex_l   = FTk(mai(1).mus.l',1);
            T_hip_flex_l    = f_T27(MA.hip_flex.l,Ft_hip_flex_l);
            g               = {g{:},Tk(jointi.hip_flex.l,1) - ...
                (T_hip_flex_l + Tau_passk.hip.flex.l - ...
                0.1*Xk_nsc(jointi.hip_flex.l*2,1))};        
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];    
            % Hip flexion, right
            Ft_hip_flex_r   = FTk(mai(1).mus.r',1);
            T_hip_flex_r    = f_T27(MA.hip_flex.r,Ft_hip_flex_r);
            g               = {g{:},Tk(jointi.hip_flex.r,1) - ...
                (T_hip_flex_r + Tau_passk.hip.flex.r - ...
                0.1*Xk_nsc(jointi.hip_flex.r*2,1))};        
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];   
            % Hip adduction, left
            Ft_hip_add_l    = FTk(mai(2).mus.l',1);
            T_hip_add_l     = f_T27(MA.hip_add.l,Ft_hip_add_l);
            g               = {g{:},Tk(jointi.hip_add.l,1) - ...
                (T_hip_add_l + Tau_passk.hip.add.l - ...
                0.1*Xk_nsc(jointi.hip_add.l*2,1))};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];    
            % Hip adduction, right
            Ft_hip_add_r    = FTk(mai(2).mus.r',1);
            T_hip_add_r     = f_T27(MA.hip_add.r,Ft_hip_add_r);
            g               = {g{:},Tk(jointi.hip_add.r,1) - ...
                (T_hip_add_r + Tau_passk.hip.add.r) - ...
                0.1*Xk_nsc(jointi.hip_add.r*2,1)};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];  
            % Hip rotation, left
            Ft_hip_rot_l    = FTk(mai(3).mus.l',1);
            T_hip_rot_l     = f_T27(MA.hip_rot.l,Ft_hip_rot_l);
            g               = {g{:},Tk(jointi.hip_rot.l,1) - ...
                (T_hip_rot_l + Tau_passk.hip.rot.l - ...
                0.1*Xk_nsc(jointi.hip_rot.l*2,1))};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];    
            % Hip rotation, right
            Ft_hip_rot_r    = FTk(mai(3).mus.r',1);
            T_hip_rot_r     = f_T27(MA.hip_rot.r,Ft_hip_rot_r);
            g               = {g{:},Tk(jointi.hip_rot.r,1) - ...
                (T_hip_rot_r + Tau_passk.hip.rot.r - ...
                0.1*Xk_nsc(jointi.hip_rot.r*2,1))};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];   
            % Knee, left
            Ft_knee_l       = FTk(mai(4).mus.l',1);
            T_knee_l        = f_T13(MA.knee.l,Ft_knee_l);
            g               = {g{:},Tk(jointi.knee.l,1) - ...
                (T_knee_l + Tau_passk.knee.l - ...
                0.1*Xk_nsc(jointi.knee.l*2,1))};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];    
            % Knee, right
            Ft_knee_r       = FTk(mai(4).mus.r',1);
            T_knee_r        = f_T13(MA.knee.r,Ft_knee_r);
            g               = {g{:},Tk(jointi.knee.r,1) - ...
                (T_knee_r + Tau_passk.knee.r - ...
                0.1*Xk_nsc(jointi.knee.r*2,1))};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];    
            % Ankle, left
            Ft_ankle_l      = FTk(mai(5).mus.l',1);
            T_ankle_l       = f_T12(MA.ankle.l,Ft_ankle_l);
            g               = {g{:},Tk(jointi.ankle.l,1) - ...
                (T_ankle_l + Tau_passk.ankle.l - ...
                0.1*Xk_nsc(jointi.ankle.l*2,1))};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];    
            % Ankle, right
            Ft_ankle_r      = FTk(mai(5).mus.r',1);
            T_ankle_r       = f_T12(MA.ankle.r,Ft_ankle_r);
            g               = {g{:},Tk(jointi.ankle.r,1) - ...
                (T_ankle_r + Tau_passk.ankle.r - ...
                0.1*Xk_nsc(jointi.ankle.r*2,1))};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];
            % Subtalar, left
            Ft_subt_l       = FTk(mai(6).mus.l',1);
            T_subt_l        = f_T12(MA.subt.l,Ft_subt_l);
            g               = {g{:},(Tk(jointi.subt.l,1) - ...
                (T_subt_l + Tau_passk.subt.l - ...
                0.1*Xk_nsc(jointi.subt.l*2,1)))};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];    
            % Subtalar, right
            Ft_subt_r       = FTk(mai(6).mus.r',1);
            T_subt_r        = f_T12(MA.subt.r,Ft_subt_r);
            g               = {g{:},(Tk(jointi.subt.r,1) - ...
                (T_subt_r + Tau_passk.subt.r - ...
                0.1*Xk_nsc(jointi.subt.r*2,1)))};
            lbg             = [lbg; 0];
            ubg             = [ubg; 0];
            % Torque-driven joint torques for the trunk
            % Trunk
            g       = ...
                {g{:},Tk(jointi.res_trunki,1)/scaling(mpi).mp.TrunkTau - a_bk};
            lbg     = [lbg; zeros(nq.res_trunk,1)];
            ubg 	= [ubg; zeros(nq.res_trunk,1)];
            % Activation dynamics (implicit formulation)            
            act1 = vAk*scaling(mpi).mp.vA + ak./(ones(size(ak,1),1)*tdeact);
            act2 = vAk*scaling(mpi).mp.vA + ak./(ones(size(ak,1),1)*tact);
            % act1
            g               = {g{:},act1};
            lbg             = [lbg; zeros(NMact,1)];
            ubg             = [ubg; inf*ones(NMact,1)]; 
            % act2
            g               = {g{:},act2};
            lbg             = [lbg; -inf*ones(NMact,1)];
            ubg             = [ubg; ones(NMact,1)./(ones(NMact,1)*tact)];
            % Contraction dynamics (implicit formulation)
            g               = {g{:},Hilldiffk};
            lbg             = [lbg; zeros(NMuscles_FLV,1)];
            ubg             = [ubg; zeros(NMuscles_FLV,1)];      
            % Total activations (feedback and feedforward) cannot exceed 1
            if spasi == 1
                g   = {g{:},a_totk(musi_Spas)};
                lbg = [lbg; zeros(NMuscles_Spas,1)];
                ubg = [ubg; ones(NMuscles_Spas,1)]; 
            end
            % Constraints to prevent parts of the skeleton to penetrate each
            % other.
            % Origins calcaneus (transv plane) at min 9 cm from each other.
            g       = {g{:},sumsqr(Tk(calcOr.r,1) - Tk(calcOr.l,1))};
            lbg     = [lbg; 0.0081];
            ubg     = [ubg; 4];   
            % Origins tibia (transv plane) at min 11 cm from each other.   
            g       = {g{:},sumsqr(Tk(tibiaOr.r,1) - Tk(tibiaOr.l,1))};
            lbg     = [lbg; 0.0121];
            ubg     = [ubg; 4]; 
            % New NLP variables for states at end of interval
            if k ~= N-1
                % Muscle activations
                ak      = MX.sym(['a_' num2str(k+1)],NMact);
                w       = {w{:}, ak};
                lbw 	= [lbw; bounds(mpi).mp.a.lower'];
                ubw 	= [ubw; bounds(mpi).mp.a.upper'];
                w0      = [w0;  guess(mpi).mp.a(k+2,:)'];
                % Muscle-tendon forces
                FTtildek	= MX.sym(['FTtilde_' num2str(k+1)],NMuscles_FLV);
                w           = {w{:}, FTtildek};
                lbw         = [lbw; bounds(mpi).mp.FTtilde.lower'];
                ubw         = [ubw; bounds(mpi).mp.FTtilde.upper'];
                w0          = [w0;  guess(mpi).mp.FTtilde(k+2,:)'];    
                % Qs and Qdots
                w = {w{:}, Xk{k+2,1}};
                lbw = [lbw; bounds(mpi).mp.QsQdots.lower(k+2,:)'];
                ubw = [ubw; bounds(mpi).mp.QsQdots.upper(k+2,:)'];
                w0 = [w0;  guess(mpi).mp.QsQdots(k+2,:)'];
                if spasi == 1
                    % Muscle activations from muscle-tendon force feedback
                    a_Ffk	= MX.sym(['a_Ff_' num2str(k+1)],NMuscles_Spas);
                    w       = {w{:}, a_Ffk};
                    lbw 	= [lbw; bounds(mpi).mp.a_Ff.lower'];
                    ubw 	= [ubw; bounds(mpi).mp.a_Ff.upper'];
                    w0  	= [w0;  guess(mpi).mp.a_Ff(k+2,:)'];
                    % Muscle activations from muscle-tendon force rate
                    % feedback
                    a_dFfk	= MX.sym(['a_dFf_' num2str(k+1)],NMuscles_Spas);
                    w       = {w{:}, a_dFfk};
                    lbw  	= [lbw; bounds(mpi).mp.a_dFf.lower'];
                    ubw  	= [ubw; bounds(mpi).mp.a_dFf.upper'];
                    w0    	= [w0;  guess(mpi).mp.a_dFf(k+2,:)'];
                end
                % Trunk activations
                a_bk	= MX.sym(['a_b_' num2str(k+1)],nq.res_trunk);
                w       = {w{:}, a_bk};
                lbw     = [lbw; bounds(mpi).mp.a_b.lower'];
                ubw     = [ubw; bounds(mpi).mp.a_b.upper'];
                w0      = [w0;  guess(mpi).mp.a_b(k+2,:)'];
            else
                % Muscle activations
                ak      = MX.sym(['a_' num2str(k+1)],NMact);
                w       = {w{:}, ak};
                lbw     = [lbw; bounds(mpi).mp.a.lower'];
                ubw     = [ubw; bounds(mpi).mp.a.upper'];
                w0      = [w0;  guess(mpi).mp.a(end,:)'];
                % Muscle-tendon forces
                FTtildek	= MX.sym(['FTtilde_' num2str(k+1)],NMuscles_FLV);
                w           = {w{:}, FTtildek};
                lbw         = [lbw; bounds(mpi).mp.FTtilde.lower'];
                ubw         = [ubw; bounds(mpi).mp.FTtilde.upper'];
                w0          = [w0;  guess(mpi).mp.FTtilde(end,:)'];    
                % Qs and Qdots
                w = {w{:}, Xk{k+2,1}};
                lbw	= [lbw; bounds(mpi).mp.QsQdots.lower(end,:)'];
                ubw = [ubw; bounds(mpi).mp.QsQdots.upper(end,:)'];
                w0 = [w0;  guess(mpi).mp.QsQdots(end,:)'];
                if spasi == 1
                    % Muscle activations from muscle-tendon force feedback
                    a_Ffk	= MX.sym(['a_Ff_' num2str(k+1)],NMuscles_Spas);
                    w       = {w{:}, a_Ffk};
                    lbw   	= [lbw; bounds(mpi).mp.a_Ff.lower'];
                    ubw    	= [ubw; bounds(mpi).mp.a_Ff.upper'];
                    w0    	= [w0;  guess(mpi).mp.a_Ff(end,:)'];
                    % Muscle activations from muscle-tendon force rate
                    % feedback
                    a_dFfk 	= MX.sym(['a_dFf_' num2str(k+1)],NMuscles_Spas);
                    w     	= {w{:}, a_dFfk};
                    lbw    	= [lbw; bounds(mpi).mp.a_dFf.lower'];
                    ubw    	= [ubw; bounds(mpi).mp.a_dFf.upper'];
                    w0     	= [w0;  guess(mpi).mp.a_dFf(end,:)'];
                end
                % Trunk activations
                a_bk	= MX.sym(['a_b_' num2str(k+1)], nq.res_trunk);
                w       = {w{:}, a_bk};
                lbw     = [lbw; bounds(mpi).mp.a_b.lower'];
                ubw     = [ubw; bounds(mpi).mp.a_b.upper'];
                w0      = [w0;  guess(mpi).mp.a_b(end,:)'];
            end
            % Rescale variables to impose equality constraints
            Xk_end = (Xk_nsc_end)./scaling(mpi).mp.QsQdots';
            FTtildek_end = (FTtildek_nsc_end)./scaling(mpi).mp.FTtilde';
            % Add equality constraints (next interval starts with end values of 
            % states from previous interval)
            g   = {g{:}, Xk_end-Xk{k+2,1}, FTtildek_end-FTtildek, ak_end-ak};
            lbg = [lbg; zeros(2*nq.res + NMact + NMuscles_FLV,1)];
            ubg = [ubg; zeros(2*nq.res + NMact + NMuscles_FLV,1)]; 
            if spasi == 1
                g   = {g{:}, a_Ffk_end-a_Ffk, a_dFfk_end-a_dFfk};
                lbg = [lbg; zeros(NMuscles_Spas + NMuscles_Spas,1)];
                ubg = [ubg; zeros(NMuscles_Spas + NMuscles_Spas,1)];   
            end
            g   = {g{:}, a_bk_end-a_bk};
            lbg = [lbg; zeros(nq.res_trunk,1)];
            ubg = [ubg; zeros(nq.res_trunk,1)];
        end   
        % Periodic constraints
        % Periodicity on joint states
        % All states but Q pelvis_tx
        idx_per = [2*jointi.pelvis.list-1:2*jointi.pelvis.tilt,...
            2*jointi.pelvis.tx:2*jointi.trunk.rot];
        g   = {g{:}, Xk_end(idx_per)-X0(idx_per,1)};
        lbg = [lbg; zeros(length(idx_per),1)];
        ubg = [ubg; zeros(length(idx_per),1)];        
        % Periodicity on muscle states
        % Muscle activations
        g   = {g{:}, ak_end-a0};
        lbg = [lbg; zeros(NMact,1)];
        ubg = [ubg; zeros(NMact,1)];
        % Muscle-tendon forces
        g   = {g{:}, FTtildek_end-FTtilde0};
        lbg = [lbg; zeros(NMuscles_FLV,1)];
        ubg = [ubg; zeros(NMuscles_FLV,1)];
        % Speed constraint
        vel_aver_tot = dist_trav_tot/tf;
        g   = {g{:}, vel_aver_tot - gaitSpeed};
        lbg = [lbg; 0];
        ubg = [ubg; 0];            
    end    
    % Assert bounds / IG
    % Lower bounds smaller than upper bounds
    assert_bw = isempty(find(lbw <= ubw == 0,1));
    assert_bg = isempty(find(lbg <= ubg == 0,1));
    % Design varibles between -1 and 1    
    assert_bwl = isempty(find(lbw < -1 == 1,1));
    assert_bwu = isempty(find(1 < ubw == 1,1));
    % Initial guess within bounds
    assert_w0_ubw = isempty(find(w0 <= ubw == 0,1));
    assert_w0_lbw = isempty(find(lbw <= w0 == 0,1)); 
                
    % Create an NLP solver
    prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy      = 'adaptive';
    options.ipopt.max_iter = 20000;
    options.ipopt.tol = 1*10^(-tol_ipopt);
    solver = nlpsol('solver', 'ipopt', prob, options);
    % Create and save diary    
    p = mfilename('fullpath');
    [~,namescript,~] = fileparts(p);
    pathresults = [pathPredictiveSimulations,'/Results'];
    if ~(exist([pathresults,'/',namescript],'dir')==7)
        mkdir(pathresults,namescript);
    end
    if (exist([pathresults,'/',namescript,'/D',savename],'file')==2)
        delete ([pathresults,'/',namescript,'/D',savename])
    end 
    diary([pathresults,'/',namescript,'/D',savename]); 
    % Solve problem
    sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
        'lbg', lbg, 'ubg', ubg);    
    diary off
    % Extract results
    w_opt = full(sol.x);
    g_opt = full(sol.g);  
    % Save results
    setup.tolerance.ipopt = tol_ipopt;
    for mpi = 1:NPhases  
        setup.bounds(mpi).mp = bounds(mpi).mp;
        setup.scaling(mpi).mp = scaling(mpi).mp;
        setup.guess(mpi).mp = guess(mpi).mp;
    end
    setup.lbw = lbw;
    setup.ubw = ubw;
    % Extract stats
    stats = solver.stats();
    % Save results and setup
    save([pathresults,'/',namescript,'/w',savename],'w_opt');
    save([pathresults,'/',namescript,'/g',savename],'g_opt');
    save([pathresults,'/',namescript,'/s',savename],'setup');
    save([pathresults,'/',namescript,'/stats',savename],'stats');
end

%% Analyze results
if analyseResults
    %% Load results
    if loadResults
        p = mfilename('fullpath');
        [~,namescript,~] = fileparts(p);
        pathresults = [pathPredictiveSimulations,'/Results'];
        load([pathresults,'/',namescript,'/w',savename]);
        load([pathresults,'/',namescript,'/g',savename]);
        load([pathresults,'/',namescript,'/s',savename]);
        load([pathresults,'/',namescript,'/stats',savename]);
    end
    
    %% Extract results
    % All optimized design variables are saved in a single column vector      
    % Number of design variables  
    NControls = NMact+NMuscles_FLV+nq.res+nq.res_trunk;
    NStates = NMact+NMuscles_FLV+2*nq.res+nq.res_trunk;
    NParameters = 1;
    if NSyn ~= 99        
        NParameters = NParameters + NMuscles*NSyn;   
    end
    if spasi == 1
        NStates = NStates+2*NMuscles_Spas;
    end    
    % In the loop
    Nwl = NControls+d*NStates+NStates;
    % In total
    Nw = NParameters+NStates+N*Nwl;
    % Before the variable corresponding to the first collocation point
    Nwm = NParameters+NStates+NControls;
    % Here we extract the results and re-organize them for analysis  
    % Static parameters  
    np_acc = 1;
    tf_opt = w_opt(1);
    np_acc = np_acc + 1;
    if NSyn ~= 99     
        syn_wl_opt = w_opt(np_acc:np_acc+NMuscles/2*NSyn-1);
        syn_wr_opt = w_opt(np_acc+NMuscles/2*NSyn:NParameters);
    end    
    Nw_acc = 0;    
    a_opt = struct('mp',[]); 
    FTtilde_opt = struct('mp',[]); 
    q_opt = struct('mp',[]); 
    qdot_opt = struct('mp',[]); 
    qqdot_opt = struct('mp',[]); 
    a_b_opt = struct('mp',[]); 
    e_b_opt = struct('mp',[]); 
    vA_opt = struct('mp',[]); 
    dFTtilde_opt = struct('mp',[]); 
    qdotdot_opt = struct('mp',[]); 
    a_opt_ext = struct('mp',[]);
    a_opt_ext_col = struct('mp',[]);
    FTtilde_opt_ext = struct('mp',[]);
    q_opt_ext = struct('mp',[]);
    q_dot_opt_ext = struct('mp',[]);
    a_Ff_opt = struct('mp',[]);
    a_dFf_opt = struct('mp',[]);
    tempi = Nw_acc+NParameters;
    for mpi = 1:NPhases 
        % Mesh points
        % Muscle activations
        a_opt(mpi).mp = zeros(N+1,NMact);
        for i = 1:NMact
            a_opt(mpi).mp(:,i) = w_opt(tempi+i:Nwl:Nw_acc+Nw);
        end
        tempi = tempi + NMact;
        % Muscle-tendon forces
        FTtilde_opt(mpi).mp = zeros(N+1,NMuscles_FLV);
        for i = 1:NMuscles_FLV
            FTtilde_opt(mpi).mp(:,i) = w_opt(tempi+i:Nwl:Nw_acc+Nw);
        end    
        tempi = tempi + NMuscles_FLV;
        % Qs and Qdots
        q_opt(mpi).mp = zeros(N+1,nq.res);
        qdot_opt(mpi).mp = zeros(N+1,nq.res);
        count = 0;
        for i = 1:2:2*nq.res
            count = count +1;
            q_opt(mpi).mp(:,count) = w_opt(tempi+i:Nwl:Nw_acc+Nw);
            qdot_opt(mpi).mp(:,count) = w_opt(tempi+i+1:Nwl:Nw_acc+Nw);
        end
        tempi = tempi + 2*nq.res;
        qqdot_opt(mpi).mp = zeros(N+1,2*nq.res); 
        qqdot_opt(mpi).mp(:,1:2:end) = q_opt(mpi).mp;
        qqdot_opt(mpi).mp(:,2:2:end) = qdot_opt(mpi).mp;
        if spasi == 1
            % Muscle activations from muscle-tendon force feedback and 
            % muscle-tendon force rate feedback
            a_Ff_opt(mpi).mp = zeros(N+1,NMuscles_Spas); 
            a_dFf_opt(mpi).mp = zeros(N+1,NMuscles_Spas);
            for i = 1:NMuscles_Spas
                a_Ff_opt(mpi).mp(:,i) = w_opt(tempi+i:Nwl:Nw_acc+Nw);
                a_dFf_opt(mpi).mp(:,i) = ...
                    w_opt(tempi+NMuscles_Spas+i:Nwl:Nw_acc+Nw);
            end 
            tempi = tempi + 2*NMuscles_Spas;            
        end
        % Trunk activations
        a_b_opt(mpi).mp = zeros(N+1,nq.res_trunk);
        for i = 1:nq.res_trunk
            a_b_opt(mpi).mp(:,i) = w_opt(tempi+i:Nwl:Nw_acc+Nw);
        end     
        tempi = tempi + nq.res_trunk; 
        % Time derivative of muscle activations
        vA_opt(mpi).mp = zeros(N,NMact);
        for i = 1:NMact
            vA_opt(mpi).mp(:,i) = w_opt(tempi+i:Nwl:Nw_acc+Nw);
        end
        tempi = tempi + NMact; 
        % Time derivative of muscle-tendon forces
        dFTtilde_opt(mpi).mp = zeros(N,NMuscles_FLV);
        for i = 1:NMuscles_FLV
            dFTtilde_opt(mpi).mp(:,i) = w_opt(tempi+i:Nwl:Nw_acc+Nw);
        end
        tempi = tempi + NMuscles_FLV; 
        % Time derivative of joint velocities
        qdotdot_opt(mpi).mp = zeros(N,nq.res);
        for i = 1:nq.res
            qdotdot_opt(mpi).mp(:,i) = w_opt(tempi+i:Nwl:Nw_acc+Nw);
        end
        tempi = tempi + nq.res;
        % Trunk excitations
        e_b_opt(mpi).mp = zeros(N,nq.res_trunk);
        for i = 1:nq.res_trunk
            e_b_opt(mpi).mp(:,i) = w_opt(tempi+i:Nwl:Nw_acc+Nw); 
        end                
        % Collocation points
        % Muscle activations
        a_opt_ext(mpi).mp=zeros(N*(d+1)+1,NMact);
        a_opt_ext(mpi).mp(1:(d+1):end,:)= a_opt(mpi).mp;
        for nmusi=1:NMact
            a_opt_ext(mpi).mp(2:(d+1):end,nmusi) = ...
                w_opt(Nw_acc+Nwm+nmusi:Nwl:Nw_acc+Nw);
            a_opt_ext(mpi).mp(3:(d+1):end,nmusi) = ...
                w_opt(Nw_acc+Nwm+NMact+nmusi:Nwl:Nw_acc+Nw);
            a_opt_ext(mpi).mp(4:(d+1):end,nmusi) = ...
                w_opt(Nw_acc+Nwm+NMact+NMact+nmusi:Nwl:Nw_acc+Nw);
        end
        % Muscle activations at collocation points only
        a_opt_ext_col(mpi).mp = zeros(N*d,NMact); 
        for nmusi=1:NMact
            a_opt_ext_col(mpi).mp(1:d:end,nmusi) = ...
                w_opt(Nw_acc+Nwm+nmusi:Nwl:Nw_acc+Nw);
            a_opt_ext_col(mpi).mp(2:d:end,nmusi) = ...
                w_opt(Nw_acc+Nwm+NMact+nmusi:Nwl:Nw_acc+Nw);
            a_opt_ext_col(mpi).mp(3:d:end,nmusi) = ...
                w_opt(Nw_acc+Nwm+NMact+NMact+nmusi:Nwl:Nw_acc+Nw);   
        end 
        % Muscle-tendon forces
        FTtilde_opt_ext(mpi).mp=zeros(N*(d+1)+1,NMuscles_FLV);
        FTtilde_opt_ext(mpi).mp(1:(d+1):end,:)= FTtilde_opt(mpi).mp;
        for nmusi=1:NMuscles_FLV
            FTtilde_opt_ext(mpi).mp(2:(d+1):end,nmusi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+nmusi:Nwl:Nw_acc+Nw);
            FTtilde_opt_ext(mpi).mp(3:(d+1):end,nmusi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+NMuscles_FLV+nmusi:Nwl:Nw_acc+Nw);
            FTtilde_opt_ext(mpi).mp(4:(d+1):end,nmusi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+NMuscles_FLV+NMuscles_FLV+nmusi:Nwl:...
                Nw_acc+Nw);
        end
        % Qs and Qdots
        q_opt_ext(mpi).mp=zeros(N*(d+1)+1,nq.res);
        q_opt_ext(mpi).mp(1:(d+1):end,:)= q_opt(mpi).mp;
        q_dot_opt_ext(mpi).mp=zeros(N*(d+1)+1,nq.res);
        q_dot_opt_ext(mpi).mp(1:(d+1):end,:)= qdot_opt(mpi).mp;
        nqi_col = 1:2:2*nq.res;
        for nqi=1:nq.res
            nqi_q = nqi_col(nqi);
            q_opt_ext(mpi).mp(2:(d+1):end,nqi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+d*NMuscles_FLV+nqi_q:Nwl:Nw_acc+Nw);   
            q_opt_ext(mpi).mp(3:(d+1):end,nqi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+d*NMuscles_FLV+2*nq.res+nqi_q:Nwl:...
                Nw_acc+Nw);  
            q_opt_ext(mpi).mp(4:(d+1):end,nqi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+d*NMuscles_FLV+2*nq.res+2*nq.res+...
                nqi_q:Nwl:Nw_acc+Nw);  
            q_dot_opt_ext(mpi).mp(2:(d+1):end,nqi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+d*NMuscles_FLV+nqi_q+1:Nwl:Nw_acc+Nw);   
            q_dot_opt_ext(mpi).mp(3:(d+1):end,nqi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+d*NMuscles_FLV+2*nq.res+nqi_q+1:Nwl:...
                Nw_acc+Nw);  
            q_dot_opt_ext(mpi).mp(4:(d+1):end,nqi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+d*NMuscles_FLV+2*nq.res+2*nq.res+...
                nqi_q+1:Nwl:Nw_acc+Nw);
        end
        Nw_acc = Nw_acc + Nw - NParameters;
    end
    
    %% Unscale results
    q_opt_unsc = struct('mp',[]);
    q_opt_unsc_all = struct('mp',[]);
    qdot_opt_unsc = struct('mp',[]);
    qdot_opt_unsc_all = struct('mp',[]);
    a_opt_unsc = struct('mp',[]);
    FTtilde_opt_unsc = struct('mp',[]);
    qdotdot_opt_unsc = struct('mp',[]);
    vA_opt_unsc = struct('mp',[]);
    e_opt_unsc = struct('mp',[]);
    dFTtilde_opt_unsc = struct('mp',[]);
    a_syn_opt = struct('mp',[]);
    a_Ff_opt_unsc = struct('mp',[]);
    a_dFf_opt_unsc = struct('mp',[]);
    a_tot_opt = struct('mp',[]);
    a_b_opt_unsc = struct('mp',[]);
    e_b_opt_unsc = struct('mp',[]);
    for mpi = 1:NPhases
        % States at mesh points
        % Qs (1:N-1)
        q_opt_unsc(mpi).mp.rad = q_opt(mpi).mp(1:end-1,:)...
            .*repmat(scaling(mpi).mp.Qs,size(q_opt(mpi).mp(1:end-1,:),1),1); 
        % Convert in degrees
        q_opt_unsc(mpi).mp.deg = q_opt_unsc(mpi).mp.rad;
        q_opt_unsc(mpi).mp.deg(:,jointi.res_roti) = ...
            q_opt_unsc(mpi).mp.deg(:,jointi.res_roti).*180/pi;
        % Qs (1:N)
        q_opt_unsc_all(mpi).mp.rad = q_opt(mpi).mp(1:end,:)...
            .*repmat(scaling(mpi).mp.Qs,size(q_opt(mpi).mp(1:end,:),1),1); 
        % Convert in degrees
        q_opt_unsc_all(mpi).mp.deg = q_opt_unsc_all(mpi).mp.rad;
        q_opt_unsc_all(mpi).mp.deg(:,jointi.res_roti) = ...
            q_opt_unsc_all(mpi).mp.deg(:,jointi.res_roti).*180/pi;
        % Qdots (1:N-1)
        qdot_opt_unsc(mpi).mp.rad = qdot_opt(mpi).mp(1:end-1,:).*repmat(...
            scaling(mpi).mp.Qdots,size(qdot_opt(mpi).mp(1:end-1,:),1),1);
        % Convert in degrees
        qdot_opt_unsc(mpi).mp.deg = qdot_opt_unsc(mpi).mp.rad;
        qdot_opt_unsc(mpi).mp.deg(:,jointi.res_roti) = ...
            qdot_opt_unsc(mpi).mp.deg(:,jointi.res_roti).*180/pi;
        % Qdots (1:N)
        qdot_opt_unsc_all(mpi).mp.rad = qdot_opt(mpi).mp(1:end,:).*repmat(...
            scaling(mpi).mp.Qdots,size(qdot_opt(mpi).mp(1:end,:),1),1); 
        % Muscle /synergy activations
        a_syn_opt(mpi).mp = a_opt(mpi).mp(1:end-1,:)...
            .*repmat(scaling(mpi).mp.a,size(a_opt(mpi).mp(1:end-1,:),1),...
            size(a_opt(mpi).mp,2));
        % Muscle-tendon forces
        FTtilde_opt_unsc(mpi).mp = FTtilde_opt(mpi).mp(1:end-1,:)...
            .*repmat(scaling(mpi).mp.FTtilde,...
            size(FTtilde_opt(mpi).mp(1:end-1,:),1),1);
        % Trunk activations
        a_b_opt_unsc(mpi).mp = a_b_opt(mpi).mp(1:end-1,:)...
            .*repmat(scaling(mpi).mp.a_b,size(a_b_opt(mpi).mp(1:end-1,:),1),...
            size(a_b_opt(mpi).mp,2));
        % Controls at mesh points
        % Time derivative of Qdots
        qdotdot_opt_unsc(mpi).mp.rad = qdotdot_opt(mpi).mp...
            .*repmat(scaling(mpi).mp.Qdotdots,size(qdotdot_opt(mpi).mp,1),1);        
        % Convert in degrees
        qdotdot_opt_unsc(mpi).mp.deg = qdotdot_opt_unsc(mpi).mp.rad;
        qdotdot_opt_unsc(mpi).mp.deg(:,jointi.res_roti) = ...
            qdotdot_opt_unsc(mpi).mp.deg(:,jointi.res_roti).*180/pi;
        % Trunk excitations
        e_b_opt_unsc(mpi).mp = e_b_opt(mpi).mp.*repmat(scaling(mpi).mp.e_b,...
            size(e_b_opt(mpi).mp,1),size(e_b_opt(mpi).mp,2));
        % Time derivative of muscle /synergy activations (states)
        vA_opt_unsc(mpi).mp = vA_opt(mpi).mp.*repmat(scaling(mpi).mp.vA,...
            size(vA_opt(mpi).mp,1),size(vA_opt(mpi).mp,2));
        % Time derivative of muscle-tendon forces (states)
        dFTtilde_opt_unsc(mpi).mp = dFTtilde_opt(mpi).mp...
            .*repmat(scaling(mpi).mp.dFTtilde,size(dFTtilde_opt(mpi).mp,1),...
            size(dFTtilde_opt(mpi).mp,2));
        
        % Reconstruct muscle activations from synergies
        if NSyn == 99 % no synergies
            a_opt_unsc(mpi).mp = a_syn_opt(mpi).mp;
        else
            a_syn_opt_l = zeros(N,NMuscles/2);
            a_syn_opt_r = zeros(N,NMuscles/2);
            for n = 1:N
                a_syn_opt_l(n,:) = full(f_SynergyProduct(syn_wl_opt,...
                    a_syn_opt(mpi).mp(n,1:NSyn)))';
                a_syn_opt_r(n,:) = full(f_SynergyProduct(syn_wr_opt,...
                    a_syn_opt(mpi).mp(n,NSyn+1:2*NSyn)))';  
            end
            a_opt_unsc(mpi).mp = [a_syn_opt_l,a_syn_opt_r];   
        end        
        
        % Muscle / synergy excitations
        e_opt_unsc(mpi).mp = computeExcitationRaasch(a_syn_opt(mpi).mp,...
            vA_opt_unsc(mpi).mp,ones(1,NMact)*tdeact,ones(1,NMact)*tact);
        
        % Combination of feedback and feedforward activations 
        a_tot_opt(mpi).mp = a_opt_unsc(mpi).mp;         
        if spasi == 1
            % Muscle activations from muscle-tendon force feedback and
            % of muscle-tendon force rate feedback
            a_Ff_opt_unsc(mpi).mp = a_Ff_opt(mpi).mp(1:end-1,:);
            a_dFf_opt_unsc(mpi).mp = a_dFf_opt(mpi).mp(1:end-1,:);        
            a_tot_opt(mpi).mp(:,musi_Spas) = ...
                a_tot_opt(mpi).mp(:,musi_Spas) + a_Ff_opt_unsc(mpi).mp + ...
                a_dFf_opt_unsc(mpi).mp; 
        end
    end
    
    %% Time grid  
    tgrid = struct('mp',[]);
    tgrid_ext = struct('mp',[]);
    for mpi = 1:NPhases
        % Mesh points
        tgrid(mpi).mp = linspace(0,tf_opt,N+1);
        dtime = zeros(1,d+1);
        for i=1:4
            dtime(i)=tau_root(i)*(tf_opt/N);
        end
        % Mesh points and collocation points
        tgrid_ext(mpi).mp = zeros(1,(d+1)*N+1);
        for i=1:N
            tgrid_ext(mpi).mp(((i-1)*4+1):1:i*4)=tgrid(mpi).mp(i)+dtime;
        end
        tgrid_ext(mpi).mp(end)=tf_opt;
    end
 
    %% Joint torques, ground reaction forces and torques at optimal solution
    Xk_Qs_Qdots_opt = struct('mp',[]);
    Xk_Qdotdots_opt = struct('mp',[]);
    out_res_opt = struct('mp',[]);
    Tauk_out = struct('mp',[]);
    GRF_opt_unsc = struct('mp',[]);
    GRM_opt_unsc = struct('mp',[]);
    StrideLength_opt = struct('mp',[]);
    StepWidth_opt_mean = struct('mp',[]);
    StepWidth_opt_std = struct('mp',[]);
    for mpi = 1:NPhases
        Xk_Qs_Qdots_opt(mpi).mp = zeros(N,2*size(q_opt_unsc(mpi).mp.rad,2));
        Xk_Qs_Qdots_opt(mpi).mp(:,1:2:end) = q_opt_unsc(mpi).mp.rad;
        Xk_Qs_Qdots_opt(mpi).mp(:,2:2:end) = qdot_opt_unsc(mpi).mp.rad;
        Xk_Qdotdots_opt(mpi).mp = qdotdot_opt_unsc(mpi).mp.rad;  
        for i = 1:N
            out2_res = F([Xk_Qs_Qdots_opt(mpi).mp(i,:)';...
                Xk_Qdotdots_opt(mpi).mp(i,:)']);        
            out_res_opt(mpi).mp(i,:) = full(out2_res);  
        end
        % Optimal joint torques, ground reaction forces and moments
        Tauk_out(mpi).mp = out_res_opt(mpi).mp(:,jointi.resi);
        GRF_opt_unsc(mpi).mp = out_res_opt(mpi).mp(:,GRFi.all);
        GRM_opt_unsc(mpi).mp = out_res_opt(mpi).mp(:,GRMi.all);
        % Stride length       
        dist = sqrt(sumsqr(out_res_opt(mpi).mp(end,calcOr.r)-...
            out_res_opt(mpi).mp(1,calcOr.r)));
        StrideLength_opt(mpi).mp = full(dist);
        % Step width 
        StepWidth_opt = full(abs(out_res_opt(mpi).mp(:,calcOr.r(2)) - ...
            out_res_opt(mpi).mp(:,calcOr.l(2))));
        StepWidth_opt_mean(mpi).mp = mean(StepWidth_opt);
        StepWidth_opt_std(mpi).mp = std(StepWidth_opt);
        % assertTrunkTmax should be < 1*10^(-tol_ipopt)
        assertTrunkTmax = ...
            max(max(abs(out_res_opt(mpi).mp(:,jointi.res_trunki)-...
            (a_b_opt_unsc(mpi).mp)*scaling(mpi).mp.TrunkTau)));  
        if assertTrunkTmax > 1*10^(-tol_ipopt)
            disp('Issue when reconstructing residual forces')
        end
    end
       
    %% Passive joint torques at optimal solution
    TPass_opt_all = struct('mp',[]);
    TPass_opt = struct('mp',[]);
    for mpi = 1:NPhases
        TPass_opt_all(mpi).mp = zeros(N,nq.res_bgp);
        for i = 1:N    
            TPass_opt(mpi).mp.hip.flex.l = f_PassiveTorques(...
                k_pass.hip.flex,theta.pass.hip.flex,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.l*2));
            TPass_opt(mpi).mp.hip.flex.r = f_PassiveTorques(...
                k_pass.hip.flex,theta.pass.hip.flex,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.r*2));
            TPass_opt(mpi).mp.hip.add.l = f_PassiveTorques(...
                k_pass.hip.add,theta.pass.hip.add,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.l*2));
            TPass_opt(mpi).mp.hip.add.r = f_PassiveTorques(...
                k_pass.hip.add,theta.pass.hip.add,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.r*2));
            TPass_opt(mpi).mp.hip.rot.l = f_PassiveTorques(...
                k_pass.hip.rot,theta.pass.hip.rot,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.l*2));
            TPass_opt(mpi).mp.hip.rot.r = f_PassiveTorques(...
                k_pass.hip.rot,theta.pass.hip.rot,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.r*2));
            TPass_opt(mpi).mp.knee.l = f_PassiveTorques(...
                k_pass.knee,theta.pass.knee,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.l*2));
            TPass_opt(mpi).mp.knee.r = f_PassiveTorques(...
                k_pass.knee,theta.pass.knee,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.r*2));
            TPass_opt(mpi).mp.ankle.l = f_PassiveTorques(...
                k_pass.ankle,theta.pass.ankle,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.l*2));
            TPass_opt(mpi).mp.ankle.r = f_PassiveTorques(...
                k_pass.ankle,theta.pass.ankle,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.r*2));
            TPass_opt(mpi).mp.subt.l = f_PassiveTorques(...
                k_pass.subt,theta.pass.subt,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.l*2));
            TPass_opt(mpi).mp.subt.r = f_PassiveTorques(...
                k_pass.subt,theta.pass.subt,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.r*2));   
            TPass_opt(mpi).mp.trunk.ext = f_PassiveTorques(...
                k_pass.trunk.ext,theta.pass.trunk.ext,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.ext*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.ext*2));
            TPass_opt(mpi).mp.trunk.ben = f_PassiveTorques(...
                k_pass.trunk.ben,theta.pass.trunk.ben,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.ben*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.ben*2));
            TPass_opt(mpi).mp.trunk.rot = f_PassiveTorques(...
                k_pass.trunk.rot,theta.pass.trunk.rot,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.rot*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.rot*2)); 
            TPass_opt_all(mpi).mp(i,:) = full([...
                TPass_opt(mpi).mp.hip.flex.l,TPass_opt(mpi).mp.hip.add.l,...
                TPass_opt(mpi).mp.hip.rot.l,TPass_opt(mpi).mp.hip.flex.r,...
                TPass_opt(mpi).mp.hip.add.r,TPass_opt(mpi).mp.hip.rot.r,...
                TPass_opt(mpi).mp.knee.l,TPass_opt(mpi).mp.knee.r,...
                TPass_opt(mpi).mp.ankle.l,TPass_opt(mpi).mp.ankle.r,...
                TPass_opt(mpi).mp.subt.l,TPass_opt(mpi).mp.subt.r,...
                TPass_opt(mpi).mp.trunk.ext,TPass_opt(mpi).mp.trunk.ben,...
                TPass_opt(mpi).mp.trunk.rot]); 
        end        
    end
    
    %% Reconstruct cost function and compute metabolic cost
    J_opt           = 0;
    Qs_cost         = 0;
    a_cost          = 0;
    Qdotdots_cost   = 0;
    vA_cost         = 0;
    dFTtilde_cost   = 0;
    mE_cost         = 0;
    passT_cost      = 0;
    trunk_cost      = 0;
    e_mo_opt        = struct('mp',[]);
    e_mo_opt_tr     = struct('mp',[]);
    COT_opt         = struct('mp',[]);
    for mpi = 1:NPhases        
        count = 1;
        h_opt = tf_opt/N;  
        % Assert speed
        dist_trav_opt = Xk_Qs_Qdots_opt(mpi).mp(end,2*jointi.pelvis.tx-1) - ...
            Xk_Qs_Qdots_opt(mpi).mp(1,2*jointi.pelvis.tx-1); % distance traveled 
        speed_opt = dist_trav_opt/tf_opt;
        assertSpeed = abs(speed_opt - gaitSpeed);
        if assertSpeed > 1e-10
            disp(['Issue when reconstructing target speed: difference is ',...
                num2str(assertSpeed)])
        end         
        for k=1:N               
            % Get muscle-tendon lengths, velocities, moment arms
            % Left leg
            qin_l_opt_all = [...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_flex.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_add.l*2-1), ...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_rot.l*2-1), ...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.knee.l*2-1), ...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.ankle.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.subt.l*2-1)];  
            qdotin_l_opt_all = [...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_flex.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_add.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_rot.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.knee.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.ankle.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.subt.l*2)];  
            [lMTk_l_opt_all,vMTk_l_opt_all,~] = ...
                f_lMT_vMT_dM_l(qin_l_opt_all,qdotin_l_opt_all);    
            % Right leg
            qin_r_opt_all = [...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_flex.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_add.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_rot.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.knee.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.ankle.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.subt.r*2-1)];  
            qdotin_r_opt_all = [...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_flex.r*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_add.r*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_rot.r*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.knee.r*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.ankle.r*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.subt.r*2)];      
            [lMTk_r_opt_all,vMTk_r_opt_all,~] = ...
                f_lMT_vMT_dM_r(qin_r_opt_all,qdotin_r_opt_all);
            % Both legs
            lMTk_lr_opt_all = [lMTk_l_opt_all;lMTk_r_opt_all];
            vMTk_lr_opt_all = [vMTk_l_opt_all;vMTk_r_opt_all];        
            % Metabolic energy rate
            [~,~,Fce_opt_all,Fiso_opt_all,~,massM_opt_all,Fpass_opt_all,~] = ...
                f_forceEquilibrium_FtildeState_all(...
                    a_tot_opt(mpi).mp(k,musi_FLV)',...
                    FTtilde_opt_unsc(mpi).mp(k,:)',...
                    dFTtilde_opt_unsc(mpi).mp(k,:)',...
                    full(lMTk_lr_opt_all(musi_FLV)),...
                    full(vMTk_lr_opt_all(musi_FLV)),tensions(musi_FLV));                  
            [~,lMtilde_opt_all] = f_FiberLength_TendonForce(...
                FTtilde_opt_unsc(mpi).mp(k,:)',full(lMTk_lr_opt_all(musi_FLV)));                
            [vM_opt_all,~] = f_FiberVelocity_TendonForce(...
                FTtilde_opt_unsc(mpi).mp(k,:)',...
                dFTtilde_opt_unsc(mpi).mp(k,:)',...
                full(lMTk_lr_opt_all(musi_FLV)),...
                full(vMTk_lr_opt_all(musi_FLV)));
            [e_tot_all,~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2004all(...
            a_tot_opt(mpi).mp(k,musi_FLV)',a_tot_opt(mpi).mp(k,musi_FLV)',...
            full(lMtilde_opt_all),full(vM_opt_all),full(Fce_opt_all)',...
            full(Fpass_opt_all)',full(massM_opt_all)',pctsts(musi_FLV),...
            full(Fiso_opt_all)',MTParameters_FLV(1,:)',body_mass,10);                 
            e_tot_opt_all = full(e_tot_all)';
            e_mo_opt(mpi).mp(k,:) = full(e_mot)'; 
            for j=1:d 
                % Tracking term
                Qs_cost_optk = W.TrackTD*B(j+1)...
                    *(f_JNq_bpty(qqdot_opt(mpi).mp(k,Qsi(jointi.res_bptyi))- ...
                    Qs(mpi).mp.allinterpfilt(k,jointi.res_bptyi+1)... 
                    ./scaling(mpi).mp.Qs(jointi.res_bptyi)))*h_opt;
                Qs_cost = Qs_cost + Qs_cost_optk;                
                trackTerms_optk = Qs_cost_optk;
                % Motor control terms   
                % Muscle activations
                a_cost_optk = W.Act*B(j+1)...
                    *(f_JNM_exp(a_opt_unsc(mpi).mp(k,:),exp_A))*h_opt;
                a_cost = a_cost + a_cost_optk;
                % Metabolic energy rate
                if W.MER == 0
                    mE_cost_optk = 0;
                else
                    mE_cost_optk = W.MER*B(j+1)...
                        *(f_JNM_FLVexp(e_tot_opt_all,exp_E))/body_mass*h_opt;
                end
                mE_cost = mE_cost + mE_cost_optk;
                % Trunk excitations
                trunk_cost_optk = W.TrunkExc*B(j+1)...
                    *(f_Jnq_trunk(e_b_opt_unsc(mpi).mp(k,:)))*h_opt; 
                trunk_cost = trunk_cost + trunk_cost_optk;
                % Joint accelerations
                Qdotdots_cost_optk = W.Acc*B(j+1)...
                    *(f_JNq_all(qdotdot_opt(mpi).mp(k,:)))*h_opt;
                Qdotdots_cost = Qdotdots_cost + Qdotdots_cost_optk;
                % Time derivative of muscle activations
                vA_cost_optk = W.u*B(j+1)*(f_JNMact(vA_opt(mpi).mp(k,:)))*h_opt;
                vA_cost = vA_cost + vA_cost_optk;
                % Time derivative of muscle-tendon forces
                dFTtilde_cost_optk = W.u*B(j+1)...
                    *(f_JNM_FLV(dFTtilde_opt(mpi).mp(k,:)))*h_opt;
                dFTtilde_cost = dFTtilde_cost + dFTtilde_cost_optk;   
                % Passive joint torques
                passT_cost_optk = W.PassT*B(j+1)...
                    *(f_JNq_act(TPass_opt_all(mpi).mp(k,:)))*h_opt;
                passT_cost = passT_cost + passT_cost_optk;               
                mcTerms_optk = a_cost_optk + mE_cost_optk +...
                    Qdotdots_cost_optk + passT_cost_optk + trunk_cost_optk + ...
                    vA_cost_optk + dFTtilde_cost_optk; 
                dist_norm_opt = dist_trav_opt;
                J_opt = J_opt + trackTerms_optk + 1/dist_norm_opt*mcTerms_optk;  
                count = count + 1;                 
            end
        end            
        e_mo_opt_tr(mpi).mp = trapz(tgrid(mpi).mp(1:end-1),e_mo_opt(mpi).mp);
        % Energy model from Bhargava et al. (2004)
        COT_opt(mpi).mp = e_mo_opt_tr(mpi).mp/body_mass/dist_trav_opt; 
    end   
    % assertCost should be ~ 0     
    assertCost = abs(full(J_opt) - stats.iterations.obj(end));
    if assertCost > 1e-8
        disp(['Issue when reconstructing optimal cost: difference is ',...
            num2str(assertCost)])
    end

    %% Extract fiber lengths and passive forces, and assert Hill-equilibrium   
    Hilldiffk_opt = struct('mp',[]);
    Fpe_opt = struct('mp',[]);
    lMtilde_opt = struct('mp',[]);
    lMTk_opt_lr = struct('mp',[]);
    vMTk_opt_lr = struct('mp',[]);
    for mpi = 1:NPhases
        for i = 1:N   
            % Get muscle-tendon lengths, velocities, and moment arms
            % Left leg
            qin_opt_l = [Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.l*2-1), ...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.l*2-1), ...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.l*2-1), ...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.l*2-1)];  
            qdotin_opt_l = [Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.l*2)];  
            [lMTk_opt_l,vMTk_opt_l,~] = f_lMT_vMT_dM_l(qin_opt_l,qdotin_opt_l); 
            % Right leg
            qin_opt_r = [Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.r*2-1)];  
            qdotin_opt_r = [Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.r*2),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.r*2),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.r*2),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.r*2),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.r*2),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.r*2)];      
            [lMTk_opt_r,vMTk_opt_r,~] = f_lMT_vMT_dM_r(qin_opt_r,qdotin_opt_r);
            % Both legs
            lMTk_opt_lr(mpi).mp(i,:) = [full(lMTk_opt_l);full(lMTk_opt_r)];
            vMTk_opt_lr(mpi).mp(i,:) = [full(vMTk_opt_l);full(vMTk_opt_r)];   
            % Get muscle-tendon forces and derive Hill-equilibrium       
            [Hilldiffk_t,~,~,~,~,~,Fpe_t,lMtilde_t] = ...
                f_forceEquilibrium_FtildeState_all(...
                a_tot_opt(mpi).mp(i,musi_FLV),FTtilde_opt_unsc(mpi).mp(i,:),...
                dFTtilde_opt_unsc(mpi).mp(i,:),...
                lMTk_opt_lr(mpi).mp(i,musi_FLV),...
                vMTk_opt_lr(mpi).mp(i,musi_FLV),tensions(musi_FLV)');
           Hilldiffk_opt(mpi).mp(i,:) = full(Hilldiffk_t');
           Fpe_opt(mpi).mp(i,:) = full(Fpe_t');     
           lMtilde_opt(mpi).mp(i,:) = full(lMtilde_t');      
        end
        % assertHill should be ~ 0  
        assertHill = max(max(Hilldiffk_opt(mpi).mp));
        if assertHill > 1*10^(-tol_ipopt)
            disp('Issue in Hill-equilibrium')
        end
    end
    
    %% Reconstruct gait cycle: starting with right heel strike
    % Identify heel strike
    threshold = 30;
    if ww == 34
        threshold = 28;
    end
    IC1i = struct('mp',[]);
    Qs_opt_GC_r = struct('mp',[]);
    a_opt_unsc_GC_r = struct('mp',[]);
    a_Ff_opt_unsc_GC_r = struct('mp',[]);
    a_dFf_opt_unsc_GC_r = struct('mp',[]);
    a_tot_opt_GC_r = struct('mp',[]);
    GRF_opt_GC_r = struct('mp',[]);
    GRM_opt_GC_r = struct('mp',[]);
    Tauk_out_GC_r = struct('mp',[]);
    Fpe_opt_GC_r = struct('mp',[]);
    lMtilde_opt_GC_r = struct('mp',[]);
    lMTk_opt_lr_GC_r = struct('mp',[]);
    vMTk_opt_lr_GC_r = struct('mp',[]);
    a_syn_opt_GC_r = struct('mp',[]);
    for mpi = 1:NPhases
        if exist('HS1','var')
            clear HS1
        end
        % Right heel strike first    
        phase_tran_tgridi = find(GRF_opt_unsc(mpi).mp(:,2)<threshold,1,'last');
        if ~isempty(phase_tran_tgridi)        
            if phase_tran_tgridi == N
                temp_idx = find(GRF_opt_unsc(mpi).mp(:,2)>threshold,1,'first');
                if ~isempty(temp_idx)
                    if temp_idx-1 ~= 0 && ...
                            find(GRF_opt_unsc(mpi).mp(temp_idx-1,2)<threshold)
                        phase_tran_tgridi_t = temp_idx;             
                        IC1i(mpi).mp = phase_tran_tgridi_t;
                        HS1 = 'r';
                    end 
                else            
                    IC1i(mpi).mp = phase_tran_tgridi + 1; 
                    HS1 = 'r';
                end
            else            
                IC1i(mpi).mp = phase_tran_tgridi + 1; 
                HS1 = 'r';
            end        
        end
        if isempty(phase_tran_tgridi)
            continue;
        end   
        % Joint angles
        Qs_opt_GC_r(mpi).mp = zeros(N,size(q_opt_unsc(mpi).mp.deg,2));    
        Qs_opt_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
            q_opt_unsc(mpi).mp.deg(IC1i(mpi).mp:end,:);
        Qs_opt_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
            q_opt_unsc(mpi).mp.deg(1:IC1i(mpi).mp-1,:);           
        Qs_opt_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,jointi.pelvis.tx) = ...
            q_opt_unsc(mpi).mp.deg(1:IC1i(mpi).mp-1,jointi.pelvis.tx) + ...
            (q_opt_unsc_all(mpi).mp.deg(end,jointi.pelvis.tx)-...
            q_opt_unsc_all(mpi).mp.deg(1,jointi.pelvis.tx));    
        temp_q_opt_GC_pelvis_tx = Qs_opt_GC_r(mpi).mp(1,jointi.pelvis.tx);
        Qs_opt_GC_r(mpi).mp(:,jointi.pelvis.tx) = ...
            Qs_opt_GC_r(mpi).mp(:,jointi.pelvis.tx) - temp_q_opt_GC_pelvis_tx;           
        % Muscle activations (excluding feedback component)
        a_opt_unsc_GC_r(mpi).mp = zeros(N,NMuscles);    
        a_opt_unsc_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
            a_opt_unsc(mpi).mp(IC1i(mpi).mp:end,:);    
        a_opt_unsc_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
            a_opt_unsc(mpi).mp(1:IC1i(mpi).mp-1,:);            
        % Muscle activations (feedback component)
        if spasi == 1
            % Force feedback
            a_Ff_opt_unsc_GC_r(mpi).mp = zeros(N,NMuscles_Spas);    
            a_Ff_opt_unsc_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
                a_Ff_opt_unsc(mpi).mp(IC1i(mpi).mp:end,:);    
            a_Ff_opt_unsc_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
                a_Ff_opt_unsc(mpi).mp(1:IC1i(mpi).mp-1,:);
            % Force rate feedback
            a_dFf_opt_unsc_GC_r(mpi).mp = zeros(N,NMuscles_Spas);    
            a_dFf_opt_unsc_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
                a_dFf_opt_unsc(mpi).mp(IC1i(mpi).mp:end,:);    
            a_dFf_opt_unsc_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
                a_dFf_opt_unsc(mpi).mp(1:IC1i(mpi).mp-1,:);
        end        
        % Muscle activations (including feedback component)
        a_tot_opt_GC_r(mpi).mp = zeros(N,NMuscles);    
        a_tot_opt_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
            a_tot_opt(mpi).mp(IC1i(mpi).mp:end,:);    
        a_tot_opt_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
            a_tot_opt(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Ground reaction forces
        GRF_opt_GC_r(mpi).mp = zeros(N,nGRF);
        GRF_opt_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
            GRF_opt_unsc(mpi).mp(IC1i(mpi).mp:end,:);
        GRF_opt_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
            GRF_opt_unsc(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Ground reaction torques
        GRM_opt_GC_r(mpi).mp = zeros(N,nGRF);
        GRM_opt_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
            GRM_opt_unsc(mpi).mp(IC1i(mpi).mp:end,:);
        GRM_opt_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
            GRM_opt_unsc(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Joint torques
        Tauk_out_GC_r(mpi).mp = zeros(N,nq.res);    
        Tauk_out_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
            Tauk_out(mpi).mp(IC1i(mpi).mp:end,:);    
        Tauk_out_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
            Tauk_out(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Passive forces
        Fpe_opt_GC_r(mpi).mp = zeros(N,NMuscles);    
        Fpe_opt_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
            Fpe_opt(mpi).mp(IC1i(mpi).mp:end,:);    
        Fpe_opt_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
            Fpe_opt(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Fiber lengths
        lMtilde_opt_GC_r(mpi).mp = zeros(N,NMuscles);    
        lMtilde_opt_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
            lMtilde_opt(mpi).mp(IC1i(mpi).mp:end,:);    
        lMtilde_opt_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
            lMtilde_opt(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Muscle-tendon lengths
        lMTk_opt_lr_GC_r(mpi).mp = zeros(N,NMuscles);    
        lMTk_opt_lr_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
            lMTk_opt_lr(mpi).mp(IC1i(mpi).mp:end,:);    
        lMTk_opt_lr_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
            lMTk_opt_lr(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Muscle-tendon velocities
        vMTk_opt_lr_GC_r(mpi).mp = zeros(N,NMuscles);    
        vMTk_opt_lr_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
            vMTk_opt_lr(mpi).mp(IC1i(mpi).mp:end,:);    
        vMTk_opt_lr_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
            vMTk_opt_lr(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Synergy activations
        if NSyn ~= 99
            a_syn_opt_GC_r(mpi).mp = zeros(N,NMact);    
            a_syn_opt_GC_r(mpi).mp(1:N-IC1i(mpi).mp+1,:) = ...
                a_syn_opt(mpi).mp(IC1i(mpi).mp:end,:);    
            a_syn_opt_GC_r(mpi).mp(N-IC1i(mpi).mp+2:N,:) = ...
                a_syn_opt(mpi).mp(1:IC1i(mpi).mp-1,:);
        end      
        % Visualization in OpenSim GUI    
        if writeMotion_r    
            % Joint angles
            q_opt_GUI_r = zeros(N,1+nq.res);
            q_opt_GUI_r(:,1) = tgrid(mpi).mp(1:end-1)';
            q_opt_GUI_r(:,2:nq.res+1)  = Qs_opt_GC_r(mpi).mp;
            % Muscle activations (to have muscles turning red when activated)
            Acts_opt_GUI_r = a_tot_opt_GC_r(mpi).mp;
            % Combine data joint angles and muscle activations
            JointAngleMuscleAct_r.data = [q_opt_GUI_r,Acts_opt_GUI_r];
            % Get muscle labels
            muscleNamesAll = cell(1,NMuscles);
            for i = 1:NMuscles/2
                muscleNamesAll{i} = [muscleNames{i}(1:end-2),'_l'];
                muscleNamesAll{i+NMuscles/2} = [muscleNames{i}(1:end-2),'_r'];
            end  
            JointAngles_r.labels = {'time','lower_torso_RX','lower_torso_RY',...
                'lower_torso_RZ','lower_torso_TX','lower_torso_TY',...
                'lower_torso_TZ','hip_flex_l','hip_add_l','hip_rot_l',...
                'hip_flex_r','hip_add_r','hip_rot_r','knee_flex_l',...
                'knee_flex_r','ankle_flex_l','ankle_flex_r','subt_angle_l',...
                'subt_angle_r','lumbar_pitch','lumbar_roll','lumbar_yaw'};
            % Combine labels joint angles and muscle activations
            JointAngleMuscleAct_r.labels = JointAngles_r.labels;
            for i = 1:NMuscles
                JointAngleMuscleAct_r.labels{i+size(q_opt_GUI_r,2)} = ...
                    [muscleNamesAll{i},'/activation'];
            end
            filenameJointAngleMuscleAct_r = [pathPredictiveSimulations,...
                '/Results/',namescript,'/IK',savenamePhase(mpi).mp,'_r.mot'];
            write_motionFile(JointAngleMuscleAct_r,...
                filenameJointAngleMuscleAct_r);
        end        
    end
    
    %% Reconstruct gait cycle: starting with left heel strike
    % Identify heel strike
    threshold = 30;
    IC1i_l = struct('mp',[]);
    Qs_opt_GC_l = struct('mp',[]);
    a_opt_unsc_GC_l = struct('mp',[]);
    a_Ff_opt_unsc_GC_l = struct('mp',[]);
    a_dFf_opt_unsc_GC_l = struct('mp',[]);
    a_tot_opt_GC_l = struct('mp',[]);
    GRF_opt_GC_l = struct('mp',[]);
    GRM_opt_GC_l = struct('mp',[]);
    Tauk_out_GC_l = struct('mp',[]);
    Fpe_opt_GC_l = struct('mp',[]);
    lMtilde_opt_GC_l = struct('mp',[]);
    lMTk_opt_lr_GC_l = struct('mp',[]);
    vMTk_opt_lr_GC_l = struct('mp',[]);
    a_syn_opt_GC_l = struct('mp',[]);
    for mpi = 1:NPhases
        if exist('HS1_l','var')
            clear HS1_l
        end
        % Left heel strike first    
        phase_tran_tgridi_l =find(GRF_opt_unsc(mpi).mp(:,5)<threshold,1,'last');
        if ~isempty(phase_tran_tgridi_l)        
            if phase_tran_tgridi_l == N
                temp_idx = find(GRF_opt_unsc(mpi).mp(:,5)>threshold,1,'first');
                if ~isempty(temp_idx)
                    if temp_idx-1 ~= 0 && ...
                            find(GRF_opt_unsc(mpi).mp(temp_idx-1,5)<threshold)
                        phase_tran_tgridi_t = temp_idx;             
                        IC1i_l(mpi).mp = phase_tran_tgridi_t;
                        HS1_l = 'l';
                    end 
                else            
                    IC1i_l(mpi).mp = phase_tran_tgridi_l + 1; 
                    HS1_l = 'l';
                end
            else            
                IC1i_l(mpi).mp = phase_tran_tgridi_l + 1; 
                HS1_l = 'l';
            end        
        end
        if isempty(phase_tran_tgridi_l)
            continue;
        end   
        % Joint angles
        Qs_opt_GC_l(mpi).mp = zeros(N,size(q_opt_unsc(mpi).mp.deg,2));    
        Qs_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
            q_opt_unsc(mpi).mp.deg(IC1i_l(mpi).mp:end,:);
        Qs_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
            q_opt_unsc(mpi).mp.deg(1:IC1i_l(mpi).mp-1,:);           
        Qs_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,jointi.pelvis.tx) = ...
            q_opt_unsc(mpi).mp.deg(1:IC1i_l(mpi).mp-1,jointi.pelvis.tx) + ...
            (q_opt_unsc_all(mpi).mp.deg(end,jointi.pelvis.tx)-...
            q_opt_unsc_all(mpi).mp.deg(1,jointi.pelvis.tx));    
        temp_q_opt_GC_pelvis_tx = Qs_opt_GC_l(mpi).mp(1,jointi.pelvis.tx);
        Qs_opt_GC_l(mpi).mp(:,jointi.pelvis.tx) = ...
            Qs_opt_GC_l(mpi).mp(:,jointi.pelvis.tx) - temp_q_opt_GC_pelvis_tx;  
        % Muscle activations (excluding feedback component)
        a_opt_unsc_GC_l(mpi).mp = zeros(N,NMuscles);    
        a_opt_unsc_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
            a_opt_unsc(mpi).mp(IC1i_l(mpi).mp:end,:);    
        a_opt_unsc_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
            a_opt_unsc(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Muscle activations (feedback component)
        a_tot_opt_GC_l(mpi).mp = zeros(N,NMuscles);    
        a_tot_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
            a_tot_opt(mpi).mp(IC1i_l(mpi).mp:end,:);    
        a_tot_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
            a_tot_opt(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        if spasi == 1
            % Force feedback
            a_Ff_opt_unsc_GC_l(mpi).mp = zeros(N,NMuscles_Spas);    
            a_Ff_opt_unsc_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
                a_Ff_opt_unsc(mpi).mp(IC1i_l(mpi).mp:end,:);    
            a_Ff_opt_unsc_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
                a_Ff_opt_unsc(mpi).mp(1:IC1i_l(mpi).mp-1,:);
            % Force rate feedback
            a_dFf_opt_unsc_GC_l(mpi).mp = zeros(N,NMuscles_Spas);    
            a_dFf_opt_unsc_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
                a_dFf_opt_unsc(mpi).mp(IC1i_l(mpi).mp:end,:);    
            a_dFf_opt_unsc_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
                a_dFf_opt_unsc(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        end
        % Ground reaction forces
        GRF_opt_GC_l(mpi).mp = zeros(N,nGRF);
        GRF_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
            GRF_opt_unsc(mpi).mp(IC1i_l(mpi).mp:end,:);
        GRF_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
            GRF_opt_unsc(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Ground reaction torques
        GRM_opt_GC_l(mpi).mp = zeros(N,nGRF);
        GRM_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
            GRM_opt_unsc(mpi).mp(IC1i_l(mpi).mp:end,:);
        GRM_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
            GRM_opt_unsc(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Joint torques
        Tauk_out_GC_l(mpi).mp = zeros(N,nq.res);    
        Tauk_out_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
            Tauk_out(mpi).mp(IC1i_l(mpi).mp:end,:);    
        Tauk_out_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
            Tauk_out(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Passive forces
        Fpe_opt_GC_l(mpi).mp = zeros(N,NMuscles);    
        Fpe_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
            Fpe_opt(mpi).mp(IC1i_l(mpi).mp:end,:);    
        Fpe_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
            Fpe_opt(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Fiber lengths
        lMtilde_opt_GC_l(mpi).mp = zeros(N,NMuscles);    
        lMtilde_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
            lMtilde_opt(mpi).mp(IC1i_l(mpi).mp:end,:);    
        lMtilde_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
            lMtilde_opt(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Muscle-tendon lengths
        lMTk_opt_lr_GC_l(mpi).mp = zeros(N,NMuscles);    
        lMTk_opt_lr_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
            lMTk_opt_lr(mpi).mp(IC1i_l(mpi).mp:end,:);    
        lMTk_opt_lr_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
            lMTk_opt_lr(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Muscle-tendon velocities
        vMTk_opt_lr_GC_l(mpi).mp = zeros(N,NMuscles);    
        vMTk_opt_lr_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
            vMTk_opt_lr(mpi).mp(IC1i_l(mpi).mp:end,:);    
        vMTk_opt_lr_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
            vMTk_opt_lr(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Synergy activations
        if NSyn ~= 99
            a_syn_opt_GC_l(mpi).mp = zeros(N,NMact);    
            a_syn_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = ...
                a_syn_opt(mpi).mp(IC1i_l(mpi).mp:end,:);    
            a_syn_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = ...
                a_syn_opt(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        end
        % Visualization in OpenSim GUI    
        if writeMotion_l 
            % Joint angles
            q_opt_GUI_l = zeros(N,1+nq.res);
            q_opt_GUI_l(:,1) = tgrid(mpi).mp(1:end-1)';
            q_opt_GUI_l(:,2:nq.res+1)  = Qs_opt_GC_l(mpi).mp;
            % Muscle activations (to have muscles turning red when activated)
            Acts_opt_GUI_l = a_tot_opt_GC_l(mpi).mp;
            % Combine data joint angles and muscle activations
            JointAngleMuscleAct_l.data = [q_opt_GUI_l,Acts_opt_GUI_l];
            % Get muscle labels
            muscleNamesAll = cell(1,NMuscles);
            for i = 1:NMuscles/2
                muscleNamesAll{i} = [muscleNames{i}(1:end-2),'_l'];
                muscleNamesAll{i+NMuscles/2} = [muscleNames{i}(1:end-2),'_r'];
            end  
            JointAngles_l.labels = {'time','lower_torso_RX','lower_torso_RY',...
                'lower_torso_RZ','lower_torso_TX','lower_torso_TY',...
                'lower_torso_TZ','hip_flex_l','hip_add_l','hip_rot_l',...
                'hip_flex_r','hip_add_r','hip_rot_r','knee_flex_l',...
                'knee_flex_r','ankle_flex_l','ankle_flex_r','subt_angle_l',...
                'subt_angle_r','lumbar_pitch','lumbar_roll','lumbar_yaw'};
            % Combine labels joint angles and muscle activations
            JointAngleMuscleAct_l.labels = JointAngles_l.labels;
            for i = 1:NMuscles
                JointAngleMuscleAct_l.labels{i+size(q_opt_GUI_l,2)} = ...
                    [muscleNamesAll{i},'/activation'];
            end
            filenameJointAngleMuscleAct_l = [pathPredictiveSimulations,...
                '/Results/',namescript,'/IK',savenamePhase(mpi).mp,'_l.mot'];
            write_motionFile(JointAngleMuscleAct_l,...
                filenameJointAngleMuscleAct_l);
        end
    end    
        
    %% Save results
    for mpi = 1:NPhases
        if saveResults
            if (exist([pathresults,'/',namescript,'/ResultsPredSim.mat'],...
                'file')==2) 
            load([pathresults,'/',namescript,'/ResultsPredSim.mat']);
            else
                ResultsPredSim.(['Case_',num2str(ww)]) = struct('Qs_opt',[]);
            end
            % Structure results
            ResultsPredSim.(['Case_',num2str(ww)]).tgrid(mpi).mp = ...
                tgrid(mpi).mp(1:end-1)';
            ResultsPredSim.(['Case_',num2str(ww)]).Qs_r(mpi).mp = ...
                Qs_opt_GC_r(mpi).mp;
            ResultsPredSim.(['Case_',num2str(ww)]).Qs_l(mpi).mp = ...
                Qs_opt_GC_l(mpi).mp;
            ResultsPredSim.(['Case_',num2str(ww)]).Acts_noSpas_r(mpi).mp = ...
                a_opt_unsc_GC_r(mpi).mp;  
            ResultsPredSim.(['Case_',num2str(ww)]).Acts_noSpas_l(mpi).mp = ...
                a_opt_unsc_GC_l(mpi).mp;  
            ResultsPredSim.(['Case_',num2str(ww)]).Acts_r(mpi).mp = ...
                a_tot_opt_GC_r(mpi).mp;   
            ResultsPredSim.(['Case_',num2str(ww)]).Acts_l(mpi).mp = ...
                a_tot_opt_GC_l(mpi).mp;   
            if spasi == 1
                ResultsPredSim.(['Case_',num2str(ww)]).Acts_Ff_r(mpi).mp = ...
                    a_Ff_opt_unsc_GC_r(mpi).mp;
                ResultsPredSim.(['Case_',num2str(ww)]).Acts_Ff_l(mpi).mp = ...
                    a_Ff_opt_unsc_GC_l(mpi).mp;
                ResultsPredSim.(['Case_',num2str(ww)]).Acts_dFf_r(mpi).mp = ...
                    a_dFf_opt_unsc_GC_r(mpi).mp;
                ResultsPredSim.(['Case_',num2str(ww)]).Acts_dFf_l(mpi).mp = ...
                    a_dFf_opt_unsc_GC_l(mpi).mp;
            else
                ResultsPredSim.(['Case_',num2str(ww)]).Acts_Ff_r(mpi).mp = NaN;
                ResultsPredSim.(['Case_',num2str(ww)]).Acts_Ff_l(mpi).mp = NaN;
                ResultsPredSim.(['Case_',num2str(ww)]).Acts_dFf_r(mpi).mp = NaN;
                ResultsPredSim.(['Case_',num2str(ww)]).Acts_dFf_l(mpi).mp = NaN;
            end
            ResultsPredSim.(['Case_',num2str(ww)]).Ts_r(mpi).mp = ...
                Tauk_out_GC_r(mpi).mp;     
            ResultsPredSim.(['Case_',num2str(ww)]).Ts_l(mpi).mp = ...
                Tauk_out_GC_l(mpi).mp; 
            ResultsPredSim.(['Case_',num2str(ww)]).GRFs_r(mpi).mp = ...
                GRF_opt_GC_r(mpi).mp;        
            ResultsPredSim.(['Case_',num2str(ww)]).GRFs_l(mpi).mp = ...
                GRF_opt_GC_l(mpi).mp;     
            ResultsPredSim.(['Case_',num2str(ww)]).GRMs_r(mpi).mp = ...
                GRM_opt_GC_r(mpi).mp;  
            ResultsPredSim.(['Case_',num2str(ww)]).GRMs_l(mpi).mp = ...
                GRM_opt_GC_l(mpi).mp; 
            ResultsPredSim.(['Case_',num2str(ww)]).Fpe_r(mpi).mp = ...
                Fpe_opt_GC_r(mpi).mp;
            ResultsPredSim.(['Case_',num2str(ww)]).Fpe_l(mpi).mp = ...
                Fpe_opt_GC_l(mpi).mp;
            ResultsPredSim.(['Case_',num2str(ww)]).lMtilde_r(mpi).mp = ...
                lMtilde_opt_GC_r(mpi).mp;    
            ResultsPredSim.(['Case_',num2str(ww)]).lMtilde_l(mpi).mp = ...
                lMtilde_opt_GC_l(mpi).mp; 
            ResultsPredSim.(['Case_',num2str(ww)]).lMTk_r(mpi).mp = ...
                lMTk_opt_lr_GC_r(mpi).mp;         
            ResultsPredSim.(['Case_',num2str(ww)]).lMTk_l(mpi).mp = ...
                lMTk_opt_lr_GC_l(mpi).mp; 
            ResultsPredSim.(['Case_',num2str(ww)]).vMTk_r(mpi).mp = ...
                vMTk_opt_lr_GC_r(mpi).mp;
            ResultsPredSim.(['Case_',num2str(ww)]).vMTk_l(mpi).mp = ...
                vMTk_opt_lr_GC_l(mpi).mp;
            if NSyn ~= 99
                ResultsPredSim.(['Case_',num2str(ww)]).Acts_syn_r(mpi).mp = ...
                    a_syn_opt_GC_r(mpi).mp; 
                ResultsPredSim.(['Case_',num2str(ww)]).Acts_syn_l(mpi).mp = ...
                    a_syn_opt_GC_l(mpi).mp; 
                ResultsPredSim.(['Case_',num2str(ww)]).w_syn(mpi).mp = ...
                    [syn_wl_opt,syn_wr_opt];     
            end
            ResultsPredSim.(['Case_',num2str(ww)]).COT(mpi).mp= COT_opt(mpi).mp;
            ResultsPredSim.(['Case_',num2str(ww)]).Qs_toTrack(mpi).mp = ...
                Qs(mpi).mp.allinterpfilt;
            ResultsPredSim.(['Case_',num2str(ww)]).Ts_toTrack(mpi).mp = ...
                ID(mpi).mp.allinterp;
            ResultsPredSim.(['Case_',num2str(ww)]).GRFs_toTrack(mpi).mp = ...
                GRF(mpi).mp.val.allinterp;
            ResultsPredSim.(['Case_',num2str(ww)]).GRMs_toTrack(mpi).mp = ...
                GRF(mpi).mp.MorGF.allinterp; 
            ResultsPredSim.(['Case_',num2str(ww)]).stats = stats;
            ResultsPredSim.(['Case_',num2str(ww)]).StrideLength = ...
                StrideLength_opt;
            ResultsPredSim.(['Case_',num2str(ww)]).StepWidthMean = ...
                StepWidth_opt_mean;
            ResultsPredSim.(['Case_',num2str(ww)]).StepWidthStd = ...
                StepWidth_opt_std;
            ResultsPredSim.colheaders.joints = joints;
            ResultsPredSim.colheaders.GRF = {'fore_aft_r','vertical_r',...
                'lateral_r','fore_aft_l','vertical_l','lateral_l'};
            for i = 1:NMuscles/2
                    ResultsPredSim.colheaders.muscles{i} = ...
                        [muscleNames{i}(1:end-2),'_l'];
                    ResultsPredSim.colheaders.muscles{i+NMuscles/2} = ...
                        [muscleNames{i}(1:end-2),'_r'];
            end
            % Save data
            save([pathresults,'/',namescript,'/ResultsPredSim.mat'],...
                'ResultsPredSim');
        end 
    end
end
end
