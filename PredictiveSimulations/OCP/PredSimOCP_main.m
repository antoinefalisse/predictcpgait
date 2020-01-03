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
% num_set(5): set to 1 to visualize guess-bounds 
% num_set(6): set to 1 to write .mot file
% num_set(7): set to 1 to write .mot file (HS on left leg)
% num_set(8): set to 1 to check for contact forces
% num_set(9): set to 1 to decompose cost function

% num_set = [1,1,0,1,0,1,1,1]; % This configuration solves the problem
num_set = [0,1,1,0,0,0,0,0,1]; % This configuration analyzes the results

% The variable settings in the following section will set some parameters 
% of the optimization problem. Through the variable idx_ww, the user can 
% select which row of parameters will be used.
idx_ww = [695]; % Index row in matrix settings

%% Settings
import casadi.*
subject = 'subject1';

solveProblem    = num_set(1); % set to 1 to solve problem
analyseResults  = num_set(2); % set to 1 to analyze results
loadResults     = num_set(3); % set to 1 to load results
saveResults     = num_set(4); % set to 1 to save sens. results
checkBoundsIG   = num_set(5); % set to 1 to visualize guess-bounds 
writeIKmotion   = num_set(6); % set to 1 to write .mot file starting at right heel strike
writeIKmotion_l = num_set(7); % set to 1 to write .mot file starting at left heel strike
testContact     = num_set(8); % set to 1 to check for contact forces
decomposeCost   = num_set(9); % set to 1 to check for contact forces 

pathmain = pwd;
[pathRepo,~,~] = fileparts(pathmain);
pathSettings = [pathRepo,'\Settings'];
addpath(genpath(pathSettings));
PredSimOCP_settings

%% Select settings
for www = 1:length(idx_ww)
ww = idx_ww(www);
% Variable parameters
W.Qs        = settings(ww,1);  % weight joint kinematics
W.a         = settings(ww,2);  % weight muscle activations
N           = settings(ww,3);  % number of mesh intervals
tol_ipopt   = settings(ww,4);  % NLP error tolerance: 1*10^(-settings(9))
NSyn        = settings(ww,5); % number of synergies
W.Syn       = settings(ww,6); % weight synergies
dev.pr      = settings(ww,7); % allowed deviation from pelvis residuals
paramsi     = settings(ww,8); % index MT-parameters
IGi         = settings(ww,9); % index initial guess
W.E         = settings(ww,10); % weight metabolic energy
exp_E       = settings(ww,11); % exponent metabolic energy
exp_A       = settings(ww,12); % exponent muscle activity
spasi       = settings(ww,13); % whether to include spasticity
W.qdotdot   = settings(ww,14);  % weight muscle activations
v_tgt       = settings(ww,15);  % imposed speed
W.Trunk     = settings(ww,16);  % minimize torque excitations
W.passT     = settings(ww,17);  % minimize passive torques
IGc         = settings(ww,18);  % case for ig

% Experimental walking trial to track
trials = {''};
nametrial = struct('mp',[]);
for tr = 1:length(trials)
    nametrial(tr).mp.id = trials{tr};
end
NPhases = length(nametrial);
% Fixed parameter
W.u = 0.001;
% Identifiers for experimental data
nametrial_all = '';
time_opt = struct('mp',[]);
for mpi = 1:NPhases
    nametrial(mpi).mp.ID = ['ID_',nametrial(mpi).mp.id];
    nametrial(mpi).mp.GRF = ['GRF_',nametrial(mpi).mp.id];
    nametrial(mpi).mp.IK = ['KS_',nametrial(mpi).mp.id];
    switch nametrial(mpi).mp.id
        case {'3DGAIT_A_W7_s','3DGAIT_A_W7_s_off'}
            time_opt(mpi).mp = [1.61,2.39];
        case {'3DGAIT_A_W6_s','3DGAIT_A_W6_s_off'}
            time_opt(mpi).mp = [1.88,2.75];
        case {'3DGAIT_A_W9_s','3DGAIT_A_W9_s_off'}
            time_opt(mpi).mp = [2.56,3.43];
        case {'3DGAIT_A_W18_s','3DGAIT_A_W18_s_off'}
            time_opt(mpi).mp = [1.45,2.31];
        case {'3DGAIT_A_W25_s','3DGAIT_A_W25_s_off'}
            time_opt(mpi).mp = [2.27,3.14];
        case {'3DGAIT_A_W27_s','3DGAIT_A_W27_s_off'}
            time_opt(mpi).mp = [1.05,1.97];
        case 'gait_10'
            time_opt(mpi).mp = [2.54,3.05];
        case 'gait_12'
            time_opt(mpi).mp = [2.63,3.13];
    end  
    if mpi == NPhases
        nametrial_all = [nametrial_all,nametrial(mpi).mp.id];
    else
        nametrial_all = [nametrial_all,nametrial(mpi).mp.id,'_'];
    end
end
% The filename used to save the results depends on the settings 
savename = ['_c',num2str(ww)];
savename_ind = struct('mp',[]); 
phasesid = {'a','b','c','d'};
for mpi = 1:NPhases
    savename_ind(mpi).mp = ['_c',num2str(ww),phasesid{mpi}];
end

%% Load external functions
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB.
% We use different external functions. A first external function extracts 
% several parameters of the bodies to which the contact spheres are attached.
% The contact forces are then computed in MATLAB and are inputs of the
% second external function in which the skeleton dynamics is described. The
% motivation for this decoupling is to limit the number of times we need to
% build the model. By defining the contact model in MATLAB, we only need to
% build the model once per external function whereas keeping the contact
% model in the external function would require re-building the model during
% the optimization.
pathExternalFunctions = [pathRepo,'\ExternalFunctions'];
% Loading external functions. 
setup.derivatives =  'AD'; % Algorithmic differentiation
switch setup.derivatives
    case 'AD'     
        cd(pathExternalFunctions);
        F = external('F',['TrackSim_',subject,'_SSCM1_pos.dll']);
end
cd(pathmain);
% This is an example of how to call an external function with some
% numerical values.
% vec2 = zeros(63,1);
% vec2(9) = 0.934;
% res2 = full(F(vec2));

%% Indices external function
% External function: F
% First, joint torques. 
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
% Vectors of indices for later use
residualsi          = jointi.pelvis.list:jointi.trunk.rot; % all 
residuals_acti      = jointi.hip_flex.l:jointi.trunk.rot; % all but gr-pelvis
residual_bptyi      = [jointi.pelvis.list:jointi.pelvis.tx,...
    jointi.pelvis.tz:jointi.trunk.rot]; % all but pelvis_ty
backi              = jointi.trunk.ext:jointi.trunk.rot; % trunk
nq.trunk            = length(backi); % trunk
jointi.radi         = [jointi.pelvis.list:jointi.pelvis.tilt,...
                      jointi.hip_flex.l:jointi.trunk.rot];  
ground_pelvisi      = jointi.pelvis.list:jointi.pelvis.tz; % ground-pelvis
jointi.degi         = jointi.pelvis.tx:jointi.pelvis.tz;
% Number of degrees of freedom for later use
nq.all              = length(residualsi); % all 
nq.abs              = length(ground_pelvisi); % ground-pelvis
nq.act              = nq.all-nq.abs;% all but ground-pelvis
nq.leg              = 6; % #joints needed for polynomials
Qsi                 = 1:2:2*nq.all; % indices Qs only 
% Second, GRFs
GRFi.r              = 22:24;
GRFi.l              = 25:27;
% Third, GRMs
GRMi.r              = 28:30;
GRMi.l              = 31:33;    
GRFi.all            = [GRFi.r,GRFi.l];
nGRF                = length(GRFi.all);
GRMi.all            = [GRMi.r,GRMi.l];
nGRM                = length(GRMi.all);
% Calcaneus
calcOr.r    = 34:35;
calcOr.l    = 36:37;        
% Tibias
tibiaOr.r   = 38:39;
tibiaOr.l   = 40:41;        
calcOr.all  = [calcOr.r,calcOr.l];
NcalcOr     = length(calcOr.all);
tibiaOr.all = [tibiaOr.r,tibiaOr.l];
NtibiaOr    = length(tibiaOr.all);

%% Model info
body_mass = 3.69727 + 4.64393 + 1.43323 + 0.01986 + 0.39058 + 0.04303 + ...
    4.64393 + 1.43323 + 0.01986 + 0.39058 + 0.04303 + 16.25541;
body_weight = body_mass*9.81;

%% Collocation scheme
% We use a pseudospectral direct collocation method, i.e. we use Lagrange
% polynomials to approximate the state derivatives at the collocation
% points in each mesh interval. We use d=3 collocation points per mesh
% interval and Radau collocation points. 
pathCollocationScheme = [pathRepo,'\CollocationScheme'];
addpath(genpath(pathCollocationScheme));
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[tau_root,C,D,B] = CollocationScheme(d,method);

%% Muscle-tendon parameters 
% Muscles from one leg and from the back
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
% Muscle indices for later use
pathmusclemodel = [pathRepo,'\MuscleModel'];
addpath(genpath(pathmusclemodel));    
musi = MuscleIndices_MRI(muscleNames);
% Total number of muscles
NMuscle = length(muscleNames)*2;
% Muscle-tendon parameters. Row 1: maximal isometric forces; Row 2: optimal
% fiber lengths; Row 3: tendon slack lengths; Row 4: optimal pennation 
% angles; Row 5: maximal contraction velocities
if paramsi == 1 % Tuned-parameters
    load([pathmusclemodel,'\MTparameters_',subject,'_L.mat']);
    MTparameters_m = eval(['MTparameters_',subject,'_L']);
elseif paramsi == 2 % Linearly-scaled parameters (as in saved in model)
    load([pathmusclemodel,'\MTparameters_',subject,'_ls.mat']);
    MTparameters_m = eval(['MTparameters_',subject,'_ls']);
elseif paramsi == 3 % Linearly-scaled parameters (new)
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_ls.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 4 % Optimized parameters (new): out 9
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_out_9_out_9_ublM14_ls_simm.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 5 % Optimized parameters (new): out 18
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_out_18_out_18_ublM14_ls_simm.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 6 % Optimized parameters (new): out 18
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_out_18_out_18_ublM15_ls_simm.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 7 % Optimized parameters (new): out 18
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_out_18_out_18_ublM16_ls_simm.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 8 % Optimized parameters (new): out 18
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_out_18_out_18_ublM16_ls_osim.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 9 % Optimized parameters (new): out 18
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_out_18_out_18_ublM15_ls_osim.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 10 % Optimized parameters (new): out 18
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_out_18_out_18_ublM14_ls_osim.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 11 % Optimized parameters (new): out 18
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_ls_osim_opt_L.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 12 % Optimized parameters (new): out 18
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_out_18_out_18_ublM16_ls_L.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 13 % Optimized parameters (new): out 18
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_out_18_out_18_ublM15_ls_old.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 14 % Optimized parameters (new): out 18
    lMtilde_ext.max = 1.4;
    ls = 'simm_gen';
    meth = 'cas';
    optMuscles = 'out_18';
    lMt_ub = num2str(10*lMtilde_ext.max);
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',optMuscles,...
        '_ublM',lMt_ub,'_ls_',ls,'_meth_',meth,'.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 15 % Optimized parameters (new): out 18 / lMopt_max=1.5 / Simm_gen / CasADi
    lMtilde_ext.max = 1.5;
    ls = 'simm_gen';
    meth = 'cas';
    optMuscles = 'out_18';
    lMt_ub = num2str(10*lMtilde_ext.max);
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',optMuscles,...
        '_ublM',lMt_ub,'_ls_',ls,'_meth_',meth,'.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 16 % Scaled generic parameters from Simm model  
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_ls_simm_gen.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 17 % Opt parameters: out 18 /lMopt_max=1.6 /Simm_gen /CasADi
    lMtilde_ext.max = 1.6;
    ls = 'simm_gen';
    meth = 'cas';
    optMuscles = 'out_18';
    lMt_ub = num2str(10*lMtilde_ext.max);
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',optMuscles,...
        '_ublM',lMt_ub,'_ls_',ls,'_meth_',meth,'.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 18 % Opt parameters: out 0 /lMopt_max=1.5 /Simm_gen /CasADi
    lMtilde_ext.max = 1.5;
    ls = 'simm_gen';
    meth = 'cas';
    optMuscles = 'out_0';
    lMt_ub = num2str(10*lMtilde_ext.max);
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',optMuscles,...
        '_ublM',lMt_ub,'_ls_',ls,'_meth_',meth,'.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 19 % Opt parameters: out 9/ lMopt_max=1.5/ Simm_gen /CasADi
    lMtilde_ext.max = 1.5;
    ls = 'simm_gen';
    meth = 'cas';
    optMuscles = 'out_9';
    lMt_ub = num2str(10*lMtilde_ext.max);
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',optMuscles,...
        '_ublM',lMt_ub,'_ls_',ls,'_meth_',meth,'.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 23 % Optimized parameters (new): out 0
    lMtilde_ext.max = 1.5;
    ls = 'simm_gen';
    meth = 'opti';
    optMuscles = 'out_0';
    app = 'DO';
    numbTrials = '4';
    cons_k = 'no';
    cost_k = 'yes';
    NISA = '6';
    lMt_ub = num2str(10*lMtilde_ext.max);
    nameParameters = [optMuscles,'_ublM',lMt_ub,'_ls_',ls,'_meth_',meth,'_app_',...
    app,'_Ntr',numbTrials,'_conK_',cons_k,'_costK_',cost_k,'_NISA',NISA];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters; 
elseif paramsi == 24
    lMtilde_ext.max = 1.5;
    ls = 'simm_gen';
    meth = 'opti';
    optMuscles = 'out_10';
    app = 'DO';
    numbTrials = '4';
    cons_k = 'no';
    cost_k = 'no';
    NISA = '6';
    dev_p = 25;
    lMt_ub = num2str(10*lMtilde_ext.max);
    nameParameters = [optMuscles,'_ublM',lMt_ub,'_ls_',ls,'_meth_',meth,'_app_',...
    app,'_Ntr',numbTrials,'_conK_',cons_k,'_costK_',cost_k,'_NISA',NISA,...
    '_devp',num2str(dev_p),'_2sides'];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters; 
elseif paramsi == 25
    lMtilde_ext.max = 1.5;
    optMuscles = 'out_10';
    numbTrials = '4';
    NISA = '6';
    dev_p = 5;
    lMt_ub = num2str(10*lMtilde_ext.max);
    W.act         = 0.00100;
    W.aT        = 0.99000;
    W.vA        = 0.00001;
    W.dF        = 0.00001;
    W.EMG.all   = 0.003500-W.vA-W.dF;
    W.lMopt_max = 0.004500;
    W.pen       = 0.00100;
    nameParameters = [optMuscles,'_ublM',lMt_ub,...
        '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
        '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
        '_EMG',num2str(round(W.EMG.all*10000)),...
        '_lM',num2str(round(W.lMopt_max*10000)),...
        '_pen',num2str(round(W.pen*10000)),'_2sides'];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters; 
elseif paramsi == 26
    lMtilde_ext.max = 1.5;
    optMuscles = 'out_10';
    numbTrials = '4';
    NISA = '6';
    dev_p = 5;
    lMt_ub = num2str(10*lMtilde_ext.max);
    W.act         = 0.00100;
    W.aT        = 0.98500;
    W.vA        = 0.00001;
    W.dF        = 0.00001;
    W.EMG.all   = 0.003500-W.vA-W.dF;
    W.lMopt_max = 0.009500;
    W.pen       = 0.00100;
    nameParameters = [optMuscles,'_ublM',lMt_ub,...
        '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
        '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
        '_EMG',num2str(round(W.EMG.all*10000)),...
        '_lM',num2str(round(W.lMopt_max*10000)),...
        '_pen',num2str(round(W.pen*10000)),'_2sides'];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters; 
elseif paramsi == 27
    lMtilde_ext.max = 1.5;
    optMuscles = 'out_10';
    numbTrials = '4';
    NISA = '6';
    dev_p = 5;
    lMt_ub = num2str(10*lMtilde_ext.max);
    W.act         = 0.00100;
    W.aT        = 0.98000;
    W.vA        = 0.00001;
    W.dF        = 0.00001;
    W.EMG.all   = 0.003500-W.vA-W.dF;
    W.lMopt_max = 0.014500;
    W.pen       = 0.00100;
    nameParameters = [optMuscles,'_ublM',lMt_ub,...
        '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
        '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
        '_EMG',num2str(round(W.EMG.all*10000)),...
        '_lM',num2str(round(W.lMopt_max*10000)),...
        '_pen',num2str(round(W.pen*10000)),'_2sides'];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters; 
elseif paramsi == 28
    lMtilde_ext.max = 1.5;
    optMuscles = 'out_10';
    numbTrials = '4';
    NISA = '6';
    dev_p = 5;
    lMt_ub = num2str(10*lMtilde_ext.max);
    W.act         = 0.00100;
    W.aT        = 0.90000;
    W.vA        = 0.00001;
    W.dF        = 0.00001;
    W.EMG.all   = 0.003500-W.vA-W.dF;
    W.lMopt_max = 0.094500;
    W.pen       = 0.00100;
    nameParameters = [optMuscles,'_ublM',lMt_ub,...
        '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
        '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
        '_EMG',num2str(round(W.EMG.all*10000)),...
        '_lM',num2str(round(W.lMopt_max*10000)),...
        '_pen',num2str(round(W.pen*10000)),'_2sides'];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters; 
elseif paramsi == 29
    lMtilde_ext.max = 1.5;
    optMuscles = 'out_0';
    numbTrials = '4';
    NISA = '6';
    dev_p = 5;
    lMt_ub = num2str(10*lMtilde_ext.max);
    W.act         = 0.00100;
    W.aT        = 0.85000;
    W.vA        = 0.00001;
    W.dF        = 0.00001;
    W.EMG.all   = 0.003500-W.vA-W.dF;
    W.lMopt_max = 0.144500;
    W.pen       = 0.00100;
    nameParameters = [optMuscles,'_ublM',lMt_ub,...
        '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
        '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
        '_EMG',num2str(round(W.EMG.all*10000)),...
        '_lM',num2str(round(W.lMopt_max*10000)),...
        '_pen',num2str(round(W.pen*10000)),'_2sides'];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters; 
elseif paramsi == 30
    lMtilde_ext.max = 1.5;
    optMuscles = 'out_10';
    numbTrials = '4';
    NISA = '6';
    dev_p = 5;
    lMt_ub = num2str(10*lMtilde_ext.max);
    W.act         = 0.00100;
    W.aT        = 0.85000;
    W.vA        = 0.00001;
    W.dF        = 0.00001;
    W.EMG.all   = 0.003500-W.vA-W.dF;
    W.lMopt_max = 0.144500;
    W.pen       = 0.00100;
    nameParameters = [optMuscles,'_ublM',lMt_ub,...
        '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
        '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
        '_EMG',num2str(round(W.EMG.all*10000)),...
        '_lM',num2str(round(W.lMopt_max*10000)),...
        '_pen',num2str(round(W.pen*10000)),'_2sides'];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters; 
elseif paramsi == 31
    lMtilde_ext.max = 1.5;
    optMuscles = 'out_0';
    numbTrials = '4';
    NISA = '6';
    dev_p = 5;
    lMt_ub = num2str(10*lMtilde_ext.max);
    W.act         = 0.00100;
    W.aT        = 0.60000;
    W.vA        = 0.00001;
    W.dF        = 0.00001;
    W.EMG.all   = 0.003500-W.vA-W.dF;
    W.lMopt_max = 0.394500;
    W.pen       = 0.00100;
    nameParameters = [optMuscles,'_ublM',lMt_ub,...
        '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
        '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
        '_EMG',num2str(round(W.EMG.all*10000)),...
        '_lM',num2str(round(W.lMopt_max*10000)),...
        '_pen',num2str(round(W.pen*10000)),'_2sides'];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters; 
elseif paramsi == 32
    ScaleMIN = 0.01;
    lMtilde_ext.max = 1.5;
    optMuscles = 'out_0';
    numbTrials = '4';
    NISA = '6';
    dev_p = 5;
    lMt_ub = num2str(10*lMtilde_ext.max);
    delta = 0.1;
    W.act       = 0.00100;
    W.aT        = 0.60000;
    W.vA        = 0.00001;
    W.dF        = 0.00001;
    W.EMG.all   = 0.003500-W.vA-W.dF;
    W.lMopt_max = 0.394500;
    W.lMopt     = 0.00100;
    W.act       = 0.00100-W.act*delta/(1-0.00100);
    W.aT        = 0.60000-W.aT*delta/(1-0.00100);
    W.vA        = 0.00001-W.vA*delta/(1-0.00100);
    W.dF        = 0.00001-W.dF*delta/(1-0.00100);
    W.EMG.all   = 0.003480-W.EMG.all*delta/(1-0.00100);
    W.lMopt_max = 0.394500-W.lMopt_max*delta/(1-0.00100);
    W.lMopt     = 0.00100+delta;
    nameParameters = [optMuscles,'_ublM',lMt_ub,...
        '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
        '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
        '_EMG',num2str(round(W.EMG.all*10000)),...
        '_lM',num2str(round(W.lMopt_max*10000)),...
        '_lMopt',num2str(round(W.lMopt*10000)),'_2sides',...
        '_scMin',num2str(ScaleMIN*100)];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters; 
elseif paramsi == 33
    lMtilde_ext.max = 1.8;
    optMuscles = 'out_0';
    numbTrials = '4';
    NISA = '6';
    dev_p = 5;
    lMt_ub = num2str(10*lMtilde_ext.max);
    delta = -0.001;
    W.act       = 0.00100;
    W.aT        = 0.60000;
    W.vA        = 0.00001;
    W.dF        = 0.00001;
    W.EMG.all   = 0.003500-W.vA-W.dF;
    W.lMopt_max = 0.394500;
    W.lMopt     = 0.00100;
    W.act       = 0.00100-W.act*delta/(1-0.00100);
    W.aT        = 0.60000-W.aT*delta/(1-0.00100);
    W.vA        = 0.00001-W.vA*delta/(1-0.00100);
    W.dF        = 0.00001-W.dF*delta/(1-0.00100);
    W.EMG.all   = 0.003480-W.EMG.all*delta/(1-0.00100);
    W.lMopt_max = 0.394500-W.lMopt_max*delta/(1-0.00100);
    W.lMopt     = 0.00100+delta;
    nameParameters = [optMuscles,'_ublM',lMt_ub,...
        '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
        '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
        '_EMG',num2str(round(W.EMG.all*10000)),...
        '_lM',num2str(round(W.lMopt_max*10000)),...
        '_lMopt',num2str(round(W.lMopt*10000)),'_2sides'];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters; 
elseif paramsi == 34
    lMtilde_ext.max = 1.5;
    optMuscles = 'out_0';
    numbTrials = '4';
    NISA = '6';
    dev_p = 5;
    lMt_ub = num2str(10*lMtilde_ext.max);
    ScaleMIN = 0.3;
    delta = 0.1;
    W.act       = 0.00100;
    W.aT        = 0.60000;
    W.vA        = 0.00001;
    W.dF        = 0.00001;
    W.EMG.all   = 0.003500-W.vA-W.dF;
    W.lMopt_max = 0.394500;
    W.lMopt     = 0.00100;
    W.act       = 0.00100-W.act*delta/(1-0.00100);
    W.aT        = 0.60000-W.aT*delta/(1-0.00100);
    W.vA        = 0.00001-W.vA*delta/(1-0.00100);
    W.dF        = 0.00001-W.dF*delta/(1-0.00100);
    W.EMG.all   = 0.003480-W.EMG.all*delta/(1-0.00100);
    W.lMopt_max = 0.394500-W.lMopt_max*delta/(1-0.00100);
    W.lMopt     = 0.00100+delta;
    nameParameters = [optMuscles,'_ublM',lMt_ub,...
        '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
        '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
        '_EMG',num2str(round(W.EMG.all*10000)),...
        '_lM',num2str(round(W.lMopt_max*10000)),...
        '_lMopt',num2str(round(W.lMopt*10000)),'_2sides',...
        '_scMin',num2str(ScaleMIN*100)];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters;
elseif paramsi == 35
    lMtilde_ext.max = 1.5;
    optMuscles = 'out_0';
    numbTrials = '4';
    NISA = '6';
    dev_p = 5;
    lMt_ub = num2str(10*lMtilde_ext.max);
    ScaleMIN = 0.3;
    delta = 0;
    W.act       = 0.00100;
    W.aT        = 0.70000;
    W.vA        = 0.00001;
    W.dF        = 0.00001;
    W.EMG.all   = 0.003500-W.vA-W.dF;
    W.lMopt_max = 0.794500;
    W.lMopt     = 0.00100;
    W.act       = 0.00100-W.act*delta/(1-0.00100);
    W.aT        = 0.70000-W.aT*delta/(1-0.00100);
    W.vA        = 0.00001-W.vA*delta/(1-0.00100);
    W.dF        = 0.00001-W.dF*delta/(1-0.00100);
    W.EMG.all   = 0.003480-W.EMG.all*delta/(1-0.00100);
    W.lMopt_max = 0.294500-W.lMopt_max*delta/(1-0.00100);
    W.lMopt     = 0.00100+delta;
    nameParameters = [optMuscles,'_ublM',lMt_ub,...
        '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
        '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
        '_EMG',num2str(round(W.EMG.all*10000)),...
        '_lM',num2str(round(W.lMopt_max*10000)),...
        '_lMopt',num2str(round(W.lMopt*10000)),'_2sides',...
        '_scMin',num2str(ScaleMIN*100)];
    pathMTParameters = [pathRepo,'\MTParameters'];
    load([pathMTParameters,'\MTParameters_',subject,'_MRI_opt_',...
        nameParameters,'.mat']);
    MTparameters_m = MTParameters;
end
temp_paper(1,:) =  round(MTparameters_m(2,:)*100,2);
temp_paper(2,:) =  round(MTparameters_m(3,:)*100,2);

temp_paper2 = zeros(size(temp_paper));
temp_paper2(:,1:2:end) = temp_paper(:,1:43);
temp_paper2(:,2:2:end) = temp_paper(:,44:end);

% Indices of the muscles actuating the different joints for later use
pathpolynomial = [pathRepo,'\Polynomials'];
addpath(genpath(pathpolynomial));
tl = load([pathpolynomial,'\muscle_spanning_joint_INFO_',subject,'_r.mat']);
[~,mai] = MomentArmIndices_MRI(muscleNames,tl.muscle_spanning_joint_INFO_r);
% Indices of muscles whose FLV relationships are not used
muscleNames_noFLV = {};
NMuscle_noFLV_r  = length(muscleNames_noFLV); % one side
musi_noFLV_r = zeros(1,NMuscle_noFLV_r);
for i = 1:NMuscle_noFLV_r
    musi_noFLV_r(i) = find(strcmp(muscleNames,muscleNames_noFLV{i}));
end
musi_noFLV = [musi_noFLV_r,musi_noFLV_r+NMuscle/2]; % both sides
NMuscle_noFLV = 2*NMuscle_noFLV_r; % both sides
NMuscle_FLV = NMuscle-NMuscle_noFLV;
musi_FLV_r = musi;
musi_FLV_r(musi_noFLV_r) = [];
musi_FLV = [musi,musi+NMuscle/2];
musi_FLV(musi_noFLV) = [];
MTparameters_m_FLV = MTparameters_m(:,musi_FLV);
MTparameters_m_noFLV = MTparameters_m(:,musi_noFLV);

%% Metabolic energy model parameters
% We extract the specific tensions and slow twitch rations.
pathMetabolicEnergy = [pathRepo,'\MetabolicEnergy'];
addpath(genpath(pathMetabolicEnergy));
tension = getSpecificTensions(muscleNames); 
tensions = [tension;tension];
pctst = getSlowTwitchRatios(muscleNames); 
pctsts = [pctst;pctst];

%% Spasticity
pathData = [pathRepo,'\OpenSimModel\',subject];
muscleNames_Spas = {'bi_fem_lh_r','semimem_r','semiten_r',...
    'gas_med_r','gas_lat_r'};
NMuscle_Spas_r  = length(muscleNames_Spas); % one side
musi_Spas_r = zeros(1,NMuscle_Spas_r);
musi_SpasInFLV_r = zeros(1,NMuscle_Spas_r);
for i = 1:NMuscle_Spas_r
    musi_Spas_r(i) = find(strcmp(muscleNames,muscleNames_Spas{i}));
    musi_SpasInFLV_r(i) = ...
        find(strcmp(muscleNames(musi_FLV_r),muscleNames_Spas{i}));    
end
musi_Spas = [musi_Spas_r,musi_Spas_r+NMuscle/2]; % both sides
musi_SpasInFLV = [musi_SpasInFLV_r,musi_SpasInFLV_r+NMuscle_FLV/2]; % both sides
NMuscle_Spas = 2*NMuscle_Spas_r; % both sides
NMuscle_noSpas = NMuscle-NMuscle_Spas;
musi_noSpas = [musi,musi+NMuscle/2];
musi_noSpas(musi_Spas) = [];
if spasi == 1
    bspas = 100;
    
    % load parameters for the spasticity models    
    pathSpasticity = [pathData,'\Spasticity\SpasticGainEstimation\forceModel_',...
        nameParameters,'_',num2str(bspas)];
    formulation_fM = 'forceModel';
    % Hamstrings
    joint = 'knee'; segment_sel = [6 7 8];
    number_save_fM = [formulation_fM '_' int2str(segment_sel)];
    load([pathSpasticity,'\output_',joint,'_',number_save_fM]);  
    % Extract data        
    NMuscles_spas = output.result.setup.auxdata.NMuscles_spas;
    gFf.ham = output.result.solution.parameter(:,1:NMuscles_spas);
    gdFf.ham = output.result.solution.parameter(:,NMuscles_spas+1:2*NMuscles_spas);
    tauFf.ham = output.result.setup.auxdata.Ff_td;
    taudFf.ham = output.result.setup.auxdata.dFf_td;         
    threshold_Ff.ham = output.result.setup.auxdata.threshold_Ff;
    threshold_dFf.ham = output.result.setup.auxdata.threshold_dFf;   
    % The thresholds during gait are the lowest thresholds during passive
    % motions (IPSA)
    threshold_dFf_gait.ham = min(threshold_dFf.ham);
    threshold_Ff_gait.ham = min(threshold_Ff.ham);
    % Gastrocs
    joint = 'ankle'; segment_sel = [19 20 21];
    number_save_fM = [formulation_fM '_' int2str(segment_sel)];
    load([pathSpasticity,'\output_',joint,'_',number_save_fM]);  
    % Extract data        
    NMuscles_spas = output.result.setup.auxdata.NMuscles_spas;
    gFf.gas = output.result.solution.parameter(:,1:NMuscles_spas);
    gdFf.gas = output.result.solution.parameter(:,NMuscles_spas+1:2*NMuscles_spas);
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
    gFf.all = [gFf.side,gFf.side]; % option to increase gains
    gdFf.all = [gdFf.side,gdFf.side]; % option to increase gains
    threshold_dFf_gait.all = [threshold_dFf_gait.side,threshold_dFf_gait.side];
    threshold_Ff_gait.all = [threshold_Ff_gait.side,threshold_Ff_gait.side];
    
    
        
end

%% Synergies
pathSynergies = [pathData,'\Synergies\'];
synData = load([pathSynergies,'\Syn_',subject,'.mat']);
synW.l = synData.Left.W{1,NSyn};
synW.r = synData.Right.W{1,NSyn};

if W.Syn == 0
    NMact = NMuscle;
else
    NMact = 2*NSyn;
end

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'\CasADiFunctions'];
addpath(genpath(pathCasADiFunctions));
pathContactModel = [pathRepo,'\Contact'];
addpath(genpath(pathContactModel));
% We load some variables for the polynomial approximations
load([pathpolynomial,'\muscle_spanning_joint_INFO_',subject,'_r.mat']);
load([pathpolynomial,'\muscle_spanning_joint_INFO_',subject,'_l.mat']);
load([pathpolynomial,'\MuscleInfo_',subject,'_r.mat']);
load([pathpolynomial,'\MuscleInfo_',subject,'_l.mat']);
% For the polynomials, we want all independent muscles. So we do not need 
% the muscles from both legs, since we assume bilateral symmetry, but want
% all muscles from the back (indices 47:49).
musi_pol = musi;
NMuscle_pol = NMuscle/2;
CasADiFunctions_tracking_TD_MRI_FLV_spas

%% Passive joint torques
% We extract the parameters for the passive torques of the lower limbs and
% the trunk
pathPassiveTorques = [pathRepo,'\PassiveTorques'];
addpath(genpath(pathPassiveTorques));
PassiveTorquesData

%% Experimental data
joints = {'lower_torso_RX','lower_torso_RY','lower_torso_RZ',...
    'lower_torso_TX','lower_torso_TY','lower_torso_TZ','hip_flex_l',...
    'hip_add_l','hip_rot_l','hip_flex_r','hip_add_r','hip_rot_r',...
    'knee_flex_l','knee_flex_r','ankle_flex_l','ankle_flex_r',...
    'subt_angle_l','subt_angle_r','lumbar_pitch','lumbar_roll',...
    'lumbar_yaw'};
pathVariousFunctions = [pathRepo,'\VariousFunctions'];
addpath(genpath(pathVariousFunctions));

Qs = struct('mp',[]); 
GRF = struct('mp',[]); 
ID = struct('mp',[]);

trackSimRes = 93;
p = mfilename('fullpath');
[~,namescript,~] = fileparts(p);
pathresults = [pathRepo,'\Results'];
load([pathresults,'\',namescript,'\Results_tracking.mat']);
for mpi = 1:NPhases
    ResultsToTrack = Results_tracking.(['Case_',num2str(trackSimRes)]);
    % Re-create the joint angles
    Qs_temp = ResultsToTrack.Qs_opt(mpi).mp;
    % TODO
    Qs_temp(:,end+1:end+3) = zeros(size(Qs_temp,1),3);
    Qs_temp(:,jointi.radi) = Qs_temp(:,jointi.radi)*pi/180;
    % Re-create the ground reaction forces
    GRF_temp = ResultsToTrack.GRFs_opt(mpi).mp;
    % Re-create the ground reaction torques
    GRM_temp = ResultsToTrack.GRMs_opt(mpi).mp;
    % Re-create the torques
    ID_temp = ResultsToTrack.Ts_opt(mpi).mp;
    % TODO
    ID_temp(:,end+1:end+3) = zeros(size(ID_temp,1),3);
    % Apply a low-pass filter to the trajectories
    order = 4;
    cutoff_low = 20;
    fs=1/mean(diff(ResultsToTrack.tgrid(mpi).mp));
    [af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
    Qs_tempfilt = filtfilt(af,bf,Qs_temp);  
    Qs_temp_filt = [ResultsToTrack.tgrid(mpi).mp,Qs_tempfilt];
    GRF_tempfilt = filtfilt(af,bf,GRF_temp);  
    GRF_temp_filt = [ResultsToTrack.tgrid(mpi).mp,GRF_tempfilt];
    GRM_tempfilt = filtfilt(af,bf,GRM_temp);  
    GRM_temp_filt = [ResultsToTrack.tgrid(mpi).mp,GRM_tempfilt];
    ID_tempfilt = filtfilt(af,bf,ID_temp);  
    ID_temp_filt = [ResultsToTrack.tgrid(mpi).mp,ID_tempfilt];       
    % Adjust the time vector (shorter because N controls and N+1 states)
    time_opt(mpi).mp(1) = ResultsToTrack.tgrid(mpi).mp(1);
    time_opt(mpi).mp(2) = ResultsToTrack.tgrid(mpi).mp(end);                
    % Interpolation to number of mesh intervals
    step = (time_opt(mpi).mp(2)-time_opt(mpi).mp(1))/(N-1);
    interval = time_opt(mpi).mp(1):step:time_opt(mpi).mp(2);        
    Qs(mpi).mp.allinterpfilt = ...
        interp1(Qs_temp_filt(:,1),Qs_temp_filt,interval);
    Qs(mpi).mp.colheaders = ['Time',Results_tracking.colheaders.joints];
    ID(mpi).mp.allinterp = interp1(ID_temp_filt(:,1),ID_temp_filt,interval);        
    GRF(mpi).mp.val.allinterp = ...
        interp1(GRF_temp_filt(:,1),GRF_temp_filt,interval);
    GRF(mpi).mp.val.all = GRF(mpi).mp.val.allinterp;
    GRF(mpi).mp.MorGF.allinterp = ...
        interp1(GRM_temp_filt(:,1),GRM_temp_filt,interval);           
end  


%% Bounds
pathBounds = [pathRepo,'\Bounds'];
addpath(genpath(pathBounds));
bounds = struct('mp',[]); 
scaling = struct('mp',[]); 
Qs_CP = struct('mp',[]);
trunki = 1;
for mpi = 1:NPhases
    % Bounds accounting for both TD and CD walking patterns
    pathIK = [pathData,'\KS\KS_walking_average_r_MRI_extROM.mot'];
    Qs_CP(mpi).mp = getIK_MRI(pathIK,joints);           
    step = (Qs_CP(mpi).mp.time(end)-Qs_CP(mpi).mp.time(1))/(N-1);
    interval = Qs_CP(mpi).mp.time(1):step:Qs_CP(mpi).mp.time(end);        
    Qs_CP(mpi).mp.allinterpfilt = interp1(Qs_CP(mpi).mp.allfilt(:,1),...
        Qs_CP(mpi).mp.allfilt,interval);  
    [bounds(mpi).mp,scaling(mpi).mp] = ...
        getBounds_tracking_TD_MRI_FLV_Syn3_spas_TDCP(Qs(mpi).mp,...
        Qs_CP(mpi).mp,NMuscle,NMuscle_FLV,nq,jointi,GRF(mpi).mp,N,NSyn,...
        synW,dev,NMuscle_Spas,trunki,W);
end

%% Initial guess
pathIG = [pathRepo,'\IG'];
addpath(genpath(pathIG));
% Data-informed initial guess
guess = struct('mp',[]); 
Qs_CP = struct('mp',[]);
Qs_ig = struct('mp',[]);
for mpi = 1:NPhases
    if IGi == 1 
        % IG is walking pattern to track (TD)
        guess(mpi).mp = getGuess_DI_tracking_TD_MRI_FLV_Syn2_spas(...
            Qs(mpi).mp,nq,N,NMuscle,NMuscle_FLV,jointi,scaling(mpi).mp,...
            NSyn,synW,NMuscle_Spas,trunki);
    elseif IGi == 2 
        guess(mpi).mp = getGuess_QR_tracking_TD_MRI_FLV_Syn_spas(nq,N,...
            NMuscle,NMuscle_FLV,jointi,scaling(mpi).mp,NSyn,v_tgt,trunki,...
            NMuscle_Spas);
    elseif IGi == 3
        % IG is CP walking pattern
        pathIK = [pathData,'\KS\KS_walking_average_r.mot'];
        Qs_CP(mpi).mp = getIK_MRI(pathIK,joints);           
        step = (Qs_CP(mpi).mp.time(end)-Qs_CP(mpi).mp.time(1))/(N-1);
        interval = Qs_CP(mpi).mp.time(1):step:Qs_CP(mpi).mp.time(end);        
        Qs_CP(mpi).mp.allinterpfilt = interp1(Qs_CP(mpi).mp.allfilt(:,1),...
            Qs_CP(mpi).mp.allfilt,interval);        
        guess(mpi).mp = getGuess_DI_tracking_TD_MRI_FLV_Syn3_spas(...
            Qs_CP(mpi).mp,nq,N,NMuscle,NMuscle_FLV,jointi,scaling(mpi).mp,...
            NSyn,synW,NMuscle_Spas,trunki,W);   
    elseif IGi == 4 
        % IG is walking pattern to track (TD)
        guess(mpi).mp = getGuess_DI_tracking_TD_MRI_FLV_Syn3_spas(...
            Qs(mpi).mp,nq,N,NMuscle,NMuscle_FLV,jointi,scaling(mpi).mp,...
            NSyn,synW,NMuscle_Spas,trunki,W);
    elseif IGi == 5
        % IG is certain simulated movement        
        savename_ig = ['_c',num2str(IGc),'a'];
        pathIK_ig = [pathRepo,'\Results\',namescript,...
            '\IK',savename_ig,'_HS.mot'];
        Qs_ig(mpi).mp = getIK_MRI(pathIK_ig,joints);           
        step = (Qs_ig(mpi).mp.time(end)-Qs_ig(mpi).mp.time(1))/(N-1);        
        interval = Qs_ig(mpi).mp.time(1):step:Qs_ig(mpi).mp.time(end);        
        Qs_ig(mpi).mp.allinterpfilt = interp1(Qs_ig(mpi).mp.allfilt(:,1),...
            Qs_ig(mpi).mp.allfilt,interval);  
        guess(mpi).mp =  getGuess_DI_tracking_TD_MRI_FLV_Syn3_spas(...
            Qs_ig(mpi).mp,nq,N,NMuscle,NMuscle_FLV,jointi,scaling(mpi).mp,...
            NSyn,synW,NMuscle_Spas,trunki,W);  
    end
end
% This allows visualizing the initial guess and the bounds
if checkBoundsIG
    pathPlots = [pathRepo,'\Plots'];
    addpath(genpath(pathPlots));
    for mpi = 1:NPhases
        plot_BoundsVSInitialGuess_tracking_TD_env_MRI_MP_FLV
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
    if W.Syn ~= 0
        % Synergy weights: left side
        syn_wl          = MX.sym('syn_wl',NMuscle/2*NSyn);
        w               = [w {syn_wl}];
        lbw             = [lbw; bounds(1).mp.synw.lower.l(:)];
        ubw             = [ubw; bounds(1).mp.synw.upper.l(:)];
        w0              = [w0;  guess(1).mp.synw.l(:)];
        % Synergy weights: right side
        syn_wr          = MX.sym('syn_wr',NMuscle/2*NSyn);
        w               = [w {syn_wr}];
        lbw             = [lbw; bounds(1).mp.synw.lower.r(:)];
        ubw             = [ubw; bounds(1).mp.synw.upper.r(:)];
        w0              = [w0;  guess(1).mp.synw.r(:)];  
    end
    % loop over phases
    for mpi = 1:NPhases  
        % Define states at first mesh point
        % Muscle activations
        a0              = MX.sym('a0',NMact);
        w               = [w {a0}];
        lbw             = [lbw; bounds(mpi).mp.a.lower'];
        ubw             = [ubw; bounds(mpi).mp.a.upper'];
        w0              = [w0;  guess(mpi).mp.a(1,:)'];
        % Muscle-tendon forces
        FTtilde0        = MX.sym('FTtilde0',NMuscle_FLV);
        w               = [w {FTtilde0}];
        lbw             = [lbw; bounds(mpi).mp.FTtilde.lower'];
        ubw             = [ubw; bounds(mpi).mp.FTtilde.upper'];
        w0              = [w0;  guess(mpi).mp.FTtilde(1,:)'];
        % Qs and Qdots
        X0              = MX.sym('X0',2*nq.all);
        w               = [w {X0}];    
        lbw             = [lbw; bounds(mpi).mp.QsQdots.lower(1,:)'];
        ubw             = [ubw; bounds(mpi).mp.QsQdots.upper(1,:)'];    
        w0              = [w0;  guess(mpi).mp.QsQdots(1,:)']; 
        if spasi == 1
            % Muscle activations from muscle-tendon force feedback
            a_Ff0           = MX.sym('a_Ff0',NMuscle_Spas);
            w               = [w {a_Ff0}];
            lbw             = [lbw; bounds(mpi).mp.a_Ff.lower'];
            ubw             = [ubw; bounds(mpi).mp.a_Ff.upper'];
            w0              = [w0;  guess(mpi).mp.a_Ff(1,:)'];    
            % Muscle activations from time derivative of muscle-tendon force feedback
            a_dFf0           = MX.sym('a_dFf0',NMuscle_Spas);
            w               = [w {a_dFf0}];
            lbw             = [lbw; bounds(mpi).mp.a_dFf.lower'];
            ubw             = [ubw; bounds(mpi).mp.a_dFf.upper'];
            w0              = [w0;  guess(mpi).mp.a_dFf(1,:)'];   
        end
        % Back activations
        a_b0            = MX.sym('a_b0',nq.trunk);
        w               = [w {a_b0}];
        lbw             = [lbw; bounds(mpi).mp.a_b.lower'];
        ubw             = [ubw; bounds(mpi).mp.a_b.upper'];
        w0              = [w0;  guess(mpi).mp.a_b(1,:)'];    
        % We pre-allocate some of the states so that we can provide an
        % expression for the distance traveled
        for k=0:N
            Xk{k+1,1} = MX.sym(['X_' num2str(k+1)], 2*nq.all);
        end   
        % "Lift" initial conditions
        ak          = a0;
        FTtildek    = FTtilde0;
        Xk{1,1}     = X0;
        if spasi == 1
            a_Ffk       = a_Ff0; 
            a_dFfk      = a_dFf0;
        end
        a_bk        = a_b0; 
        % Provide expression for the distance traveled
        pelvis_tx0 = Xk{1,1}(2*jointi.pelvis.tx-1,1).*...
            scaling(mpi).mp.QsQdots(2*jointi.pelvis.tx-1); % initial position pelvis_tx  
         % final position pelvis_tx: N and not N+1 to match experimental data 
        pelvis_txf = Xk{N,1}(2*jointi.pelvis.tx-1,1).*...
            scaling(mpi).mp.QsQdots(2*jointi.pelvis.tx-1);   
        dist_trav_tot = pelvis_txf-pelvis_tx0;% distance traveled  
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
            dFTtildek           = MX.sym(['dFTtilde_' num2str(k)],NMuscle_FLV);
            w                   = [w {dFTtildek}];
            lbw                 = [lbw; bounds(mpi).mp.dFTtilde.lower'];
            ubw                 = [ubw; bounds(mpi).mp.dFTtilde.upper'];
            w0                  = [w0; guess(mpi).mp.dFTtilde(k+1,:)'];  
            % Time derivative of Qdots (states) 
            Ak                  = MX.sym(['A_' num2str(k)],nq.all);
            w                   = [w {Ak}];
            lbw                 = [lbw; bounds(mpi).mp.Qdotdots.lower'];
            ubw                 = [ubw; bounds(mpi).mp.Qdotdots.upper'];
            w0                  = [w0; guess(mpi).mp.Qdotdots(k+1,:)'];
            % Back excitations
            e_bk                = MX.sym(['e_b_' num2str(k)],nq.trunk);
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
                    MX.sym(['FTtilde_' num2str(k) '_' num2str(j)],NMuscle_FLV);
                w            = {w{:}, FTtildekj{j}};
                lbw          = [lbw; bounds(mpi).mp.FTtilde.lower'];
                ubw          = [ubw; bounds(mpi).mp.FTtilde.upper'];
                w0           = [w0;  guess(mpi).mp.FTtilde(k+1,:)'];
            end
            % Qs and Qdots        
            Xkj = {};
            for j=1:d
                Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)],2*nq.all);
                w      = {w{:}, Xkj{j}};
                if boundsi == 1
                    lbw    = [lbw; bounds(mpi).mp.QsQdots.lower'];
                    ubw    = [ubw; bounds(mpi).mp.QsQdots.upper'];
                elseif boundsi == 2 ||  boundsi == 3 || boundsi == 4 || boundsi == 5 || boundsi == 6 
                    lbw    = [lbw; bounds(mpi).mp.QsQdots.lower(k+1,:)'];
                    ubw    = [ubw; bounds(mpi).mp.QsQdots.upper(k+1,:)'];
                end
                w0     = [w0;  guess(mpi).mp.QsQdots(k+1,:)'];
            end   
            if spasi == 1
                % Muscle activations from muscle-tendon force feedback
                a_Ffkj = {};
                for j=1:d
                    a_Ffkj{j}= MX.sym(['	a_Ff_' num2str(k) '_' num2str(j)],NMuscle_Spas);
                    w       = {w{:}, a_Ffkj{j}};
                    lbw     = [lbw; bounds(mpi).mp.a_Ff.lower'];
                    ubw     = [ubw; bounds(mpi).mp.a_Ff.upper'];
                    w0      = [w0;  guess(mpi).mp.a_Ff(k+1,:)'];
                end
                % Muscle activations from time derivative of muscle-tendon force feedback
                a_dFfkj = {};
                for j=1:d
                    a_dFfkj{j}= MX.sym(['	a_dFf_' num2str(k) '_' num2str(j)],NMuscle_Spas);
                    w       = {w{:}, a_dFfkj{j}};
                    lbw     = [lbw; bounds(mpi).mp.a_dFf.lower'];
                    ubw     = [ubw; bounds(mpi).mp.a_dFf.upper'];
                    w0      = [w0;  guess(mpi).mp.a_dFf(k+1,:)'];
                end
            end
            % Back activations
            a_bkj = {};
            for j=1:d
                a_bkj{j}= MX.sym(['	a_b_' num2str(k) '_' num2str(j)],nq.trunk);
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
            Tauk_abs = Tk(ground_pelvisi,1);
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
            if W.Syn == 0
                asynk = ak;
            else
                alk = f_SynergyProduct(syn_wl,ak(1:NSyn,1));
                ark = f_SynergyProduct(syn_wr,ak(NSyn+1:2*NSyn,1));            
                asynk = [alk;ark];
            end                   
            % Muscle activations come from supra-spinal input and reflex feedback
            a_totk = MX(NMuscle,1); 
            if spasi == 1
                a_totk(musi_noSpas) = asynk(musi_noSpas); % supra-spinal inputs
                a_totk(musi_Spas) = asynk(musi_Spas) + a_Ffk + a_dFfk; % reflexes     
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
            FTk = MX(NMuscle,1);
            FTk(musi_FLV) = FTk_FLV;
            if ~isempty(musi_noFLV)
                FTk(musi_noFLV) = FTk_noFLV;
            end
            % Get metabolic energy rate if in the cost function   
            if W.E ~= 0    
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
                    MTparameters_m_FLV(1,:)',body_mass,10);
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
                a_Ffk_end           = D(1)*a_Ffk;
                a_dFfk_end          = D(1)*a_dFfk;
            end
            a_bk_end            = D(1)*a_bk; 
            for j=1:d
                % Expression for the state derivatives at the collocation point
                xp_nsc          = C(1,j+1)*Xk_nsc;
                FTtildep_nsc    = C(1,j+1)*FTtildek_nsc;
                ap              = C(1,j+1)*ak;
                if spasi == 1
                    a_Ffp           = C(1,j+1)*a_Ffk;
                    a_dFfp          = C(1,j+1)*a_dFfk;
                end
                a_bp            = C(1,j+1)*a_bk;
                for r=1:d
                    xp_nsc       = xp_nsc + C(r+1,j+1)*Xkj_nsc{r};
                    FTtildep_nsc = FTtildep_nsc + C(r+1,j+1)*FTtildekj_nsc{r};
                    ap           = ap + C(r+1,j+1)*akj{r};
                    if spasi == 1
                        a_Ffp        = a_Ffp + C(r+1,j+1)*a_Ffkj{r};
                        a_dFfp       = a_dFfp + C(r+1,j+1)*a_dFfkj{r};
                    end
                    a_bp         = a_bp + C(r+1,j+1)*a_bkj{r};
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
                lbg     = [lbg; zeros(NMuscle_FLV,1)];
                ubg     = [ubg; zeros(NMuscle_FLV,1)];
                if spasi == 1
                    % Spindle dynamics (explicit formulation)
                    % Force feedback            
                    da_Ffdt = f_spindleDynamics(a_Ffkj{j},...
                        FTtildek_nsc(musi_SpasInFLV),tauFf.all,gFf.all,bspas,...
                        threshold_Ff_gait.all);
                    g       = {g{:}, (h*da_Ffdt - a_Ffp)};
                    lbg     = [lbg; zeros(NMuscle_Spas,1)];
                    ubg     = [ubg; zeros(NMuscle_Spas,1)]; 
                    % Time derivative of force feedback feedback            
                    da_dFfdt = f_spindleDynamics(a_dFfkj{j},...
                        dFTtildek(musi_SpasInFLV).*scaling(mpi).mp.dFTtilde,...
                        taudFf.all,gdFf.all,bspas,threshold_dFf_gait.all);
                    g       = {g{:}, (h*da_dFfdt - a_dFfp)};
                    lbg     = [lbg; zeros(NMuscle_Spas,1)];
                    ubg     = [ubg; zeros(NMuscle_Spas,1)];
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
                g       = {g{:}, (h*xj_nsc - xp_nsc)./(scaling(mpi).mp.QsQdots')};
                lbg     = [lbg; zeros(2*nq.all,1)];
                ubg     = [ubg; zeros(2*nq.all,1)];   
                % Back activation dynamics (explicit formulation)   
                dadt    = f_TrunkActivationDynamics(e_bk,a_bkj{j});
                g       = {g{:}, (h*dadt - a_bp)./scaling(mpi).mp.a_b};
                lbg     = [lbg; zeros(nq.trunk,1)];
                ubg     = [ubg; zeros(nq.trunk,1)];  
                % Add contribution to the end state
                Xk_nsc_end       = Xk_nsc_end + D(j+1)*Xkj_nsc{j};
                FTtildek_nsc_end = FTtildek_nsc_end + D(j+1)*FTtildekj_nsc{j};
                ak_end           = ak_end + D(j+1)*akj{j};
                if spasi == 1
                    a_Ffk_end   = a_Ffk_end + D(j+1)*a_Ffkj{j};   
                    a_dFfk_end  = a_dFfk_end + D(j+1)*a_dFfkj{j}; 
                end                 
               a_bk_end     = a_bk_end + D(j+1)*a_bkj{j};  
                % Add contribution to quadrature function
                % Tracking terms
                % Qs
                Qs_costk = B(j+1)*(f_JNq_bpty(Xk{k+1,1}(Qsi(residual_bptyi))-...
                        Qs(mpi).mp.allinterpfilt(k+1,residual_bptyi+1)'... 
                        ./scaling(mpi).mp.Qs(residual_bptyi)'))*h;
                % GRF
                GRF_costk = B(j+1)*(f_J6((Tk(GRFi.all,1)./...
                    scaling(mpi).mp.GRF')-...
                    GRF(mpi).mp.val.allinterp(k+1,2:end)'./...
                    scaling(mpi).mp.GRF'))*h;
                % GRT
                GRT_costk = B(j+1)*(f_J6((Tk(GRMi.all,1)./...
                    scaling(mpi).mp.GRM')-...
                    GRF(mpi).mp.MorGF.allinterp(k+1,2:end)'./...
                    scaling(mpi).mp.GRM'))*h;
                % Joint torques
                ID_costk = B(j+1)*...
                    (f_JNq_act((Tk(residuals_acti,1)./scaling(mpi).mp.T(1)')-...
                    ID(mpi).mp.allinterp(k+1,2+nq.abs:end)'./...
                    scaling(mpi).mp.T(1)))*h;
                % All combined
                track_terms = W.Qs*Qs_costk + W.GRF*GRF_costk +...
                    W.GRM*GRT_costk + W.ID_act*ID_costk;
                % Motor control terms
                % Activations
                 a_costk = B(j+1)*(f_JNM_exp(asynk,exp_A))*h;
                % Metabolic energy rate
                if W.E == 0
                    mE_costk = 0;
                else
                    mE_costk = B(j+1)*(f_JNM_FLVexp(e_tot,exp_E))/body_mass*h;
                end
                trunk_costk = B(j+1)*(f_Jnq_trunk(e_bk))*h;     
                % Joint accelerations
                Qdotdot_costk = B(j+1)*(f_JNq_all(Ak))*h;
                % Passive torques
                passT_costk = B(j+1)*(f_JNq_act(Tau_passk_all))*h;
                % Time derivative of muscle activations / muscle-tendon
                % forces
                vA_costk = B(j+1)*(f_JNMact(vAk))*h;
                dFTtilde_costk = B(j+1)*(f_JNM_FLV(dFTtildek))*h;
                % All combined
                mc_terms = W.a*a_costk + W.E*mE_costk +...
                    W.qdotdot*Qdotdot_costk + W.passT*passT_costk + ...
                    W.Trunk*trunk_costk + W.u*vA_costk + W.u*dFTtilde_costk;   
                % Quadrature function   
                dist_norm = dist_trav_tot;
                J = J + track_terms + 1/dist_norm*mc_terms;                  
            end                              
            % Add path constraints
            % Pelvis residuals
            Tauk_exp = (ID(mpi).mp.allinterp(k+1,2:7)')./scaling(mpi).mp.T(1);
            g   = {g{:},(Tauk_abs)./scaling(mpi).mp.T(1)};
            if dev.pr == 1000
                % Null pelvis residuals
                lbg = [lbg; zeros(size(Tauk_abs,1),1)];        
                ubg = [ubg; zeros(size(Tauk_abs,1),1)];
            else
                % Constrained deviation from experimental pelvis residuals
                lbg = [lbg; Tauk_exp - dev.pr/100*abs(Tauk_exp)];        
                ubg = [ubg; Tauk_exp + dev.pr/100*abs(Tauk_exp)];      
            end            
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
            g       = {g{:},Tk(backi,1)/scaling(mpi).mp.TrunkTau - a_bk};
            lbg     = [lbg; zeros(nq.trunk,1)];
            ubg 	= [ubg; zeros(nq.trunk,1)];
            % Activation dynamics (implicit formulation)
            tact = 0.015;
            tdeact = 0.06;
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
            lbg             = [lbg; zeros(NMuscle_FLV,1)];
            ubg             = [ubg; zeros(NMuscle_FLV,1)];      
            % Total activations (supra-spinal and reflexes) cannot exceed 1
            if spasi == 1
                g   = {g{:},a_totk(musi_Spas)};
                lbg = [lbg; zeros(NMuscle_Spas,1)];
                ubg = [ubg; ones(NMuscle_Spas,1)]; 
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
                ak              = MX.sym(['a_' num2str(k+1)],NMact);
                w               = {w{:}, ak};
                lbw             = [lbw; bounds(mpi).mp.a.lower'];
                ubw             = [ubw; bounds(mpi).mp.a.upper'];
                w0              = [w0;  guess(mpi).mp.a(k+2,:)'];
                % Muscle-tendon forces
                FTtildek        = MX.sym(['FTtilde_' num2str(k+1)],NMuscle_FLV);
                w               = {w{:}, FTtildek};
                lbw             = [lbw; bounds(mpi).mp.FTtilde.lower'];
                ubw             = [ubw; bounds(mpi).mp.FTtilde.upper'];
                w0              = [w0;  guess(mpi).mp.FTtilde(k+2,:)'];    
                % Qs and Qdots
                w = {w{:}, Xk{k+2,1}};
                if boundsi == 1
                    lbw = [lbw; bounds(mpi).mp.QsQdots.lower'];
                    ubw = [ubw; bounds(mpi).mp.QsQdots.upper'];
                elseif boundsi == 2 ||  boundsi == 3 || boundsi == 4 || boundsi == 5 || boundsi == 6
                    lbw = [lbw; bounds(mpi).mp.QsQdots.lower(k+2,:)'];
                    ubw = [ubw; bounds(mpi).mp.QsQdots.upper(k+2,:)'];
                end 
                w0 = [w0;  guess(mpi).mp.QsQdots(k+2,:)'];
                if spasi == 1
                    % Muscle activations from muscle-tendon force feedback
                    a_Ffk	= MX.sym(['a_Ff_' num2str(k+1)],NMuscle_Spas);
                    w       = {w{:}, a_Ffk};
                    lbw 	= [lbw; bounds(mpi).mp.a_Ff.lower'];
                    ubw 	= [ubw; bounds(mpi).mp.a_Ff.upper'];
                    w0  	= [w0;  guess(mpi).mp.a_Ff(k+2,:)'];
                    % Muscle activations from time derivative of muscle-tendon
                    % force feedback
                    a_dFfk	= MX.sym(['a_dFf_' num2str(k+1)],NMuscle_Spas);
                    w       = {w{:}, a_dFfk};
                    lbw  	= [lbw; bounds(mpi).mp.a_dFf.lower'];
                    ubw  	= [ubw; bounds(mpi).mp.a_dFf.upper'];
                    w0    	= [w0;  guess(mpi).mp.a_dFf(k+2,:)'];
                end
                % Back activations
                a_bk            = MX.sym(['a_b_' num2str(k+1)],nq.trunk);
                w               = {w{:}, a_bk};
                lbw             = [lbw; bounds(mpi).mp.a_b.lower'];
                ubw             = [ubw; bounds(mpi).mp.a_b.upper'];
                w0              = [w0;  guess(mpi).mp.a_b(k+2,:)'];
            else
                % Muscle activations
                ak              = MX.sym(['a_' num2str(k+1)],NMact);
                w               = {w{:}, ak};
                lbw             = [lbw; bounds(mpi).mp.a.lower'];
                ubw             = [ubw; bounds(mpi).mp.a.upper'];
                w0              = [w0;  guess(mpi).mp.a(end,:)'];
                % Muscle-tendon forces
                FTtildek        = MX.sym(['FTtilde_' num2str(k+1)],NMuscle_FLV);
                w               = {w{:}, FTtildek};
                lbw             = [lbw; bounds(mpi).mp.FTtilde.lower'];
                ubw             = [ubw; bounds(mpi).mp.FTtilde.upper'];
                w0              = [w0;  guess(mpi).mp.FTtilde(end,:)'];    
                % Qs and Qdots
                w = {w{:}, Xk{k+2,1}};
                if boundsi == 1
                    lbw	= [lbw; bounds(mpi).mp.QsQdots.lower'];
                    ubw	= [ubw; bounds(mpi).mp.QsQdots.upper'];
                elseif boundsi == 2 ||  boundsi == 3 || boundsi == 4 || boundsi == 5 || boundsi == 6 
                    lbw	= [lbw; bounds(mpi).mp.QsQdots.lower(end,:)'];
                    ubw = [ubw; bounds(mpi).mp.QsQdots.upper(end,:)'];
                end
                w0 = [w0;  guess(mpi).mp.QsQdots(end,:)'];
                if spasi == 1
                    % Muscle activations from muscle-tendon force feedback
                    a_Ffk	= MX.sym(['a_Ff_' num2str(k+1)],NMuscle_Spas);
                    w       = {w{:}, a_Ffk};
                    lbw   	= [lbw; bounds(mpi).mp.a_Ff.lower'];
                    ubw    	= [ubw; bounds(mpi).mp.a_Ff.upper'];
                    w0    	= [w0;  guess(mpi).mp.a_Ff(end,:)'];
                    % Muscle activations from time derivative of muscle-tendon
                    % force feedback
                    a_dFfk 	= MX.sym(['a_dFf_' num2str(k+1)],NMuscle_Spas);
                    w     	= {w{:}, a_dFfk};
                    lbw    	= [lbw; bounds(mpi).mp.a_dFf.lower'];
                    ubw    	= [ubw; bounds(mpi).mp.a_dFf.upper'];
                    w0     	= [w0;  guess(mpi).mp.a_dFf(end,:)'];
                end
                % Back activations
                a_bk            = MX.sym(['a_b_' num2str(k+1)], nq.trunk);
                w               = {w{:}, a_bk};
                lbw             = [lbw; bounds(mpi).mp.a_b.lower'];
                ubw             = [ubw; bounds(mpi).mp.a_b.upper'];
                w0              = [w0;  guess(mpi).mp.a_b(end,:)'];
            end
            % Rescale variables to impose equality constraints
            Xk_end = (Xk_nsc_end)./scaling(mpi).mp.QsQdots';
            FTtildek_end = (FTtildek_nsc_end)./scaling(mpi).mp.FTtilde';
            % Add equality constraints (next interval starts with end values of 
            % states from previous interval)
            g   = {g{:}, Xk_end-Xk{k+2,1}, FTtildek_end-FTtildek, ak_end-ak};
            lbg = [lbg; zeros(2*nq.all + NMact + NMuscle_FLV,1)];
            ubg = [ubg; zeros(2*nq.all + NMact + NMuscle_FLV,1)]; 
            if spasi == 1
                g   = {g{:}, a_Ffk_end-a_Ffk, a_dFfk_end-a_dFfk};
                lbg = [lbg; zeros(NMuscle_Spas + NMuscle_Spas,1)];
                ubg = [ubg; zeros(NMuscle_Spas + NMuscle_Spas,1)];   
            end
            g   = {g{:}, a_bk_end-a_bk};
            lbg = [lbg; zeros(nq.trunk,1)];
            ubg = [ubg; zeros(nq.trunk,1)];
        end   
        % Periodic constraints
        % Periodicity on joint states
        % All states but Q pelvis_tx, trunk included
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
        lbg = [lbg; zeros(NMuscle_FLV,1)];
        ubg = [ubg; zeros(NMuscle_FLV,1)];
        % Speed constraint
        vel_aver_tot = dist_trav_tot/tf;
        g   = {g{:}, vel_aver_tot - v_tgt};
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
    pathresults = [pathRepo,'\Results'];
    if ~(exist([pathresults,'\',namescript],'dir')==7)
        mkdir(pathresults,namescript);
    end
    if (exist([pathresults,'\',namescript,'\D',savename],'file')==2)
        delete ([pathresults,'\',namescript,'\D',savename])
    end 
    diary([pathresults,'\',namescript,'\D',savename]); 
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
    save([pathresults,'\',namescript,'\w',savename],'w_opt');
    save([pathresults,'\',namescript,'\g',savename],'g_opt');
    save([pathresults,'\',namescript,'\s',savename],'setup');
    save([pathresults,'\',namescript,'\stats',savename],'stats');
end

%% Analyze results
if analyseResults
    %% Load results
    if loadResults
        p = mfilename('fullpath');
        [~,namescript,~] = fileparts(p);
        pathresults = [pathRepo,'\Results'];
        load([pathresults,'\',namescript,'\w',savename]);
        load([pathresults,'\',namescript,'\g',savename]);
        load([pathresults,'\',namescript,'\s',savename]);
        load([pathresults,'\',namescript,'\stats',savename]);
    end
    
    %% Extract results
    % All optimized design variables are saved in a single column vector      
    % Number of design variables  
    NControls = NMact+NMuscle_FLV+nq.all;
    NStates = NMact+NMuscle_FLV+2*nq.all;
    NParameters = 0;
    NParameters = NParameters + 1;        
    if W.Syn ~= 0        
        NParameters = NParameters + NMuscle*NSyn;   
    end
    if spasi == 1
        NStates = NStates+2*NMuscle_Spas;
    end    
    NControls = NControls+nq.trunk;
    NStates = NStates+nq.trunk;
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
    if W.Syn ~= 0     
        syn_wl_opt = w_opt(np_acc:np_acc+NMuscle/2*NSyn-1);
        syn_wr_opt = w_opt(np_acc+NMuscle/2*NSyn:NParameters);
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
    syn_a_opt_l = struct('mp',[]);
    syn_a_opt_r = struct('mp',[]);
    a_Ff_opt = struct('mp',[]);
    a_dFf_opt = struct('mp',[]);
    for mpi = 1:NPhases 
        % Mesh points
        % Muscle activations and muscle-tendon forces
        a_opt(mpi).mp = zeros(N+1,NMact);
        for i = 1:NMact
            a_opt(mpi).mp(:,i) = w_opt(Nw_acc+NParameters+i:Nwl:Nw_acc+Nw);
        end
        FTtilde_opt(mpi).mp = zeros(N+1,NMuscle_FLV);
        for i = 1:NMuscle_FLV
            FTtilde_opt(mpi).mp(:,i) = w_opt(Nw_acc+NParameters+NMact+i:Nwl:Nw_acc+Nw);
        end        
        % Qs and Qdots
        q_opt(mpi).mp = zeros(N+1,nq.all);
        qdot_opt(mpi).mp = zeros(N+1,nq.all);
        count = 0;
        for i = 1:2:2*nq.all
            count = count +1;
            q_opt(mpi).mp(:,count)      = w_opt(Nw_acc+NParameters+NMact+NMuscle_FLV+i:Nwl:Nw_acc+Nw);
            qdot_opt(mpi).mp(:,count)   = w_opt(Nw_acc+NParameters+NMact+NMuscle_FLV+i+1:Nwl:Nw_acc+Nw);
        end
        qqdot_opt(mpi).mp = zeros(N+1,2*nq.all); 
        qqdot_opt(mpi).mp(:,1:2:end) = q_opt(mpi).mp;
        qqdot_opt(mpi).mp(:,2:2:end) = qdot_opt(mpi).mp;
        tempi = Nw_acc+NParameters+NMact+NMuscle_FLV+2*nq.all;
        if spasi == 1
            % Muscle activations from muscle-tendon force feedback and time derivative 
            % of muscle-tendon force feedback
            a_Ff_opt(mpi).mp = zeros(N+1,NMuscle_Spas); 
            a_dFf_opt(mpi).mp = zeros(N+1,NMuscle_Spas);
            for i = 1:NMuscle_Spas
                a_Ff_opt(mpi).mp(:,i) =  w_opt(tempi+i:Nwl:Nw_acc+Nw);
                a_dFf_opt(mpi).mp(:,i) = w_opt(tempi+NMuscle_Spas+i:Nwl:Nw_acc+Nw);
            end 
            tempi = tempi + 2*NMuscle_Spas;            
        end
        % Back activations
        a_b_opt(mpi).mp = zeros(N+1,nq.trunk);
        for i = 1:nq.trunk
            a_b_opt(mpi).mp(:,i) = w_opt(tempi+i:Nwl:Nw_acc+Nw);
        end           
        % Time derivative of muscle activations and muscle-tendon forces
        vA_opt(mpi).mp = zeros(N,NMact);
        for i = 1:NMact
            vA_opt(mpi).mp(:,i) = w_opt(Nw_acc+NParameters+NStates+i:Nwl:Nw_acc+Nw);
        end
        dFTtilde_opt(mpi).mp = zeros(N,NMuscle_FLV);
        for i = 1:NMuscle_FLV
            dFTtilde_opt(mpi).mp(:,i)   = w_opt(Nw_acc+NParameters+NStates+...
                NMact+i:Nwl:Nw_acc+Nw);
        end
        % Time derivative of joint velocities
        qdotdot_opt(mpi).mp = zeros(N,nq.all);
        for i = 1:nq.all
            qdotdot_opt(mpi).mp(:,i)    = w_opt(Nw_acc+NParameters+NStates+...
                NMact+NMuscle_FLV+i:Nwl:Nw_acc+Nw);
        end
        tempi = Nw_acc+NParameters+NStates+NMact+NMuscle_FLV+nq.all;
        % Back excitations
        e_b_opt(mpi).mp = zeros(N,nq.trunk);
        for i = 1:nq.trunk
            e_b_opt(mpi).mp(:,i) = w_opt(tempi+i:Nwl:Nw_acc+Nw); 
        end                
        % Collocation points
        % Muscle activations
        a_opt_ext(mpi).mp=zeros(N*(d+1)+1,NMact);
        a_opt_ext(mpi).mp(1:(d+1):end,:)= a_opt(mpi).mp;
        for nmusi=1:NMact
            a_opt_ext(mpi).mp(2:(d+1):end,nmusi) = w_opt(Nw_acc+Nwm+nmusi:Nwl:Nw_acc+Nw);
            a_opt_ext(mpi).mp(3:(d+1):end,nmusi) = w_opt(Nw_acc+Nwm+NMact+nmusi:Nwl:Nw_acc+Nw);
            a_opt_ext(mpi).mp(4:(d+1):end,nmusi) = w_opt(Nw_acc+Nwm+NMact+NMact+nmusi:Nwl:Nw_acc+Nw);
        end
        % Muscle activations at collocation points only
        a_opt_ext_col(mpi).mp = zeros(N*d,NMact); 
        for nmusi=1:NMact
            a_opt_ext_col(mpi).mp(1:d:end,nmusi) = w_opt(Nw_acc+Nwm+nmusi:Nwl:Nw_acc+Nw);
            a_opt_ext_col(mpi).mp(2:d:end,nmusi) = w_opt(Nw_acc+Nwm+NMact+nmusi:Nwl:Nw_acc+Nw);
            a_opt_ext_col(mpi).mp(3:d:end,nmusi) = w_opt(Nw_acc+Nwm+NMact+NMact+nmusi:Nwl:Nw_acc+Nw);   
        end 
        % Muscle-tendon forces
        FTtilde_opt_ext(mpi).mp=zeros(N*(d+1)+1,NMuscle_FLV);
        FTtilde_opt_ext(mpi).mp(1:(d+1):end,:)= FTtilde_opt(mpi).mp;
        for nmusi=1:NMuscle_FLV
            FTtilde_opt_ext(mpi).mp(2:(d+1):end,nmusi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+nmusi:Nwl:Nw_acc+Nw);
            FTtilde_opt_ext(mpi).mp(3:(d+1):end,nmusi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+NMuscle_FLV+nmusi:Nwl:Nw_acc+Nw);
            FTtilde_opt_ext(mpi).mp(4:(d+1):end,nmusi) = ...
                w_opt(Nw_acc+Nwm+d*NMact+NMuscle_FLV+NMuscle_FLV+nmusi:Nwl:Nw_acc+Nw);
        end
        % Qs and Qdots
        q_opt_ext(mpi).mp=zeros(N*(d+1)+1,nq.all);
        q_opt_ext(mpi).mp(1:(d+1):end,:)= q_opt(mpi).mp;
        q_dot_opt_ext(mpi).mp=zeros(N*(d+1)+1,nq.all);
        q_dot_opt_ext(mpi).mp(1:(d+1):end,:)= qdot_opt(mpi).mp;
        nqi_col = 1:2:2*nq.all;
        for nqi=1:nq.all
            nqi_q = nqi_col(nqi);
            q_opt_ext(mpi).mp(2:(d+1):end,nqi) = w_opt(Nw_acc+Nwm+d*NMact+d*NMuscle_FLV+nqi_q:Nwl:Nw_acc+Nw);   
            q_opt_ext(mpi).mp(3:(d+1):end,nqi) = w_opt(Nw_acc+Nwm+d*NMact+d*NMuscle_FLV+2*nq.all+nqi_q:Nwl:Nw_acc+Nw);  
            q_opt_ext(mpi).mp(4:(d+1):end,nqi) = w_opt(Nw_acc+Nwm+d*NMact+d*NMuscle_FLV+2*nq.all+2*nq.all+nqi_q:Nwl:Nw_acc+Nw);  
            q_dot_opt_ext(mpi).mp(2:(d+1):end,nqi) = w_opt(Nw_acc+Nwm+d*NMact+d*NMuscle_FLV+nqi_q+1:Nwl:Nw_acc+Nw);   
            q_dot_opt_ext(mpi).mp(3:(d+1):end,nqi) = w_opt(Nw_acc+Nwm+d*NMact+d*NMuscle_FLV+2*nq.all+nqi_q+1:Nwl:Nw_acc+Nw);  
            q_dot_opt_ext(mpi).mp(4:(d+1):end,nqi) = w_opt(Nw_acc+Nwm+d*NMact+d*NMuscle_FLV+2*nq.all+2*nq.all+nqi_q+1:Nwl:Nw_acc+Nw);
        end
        Nw_acc = Nw_acc + Nw - NParameters;
    end
    
    %% Unscale results
    % Parameters
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
        q_opt_unsc(mpi).mp.rad = q_opt(mpi).mp(1:end-1,:).*repmat(scaling(mpi).mp.Qs,size(q_opt(mpi).mp(1:end-1,:),1),1); 
        % Convert in degrees
        q_opt_unsc(mpi).mp.deg = q_opt_unsc(mpi).mp.rad;
        q_opt_unsc(mpi).mp.deg(:,[1:3,7:end]) = q_opt_unsc(mpi).mp.deg(:,[1:3,7:end]).*180/pi;
        % Qs (1:N)
        q_opt_unsc_all(mpi).mp.rad = q_opt(mpi).mp(1:end,:).*repmat(scaling(mpi).mp.Qs,size(q_opt(mpi).mp(1:end,:),1),1); 
        % Convert in degrees
        q_opt_unsc_all(mpi).mp.deg = q_opt_unsc_all(mpi).mp.rad;
        q_opt_unsc_all(mpi).mp.deg(:,[1:3,7:end]) = q_opt_unsc_all(mpi).mp.deg(:,[1:3,7:end]).*180/pi;
        % Qdots (1:N-1)
        qdot_opt_unsc(mpi).mp.rad = qdot_opt(mpi).mp(1:end-1,:).*repmat(scaling(mpi).mp.Qdots,size(qdot_opt(mpi).mp(1:end-1,:),1),1);
        % Convert in degrees
        qdot_opt_unsc(mpi).mp.deg = qdot_opt_unsc(mpi).mp.rad;
        qdot_opt_unsc(mpi).mp.deg(:,[1:3,7:end]) = qdot_opt_unsc(mpi).mp.deg(:,[1:3,7:end]).*180/pi;
        % Qdots (1:N)
        qdot_opt_unsc_all(mpi).mp.rad = qdot_opt(mpi).mp(1:end,:).*repmat(scaling(mpi).mp.Qdots,size(qdot_opt(mpi).mp(1:end,:),1),1); 
        % Muscle activations
        a_syn_opt(mpi).mp = a_opt(mpi).mp(1:end-1,:).*repmat(scaling(mpi).mp.a,size(a_opt(mpi).mp(1:end-1,:),1),size(a_opt(mpi).mp,2));
        % Muscle-tendon forces
        FTtilde_opt_unsc(mpi).mp = FTtilde_opt(mpi).mp(1:end-1,:).*repmat(scaling(mpi).mp.FTtilde,size(FTtilde_opt(mpi).mp(1:end-1,:),1),1);
        % Controls at mesh points
        % Time derivative of Qdots
        qdotdot_opt_unsc(mpi).mp.rad = qdotdot_opt(mpi).mp.*repmat(scaling(mpi).mp.Qdotdots,size(qdotdot_opt(mpi).mp,1),1);
        % Back activations
        a_b_opt_unsc(mpi).mp = a_b_opt(mpi).mp(1:end-1,:).*repmat(scaling(mpi).mp.a_b,size(a_b_opt(mpi).mp(1:end-1,:),1),size(a_b_opt(mpi).mp,2));
        e_b_opt_unsc(mpi).mp = e_b_opt(mpi).mp.*repmat(scaling(mpi).mp.e_b,size(e_b_opt(mpi).mp,1),size(e_b_opt(mpi).mp,2));
        % Convert in degrees
        qdotdot_opt_unsc(mpi).mp.deg = qdotdot_opt_unsc(mpi).mp.rad;
        qdotdot_opt_unsc(mpi).mp.deg(:,[1:3,7:end]) = qdotdot_opt_unsc(mpi).mp.deg(:,[1:3,7:end]).*180/pi;
        % Time derivative of muscle activations (states)
        vA_opt_unsc(mpi).mp = vA_opt(mpi).mp.*repmat(scaling(mpi).mp.vA,size(vA_opt(mpi).mp,1),size(vA_opt(mpi).mp,2));
        tact = 0.015;
        tdeact = 0.06;
        % Get muscle excitations from time derivative of muscle activations
        % Time derivative of muscle-tendon forces
        dFTtilde_opt_unsc(mpi).mp = dFTtilde_opt(mpi).mp.*repmat(scaling(mpi).mp.dFTtilde,size(dFTtilde_opt(mpi).mp,1),size(dFTtilde_opt(mpi).mp,2));
        % Reconstruct synergies
        if W.Syn == 0
            a_opt_unsc(mpi).mp = a_syn_opt(mpi).mp;
        else
            a_syn_opt_l = zeros(N,NMuscle/2);
            a_syn_opt_r = zeros(N,NMuscle/2);
            for n = 1:N
                a_syn_opt_l(n,:) = full(f_SynergyProduct(syn_wl_opt,a_syn_opt(mpi).mp(n,1:NSyn)))';
                a_syn_opt_r(n,:) = full(f_SynergyProduct(syn_wr_opt,a_syn_opt(mpi).mp(n,NSyn+1:2*NSyn)))';  
            end
            a_opt_unsc(mpi).mp = [a_syn_opt_l,a_syn_opt_r];   
        end
        e_opt_unsc(mpi).mp = computeExcitationRaasch(a_syn_opt(mpi).mp,vA_opt_unsc(mpi).mp,ones(1,NMact)*tdeact,ones(1,NMact)*tact);
        % Combination of supra-spinal inputs and reflex feedback
        a_tot_opt(mpi).mp = a_opt_unsc(mpi).mp;         
        if spasi == 1
            % Muscle activations from muscle-tendon force feedback and time derivative
            % of muscle-tendon force feedback
            a_Ff_opt_unsc(mpi).mp = a_Ff_opt(mpi).mp(1:end-1,:);
            a_dFf_opt_unsc(mpi).mp = a_dFf_opt(mpi).mp(1:end-1,:);        
            a_tot_opt(mpi).mp(:,musi_Spas) = a_tot_opt(mpi).mp(:,musi_Spas) + a_Ff_opt_unsc(mpi).mp + a_dFf_opt_unsc(mpi).mp; 
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
    for mpi = 1:NPhases
        Xk_Qs_Qdots_opt(mpi).mp             = zeros(N,2*size(q_opt_unsc(mpi).mp.rad,2));
        Xk_Qs_Qdots_opt(mpi).mp(:,1:2:end)  = q_opt_unsc(mpi).mp.rad;
        Xk_Qs_Qdots_opt(mpi).mp(:,2:2:end)  = qdot_opt_unsc(mpi).mp.rad;
        Xk_Qdotdots_opt(mpi).mp             = qdotdot_opt_unsc(mpi).mp.rad;  
        for i = 1:N
            out2_res = F([Xk_Qs_Qdots_opt(mpi).mp(i,:)';Xk_Qdotdots_opt(mpi).mp(i,:)']);        
            out_res_opt(mpi).mp(i,:) = full(out2_res);  
        end
        % Optimal joint torques, ground reaction forces and moments
        Tauk_out(mpi).mp        = out_res_opt(mpi).mp(:,residualsi);
        GRF_opt_unsc(mpi).mp    = out_res_opt(mpi).mp(:,GRFi.all);
        GRM_opt_unsc(mpi).mp    = out_res_opt(mpi).mp(:,GRMi.all);
        % Stride length (transversal plane)        
        dist = sqrt(sumsqr(out_res_opt(mpi).mp(end,calcOr.r)-...
            out_res_opt(mpi).mp(1,calcOr.r)));
        StrideLength_opt(mpi).mp = full(dist);
        % Step width 
        StepWidth_opt = full(abs(out_res_opt(mpi).mp(:,calcOr.r(2)) - ...
            out_res_opt(mpi).mp(:,calcOr.l(2))));
        StepWidth_opt_mean(mpi).mp = mean(StepWidth_opt);
        StepWidth_opt_std = std(StepWidth_opt);
        % assertArmTmax should be 0
        assertTrunkTmax = max(max(abs(out_res_opt(mpi).mp(:,backi)-...
            (a_b_opt_unsc(mpi).mp)*scaling(mpi).mp.TrunkTau)));  
        if assertTrunkTmax > 1*10^(-tol_ipopt)
            disp('Issue when reconstructing residual forces')
        end
    end
       
    %% Passive joint torques at optimal solution
    Tau_pass_opt_all = struct('mp',[]);
    Tau_pass_opt = struct('mp',[]);
    for mpi = 1:NPhases
        Tau_pass_opt_all(mpi).mp = zeros(N,nq.act);
        for i = 1:N    
            Tau_pass_opt(mpi).mp.hip.flex.l = f_PassiveTorques(k_pass.hip.flex,...
                theta.pass.hip.flex,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.l*2));
            Tau_pass_opt(mpi).mp.hip.flex.r = f_PassiveTorques(k_pass.hip.flex,...
                theta.pass.hip.flex,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_flex.r*2));
            Tau_pass_opt(mpi).mp.hip.add.l = f_PassiveTorques(k_pass.hip.add,...
                theta.pass.hip.add,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.l*2));
            Tau_pass_opt(mpi).mp.hip.add.r = f_PassiveTorques(k_pass.hip.add,...
                theta.pass.hip.add,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_add.r*2));
            Tau_pass_opt(mpi).mp.hip.rot.l = f_PassiveTorques(k_pass.hip.rot,...
                theta.pass.hip.rot,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.l*2));
            Tau_pass_opt(mpi).mp.hip.rot.r = f_PassiveTorques(k_pass.hip.rot,...
                theta.pass.hip.rot,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.hip_rot.r*2));
            Tau_pass_opt(mpi).mp.knee.l = f_PassiveTorques(k_pass.knee,...
                theta.pass.knee,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.l*2));
            Tau_pass_opt(mpi).mp.knee.r = f_PassiveTorques(k_pass.knee,...
                theta.pass.knee,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.knee.r*2));
            Tau_pass_opt(mpi).mp.ankle.l = f_PassiveTorques(k_pass.ankle,...
                theta.pass.ankle,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.l*2));
            Tau_pass_opt(mpi).mp.ankle.r = f_PassiveTorques(k_pass.ankle,...
                theta.pass.ankle,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.ankle.r*2));
            Tau_pass_opt(mpi).mp.subt.l = f_PassiveTorques(k_pass.subt,...
                theta.pass.subt,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.l*2));
            Tau_pass_opt(mpi).mp.subt.r = f_PassiveTorques(k_pass.subt,...
                theta.pass.subt,Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.subt.r*2));   
            Tau_pass_opt(mpi).mp.trunk.ext = f_PassiveTorques(...
                k_pass.trunk.ext,theta.pass.trunk.ext,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.ext*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.ext*2));
            Tau_pass_opt(mpi).mp.trunk.ben = f_PassiveTorques(...
                k_pass.trunk.ben,theta.pass.trunk.ben,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.ben*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.ben*2));
            Tau_pass_opt(mpi).mp.trunk.rot = f_PassiveTorques(...
                k_pass.trunk.rot,theta.pass.trunk.rot,...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.rot*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(i,jointi.trunk.rot*2)); 
            Tau_pass_opt_all(mpi).mp(i,:) = full([...
                Tau_pass_opt(mpi).mp.hip.flex.l,Tau_pass_opt(mpi).mp.hip.add.l,...
                Tau_pass_opt(mpi).mp.hip.rot.l,Tau_pass_opt(mpi).mp.hip.flex.r,...
                Tau_pass_opt(mpi).mp.hip.add.r,Tau_pass_opt(mpi).mp.hip.rot.r,...
                Tau_pass_opt(mpi).mp.knee.l,Tau_pass_opt(mpi).mp.knee.r,...
                Tau_pass_opt(mpi).mp.ankle.l,Tau_pass_opt(mpi).mp.ankle.r,...
                Tau_pass_opt(mpi).mp.subt.l,Tau_pass_opt(mpi).mp.subt.r,...
                Tau_pass_opt(mpi).mp.trunk.ext,Tau_pass_opt(mpi).mp.trunk.ben,...
                Tau_pass_opt(mpi).mp.trunk.rot]); 
        end        
    end
    
    %% Reconstruct cost function and compute metabolic cost
if decomposeCost
    J_opt           = 0;
    Qs_cost         = 0;
    GRF_cost        = 0;
    GRM_cost        = 0;
    ID_cost         = 0;
    a_cost          = 0;
    a_cost_temp     = 0;
    Qdotdots_cost   = 0;
    vA_cost         = 0;
    dFTtilde_cost   = 0;
    syn_cost        = 0;
    mE_cost         = 0;
    passT_cost      = 0;
    trunk_cost      = 0;
    a_effort        = 0;
    a_fatigue       = 0;
    e_mo_opt        = struct('mp',[]);
    e_mo_opt_tr     = struct('mp',[]);
    COT_opt         = struct('mp',[]);
    for mpi = 1:NPhases        
        count           = 1;
        h_opt	= tf_opt/N;  
        % Get COT
        dist_trav_opt = Xk_Qs_Qdots_opt(mpi).mp(end,2*jointi.pelvis.tx-1) - ...
        Xk_Qs_Qdots_opt(mpi).mp(1,2*jointi.pelvis.tx-1); % distance traveled 
        speed_opt = dist_trav_opt/tf_opt;
        assertSpeed = abs(speed_opt - v_tgt);
        if assertSpeed > 1e-10
            disp(['Issue when reconstructing speed: difference is ',...
                num2str(assertSpeed)])
        end         
        for k=1:N               
            % Get muscle-tendon lengths, velocities, moment arms
            % Left leg
            qin_l_opt_all = [Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_flex.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_add.l*2-1), ...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_rot.l*2-1), ...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.knee.l*2-1), ...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.ankle.l*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.subt.l*2-1)];  
            qdotin_l_opt_all = [Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_flex.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_add.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_rot.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.knee.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.ankle.l*2),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.subt.l*2)];  
            [lMTk_l_opt_all,vMTk_l_opt_all,~] = ...
                f_lMT_vMT_dM_l(qin_l_opt_all,qdotin_l_opt_all);    
            % Right leg
            qin_r_opt_all = [Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_flex.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_add.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_rot.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.knee.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.ankle.r*2-1),...
                Xk_Qs_Qdots_opt(mpi).mp(k,jointi.subt.r*2-1)];  
            qdotin_r_opt_all = [Xk_Qs_Qdots_opt(mpi).mp(k,jointi.hip_flex.r*2),...
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
                    a_tot_opt(mpi).mp(k,musi_FLV)',FTtilde_opt_unsc(mpi).mp(k,:)',...
                    dFTtilde_opt_unsc(mpi).mp(k,:)',full(lMTk_lr_opt_all(musi_FLV)),...
                    full(vMTk_lr_opt_all(musi_FLV)),tensions(musi_FLV));                  
            [~,lMtilde_opt_all] = f_FiberLength_TendonForce(...
                FTtilde_opt_unsc(mpi).mp(k,:)',full(lMTk_lr_opt_all(musi_FLV)));                
            [vM_opt_all,~] = ...
                f_FiberVelocity_TendonForce(FTtilde_opt_unsc(mpi).mp(k,:)',...
                dFTtilde_opt_unsc(mpi).mp(k,:)',full(lMTk_lr_opt_all(musi_FLV)),...
                full(vMTk_lr_opt_all(musi_FLV)));
            [e_tot_all,~,~,~,~,e_mot] = fgetMetabolicEnergySmooth2004all(...
            a_tot_opt(mpi).mp(k,musi_FLV)',a_tot_opt(mpi).mp(k,musi_FLV)',...
            full(lMtilde_opt_all),full(vM_opt_all),full(Fce_opt_all)',...
            full(Fpass_opt_all)',full(massM_opt_all)',pctsts(musi_FLV),...
            full(Fiso_opt_all)',MTparameters_m_FLV(1,:)',body_mass,10);                 
            e_tot_opt_all = full(e_tot_all)';
            e_mo_opt(mpi).mp(k,:) = full(e_mot)'; 
            % Muscle effort and fatigue
            a_effort = a_effort + (f_JNM_exp(a_opt_unsc(mpi).mp(k,:),2));
            a_fatigue = a_fatigue + (f_JNM_exp(a_opt_unsc(mpi).mp(k,:),10));
            for j=1:d 
                % Tracking terms
                Qs_cost_optk = W.Qs*B(j+1)*(f_JNq_bpty(qqdot_opt(mpi).mp(k,Qsi(residual_bptyi))-...
                    Qs(mpi).mp.allinterpfilt(k,residual_bptyi+1)... 
                    ./scaling(mpi).mp.Qs(residual_bptyi)))*h_opt;
                Qs_cost = Qs_cost + Qs_cost_optk;
                GRF_cost_optk = W.GRF*B(j+1)*(f_J6((out_res_opt(mpi).mp(k,GRFi.all)./scaling(mpi).mp.GRF)-...
                    GRF(mpi).mp.val.allinterp(k,2:end)./scaling(mpi).mp.GRF))*h_opt;
                GRF_cost = GRF_cost + GRF_cost_optk;
                GRM_cost_optk = W.GRM*B(j+1)*(f_J6((out_res_opt(mpi).mp(k,GRMi.all)./scaling(mpi).mp.GRM)-...
                    GRF(mpi).mp.MorGF.allinterp(k,2:end)./scaling(mpi).mp.GRM))*h_opt;
                GRM_cost = GRM_cost + GRM_cost_optk;
                ID_cost_optk = W.ID_act*B(j+1)*(f_JNq_act((out_res_opt(mpi).mp(k,residuals_acti)./scaling(mpi).mp.T(1))-...
                    ID(mpi).mp.allinterp(k,2+nq.abs:end)./scaling(mpi).mp.T(1)))*h_opt;
                ID_cost = ID_cost + ID_cost_optk;
                
                track_terms_optk = Qs_cost_optk + GRF_cost_optk +...
                    GRM_cost_optk + ID_cost_optk;
                % Motor control terms
                if acti == 0
                    a_cost_optk = W.a*B(j+1)*(f_JNM_exp(a_tot_opt(mpi).mp(k,:),exp_A))*h_opt;
                elseif acti == 1                
                    a_cost_optk = W.a*B(j+1)*(f_JNM_exp(a_opt_unsc(mpi).mp(k,:),exp_A))*h_opt;
                    a_cost_optk_temp = W.a*B(j+1)*(f_JNM_exp(a_opt_unsc(mpi).mp(k,:),1))*h_opt;
                end
                if W.E == 0
                    mE_cost_optk = 0;
                else
                    mE_cost_optk = W.E*B(j+1)*(f_JNM_FLVexp(e_tot_opt_all,exp_E))/body_mass*h_opt;
                end
                trunk_cost_optk = W.Trunk*B(j+1)*(f_Jnq_trunk(e_b_opt_unsc(mpi).mp(k,:)))*h_opt; 
                a_cost = a_cost + a_cost_optk;
                a_cost_temp = a_cost_temp + a_cost_optk_temp;
                mE_cost = mE_cost + mE_cost_optk;
                Qdotdots_cost_optk = W.qdotdot*B(j+1)*(f_JNq_all(qdotdot_opt(mpi).mp(k,:)))*h_opt;
                Qdotdots_cost = Qdotdots_cost + Qdotdots_cost_optk;
                vA_cost_optk = W.u*B(j+1)*(f_JNMact(vA_opt(mpi).mp(k,:)))*h_opt;
                vA_cost = vA_cost + vA_cost_optk;
                dFTtilde_cost_optk = W.u*B(j+1)*(f_JNM_FLV(dFTtilde_opt(mpi).mp(k,:)))*h_opt;
                dFTtilde_cost = dFTtilde_cost + dFTtilde_cost_optk;
                trunk_cost = trunk_cost + trunk_cost_optk;
                passT_cost_optk = W.passT*B(j+1)*(f_JNq_act(Tau_pass_opt_all(mpi).mp(k,:)))*h_opt;
                passT_cost = passT_cost + passT_cost_optk;                
                mc_terms_optk = a_cost_optk + mE_cost_optk +...
                    Qdotdots_cost_optk + passT_cost_optk + trunk_cost_optk + ...
                    vA_cost_optk + dFTtilde_cost_optk; 
                dist_norm_opt = dist_trav_opt;
                J_opt = J_opt + track_terms_optk + ...
                    1/dist_norm_opt*mc_terms_optk;  
                count = count + 1;                 
            end
        end            
        e_mo_opt_tr(mpi).mp = trapz(tgrid(mpi).mp(1:end-1),e_mo_opt(mpi).mp);
        % Energy model from Bhargava et al. (2004)
        COT_opt(mpi).mp = e_mo_opt_tr(mpi).mp/body_mass/dist_trav_opt; 
    end
    J_optf = full(J_opt);     
    Qs_costf = full(Qs_cost);
    GRF_costf = full(GRF_cost);
    GRM_costf = full(GRM_cost);
    ID_costf = full(ID_cost);
    Act_costf = full(a_cost);
    Act_costf_temp = full(a_cost_temp);
    Qdotdots_costf = full(Qdotdots_cost);
    vA_costf = full(vA_cost);
    dFTtilde_costf = full(dFTtilde_cost);
    syn_costf = full(syn_cost);
    trunk_costf = full(trunk_cost);
    passT_costf = full(passT_cost);
    mE_costf = full(mE_cost);
    % assertCost should be 0     
    assertCost = abs(J_optf - stats.iterations.obj(end));
    if assertCost > 1e-10
        disp(['Issue when reconstructing optimal cost: difference is ',...
            num2str(assertCost)])
    end 
    cont_ob = [Qs_costf,GRF_costf,GRM_costf,ID_costf,...
        1/dist_norm_opt*Act_costf,...
        1/dist_norm_opt*Qdotdots_costf,1/dist_norm_opt*vA_costf,...
        1/dist_norm_opt*dFTtilde_costf,1/dist_norm_opt*syn_costf,...
        1/dist_norm_opt*mE_costf,1/dist_norm_opt*trunk_costf,...
        1/dist_norm_opt*passT_costf]./J_optf*100;
    cont_ob_abs = [Qs_costf,GRF_costf,GRM_costf,ID_costf,...
        1/dist_norm_opt*Act_costf,...
        1/dist_norm_opt*Qdotdots_costf,1/dist_norm_opt*vA_costf,...
        1/dist_norm_opt*dFTtilde_costf,1/dist_norm_opt*syn_costf,...
        1/dist_norm_opt*mE_costf,1/dist_norm_opt*trunk_costf,...
        1/dist_norm_opt*passT_costf];
    cont_ob_abs_notNorm = [Qs_costf/W.Qs,GRF_costf/W.GRF,GRM_costf/W.GRM,...
        ID_costf/W.ID_act,...
        1/dist_norm_opt*Act_costf_temp/W.a,...
        1/dist_norm_opt*Qdotdots_costf/W.qdotdot,...
        1/dist_norm_opt*vA_costf/W.u,...
        1/dist_norm_opt*dFTtilde_costf/W.u,...
        1/dist_norm_opt*syn_costf/W.Syn,...
        1/dist_norm_opt*mE_costf/W.E,...
        1/dist_norm_opt*trunk_costf/W.Trunk,...
        1/dist_norm_opt*passT_costf/W.passT,full(a_effort),full(a_fatigue)];
end

    %% Extract passive forces and assert Hill-equilibrium   
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
            lMTk_opt_lr(mpi).mp(i,:)     = [full(lMTk_opt_l);full(lMTk_opt_r)];
            vMTk_opt_lr(mpi).mp(i,:)     = [full(vMTk_opt_l);full(vMTk_opt_r)];   
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
        assertHill = max(max(Hilldiffk_opt(mpi).mp));
        if assertHill > 1*10^(-tol_ipopt)
            disp('Issue in Hill-equilibrium')
        end
    end
    
    %% Reconstruct gait cycle: starting with right heel strike
    % Identify heel strike
    threshold = 30;
    if ww == 729
        threshold = 28;
    end
    if ww == 757
        threshold = 22;
    end
    if ww == 762
        threshold = 21;
    end
    IC1i = struct('mp',[]);
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
    
        Qs_opt_GC(mpi).mp       = zeros(N,size(q_opt_unsc(mpi).mp.deg,2));    
        Qs_opt_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:)  = q_opt_unsc(mpi).mp.deg(IC1i(mpi).mp:end,:);
        Qs_opt_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:)  = q_opt_unsc(mpi).mp.deg(1:IC1i(mpi).mp-1,:);           
        Qs_opt_GC(mpi).mp(N-IC1i(mpi).mp+2:N,jointi.pelvis.tx) = ...
            q_opt_unsc(mpi).mp.deg(1:IC1i(mpi).mp-1,jointi.pelvis.tx) + ...
            (q_opt_unsc_all(mpi).mp.deg(end,jointi.pelvis.tx)-...
            q_opt_unsc_all(mpi).mp.deg(1,jointi.pelvis.tx));    
        temp_q_opt_GC_pelvis_tx = Qs_opt_GC(mpi).mp(1,jointi.pelvis.tx);
        Qs_opt_GC(mpi).mp(:,jointi.pelvis.tx) = Qs_opt_GC(mpi).mp(:,jointi.pelvis.tx)-...
            temp_q_opt_GC_pelvis_tx;    
        % Visualization in OpenSim GUI    
        % Create .mot file for OpenSim GUI
        if writeIKmotion 
            pathOpenSim = [pathRepo,'\OpenSim'];
            addpath(genpath(pathOpenSim));
%             q_opt_GUI = zeros(N+1,1+nq.all);
%             if trunki ~= 1
%                 q_opt_GUI = zeros(N+1,1+nq.all+3);
%             end
%             q_opt_GUI(:,1) = tgrid(mpi).mp';
%             q_opt_GUI(:,2:nq.all+1)  = q_opt_unsc_all(mpi).mp.deg;
%             JointAngle.labels = {'time','lower_torso_RX','lower_torso_RY',...
%                 'lower_torso_RZ','lower_torso_TX','lower_torso_TY',...
%                 'lower_torso_TZ',...
%                 'hip_flex_l','hip_add_l','hip_rot_l',...
%                 'hip_flex_r','hip_add_r','hip_rot_r',...
%                 'knee_flex_l','knee_flex_r','ankle_flex_l',...
%                 'ankle_flex_r','subt_angle_l','subt_angle_r',...
%                 'lumbar_pitch','lumbar_roll','lumbar_yaw'};
%             JointAngle.data = q_opt_GUI;
%             filenameJointAngles = [pathRepo,'\Results\',namescript,...
%                     '\IK',savename_ind(mpi).mp,'.mot'];
%             write_motionFile(JointAngle, filenameJointAngles)            
            q_opt_GUI_HS = zeros(N,1+nq.all);
            q_opt_GUI_HS = zeros(N+1,1+nq.all+3);
            q_opt_GUI_HS(:,1) = tgrid(mpi).mp(1:end-1)';
            q_opt_GUI_HS(:,2:nq.all+1)  = Qs_opt_GC(mpi).mp;
            JointAngle_HS.data = q_opt_GUI_HS;
            JointAngle_HS.labels = {'time','lower_torso_RX','lower_torso_RY',...
                'lower_torso_RZ','lower_torso_TX','lower_torso_TY',...
                'lower_torso_TZ',...
                'hip_flex_l','hip_add_l','hip_rot_l',...
                'hip_flex_r','hip_add_r','hip_rot_r',...
                'knee_flex_l','knee_flex_r','ankle_flex_l',...
                'ankle_flex_r','subt_angle_l','subt_angle_r',...
                'lumbar_pitch','lumbar_roll','lumbar_yaw'};
            filenameJointAngles = [pathRepo,'\Results\',namescript,...
                    '\IK',savename_ind(mpi).mp,'_HS.mot'];
            write_motionFile(JointAngle_HS, filenameJointAngles)
        end
        % Ground reaction forces
        GRF_opt_GC(mpi).mp = zeros(N,nGRF);
        GRF_opt_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = GRF_opt_unsc(mpi).mp(IC1i(mpi).mp:end,:);
        GRF_opt_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = GRF_opt_unsc(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Ground reaction torques
        GRM_opt_GC(mpi).mp = zeros(N,nGRF);
        GRM_opt_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = GRM_opt_unsc(mpi).mp(IC1i(mpi).mp:end,:);
        GRM_opt_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = GRM_opt_unsc(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Joint torques
        Tauk_out_GC(mpi).mp = zeros(N,nq.all);    
        Tauk_out_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = Tauk_out(mpi).mp(IC1i(mpi).mp:end,:);    
        Tauk_out_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = Tauk_out(mpi).mp(1:IC1i(mpi).mp-1,:);

        % Muscle activations
        a_opt_unsc_GC(mpi).mp = zeros(N,NMuscle);    
        a_opt_unsc_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = a_opt_unsc(mpi).mp(IC1i(mpi).mp:end,:);    
        a_opt_unsc_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = a_opt_unsc(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Muscle activations
        a_tot_opt_GC(mpi).mp = zeros(N,NMuscle);    
        a_tot_opt_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = a_tot_opt(mpi).mp(IC1i(mpi).mp:end,:);    
        a_tot_opt_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = a_tot_opt(mpi).mp(1:IC1i(mpi).mp-1,:);
        % Muscle activations
        if spasi == 1
            a_Ff_opt_unsc_GC(mpi).mp = zeros(N,NMuscle_Spas);    
            a_Ff_opt_unsc_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = a_Ff_opt_unsc(mpi).mp(IC1i(mpi).mp:end,:);    
            a_Ff_opt_unsc_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = a_Ff_opt_unsc(mpi).mp(1:IC1i(mpi).mp-1,:);
            % Muscle activations
            a_dFf_opt_unsc_GC(mpi).mp = zeros(N,NMuscle_Spas);    
            a_dFf_opt_unsc_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = a_dFf_opt_unsc(mpi).mp(IC1i(mpi).mp:end,:);    
            a_dFf_opt_unsc_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = a_dFf_opt_unsc(mpi).mp(1:IC1i(mpi).mp-1,:);
        end
        % Passive forces
        Fpe_opt_GC(mpi).mp = zeros(N,NMuscle);    
        Fpe_opt_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = Fpe_opt(mpi).mp(IC1i(mpi).mp:end,:);    
        Fpe_opt_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = Fpe_opt(mpi).mp(1:IC1i(mpi).mp-1,:);

        % Fiber lengths
        lMtilde_opt_GC(mpi).mp = zeros(N,NMuscle);    
        lMtilde_opt_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = lMtilde_opt(mpi).mp(IC1i(mpi).mp:end,:);    
        lMtilde_opt_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = lMtilde_opt(mpi).mp(1:IC1i(mpi).mp-1,:);

        % muscle-tendon lengths
        lMTk_opt_lr_GC(mpi).mp = zeros(N,NMuscle);    
        lMTk_opt_lr_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = lMTk_opt_lr(mpi).mp(IC1i(mpi).mp:end,:);    
        lMTk_opt_lr_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = lMTk_opt_lr(mpi).mp(1:IC1i(mpi).mp-1,:);

        % muscle-tendon velocities
        vMTk_opt_lr_GC(mpi).mp = zeros(N,NMuscle);    
        vMTk_opt_lr_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = vMTk_opt_lr(mpi).mp(IC1i(mpi).mp:end,:);    
        vMTk_opt_lr_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = vMTk_opt_lr(mpi).mp(1:IC1i(mpi).mp-1,:);

        % synergy activations
        if W.Syn ~= 0
            a_syn_opt_GC(mpi).mp = zeros(N,NMact);    
            a_syn_opt_GC(mpi).mp(1:N-IC1i(mpi).mp+1,:) = a_syn_opt(mpi).mp(IC1i(mpi).mp:end,:);    
            a_syn_opt_GC(mpi).mp(N-IC1i(mpi).mp+2:N,:) = a_syn_opt(mpi).mp(1:IC1i(mpi).mp-1,:);
        end
    end
    
    %% Reconstruct gait cycle: starting with left heel strike
    % Identify heel strike
    threshold = 30;
    IC1i_l = struct('mp',[]);
    for mpi = 1:NPhases
        if exist('HS1_l','var')
            clear HS1_l
        end
        % Right heel strike first    
        phase_tran_tgridi_l = find(GRF_opt_unsc(mpi).mp(:,5)<threshold,1,'last');
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
        % Joint kinematics
        % Qs   
        Qs_opt_GC_l(mpi).mp = zeros(N,size(q_opt_unsc(mpi).mp.deg,2));    
        Qs_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:)  = q_opt_unsc(mpi).mp.deg(IC1i_l(mpi).mp:end,:);
        Qs_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:)  = q_opt_unsc(mpi).mp.deg(1:IC1i_l(mpi).mp-1,:);           
        Qs_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,jointi.pelvis.tx) = ...
            q_opt_unsc(mpi).mp.deg(1:IC1i_l(mpi).mp-1,jointi.pelvis.tx) + ...
            (q_opt_unsc_all(mpi).mp.deg(end,jointi.pelvis.tx)-...
            q_opt_unsc_all(mpi).mp.deg(1,jointi.pelvis.tx));    
        temp_q_opt_GC_pelvis_tx = Qs_opt_GC_l(mpi).mp(1,jointi.pelvis.tx);
        Qs_opt_GC_l(mpi).mp(:,jointi.pelvis.tx) = Qs_opt_GC_l(mpi).mp(:,jointi.pelvis.tx)-...
            temp_q_opt_GC_pelvis_tx;    
        % Visualization in OpenSim GUI    
        % Create .mot file for OpenSim GUI
        if writeIKmotion_l 
            pathOpenSim = [pathRepo,'\OpenSim'];
            addpath(genpath(pathOpenSim));        
            q_opt_GUI_HS_l = zeros(N,1+nq.all);
            q_opt_GUI_HS_l(:,1) = tgrid(mpi).mp(1:end-1)';
            q_opt_GUI_HS_l(:,2:nq.all+1)  = Qs_opt_GC_l(mpi).mp;
            JointAngle_HS.data = q_opt_GUI_HS_l;
            JointAngle_HS.labels = {'time','lower_torso_RX','lower_torso_RY',...
                'lower_torso_RZ','lower_torso_TX','lower_torso_TY',...
                'lower_torso_TZ',...
                'hip_flex_l','hip_add_l','hip_rot_l',...
                'hip_flex_r','hip_add_r','hip_rot_r',...
                'knee_flex_l','knee_flex_r','ankle_flex_l',...
                'ankle_flex_r','subt_angle_l','subt_angle_r',...
                'lumbar_pitch','lumbar_roll','lumbar_yaw'};
            filenameJointAngles = [pathRepo,'\Results\',namescript,...
                    '\IK',savename_ind(mpi).mp,'_HS_l.mot'];
            write_motionFile(JointAngle_HS, filenameJointAngles)
        end
        % Ground reaction forces
        GRF_opt_GC_l(mpi).mp = zeros(N,nGRF);
        GRF_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = GRF_opt_unsc(mpi).mp(IC1i_l(mpi).mp:end,:);
        GRF_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = GRF_opt_unsc(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Ground reaction torques
        GRM_opt_GC_l(mpi).mp = zeros(N,nGRF);
        GRM_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = GRM_opt_unsc(mpi).mp(IC1i_l(mpi).mp:end,:);
        GRM_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = GRM_opt_unsc(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Joint torques
        Tauk_out_GC_l(mpi).mp = zeros(N,nq.all);    
        Tauk_out_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = Tauk_out(mpi).mp(IC1i_l(mpi).mp:end,:);    
        Tauk_out_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = Tauk_out(mpi).mp(1:IC1i_l(mpi).mp-1,:);

        % Muscle activations
        a_opt_unsc_GC_l(mpi).mp = zeros(N,NMuscle);    
        a_opt_unsc_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = a_opt_unsc(mpi).mp(IC1i_l(mpi).mp:end,:);    
        a_opt_unsc_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = a_opt_unsc(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Muscle activations
        a_tot_opt_GC_l(mpi).mp = zeros(N,NMuscle);    
        a_tot_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = a_tot_opt(mpi).mp(IC1i_l(mpi).mp:end,:);    
        a_tot_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = a_tot_opt(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        % Muscle activations
        if spasi == 1
            a_Ff_opt_unsc_GC_l(mpi).mp = zeros(N,NMuscle_Spas);    
            a_Ff_opt_unsc_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = a_Ff_opt_unsc(mpi).mp(IC1i_l(mpi).mp:end,:);    
            a_Ff_opt_unsc_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = a_Ff_opt_unsc(mpi).mp(1:IC1i_l(mpi).mp-1,:);
            % Muscle activations
            a_dFf_opt_unsc_GC_l(mpi).mp = zeros(N,NMuscle_Spas);    
            a_dFf_opt_unsc_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = a_dFf_opt_unsc(mpi).mp(IC1i_l(mpi).mp:end,:);    
            a_dFf_opt_unsc_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = a_dFf_opt_unsc(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        end
        % Passive forces
        Fpe_opt_GC_l(mpi).mp = zeros(N,NMuscle);    
        Fpe_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = Fpe_opt(mpi).mp(IC1i_l(mpi).mp:end,:);    
        Fpe_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = Fpe_opt(mpi).mp(1:IC1i_l(mpi).mp-1,:);

        % Fiber lengths
        lMtilde_opt_GC_l(mpi).mp = zeros(N,NMuscle);    
        lMtilde_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = lMtilde_opt(mpi).mp(IC1i_l(mpi).mp:end,:);    
        lMtilde_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = lMtilde_opt(mpi).mp(1:IC1i_l(mpi).mp-1,:);

        % muscle-tendon lengths
        lMTk_opt_lr_GC_l(mpi).mp = zeros(N,NMuscle);    
        lMTk_opt_lr_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = lMTk_opt_lr(mpi).mp(IC1i_l(mpi).mp:end,:);    
        lMTk_opt_lr_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = lMTk_opt_lr(mpi).mp(1:IC1i_l(mpi).mp-1,:);

        % muscle-tendon velocities
        vMTk_opt_lr_GC_l(mpi).mp = zeros(N,NMuscle);    
        vMTk_opt_lr_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = vMTk_opt_lr(mpi).mp(IC1i_l(mpi).mp:end,:);    
        vMTk_opt_lr_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = vMTk_opt_lr(mpi).mp(1:IC1i_l(mpi).mp-1,:);

        % synergy activations
        if W.Syn ~= 0
            a_syn_opt_GC_l(mpi).mp = zeros(N,NMact);    
            a_syn_opt_GC_l(mpi).mp(1:N-IC1i_l(mpi).mp+1,:) = a_syn_opt(mpi).mp(IC1i_l(mpi).mp:end,:);    
            a_syn_opt_GC_l(mpi).mp(N-IC1i_l(mpi).mp+2:N,:) = a_syn_opt(mpi).mp(1:IC1i_l(mpi).mp-1,:);
        end
    end
    
        
    %% Save results
    for mpi = 1:NPhases
        if saveResults
            if (exist([pathresults,'\',namescript,'\Results_tracking.mat'],...
                'file')==2) 
            load([pathresults,'\',namescript,'\Results_tracking.mat']);
            else
                Results_tracking.(['Case_',num2str(ww)]) = struct('Qs_opt',[]);
            end
        % Structure results
        Results_tracking.(['Case_',num2str(ww)]).tgrid(mpi).mp = ...
            tgrid(mpi).mp(1:end-1)';
        Results_tracking.(['Case_',num2str(ww)]).Qs_opt(mpi).mp = ...
            q_opt_unsc(mpi).mp.deg;
        Results_tracking.(['Case_',num2str(ww)]).Qs_opt_GC(mpi).mp = ...
            Qs_opt_GC(mpi).mp;
        Results_tracking.(['Case_',num2str(ww)]).Qs_opt_GC_l(mpi).mp = ...
            Qs_opt_GC_l(mpi).mp;
        Results_tracking.(['Case_',num2str(ww)]).Acts_opt(mpi).mp = ...
            a_opt_unsc(mpi).mp;
        Results_tracking.(['Case_',num2str(ww)]).Acts_opt_GC(mpi).mp = ...
            a_opt_unsc_GC(mpi).mp;  
        Results_tracking.(['Case_',num2str(ww)]).Acts_opt_GC_l(mpi).mp = ...
            a_opt_unsc_GC_l(mpi).mp;  
        Results_tracking.(['Case_',num2str(ww)]).Acts_tot_opt(mpi).mp = ...
            a_tot_opt(mpi).mp;
        Results_tracking.(['Case_',num2str(ww)]).Acts_tot_opt_GC(mpi).mp = ...
            a_tot_opt_GC(mpi).mp;   
        Results_tracking.(['Case_',num2str(ww)]).Acts_tot_opt_GC_l(mpi).mp = ...
            a_tot_opt_GC_l(mpi).mp;   
        if spasi == 1
            Results_tracking.(['Case_',num2str(ww)]).a_Ff_opt_unsc(mpi).mp = ...
                a_Ff_opt_unsc(mpi).mp;
            Results_tracking.(['Case_',num2str(ww)]).a_Ff_opt_unsc_GC(mpi).mp = ...
                a_Ff_opt_unsc_GC(mpi).mp;
            Results_tracking.(['Case_',num2str(ww)]).a_Ff_opt_unsc_GC_l(mpi).mp = ...
                a_Ff_opt_unsc_GC_l(mpi).mp;
            Results_tracking.(['Case_',num2str(ww)]).a_dFf_opt_unsc(mpi).mp = ...
                a_dFf_opt_unsc(mpi).mp;
            Results_tracking.(['Case_',num2str(ww)]).a_dFf_opt_unsc_GC(mpi).mp = ...
                a_dFf_opt_unsc_GC(mpi).mp;
            Results_tracking.(['Case_',num2str(ww)]).a_dFf_opt_unsc_GC_l(mpi).mp = ...
                a_dFf_opt_unsc_GC_l(mpi).mp;
        else
            Results_tracking.(['Case_',num2str(ww)]).a_Ff_opt_unsc(mpi).mp = ...
            NaN;
            Results_tracking.(['Case_',num2str(ww)]).a_Ff_opt_unsc_GC(mpi).mp = ...
            NaN;
            Results_tracking.(['Case_',num2str(ww)]).a_Ff_opt_unsc_GC_l(mpi).mp = ...
            NaN;
            Results_tracking.(['Case_',num2str(ww)]).a_dFf_opt_unsc(mpi).mp =...
            NaN;
            Results_tracking.(['Case_',num2str(ww)]).a_dFf_opt_unsc_GC(mpi).mp =...
            NaN;
            Results_tracking.(['Case_',num2str(ww)]).a_dFf_opt_unsc_GC_l(mpi).mp =...
            NaN;
        end
        Results_tracking.(['Case_',num2str(ww)]).Ts_opt(mpi).mp = ...
            Tauk_out(mpi).mp;
        Results_tracking.(['Case_',num2str(ww)]).Ts_opt_GC(mpi).mp = ...
            Tauk_out_GC(mpi).mp;     
        Results_tracking.(['Case_',num2str(ww)]).Ts_opt_GC_l(mpi).mp = ...
            Tauk_out_GC_l(mpi).mp; 
        Results_tracking.(['Case_',num2str(ww)]).GRFs_opt(mpi).mp = ...
            GRF_opt_unsc(mpi).mp;
        Results_tracking.(['Case_',num2str(ww)]).GRFs_opt_GC(mpi).mp = ...
            GRF_opt_GC(mpi).mp;        
        Results_tracking.(['Case_',num2str(ww)]).GRFs_opt_GC_l(mpi).mp = ...
            GRF_opt_GC_l(mpi).mp;     
        Results_tracking.(['Case_',num2str(ww)]).GRMs_opt(mpi).mp = ...
            GRM_opt_unsc(mpi).mp; 
        Results_tracking.(['Case_',num2str(ww)]).GRMs_opt_GC(mpi).mp = ...
            GRM_opt_GC(mpi).mp;  
        Results_tracking.(['Case_',num2str(ww)]).GRMs_opt_GC_l(mpi).mp = ...
            GRM_opt_GC_l(mpi).mp; 
        Results_tracking.(['Case_',num2str(ww)]).Fpe_opt(mpi).mp = ...
            Fpe_opt(mpi).mp; 
        Results_tracking.(['Case_',num2str(ww)]).Fpe_opt_GC(mpi).mp = ...
            Fpe_opt_GC(mpi).mp;
        Results_tracking.(['Case_',num2str(ww)]).Fpe_opt_GC_l(mpi).mp = ...
            Fpe_opt_GC_l(mpi).mp;
        Results_tracking.(['Case_',num2str(ww)]).lMtilde_opt(mpi).mp = ...
            lMtilde_opt(mpi).mp; 
        Results_tracking.(['Case_',num2str(ww)]).lMtilde_opt_GC(mpi).mp = ...
            lMtilde_opt_GC(mpi).mp;    
        Results_tracking.(['Case_',num2str(ww)]).lMtilde_opt_GC_l(mpi).mp = ...
            lMtilde_opt_GC_l(mpi).mp; 
        Results_tracking.(['Case_',num2str(ww)]).lMTk_opt(mpi).mp = ...
            lMTk_opt_lr(mpi).mp; 
        Results_tracking.(['Case_',num2str(ww)]).lMTk_opt_GC(mpi).mp = ...
            lMTk_opt_lr_GC(mpi).mp;         
        Results_tracking.(['Case_',num2str(ww)]).lMTk_opt_GC_l(mpi).mp = ...
            lMTk_opt_lr_GC_l(mpi).mp; 
        Results_tracking.(['Case_',num2str(ww)]).vMTk_opt(mpi).mp = ...
            vMTk_opt_lr(mpi).mp; 
        Results_tracking.(['Case_',num2str(ww)]).vMTk_opt_GC(mpi).mp = ...
            vMTk_opt_lr_GC(mpi).mp;
        Results_tracking.(['Case_',num2str(ww)]).vMTk_opt_GC_l(mpi).mp = ...
            vMTk_opt_lr_GC_l(mpi).mp;
        
        if W.Syn ~= 0
            Results_tracking.(['Case_',num2str(ww)]).a_syn_opt(mpi).mp = ...
                a_syn_opt(mpi).mp;  
            Results_tracking.(['Case_',num2str(ww)]).a_syn_opt_GC(mpi).mp = ...
                a_syn_opt_GC(mpi).mp; 
            Results_tracking.(['Case_',num2str(ww)]).a_syn_opt_GC_l(mpi).mp = ...
                a_syn_opt_GC_l(mpi).mp; 
            Results_tracking.(['Case_',num2str(ww)]).w_syn_opt_GC(mpi).mp = ...
                [syn_wl_opt,syn_wr_opt];     
        end
        
        Results_tracking.(['Case_',num2str(ww)]).COT_opt(mpi).mp = ...
            COT_opt(mpi).mp;
        Results_tracking.(['Case_',num2str(ww)]).Qs_toTrack(mpi).mp = ...
            Qs(mpi).mp.allinterpfilt;
        Results_tracking.(['Case_',num2str(ww)]).Ts_toTrack(mpi).mp = ...
            ID(mpi).mp.allinterp;
        Results_tracking.(['Case_',num2str(ww)]).GRFs_toTrack(mpi).mp = ...
            GRF(mpi).mp.val.allinterp;
        Results_tracking.(['Case_',num2str(ww)]).GRMs_toTrack(mpi).mp = ...
            GRF(mpi).mp.MorGF.allinterp; 
        Results_tracking.(['Case_',num2str(ww)]).stats = stats;
        Results_tracking.(['Case_',num2str(ww)]).cont_ob = cont_ob;
        Results_tracking.(['Case_',num2str(ww)]).cont_ob_abs = cont_ob_abs;
        Results_tracking.(['Case_',num2str(ww)]).cont_ob_abs_notNorm = cont_ob_abs_notNorm;
        Results_tracking.(['Case_',num2str(ww)]).StrideLength_opt = ...
            StrideLength_opt;
        Results_tracking.(['Case_',num2str(ww)]).StepWidth_opt_mean = ...
            StepWidth_opt_mean;
        Results_tracking.colheaders.joints = joints;
        Results_tracking.colheaders.GRF = {'fore_aft_r','vertical_r',...
            'lateral_r','fore_aft_l','vertical_l','lateral_l'};
        for i = 1:NMuscle/2
                Results_tracking.colheaders.muscles{i} = ...
                    [muscleNames{i}(1:end-2),'_l'];
                Results_tracking.colheaders.muscles{i+NMuscle/2} = ...
                    [muscleNames{i}(1:end-2),'_r'];
        end
        % Save data
        save([pathresults,'\',namescript,'\Results_tracking.mat'],...
            'Results_tracking');
        end 
    end
end
end