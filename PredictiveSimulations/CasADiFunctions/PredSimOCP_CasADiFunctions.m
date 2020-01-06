% This script contains several CasADi-based functions that are
% used when solving the OCPs
%
% Author: Antoine Falisse
% Date: 12/19/2018
%
import casadi.*

%% Polynomial approximation
pathPolynomials = [pathPredictiveSimulations,'/Polynomials'];
addpath(genpath(pathPolynomials));
% Since the model is asymmetric, we use different functions for both legs
% Right leg
muscle_spanning_info_m_r = muscle_spanning_joint_INFO_r(musi_pol,:);
MuscleInfo_m_r.muscle    = MuscleInfo_r.muscle(musi_pol);                  
qin     = SX.sym('qin',1,nq.leg);
qdotin  = SX.sym('qdotin',1,nq.leg);
lMT     = SX(NMuscles_pol,1);
vMT     = SX(NMuscles_pol,1);
dM      = SX(NMuscles_pol,nq.leg);
for i=1:NMuscles_pol      
    index_dof_crossing  = find(muscle_spanning_info_m_r(i,:)==1);
    order               = MuscleInfo_m_r.muscle(i).order;
    [mat,diff_mat_q]    = n_art_mat_3_cas_SX(qin(1,index_dof_crossing),order);
    lMT(i,1)            = mat*MuscleInfo_m_r.muscle(i).coeff;
    vMT(i,1)            = 0;
    dM(i,1:nq.leg)      = 0;
    nr_dof_crossing     = length(index_dof_crossing); 
    for dof_nr = 1:nr_dof_crossing
        dM(i,index_dof_crossing(dof_nr)) = ...
            (-(diff_mat_q(:,dof_nr)))'*MuscleInfo_m_r.muscle(i).coeff;
        vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing(dof_nr))*...
            qdotin(1,index_dof_crossing(dof_nr)));
    end 
end
f_lMT_vMT_dM_r = Function('f_lMT_vMT_dM_r',{qin,qdotin},{lMT,vMT,dM});
% Left leg
muscle_spanning_info_m_l = muscle_spanning_joint_INFO_l(musi_pol,:);
MuscleInfo_m_l.muscle    = MuscleInfo_l.muscle(musi_pol);                  
qin     = SX.sym('qin',1,nq.leg);
qdotin  = SX.sym('qdotin',1,nq.leg);
lMT     = SX(NMuscles_pol,1);
vMT     = SX(NMuscles_pol,1);
dM      = SX(NMuscles_pol,nq.leg);
for i=1:NMuscles_pol      
    index_dof_crossing  = find(muscle_spanning_info_m_l(i,:)==1);
    order               = MuscleInfo_m_l.muscle(i).order;
    [mat,diff_mat_q]    = n_art_mat_3_cas_SX(qin(1,index_dof_crossing),order);
    lMT(i,1)            = mat*MuscleInfo_m_l.muscle(i).coeff;
    vMT(i,1)            = 0;
    dM(i,1:nq.leg)      = 0;
    nr_dof_crossing     = length(index_dof_crossing); 
    for dof_nr = 1:nr_dof_crossing
        dM(i,index_dof_crossing(dof_nr)) = ...
            (-(diff_mat_q(:,dof_nr)))'*MuscleInfo_m_l.muscle(i).coeff;
        vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing(dof_nr))*...
            qdotin(1,index_dof_crossing(dof_nr)));
    end 
end
f_lMT_vMT_dM_l = Function('f_lMT_vMT_dM_l',{qin,qdotin},{lMT,vMT,dM});

%% Normalized sum of squared values
% Function for 6 elements 
etemp6 = SX.sym('etemp6',6);
Jtemp6 = 0;
for i=1:length(etemp6)
    Jtemp6 = Jtemp6 + etemp6(i).^2;
end
Jtemp6 = Jtemp6/6;
f_J6 = Function('f_J6',{etemp6},{Jtemp6});
% Function for NMact elements
etempNMact = SX.sym('etempNMact',NMact);
JtempNMact = 0;
for i=1:length(etempNMact)
    JtempNMact = JtempNMact + etempNMact(i).^2;
end
JtempNMact = JtempNMact/NMact;
f_JNMact = Function('f_JNMact',{etempNMact},{JtempNMact});
% Function for nq.res-1 elements
etempNq_bpty = SX.sym('etempNq_bpty',nq.res-1);
JtempNq_bpty = 0;
for i=1:length(etempNq_bpty)
    JtempNq_bpty = JtempNq_bpty + etempNq_bpty(i).^2;
end
JtempNq_bpty = JtempNq_bpty/(nq.res-1);
f_JNq_bpty = Function('f_JNq_bpty',{etempNq_bpty},{JtempNq_bpty});
% Function for nq.res_bg elements
etempNq_act = SX.sym('etempNq_act',nq.res_bgp);
JtempNq_act = 0;
for i=1:length(etempNq_act)
    JtempNq_act = JtempNq_act + etempNq_act(i).^2;
end
JtempNq_act = JtempNq_act/nq.res_bgp;
f_JNq_act = Function('f_JNq_act',{etempNq_act},{JtempNq_act});
% Function for nq.res elements
etempNq_all = SX.sym('etempNq_all',nq.res);
JtempNq_all = 0;
for i=1:length(etempNq_all)
    JtempNq_all = JtempNq_all + etempNq_all(i).^2;
end
JtempNq_all = JtempNq_all/nq.res;
f_JNq_all = Function('f_JNq_all',{etempNq_all},{JtempNq_all});
% Function for NMuscles_FLV elements
etempNM_FLV = SX.sym('etempNM_FLV',NMuscles_FLV);
JtempNM_FLV = 0;
for i=1:length(etempNM_FLV)
    JtempNM_FLV = JtempNM_FLV + etempNM_FLV(i).^2;
end
JtempNM_FLV = JtempNM_FLV/NMuscles_FLV;
f_JNM_FLV = Function('f_JNM_FLV',{etempNM_FLV},{JtempNM_FLV});
% Function for nq.res_trunk elements 
etempnq_trunk = SX.sym('etempnq_trunk',nq.res_trunk);
Jtempnq_trunk = 0;
for i=1:length(etempnq_trunk)
    Jtempnq_trunk = Jtempnq_trunk + etempnq_trunk(i).^2;
end
Jtempnq_trunk = Jtempnq_trunk/nq.res_trunk;
f_Jnq_trunk = Function('f_Jnq_trunk',{etempnq_trunk},{Jtempnq_trunk});

%% Normalized sum of values to a certain power
% Function for NMuscles_FLV elements 
etempNM_FLVexp  = SX.sym('etempNM_FLVexp',NMuscles_FLV);
expo            = SX.sym('exp',1);
JtempNM_FLVexp  = 0;
for i=1:length(etempNM_FLVexp)
    JtempNM_FLVexp = JtempNM_FLVexp + etempNM_FLVexp(i).^expo;
end
JtempNM_FLVexp = JtempNM_FLVexp/NMuscles_FLV;
f_JNM_FLVexp = Function('f_JNM_FLVexp',{etempNM_FLVexp,expo},{JtempNM_FLVexp});

% Function for NMuscles elements 
etempNM_exp = SX.sym('etempNM_exp',NMuscles);
expo        = SX.sym('exp',1);
JtempNM_exp  = 0;
for i=1:length(etempNM_exp)
    JtempNM_exp = JtempNM_exp + etempNM_exp(i).^expo;
end
JtempNM_exp = JtempNM_exp/NMuscles;
f_JNM_exp = Function('f_JNM_exp',{etempNM_exp,expo},{JtempNM_exp});

%% Sum of products 
% Function for 27 elements 
ma_temp27 = SX.sym('ma_temp27',27);
ft_temp27 = SX.sym('ft_temp27',27);
J_sptemp27 = 0;
for i=1:length(ma_temp27)
    J_sptemp27 = J_sptemp27 + ma_temp27(i,1)*ft_temp27(i,1);    
end
f_T27 = Function('f_T27',{ma_temp27,ft_temp27},{J_sptemp27});
% Function for 13 elements 
ma_temp13 = SX.sym('ma_temp13',13);
ft_temp13 = SX.sym('ft_temp13',13);
J_sptemp13 = 0;
for i=1:length(ma_temp13)
    J_sptemp13 = J_sptemp13 + ma_temp13(i,1)*ft_temp13(i,1);    
end
f_T13 = Function('f_T13',{ma_temp13,ft_temp13},{J_sptemp13});
% Function for 12 elements 
ma_temp12 = SX.sym('ma_temp12',12);
ft_temp12 = SX.sym('ft_temp12',12);
J_sptemp12 = 0;
for i=1:length(ma_temp12)
    J_sptemp12 = J_sptemp12 + ma_temp12(i,1)*ft_temp12(i,1);    
end
f_T12 = Function('f_T12',{ma_temp12,ft_temp12},{J_sptemp12});

%% Trunk activation dynamics
e_b = SX.sym('e_b',nq.res_trunk); % trunk excitations
a_b = SX.sym('a_b',nq.res_trunk); % trunk activations
dadt = TrunkActivationDynamics(e_b,a_b);
f_TrunkActivationDynamics = ...
    Function('f_TrunkActivationDynamics',{e_b,a_b},{dadt});

%% Muscle contraction dynamics
pathmusclemodel = [pathRepo,'\MuscleModel'];
addpath(genpath(pathmusclemodel));
% Function for Hill-equilibrium
FTtilde     = SX.sym('FTtilde',NMuscles_FLV); % Normalized tendon forces
a           = SX.sym('a',NMuscles_FLV); % Muscle activations
dFTtilde    = SX.sym('dFTtilde',NMuscles_FLV); % Time derivative tendon forces
lMT         = SX.sym('lMT',NMuscles_FLV); % Muscle-tendon lengths
vMT         = SX.sym('vMT',NMuscles_FLV); % Muscle-tendon velocities
tension_SX  = SX.sym('tension',NMuscles_FLV); % Tensions
Hilldiff    = SX(NMuscles_FLV,1); % Hill-equilibrium
FT          = SX(NMuscles_FLV,1); % Tendon forces
Fce         = SX(NMuscles_FLV,1); % Contractile element forces
Fiso        = SX(NMuscles_FLV,1); % Normalized forces from force-length curve
vMmax       = SX(NMuscles_FLV,1); % Maximum contraction velocities
massM       = SX(NMuscles_FLV,1); % Muscle mass
Fpe         = SX(NMuscles_FLV,1); % Passive forces
lMtilde     = SX(NMuscles_FLV,1); % Muscle fiber lengths
% Parameters of force-length-velocity curves
load Fvparam
load Fpparam
load Faparam
for m = 1:NMuscles_FLV
    [Hilldiff(m),FT(m),Fce(m),Fiso(m),vMmax(m),massM(m),Fpe(m),lMtilde(m)] = ...
        ForceEquilibrium_FtildeState_all(a(m),FTtilde(m),dFTtilde(m),...
        lMT(m),vMT(m),MTParameters_FLV(:,m),Fvparam,Fpparam,Faparam,...
        tension_SX(m));
end
f_forceEquilibrium_FtildeState_all = ...
    Function('f_forceEquilibrium_FtildeState_all',{a,FTtilde,dFTtilde,...
    lMT,vMT,tension_SX},{Hilldiff,FT,Fce,Fiso,vMmax,massM,Fpe,lMtilde});

% Function to get (normalized) muscle fiber lengths
lM      = SX(NMuscles_FLV,1);
lMtilde = SX(NMuscles_FLV,1);
for m = 1:NMuscles_FLV
    [lM(m),lMtilde(m)] = FiberLength_TendonForce(FTtilde(m),...
        MTParameters_FLV(:,m),lMT(m));
end
f_FiberLength_TendonForce = Function('f_FiberLength_Ftilde',...
    {FTtilde,lMT},{lM,lMtilde});

% Function to get (normalized) muscle fiber velocities
vM      = SX(NMuscles_FLV,1);
vMtilde = SX(NMuscles_FLV,1);
for m = 1:NMuscles_FLV
    [vM(m),vMtilde(m)] = FiberVelocity_TendonForce(FTtilde(m),...
        dFTtilde(m),MTParameters_FLV(:,m),lMT(m),vMT(m));
end
f_FiberVelocity_TendonForce = Function('f_FiberVelocity_Ftilde',...
    {FTtilde,dFTtilde,lMT,vMT},{vM,vMtilde});

%% Muscle contraction: no force length velocity
a_noFLV     = SX.sym('a_noFLV',NMuscles_noFLV); % Muscle activations
FT_noFLV    = a_noFLV.*MTParameters_noFLV(1,:)';
f_force_noFLV = Function('f_force_noFLV',{a_noFLV},{FT_noFLV});

%% Passive joint torques
K_pass      = SX.sym('K_pass',4);
theta_pass  = SX.sym('theta_pass',2);
qin_pass    = SX.sym('qin_pass',1);
qdotin_pass = SX.sym('qdotin_pass',1);
% theta_pass 1 and 2 are inverted on purpose.
Tau_pass = K_pass(1,1)*exp(K_pass(2,1)*(qin_pass-theta_pass(2,1))) + ...
    K_pass(3,1)*exp(K_pass(4,1)*(qin_pass-theta_pass(1,1))) ...
    - 0.001*qdotin_pass;
f_PassiveTorques = Function('f_PassiveTorques',{K_pass,theta_pass,...
    qin_pass,qdotin_pass},{Tau_pass});

%% Synergy product
synW_test  = SX.sym('synW_test',NMuscles/2*NSyn);
synA_test  = SX.sym('synA_test',NSyn);
synW_test_res = reshape(synW_test,NMuscles/2,NSyn);
asyn_test = synW_test_res*synA_test;
f_SynergyProduct = Function('f_SynergyProduct',...
  {synW_test,synA_test},{asyn_test});

%% Metabolic energy models
act_SX          = SX.sym('act_SX',NMuscles_FLV,1); % Muscle activations
exc_SX          = SX.sym('exc_SX',NMuscles_FLV,1); % Muscle excitations
lMtilde_SX      = SX.sym('lMtilde_SX',NMuscles_FLV,1); % N muscle fiber lengths
vMtilde_SX      = SX.sym('vMtilde_SX',NMuscles_FLV,1); % N muscle fiber vel
vM_SX           = SX.sym('vM_SX',NMuscles_FLV,1); % Muscle fiber velocities
Fce_SX          = SX.sym('FT_SX',NMuscles_FLV,1); % Contractile element forces
Fpass_SX        = SX.sym('FT_SX',NMuscles_FLV,1); % Passive element forces
Fiso_SX         = SX.sym('Fiso_SX',NMuscles_FLV,1); % N forces (F-L curve)
musclemass_SX   = SX.sym('musclemass_SX',NMuscles_FLV,1); % Muscle mass 
vcemax_SX       = SX.sym('vcemax_SX',NMuscles_FLV,1); % Max contraction vel
pctst_SX        = SX.sym('pctst_SX',NMuscles_FLV,1); % Slow twitch ratio 
Fmax_SX         = SX.sym('Fmax_SX',NMuscles_FLV,1); % Max iso forces
modelmass_SX    = SX.sym('modelmass_SX',1); % Model mass
b_SX            = SX.sym('b_SX',1); % Parameter determining tanh smoothness
% Bhargava et al. (2004)
[energy_total_sm_SX,Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,...
    energy_model_sm_SX] = getMetabolicEnergySmooth2004all(exc_SX,act_SX,...
    lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass_SX,pctst_SX,Fiso_SX,...
    Fmax_SX,modelmass_SX,b_SX);
fgetMetabolicEnergySmooth2004all = ...
    Function('fgetMetabolicEnergySmooth2004all',...
    {exc_SX,act_SX,lMtilde_SX,vM_SX,Fce_SX,Fpass_SX,musclemass_SX,...
    pctst_SX,Fiso_SX,Fmax_SX,modelmass_SX,b_SX},{energy_total_sm_SX,...
    Adot_sm_SX,Mdot_sm_SX,Sdot_sm_SX,Wdot_sm_SX,energy_model_sm_SX});

%% Spindle dynamics
if spasi == 1
    a_sf_SX         = SX.sym('a_sf_SX',NMuscles_Spas,1); % activation from feedback
    s_SX            = SX.sym('s_SX',NMuscles_Spas,1); % sensory information
    tau_SX          = SX.sym('tau_SX',NMuscles_Spas,1); % time delay
    g_SX            = SX.sym('g_SX',NMuscles_Spas,1); % feedback gains
    threshold_SX    = SX.sym('threshold_SX',NMuscles_Spas,1); % feddback threshold
    daSfdt_SX = spindleDynamics(a_sf_SX,s_SX,tau_SX,g_SX,b_SX,threshold_SX);
    f_spindleDynamics = Function('f_spindleDynamics',{a_sf_SX,s_SX,tau_SX,g_SX,...
        b_SX,threshold_SX},{daSfdt_SX});
end
