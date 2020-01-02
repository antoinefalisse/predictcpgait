% % Function for Hill-equilibrium
% import casadi.*
% % States 
% NMuscle_SX = SX.sym('NMuscle_SX',1);
% a_SX = SX.sym('a_SX',NMuscle_SX); % Muscle activations
% fse_SX = SX.sym('fse_SX',NMuscle_SX); % Normalized tendon forces
% % Controls 
% dfse_SX = SX.sym('dfse_SX',NMuscle_SX); % Derivative tendon forces
% % Parameters
% Fmax_SX = SX.sym('Fmax_SX',NMuscle_SX);
% lMopt_SX = SX.sym('lMopt_SX',NMuscle_SX);
% lTs_SX = SX.sym('lTs_SX',NMuscle_SX);
% alpha_opt_SX = SX.sym('alpha_opt_SX',NMuscle_SX);
% aTendon_SX = SX.sym('aTendon_SX',NMuscle_SX);
% shift_SX = SX.sym('shift_SX',NMuscle_SX);
% lMT_SX = SX.sym('lMT',NMuscle_SX); % Muscle-tendon lengths
% vMT_SX = SX.sym('vMT',NMuscle_SX); % Muscle-tendon velocities
% % Parameters of force-length-velocity curves
% load Fvparam
% load Fpparam
% load Faparam
% % Results
% Hilldiff_SX = SX(NMuscle_SX,1); % Hill equilibrium
% FT_SX = SX(NMuscle_SX,1); % Tendon forces
% for m = 1:NMuscle_SX 
%     [Hilldiff_SX, FT_SX] = ...
%         ForceEquilibrium_FtildeState_ParamEst(a_SX,fse_SX,dfse_SX,lMT_SX,...
%         vMT_SX,Fmax_SX,lMopt_SX,lTs_SX,alpha_opt_SX,Fvparam,Fpparam,Faparam,...
%         aTendon_SX,shift_SX);       
% end
% f_ForceEquilibrium_FtildeState_ParamEst = ...
%     Function('f_ForceEquilibrium_FtildeState_ParamEst',{a_SX,fse_SX,dfse_SX,...
%     lMT_SX,vMT_SX,Fmax_SX,lMopt_SX,lTs_SX,alpha_opt_SX,aTendon_SX,shift_SX},...
%     {Hilldiff_SX,FT_SX});


% Function for Hill-equilibrium
import casadi.*
a_SX           = SX.sym('a_SX',1,NselMusc_all);
lMT_SX         = SX.sym('lMT_SX',1,NselMusc_all);
vMT_SX         = SX.sym('vMT_SX',1,NselMusc_all);
lTs_SX         = SX.sym('lTs_SX',NselMusc_all);
lMo_SX         = SX.sym('lMo_SX',NselMusc_all);
% Parameters of force-length-velocity curves
load Fvparam
load Fpparam
load Faparam

% it does not matter which side we pick for the pennation angle since they
% are the same for both sides.
[FM_SX,Fce_SX,Fpe_SX,lMtilde_SX] = ...
    HillModel_RigidTendon(a_SX,lMT_SX,vMT_SX,ParametersInit.l.PennAng(...
    :,SelectedMusclesInfo{1}.Index),lTs_SX',lMo_SX',Fvparam,Fpparam,Faparam);

f_HillModel_RigidTendon = ...
    Function('f_HillModel_RigidTendon',{a_SX,lMT_SX,vMT_SX,...
    lTs_SX,lMo_SX},{FM_SX,Fce_SX,Fpe_SX,lMtilde_SX});
