% This script provides bounds and scaling factors for the design variables.
% The bounds on the joint variables are informed by experimental data.
% The bounds on the remaining variables are fixed.
% The bounds are scaled such that the upper/lower bounds cannot be
% larger/smaller than 1/-1.
%
% Author: Antoine Falisse
% Date: 12/19/2018
%
function [bounds,scaling] = getBounds(...
    Qs,Qs_CP,NMuscles,NMuscles_FLV,nq,jointi,GRF,N,NSyn,NMuscles_Spas)

%% Spline approximation of Qs to get Qdots and Qdotdots
% TD data
Qs_spline.data = zeros(size(Qs.allinterpfilt));
Qs_spline.data(:,1) = Qs.allinterpfilt(:,1);
Qdots_spline.data = zeros(size(Qs.allinterpfilt));
Qdots_spline.data(:,1) = Qs.allinterpfilt(:,1);
Qdotdots_spline.data = zeros(size(Qs.allinterpfilt));
Qdotdots_spline.data(:,1) = Qs.allinterpfilt(:,1);
for i = 2:size(Qs.allinterpfilt,2)
    Qs.datafiltspline(i) = spline(Qs.allinterpfilt(:,1),Qs.allinterpfilt(:,i));
    [Qs_spline.data(:,i),Qdots_spline.data(:,i),...
        Qdotdots_spline.data(:,i)] = ...
        SplineEval_ppuval(Qs.datafiltspline(i),Qs.allinterpfilt(:,1),1);
end
% Filter the accelerations
order = 4;
cutoff_low = 10;
fs=1/mean(diff(Qs_spline.data(:,1)));
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
Qdotdots_spline.data(:,2:end) = filtfilt(af,bf,Qdotdots_spline.data(:,2:end));  

% CP data
Qs_CP_spline.data               = zeros(size(Qs_CP.allinterpfilt));
Qs_CP_spline.data(:,1)          = Qs_CP.allinterpfilt(:,1);
Qdots_CP_spline.data            = zeros(size(Qs_CP.allinterpfilt));
Qdots_CP_spline.data(:,1)       = Qs_CP.allinterpfilt(:,1);
Qdotdots_CP_spline.data         = zeros(size(Qs_CP.allinterpfilt));
Qdotdots_CP_spline.data(:,1)    = Qs_CP.allinterpfilt(:,1);
for i = 2:size(Qs_CP.allinterpfilt,2)
    Qs_CP.datafiltspline(i) = spline(Qs_CP.allinterpfilt(:,1),...
        Qs_CP.allinterpfilt(:,i));
    [Qs_CP_spline.data(:,i),Qdots_CP_spline.data(:,i),...
        Qdotdots_CP_spline.data(:,i)] = ...
        SplineEval_ppuval(Qs_CP.datafiltspline(i),Qs_CP.allinterpfilt(:,1),1);
end
% Filter the accelerations
order = 4;
cutoff_low = 10;
fs=1/mean(diff(Qs_CP_spline.data(:,1)));
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
Qdotdots_CP_spline.data(:,2:end) = ...
    filtfilt(af,bf,Qdotdots_CP_spline.data(:,2:end)); 

%% Qs: based on filtered experimental data (inverse kinematics)
% Pelvis tilt
bounds.Qs.upper(jointi.pelvis.list) = max(max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RX')))),...
    max((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RX')))));
bounds.Qs.lower(jointi.pelvis.list) = min(min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RX')))),...
    min((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RX')))));
% bounds.Qs.data_all(:,jointi.pelvis.list) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RX'));
% Pelvis list
bounds.Qs.upper(jointi.pelvis.rot) = max(max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RY')))),...
    max((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RY')))));
bounds.Qs.lower(jointi.pelvis.rot) = min(min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RY')))),...
    min((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RY')))));
% bounds.Qs.data_all(:,jointi.pelvis.rot) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RY'));
% Pelvis rot
bounds.Qs.upper(jointi.pelvis.tilt) = max(max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RZ')))),...
    max((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RZ')))));
bounds.Qs.lower(jointi.pelvis.tilt) = min(min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RZ')))),...
    min((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RZ')))));
% bounds.Qs.data_all(:,jointi.pelvis.tilt) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RZ'));
% Pelvis tx
bounds.Qs.upper(jointi.pelvis.tx) = max(max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TX')))),...
    max((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TX')))));
bounds.Qs.lower(jointi.pelvis.tx) = min(min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TX')))),...
    min((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TX')))));
% bounds.Qs.data_all(:,jointi.pelvis.tx) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TX'));
% Pelvis ty
bounds.Qs.upper(jointi.pelvis.ty) = max(max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TY')))),...
    max((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TY')))));
bounds.Qs.lower(jointi.pelvis.ty) = min(min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TY')))),...
    min((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TY')))));
% bounds.Qs.data_all(:,jointi.pelvis.ty) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TY'));
% Pelvis tz
bounds.Qs.upper(jointi.pelvis.tz) = max(max((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TZ')))),...
    max((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TZ')))));
bounds.Qs.lower(jointi.pelvis.tz) = min(min((Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TZ')))),...
    min((Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TZ')))));
% bounds.Qs.data_all(:,jointi.pelvis.tz) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TZ'));
% Hip flexion
bounds.Qs.upper(jointi.hip_flex.l) = max([max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_l'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_r'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_l')))]);
bounds.Qs.lower(jointi.hip_flex.l) = min([min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_l'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_r'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_l')))]);
bounds.Qs.upper(jointi.hip_flex.r) = bounds.Qs.upper(jointi.hip_flex.l);
bounds.Qs.lower(jointi.hip_flex.r) = bounds.Qs.lower(jointi.hip_flex.l);
% bounds.Qs.data_all(:,jointi.hip_flex.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_l'));
% bounds.Qs.data_all(:,jointi.hip_flex.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_r'));
% Hip adduction
bounds.Qs.upper(jointi.hip_add.l) = max([max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_l'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_r'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_l')))]);
bounds.Qs.lower(jointi.hip_add.l) = min([min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_l'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_r'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_l')))]);
bounds.Qs.upper(jointi.hip_add.r) = bounds.Qs.upper(jointi.hip_add.l);
bounds.Qs.lower(jointi.hip_add.r) = bounds.Qs.lower(jointi.hip_add.l);
% bounds.Qs.data_all(:,jointi.hip_add.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_l'));
% bounds.Qs.data_all(:,jointi.hip_add.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_r'));
% Hip rotation
bounds.Qs.upper(jointi.hip_rot.l) = max([max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_l'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_r'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_l')))]);
bounds.Qs.lower(jointi.hip_rot.l) = min([min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_l'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_r'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_l')))]);
bounds.Qs.upper(jointi.hip_rot.r) = bounds.Qs.upper(jointi.hip_rot.l);
bounds.Qs.lower(jointi.hip_rot.r) = bounds.Qs.lower(jointi.hip_rot.l);
% bounds.Qs.data_all(:,jointi.hip_rot.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_l'));
% bounds.Qs.data_all(:,jointi.hip_rot.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_r'));
% Knee
bounds.Qs.upper(jointi.knee.l) = max([max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_l'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_r'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_l')))]);
bounds.Qs.lower(jointi.knee.l) = min([min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_l'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_r'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_l')))]);
bounds.Qs.upper(jointi.knee.r) = bounds.Qs.upper(jointi.knee.l);
bounds.Qs.lower(jointi.knee.r) = bounds.Qs.lower(jointi.knee.l);
% bounds.Qs.data_all(:,jointi.knee.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_l'));
% bounds.Qs.data_all(:,jointi.knee.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_r'));
% Ankle
bounds.Qs.upper(jointi.ankle.l) = max([max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_l'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_r'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_l')))]);
bounds.Qs.lower(jointi.ankle.l) = min([min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_l'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_r'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_l')))]);
bounds.Qs.upper(jointi.ankle.r) = bounds.Qs.upper(jointi.ankle.l);
bounds.Qs.lower(jointi.ankle.r) = bounds.Qs.lower(jointi.ankle.l);
% bounds.Qs.data_all(:,jointi.ankle.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_l'));
% bounds.Qs.data_all(:,jointi.ankle.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_r'));
% Subtalar
bounds.Qs.upper(jointi.subt.l) = max([max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_r'))),...
    max(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_l'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_r'))),...
    max(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_l')))]);
bounds.Qs.lower(jointi.subt.l) = min([min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_r'))),...
    min(Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_l'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_r'))),...
    min(Qs_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_l')))]);
bounds.Qs.upper(jointi.subt.r) = bounds.Qs.upper(jointi.subt.l);
bounds.Qs.lower(jointi.subt.r) = bounds.Qs.lower(jointi.subt.l);

% The bounds are extended by twice the absolute difference between upper
% and lower bounds.
Qs_range = abs(bounds.Qs.upper - bounds.Qs.lower);
bounds.Qs.lower = bounds.Qs.lower - 3*Qs_range;
bounds.Qs.upper = bounds.Qs.upper + 3*Qs_range;
% For several joints, we manually adjust the bounds
% lower_torso_TX
bounds.Qs.upper(jointi.pelvis.tx) = 2;  
bounds.Qs.lower(jointi.pelvis.tx) = -2;
% lower_torso_TY
bounds.Qs.upper(jointi.pelvis.ty) = 1.1;  
bounds.Qs.lower(jointi.pelvis.ty) = 0.55;
% lower_torso_TZ
bounds.Qs.upper(jointi.pelvis.tz) = 5;
bounds.Qs.lower(jointi.pelvis.tz) = -5;
bounds.Qs.upper(jointi.trunk.ext) = 45*pi/180;
bounds.Qs.lower(jointi.trunk.ext) = -45*pi/180;    
bounds.Qs.upper(jointi.trunk.ben) = 45*pi/180;
bounds.Qs.lower(jointi.trunk.ben) = -45*pi/180;    
bounds.Qs.upper(jointi.trunk.rot) = 60*pi/180;
bounds.Qs.lower(jointi.trunk.rot) = -60*pi/180;

%% Qdots
% The extreme values are selected as upper/lower bounds, which are then
% further extended.
% Pelvis tilt
bounds.Qdots.upper(jointi.pelvis.list) = max(max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RX')))),...
    max((Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RX')))));
bounds.Qdots.lower(jointi.pelvis.list) = min(min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RX')))),...
    min((Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RX')))));
% Pelvis list
bounds.Qdots.upper(jointi.pelvis.rot) = max(max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RY')))),...
    max((Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RY')))));
bounds.Qdots.lower(jointi.pelvis.rot) = min(min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RY')))),...
    min((Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RY')))));
% Pelvis rotation
bounds.Qdots.upper(jointi.pelvis.tilt) = max(max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RZ')))),...
    max((Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RZ')))));
bounds.Qdots.lower(jointi.pelvis.tilt) = min(min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RZ')))),...
    min((Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RZ')))));
% Pelvis_tx
bounds.Qdots.upper(jointi.pelvis.tx) = max(max((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TX')))),...
    max((Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TX'))))); 
bounds.Qdots.lower(jointi.pelvis.tx) = min(min((Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TX')))),...
    min((Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TX')))));
% Pelvis_ty
bounds.Qdots.upper(jointi.pelvis.ty) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TY'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TY')))); 
bounds.Qdots.lower(jointi.pelvis.ty) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TY'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TY')))); 
% Pelvis_tz
bounds.Qdots.upper(jointi.pelvis.tz) = max(max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TZ'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TZ')))); 
bounds.Qdots.lower(jointi.pelvis.tz) = min(min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TZ'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TZ'))));
% Hip flexion
bounds.Qdots.upper(jointi.hip_flex.l) = max([max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_l'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_r'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_l')))]);
bounds.Qdots.lower(jointi.hip_flex.l) = min([min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_l'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_r'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_l')))]);
bounds.Qdots.upper(jointi.hip_flex.r) = bounds.Qdots.upper(jointi.hip_flex.l);
bounds.Qdots.lower(jointi.hip_flex.r) = bounds.Qdots.lower(jointi.hip_flex.l);
% Hip adduction
bounds.Qdots.upper(jointi.hip_add.l) = max([max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_l'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_r'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_l')))]);
bounds.Qdots.lower(jointi.hip_add.l) = min([min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_l'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_r'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_l')))]);
bounds.Qdots.upper(jointi.hip_add.r) = bounds.Qdots.upper(jointi.hip_add.l);
bounds.Qdots.lower(jointi.hip_add.r) = bounds.Qdots.lower(jointi.hip_add.l);
% Hip rotation
bounds.Qdots.upper(jointi.hip_rot.l) = max([max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_l'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_r'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_l')))]);
bounds.Qdots.lower(jointi.hip_rot.l) = min([min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_l'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_r'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_l')))]);
bounds.Qdots.upper(jointi.hip_rot.r) = bounds.Qdots.upper(jointi.hip_rot.l);
bounds.Qdots.lower(jointi.hip_rot.r) = bounds.Qdots.lower(jointi.hip_rot.l);
% Knee
bounds.Qdots.upper(jointi.knee.l) = max([max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_l'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_r'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_l')))]);
bounds.Qdots.lower(jointi.knee.l) = min([min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_l'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_r'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_l')))]);
bounds.Qdots.upper(jointi.knee.r) = bounds.Qdots.upper(jointi.knee.l);
bounds.Qdots.lower(jointi.knee.r) = bounds.Qdots.lower(jointi.knee.l);
% Ankle
bounds.Qdots.upper(jointi.ankle.l) = max([max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_l'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_r'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_l')))]);
bounds.Qdots.lower(jointi.ankle.l) = min([min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_l'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_r'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_l')))]);
bounds.Qdots.upper(jointi.ankle.r) = bounds.Qdots.upper(jointi.ankle.l);
bounds.Qdots.lower(jointi.ankle.r) = bounds.Qdots.lower(jointi.ankle.l);
% Subtalar
bounds.Qdots.upper(jointi.subt.l) = max([max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_r'))),...
    max(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_l'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_r'))),...
    max(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_l')))]);
bounds.Qdots.lower(jointi.subt.l) = min([min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_r'))),...
    min(Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_l'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_r'))),...
    min(Qdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_l')))]);
bounds.Qdots.upper(jointi.subt.r) = bounds.Qdots.upper(jointi.subt.l);
bounds.Qdots.lower(jointi.subt.r) = bounds.Qdots.lower(jointi.subt.l);
% The bounds are extended by the absolute difference between upper
% and lower bounds.
Qdots_range = abs(bounds.Qdots.upper - bounds.Qdots.lower);
bounds.Qdots.lower = bounds.Qdots.lower - 3*Qdots_range;
bounds.Qdots.upper = bounds.Qdots.upper + 3*Qdots_range;
bounds.Qdots.upper(jointi.trunk.ext) = 1;
bounds.Qdots.lower(jointi.trunk.ext) = -1;    
bounds.Qdots.upper(jointi.trunk.ben) = 1;
bounds.Qdots.lower(jointi.trunk.ben) = -1;    
bounds.Qdots.upper(jointi.trunk.rot) = 1;
bounds.Qdots.lower(jointi.trunk.rot) = -1;    

%% Qdotdots
% The extreme values are selected as upper/lower bounds, which are then
% further extended.
% Pelvis tilt
bounds.Qdotdots.upper(jointi.pelvis.list) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RX'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RX'))));
bounds.Qdotdots.lower(jointi.pelvis.list) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RX'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RX'))));
% Pelvis list
bounds.Qdotdots.upper(jointi.pelvis.rot) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RY'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RY'))));
bounds.Qdotdots.lower(jointi.pelvis.rot) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RY'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RY'))));
% Pelvis rotation
bounds.Qdotdots.upper(jointi.pelvis.tilt) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RZ'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RZ'))));
bounds.Qdotdots.lower(jointi.pelvis.tilt) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RZ'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_RZ'))));
% Pelvis_tx
bounds.Qdotdots.upper(jointi.pelvis.tx) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TX'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TX')))); 
bounds.Qdotdots.lower(jointi.pelvis.tx) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TX'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TX'))));
% Pelvis_ty
bounds.Qdotdots.upper(jointi.pelvis.ty) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TY'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TY'))));
bounds.Qdotdots.lower(jointi.pelvis.ty) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TY'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TY'))));
% Pelvis_tz
bounds.Qdotdots.upper(jointi.pelvis.tz) = max(max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TZ'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TZ'))));
bounds.Qdotdots.lower(jointi.pelvis.tz) = min(min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TZ'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'lower_torso_TZ'))));
% Hip flexion
bounds.Qdotdots.upper(jointi.hip_flex.l) = max([max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_l'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_r'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_l')))]);
bounds.Qdotdots.lower(jointi.hip_flex.l) = min([min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_l'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_r'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_flex_l')))]);
bounds.Qdotdots.upper(jointi.hip_flex.r) = bounds.Qdotdots.upper(jointi.hip_flex.l);
bounds.Qdotdots.lower(jointi.hip_flex.r) = bounds.Qdotdots.lower(jointi.hip_flex.l);
% Hip adduction
bounds.Qdotdots.upper(jointi.hip_add.l) = max([max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_l'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_r'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_l')))]);
bounds.Qdotdots.lower(jointi.hip_add.l) = min([min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_l'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_r'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_add_l')))]);
bounds.Qdotdots.upper(jointi.hip_add.r) = bounds.Qdotdots.upper(jointi.hip_add.l);
bounds.Qdotdots.lower(jointi.hip_add.r) = bounds.Qdotdots.lower(jointi.hip_add.l);
% Hip rotation
bounds.Qdotdots.upper(jointi.hip_rot.l) = max([max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_l'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_r'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_l')))]);
bounds.Qdotdots.lower(jointi.hip_rot.l) = min([min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_l'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_r'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'hip_rot_l')))]);
bounds.Qdotdots.upper(jointi.hip_rot.r) = bounds.Qdotdots.upper(jointi.hip_rot.l);
bounds.Qdotdots.lower(jointi.hip_rot.r) = bounds.Qdotdots.lower(jointi.hip_rot.l);
% Knee
bounds.Qdotdots.upper(jointi.knee.l) = max([max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_l'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_r'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_l')))]);
bounds.Qdotdots.lower(jointi.knee.l) = min([min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_l'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_r'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'knee_flex_l')))]);
bounds.Qdotdots.upper(jointi.knee.r) = bounds.Qdotdots.upper(jointi.knee.l);
bounds.Qdotdots.lower(jointi.knee.r) = bounds.Qdotdots.lower(jointi.knee.l);
% Ankle
bounds.Qdotdots.upper(jointi.ankle.l) = max([max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_l'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_r'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_l')))]);
bounds.Qdotdots.lower(jointi.ankle.l) = min([min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_l'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_r'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'ankle_flex_l')))]);
bounds.Qdotdots.upper(jointi.ankle.r) = bounds.Qdotdots.upper(jointi.ankle.l);
bounds.Qdotdots.lower(jointi.ankle.r) = bounds.Qdotdots.lower(jointi.ankle.l);
% Subtalar
bounds.Qdotdots.upper(jointi.subt.l) = max([max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_r'))),...
    max(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_l'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_r'))),...
    max(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_l')))]);
bounds.Qdotdots.lower(jointi.subt.l) = min([min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_r'))),...
    min(Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_l'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_r'))),...
    min(Qdotdots_CP_spline.data(:,strcmp(Qs_CP.colheaders(1,:),'subt_angle_l')))]);
bounds.Qdotdots.upper(jointi.subt.r) = bounds.Qdotdots.upper(jointi.subt.l);
bounds.Qdotdots.lower(jointi.subt.r) = bounds.Qdotdots.lower(jointi.subt.l);
% The bounds are extended by the absolute difference between upper
% and lower bounds.
Qdotdots_range = abs(bounds.Qdotdots.upper - bounds.Qdotdots.lower);
bounds.Qdotdots.lower = bounds.Qdotdots.lower - Qdotdots_range;
bounds.Qdotdots.upper = bounds.Qdotdots.upper + Qdotdots_range;
bounds.Qdotdots.upper(jointi.trunk.ext) = 10;
bounds.Qdotdots.lower(jointi.trunk.ext) = -10;    
bounds.Qdotdots.upper(jointi.trunk.ben) = 10;
bounds.Qdotdots.lower(jointi.trunk.ben) = -10;    
bounds.Qdotdots.upper(jointi.trunk.rot) = 10;
bounds.Qdotdots.lower(jointi.trunk.rot) = -10;     
    
%% Muscle activations
if NSyn == 99
    bounds.a.lower = zeros(1,NMuscles);
    bounds.a.upper = ones(1,NMuscles);
else
    bounds.a.lower = zeros(1,2*NSyn);
    bounds.a.upper = ones(1,2*NSyn);
end

%% Muscle-tendon forces
bounds.FTtilde.lower = zeros(1,NMuscles_FLV);
bounds.FTtilde.upper = 5*ones(1,NMuscles_FLV);

%% Time derivative of muscle activations
tact = 0.015;
tdeact = 0.06;
if NSyn == 99
    bounds.vA.lower = (-1/100*ones(1,NMuscles))./(ones(1,NMuscles)*tdeact);
    bounds.vA.upper = (1/100*ones(1,NMuscles))./(ones(1,NMuscles)*tact);
else
    bounds.vA.lower = (-1/100*ones(1,2*NSyn))./(ones(1,2*NSyn)*tdeact);
    bounds.vA.upper = (1/100*ones(1,2*NSyn))./(ones(1,2*NSyn)*tact);
end
%% Time derivative of muscle-tendon forces
bounds.dFTtilde.lower = -1*ones(1,NMuscles_FLV);
bounds.dFTtilde.upper = 1*ones(1,NMuscles_FLV);

%% Trunk activations
bounds.a_b.lower = -ones(1,nq.res_trunk);
bounds.a_b.upper = ones(1,nq.res_trunk);

%% Trunk excitations
bounds.e_b.lower = -ones(1,nq.res_trunk);
bounds.e_b.upper = ones(1,nq.res_trunk);

%% GRF
bounds.GRF.lower = min(GRF.val.all(:,2:end));
bounds.GRF.upper = max(GRF.val.all(:,2:end));
% Extend bounds to give some flexibility
GRF_ROM = abs(bounds.GRF.upper - bounds.GRF.lower);
bounds.GRF.lower = bounds.GRF.lower - GRF_ROM;
bounds.GRF.upper = bounds.GRF.upper + GRF_ROM;

%% GRM
bounds.GRM.lower = min(GRF.MorGF.allinterp(:,2:end));
bounds.GRM.upper = max(GRF.MorGF.allinterp(:,2:end));
% Extend bounds to give some flexibility
GRM_ROM = abs(bounds.GRM.upper - bounds.GRM.lower);
bounds.GRM.lower = bounds.GRM.lower - GRM_ROM;
bounds.GRM.upper = bounds.GRM.upper + GRM_ROM;

%% Synergy activations
bounds.syna.lower = zeros(1,NSyn);
bounds.syna.upper = ones(1,NSyn);

%% Synergy weights
bounds.synw.lower.l(:,:)=0; 
bounds.synw.upper.l(:,:)=1;  
bounds.synw.lower.r(:,:)=0;  
bounds.synw.upper.r(:,:)=1;   

%% Muscle activations from muscle-tendon force feedback 
bounds.a_Ff.lower = zeros(1,NMuscles_Spas);
bounds.a_Ff.upper = ones(1,NMuscles_Spas);

%% Muscle activations from time derivative of muscle-tendon force feedback 
bounds.a_dFf.lower = zeros(1,NMuscles_Spas);
bounds.a_dFf.upper = ones(1,NMuscles_Spas);

%% Final time
bounds.tf.lower = 0.1;
bounds.tf.upper = 2;

%% Scaling
% Qs
scaling.Qs      = max(abs(bounds.Qs.lower),abs(bounds.Qs.upper));
bounds.Qs.lower = (bounds.Qs.lower)./repmat(scaling.Qs,N,1);
bounds.Qs.upper = (bounds.Qs.upper)./repmat(scaling.Qs,N,1);
% Qdots
scaling.Qdots      = max(abs(bounds.Qdots.lower),abs(bounds.Qdots.upper));
bounds.Qdots.lower = (bounds.Qdots.lower)./scaling.Qdots;
bounds.Qdots.upper = (bounds.Qdots.upper)./scaling.Qdots;
% Qs and Qdots are intertwined
bounds.QsQdots.lower = zeros(N,2*nq.res);
bounds.QsQdots.upper = zeros(N,2*nq.res);
bounds.QsQdots.lower(:,1:2:end) = bounds.Qs.lower;
bounds.QsQdots.upper(:,1:2:end) = bounds.Qs.upper;
bounds.QsQdots.lower(:,2:2:end) = repmat(bounds.Qdots.lower,N,1);
bounds.QsQdots.upper(:,2:2:end) = repmat(bounds.Qdots.upper,N,1);
scaling.QsQdots                 = zeros(1,2*nq.res);
scaling.QsQdots(1,1:2:end)      = scaling.Qs ;
scaling.QsQdots(1,2:2:end)      = scaling.Qdots ;
% Qdotdots
scaling.Qdotdots = max(abs(bounds.Qdotdots.lower),...
    abs(bounds.Qdotdots.upper));
bounds.Qdotdots.lower = (bounds.Qdotdots.lower)./scaling.Qdotdots;
bounds.Qdotdots.upper = (bounds.Qdotdots.upper)./scaling.Qdotdots;
bounds.Qdotdots.lower(isnan(bounds.Qdotdots.lower)) = 0;
bounds.Qdotdots.upper(isnan(bounds.Qdotdots.upper)) = 0;
% Time derivative of muscle activations
% Fixed scaling factor
scaling.vA = 100;
% Muscle activations
scaling.a = 1;
% Trunk activations
scaling.a_b = 1;
% Trunk excitations
scaling.e_b = 1;
% Muscle activations from muscle-tendon force feedback 
scaling.a_Ff = 1;
% Muscle activations from time derivative of muscle-tendon force feedback
scaling.a_dFf = 1;
% Time derivative of muscle-tendon forces
% Fixed scaling factor
scaling.dFTtilde = 100;
% Trunk torque actuators
% Fixed scaling factor
scaling.TrunkTau = 150;
% Muscle-tendon forces
scaling.FTtilde         = max(...
    abs(bounds.FTtilde.lower),abs(bounds.FTtilde.upper)); 
bounds.FTtilde.lower    = (bounds.FTtilde.lower)./scaling.FTtilde;
bounds.FTtilde.upper    = (bounds.FTtilde.upper)./scaling.FTtilde;
% GRF
scaling.GRF = max(abs(bounds.GRF.lower),abs(bounds.GRF.upper));
bounds.GRF.lower = (bounds.GRF.lower)./scaling.GRF;
bounds.GRF.upper = (bounds.GRF.upper)./scaling.GRF;
% GRM
scaling.GRM      = max(abs(bounds.GRM.lower),abs(bounds.GRM.upper));
bounds.GRM.lower = (bounds.GRM.lower)./scaling.GRM;
bounds.GRM.upper = (bounds.GRM.upper)./scaling.GRM;

end
