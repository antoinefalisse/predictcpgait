% This script provides an inital guess for the design variables.
% The guess is data-informed (DI). We use experimental data to provide an
% initial guess of the joint variables but set constant values to the 
% muscle variable and the arm variables. We use a pre-defined final time 
% that is function of the imposed speed.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function guess = ...
    getGuess(Qs,nq,N,NMuscles,NMuscles_FLV,jointi,scaling,NSyn,NMuscles_Spas)

%% Spline approximation of Qs to get Qdots and Qdotdots
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

%% Qs: data-informed
% Pelvis tilt
guess.Qs(:,jointi.pelvis.list) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RX'));
% Pelvis list
guess.Qs(:,jointi.pelvis.rot) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RY'));
% Pelvis rotation
guess.Qs(:,jointi.pelvis.tilt) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RZ'));
% Pelvis_tx
guess.Qs(:,jointi.pelvis.tx) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TX'));
% Pelvis_ty
guess.Qs(:,jointi.pelvis.ty) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TY'));
% Pelvis_tz
guess.Qs(:,jointi.pelvis.tz) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TZ'));
% Hip flexion
guess.Qs(:,jointi.hip_flex.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_l'));
guess.Qs(:,jointi.hip_flex.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_r'));
% Hip adduction
guess.Qs(:,jointi.hip_add.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_l'));
guess.Qs(:,jointi.hip_add.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_r'));
% Hip rotation
guess.Qs(:,jointi.hip_rot.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_l'));
guess.Qs(:,jointi.hip_rot.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_r'));
% Knee
guess.Qs(:,jointi.knee.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_l'));
guess.Qs(:,jointi.knee.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_r'));
% Ankle
guess.Qs(:,jointi.ankle.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_l'));
guess.Qs(:,jointi.ankle.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_r'));
% Subtalar
guess.Qs(:,jointi.subt.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_l'));
guess.Qs(:,jointi.subt.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_r'));
% Lumbar extension
guess.Qs(:,jointi.trunk.ext) = 0;
% Lumbar bending
guess.Qs(:,jointi.trunk.ben) = 0;
% Lumbar rotation
guess.Qs(:,jointi.trunk.rot) = 0; 

%% Qdots: data-informed
% Pelvis tilt
guess.Qdots(:,jointi.pelvis.list) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RX'));
% Pelvis list
guess.Qdots(:,jointi.pelvis.rot) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RY'));
% Pelvis rotation
guess.Qdots(:,jointi.pelvis.tilt) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RZ'));
% Pelvis_tx
guess.Qdots(:,jointi.pelvis.tx) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TX'));
% Pelvis_ty
guess.Qdots(:,jointi.pelvis.ty) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TY'));
% Pelvis_tz
guess.Qdots(:,jointi.pelvis.tz) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TZ'));
% Hip flexion
guess.Qdots(:,jointi.hip_flex.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_l'));
guess.Qdots(:,jointi.hip_flex.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_r'));
% Hip adduction
guess.Qdots(:,jointi.hip_add.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_l'));
guess.Qdots(:,jointi.hip_add.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_r'));
% Hip rotation
guess.Qdots(:,jointi.hip_rot.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_l'));
guess.Qdots(:,jointi.hip_rot.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_r'));
% Knee
guess.Qdots(:,jointi.knee.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_l'));
guess.Qdots(:,jointi.knee.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_r'));
% Ankle
guess.Qdots(:,jointi.ankle.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_l'));
guess.Qdots(:,jointi.ankle.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_r'));
% Subtalar
guess.Qdots(:,jointi.subt.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_l'));
guess.Qdots(:,jointi.subt.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_r'));
% Lumbar extension
guess.Qdots(:,jointi.trunk.ext) = 0;
% Lumbar bending
guess.Qdots(:,jointi.trunk.ben) = 0;
% Lumbar rotation
guess.Qdots(:,jointi.trunk.rot) = 0; 

%% Qs and Qdots are intertwined
guess.QsQdots = zeros(N,2*nq.res);
guess.QsQdots(:,1:2:end) = guess.Qs;
guess.QsQdots(:,2:2:end) = guess.Qdots;

%% Qdotdots: data-informed
% Pelvis tilt
guess.Qdotdots(:,jointi.pelvis.list) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RX'));
% Pelvis list
guess.Qdotdots(:,jointi.pelvis.rot) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RY'));
% Pelvis rotation
guess.Qdotdots(:,jointi.pelvis.tilt) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_RZ'));
% Pelvis_tx
guess.Qdotdots(:,jointi.pelvis.tx) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TX'));
% Pelvis_ty
guess.Qdotdots(:,jointi.pelvis.ty) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TY'));
% Pelvis_tz
guess.Qdotdots(:,jointi.pelvis.tz) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lower_torso_TZ'));
% Hip flexion
guess.Qdotdots(:,jointi.hip_flex.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_l'));
guess.Qdotdots(:,jointi.hip_flex.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flex_r'));
% Hip adduction
guess.Qdotdots(:,jointi.hip_add.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_l'));
guess.Qdotdots(:,jointi.hip_add.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_add_r'));
% Hip rotation
guess.Qdotdots(:,jointi.hip_rot.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_l'));
guess.Qdotdots(:,jointi.hip_rot.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_rot_r'));
% Knee
guess.Qdotdots(:,jointi.knee.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_l'));
guess.Qdotdots(:,jointi.knee.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_flex_r'));
% Ankle
guess.Qdotdots(:,jointi.ankle.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_l'));
guess.Qdotdots(:,jointi.ankle.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_flex_r'));
% Subtalar
guess.Qdotdots(:,jointi.subt.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_l'));
guess.Qdotdots(:,jointi.subt.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'subt_angle_r'));
% Lumbar extension
guess.Qdotdots(:,jointi.trunk.ext) = 0;
% Lumbar bending
guess.Qdotdots(:,jointi.trunk.ben) = 0;
% Lumbar rotation
guess.Qdotdots(:,jointi.trunk.rot) = 0;

%% Muscle variables
if NSyn == 99
    guess.a = 0.1*ones(N,NMuscles);
    guess.vA = 0.01*ones(N,NMuscles);
else
    guess.a = 0.1*ones(N,2*NSyn);
    guess.vA = 0.01*ones(N,2*NSyn);
end
guess.FTtilde = 0.1*ones(N,NMuscles_FLV);
guess.dFTtilde = 0.01*ones(N,NMuscles_FLV);

%% Back activations
guess.a_b = 0.1*ones(N,nq.res_trunk);
guess.e_b = 0.1*ones(N,nq.res_trunk);

%% Synergy activations
guess.syna = 0.01*ones(N,NSyn);

%% Synergy weights
guess.synw.l = 0.1*ones(NMuscles/2,NSyn);
guess.synw.r = 0.1*ones(NMuscles/2,NSyn);

%% Muscle activations from muscle-tendon force feedback 
guess.a_Ff = 0.1*ones(N,NMuscles_Spas);

%% Muscle activations from time derivative of muscle-tendon force feedback 
guess.a_dFf = 0.1*ones(N,NMuscles_Spas);

%% Final time
guess.tf = 1;

%% Scaling
guess.QsQdots   = guess.QsQdots./repmat(scaling.QsQdots,N,1);
guess.Qdotdots  = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,...
    size(guess.dFTtilde,2));

end