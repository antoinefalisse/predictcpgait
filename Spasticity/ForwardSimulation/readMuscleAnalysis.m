% Read results from muscle analysis

function [lMT,MA,time,knee_angle] = ...
    readMuscleAnalysis(path_MuscleAnalysis,MKnee,path_JointAngles)

lMT_raw = importdata([path_MuscleAnalysis,'MuscleAnalysis_Length.sto']);
MA_raw_l = importdata([path_MuscleAnalysis,...
    'MuscleAnalysis_MomentArm_knee_flex_l.sto']);
MA_raw_r = importdata([path_MuscleAnalysis,...
    'MuscleAnalysis_MomentArm_knee_flex_r.sto']);

time = lMT_raw.data(:,1);

[lMT,MA_l] = getlMTMA(MKnee,lMT_raw,MA_raw_l);
[~,MA_r] = getlMTMA(MKnee,lMT_raw,MA_raw_r);

NMuscles = length(MKnee);

MA(:,1:NMuscles/2) = MA_r(:,1:NMuscles/2);
MA(:,NMuscles/2+1:NMuscles) = MA_l(:,NMuscles/2+1:NMuscles);

if nargin <3
    knee_angle = [];
else
    [~, ~, data, ~] = SIMM_ReadMotion(path_JointAngles);
    knee_angle = data(:,4);
end

end
