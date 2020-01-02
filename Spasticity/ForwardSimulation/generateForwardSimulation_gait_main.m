% This script generates EMG-driven forward simulations for gait trials. The
% code is inspired from Falisse et al, PLoS ONE (2018):
% https://doi.org/10.1371/journal.pone.0208811
%
% Author: Antoine Falisse
% Date: 1/2/2020
%
clear all
close all
clc

%% Paths 
pathmain = pwd;
[pathSpasticity,~,~] = fileparts(pathmain);
addpath(genpath(pathSpasticity));
[pathRepo,~,~] = fileparts(pathSpasticity);
pathOpenSimModel = [pathRepo,'\OpenSimModel\'];
pathParameterEstimation = [pathRepo,'\ParameterEstimation\'];
pathMuscleModel = [pathRepo,'\MuscleModel\'];

%% Loop over cases
for iii = 1
cond = num2str(iii);
switch cond
    case '1'
        subject = 'subject1';
        side = 'r';
        joint = 'knee';
        segment_sel = [6 7 8];       
end
% Load data from subjects
switch subject 
    case 'subject1'
        load([pathOpenSimModel,subject,'\IPSA\IPSA_data.mat']);
        % Available gait trials
        trialnumbers = {'10','12','13','22'};
        % Affected side
        sidename = 'RIGHT';   
end

%% EMG-driven forward simulations
savefolder = [pathOpenSimModel,subject,'\Spasticity\ForwardSimulations\Gait'];
if ~(exist(savefolder,'dir')==7)
    mkdir(savefolder);
end
% Load already existing data if existing
if (exist([savefolder,'\outputGPOPSFTtilde.mat'],'file')==2)
    load([savefolder,'\outputGPOPSFTtilde']);
    load([savefolder,'\HilldiffGPOPSFTtilde']);
    load([savefolder,'\M_muscles_GPOPSFTtilde']);
    load([savefolder,'\M_knee_GPOPSFTtilde']);
    load([savefolder,'\FTGPOPSFTtilde']);
    load([savefolder,'\FtildeGPOPSFTtilde']);
    load([savefolder,'\dFtildeGPOPSFTtilde']);
    load([savefolder,'\MAGPOPSFTtilde']);
    load([savefolder,'\lMtildeGPOPSFTtilde']);
    load([savefolder,'\vMtildeGPOPSFTtilde']);
end
% Load EMG data (non-normalized)
pathEMG = [pathOpenSimModel,subject,'\EMG\Gait\'];
load([pathEMG,'EMG_filt'],'EMG_filt')
% Loop over gait trials
for tn = 1:length(trialnumbers)       
    % Various inputs 
    trialnumber = trialnumbers{tn};
    trialnumbername = ['Gait_',trialnumber]; 
    trialnumbername_lc = ['gait_',trialnumber]; 
    switch sidename
        case 'LEFT'
            vec = 13:24;
            side = 'l'; Uletter_side = 'L';
            vecAll = 1:43;
        case 'RIGHT'
            vec = 1:12;
            side = 'r'; Uletter_side = 'R';
            vecAll = 44:86;
    end
    % We only perform the simulation for certain muscles for which we have
    % EMG scale factors and for which we want to estimate the spastic
    % response
    muscleNamesSel = {'bi_fem_lh';'gas_lat';'gas_med';'rectus_fem';...
        'semimem';'semiten'};
    getMuscleNamesIndices
    % Extract muscle-tendon lengths and moment arms
    path_MuscleAnalysis = [pathOpenSimModel,subject,...
        '\MuscleAnalysis\Gait\Gait_',trialnumber,side,...
        '\Gait_',trialnumber,side,'_'];      
    [lMT,MA,time_MA,~] = readMuscleAnalysis(path_MuscleAnalysis,MKnee);
    lMT_side = lMT(:,vec);
    MA_side = MA(:,vec);
    idx_MselInMKnee_side = zeros(1,length(muscleNamesSel));
    for ii = 1:length(muscleNamesSel)
        idx_MselInMKnee_side(ii) = ...
            find(strcmp(MKnee_side,[muscleNamesSel{ii},'_',side]));
    end
    lMT_sel = lMT_side(:,idx_MselInMKnee_side);
    MA_sel = MA_side(:,idx_MselInMKnee_side);    
    % Normalize EMG
    muscleChannels = {'BIF','GAS','GAS','REF','MEH','MEH'};    
    EMG_ordered = zeros(size(EMG_filt.(trialnumbername_lc).data,1),...
        length(muscleChannels));
    EMG_ordered_norm1 = zeros(size(EMG_filt.(trialnumbername_lc).data,1),...
        length(muscleChannels));
    EMG_ordered_norm2 = zeros(size(EMG_filt.(trialnumbername_lc).data,1),...
        length(muscleChannels));
    load([pathParameterEstimation,'\Results\',...
        subject,'\ParameterEstimationResults.mat']); 
    for j = 1:length(muscleChannels)
        % Select filtered EMG data
        EMG_ordered(:,j) = EMG_filt.(trialnumbername_lc).data(:,strcmp(...
            EMG_filt.(trialnumbername_lc).colheaders,...
            [Uletter_side,muscleChannels{j}]));  
        % We first need to normalize the EMGs to the max value
        % during gait as was done during the parameter optimization
        EMG_ordered_norm1(:,j) = EMG_ordered(:,j)...
            ./ParameterEstimationResults.Scale.maxEMG_s.(side)(strcmp(...
            ParameterEstimationResults.Scale.EMGchannels,...
            muscleChannels{j}));
        % We then normalize
        EMG_ordered_norm2(:,j) = EMG_ordered_norm1(:,j)...
            .*ParameterEstimationResults.Scale.values_s.(side)(strcmp(...
            ParameterEstimationResults.Scale.muscles,...
            [muscleNamesSel{j},'_']));    
    end
    % Interpolate
    EMG_ordered_norm2_interp = ...
        interp1(round(EMG_filt.(trialnumbername_lc).data(:,1),4),...
        EMG_ordered_norm2,round(time_MA,4));    
    % MT-parameters optimized during parameter optimization
    load([pathParameterEstimation,'\Results\',subject,...
        '\MTParameters_personalized.mat']);    
    MTParameters_side = MTParameters(:,vecAll);  
    idx_MselInMall_side = zeros(1,length(muscleNamesSel));
    for ii = 1:length(muscleNamesSel)
        idx_MselInMall_side(ii) = ...
            find(strcmp(Mall_side,[muscleNamesSel{ii},'_',side]));
    end
    params_scaled = MTParameters_side(:,idx_MselInMall_side);    
    % Perform EMG-driven forward simulations
    [outputGPOPS.(trialnumbername),HilldiffGPOPS.(trialnumbername),...
        M_muscles_GPOPS.(trialnumbername),...
        M_knee_GPOPS.(trialnumbername),FTGPOPS.(trialnumbername),...
        FtildeGPOPS.(trialnumbername),dFtildeGPOPS.(trialnumbername),...
        MAGPOPS.(trialnumbername),lMtildeGPOPS.(trialnumbername),...
        vMtildeGPOPS.(trialnumbername)] = ForwardSimulation_main(...
        time_MA,EMG_ordered_norm2_interp,params_scaled,lMT_sel,MA_sel,...
        pathMuscleModel);
end
% Save data
save([savefolder,'\outputGPOPSFTtilde'],'outputGPOPS');
save([savefolder,'\HilldiffGPOPSFTtilde'],'HilldiffGPOPS');
save([savefolder,'\M_muscles_GPOPSFTtilde'],'M_muscles_GPOPS');
save([savefolder,'\M_knee_GPOPSFTtilde'],'M_knee_GPOPS');
save([savefolder,'\FTGPOPSFTtilde'],'FTGPOPS');
save([savefolder,'\FtildeGPOPSFTtilde'],'FtildeGPOPS');
save([savefolder,'\dFtildeGPOPSFTtilde'],'dFtildeGPOPS');
save([savefolder,'\MAGPOPSFTtilde'],'MAGPOPS');
save([savefolder,'\lMtildeGPOPSFTtilde'],'lMtildeGPOPS');
save([savefolder,'\vMtildeGPOPSFTtilde'],'vMtildeGPOPS');

end
