% This script applies the spasticity models calibrated based on sensory
% information from IPSA during gait. The code is inspired from Falisse et al,
% PLoS ONE (2018): https://doi.org/10.1371/journal.pone.0208811
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
pathOpenSimModel = [pathRepo,'/OpenSimModel/'];

% User setting
savecorrelation = 1; % 1 to save the correlation results, 0 otherwise

for iii = 1:2
% Different cases
cond = num2str(iii);
switch cond
    case '1'
        subject_name = 'subject1';
        side = 'r';
        joint = 'knee';
        segment_sel = [6 7 8];
    case '2'
        subject_name = 'subject1';
        side = 'r';
        joint = 'ankle';
        segment_sel = [19 20 21];
end
% Load data from subjects
switch subject_name 
    case 'subject1'
        subject = subject_name;
        load([pathOpenSimModel,subject,'/IPSA/IPSA_data.mat']);
        % Available gait trials
        trialnumbers = {'10','12','13','22'};
        % Affected side
        sidename = 'RIGHT';
end
% Paths to saved data
savefolder = [pathOpenSimModel,subject,'/Spasticity/',...
    'ForwardSimulations/Gait'];
savefolder_sge = [pathOpenSimModel,subject,'/Spasticity/',...
    'SpasticGainEstimation'];
pathParameterEstimation = ...
        [pathOpenSimModel,'/',subject,'/ParameterEstimation/'];
% Pre-allocation
trialnumber = trialnumbers{1};
trialnumbername = ['Gait_',trialnumber];
ulMf_v = struct(trialnumbername,[]);
uvMf_v = struct(trialnumbername,[]);
ulMf_a = struct(trialnumbername,[]);
uvMf_a = struct(trialnumbername,[]);
uaMf_a = struct(trialnumbername,[]);
uFf = struct(trialnumbername,[]);
udFf = struct(trialnumbername,[]);
EMG_interp = struct(trialnumbername,[]);
ulMf_v_interp = struct(trialnumbername,[]);
uvMf_v_interp = struct(trialnumbername,[]);
ulMf_a_interp = struct(trialnumbername,[]);
uvMf_a_interp = struct(trialnumbername,[]);
uaMf_a_interp = struct(trialnumbername,[]);
uFf_interp = struct(trialnumbername,[]);
udFf_interp = struct(trialnumbername,[]); 
ind_TO = zeros(1,length(trialnumbers));    
% Loop over gait trials
% Load EMG data (non-normalized)
pathEMG = [pathOpenSimModel,subject,'/EMG/Gait/'];
load([pathEMG,'EMG_filt'],'EMG_filt')
for tn = 1:length(trialnumbers)
    % Various inputs 
    trialnumber = trialnumbers{tn};
    trialnumbername = ['Gait_',trialnumber];
    trialnumbername_lc = ['gait_',trialnumber];
    threshold = 20;    
    % Extract indices of initial contact and toe off of both legs    
    path_GRF = ...
        [pathOpenSimModel,subject,'/GRF/Gait/GRF_',trialnumbername,'.mot'];        
    [IC_L,IC_R,TO_L,TO_R] =extractStancePhases_MRI(path_GRF,threshold,sidename);
    switch sidename
        case 'LEFT'
            vec = 13:24;
            side = 'l'; Uletter_side = 'L';
            IC = IC_L; TO = TO_L;
        case 'RIGHT'
            vec = 1:12;
            side = 'r'; Uletter_side = 'R';
            IC = IC_R; TO = TO_R;
    end
    getMuscleNamesIndices_MRI
    % Extract muscle-tendon lengths and moment arms
    path_MuscleAnalysis = [pathOpenSimModel,subject,'/MuscleAnalysis/Gait/',...
        'Gait_',trialnumber,side,'/Gait_',trialnumber,side,'_'];
    [~,~,time_MA,~] = extractOSInfo_MRI(path_MuscleAnalysis,MKnee);    
    % Load and normalize EMG data
    muscleNamesSel = {'bi_fem_lh';'gas_lat';'gas_med';'rectus_fem';...
        'semimem';'semiten'};
    muscleChannels = {'BIF','GAS','GAS','REF','MEH','MEH'};    
    EMG_ordered = zeros(size(EMG_filt.(trialnumbername_lc).data,1),...
        length(muscleChannels));
    EMG_ordered_norm1 = zeros(size(EMG_filt.(trialnumbername_lc).data,1),...
        length(muscleChannels));
    EMG_ordered_norm2 = zeros(size(EMG_filt.(trialnumbername_lc).data,1),...
        length(muscleChannels));
    load([pathParameterEstimation,'/ParameterEstimationResults.mat']); 
    for j = 1:length(muscleChannels)
        % Select filtered EMG data
        EMG_ordered(:,j) = EMG_filt.(trialnumbername_lc).data(:,strcmp(...
            EMG_filt.(trialnumbername_lc).colheaders,...
            [Uletter_side,muscleChannels{j}]));  
        % We first need to normalize the EMGs to the max value
        % during gait as was done during the parameter optimization
        EMG_ordered_norm1(:,j) = EMG_ordered(:,j)...
            ./ParameterEstimationResults.Scale.maxEMG_s.(side)(strcmp(...
            ParameterEstimationResults.Scale.EMGchannels,muscleChannels{j}));
        % We then normalize
        EMG_ordered_norm2(:,j) = EMG_ordered_norm1(:,j)...
            .*ParameterEstimationResults.Scale.values_s.(side)(strcmp(...
            ParameterEstimationResults.Scale.muscles,[muscleNamesSel{j},'_']));
    end   
    % Load EMG-driven forward simulations
    load([savefolder,'/outputGPOPSFTtilde']);
    load([savefolder,'/FtildeGPOPSFTtilde']);
    load([savefolder,'/dFtildeGPOPSFTtilde']);
    % Get acceleration    
    % 1) Get continuous time vector
    time_GPOPS  = outputGPOPS.(trialnumbername).result.solution.phase.time;   
    time_GPOPScont = time_GPOPS(1):(time_GPOPS(end)- time_GPOPS(1))/...
        (length(time_GPOPS)-1):time_GPOPS(end);
    time_GPOPScont = time_GPOPScont';
    
    %% Spastic contribution from force model
    % Load data from spastic gain estimation (sge)  
    formulation_fM = 'forceModel';
    number_save_fM = [formulation_fM '_' int2str(segment_sel)];
    load([savefolder_sge,'/',formulation_fM,'/output_',joint,'_',...
        number_save_fM]);  
    % Extract data        
    NMuscles_spas = output.result.setup.auxdata.NMuscles_spas;
    Muscles_spas = output.result.setup.auxdata.Muscles_spas;
    gFf = output.result.solution.parameter(:,1:NMuscles_spas);
    gdFf = output.result.solution.parameter(:,NMuscles_spas+1:2*NMuscles_spas);
    tauFf = output.result.setup.auxdata.Ff_td;
    taudFf = output.result.setup.auxdata.dFf_td;         
    threshold_Ff = output.result.setup.auxdata.threshold_Ff;
    threshold_dFf = output.result.setup.auxdata.threshold_dFf;   
    % The thresholds during gait are the lowest thresholds during passive
    % motions (IPSA)
    threshold_dFf_gait = min(threshold_dFf);
    threshold_Ff_gait = min(threshold_Ff);
    % Gastrocnemii have different indices when considering the ankle (GM is
    % 1 and GL is 2) or the ankle (GM is 5 and GL is 4)  
    switch joint
        case {'knee'}
            for m = 1:NMuscles_spas
                Muscles_spas = [1,5,6]; % BFLH, SM, and ST
                mm = Muscles_spas(m);
                uFf.(trialnumbername)(:,m) = forceFDynamicsExplicit(...
                    time_GPOPS,FtildeGPOPS.(trialnumbername)(:,mm),...
                    tauFf,gFf(m),100,threshold_Ff_gait(m));
                udFf.(trialnumbername)(:,m) = dFdtFDynamicsExplicit(...
                    time_GPOPS,dFtildeGPOPS.(trialnumbername)(:,mm),...
                    taudFf,gdFf(m),100,threshold_dFf_gait(m));
            end
        case 'ankle'
            for m = 1:NMuscles_spas
                Muscles_spas = [3,2]; % GM and GL
                mm = Muscles_spas(m);
                uFf.(trialnumbername)(:,m) = forceFDynamicsExplicit(...
                    time_GPOPS,FtildeGPOPS.(trialnumbername)(:,mm),...
                    tauFf,gFf(m),100,threshold_Ff_gait(m));
                udFf.(trialnumbername)(:,m) = dFdtFDynamicsExplicit(...
                    time_GPOPS,...
                    dFtildeGPOPS.(trialnumbername)(:,mm),taudFf,...
                    gdFf(m),100,threshold_dFf_gait(m));
            end
    end
    %% Interpolation over gait cycle 
    % Index initial contact
    roundv = 2;
    ind_IC = find(round(time_MA,roundv) == round(IC,roundv),1,'first');
    ind_IC_GPOPS = find(round(time_GPOPS,roundv) == round(IC,roundv),1,'first');
    ind_IC_GPOPScont = ...
        find(round(time_GPOPScont,roundv) == round(IC,roundv),1,'first');
    % Sampling frequency
    fs = 1/mean(diff(time_MA));
    fs_GPOPS = 1/mean(diff(time_GPOPS));
    fs_GPOPScont = 1/mean(diff(time_GPOPScont));
    % Here we correct for some margin (200ms) manually taken when selecting
    % the gait cycle
    shift = round(fs*0.2);    
    shift_GPOPS = round(fs_GPOPS*0.2);
    shift_GPOPScont = round(fs_GPOPScont*0.2);
    % Interpolation
    step = (time_MA(end-shift)-time_MA(ind_IC))/999;
    interval = time_MA(ind_IC):step:time_MA(end-shift);
    
    interval_GPOPS = round(time_GPOPS(ind_IC_GPOPS),roundv):step:...
        round(time_GPOPS(end-shift_GPOPS),roundv);
    
    interval_GPOPScont = round(time_GPOPScont(ind_IC_GPOPScont),roundv):...
        step:round(time_GPOPScont(end-shift_GPOPScont),roundv);
    time_cycle = time_MA(ind_IC:end-shift);
    time_cycle_GPOPS = round(time_GPOPS(ind_IC_GPOPS:end-shift_GPOPS),4);
    time_cycle_GPOPScont = round(time_GPOPScont(ind_IC_GPOPScont:...
        end-shift_GPOPScont),4);
    ind_TO(tn) = find(round(interval,2) == round(TO,2),1,'first');    
    % Interpolation EMG
    EMG_interp.(trialnumbername) = interp1(...
        round(EMG_filt.(trialnumbername_lc).data(:,1),4),...
        EMG_ordered_norm2,round(interval,4));
    % Interpolation spastic contribution: force model
    for m = 1:NMuscles_spas
        uFf_interp.(trialnumbername)(:,m) = interp1(time_cycle_GPOPS,...
            uFf.(trialnumbername)(ind_IC_GPOPS:end-shift_GPOPS,m),...
            interval_GPOPS);
        udFf_interp.(trialnumbername)(:,m) = interp1(time_cycle_GPOPS,...
            udFf.(trialnumbername)(ind_IC_GPOPS:end-shift_GPOPS,m),...
            interval_GPOPS);
    end
end   

%% Plot spastic contribution during gait
switch joint
    case 'knee'
        name_muscles_spastic = {'Biceps femoris lh',...
                'Semimembranosus','Semitendinosus'};
        Muscles_spas_in_knee = [1,5,6]; % BFLH, SM, and ST
    case 'ankle'
        name_muscles_spastic = {'Gas med','Gas lat'};
        Muscles_spas_in_knee = [3,2]; % GM and GL
end
% Mean index TO
ind_TO_mean_round = round(mean(ind_TO(:,1)/10));
% Pre-allocation
trialnumber = trialnumbers{1};
trialnumbername = ['Gait_',trialnumber];
uFf_dFf_interp = struct(trialnumbername,[]);
ulf_vf_interp = struct(trialnumbername,[]);
ulf_vf_af_interp = struct(trialnumbername,[]);
figure()
x = linspace(0,100,1000);
count = 0;
for tn = 1:length(trialnumbers)
    for m = 1:NMuscles_spas
        count = count +1;
        trialnumber = trialnumbers{tn};
        trialnumbername = ['Gait_',trialnumber];
        uFf_dFf_interp.(trialnumbername)(:,m) = ...
            udFf_interp.(trialnumbername)(:,m) + ...
            uFf_interp.(trialnumbername)(:,m);
        uFf_dFf_interp.(trialnumbername)(...
            uFf_dFf_interp.(trialnumbername)<0) = 0;
        % Plots
        subplot(length(trialnumbers),NMuscles_spas,count)
        rr(1) = plot(x,EMG_interp.(trialnumbername)...
            (:,Muscles_spas_in_knee(m)),'k','linewidth',3); hold on;
        rr(2) = plot(x,uFf_dFf_interp.(trialnumbername)(:,m),'m',...
            'linewidth',2);
        ylim([0 1]);
        set(gca,'Fontsize',16);
        set(gca, 'XTickLabel', [])
        if count > length(trialnumbers)*(NMuscles_spas)-(NMuscles_spas)
            NumTicks = 3;
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks))
            xlabel('Gait cycle (%)','Fontsize',16);
            set(gca, 'XTickLabel',[0 50 100])
        end
        NumTicks = 2;
        L = get(gca,'YLim');
        set(gca,'YTick',linspace(L(1),L(2),NumTicks))            
        plot([ind_TO_mean_round ind_TO_mean_round],ylim,'k','linestyle',...
            ':','linewidth',3);
        ylabel('m. exc','Fontsize',16);
        if count < NMuscles_spas+1
            title(name_muscles_spastic{m},'Fontsize',16);
        end

    end
end
l = legend('EMG','Force model');
set(l,'Fontsize',16);
sp = suptitle('Comparison EMG with sensory feedback during gait');
set(sp,'Fontsize',20);
    
%% Compute correlations between EMG and spastic contributions during gait
if strcmp(joint,'knee')
    muscle_names = {'BFLH','SM','ST'};
elseif strcmp(joint,'ankle')
    muscle_names = {'GM','GL'};
end
% Pre-allocation
rcorr_gait = struct('uFf_dFf',[]);
for m = 1:NMuscles_spas
    for tn = 1:length(trialnumbers)
        trialnumber = trialnumbers{tn};
        trialnumbername = ['Gait_',trialnumber]; 
        rcorr_gait.uFf_dFf.(trialnumbername).(muscle_names{m}) = ...
            xcorr(EMG_interp.(trialnumbername)...
            (:,Muscles_spas_in_knee(m)),...
            uFf_dFf_interp.(trialnumbername)(:,m),0,'coeff');                   
        rcorr_gait.uFf_dFf.all.(muscle_names{m})(tn) = ...
            rcorr_gait.uFf_dFf.(trialnumbername).(muscle_names{m});
    end
    rcorr_gait.uFf_dFf.all_mean.(muscle_names{m}) = ...
        mean(rcorr_gait.uFf_dFf.all.(muscle_names{m}));
    rcorr_gait.uFf_dFf.all_std.(muscle_names{m}) = ...
        std(rcorr_gait.uFf_dFf.all.(muscle_names{m}));
end
% Save correlation for statistical analysis
if savecorrelation
    saveCorrGait = [pathOpenSimModel,subject,'/Spasticity/',...
        'SpasticityGait/']; 
    if ~(exist(saveCorrGait,'dir')==7)
        mkdir(saveCorrGait);
    end
    save([saveCorrGait,'/rcorr_gait_',joint,'_',side,'_',...
        int2str(segment_sel)],'rcorr_gait');
end

end
