% This scripts allows the user to manually adjust the EMG onset, which is
% sometimes misidentified, likely due to filtering.

clc
close all
clear all

lMtilde_ext.max = 1.5;
ls = 'simm_gen';
meth = 'cas';
optMuscles = 'out_18';
lMt_ub = num2str(10*lMtilde_ext.max);

% Add folder to MATLAB search path 
pathmain = pwd;
[pathRepo,~,~] = fileparts(pathmain);
pathRepo_simcpspasticity = 'C:\Users\u0101727\Documents\MyRepositories\simcpspasticity_cases';
addpath(genpath(pathRepo));
pathOpenSimModel = [pathRepo,'\OpenSimModel\'];


subjects_names = {'EF_r'};

for sn = 1%length(subjects_names) 
    subject_name = subjects_names{sn};
    % Load data from subjects
    switch subject_name
        case 'EF_r'
            real_subject_name = 'EF';
            model = 'FrEu_pre_MRI_scaledtorso'; 
            path_folder_HD = [pathRepo_simcpspasticity,'\Alldata\'];
            load([path_folder_HD,'\IPSA\EF\001217b123_segments_deleted.mat']);
            segment_sel = 1:22;            
    end      
    savefolder = [pathOpenSimModel,real_subject_name,...
        '\Spasticity\SpasticGainEstimation'];
    Nsegment = length(segment_sel);              
    onset_man = zeros(Nsegment,1);                
    for ms = 1:Nsegment     
        % They changed the way ISA are structured, now GAS corresponds to 1
        % MEH corresponds to 2, REF corresponds to 3
        if strcmp(real_subject_name,'EF') && ...
                strcmp(allsegments(ms).joint_name,'knee_ext')
            onset_idx = 'onset2';
        elseif strcmp(real_subject_name,'EF') && ...
                strcmp(allsegments(ms).joint_name,'knee_flex')
            onset_idx = 'onset3';
        else
            onset_idx = 'onset1';
        end           
            
        if isnan(allsegments(ms).(onset_idx))
            onset_man(ms) = NaN;
        % No low velocity trials since not used for estimation
        elseif (strcmp(allsegments(ms).meas_type,'low') == 1)
            onset_man(ms) = NaN;
        % No medium velocity trials since not used for estimation
        elseif (strcmp(allsegments(ms).meas_type,'medium') == 1) 
            onset_man(ms) = NaN;    
        % No adductors since not included in study
        elseif (strcmp(allsegments(ms).joint_name,'adduct') == 1)
            onset_man(ms) = NaN;
         % special case for LP due to normalization issues with hamstrings
        elseif (strcmp(subject_name,'LP') == 1 && ...
                strcmp(allsegments(ms).joint_name,'hamstrings') == 1)
            onset_man(ms) = NaN;
        else
        trial = ['segment_' int2str(segment_sel(ms))];
        % Load EMG
        load([pathOpenSimModel,real_subject_name,...
            '\EMG\',['Stretch_',trial],'\emg1_norm_',optMuscles,'_ublM',...
            lMt_ub,'_ls_',ls,'_meth_',meth,'.mat']);
        figure(segment_sel(ms))
        plot(emg1processednorm(:,1),'linewidth',2); hold on;
        pl = plot([allsegments(ms).(onset_idx) allsegments(ms).(onset_idx)],ylim,...
            'k','linewidth',2);
        l = legend(pl,'Automatically identified EMG onset');
        set(l,'Fontsize',16);
        title(['Select EMG onset, if difference with automatically', ...
            ' identified is > 500 x-frames then automatically ', ...
            'identified is kept'],'Fontsize',16);
        [startstop, ~]=ginput(1);
        % if this is already ok, then keep
        if startstop < (allsegments(ms).(onset_idx) - 500) 
            onset_man(ms) = allsegments(ms).(onset_idx);
        else
        % if this seems bad, then adjust
            onset_man(ms) = startstop; 
        end
        plot([onset_man(ms) onset_man(ms)],ylim,'g',...
            'linewidth',3);
        end
    end
    onset_man = round(onset_man);
    if ~(exist(savefolder,'dir')==7)
        mkdir(savefolder);
    end
    save([savefolder,'\onset_man_user'],'onset_man');
end
    