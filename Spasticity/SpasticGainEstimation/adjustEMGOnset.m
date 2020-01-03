% This scripts allows users to manually adjust the EMG onset, which is
% sometimes misidentified, likely due to filtering.

clc
close all
clear all

%% Paths 
pathmain = pwd;
[pathSpasticity,~,~] = fileparts(pathmain);
addpath(genpath(pathSpasticity));
[pathRepo,~,~] = fileparts(pathSpasticity);
pathOpenSimModel = [pathRepo,'\OpenSimModel\'];

subjects_names = {'subject1'};
for sn = 1:length(subjects_names) 
    subject_name = subjects_names{sn};
    % Load data from subjects
    switch subject_name
        case 'subject1'
            subject = subject_name;
            load([pathOpenSimModel,subject,'\IPSA\IPSA_data.mat']);
            segment_sel = 1:22;            
    end      
    savefolder = [pathOpenSimModel,subject,'\Spasticity\SpasticGainEstimation'];
    Nsegment = length(segment_sel);              
    onset_man = zeros(Nsegment,1);                
    for ms = 1:Nsegment     
        % They changed the way ISA are structured, now GAS corresponds to 1
        % MEH corresponds to 2, REF corresponds to 3
        if strcmp(subject,'subject1') && ...
                strcmp(allsegments(ms).joint_name,'knee_ext')
            onset_idx = 'onset2';
        elseif strcmp(subject,'subject1') && ...
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
        load([pathOpenSimModel,subject,...
            '\EMG\IPSA\',['Stretch_',trial],'\emg1_norm_personalized.mat']);
        figure(segment_sel(ms))
        plot(emg1processednorm(:,1),'linewidth',2); hold on;
        pl = plot([allsegments(ms).(onset_idx) allsegments(ms).(onset_idx)],...
            ylim,'k','linewidth',2);
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
        plot([onset_man(ms) onset_man(ms)],ylim,'g','linewidth',3);
        end
    end
    onset_man = round(onset_man);
    if ~(exist(savefolder,'dir')==7)
        mkdir(savefolder);
    end
    save([savefolder,'\onset_man_user'],'onset_man');
end
    