% This scripts allows the user to manually adjust the optimization time
% interval

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
addpath(genpath(pathRepo));
pathRepo_simcpspasticity = 'C:\Users\u0101727\Documents\MyRepositories\simcpspasticity_cases';
pathOpenSimModel = [pathRepo,'\OpenSimModel\'];

subjects_names = {'EF_r'};
joints = {'knee','ankle'};
for sn = length(subjects_names)     
    subject_name = subjects_names{sn};
    for j = length(joints)
    joint = joints{j};
    if exist('segment_sel','var')
        clear segment_sel
    end
    switch subject_name
        case 'EF_r'
            real_subject_name = 'EF';
            path_folder_HD = [pathRepo_simcpspasticity,'\Alldata\'];
            load([path_folder_HD,'\IPSA\EF\001217b123_segments_deleted.mat']);
            switch joint
                case 'knee'
                    segment_sel = 1:8;
                case 'ankle'
                    segment_sel = 15:22;
            end               
            segment_sel_all = 1:22;
    end
   savefolder = [pathOpenSimModel,real_subject_name,...
        '\Spasticity\SpasticGainEstimation'];
    if exist('segment_sel','var') == 0
        continue
    end  
    Nsegment = length(segment_sel);              
    SelectRangeOnset = zeros(Nsegment,2);  
    load([savefolder,'\onset_man_user']);
    for ms = 1:Nsegment      
        % They changed the way ISA are structured, now GAS corresponds to 1
        % MEH corresponds to 2, REF corresponds to 3
        if strcmp(real_subject_name,'EF') && ...
                strcmp(allsegments(segment_sel(ms)).joint_name,'knee_ext')
            onset_idx = 'onset2';
        elseif strcmp(real_subject_name,'EF') && ...
                strcmp(allsegments(segment_sel(ms)).joint_name,'knee_flex')
            onset_idx = 'onset3';
        else
            onset_idx = 'onset1';
        end  
        
        trial = ['segment_' int2str(segment_sel(ms))];
        load([pathOpenSimModel,real_subject_name,...
            '\EMG\',['Stretch_',trial],'\emg1_norm_',optMuscles,'_ublM',...
            lMt_ub,'_ls_',ls,'_meth_',meth,'.mat']);
        load([path_folder_HD,'OpenSim\',real_subject_name,'\JointAngles\',...
            'IPSA\',['JointAngles_',trial,'_range.mat']]);
        if isnan(allsegments(segment_sel(ms)).(onset_idx))
            SelectRangeOnset(ms,:) = 100;
        % No low velocity trials since not used for estimation
        elseif (strcmp(allsegments(segment_sel(ms)).meas_type,'low') == 1)
            SelectRangeOnset(ms,:) = 100;
        % No medium velocity trials since not used for estimation
        elseif (strcmp(allsegments(segment_sel(ms)).meas_type,'medium') == 1)
            SelectRangeOnset(ms,:) = 100;
        else
        figure(segment_sel(ms))    
        plot(emg1processednorm(range,1),'linewidth',2); hold on;
        plot([onset_man(segment_sel(ms))-range(1)-1,...
            onset_man(segment_sel(ms))-range(1)-1],ylim,'k','linewidth',2);
        pl(1) = plot([onset_man(segment_sel(ms))-50-range(1)-1,...
            onset_man(segment_sel(ms))-50-range(1)-1],ylim,'k:',...
            'linewidth',2);
        plot([onset_man(segment_sel(ms))+120-range(1)-1,...
            onset_man(segment_sel(ms))+120-range(1)-1],ylim,'k:',...
            'linewidth',2);
        if exist([savefolder,'\SelectRangeOnset_',joint,'.mat'],'file')==2
            load([savefolder,'\SelectRangeOnset_',joint,'.mat']);
            pl(2) = plot([SelectRangeOnset3(segment_sel(ms),1),...
                SelectRangeOnset3(segment_sel(ms),1)],ylim,'b:',...
                'linewidth',2);
            plot([SelectRangeOnset3(segment_sel(ms),2),...
                SelectRangeOnset3(segment_sel(ms),2)],ylim,'b:',...
                'linewidth',2);
        end
        l = legend(pl,'Automatically selected range (-50, +120)',...
            'Already manually adjusted');
        title(['Select optimization range: first lower bound, second ',...
            'upper bound. If selected value < 0 then "automatically ',...
            'selected range" is selected'],'Fontsize',16);
        set(l,'Fontsize',16);
        ylim1 = ylim;
        ylim([-0.05 inf]);
        [startstop, yval]=ginput(2);
        if yval(1) < 0 
            % if this is already ok, then keep
            SelectRangeOnset(ms,1) = ...
                onset_man(segment_sel(ms))-50-range(1)-1;
        else
            % if this seems bad, then adjust
            SelectRangeOnset(ms,1) = startstop(1); 
        end

        if yval(2) < 0 
            % if this is already ok, then keep
            SelectRangeOnset(ms,2) = ...
                onset_man(segment_sel(ms))+120-range(1)-1;
        else
            % if this seems bad, then adjust
            SelectRangeOnset(ms,2) = startstop(2);
        end
        if SelectRangeOnset(ms,2)>length(emg1processednorm(range))
            SelectRangeOnset(ms,2) = length(emg1processednorm(range));
        end
        if SelectRangeOnset(ms,1)<1
            SelectRangeOnset(ms,1) = 1;
        end
        plot([SelectRangeOnset(ms,1) SelectRangeOnset(ms,1)],ylim1,'r',...
            'linewidth',1);
        plot([SelectRangeOnset(ms,2) SelectRangeOnset(ms,2)],ylim1,'r',...
            'linewidth',1);
        end    
    end

    SelectRangeOnset_temp = round(SelectRangeOnset);
    if ~(exist([savefolder,'\SelectRangeOnset_',joint,'_user.mat'],...
            'file')==2) 
        SelectRangeOnset3 = 10*ones(size(segment_sel_all,2),2);
    end
    SelectRangeOnset3(segment_sel,:) = SelectRangeOnset_temp;
    save([savefolder,'\SelectRangeOnset_',joint,'_user'],...
        'SelectRangeOnset3');
    end
end
