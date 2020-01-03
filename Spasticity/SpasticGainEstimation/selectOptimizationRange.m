% This scripts allows users to manually adjust the optimization time interval

clc
close all
clear all

%% Paths 
pathmain = pwd;
[pathSpasticity,~,~] = fileparts(pathmain);
addpath(genpath(pathSpasticity));
[pathRepo,~,~] = fileparts(pathSpasticity);
pathOpenSimModel = [pathRepo,'/OpenSimModel/'];

subjects_names = {'subject1'};
joints = {'knee','ankle'};
for sn = 1:length(subjects_names)     
    subject_name = subjects_names{sn};
    for j = 1:length(joints)
    joint = joints{j};
    if exist('segment_sel','var')
        clear segment_sel
    end
    switch subject_name
        case 'subject1'
            subject = subject_name;
            load([pathOpenSimModel,subject,'/IPSA/IPSA_data.mat']);
            switch joint
                case 'knee'
                    segment_sel = 6:8;
                case 'ankle'
                    segment_sel = 19:22;
            end               
            segment_sel_all = 1:22;
    end
   savefolder = [pathOpenSimModel,subject,'/Spasticity/SpasticGainEstimation'];
    if exist('segment_sel','var') == 0
        continue
    end  
    Nsegment = length(segment_sel);              
    SelectRangeOnset = zeros(Nsegment,2);  
    load([savefolder,'/onset_man_user']);
    for ms = 1:Nsegment      
        % They changed the way ISA are structured, now GAS corresponds to 1
        % MEH corresponds to 2, REF corresponds to 3
        if strcmp(subject,'subject1') && ...
                strcmp(allsegments(segment_sel(ms)).joint_name,'knee_ext')
            onset_idx = 'onset2';
        elseif strcmp(subject,'subject1') && ...
                strcmp(allsegments(segment_sel(ms)).joint_name,'knee_flex')
            onset_idx = 'onset3';
        else
            onset_idx = 'onset1';
        end  
        
        trial = ['segment_' int2str(segment_sel(ms))];
        load([pathOpenSimModel,subject,...
            '/EMG/IPSA/',['Stretch_',trial],'/emg1_norm_personalized.mat']);
        load([pathOpenSimModel,subject,...
            '/IK/','IPSA/',['JointAngles_',trial,'_range.mat']]);
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
        if exist([savefolder,'/SelectRangeOnset_',joint,'.mat'],'file')==2
            load([savefolder,'/SelectRangeOnset_',joint,'.mat']);
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
    if ~(exist([savefolder,'/SelectRangeOnset_',joint,'_user.mat'],...
            'file')==2) 
        SelectRangeOnset3 = 10*ones(size(segment_sel_all,2),2);
    end
    SelectRangeOnset3(segment_sel,:) = SelectRangeOnset_temp;
    save([savefolder,'/SelectRangeOnset_',joint,'_user'],...
        'SelectRangeOnset3');
    end
end
