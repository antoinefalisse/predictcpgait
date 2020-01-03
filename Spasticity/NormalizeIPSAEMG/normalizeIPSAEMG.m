% This script processes the EMG from ISA: filtering and normalizing.
%
% Author: Antoine Falisse
% Date: 1/2/2020
%
clear all
close all
clc

%% User settings
saveNorm = 0;
saveFilt = 0;

%% Paths 
pathmain = pwd;
[pathSpasticity,~,~] = fileparts(pathmain);
addpath(genpath(pathSpasticity));
[pathRepo,~,~] = fileparts(pathSpasticity);
pathOpenSimModel = [pathRepo,'/OpenSimModel/'];

plotSignals = 0;
subjects_names = {'subject1'};
muscles = {'hamstrings','rectusfemoris','gastrocnemius'};

for sn = length(subjects_names)
    subject = subjects_names{sn};
    for mu = 1:length(muscles)
        muscle = muscles{mu};     
        switch subject
            case 'subject1'
                load([pathOpenSimModel,subject,'/IPSA/IPSA_data.mat']);
                Nslowsegment_ham = 4;   Nmedsegment_ham = 0;    Nfastsegment_ham = 4;   % hamstrings (1:8)  
                Nhamstrings = Nslowsegment_ham + Nmedsegment_ham + Nfastsegment_ham;
                Nslowsegment_rf = 2;    Nmedsegment_rf = 0;     Nfastsegment_rf = 4;    % rectus femoris (9:14)
                Nrectus = Nslowsegment_rf + Nmedsegment_rf + Nfastsegment_rf;
                Nslowsegment_gas = 4;   Nmedsegment_gas = 0;    Nfastsegment_gas = 4;   % gastrocnemius (15:22)
                Ngastroc = Nslowsegment_gas + Nmedsegment_gas + Nfastsegment_gas;
                side = 'right';
                offsetleft = 0;                      
        end        
        Misc.savefolder = [pathOpenSimModel,subject,'/EMG/IPSA'];       
        pathParameterEstimation = ...
            [pathOpenSimModel,'/',subject,'/ParameterEstimation/'];
        if strcmp(side,'right')
            side_lc = 'r';
            side_LC = 'R';
        else
            side_lc = 'l';
            side_LC = 'L';
        end        
        % Agonists
        switch muscle
            case 'hamstrings'
                 ms = 1+offsetleft:Nhamstrings+offsetleft;
            case 'rectusfemoris'
                ms = Nhamstrings+1+offsetleft:Nhamstrings+Nrectus+...
                    offsetleft;
            case 'gastrocnemius'
                ms = Nhamstrings+Nrectus+1+offsetleft:Nhamstrings+...
                    Nrectus+Ngastroc+offsetleft;
        end        
        for ms = ms   
            fs_mvt = allsegments(ms).info.movement_rate;
            % 1. High-pass, 20 Hz with dual-pass 6th order butterworth
            % filter. Already done => allsegments.emg1raw
            % 2a. Remove artifacts from 50Hz            
            if strcmp(subject,'LP_l')&&strcmp(muscle,'gastrocnemius')
                % Large artifact around 50Hz, need larger window
                [b,a] = butter(6,[48.5 51]/fs_mvt*2,'stop'); 
            else
                [b,a] = butter(6,[49.7 50.1]/fs_mvt*2,'stop');
            end
            % The structure of the data has changed with 8 emg cells
            % corresponding to different muscles
            if strcmp(subject,'subject1')
                switch muscle
                    case 'hamstrings'
                        emg1raw_s = filtfilt(b,a,allsegments(ms).emg2raw);                        
                    case 'rectusfemoris'
                        emg1raw_s = filtfilt(b,a,allsegments(ms).emg3raw); 
                    case 'gastrocnemius'
                        emg1raw_s = filtfilt(b,a,allsegments(ms).emg1raw); 
                end 
            else                   
                emg1raw_s = filtfilt(b,a,allsegments(ms).emg1raw);
            end
            % 2b. Remove artifacts from 46Hz
            [b,a] = butter(6,[45.95 46.2]/fs_mvt*2,'stop');
            emg1raw_ss = filtfilt(b,a,emg1raw_s);
            % 3. Rectification
            emg1_ss_hp_rec = abs(emg1raw_ss);
            % 4. Low-pass, 10 Hz with dual-pass 6th order butterworth 
            % filter
            [b,a] = butter(6,10/fs_mvt*2);
            emg1_ss_hp_rec_low6_10 = filtfilt(b,a,emg1_ss_hp_rec);          
            % 5. Normalize        
            emg1processed = emg1_ss_hp_rec_low6_10;
            % Load scaling factors  
            load([pathParameterEstimation,'ParameterEstimationResults.mat']);      
            switch muscle
                case 'hamstrings'
                    muscle_spa_names = {'bi_fem_lh_','semimem_','semiten_'}; 
                    indEMGISAinGait = find(strcmp(ParameterEstimationResults.Scale.EMGchannels,'MEH'));
                case 'rectusfemoris'
                    muscle_spa_names = {'rectus_fem_'};
                    indEMGISAinGait = find(strcmp(ParameterEstimationResults.Scale.EMGchannels,'REF'));
                case 'gastrocnemius'
                    muscle_spa_names = {'gas_med_','gas_lat_'};
                    indEMGISAinGait = find(strcmp(ParameterEstimationResults.Scale.EMGchannels,'GAS'));
            end    
            emg1processednorm = ...
                zeros(size(emg1processed,1),length(muscle_spa_names));
            for j = 1:length(muscle_spa_names)
                % we first need to normalize the EMGs to the max value
                % during gait as was done during the parameter optimization
                emg1processed_norm1 = emg1processed./ParameterEstimationResults.Scale.maxEMG_s.(side_lc)(indEMGISAinGait);
                % we then normalize
                normalizingvalue = ParameterEstimationResults.Scale.values_s.(side_lc)(strcmp(ParameterEstimationResults.Scale.muscles,muscle_spa_names{j}));                    
                emg1processednorm(:,j) = emg1processed_norm1.*normalizingvalue;                
            end                 
            if plotSignals   
                figure()
                subplot(4,1,1)
                plot(allsegments(ms).emg1raw)
                set(gca,'Fontsize',16);
                title('High-pass','Fontsize',16);
                subplot(4,1,2)
                plot(emg1_ss_hp_rec)
                set(gca,'Fontsize',16);
                title('Rectified','Fontsize',16);
                subplot(4,1,3)
                plot(emg1_ss_hp_rec_low6_10)
                set(gca,'Fontsize',16);
                title('Low-pass','Fontsize',16);
                subplot(4,1,4)
                for j = 1:length(muscle_spa_names)
                    plot(emg1processednorm(:,j)); hold on;
                end
                set(gca,'Fontsize',16);
                title('Normalized','Fontsize',16);
                sp = suptitle('Agonists');
                set(sp,'Fontsize',20);
            end
            trial = ['segment_' int2str(ms)];  
            path_filename = [Misc.savefolder,'/Stretch_',trial];
            if ~(exist(path_filename,'dir')==7)
                mkdir(path_filename);
            end
            if saveNorm
                save([path_filename,'/emg1_norm_personalized'],...
                    'emg1processednorm');
            end
            if saveFilt
                save([path_filename,'/emg1_filt'],'emg1processed'); 
            end
            clear emg1processednorm
            clear emg1processed
        end

        % Antagonists
        switch muscle
            case 'hamstrings'
                 ms = 1+offsetleft:Nhamstrings+offsetleft;
            case 'rectusfemoris'
                ms = Nhamstrings+1+offsetleft:Nhamstrings+Nrectus+...
                    offsetleft;
            case 'gastrocnemius'
                ms = Nhamstrings+Nrectus+1+offsetleft:Nhamstrings+...
                    Nrectus+Ngastroc+offsetleft;
        end
        close all
        for ms = ms 
            fs_mvt = allsegments(ms).info.movement_rate;
            % 1. High-pass, 20 Hz with dual-pass 6th order butterworth
            % filter. Already done => allsegments.emg2raw
            % 2a. Remove artificats from 50Hz
            if strcmp(subject,'TE_r')&&strcmp(muscle,'hamstrings')...
                    || strcmp(subject,'TE_l') && ...
                    strcmp(muscle,'hamstrings') || ...
                    strcmp(subject,'LP_r') && ...
                    strcmp(muscle,'hamstrings') || ...
                    strcmp(subject,'LP_r') && ...
                    strcmp(muscle,'gastrocnemius') || ...
                    strcmp(subject,'LP_l') && ...
                    strcmp(muscle,'gastrocnemius') || ...
                    strcmp(subject,'BS_postBTX') && ...
                    strcmp(muscle,'rectusfemoris') || ...
                    strcmp(subject,'LVV_postBTX') && ...
                    strcmp(muscle,'gastrocnemius') || ...
                    strcmp(subject,'LP_postBTX_L') && ...
                    strcmp(muscle,'gastrocnemius')
                % Large artifact around 50Hz, need larger window
                [b,a] = butter(6,[48.5 51]/fs_mvt*2,'stop'); 
            else
                [b,a] = butter(6,[49.7 50.1]/fs_mvt*2,'stop');
            end
            if strcmp(subject,'subject1')
                switch muscle
                    case 'hamstrings'
                        emg2raw_s = filtfilt(b,a,allsegments(ms).emg3raw);  % antagonist is rectus                      
                    case 'rectusfemoris'
                        emg2raw_s = filtfilt(b,a,allsegments(ms).emg2raw);  % antagonists are hamstrings.
                    case 'gastrocnemius'
                        emg2raw_s = filtfilt(b,a,allsegments(ms).emg5raw);  % antagonist is tibialis anterior
                end 
            else
                emg2raw_s = filtfilt(b,a,allsegments(ms).emg2raw);
            end
            % 2b. Remove artifacts from 46Hz
            [b,a] = butter(6,[45.95 46.2]/fs_mvt*2,'stop');
            emg2raw_ss = filtfilt(b,a,emg2raw_s);
            % 3. Rectification
            emg2_ss_hp_rec = abs(emg2raw_ss);
            % 4. Low-pass, 10 Hz with dual-pass 6th order butterworth 
            % filter
            [b,a] = butter(6,10/fs_mvt*2);
            emg2_ss_hp_rec_low6_10 = filtfilt(b,a,emg2_ss_hp_rec);
            % 5. Normalize        
            emg2processed = emg2_ss_hp_rec_low6_10;
            switch muscle
                case 'hamstrings' % antagonist is rectus
                    muscle_spa_names = {'rectus_fem_'};
                    indEMGISAinGait = find(strcmp(ParameterEstimationResults.Scale.EMGchannels,'REF'));
                case 'rectusfemoris' % antagonists are hamstrings.
                    muscle_spa_names = {'bi_fem_lh_','semimem_','semiten_'};
                    indEMGISAinGait = find(strcmp(ParameterEstimationResults.Scale.EMGchannels,'MEH'));
                case 'gastrocnemius' % antagonist is tibialis anterior
                    muscle_spa_names = {'tib_ant_'};
                    indEMGISAinGait = find(strcmp(ParameterEstimationResults.Scale.EMGchannels,'GAS'));
            end
            emg2processednorm = zeros(size(emg2processed,1),...
                length(muscle_spa_names));                
            for j = 1:length(muscle_spa_names)
                % we first need to normalize the EMGs to the max value
                % during gait as was done during the parameter optimization
                emg2processed_norm1 = emg2processed./ParameterEstimationResults.Scale.maxEMG_s.(side_lc)(indEMGISAinGait);
                % we then normalize
                normalizingvalue = ParameterEstimationResults.Scale.values_s.(side_lc)(strcmp(ParameterEstimationResults.Scale.muscles,muscle_spa_names{j}));                    
                emg2processednorm(:,j) = emg2processed_norm1.*normalizingvalue;     
            end
            if plotSignals    
                figure()
                subplot(4,1,1)
                plot(allsegments(ms).emg2raw)
                set(gca,'Fontsize',16);
                title('High-pass','Fontsize',16);
                subplot(4,1,2)
                plot(emg2_ss_hp_rec)
                set(gca,'Fontsize',16);
                title('Rectified','Fontsize',16);
                subplot(4,1,3)
                plot(emg2_ss_hp_rec_low6_10)
                set(gca,'Fontsize',16);
                title('Low-pass','Fontsize',16);
                subplot(4,1,4)
                for j = 1:length(muscle_spa_names)
                    plot(emg2processednorm(:,j)); hold on;
                end
                set(gca,'Fontsize',16);
                title('Normalized','Fontsize',16);
                sp = suptitle('Antagonists');
                set(sp,'Fontsize',20);
            end
            trial = ['segment_' int2str(ms)];            
            path_filename = [Misc.savefolder,'/Stretch_',trial]; 
            if ~(exist(path_filename,'dir')==7)
                mkdir(path_filename);
            end
            if saveNorm
                save([path_filename,'/emg2_norm_personalized'],...
                    'emg2processednorm');
            end
            if saveFilt
                save([path_filename,'/emg2_filt'],'emg2processed'); 
            end
            clear emg2processednorm
            clear emg2processed
        end
    end
end
