% This script plots EMG against simulated muscle excitations resulting drom
% the spastic models using the optimized feedback gains
clear all
close all
clc

% lMtilde_ext.max = 1.5;
% optMuscles = 'out_0';
% numbTrials = '4';
% NISA = '6';
% dev_p = 5;
% lMt_ub = num2str(10*lMtilde_ext.max);
% W.act         = 0.00100;
% W.aT        = 0.60000;
% W.vA        = 0.00001;
% W.dF        = 0.00001;
% W.EMG.all   = 0.003500-W.vA-W.dF;
% W.lMopt_max = 0.394500;
% W.pen       = 0.00100;
% nameParameters = [optMuscles,'_ublM',lMt_ub,...
%     '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
%     '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
%     '_EMG',num2str(round(W.EMG.all*10000)),...
%     '_lM',num2str(round(W.lMopt_max*10000)),...
%     '_pen',num2str(round(W.pen*10000)),'_2sides'];
testCase = '32';
Misc.bspas = 100;
lMtilde_ext.max = 1.5;
optMuscles = 'out_0';
numbTrials = '4';
NISA = '6';
dev_p = 5;
lMt_ub = num2str(10*lMtilde_ext.max);

switch testCase
    case '32'
        ScaleMIN = 0.01;
        delta = 0.1;
        W.act       = 0.00100;
        W.aT        = 0.60000;
        W.vA        = 0.00001;
        W.dF        = 0.00001;
        W.EMG.all   = 0.003500-W.vA-W.dF;
        W.lMopt_max = 0.394500;
        W.lMopt     = 0.00100;
        W.act       = 0.00100-W.act*delta/(1-0.00100);
        W.aT        = 0.60000-W.aT*delta/(1-0.00100);
        W.vA        = 0.00001-W.vA*delta/(1-0.00100);
        W.dF        = 0.00001-W.dF*delta/(1-0.00100);
        W.EMG.all   = 0.003480-W.EMG.all*delta/(1-0.00100);
        W.lMopt_max = 0.394500-W.lMopt_max*delta/(1-0.00100);
        W.lMopt     = 0.00100+delta;    
    case '34'
        ScaleMIN = 0.3;
        delta = 0.1;
        W.act       = 0.00100;
        W.aT        = 0.60000;
        W.vA        = 0.00001;
        W.dF        = 0.00001;
        W.EMG.all   = 0.003500-W.vA-W.dF;
        W.lMopt_max = 0.394500;
        W.lMopt     = 0.00100;
        W.act       = 0.00100-W.act*delta/(1-0.00100);
        W.aT        = 0.60000-W.aT*delta/(1-0.00100);
        W.vA        = 0.00001-W.vA*delta/(1-0.00100);
        W.dF        = 0.00001-W.dF*delta/(1-0.00100);
        W.EMG.all   = 0.003480-W.EMG.all*delta/(1-0.00100);
        W.lMopt_max = 0.394500-W.lMopt_max*delta/(1-0.00100);
        W.lMopt     = 0.00100+delta;        
end
nameParameters = [optMuscles,'_ublM',lMt_ub,...
    '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
    '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
    '_EMG',num2str(round(W.EMG.all*10000)),...
    '_lM',num2str(round(W.lMopt_max*10000)),...
    '_lMopt',num2str(round(W.lMopt*10000)),'_2sides',...
    '_scMin',num2str(ScaleMIN*100)];

% Add folder to MATLAB search path 
pathmain = pwd;
[pathRepo,~,~] = fileparts(pathmain);
addpath(genpath(pathRepo));
pathOpenSimModel = [pathRepo,'\OpenSimModel\'];
pathRepo_simcpspasticity = 'C:\Users\u0101727\Documents\MyRepositories\simcpspasticity_cases';


% User setting
savefits = 0; % 1 to save the fitting results (RMSE and R²), 0 otherwise

% Three spasticity models are proposed: force-related model (muscle-tendon
% force and dF/dt feedback), velocity-related model (muscle fiber length
% and velocity feedback), and acceleration-related model (muscle fiber
% length, velocity, and acceleration feedback)
formulations = {'forceModel'};

for ff = 1:length(formulations)
    formulation = formulations{ff};
for iii = 2%:22%:20
% Different cases
cond = num2str(iii);
switch cond    
    case '1'
        subject_name = 'EF_r';
        side = 'r';
        joint = 'knee';
        segment_sel = [6 7 8];
    case '2'
        subject_name = 'EF_r';
        side = 'r';
        joint = 'ankle';
        segment_sel = [19 20 21];
    case '3'
        subject_name = 'EF_r';
        side = 'r';
        joint = 'ankle';
        segment_sel = [19 20 22];
    case '4'
        subject_name = 'EF_r';
        side = 'r';
        joint = 'ankle';
        segment_sel = [19 21 22];
    case '5'
        subject_name = 'EF_r';
        side = 'r';
        joint = 'ankle';
        segment_sel = [20 21 22];
        
end

switch subject_name
    case 'EF_r'
        real_subject_name = 'EF';
        path_folder_HD = [pathRepo_simcpspasticity,'\Alldata\'];
end

% Load data from spastic gain estimation (sge)
number_save = [formulation '_' int2str(segment_sel)];
savefolder_sge = [pathOpenSimModel,real_subject_name,...
    '\Spasticity\SpasticGainEstimation'];
savefolder_ForwardSimulation = [pathOpenSimModel,real_subject_name,...
    '\Spasticity\ForwardSimulations\IPSA_',nameParameters];
load([savefolder_sge,'\',formulation,'_',nameParameters,'_',num2str(Misc.bspas),...
    '\output_',joint,'_',number_save]);
load([savefolder_sge,'\',formulation,'_',nameParameters,'_',num2str(Misc.bspas),...
    '\DatStore_',joint,'_',number_save]);
% Load data from forward simulations
load([savefolder_ForwardSimulation,'\FtildeGPOPSFTtilde_',joint]);
load([savefolder_ForwardSimulation,'\dFtildeGPOPSFTtilde_',joint]);
load([savefolder_ForwardSimulation,'\lMtildeGPOPSFTtilde_',joint]);
load([savefolder_ForwardSimulation,'\vMtildeGPOPSFTtilde_',joint]);
load([savefolder_ForwardSimulation,'\outputGPOPSFTtilde_',joint]);

switch formulation
    %% Velocity model
    case 'velocityModel'
        % Extract data
        auxdata         = output.result.setup.auxdata;
        NMuscles        = auxdata.NMuscles; 
        NMuscles_spas   = auxdata.NMuscles_spas;
        Muscles_spas    = auxdata.Muscles_spas;
        Nsegment        = auxdata.Nsegment;
        joint           = auxdata.joint;   
        % Pre-allocation
        time = struct('MS',[]);
        elMf = struct('MS',[]);
        evMf = struct('MS',[]);
        glMf = struct('MS',[]);
        gvMf = struct('MS',[]);
        bc = struct('MS',[]);
        e = struct('MS',[]);
        EMGSpline = struct('MS',[]);
        % Results from optimal control problem
        for ms = 1:Nsegment
            % Get time
            time(ms).MS = output.result.solution.phase(ms).time;
            numcolpoints = size(time(ms).MS,1);
            % Get states
            elMf(ms).MS = output.result.solution.phase(ms).state(...
                :,1:NMuscles_spas);
            evMf(ms).MS = output.result.solution.phase(ms).state(...
                :,NMuscles_spas+1:2*NMuscles_spas);                 
            % Get parameters
            glMf(ms).MS = output.result.solution.parameter(...
                :,1:NMuscles_spas);
            gvMf(ms).MS = output.result.solution.parameter(...
                :,NMuscles_spas+1:2*NMuscles_spas);
            % Total excitation is the combination of baseline activation 
            % and muscle excitation from the spasticity/feedback models
            bc(ms).MS = repmat(auxdata.min_EMG_ordered_beg_mat(ms,:),...
                numcolpoints,1);                  
            e(ms).MS = bc(ms).MS + elMf(ms).MS + evMf(ms).MS;
            % Get splines
            for m = 1:NMuscles_spas                 
                EMGSpline(ms).MS(:,m) = ppval(...
                    auxdata.EMGSpline(ms).MS(m),time(ms).MS);
            end            
        end
        
        % Plot results
        col = hsv(NMuscles_spas);
        switch joint
            case 'knee'
                name_muscles_spastic = {'Biceps femoris lh',...
                    'Semimembranosus','Semitendinosus'};
            case 'ankle'
                name_muscles_spastic = {'Gas med','Gas lat'};
        end    
        % Optimization time window
        figure()
        for ms = 1:Nsegment            
            for m = 1:NMuscles_spas
                subplot(Nsegment,NMuscles_spas,(ms-1)*NMuscles_spas+m)
                plot(time(ms).MS,EMGSpline(ms).MS(:,m),'k','LineWidth',3); 
                hold on;
                plot(time(ms).MS,evMf(ms).MS(:,m),'b','LineWidth',2);
                plot(time(ms).MS,elMf(ms).MS(:,m),'r','LineWidth',2);
                plot(time(ms).MS,bc(ms).MS(:,m),'g','LineWidth',2);
                plot(time(ms).MS,e(ms).MS(:,m),'m','LineWidth',2);             
                if ms == Nsegment
                    xlabel('Time (s)','Fontsize',16);
                end        
                if m == 1
                    ylabel('m. exc. ()','Fontsize',16);
                end        
                if ms == 1
                    title(name_muscles_spastic{m},'Fontsize',16);
                end
                ylim([0 inf]);
                set(gca,'Fontsize',16);
            end
        end
        legendstr = {'EMG','vM f.','lM f.','bc','all f.'};
        l = legend(legendstr);
        set(l,'Fontsize',16);
        sp = suptitle(['Comparison EMG with sensory feedback: velocity',...
            'model (optimization time window)']);
        set(sp,'Fontsize',20);
        
        % Full time range      
        ulMf = struct('MS',[]);
        uvMf = struct('MS',[]);
        ulMfinterp = struct('MS',[]);
        uvMfinterp = struct('MS',[]);
        EMG_path_agonist = struct('MS',[]);
        EMG_path_antagonist = struct('MS',[]);
        emg1processedstr = struct('MS',[]);
        emg2processedstr = struct('MS',[]);
        EMG_ordered = struct('MS',[]);
        bc_adj = struct('MS',[]);
        fit_IPSA = struct('r2',[]);
        figure()
        for ms = 1:Nsegment            
            for m = 1:NMuscles_spas                
                time_ref = round(outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time(1):...
                    0.005:outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time(end),4);
                time_temp = round(outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time,4);  
                % Here we chop off beginning and end of time range
                time_sel = time_ref(41:end-21);   
                % Forward integration              
                ulMf(ms).MS(:,m) = lengthFDynamicsExplicit(...
                    outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time,...
                    lMtildeGPOPS.(['segment_',num2str(segment_sel(ms))])...
                    (:,auxdata.Muscles_spas(m)),auxdata.lMf_td,...
                    glMf(ms).MS(m),100,auxdata.threshold_lMf(ms,m));
                uvMf(ms).MS(:,m) = velocityFDynamicsExplicit(...
                    outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time,...
                    vMtildeGPOPS.(['segment_',num2str(segment_sel(ms))])...
                    (:,auxdata.Muscles_spas(m)),auxdata.vMf_td,...
                    gvMf(ms).MS(m),100,auxdata.threshold_vMf(ms,m));
                % Interpolation                    
                ulMfinterp(ms).MS(:,m) = interp1(...
                    time_temp,ulMf(ms).MS(:,m),time_sel');
                uvMfinterp(ms).MS(:,m) = interp1(...
                    time_temp,uvMf(ms).MS(:,m),time_sel');
                % Load EMG data                
                trial = ['segment_' int2str(segment_sel(ms))];  
                if strcmp(subject_name,'EF_r') || strcmp(subject_name,'EF_l')
                    EMG_path_agonist(ms).MS = fullfile(path_folder_HD,...
                    'OptimizedData',real_subject_name,'EMGProcessing',...
                    ['Stretch_',trial],...
                    'emg1processednorm_kneemuscles.mat');                
                    EMG_path_antagonist(ms).MS =  fullfile(path_folder_HD,...
                    'OptimizedData',real_subject_name,'EMGProcessing',...
                    ['Stretch_',trial],...
                    'emg2processednorm_kneemuscles.mat'); 
                 % Load motion range
                    Misc.range(ms).MS = load([path_folder_HD,'OpenSim\',...
                        real_subject_name,'\JointAngles\','IPSA\',...
                        ['JointAngles_',trial,'_range.mat']]);
                else
                    EMG_path_agonist(ms).MS = fullfile(path_folder_HD,...
                        'EMG','ProcessedStretch','MS',['Stretch_',trial],...
                        'emg1processednorm_kneemuscles.mat');
                    EMG_path_antagonist(ms).MS = fullfile(path_folder_HD,...
                        'EMG','ProcessedStretch','MS',['Stretch_',trial],...
                        'emg2processednorm_kneemuscles.mat');                    
                    % Load motion range
                    Misc.range(ms).MS = load([path_folder_HD '\OpenSim\',...
                        'SymModel\JointAngles\DataPellenberg\JointAngles_',...
                        trial,'_range.mat']);
                end
                load(EMG_path_agonist(ms).MS);
                load(EMG_path_antagonist(ms).MS);
                emg1processedstr(ms).MS = ...
                    emg1processednorm(Misc.range(ms).MS.range,:);
                emg2processedstr(ms).MS = ...
                    emg2processednorm(Misc.range(ms).MS.range,:);
                % Re-order EMG
                EMG_ordered(ms).MS = ...
                    zeros(size(emg1processedstr(ms).MS,1),NMuscles);   
                EMG_ordered(ms).MS(:,auxdata.Muscles_spas) = ...
                    emg1processedstr(ms).MS;
                EMG_ordered(ms).MS(:,auxdata.Muscles_antagonist) = ...
                    emg2processedstr(ms).MS;
                % Adjust EMG based on choped off bits during spastic 
                % gain estimation and here above
                EMG_ordered(ms).MS = EMG_ordered(ms).MS(10:end-10,:);
                EMG_ordered(ms).MS = EMG_ordered(ms).MS(41:end-21,:);               
                % Baseline activation
                bc_adj(ms).MS = repmat(bc(ms).MS(1,:),size(time_sel,2),1); 
                % Compare measured EMG with simulated EMG (ie, muscle
                % activity resulting from the spasticy model + baseline
                % activation)
                EMG_meas = EMG_ordered(ms).MS(:,auxdata.Muscles_spas(m));
                EMG_sim = uvMfinterp(ms).MS(:,m)+ulMfinterp(ms).MS(:,m)+...
                    bc_adj(ms).MS(:,m);
                % Compute RMSE and R²
                [fit_IPSA.r2(ms,m),fit_IPSA.rmse(ms,m)] = ...
                    rsquare(EMG_meas,EMG_sim);
                subplot(Nsegment,NMuscles_spas,(ms-1)*NMuscles_spas+m)
                plot(time_sel, EMG_meas,'k','LineWidth',3); hold on;
                plot(time_sel, uvMfinterp(ms).MS(:,m),'b','LineWidth',2);
                plot(time_sel, ulMfinterp(ms).MS(:,m),'r','LineWidth',2);
                plot(time_sel, bc_adj(ms).MS(:,m),'g','LineWidth', 2);
                plot(time_sel, EMG_sim,'m','LineWidth',2);                
                if ms == Nsegment
                    xlabel('Time (s)','Fontsize',16);
                end        
                if m == 1
                    ylabel('m. exc. ()','Fontsize',16);
                end        
                if ms == 1
                    title(name_muscles_spastic{m},'Fontsize',16);
                end            
                ylim([0 inf]);
                set(gca,'Fontsize',16);  
            end
        end
        legendstr = {'EMG','vM f.','lM f.','bc','all f.'};
        l = legend(legendstr);
        set(l,'Fontsize',16);         
        sp = suptitle(['Comparison EMG with sensory feedback: velocity'...
            'model (full time interval)']);     
        set(sp,'Fontsize',20);
        % Save fit for statistical analysis
        if savefits
            fitpath = [savefolder_sge,'\',formulation,'\fit_IPSA'];
            save([fitpath,'_joint_',number_save],'fit_IPSA');
        end     

    %% Force model
    case 'forceModel'
        % Extract data
        auxdata         = output.result.setup.auxdata;
        NMuscles        = auxdata.NMuscles; 
        NMuscles_spas   = auxdata.NMuscles_spas;
        Muscles_spas    = auxdata.Muscles_spas;
        Nsegment        = auxdata.Nsegment;
        joint           = auxdata.joint;       
        % Pre-allocation
        time = struct('MS',[]);
        eFf = struct('MS',[]);
        edFf = struct('MS',[]);
        gFf = struct('MS',[]);
        gdFf = struct('MS',[]);
        bc = struct('MS',[]);
        e = struct('MS',[]);
        EMGSpline = struct('MS',[]);
        % Results from optimal control problem
        for ms = 1:Nsegment
            % Get time
            time(ms).MS = output.result.solution.phase(ms).time;
            numcolpoints = size(time(ms).MS,1);
            % Get states
            eFf(ms).MS = output.result.solution.phase(ms).state(...
                :,1:NMuscles_spas);
            edFf(ms).MS = output.result.solution.phase(ms).state(...
                :,NMuscles_spas+1:2*NMuscles_spas);                 
            % Get parameters
            gFf(ms).MS = output.result.solution.parameter(...
                :,1:NMuscles_spas);
            gdFf(ms).MS = output.result.solution.parameter(...
                :,NMuscles_spas+1:2*NMuscles_spas);
            % Total excitation is the combination of baseline activation 
            % and muscle excitation from the spasticity/feedback models
            bc(ms).MS = repmat(auxdata.min_EMG_ordered_beg_mat(ms,:),...
                numcolpoints,1);          
            e(ms).MS =  eFf(ms).MS + edFf(ms).MS + bc(ms).MS;
            % Get splines
            for m = 1:NMuscles_spas                 
                EMGSpline(ms).MS(:,m) = ppval(...
                    auxdata.EMGSpline(ms).MS(m),time(ms).MS);
            end            
        end

        % Plot results
        col = hsv(NMuscles_spas);
        switch joint
            case 'knee'
                name_muscles_spastic = {'Biceps femoris lh',...
                    'Semimembranosus','Semitendinosus'};
            case 'ankle'
                name_muscles_spastic = {'Gas med','Gas lat'};
        end 
        
        % Optimization time window
        figure()
        for ms = 1:Nsegment            
            for m = 1:NMuscles_spas
                subplot(Nsegment,NMuscles_spas,(ms-1)*NMuscles_spas+m)
                plot(time(ms).MS, EMGSpline(ms).MS(:,m),'k','LineWidth',3); 
                hold on;
                plot(time(ms).MS, edFf(ms).MS(:,m),'b','LineWidth',2);
                plot(time(ms).MS, eFf(ms).MS(:,m),'r','LineWidth',2);
                plot(time(ms).MS, bc(ms).MS(:,m),'g','LineWidth',2);
                plot(time(ms).MS, e(ms).MS(:,m),'m','LineWidth',2);             
                if ms == Nsegment
                    xlabel('Time (s)','Fontsize',16);
                end        
                if m == 1
                    ylabel('m. exc. ()','Fontsize',16);
                end        
                if ms == 1
                    title(name_muscles_spastic{m},'Fontsize',16);
                end
                ylim([0 inf]);
                set(gca,'Fontsize',16);
            end
        end
        legendstr = {'EMG','dF/dt f.','Force f.','bc','all f.'};
        l = legend(legendstr);
        set(l,'Fontsize',16);
        sp = suptitle(['Comparison EMG with sensory feedback: force',...
            'model (optimization time window)']);
        set(sp,'Fontsize',20);
        
        % Full time range
        uFf = struct('MS',[]);
        udFf = struct('MS',[]);
        uFfinterp = struct('MS',[]);
        udFfinterp = struct('MS',[]);
        EMG_path_agonist = struct('MS',[]);
        EMG_path_antagonist = struct('MS',[]);
        emg1processedstr = struct('MS',[]);
        emg2processedstr = struct('MS',[]);
        EMG_ordered = struct('MS',[]);
        bc_adj = struct('MS',[]);
        fit_IPSA = struct('r2',[]);
        figure()
        for ms = 1:Nsegment            
            for m = 1:NMuscles_spas                
                time_ref = round(outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time(1):...
                    0.005:outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time(end),4);
                time_temp = round(outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time,4);    
                % Here we chop off beginning and end of time range
                time_sel = time_ref(41:end-21);
                % Forward integration   
                uFf(ms).MS(:,m) = forceFDynamicsExplicit(...
                    outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time,...
                    FtildeGPOPS.(['segment_',num2str(segment_sel(ms))])...
                    (:,auxdata.Muscles_spas(m)),auxdata.Ff_td,...
                    gFf(ms).MS(m),100,auxdata.threshold_Ff(ms,m));
                udFf(ms).MS(:,m) = dFdtFDynamicsExplicit(...
                    outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time,...
                    dFtildeGPOPS.(['segment_',num2str(segment_sel(ms))])...
                    (:,auxdata.Muscles_spas(m)),auxdata.dFf_td,...
                    gdFf(ms).MS(m),100,auxdata.threshold_dFf(ms,m));
                % Interpolation
                uFfinterp(ms).MS(:,m) = interp1(...
                    time_temp,uFf(ms).MS(:,m),time_sel');
                udFfinterp(ms).MS(:,m) = interp1(...
                    time_temp,udFf(ms).MS(:,m),time_sel');
                % Load EMG data                
                trial = ['segment_' int2str(segment_sel(ms))];  
                if strcmp(subject_name,'EF_r') || strcmp(subject_name,'EF_l')
                    EMG_path_agonist(ms).MS = [pathOpenSimModel,real_subject_name,...
                        '\EMG\',['Stretch_',trial],'\emg1_norm_',nameParameters,'.mat'];                
                    EMG_path_antagonist(ms).MS = [pathOpenSimModel,real_subject_name,...
                        '\EMG\',['Stretch_',trial],'\emg2_norm_',nameParameters,'.mat']; 
                 % Load motion range
                    Misc.range(ms).MS = load([path_folder_HD,'OpenSim\',...
                        real_subject_name,'\JointAngles\','IPSA\',...
                        ['JointAngles_',trial,'_range.mat']]);
                else
                    EMG_path_agonist(ms).MS = fullfile(path_folder_HD,...
                        'EMG','ProcessedStretch','MS',['Stretch_',trial],...
                        'emg1processednorm_kneemuscles.mat');
                    EMG_path_antagonist(ms).MS = fullfile(path_folder_HD,...
                        'EMG','ProcessedStretch','MS',['Stretch_',trial],...
                        'emg2processednorm_kneemuscles.mat');                    
                    % Load motion range
                    Misc.range(ms).MS = load([path_folder_HD '\OpenSim\',...
                        'SymModel\JointAngles\DataPellenberg\JointAngles_',...
                        trial,'_range.mat']);
                end
                load(EMG_path_agonist(ms).MS);
                load(EMG_path_antagonist(ms).MS);
                emg1processedstr(ms).MS = ...
                    emg1processednorm(Misc.range(ms).MS.range,:);
                emg2processedstr(ms).MS = ...
                    emg2processednorm(Misc.range(ms).MS.range,:);
                % Re-order EMG
                EMG_ordered(ms).MS = ...
                    zeros(size(emg1processedstr(ms).MS,1),NMuscles);   
                EMG_ordered(ms).MS(:,auxdata.Muscles_spas) = ...
                    emg1processedstr(ms).MS;
                EMG_ordered(ms).MS(:,auxdata.Muscles_antagonist) = ...
                    emg2processedstr(ms).MS;
                % Adjust EMG based on choped off bits during spastic 
                % gain estimation and here above
                EMG_ordered(ms).MS = EMG_ordered(ms).MS(10:end-10,:); 
                EMG_ordered(ms).MS = EMG_ordered(ms).MS(41:end-21,:); 
                % Baseline activation
                bc_adj(ms).MS = repmat(bc(ms).MS(1,:),size(time_sel,2),1); 
                % Compare measured EMG with simulated EMG (ie, muscle
                % activity resulting from the spasticy model + baseline
                % activation)
                EMG_meas = EMG_ordered(ms).MS(:,auxdata.Muscles_spas(m));
                EMG_sim = udFfinterp(ms).MS(:,m)+uFfinterp(ms).MS(:,m)+...
                    bc_adj(ms).MS(:,m);
                % Compute RMSE and R²
                [fit_IPSA.r2(ms,m),fit_IPSA.rmse(ms,m)] = ...
                    rsquare(EMG_meas,EMG_sim);                
                subplot(Nsegment,NMuscles_spas,(ms-1)*NMuscles_spas+m)
                plot(time_sel, EMG_meas,'k','LineWidth',3); hold on;
                plot(time_sel, udFfinterp(ms).MS(:,m),'b','LineWidth',2);
                plot(time_sel, uFfinterp(ms).MS(:,m),'r','LineWidth',2);
                plot(time_sel, bc_adj(ms).MS(:,m),'g','LineWidth',2);
                plot(time_sel, EMG_sim,'m','LineWidth',2);                
                if ms == Nsegment
                    xlabel('Time (s)','Fontsize',16);
                end        
                if m == 1
                    ylabel('m. exc. ()','Fontsize',16);
                end        
                if ms == 1
                    title(name_muscles_spastic{m},'Fontsize',16);
                end
                ylim([0 inf]);
                set(gca,'Fontsize',16);               
            end
        end
        legendstr = {'EMG','dF/dt f.','Force f.','bc','all f.'};
        l = legend(legendstr);
        set(l,'Fontsize',16);         
        sp = suptitle(['Comparison EMG with sensory feedback: force'...
            'model (full time interval)']);      
        set(sp,'Fontsize',20);
        % Save fit for statistical analysis
        if savefits
            fitpath = [savefolder_sge,'\',formulation,'_',nameParameters,'\fit_IPSA'];
            save([fitpath,'_joint_',number_save],'fit_IPSA');
        end     
        
    %% Acceleration model
    case 'accelerationModel'
        % Extract data
        auxdata         = output.result.setup.auxdata;
        NMuscles        = auxdata.NMuscles; 
        NMuscles_spas   = auxdata.NMuscles_spas;
        Muscles_spas    = auxdata.Muscles_spas;
        Nsegment        = auxdata.Nsegment;
        joint           = auxdata.joint;   
        % Pre-allocation
        time = struct('MS',[]);
        elf = struct('MS',[]);
        evf = struct('MS',[]);
        eaf = struct('MS',[]);
        glf = struct('MS',[]);
        gvf = struct('MS',[]);
        gaf = struct('MS',[]);
        bc = struct('MS',[]);
        e = struct('MS',[]);
        EMGSpline = struct('MS',[]);
        for ms = 1:Nsegment
            % Get time
            time(ms).MS = output.result.solution.phase(ms).time;
            numcolpoints = size(time(ms).MS,1);
            % Get states
            elf(ms).MS = output.result.solution.phase(ms).state(...
                :,1:NMuscles_spas);
            evf(ms).MS = output.result.solution.phase(ms).state(...
                :,NMuscles_spas+1:2*NMuscles_spas);    
            eaf(ms).MS = output.result.solution.phase(ms).state(...
                :,2*NMuscles_spas+1:3*NMuscles_spas); 
            % Get parameters
            glf(ms).MS = output.result.solution.parameter(...
                :,1:NMuscles_spas);
            gvf(ms).MS = output.result.solution.parameter(...
                :,NMuscles_spas+1:2*NMuscles_spas);
            gaf(ms).MS = output.result.solution.parameter(...
                :,2*NMuscles_spas+1:3*NMuscles_spas);
            % Total excitation is the combination of baseline activation 
            % and muscle excitation from the spasticity/feedback models
            bc(ms).MS = repmat(auxdata.min_EMG_ordered_beg_mat(ms,:),...
                numcolpoints,1);        
            e(ms).MS = elf(ms).MS + evf(ms).MS + eaf(ms).MS + bc(ms).MS;
            % Get splines
            for m = 1:NMuscles_spas                 
                EMGSpline(ms).MS(:,m) = ...
                    ppval(auxdata.EMGSpline(ms).MS(m),time(ms).MS);
            end            
        end

        % Plot results
        col = hsv(NMuscles_spas);
        switch joint
            case 'knee'
                name_muscles_spastic = {'Biceps femoris lh',...
                    'Semimembranosus','Semitendinosus'};
            case 'ankle'
                name_muscles_spastic = {'Gas med','Gas lat'};
        end 
        
        % Optimization time window
        figure()
        for ms = 1:Nsegment            
            for m = 1:NMuscles_spas
                subplot(Nsegment,NMuscles_spas,(ms-1)*NMuscles_spas+m)
                plot(time(ms).MS, EMGSpline(ms).MS(:,m),'k','LineWidth',3);
                hold on;
                plot(time(ms).MS, eaf(ms).MS(:,m),'b','LineWidth',2);
                plot(time(ms).MS, evf(ms).MS(:,m),'r','LineWidth',2);
                plot(time(ms).MS, elf(ms).MS(:,m),'y','LineWidth',2);
                plot(time(ms).MS, bc(ms).MS(:,m),'g','LineWidth',2);
                plot(time(ms).MS, e(ms).MS(:,m),'m','LineWidth',2);              
                if ms == Nsegment
                    xlabel('Time (s)','Fontsize',16);
                end        
                if m == 1
                    ylabel('m. exc. ()','Fontsize',16);
                end        
                if ms == 1
                    title(name_muscles_spastic{m},'Fontsize',16);
                end
                ylim([0 inf]);
                set(gca,'Fontsize',16);
            end
        end
        legendstr = {'EMG','aM f.','vM f.','lM f.','bc','all f.'};
        l = legend(legendstr);
        set(l,'Fontsize',16);
        sp = suptitle(['Comparison EMG with sensory feedback: '...
            'acceleration model (optimization time window)']);
        set(sp,'Fontsize',20);
        
        % Full time range
        ulf = struct('MS',[]);
        uvf = struct('MS',[]);
        uaf = struct('MS',[]);
        ulfinterp = struct('MS',[]);
        uvfinterp = struct('MS',[]);
        uafinterp = struct('MS',[]);
        EMG_path_agonist = struct('MS',[]);
        EMG_path_antagonist = struct('MS',[]);
        emg1processedstr = struct('MS',[]);
        emg2processedstr = struct('MS',[]);
        EMG_ordered = struct('MS',[]);
        bc_adj = struct('MS',[]);
        fit_IPSA = struct('r2',[]);
        figure()
        for ms = 1:Nsegment            
            for m = 1:NMuscles_spas                
                time_ref = round(outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time(1):...
                    0.005:outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time(end),4);
                time_temp = round(outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time,4);        
                % Here we chop off beginning and end of time range
                time_sel = time_ref(41:end-21);            
                % Get acceleration
                % 1) Get continuous time vector
                time_GPOPS = outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time;
                time_GPOPScont = time_GPOPS(1):(time_GPOPS(end) - ...
                    time_GPOPS(1))/(length(time_GPOPS)-1):time_GPOPS(end);
                time_GPOPScont = time_GPOPScont';
                time_GPOPScontround = round(time_GPOPScont,4);
                % 2) Interpolate vMtilde at evenly distributed points
                vMtildeInterp = interp1(time_GPOPS,vMtildeGPOPS.(...
                    ['segment_',num2str(segment_sel(ms))]),time_GPOPScont);
                % 3) Filter signal
                order = 4;
                cutoff_low = 10;
                fs=1/mean(diff(time_GPOPScont));
                [af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
                vMtildeInterp_filt = filtfilt(af,bf,vMtildeInterp);
                % 3) Spline approximation
                aMtildeInterp = zeros(size(time_GPOPScont,1),...
                    DatStore(1).MS.nMuscles);
                for mm = 1:DatStore(1).MS.nMuscles
                    ppy = spline(time_GPOPScont,vMtildeInterp_filt(:,mm));
                    [~,aMtildeInterp(:,mm),~] = ...
                        SplineEval_ppuval(ppy,time_GPOPScont,1);
                end
                % Forward integration  
                ulf(ms).MS(:,m) = lengthFDynamicsExplicit(...
                    outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time,...
                    lMtildeGPOPS.(['segment_',num2str(segment_sel(ms))])...
                    (:,auxdata.Muscles_spas(m)),auxdata.lMf_td,...
                    glf(ms).MS(m),100,auxdata.threshold_lMf(ms,m));
                uvf(ms).MS(:,m) = velocityFDynamicsExplicit(...
                    outputGPOPS.(['segment_',num2str(...
                    segment_sel(ms))]).result.solution.phase.time,...
                    vMtildeGPOPS.(['segment_',num2str(segment_sel(ms))])...
                    (:,auxdata.Muscles_spas(m)),auxdata.vMf_td,...
                    gvf(ms).MS(m),100,auxdata.threshold_vMf(ms,m));
                uaf(ms).MS(:,m) = accelerationFDynamicsExplicit(...
                    time_GPOPScont,aMtildeInterp(:,...
                    auxdata.Muscles_spas(m)),auxdata.aMf_td,...
                    gaf(ms).MS(m),100,auxdata.threshold_aMf(ms,m));
                % Interpolation 
                ulfinterp(ms).MS(:,m) = interp1(...
                    time_temp,ulf(ms).MS(:,m),time_sel');
                uvfinterp(ms).MS(:,m) = interp1(...
                    time_temp,uvf(ms).MS(:,m),time_sel');
                uafinterp(ms).MS(:,m) = interp1(...
                    time_GPOPScontround,uaf(ms).MS(:,m),time_sel');
                % Load EMG data                
                trial = ['segment_' int2str(segment_sel(ms))];  
                if strcmp(subject_name,'EF_r') || strcmp(subject_name,'EF_l')
                    EMG_path_agonist(ms).MS = fullfile(path_folder_HD,...
                    'OptimizedData',real_subject_name,'EMGProcessing',...
                    ['Stretch_',trial],...
                    'emg1processednorm_kneemuscles.mat');                
                    EMG_path_antagonist(ms).MS =  fullfile(path_folder_HD,...
                    'OptimizedData',real_subject_name,'EMGProcessing',...
                    ['Stretch_',trial],...
                    'emg2processednorm_kneemuscles.mat'); 
                 % Load motion range
                    Misc.range(ms).MS = load([path_folder_HD,'OpenSim\',...
                        real_subject_name,'\JointAngles\','IPSA\',...
                        ['JointAngles_',trial,'_range.mat']]);
                else
                    EMG_path_agonist(ms).MS = fullfile(path_folder_HD,...
                        'EMG','ProcessedStretch','MS',['Stretch_',trial],...
                        'emg1processednorm_kneemuscles.mat');
                    EMG_path_antagonist(ms).MS = fullfile(path_folder_HD,...
                        'EMG','ProcessedStretch','MS',['Stretch_',trial],...
                        'emg2processednorm_kneemuscles.mat');                    
                    % Load motion range
                    Misc.range(ms).MS = load([path_folder_HD '\OpenSim\',...
                        'SymModel\JointAngles\DataPellenberg\JointAngles_',...
                        trial,'_range.mat']);
                end
                load(EMG_path_agonist(ms).MS);
                load(EMG_path_antagonist(ms).MS);
                emg1processedstr(ms).MS = ...
                    emg1processednorm(Misc.range(ms).MS.range,:);
                emg2processedstr(ms).MS = ...
                    emg2processednorm(Misc.range(ms).MS.range,:);
                % Re-order EMG
                EMG_ordered(ms).MS = ...
                    zeros(size(emg1processedstr(ms).MS,1),NMuscles);   
                EMG_ordered(ms).MS(:,auxdata.Muscles_spas) = ...
                    emg1processedstr(ms).MS;
                EMG_ordered(ms).MS(:,auxdata.Muscles_antagonist) = ...
                    emg2processedstr(ms).MS;
                % Adjust EMG based on choped off bits during spastic 
                % gain estimation and here above
                EMG_ordered(ms).MS = EMG_ordered(ms).MS(10:end-10,:); 
                EMG_ordered(ms).MS = EMG_ordered(ms).MS(41:end-21,:); 
                % Baseline activation
                bc_adj(ms).MS = repmat(bc(ms).MS(1,:),size(time_sel,2),1); 
                % Compare measured EMG with simulated EMG (ie, muscle
                % activity resulting from the spasticy model + baseline
                % activation)
                EMG_meas = EMG_ordered(ms).MS(:,auxdata.Muscles_spas(m));
                EMG_sim = ulfinterp(ms).MS(:,m) + uvfinterp(ms).MS(:,m)+...
                    uafinterp(ms).MS(:,m) + bc_adj(ms).MS(:,m);
                % Compute RMSE and R²
                [fit_IPSA.r2(ms,m),fit_IPSA.rmse(ms,m)] = ...
                    rsquare(EMG_meas,EMG_sim);                
                subplot(Nsegment,NMuscles_spas,(ms-1)*NMuscles_spas+m)
                plot(time_sel, EMG_meas,'k','LineWidth', 3); hold on;
                plot(time_sel, uafinterp(ms).MS(:,m),'c','LineWidth',2);
                plot(time_sel, uvfinterp(ms).MS(:,m),'b','LineWidth',2);
                plot(time_sel, ulfinterp(ms).MS(:,m),'r','LineWidth',2);
                plot(time_sel, bc_adj(ms).MS(:,m),'g','LineWidth',2);
                plot(time_sel, EMG_sim,'m','LineWidth',2);                
                if ms == Nsegment
                    xlabel('Time (s)','Fontsize',16);
                end        
                if m == 1
                    ylabel('m. exc. ()','Fontsize',16);
                end        
                if ms == 1
                    title(name_muscles_spastic{m},'Fontsize',16);
                end
                ylim([0 inf]);
                set(gca,'Fontsize',16);     
            end
        end
        legendstr = {'EMG','aM f.','vM f.','lM f.','bc','all f.'};
        l = legend(legendstr);
        set(l,'Fontsize',16);
        sp = suptitle(['Comparison EMG with sensory feedback: '...
            'acceleration model (full time interval)']);     
        set(sp,'Fontsize',20);
        % Save fit for statistical analysis
        if savefits
            fitpath = [savefolder_sge,'\',formulation,'\fit_IPSA'];
            save([fitpath,'_joint_',number_save],'fit_IPSA');
        end   
end
end
end