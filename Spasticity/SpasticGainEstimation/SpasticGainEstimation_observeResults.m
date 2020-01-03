% This script plots EMG against simulated feedback muscle excitations.

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
savefits = 0; % 1 to save the fitting results (RMSE and R²), 0 otherwise

% Three spasticity models are proposed: force-related model (muscle-tendon
% force and dF/dt feedback), velocity-related model (muscle fiber length
% and velocity feedback), and acceleration-related model (muscle fiber
% length, velocity, and acceleration feedback)
formulations = {'forceModel'};

for ff = 1:length(formulations)
    formulation = formulations{ff};
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

switch subject_name
    case 'subject1'
        subject = subject_name;
end

% Load data from spastic gain estimation (sge)
number_save = [formulation '_' int2str(segment_sel)];
savefolder_sge = [pathOpenSimModel,subject,...
    '/Spasticity/SpasticGainEstimation'];
savefolder_ForwardSimulation = [pathOpenSimModel,subject,...
    '/Spasticity/ForwardSimulations/IPSA'];
load([savefolder_sge,'/',formulation,'/output_',joint,'_',number_save]);
load([savefolder_sge,'/',formulation,'/DatStore_',joint,'_',number_save]);
% Load data from forward simulations
load([savefolder_ForwardSimulation,'/FtildeGPOPSFTtilde_',joint]);
load([savefolder_ForwardSimulation,'/dFtildeGPOPSFTtilde_',joint]);
load([savefolder_ForwardSimulation,'/lMtildeGPOPSFTtilde_',joint]);
load([savefolder_ForwardSimulation,'/vMtildeGPOPSFTtilde_',joint]);
load([savefolder_ForwardSimulation,'/outputGPOPSFTtilde_',joint]);

switch formulation
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
                EMG_path_agonist(ms).MS = fullfile(pathOpenSimModel,...
                    subject,'/EMG/IPSA/',['Stretch_',trial],...
                    '/emg1_norm_personalized.mat');
                EMG_path_antagonist(ms).MS = fullfile(pathOpenSimModel,...
                    subject,'/EMG/IPSA/',['Stretch_',trial],...
                    '/emg2_norm_personalized.mat');                    
                % Load motion range
                Misc.range(ms).MS = load([pathOpenSimModel,subject,...
                    '/IK/','IPSA/',['JointAngles_',trial,'_range.mat']]);
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
            fitpath = [savefolder_sge,'/',formulation,'/fit_IPSA'];
            save([fitpath,'_joint_',number_save],'fit_IPSA');
        end     
end
end
end