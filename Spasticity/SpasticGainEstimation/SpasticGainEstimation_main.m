% This script estimates the spastic gains based on sensory information from
% passive fast stretch motions (IPSA)
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
% delta = 0.1;
% W.act       = 0.00100;
% W.aT        = 0.60000;
% W.vA        = 0.00001;
% W.dF        = 0.00001;
% W.EMG.all   = 0.003500-W.vA-W.dF;
% W.lMopt_max = 0.394500;
% W.lMopt     = 0.00100;
% W.act       = 0.00100-W.act*delta/(1-0.00100);
% W.aT        = 0.60000-W.aT*delta/(1-0.00100);
% W.vA        = 0.00001-W.vA*delta/(1-0.00100);
% W.dF        = 0.00001-W.dF*delta/(1-0.00100);
% W.EMG.all   = 0.003480-W.EMG.all*delta/(1-0.00100);
% W.lMopt_max = 0.394500-W.lMopt_max*delta/(1-0.00100);
% W.lMopt     = 0.00100+delta;
% nameParameters = [optMuscles,'_ublM',lMt_ub,...
%     '_Ntr',numbTrials,'_NISA',NISA,'_devp',num2str(dev_p),...
%     '_a',num2str(round(W.act*10000)),'_aT',num2str(round(W.aT*10000)),...
%     '_EMG',num2str(round(W.EMG.all*10000)),...
%     '_lM',num2str(round(W.lMopt_max*10000)),...
%     '_lMopt',num2str(round(W.lMopt*10000)),'_2sides'];

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
pathRepo_simcpspasticity = 'C:\Users\u0101727\Documents\MyRepositories\simcpspasticity_cases';
addpath(genpath(pathRepo));
pathOpenSimModel = [pathRepo,'\OpenSimModel\'];

% Three spasticity models are proposed: force-related model (muscle-tendon
% force and dF/dt feedback), velocity-related model (muscle fiber length
% and velocity feedback), and acceleration-related model (muscle fiber
% length, velocity, and acceleration feedback)
formulations = {'forceModel'};

% User settings
% If the user wants to use different EMG onsets, he should run 
% adjustEMGOnset first and set user-based as EMGonset
EMGonset = 'user-based'; % option: user-based / pre-selected
% If the user wants to use a different optimization time interval, he 
% should run selectOptimizationRange first and set user-based as OptRange
OptRange = 'user-based'; % option: user-based / pre-selected
% Different cases being considered
% Knee extension: IPSA of the hamstrings (1)
% Ankle dorsiflexion: IPSA of the gastrocnemii (2)
joints = {'knee','ankle'};
% ind: joint available for subject (1: knee, 2: ankle)
% comb: combination available among motions (cfr nc in code)
namesjoints.EF_r.ind = [1 2];
namesjoints.EF_r.comb = [1 1]; % 1, 2, 3, 4 tested for ankle comb

% Loop over formulations
for ff = 1:length(formulations)
    formulation = formulations{ff};
% Loop over subjects
for nn = 1:length(fieldnames(namesjoints))
    tempnames = fieldnames(namesjoints);
    Misc.subject_name = tempnames{nn};    
% Loop of spastic muscles / joints   
for jjj = 1:length(namesjoints.(Misc.subject_name).ind)
    jj = namesjoints.(Misc.subject_name).ind(jjj);
    % Indices of the fast trials
    switch Misc.subject_name
        case 'EF_r'
            switch joints{jj}
                case 'knee'
                    Misc.segment_sel_all = 6:8; 
                case 'ankle'
                    Misc.segment_sel_all = 19:22; 
            end        
    end            
    % Available combinations of 3 fast pasive stretch motions
    C3 = nchoosek(Misc.segment_sel_all,3);  
    % Pre-selected combination index
    nc = namesjoints.(Misc.subject_name).comb(jjj);
    % Combination selection
    segment_sel = C3(nc,:);
    Misc.segment_sel = segment_sel;

    % Load data from subjects
    switch Misc.subject_name
        case 'EF_r'
            real_subject_name = 'EF';
            model = 'EF_MRI_extROM'; 
            path_folder_HD = [pathRepo_simcpspasticity,'\Alldata\'];
            load([path_folder_HD,'\IPSA\EF\001217b123_segments_deleted.mat']);
            % Affected side
            sidename = 'RIGHT';        
    end
    Misc.savefolder = [pathOpenSimModel,real_subject_name,...
        '\Spasticity\SpasticGainEstimation'];
    Misc.savefolderForwardSimulation = [pathOpenSimModel,real_subject_name,...
        '\Spasticity\ForwardSimulations\IPSA_',nameParameters];
%     Misc.savefolderParameterEstimation = [path_folder_HD,...
%         'OptimizedData\',real_subject_name,'\ParameterEstimation'];
    % TODO
%     Misc.path_folder_HD = path_folder_HD;

    %% Various inputs 
    % Muscles and DOFs
    Misc.Nsegment = length(Misc.segment_sel);
    allsegments(1).side_name = sidename;
    switch allsegments(1).side_name
        case 'LEFT'
            switch joints{jj}
                case {'knee'}
                    vec = 13:24;
                case 'ankle'
                    vec = 13:24;
            end
            side = 'l';
            Uletter_side = 'L';
            vecAll = 1:43;
        case 'RIGHT'
            switch joints{jj}
                case {'knee'}
                    vec = 1:12;
                case 'ankle'
                    vec = 1:12;
            end
            side = 'r';
            Uletter_side = 'R';
            vecAll = 44:86;
    end
    getMuscleNamesIndices_MRI
    switch joints{jj}
        case {'knee','knee_flex'}
            Misc.MuscleNames_Input = MKnee_side;
            Misc.DofNames_Input={['knee_flex_',side]};
            Misc.DofNames_Input_moment=['knee_flex_',side,'_moment'];
        case 'ankle'
            Misc.MuscleNames_Input = MAnkle_side;
            Misc.DofNames_Input={['ankle_flex_',side]};
            Misc.DofNames_Input_moment=['ankle_flex_',side,'_moment'];
    end
    Misc.side = side;
    Misc.joint = joints{jj};        

    %% Paths
    IK_path = struct('MS',[]);
    ID_path = struct('MS',[]);
    MuscleAnalysis_path = struct('MS',[]);
    EMG_path_agonist = struct('MS',[]);
    EMG_path_antagonist = struct('MS',[]);
    for ms = 1:Misc.Nsegment
        trial = ['segment_' int2str(Misc.segment_sel(ms))];
        IK_path(ms).MS=fullfile(path_folder_HD,'OpenSim',...
            real_subject_name,'JointAngles','IPSA',['JointAngles_',...
            trial,'.mot']);
        ID_path(ms).MS=fullfile(path_folder_HD,'OpenSim',...
                real_subject_name,'ID','IPSA',['ID_',trial,'.sto']);
        MuscleAnalysis_path(ms).MS = fullfile(path_folder_HD,...
            'OpenSim',real_subject_name,'MuscleAnalysis','IPSA',...
            ['Stretch_',trial],['MuscleAnalysis_',...
            trial,'_MuscleAnalysis_']);
        EMG_path_agonist(ms).MS = [pathOpenSimModel,real_subject_name,...
            '\EMG\',['Stretch_',trial],'\emg1_norm_',nameParameters,'.mat'];
        EMG_path_antagonist(ms).MS = [pathOpenSimModel,real_subject_name,...
            '\EMG\',['Stretch_',trial],'\emg2_norm_',nameParameters,'.mat'];           
    end
    model_path=[pathOpenSimModel,real_subject_name,'\',model,'.osim'];
    
    %% Optional Input Arguments
    Misc.Atendon = [];              % Tendon Stiffness for selected muscles
    Misc.f_cutoff_ID = 8;           % cutoff frequency filtering ID
    Misc.f_order_ID = 4;            % order frequency filtering ID
    Misc.f_cutoff_lMT = 8;          % cutoff frequency filtering lMT
    Misc.f_order_lMT = 4;           % order frequency filtering lMT
    Misc.f_cutoff_dM= 8;            % cutoff frequency filtering MA
    Misc.f_order_dM = 4;            % order frequency filtering MA
    Misc.f_cutoff_IK= 8;            % cutoff frequency filtering IK
    Misc.f_order_IK = 4;            % order frequency filtering IK
            
    %% MT-parameters optimized during parameter optimization
    load([pathRepo,'\MTParameters\MTParameters_',real_subject_name,...
        '_MRI_opt_',nameParameters,'.mat']);
    MTParameters_side = MTParameters(:,vecAll);    
    for m = 1:length(Misc.MuscleNames_Input)
        Misc.params_scaled(:,m) = MTParameters_side(:,strcmp(Mall_side,Misc.MuscleNames_Input{m}));
    end

    %% Range of motion
    % Load EMG onset (possible manually adjusted). If the user wants to use
    % different EMG onsets, he should run adjustEMGOnset first and set
    % user-based as EMGonset
    switch EMGonset
        case 'pre-selected'
            load([Misc.savefolder,'\onset_man']);
        case 'user-based'
            load([Misc.savefolder,'\onset_man_user']);
    end             
    for ms = 1:Misc.Nsegment
        trial = ['segment_' int2str(segment_sel(ms))];
        Misc.allsegments(ms).MS = allsegments(segment_sel(ms));
        if strcmp(Misc.subject_name,'EF_r') || strcmp(Misc.subject_name,'EF_l')
            Misc.range(ms).MS = load([path_folder_HD,'OpenSim\',...
                    real_subject_name,'\JointAngles\','IPSA\',...
                    ['JointAngles_',trial,'_range.mat']]);
        else
            Misc.range(ms).MS = load([path_folder_HD '\OpenSim\',...
                'SymModel\JointAngles\DataPellenberg\JointAngles_',...
                trial,'_range.mat']); 
        end
        
        % They changed the way ISA are structured, now GAS corresponds to 1
        % MEH corresponds to 2, REF corresponds to 3
        if strcmp(real_subject_name,'EF') && ...
                strcmp(joints{jj},'knee')
            onoff_idx = 'onoff2';
        elseif strcmp(real_subject_name,'EF') && ...
                strcmp(joints{jj},'knee_flex')
            onoff_idx = 'onoff3';
        else
            onoff_idx = 'onoff1';
        end
        
        Misc.onoff1man(ms).MS = ...
            zeros(size(allsegments(segment_sel(ms)).(onoff_idx),1),1);    
        if ~isnan(onset_man(segment_sel(ms)))
            temp_first = onset_man(segment_sel(ms));
            temp_last = find(allsegments(segment_sel(ms)).(onoff_idx),1,'last');
            Misc.onoff1man(ms).MS(temp_first:temp_last,1) = 1;
        end    
    end
    % Load optimization time window (possible manually adjusted). If the 
    % user wants to use a different optimization time window, he should run 
    % selectOptimizationRange first and set user-based as EMGonset
    switch OptRange
        case 'pre-selected'
            load([Misc.savefolder,'\SelectRangeOnset_',joints{jj}]);
        case 'user-based'
            load([Misc.savefolder,'\SelectRangeOnset_',joints{jj},...
                '_user']);
    end   
    range_optimization = SelectRangeOnset3';
    % The intervals have been manually restrained for the forward 
    % simulations to account for bad data at beginning/end signals. We
    % adjust accordingly here with the '-9'.
    for ms = 1:Misc.Nsegment
        Misc.range_optimization(ms).MS = ...
            range_optimization(1,segment_sel(ms))-9 : ...
            range_optimization(2,segment_sel(ms))-9;
        % Duration interval of optimization (used in paper)
        % timeInterOpt(nn).(joints{jj})(ms) = ...
        %    length(Misc.range_optimization(ms).MS);
    end           
    
    %% Solve the problem 
    switch formulation                                          
        case 'forceModel'
            [output,DatStore] = forceModel_main(IK_path,ID_path,...
                MuscleAnalysis_path,EMG_path_agonist,...
                EMG_path_antagonist,Misc);  
    end
    number_save = [formulation '_' int2str(segment_sel)];
    if ~(exist([Misc.savefolder,'\',formulation,'_',nameParameters,'_',num2str(Misc.bspas)],'dir')==7)
        mkdir([Misc.savefolder,'\',formulation,'_',nameParameters,'_',num2str(Misc.bspas)]);
    end
    save([Misc.savefolder,'\',formulation,'_',nameParameters,'_',num2str(Misc.bspas),'\output_',joints{jj},'_',...
        number_save],'output');
    save([Misc.savefolder,'\',formulation,'_',nameParameters,'_',num2str(Misc.bspas),'\DatStore_',joints{jj},'_',...
        number_save],'DatStore'); 
    
end
end
end
