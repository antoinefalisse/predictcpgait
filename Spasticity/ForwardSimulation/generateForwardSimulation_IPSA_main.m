% This script generates EMG-driven forward simulations for passive stretch
% trials. The code is inspired from Falisse et al, PLoS ONE (2018):
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
Misc.pathMuscleModel = [pathRepo,'\MuscleModel\'];
addpath(genpath(Misc.pathMuscleModel));
pathVariousFunctions = [pathRepo,'\VariousFunctions\'];
addpath(genpath(pathVariousFunctions));

% Different cases being considered
% Knee extension: IPSA of the hamstrings (1)
% Ankle dorsiflexion: IPSA of the gastrocnemii (2)
joints = {'knee','ankle'};
namesjoints.subject1.ind = [1 2];
stretchvelocity = 'fast'; % other option is 'medium'
% Loop over subjects
for nn = 1:length(fieldnames(namesjoints))
    tempnames = fieldnames(namesjoints);
    Misc.subject_name = tempnames{nn};   
% Loop of spastic muscles / joints    
for jjj = 1:length(namesjoints.(Misc.subject_name).ind)      
    jj = namesjoints.(Misc.subject_name).ind(jjj);    
    % Indices of the fast trials
    switch Misc.subject_name
        case 'subject1'
            switch joints{jj}
                case 'knee'
                    switch stretchvelocity
                        case 'fast'
                            Misc.segment_sel_all = 6:8; % no onset 5
                        % no medium velocity
                    end                            
                case 'ankle'
                    switch stretchvelocity
                        case 'fast'
                            Misc.segment_sel_all = 19:22; 
                        % no medium velocity
                    end
            end
    end        
    % Load data from subjects
    switch Misc.subject_name
        case 'subject1'
            subject = Misc.subject_name;
            load([pathOpenSimModel,subject,'\IPSA\IPSA_data.mat']);
            % Affected side
            sidename = 'RIGHT';
    end        
    Misc.savefolder = [pathOpenSimModel,subject,...
        '\Spasticity\ForwardSimulations\IPSA'];
    
    %% Various inputs 
    % Muscles and DOFs
    Misc.Nsegment = length(Misc.segment_sel_all);
    allsegments(1).side_name = sidename;
    switch allsegments(1).side_name
        case 'LEFT'
            switch joints{jj}
                case {'knee','knee_flex'}
                    vec = 13:24;
                case 'ankle'
                    vec = 13:24;
            end
            side = 'l';
            Uletter_side = 'L';
            vecAll = 1:43;
        case 'RIGHT'
            switch joints{jj}
                case {'knee','knee_flex'}
                    vec = 1:12;
                case 'ankle'
                    vec = 1:12;
            end
            side = 'r';
            Uletter_side = 'R';
            vecAll = 44:86;
    end
    getMuscleNamesIndices
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
        trial = ['segment_' int2str(Misc.segment_sel_all(ms))];
        IK_path(ms).MS=fullfile(pathOpenSimModel,...
            subject,'IK','IPSA',['JointAngles_',trial,'.mot']);
        ID_path(ms).MS=fullfile(pathOpenSimModel,...
            subject,'ID','IPSA',['ID_',trial,'.sto']);
        MuscleAnalysis_path(ms).MS = fullfile(pathOpenSimModel,...
            subject,'MuscleAnalysis','IPSA',...
            ['Stretch_',trial],['MuscleAnalysis_',...
            trial,'_MuscleAnalysis_']);
        EMG_path_agonist(ms).MS = [pathOpenSimModel,subject,...
            '\EMG\IPSA\',['Stretch_',trial],'\emg1_norm_personalized.mat'];
        EMG_path_antagonist(ms).MS = [pathOpenSimModel,subject,...
            '\EMG\IPSA\',['Stretch_',trial],'\emg2_norm_personalized.mat'];            
    end

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
    load([pathParameterEstimation,'\Results\',subject,...
        '\MTParameters_personalized.mat']); 
    MTParameters_side = MTParameters(:,vecAll);    
    for m = 1:length(Misc.MuscleNames_Input)
        Misc.params_scaled(:,m) = ...
            MTParameters_side(:,strcmp(Mall_side,Misc.MuscleNames_Input{m}));
    end

    %% Range of motion  
    for ms = 1:Misc.Nsegment
        trial = ['segment_' int2str(Misc.segment_sel_all(ms))];
        Misc.allsegments(ms).MS = allsegments(Misc.segment_sel_all(ms));
        Misc.range(ms).MS = load([pathOpenSimModel,subject,...
            '\IK\','IPSA\',['JointAngles_',trial,'_range.mat']]);
    end

    %% Solve the problem 
    generateForwardSimulation_IPSA(IK_path,ID_path,MuscleAnalysis_path,...
        EMG_path_agonist,EMG_path_antagonist,Misc);

end
end
