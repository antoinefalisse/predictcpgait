%%  Muscle-tendon parameter optimal estimation
% This script estimates optimal fiber lengths and tendon slack lengths
% while solving the muscle rendundancy problem simultaneously for gait
% trials and passive stretches. The problem is solved via direct
% collocation. Details of the problem formualtion are in the associated
% publication.
%
% The code is not perfectly suited for use with data from different
% subjects although only minimal adjustment should be made. The code is
% also pretty messy and likely difficult to use with different data
% structure.
%
% Author: Antoine Falisse
% Date: 1/2/2020
%
clear all
close all
clc

%% User settings
solveProblem = 0;   % set to 1 to solve problem
saveResults = 0;    % set to 1 to save results
saveParameters = 1; % set to 1 to save MT-parameters
showResults = 1;    % set to 1 to plot results

%% Problem settings
subject = 'subject1';  
lMtilde_ext.max = 1.5;  % reference maximum normalized fiber length
lMtilde_ext.min = 0.5;  % reference minimum normalized fiber length
dev_p = 5; % parameters cannot deviate by more than +/- dev_p/100% between legs
EMGScale_min = 0.01; % minimum value for EMG scale factors
NMeshes = 50; % number of mesh intervals
NGaitTrials = '4'; % number of gait trials
NISATrials = '6'; % number of ISA trials
sides = {'l','r'}; NSides = num2str(length(sides));
% Weight factors in cost function
W.a         = 0.00100; % muscle activations
W.aT        = 0.60000; % reserve actuators
W.vA        = 0.00001; % time derivative of muscle activations
W.dF        = 0.00001; % time derivative of muscle forces
W.EMG.all   = 0.003500-W.vA-W.dF; % deviation from EMG
W.lMopt_max = 0.394500; % deviation from maximum normalized fiber length
W.lMopt     = 0.00100; % optimal fiber lengths (short fibers)
delta       = 0.1; % parameter to adjust weight factors (see below)

%% Helper variables / paths
% Paths
pathMain = pwd;
[pathRepo,~,~] = fileparts(pathMain);
pathOpenSimModel = [pathRepo,'/OpenSimModel'];
pathCollocationScheme = [pathRepo,'/CollocationScheme'];
addpath(genpath(pathCollocationScheme));
pathResults = [pathOpenSimModel,'/',subject,'/ParameterEstimation/'];
pathMuscleModel = [pathRepo,'/MuscleModel/'];
% Adjusted weights to preserve sum=1 when increasing weight on short
% fibers. The adjustment accounts for the magnitude of the weight factors.
W.a         = W.a-W.a*delta/(1-0.00100);
W.aT        = W.aT-W.aT*delta/(1-0.00100);
W.vA        = W.vA-W.vA*delta/(1-0.00100);
W.dF        = W.dF-W.dF*delta/(1-0.00100);
W.EMG.all   = W.EMG.all-W.EMG.all*delta/(1-0.00100);
W.lMopt_max = W.lMopt_max-W.lMopt_max*delta/(1-0.00100);
W.lMopt     = W.lMopt+delta;
W.all = [W.a,W.aT,W.vA,W.dF,W.EMG.all,W.lMopt_max,W.lMopt];
if round(sum(W.all),10)~=1
    disp('WARNING: the sum of the weights is not 1');
end
% Weight factors for EMG tracking, can differ between muscles.
W.EMG.RF = 0.1;
W.EMG.BFLH = 1;
W.EMG.ST = 1;
W.EMG.SM = 1;
W.EMG.TA = 1;
W.EMG.GL = 1;
W.EMG.GM = 1;
W.EMG.SOL = 1;
W.EMG.GluM = 1;
W.EMG.ind = [W.EMG.RF,W.EMG.BFLH,W.EMG.ST, W.EMG.SM,W.EMG.TA,W.EMG.GL,...
    W.EMG.GM,W.EMG.SOL,W.EMG.GluM];

%% Solve optimal control problem
if solveProblem
% Degrees of freedom taken into account when solving the muscle redundancy
% problem for gait trials
gaitDOFNames = {'hip_flex_','hip_add_','hip_rot_','knee_flex_',...
    'ankle_flex_','subt_angle_'};
% Muscles
MuscleNames = {'glut_max1_','glut_max2_','glut_max3_','glut_med1_',...
    'glut_med2_','glut_med3_','glut_min1_','glut_min2_',...
    'glut_min3_','add_long_','add_brev_','add_mag1_','add_mag2_',...
    'add_mag3_','pectineus_','iliacus_','psoas_','quad_fem_',...
    'gemellus_','piri_','TFL_','gracilis_','semimem_','semiten_',...
    'bi_fem_lh_','bi_fem_sh_','sartorius_','rectus_fem_',...
    'vas_med_','vas_int_','vas_lat_','gas_med_',...
    'gas_lat_','soleus_','tib_post_','tib_ant_','ext_dig_',...
    'ext_hal_','flex_dig_','flex_hal_','per_brev_','per_long_','per_tert_'};
% Muscles for which we optimize optimal fiber lengths and tendon slack
% lengths
MuscleNames_opt = {'glut_max1_','glut_max2_','glut_max3_','glut_med1_',...
    'glut_med2_','glut_med3_','glut_min1_','glut_min2_',...
    'glut_min3_','add_long_','add_brev_','add_mag1_','add_mag2_',...
    'add_mag3_','pectineus_','iliacus_','psoas_','quad_fem_',...
    'gemellus_','piri_','TFL_','gracilis_','semimem_','semiten_',...
    'bi_fem_lh_','bi_fem_sh_','sartorius_','rectus_fem_',...
    'vas_med_','vas_int_','vas_lat_','gas_med_',...
    'gas_lat_','soleus_','tib_post_','tib_ant_','ext_dig_',...
    'ext_hal_','flex_dig_','flex_hal_','per_brev_','per_long_','per_tert_'};
% Muscles for which we have EMGs
EMGchannels = {'REF','VAL','BIF','MEH','TIA','GAS','SOL','GLU'};
% Load parameters of FLV relationships
load([pathMuscleModel,'Faparam.mat']);
load([pathMuscleModel,'Fpparam.mat']);
load([pathMuscleModel,'Fvparam.mat']);

%% Extract and arrange data
NTr = 0;
for i = 1:length(sides)
    side_tr = sides{i};
    % For subject 1, there is no ISA for the left leg. Might different for
    % othe subjects.
    if strcmp(side_tr,'r')
        if strcmp(NISATrials,'0')
            ISA = 'no';
        else
            ISA = 'yes';
        end
    else         
        ISA = 'no';
    end        
    switch subject
        case 'subject1'
            modelName = ...
                [pathOpenSimModel,'/',subject,'/',subject,'_MRI_extROM.osim'];
            % Available trials in database. Only trials used in results of
            % publication are shared.
            gait_trials_id = {'10','10','11','12','12','13','13','15','16',...
                '17','20','20','22','22','23','25','27','28','29','30'};
            gait_trials_leg = {'r','l','r','r','l','r','l','r','r','l','r',...
                'l','r','l','l','r','r','r','r','l'};
            gait_trials_start_time = [2.4,2.92,4.43,3,2.45,1.65,1.16,2.51,...
                2.36,2.98,2.2,1.78,1.28,1.74,2.8,1.64,3.23,4.54,3.09,2.52];            
            gait_trials_end_time = [3.45,3.94,5.48,4.01,3.52,2.64,2.17,3.51,...
                3.44,4.07,3.01,2.58,2.14,2.48,3.53,2.29,4.68,5.96,4.16,3.52];         
            % Select trials
            if strcmp(side_tr,'r')                
                switch NGaitTrials
                    case '4'
                        trials_number = {'10','12','13','22'};
                    case '3'
                        trials_number = {'10','12','13'};
                    case '2'
                        trials_number = {'10','12'};
                    case '1'
                        trials_number = {'10'};
                end
                if strcmp(ISA,'yes')
                    if strcmp(NISATrials,'6')
                        ISA_trials_id = {'2','10','15','6','14','19'};
                        ISADOFNames = {'knee_flex_','knee_flex_',...
                            'ankle_flex_','knee_flex_','knee_flex_',...
                            'ankle_flex_'};
                        ISADOFTypes = {'knee_ext','knee_flex','ankle_dor',...
                            'knee_ext','knee_flex','ankle_dor'};
                        ISA_trials_time(1,:) = [5,4,1,5.8,5.1,3];
                        ISA_trials_time(2,:) = [12,9,5,6.8,6.1,4];
                    elseif strcmp(NISATrials,'3')
                        ISA_trials_id = {'6','14','19'};
                        ISADOFNames = {'knee_flex_','knee_flex_','ankle_flex_'};
                        ISADOFTypes = {'knee_ext','knee_flex','ankle_dor'};
                        ISA_trials_time(1,:) = [5.8,5.1,3];
                        ISA_trials_time(2,:) = [6.8,6.1,4];
                        
                    end
                end                
            elseif strcmp(side_tr,'l')         
                switch NGaitTrials
                    case '4'
                        trials_number = {'10','12','13','22'};                         
                    case '3'
                        trials_number = {'10','12','13'};                        
                    case '2'             
                        trials_number = {'10','12'};
                    case '1'             
                        trials_number = {'10'};
                end                  
            end
            count = 1;
            for tt = 1:length(gait_trials_leg)
                if strcmp(gait_trials_leg{tt},side_tr)
                    idx_leg(count) = tt;                        
                    count = count + 1;
                end
            end 
            gait_trials_start_time_sel = gait_trials_start_time(idx_leg);
            gait_trials_end_time_sel = gait_trials_end_time(idx_leg);
            gait_trials_id_sel = gait_trials_id(idx_leg);  
            count = 1;
            for tt = 1:length(trials_number)
                trials{tt} = ['Gait_',trials_number{tt}];   
                if find(strcmp(trials_number{tt},gait_trials_id_sel))
                    idx_leg_tr(count) = ...
                        find(strcmp(trials_number{tt},gait_trials_id_sel));                        
                    count = count + 1;
                end
            end  
            time(1,:) = gait_trials_start_time_sel(idx_leg_tr);
            time(2,:) = gait_trials_end_time_sel(idx_leg_tr)-0.05; 
            % Muscles driven by EMGs
            % The VAL EMG is pretty bad so it is ignored
            EMGchannelsOpt = ...
                {'REF','BIF','MEH','MEH','TIA','GAS','GAS','SOL','GLU'};
            Opt_channels = {'rectus_fem_','bi_fem_lh_','semiten_','semimem_',...
                'tib_ant_','gas_lat_','gas_med_','soleus_','glut_med2_'};
            idx_EMGscale = 1:length(EMGchannelsOpt);
            clinicalReport.xlsF = 'clinicalExam.xlsx';
            pathClinicalReport = ...
                [pathOpenSimModel,'/',subject,'/ClinicalExam/'];
            clinicalReport.xlsD = pathClinicalReport;
    end 
    idxEMG = zeros(length(EMGchannelsOpt),2);
    for ii = 1:length(EMGchannelsOpt)
        idxEMG(ii,1) = find(strcmp(EMGchannels,EMGchannelsOpt{ii}));
        idxEMG(ii,2) = find(strcmp(MuscleNames,Opt_channels{ii}));    
        idxEMG(ii,3) = find(strcmp(MuscleNames,Opt_channels{ii})); 
    end
   % Load linearly-scaled MT-parameters
   load([pathResults,'/MTParameters_generic.mat']);
   NMuscles = size(MTParameters,2);
   % Load EMG
   EMGpath = [pathOpenSimModel,'/',subject,'/EMG/Gait/EMG_filt.mat'];
   EMGall = load(EMGpath);
   
   %% Extract data from gait trials
   indEMG2misc_all = [];
   for j = 1:length(trials)
       trial = trials{j};
       trial_number = trials_number{j};
       leg = side_tr;
       leg_up = 'L';
       if strcmp(leg,'r')
           leg_up = 'R';
       end
       side_t{j+NTr} = sides{i};
       Fmax{j+NTr} = MTParameters(1,1:NMuscles/2);
       ParametersInit.(sides{i}).OptFibLen = MTParameters(2,1:NMuscles/2);
       ParametersInit.(sides{i}).TenSlackLen = MTParameters(3,1:NMuscles/2);
       ParametersInit.(sides{i}).PennAng = MTParameters(4,1:NMuscles/2);
       if strcmp(leg,'r')
           Fmax{j+NTr} = MTParameters(1,NMuscles/2+1:NMuscles);
           ParametersInit.(sides{i}).OptFibLen = ...
               MTParameters(2,NMuscles/2+1:NMuscles);
           ParametersInit.(sides{i}).TenSlackLen = ...
               MTParameters(3,NMuscles/2+1:NMuscles);
           ParametersInit.(sides{i}).PennAng = ...
               MTParameters(4,NMuscles/2+1:NMuscles);
       end
       % get ID
       pathID = [pathOpenSimModel,'/',subject,'/ID/Gait/ID_',trial,leg,...
           '_MRI_extROM.sto'];
       ID = importdata(pathID);
       for k = 1:length(gaitDOFNames)
           JointMom{j+NTr}(:,k) = ...
               ID.data(:,strcmp(ID.colheaders,[gaitDOFNames{k},leg,'_moment']));           
       end
       % select interval and interpolate
       step = (round(time(2,j),2) - round(time(1,j),2))/(NMeshes-1);
       interval = round(time(1,j),2):step:round(time(2,j),2);
       timeOpt{j+NTr}(1) = interval(1); timeOpt{j+NTr}(2) = interval(end); 
       JointMom{j+NTr} = interp1(ID.data(:,1),JointMom{j+NTr},interval);      
       % get lMT
       lMTpath = [pathOpenSimModel,'/',subject,'/MuscleAnalysis/Gait/',...
           trial,leg,'/',trial,leg,'_MuscleAnalysis_Length.sto'];
       lMT = importdata(lMTpath);
       for k = 1:length(MuscleNames)
           MuscTenLen{j+NTr}(:,k) = ...
               lMT.data(:,strcmp(lMT.colheaders,[MuscleNames{k},leg]));      
       end 
       % select interval and interpolate
       MuscTenLen{j+NTr} = interp1(lMT.data(:,1),MuscTenLen{j+NTr},interval); 
       % compute muscle tendon velocities (vMT)
       for m = 1:length(MuscleNames)
           pp_y = spline(interval,MuscTenLen{j+NTr}(:,m));
           [~,MuscTenVel{j+NTr}(:,m),~] = SplineEval_ppuval(pp_y,interval,1);
       end      
       % get MA
       MApath = [pathOpenSimModel,'/',subject,'/MuscleAnalysis/Gait/',...
           trial,leg,'/',trial,leg,'_MuscleAnalysis_MomentArm_'];
       for d = 1:length(gaitDOFNames)
           MApathDOF = [MApath,gaitDOFNames{d},leg,'.sto'];
           MA = importdata(MApathDOF);           
           for k = 1:length(MuscleNames)
               MuscMomArm{j+NTr}(:,k,d) = ...
                   MA.data(:,strcmp(MA.colheaders,[MuscleNames{k},leg])); 
           end           
       end
       idx_spanning{j+NTr} = sum(squeeze(sum(MuscMomArm{j+NTr},1))',1);
       idx_spanning{j+NTr}(idx_spanning{j+NTr}<=0.0001 & ...
           idx_spanning{j+NTr}>=-0.0001) = 0;
       idx_spanning{j+NTr}(idx_spanning{j+NTr}~=0) = 1;
       Fmax{j+NTr} = Fmax{j+NTr}(idx_spanning{j+NTr}==1);
       % select interval and interpolate
       MuscMomArm{j+NTr} = interp1(MA.data(:,1),MuscMomArm{j+NTr},interval);
       % get EMG
       EMG = EMGall.EMG_filt.(['gait_',trial_number]);       
       % select interval and interpolate
       idx_EMG(1) = find(round(EMG.data(:,1),3)==round(time(1,j),2));
       idx_EMG(2) = find(round(EMG.data(:,1),3)==round(time(2,j),2));
       step = (EMG.data(idx_EMG(2),1)-EMG.data(idx_EMG(1),1))/(NMeshes-1);
       intervalEMG = EMG.data(idx_EMG(1),1):step:EMG.data(idx_EMG(2),1);
       EMG.datainterp = interp1(EMG.data(idx_EMG(1):idx_EMG(2),1),...
           EMG.data(idx_EMG(1):idx_EMG(2),:),intervalEMG);        
       for k = 1:length(EMGchannels)
           EMGopt{j+NTr}(:,k) = ...
               EMG.datainterp(:,strcmp(EMG.colheaders,[leg_up,EMGchannels{k}]));
           maxEMG{j}(:,k) = max(EMGopt{j+NTr}(:,k));
       end    
       indEMG2misc{j+NTr} = idxEMG;
       idx_EMGscale2misc{j+NTr} = idx_EMGscale;
       indEMG2misc_all = [indEMG2misc_all;idxEMG];
       % get muscleNames       
       for k = 1:length(MuscleNames)
        muscleNames.(leg).s{k} = [MuscleNames{k},leg];
       end
       MuscleNames_Spanning = MuscleNames(idx_spanning{j+NTr}==1);
       % get Side
       Side = leg_up;
       % Helper vector
       % SelectedMusclesInfo contains the indices of the muscles selected
       % for the optimization in the vector with all the muscles
        count = 1;
        for ii = 1:length(MuscleNames)
            if find(strcmp(MuscleNames_opt,MuscleNames{ii}))
                SelectedMusclesInfo{j+NTr}.Bool(ii) = true;       
                SelectedMusclesInfo{j+NTr}.Index(count) = ii;
                count = count+1;
            else
                SelectedMusclesInfo{j+NTr}.Bool(ii) = false; 
            end
        end
       % Helper vector
       % SelectedMusclesInfoSelSpanning contains the indices of the muscles
       % spanning the degrees of freedom of that trial in the vector of
       % muscles selected for the optimization
        count = 1;
        for ii = 1:length(MuscleNames_opt)
            if find(strcmp(MuscleNames_Spanning,MuscleNames_opt{ii}))
                SelectedMusclesInfoSelSpanning{j+NTr}.Bool(ii) = true;       
                SelectedMusclesInfoSelSpanning{j+NTr}.Index(count) = ii;
                count = count+1;
            else
                SelectedMusclesInfoSelSpanning{j+NTr}.Bool(ii) = false; 
            end
        end     
        % Helper vector
        % SelectedMusclesInfoSpanningSel contains the indices of the
        % muscles selected for the optimization in the vector of
        % muscles spanning the degrees of freedom
        count = 1;
        for ii = 1:length(MuscleNames_Spanning)
            if find(strcmp(MuscleNames_opt,MuscleNames_Spanning{ii}))
                SelectedMusclesInfoSpanningSel{j+NTr}.Bool(ii) = true;       
                SelectedMusclesInfoSpanningSel{j+NTr}.Index(count) = ii;
                count = count+1;
            else
                SelectedMusclesInfoSpanningSel{j+NTr}.Bool(ii) = false; 
            end
        end 
        % Helper vector
        % SelectedMusclesInfoOptSpanning contains the indices of the muscles
        % selected for the optimization and spanning the degrees of freedom
        count = 1;
        for ii = 1:length(MuscleNames)
            if (find(strcmp(MuscleNames_opt,MuscleNames{ii}))) ...
                    & (find(strcmp(MuscleNames_Spanning,MuscleNames{ii})))
                SelectedMusclesInfoOptSpanning{j+NTr}.Bool(ii) = true;       
                SelectedMusclesInfoOptSpanning{j+NTr}.Index(count) = ii;
                count = count+1;
            else
                SelectedMusclesInfoOptSpanning{j+NTr}.Bool(ii) = false; 
            end
        end        
    end
    
    %% Extract data from ISA trials (if available)
    if strcmp(ISA,'yes')
    for j = 1:length(ISA_trials_id)
        side_t{j+NTr+length(trials)} = sides{i};
        Fmax{j+NTr+length(trials)} = MTParameters(1,1:NMuscles/2);
        % get ID
        pathID = [pathOpenSimModel,'/',subject,'/ID/IPSA/ID_segment_',...
            ISA_trials_id{j},'.sto'];
        ID = importdata(pathID);
        JointMom{j+NTr+length(trials)} = ID.data(:,strcmp(ID.colheaders,...
            [ISADOFNames{j},leg,'_moment']));
        % select interval and interpolate
        step = (round(ISA_trials_time(2,j),2) - ...
            round(ISA_trials_time(1,j),2)) / (NMeshes-1);
        interval = ...
            round(ISA_trials_time(1,j),2):step:round(ISA_trials_time(2,j),2);
        timeOpt{j+NTr+length(trials)}(1) = interval(1); 
        timeOpt{j+NTr+length(trials)}(2) = interval(end); 
        JointMom{j+NTr+length(trials)} = ...
            interp1(ID.data(:,1),JointMom{j+NTr+length(trials)},interval)';           
        % get MA
        MApath = [pathOpenSimModel,'/',subject,...
            '/MuscleAnalysis/IPSA/Stretch_segment_',ISA_trials_id{j},...
            '/MuscleAnalysis_segment_',ISA_trials_id{j},...
            '_MuscleAnalysis_MomentArm_'];
        MApathDOF = [MApath,ISADOFNames{j},leg,'.sto'];
        MA = importdata(MApathDOF);           
        for k = 1:length(MuscleNames)
           MuscMomArm_all{j+NTr+length(trials)}(:,k,1) = ...
               MA.data(:,strcmp(MA.colheaders,[MuscleNames{k},leg])); 
        end           
        % select only muscles actuating the joints
        idx_spanning{j+NTr+length(trials)} = ...
            sum(MuscMomArm_all{j+NTr+length(trials)},1);
        idx_spanning{j+NTr+length(trials)}(idx_spanning{j+NTr+length(trials)}...
            <= 0.0001 & idx_spanning{j+NTr+length(trials)}>=-0.0001) = 0;
        idx_spanning{j+NTr+length(trials)}...
            (idx_spanning{j+NTr+length(trials)}~=0) = 1;
        temp = MuscMomArm_all{j+NTr+length(trials)};
        temp(:,idx_spanning{j+NTr+length(trials)}==0)=[];
        MuscMomArm{j+NTr+length(trials)}(:,:,1) = temp; clear temp; 
        MuscleNames_Spanning = ...
            MuscleNames(idx_spanning{j+NTr+length(trials)}==1);
        Fmax{j+NTr+length(trials)} = Fmax{j+NTr+length(trials)}...
            (idx_spanning{j+NTr+length(trials)}==1);           
        % select interval and interpolate
        MuscMomArm{j+NTr+length(trials)} = ...
            interp1(MA.data(:,1),MuscMomArm{j+NTr+length(trials)},interval);
        % get lMT
        lMTpath = [pathOpenSimModel,'/',subject,...
            '/MuscleAnalysis/IPSA/Stretch_segment_',ISA_trials_id{j},...
            '/MuscleAnalysis_segment_',ISA_trials_id{j},...
            '_MuscleAnalysis_Length.sto'];
        lMT = importdata(lMTpath);
        for k = 1:length(MuscleNames)
           MuscTenLen_all{j+NTr+length(trials)}(:,k) = ...
               lMT.data(:,strcmp(lMT.colheaders,[MuscleNames{k},leg]));      
        end 
        % select only muscles actuating the joints
        temp = MuscTenLen_all{j+NTr+length(trials)};
        temp(:,idx_spanning{j+NTr+length(trials)}==0)=[];
        MuscTenLen{j+NTr+length(trials)}(:,:) = temp;          
        clear temp;
        % select interval and interpolate
        MuscTenLen{j+NTr+length(trials)} = ...
            interp1(lMT.data(:,1),MuscTenLen{j+NTr+length(trials)},interval); 
        % compute muscle tendon velocities
        for m = 1:size(MuscTenLen{j+NTr+length(trials)},2)
           pp_y = spline(interval,MuscTenLen{j+NTr+length(trials)}(:,m));
           [~,MuscTenVel{j+NTr+length(trials)}(:,m),~] = ...
               SplineEval_ppuval(pp_y,interval,1);
        end             
        % get EMG
        EMGpath1 = [pathOpenSimModel,'/',subject,...
            '/EMG/IPSA/Stretch_segment_',ISA_trials_id{j},'/emg1_filt.mat'];
        EMGall1 = load(EMGpath1);
        EMGpath2 = [pathOpenSimModel,'/',subject,...
            '/EMG/IPSA/Stretch_segment_',ISA_trials_id{j},'/emg2_filt.mat'];
        EMGall2 = load(EMGpath2);
        % select range in EMG corresponding to motion
        rangepath = [pathOpenSimModel,'/',subject,...
            '/IK/IPSA/JointAngles_segment_',ISA_trials_id{j},'_range.mat'];
        load(rangepath);
        EMGall1.emg1processed = EMGall1.emg1processed(range);
        EMGall2.emg2processed = EMGall2.emg2processed(range);
        % create corresponding time vector
        step = (ID.data(end,1)-ID.data(1,1))/(size(EMGall1.emg1processed,1)-1);
        EMGtime = ID.data(1,1):step:ID.data(end,1);
        EMGopt{j+NTr+length(trials)}(:,1) = ...
            interp1(EMGtime,EMGall1.emg1processed,interval);
        EMGopt{j+NTr+length(trials)}(:,2) = ...
            interp1(EMGtime,EMGall2.emg2processed,interval);           
        maxEMG{j+length(trials)} = max(EMGopt{j+NTr+length(trials)},[],1);
        switch ISADOFTypes{j}
            case 'knee_ext'
               indEMG2misc{1,j+NTr+length(trials)}(1:2,1)   = 1;
               indEMG2misc{1,j+NTr+length(trials)}(3,1)     = 2;
               indEMG2misc{1,j+NTr+length(trials)}(1,2)     = ...
                   find(strcmp(MuscleNames,'semimem_'));
               indEMG2misc{1,j+NTr+length(trials)}(2,2)     = ...
                   find(strcmp(MuscleNames,'semiten_'));
               indEMG2misc{1,j+NTr+length(trials)}(3,2)     = ...
                   find(strcmp(MuscleNames,'rectus_fem_'));                      
               indEMG2misc{1,j+NTr+length(trials)}(1,3)     = ...
                   find(strcmp(MuscleNames_Spanning,'semimem_'));
               indEMG2misc{1,j+NTr+length(trials)}(2,3)     = ...
                   find(strcmp(MuscleNames_Spanning,'semiten_'));
               indEMG2misc{1,j+NTr+length(trials)}(3,3)     = ...
                   find(strcmp(MuscleNames_Spanning,'rectus_fem_'));                    
               temp = indEMG2misc{1,j+NTr+length(trials)}(:,2);
               for n = 1:length(temp)
                   idx_EMGscale2misc{j+NTr+length(trials)}(n) = ...
                       find(indEMG2misc{1,1}(:,2) == temp(n));
               end  
               indEMGISAinGait{j}(1) = find(strcmp(EMGchannels,'MEH'));
               indEMGISAinGait{j}(2) = find(strcmp(EMGchannels,'REF'));
            case 'knee_flex'
                indEMG2misc{1,j+NTr+length(trials)}(1,1)     = 1;
                indEMG2misc{1,j+NTr+length(trials)}(2:3,1)   = 2;
                indEMG2misc{1,j+NTr+length(trials)}(1,2)     = ...
                    find(strcmp(MuscleNames,'rectus_fem_'));
                indEMG2misc{1,j+NTr+length(trials)}(2,2)     = ...
                    find(strcmp(MuscleNames,'semimem_'));
                indEMG2misc{1,j+NTr+length(trials)}(3,2)     = ...
                    find(strcmp(MuscleNames,'semiten_'));    
                indEMG2misc{1,j+NTr+length(trials)}(1,3)     = ...
                    find(strcmp(MuscleNames_Spanning,'rectus_fem_'));
                indEMG2misc{1,j+NTr+length(trials)}(2,3)     = ...
                    find(strcmp(MuscleNames_Spanning,'semimem_'));
                indEMG2misc{1,j+NTr+length(trials)}(3,3)     = ...
                    find(strcmp(MuscleNames_Spanning,'semiten_'));  
                temp = indEMG2misc{1,j+NTr+length(trials)}(:,2);
                for n = 1:length(temp)
                    idx_EMGscale2misc{j+NTr+length(trials)}(n) = ...
                        find(indEMG2misc{1,1}(:,2) == temp(n));
                end       
                indEMGISAinGait{j}(1) = find(strcmp(EMGchannels,'REF'));
                indEMGISAinGait{j}(2) = find(strcmp(EMGchannels,'MEH'));
            case 'ankle_dor'
                indEMG2misc{1,j+NTr+length(trials)}(1:2,1)   = 1;
                indEMG2misc{1,j+NTr+length(trials)}(3,1)     = 2;
                indEMG2misc{1,j+NTr+length(trials)}(1,2)     = ...
                    find(strcmp(MuscleNames,'gas_med_'));
                indEMG2misc{1,j+NTr+length(trials)}(2,2)     = ...
                    find(strcmp(MuscleNames,'gas_lat_'));
                indEMG2misc{1,j+NTr+length(trials)}(3,2)     = ...
                    find(strcmp(MuscleNames,'tib_ant_'));     
                indEMG2misc{1,j+NTr+length(trials)}(1,3)     = ...
                    find(strcmp(MuscleNames_Spanning,'gas_med_'));
                indEMG2misc{1,j+NTr+length(trials)}(2,3)     = ...
                    find(strcmp(MuscleNames_Spanning,'gas_lat_'));
                indEMG2misc{1,j+NTr+length(trials)}(3,3)     = ...
                    find(strcmp(MuscleNames_Spanning,'tib_ant_'));  
                temp = indEMG2misc{1,j+NTr+length(trials)}(:,2);
                for n = 1:length(temp)
                    idx_EMGscale2misc{j+NTr+length(trials)}(n) = ...
                        find(indEMG2misc{1,1}(:,2) == temp(n));
                end       
                indEMGISAinGait{j}(1) = find(strcmp(EMGchannels,'GAS'));
                indEMGISAinGait{j}(2) = find(strcmp(EMGchannels,'TIA'));
        end 
        indEMG2misc_all = [indEMG2misc_all;indEMG2misc{1,j+NTr+length(trials)}];        
        % Helper vector
        % SelectedMusclesInfo contains the indices of the muscles selected
        % for the optimization in the vector with all the muscles
        count = 1;
        for ii = 1:length(MuscleNames)
            if find(strcmp(MuscleNames_opt,MuscleNames{ii}))
                SelectedMusclesInfo{j+NTr+length(trials)}.Bool(ii) = true;       
                SelectedMusclesInfo{j+NTr+length(trials)}.Index(count) = ii;
                count = count+1;
            else
                SelectedMusclesInfo{j+NTr+length(trials)}.Bool(ii) = false; 
            end
        end
        % Helper vector
        % SelectedMusclesInfoSelSpanning contains the indices of the muscles
        % spanning the degrees of freedom of that trial in the vector of
        % muscles selected for the optimization
        count = 1;
        for ii = 1:length(MuscleNames_opt)
            if find(strcmp(MuscleNames_Spanning,MuscleNames_opt{ii}))
                SelectedMusclesInfoSelSpanning{...
                    j+NTr+length(trials)}.Bool(ii) = true;       
                SelectedMusclesInfoSelSpanning{...
                    j+NTr+length(trials)}.Index(count) = ii;
                count = count+1;
            else
                SelectedMusclesInfoSelSpanning{...
                    j+NTr+length(trials)}.Bool(ii) = false; 
            end
        end  
        % Helper vector
        % SelectedMusclesInfoSpanningSel contains the indices of the
        % muscles selected for the optimization in the vector of
        % muscles spanning the degrees of freedom
        count = 1;
        for ii = 1:length(MuscleNames_Spanning)
            if find(strcmp(MuscleNames_opt,MuscleNames_Spanning{ii}))
                SelectedMusclesInfoSpanningSel{...
                    j+NTr+length(trials)}.Bool(ii) = true;       
                SelectedMusclesInfoSpanningSel{...
                    j+NTr+length(trials)}.Index(count) = ii;
                count = count+1;
            else
                SelectedMusclesInfoSpanningSel{...
                    j+NTr+length(trials)}.Bool(ii) = false; 
            end
        end  
        % Helper vector
        % SelectedMusclesInfoOptSpanning contains the indices of the muscles
        % selected for the optimization and spanning the degrees of freedom
        count = 1;
        for ii = 1:length(MuscleNames)
            if (find(strcmp(MuscleNames_opt,MuscleNames{ii}))) & ...
                    (find(strcmp(MuscleNames_Spanning,MuscleNames{ii})))
                SelectedMusclesInfoOptSpanning{...
                    j+NTr+length(trials)}.Bool(ii) = true;       
                SelectedMusclesInfoOptSpanning{...
                    j+NTr+length(trials)}.Index(count) = ii;
                count = count+1;
            else
                SelectedMusclesInfoOptSpanning{...
                    j+NTr+length(trials)}.Bool(ii) = false; 
            end
        end
    end    
    end  
    % extract number of independent muscles for EMG scaling
    indEMG2misc_all = indEMG2misc_all(:,2)';
    NindEMG2misc = length(unique(indEMG2misc_all));
    % max EMG value
    if strcmp(ISA,'yes')        
        maxEMGall = ...
            NaN(length(trials)+length(ISA_trials_id),length(EMGchannels));
    else
        maxEMGall = NaN(length(trials),length(EMGchannels));    
    end        
    for j = 1:length(trials)
        maxEMGall(j,:) = maxEMG{j};        
    end
    if strcmp(ISA,'yes')
        for j = 1:length(ISA_trials_id)
            for k = 1:size(maxEMG{j+length(trials)},2)
                maxEMGall(j+length(trials),indEMGISAinGait{j}(k)) = ...
                    maxEMG{j+length(trials)}(k);
            end
        end
    end
    maxEMGall_max.(sides{i}) = nanmax(maxEMGall,[],1);
    % normalize EMGs
    for j = 1:length(trials)
        EMGoptNorm{j+NTr} = EMGopt{j+NTr}./repmat(maxEMGall_max.(sides{i}),...
            size(EMGopt{j+NTr},1),1);        
    end
    if strcmp(ISA,'yes')
        for j = 1:length(ISA_trials_id)
            for k = 1:size(maxEMG{j+length(trials)},2)
                EMGoptNorm{j+NTr+length(trials)}(:,k) = ...
                    EMGopt{j+NTr+length(trials)}(:,k)./repmat(...
                    maxEMGall_max.(sides{i})(indEMGISAinGait{j}(k)),...
                    size(EMGopt{j+NTr+length(trials)},1),1);
            end
        end   
    end
    if strcmp(ISA,'yes')        
        NTr = length(trials)+length(ISA_trials_id);
    else
        NTr = length(trials);   
    end
end 
    
%% Solve optimal control problem and estimate muscle-tendon parameters
[params,activations,residuals,forces,scale,lMtilde,time,stats,extlMT] = ...
    parameterEstimation_ocp(JointMom,MuscMomArm,Fmax,EMGoptNorm,indEMG2misc,...
    ParametersInit,MuscTenLen,MuscTenVel,SelectedMusclesInfo,modelName,...
    muscleNames,sides,Fvparam,Fpparam,Faparam,lMtilde_ext,W,...
    clinicalReport,idx_spanning,NindEMG2misc,idx_EMGscale2misc,...
    SelectedMusclesInfoSelSpanning,SelectedMusclesInfoOptSpanning,...
    SelectedMusclesInfoSpanningSel,timeOpt,side_t,dev_p,EMGScale_min);
for i = 1:length(sides)
    side_tr = sides{i};        
    % Re-arrange parameters to includes the ones that were not optimized
    PARAMETERSall.(sides{i}).PennAng = ParametersInit.(sides{i}).PennAng;    
    PARAMETERSall.(sides{i}).OptFibLen = params.(sides{i}).OptFibLen;
    PARAMETERSall.(sides{i}).TenSlackLen = params.(sides{i}).TenSlackLen;
    PARAMETERSall.(sides{i}).MaxIsomForce = MTParameters(1,1:NMuscles/2);
    if strcmp(side_tr,'r')
        PARAMETERSall.(sides{i}).MaxIsomForce = ...
            MTParameters(1,NMuscles/2+1:NMuscles);
    end
    PARAMETERSall.(sides{i}).MaxContVelo = ...
        10*PARAMETERSall.(sides{i}).OptFibLen; 
    ParameterEstimationResults.PARAMETERS.(sides{i}) = PARAMETERSall.(sides{i});             
    ParameterEstimationResults.ParametersInit.(sides{i}) = ...
        ParametersInit.(sides{i}); 
    ParameterEstimationResults.Scale.values_s.(sides{i}) = scale.(sides{i});
end
ParameterEstimationResults.extlMT = extlMT;
for j = 1:length(side_t)
    ParameterEstimationResults.Scale.values{j} = scale.(side_t{j});  
    ParameterEstimationResults.Scale.maxEMG{j} = maxEMGall_max.(side_t{j});         
end
ParameterEstimationResults.Scale.maxEMG_s = maxEMGall_max;
for ii = 1:length(scale)
    ParameterEstimationResults.Scale.EMGchannelsOpt{ii} = ...
        EMGchannels(indEMG2misc{1}(ii,1));
end  
ParameterEstimationResults.Activations = activations;
ParameterEstimationResults.Forces = forces;
ParameterEstimationResults.Residuals = residuals;   
ParameterEstimationResults.Scale.indEMG2misc = indEMG2misc;
ParameterEstimationResults.lMtilde = lMtilde;
ParameterEstimationResults.Scale.muscles = Opt_channels;
ParameterEstimationResults.Scale.EMGchannels = EMGchannels;  
ParameterEstimationResults.Scale.EMGoptNorm = EMGoptNorm;    
ParameterEstimationResults.Scale.idx_EMGscale2misc = idx_EMGscale2misc;
ParameterEstimationResults.Time = time;
ParameterEstimationResults.idx_spanning = idx_spanning;
ParameterEstimationResults.stats = stats;    
if saveResults
    if ~(exist(pathResults,'dir')==7)
        mkdir(pathResults);
    end
    save([pathResults,'/ParameterEstimationResults'],....
        'ParameterEstimationResults');
end
end

%% Save parameters
if saveParameters
% load MT-parameters
load([pathResults,'/ParameterEstimationResults']);
% Left side
MTParameters(1,1:43) = ParameterEstimationResults.PARAMETERS.l.MaxIsomForce;
MTParameters(2,1:43) = ParameterEstimationResults.PARAMETERS.l.OptFibLen;
MTParameters(3,1:43) = ParameterEstimationResults.PARAMETERS.l.TenSlackLen;
MTParameters(4,1:43) = ParameterEstimationResults.PARAMETERS.l.PennAng;
MTParameters(5,1:43) = ParameterEstimationResults.PARAMETERS.l.MaxContVelo;
% Right side
MTParameters(1,44:86) = ParameterEstimationResults.PARAMETERS.r.MaxIsomForce;
MTParameters(2,44:86) = ParameterEstimationResults.PARAMETERS.r.OptFibLen;
MTParameters(3,44:86) = ParameterEstimationResults.PARAMETERS.r.TenSlackLen;
MTParameters(4,44:86) = ParameterEstimationResults.PARAMETERS.r.PennAng;
MTParameters(5,44:86) = ParameterEstimationResults.PARAMETERS.r.MaxContVelo;
% save parameters
save([pathResults,'/MTParameters_personalized.mat'],'MTParameters');
end

%% Plot results
if showResults 
if ~solveProblem
    load([pathResults,'ParameterEstimationResults.mat']);
    PARAMETERSall = ParameterEstimationResults.PARAMETERS;
    activations = ParameterEstimationResults.Activations;
    scale = ParameterEstimationResults.Scale.values;
    indEMG2misc = ParameterEstimationResults.Scale.indEMG2misc;
    idx_EMGscale2misc = ParameterEstimationResults.Scale.idx_EMGscale2misc;
    EMGoptNorm  = ParameterEstimationResults.Scale.EMGoptNorm;
    ParametersInit = ParameterEstimationResults.ParametersInit;
    residuals = ParameterEstimationResults.Residuals;
    lMtilde = ParameterEstimationResults.lMtilde;
    time = ParameterEstimationResults.Time;
    idx_spanning = ParameterEstimationResults.idx_spanning;
end
MuscleNames_Input_tit = {'glut-max1','glut-max2','glut-max3','glut-med1',...
    'glut-med2','glut-med3','glut-min1','glut-min2',...
    'glut-min3','add-long','add-brev','add-mag1','add-mag2',...
    'add-mag3','pectineus','iliacus','psoas','quad-fem',...
    'gemellus','piri','TFL','gracilis','semimem','semiten',...
    'bi-fem-lh','bi-fem-sh','sartorius','rectus-fem',...
    'vas-med','vas-int','vas-lat','gas-med',...
    'gas-lat','soleus','tib-post','tib-ant','ext-dig',...
    'ext-hal','flex-dig','flex-hal','per-brev','per-long','per-tert'};
% MT-parameters
figure()
subplot(2,1,1)
count = 1;
s_col = {'b','r'};
for i = 1:length(sides)
    s(count) = scatter(1:size(PARAMETERSall.(sides{i}).TenSlackLen,2),...
        PARAMETERSall.(sides{i}).TenSlackLen,s_col{i}); 
    s_leg{count} = ['Personalized: ',sides{i}];
    count = count + 1;
    hold on;
    s(count) = scatter(1:size(ParametersInit.(sides{i}).TenSlackLen,2),...
        ParametersInit.(sides{i}).TenSlackLen,s_col{i},'^');
    s_leg{count} = ['Generic: ',sides{i}];
    count = count + 1;
end
title('Tendon slack length')
legend(s,s_leg)
set(gca,'XTickLabels',MuscleNames_Input_tit)
set(gca,'XTick',1:length(MuscleNames_Input_tit))
set(gca,'XTickLabelRotation',60)
subplot(2,1,2)
count = 1;
for i = 1:length(sides)
    s(count) = scatter(1:size(PARAMETERSall.(sides{i}).OptFibLen,2),...
        PARAMETERSall.(sides{i}).OptFibLen,s_col{i}); 
    s_leg{count} = ['Personalized: ',sides{i}];
    count = count + 1;
    hold on;
    s(count) = scatter(1:size(ParametersInit.(sides{i}).OptFibLen,2),...
        ParametersInit.(sides{i}).OptFibLen,s_col{i},'^');
    s_leg{count} = ['Generic: ',sides{i}];
    count = count + 1;
end
title('Optimal fiber length')
legend(s,s_leg)
set(gca,'XTickLabels',MuscleNames_Input_tit)
set(gca,'XTick',1:length(MuscleNames_Input_tit))
set(gca,'XTickLabelRotation',60)
% Activations
for j = 1:length(activations)
    figure()
    for i = 1:size(activations{j},2)
        subplot(7,7,i)
        plot(time{j}(1:end-1),activations{j}(1:end-1,i));
        hold on
        if find(indEMG2misc{j}(:,3)==i)
            idx_channel = find(indEMG2misc{j}(:,3)==i);
            plot(time{j}(1:end-1),EMGoptNorm{j}(:,indEMG2misc{j}...
                (idx_channel,1)).*scale{j}(idx_EMGscale2misc{j}...
                (idx_channel)),'r')        
        end
        set(gca,'XTick',[]);
        temp_tit = MuscleNames_Input_tit(idx_spanning{j}==1);
        title(temp_tit{i});
    end
    suptitle('Muscle activations')
end
% Residual torques
for j = 1:length(residuals)
    figure()
    for i = 1:size(residuals{j},2)
        subplot(2,3,i)
        plot(residuals{j}(:,i));
        set(gca,'XTick',[]);
    end
    suptitle('Residuals')
end

% Fiber lengths and passive forces at extremes of the range of motion
for s = 1:length(sides)
selMInfo = ParameterEstimationResults.extlMT.SelectedMusclesInfo.(sides{s});
LMTMAX = ParameterEstimationResults.extlMT.LMTMAX.(sides{s})(:,selMInfo);
LMTMIN = ParameterEstimationResults.extlMT.LMTMIN.(sides{s})(:,selMInfo);
idx_lMopt_max = ParameterEstimationResults.extlMT.idx_lMopt_max.(sides{s});
F_lMtilde_max_opt ...
    = ParameterEstimationResults.extlMT.F_lMtilde_max_opt.(sides{s});
F_lMtilde_min_opt ...
    = ParameterEstimationResults.extlMT.F_lMtilde_min_opt.(sides{s});
NMuscleAll = length(MuscleNames_Input_tit);
NselMusc_all = length(selMInfo);
Lfibmax.(sides{s}).s = zeros(1,NMuscleAll); % all for indices
Lfibmin.(sides{s}).s = zeros(1,NselMusc_all); % only selected ones
FpetildeMax.(sides{s}).s = zeros(1,NMuscleAll); % all for indices
load Fvparam
load Fpparam
load Faparam
for m = 1:NselMusc_all
    % Upper extreme
    [Hilldiff_lMtilde_max.(sides{s}).s(m), ~, ...
        Lfibmax.(sides{s}).s(1,selMInfo(m))] = ...
        ForceEquilibrium_FtildeState_ParamEst(0,F_lMtilde_max_opt(m),0,...
        LMTMAX(m),0,PARAMETERSall.(sides{s}).MaxIsomForce(1,selMInfo(m)),...
        PARAMETERSall.(sides{s}).OptFibLen(1,selMInfo(m)),...
        PARAMETERSall.(sides{s}).TenSlackLen(1,selMInfo(m)),...
        PARAMETERSall.(sides{s}).PennAng(1,selMInfo(m)),Fvparam,Fpparam,...
        Faparam,35,0); 
    % Passive muscle force
    e0 = 0.6; kpe = 4;
    t5 = exp(kpe * (Lfibmax.(sides{s}).s(1,selMInfo(m))-0.10e1)/e0);
    FpetildeMax.(sides{s}).s(1,selMInfo(m))=((t5-0.10e1)-Fpparam(1))/Fpparam(2);
    % Lower extreme
    [Hilldiff_lMtilde_min.(sides{s}).s(m), ~, Lfibmin.(sides{s}).s(m)] = ...
        ForceEquilibrium_FtildeState_ParamEst(0,F_lMtilde_min_opt(m),0,...
        LMTMIN(m),0,PARAMETERSall.(sides{s}).MaxIsomForce(1,selMInfo(m)),...
        PARAMETERSall.(sides{s}).OptFibLen(1,selMInfo(m)),...
        PARAMETERSall.(sides{s}).TenSlackLen(1,selMInfo(m)),...
        PARAMETERSall.(sides{s}).PennAng(1,selMInfo(m)),Fvparam,Fpparam,...
        Faparam,35,0);
    if (abs(Hilldiff_lMtilde_max.(sides{s}).s(m)) > 1*10^(-6))
        disp('Issue in Hill model in extreme positions')
    end
    if (abs(Hilldiff_lMtilde_min.(sides{s}).s(m)) > 1*10^(-6))
        disp('Issue in Hill model in extreme positions')
    end
end 
end
figure()
temp = 1:length(MuscleNames_Input_tit);
s_mar = {'o','+'};
s_col_min = {'c','y'};
count = 1;
for i = 1:length(sides)
    scatter(temp(:,selMInfo),Lfibmax.(sides{i}).s(:,selMInfo),'k',s_mar{i});
    hold on;
    slM(count) = scatter(temp(idx_lMopt_max),Lfibmax.(sides{i}).s(...
        idx_lMopt_max),s_col{i},s_mar{i});
    s_lM{count} = ['lM-max: ',sides{i}];
    count = count + 1;
    slM(count) = scatter(temp(:,selMInfo),Lfibmin.(sides{i}).s,...
        s_col_min{i},s_mar{i});
    s_lM{count} = ['lM-min: ',sides{i}];
    count = count + 1;
end
scatter(temp(:,selMInfo),lMtilde_ext.max*ones(1,...
    size(temp(:,selMInfo),2)),'^','m') 
scatter(temp(:,selMInfo),lMtilde_ext.min*ones(1,...
    size(temp(:,selMInfo),2)),'v','m') 
plot([1,43],[1,1],'k--')
legend(slM,s_lM)
title('Fiber lengths in extreme positions')
set(gca,'XTickLabels',MuscleNames_Input_tit)
set(gca,'XTick',1:length(MuscleNames_Input_tit))
set(gca,'XTickLabelRotation',60)
figure()
count = 1;
for i = 1:length(sides)
    sFpe(count) = scatter(temp(:,selMInfo),...
        FpetildeMax.(sides{i}).s(:,selMInfo),s_col{i},s_mar{i}); hold on
    s_Fpe{count} = sides{i};
    count = count + 1;
end
legend(sFpe,s_Fpe)
title('Passive forces in extreme positions')
set(gca,'XTickLabels',MuscleNames_Input_tit)
set(gca,'XTick',1:length(MuscleNames_Input_tit))
set(gca,'XTickLabelRotation',60)
end
