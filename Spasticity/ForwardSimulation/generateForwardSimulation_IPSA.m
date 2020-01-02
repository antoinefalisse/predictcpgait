% This script generates the EMG-driven forward simulations for the passive
% motions (IPSA)    

function generateForwardSimulation_IPSA(IK_path,ID_path,...
    MuscleAnalysis_path,EMG_path_agonist,EMG_path_antagonist,Misc)

%% Filter parameters 
% Default low-pass filter: Butterworth order: 6, Cutoff frequency: 6Hz
% Inverse Dynamics
if ~isfield(Misc,'f_cutoff_ID') || isempty(Misc.f_cutoff_ID)
    Misc.f_cutoff_ID=6;
end
if ~isfield(Misc,'f_order_ID') || isempty(Misc.f_order_ID)
    Misc.f_order_ID=6;
end
% Muscle-tendon lengths
if ~isfield(Misc,'f_cutoff_lMT') || isempty(Misc.f_cutoff_lMT)
    Misc.f_cutoff_lMT=6;
end
if ~isfield(Misc,'f_order_lMT') || isempty(Misc.f_order_lMT)
    Misc.f_order_lMT=6;
end
% Moment arms
if ~isfield(Misc,'f_cutoff_dM') || isempty(Misc.f_cutoff_dM)
    Misc.f_cutoff_dM=6;
end
if ~isfield(Misc,'f_order_dM') || isempty(Misc.f_order_dM)
    Misc.f_order_dM=6;
end
% Inverse Kinematics
if ~isfield(Misc,'f_cutoff_IK') || isempty(Misc.f_cutoff_IK)
    Misc.f_cutoff_IK=6;
end
if ~isfield(Misc,'f_order_IK') || isempty(Misc.f_order_IK)
    Misc.f_order_IK=6;
end
% Mesh Frequency
if ~isfield(Misc,'Mesh_Frequency') || isempty(Misc.Mesh_Frequency)
    Misc.Mesh_Frequency=300;
end

%%  Get muscle information 
Nsegment = Misc.Nsegment;
% Constant parameters
switch Misc.joint
    case 'knee'
        Muscles_spas = [1 8 9]; % hamstrings (BFLH, SM, ST)
        Muscles_antagonist = 6;
    case 'ankle'
        Muscles_spas = [1 2];   % gastrocnemii (GM,GL)
        Muscles_antagonist = 7;
end
% Path to muscle analysis results
Misc.MuscleAnalysisPath=MuscleAnalysis_path;
% Pre-allocation
DatStore = struct('MS',[]);
emg1processedstr = struct('MS',[]);
emg2processedstr = struct('MS',[]);
EMG_ordered = struct('MS',[]);
for ms = 1:Nsegment    
    % Get # dofs, lMT, MA
    [~,Misc.trialName(ms).MS,~]=fileparts(IK_path(ms).MS);
    if ~isfield(Misc,'MuscleNames_Input')||isempty(Misc.MuscleNames_Input)    
        Misc=getMusclesDofs_IPSA(Misc);
    end
    DatStore(ms).MS = getMuscleInfo_IPSA(IK_path(ms).MS,ID_path(ms).MS,...
        Misc,ms);
    % Get MT-parameters
    DatStore(ms).MS.params  = Misc.params_scaled;
    % Get EMG
    temp1 = load(EMG_path_agonist(ms).MS);
    temp2 = load(EMG_path_antagonist(ms).MS);
    emg1processedstr(ms).MS = ...
        temp1.emg1processednorm(Misc.range(ms).MS.range,:);
    emg2processedstr(ms).MS = ...
        temp2.emg2processednorm(Misc.range(ms).MS.range,:);  
    EMG_ordered(ms).MS = zeros(...
        size(emg1processedstr(ms).MS,1),DatStore(ms).MS.nMuscles);   
    EMG_ordered(ms).MS(:,Muscles_spas) = emg1processedstr(ms).MS;
    EMG_ordered(ms).MS(:,Muscles_antagonist) = emg2processedstr(ms).MS;
    % Remove beginning and end of motion to account for bad EMG 
    DatStore(ms).MS.time = DatStore(ms).MS.time(10:end-10,:);
    DatStore(ms).MS.dM = DatStore(ms).MS.dM(10:end-10,:,:);
    DatStore(ms).MS.LMT = DatStore(ms).MS.LMT(10:end-10,:);
    DatStore(ms).MS.T_exp = DatStore(ms).MS.T_exp(10:end-10,:);
    EMG_ordered(ms).MS = EMG_ordered(ms).MS(10:end-10,:);  
    % Remove potential negative EMG due to low-pass filtering
    EMG_ordered(ms).MS(EMG_ordered(ms).MS<0.005)=0.005;
end

%% Perform EMG-driven forward simulations
if ~(exist([Misc.savefolder],'dir')==7)
    mkdir(Misc.savefolder);
end
if (exist([Misc.savefolder,'\outputGPOPSFTtilde_',Misc.joint,'.mat'],...
        'file')==2)
    load([Misc.savefolder,'\outputGPOPSFTtilde_',Misc.joint]);
    load([Misc.savefolder,'\HillGPOPSFTtilde_',Misc.joint]);
    load([Misc.savefolder,'\M_musclesGPOPSFTtilde_',Misc.joint]);
    load([Misc.savefolder,'\M_kneeGPOPSFTtilde_',Misc.joint]);
    load([Misc.savefolder,'\FTGPOPSFTtilde_',Misc.joint]);
    load([Misc.savefolder,'\FtildeGPOPSFTtilde_',Misc.joint]);
    load([Misc.savefolder,'\dFtildeGPOPSFTtilde_',Misc.joint]);
    load([Misc.savefolder,'\lMtildeGPOPSFTtilde_',Misc.joint]);
    load([Misc.savefolder,'\vMtildeGPOPSFTtilde_',Misc.joint]);
else
    % Pre-allocation
    trialtemp = ['segment_' int2str(Misc.segment_sel_all(1))];
    outputGPOPS = struct(trialtemp,[]);
    HillGPOPS = struct(trialtemp,[]);
    M_musclesGPOPS = struct(trialtemp,[]);
    M_kneeGPOPS = struct(trialtemp,[]);
    FTGPOPS = struct(trialtemp,[]);
    FtildeGPOPS = struct(trialtemp,[]);
    dFtildeGPOPS = struct(trialtemp,[]);
    lMtildeGPOPS = struct(trialtemp,[]);
    vMtildeGPOPS = struct(trialtemp,[]);    
end
for ms = 1:Nsegment
    trial = ['segment_' int2str(Misc.segment_sel_all(ms))];
    % EMG-driven forward simulation
    [outputGPOPS.(trial),HillGPOPS.(trial),M_musclesGPOPS.(trial),...
        M_kneeGPOPS.(trial),FTGPOPS.(trial),FtildeGPOPS.(trial),...
        dFtildeGPOPS.(trial),~,lMtildeGPOPS.(trial),...
        vMtildeGPOPS.(trial)] = ...
        ForwardSimulation_main(DatStore(ms).MS.time,EMG_ordered(ms).MS,...
        DatStore(ms).MS.params,DatStore(ms).MS.LMT,...
        DatStore(ms).MS.dM(:,1,:),Misc.pathMuscleModel);
end    
save([Misc.savefolder,'\outputGPOPSFTtilde_',Misc.joint],'outputGPOPS');
save([Misc.savefolder,'\HillGPOPSFTtilde_',Misc.joint],'HillGPOPS');
save([Misc.savefolder,'\M_musclesGPOPSFTtilde_',Misc.joint],...
    'M_musclesGPOPS');
save([Misc.savefolder,'\M_kneeGPOPSFTtilde_',Misc.joint],'M_kneeGPOPS');
save([Misc.savefolder,'\FTGPOPSFTtilde_',Misc.joint],'FTGPOPS');
save([Misc.savefolder,'\FtildeGPOPSFTtilde_',Misc.joint],'FtildeGPOPS');
save([Misc.savefolder,'\dFtildeGPOPSFTtilde_',Misc.joint],'dFtildeGPOPS');
save([Misc.savefolder,'\lMtildeGPOPSFTtilde_',Misc.joint],'lMtildeGPOPS');
save([Misc.savefolder,'\vMtildeGPOPSFTtilde_',Misc.joint],'vMtildeGPOPS');   

end