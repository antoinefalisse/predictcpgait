% This script extracts experimental data to provide reference data for
% comparison with results of the simulations.

clear all
close all
clc

%% User settings
plotGRF = 0;
plotKinematics = 0;
plotKinetics = 0;
plotEMG = 0;
createAverageMotion = 0;
saveExperimentalData = 0;
leg = 'l'; % 'r'

subjects = {'subject1'};
subject = subjects{1};

%% Paths
pathmain = pwd;
[pathPredictiveSimulations,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathPredictiveSimulations);
pathOpenSimModel = [pathRepo,'/OpenSimModel/',subject,'/'];
pathVariousFunctions = [pathRepo,'/VariousFunctions'];
addpath(genpath(pathVariousFunctions));

%% Extract data
switch subject
    case 'subject1'
        body_mass = 3.69727 + 4.64393 + 1.43323 + 0.01986 + 0.39058 + ...
            0.04303 + 4.64393 + 1.43323 + 0.01986 + 0.39058 + 0.04303 + ...
            16.25541;        
        if strcmp(leg,'r')
            Trials.subject1 = {'Gait_10','Gait_11','Gait_12','Gait_13'};
            trials.subject1 = {'gait_10','gait_11','gait_12','gait_13'};
            timeICon.subject1 =  [2.4 ,4.43, 3   ,1.65];
            timeICoff.subject1 = [3.45,5.48, 4.01,2.64];
        else
            Trials.subject1 = {'Gait_10','Gait_12','Gait_13'};
            trials.subject1 = {'gait_10','gait_12','gait_13'};
            timeICon.subject1 =  [2.92 ,2.45 ,1.16];
            timeICoff.subject1 = [3.94 ,3.52 ,2.17];
        end
end
body_weight = body_mass*9.81;

% Kinematics
joints = {'lower_torso_RX','lower_torso_RY','lower_torso_RZ',...
    'lower_torso_TX','lower_torso_TY','lower_torso_TZ','hip_flex_l',...
    'hip_add_l','hip_rot_l','hip_flex_r','hip_add_r','hip_rot_r',...
    'knee_flex_l','knee_flex_r','ankle_flex_l','ankle_flex_r',...
    'subt_angle_l','subt_angle_r',...
    'lumbar_pitch','lumbar_roll','lumbar_yaw'};
joints_tit = {'lower-torso-RX','lower-torso-RY','lower-torso-RZ',...
    'lower-torso-TX','lower-torso-TY','lower-torso-TZ','hip-flex-l',...
    'hip-add-l','hip-rot-l','hip-flex-r','hip-add-r','hip-rot-r',...
    'knee-flex-l','knee-flex-r','ankle-flex-l','ankle-flex-r',...
    'subt-angle-l','subt-angle-r',...
    'lumbar-pitch','lumbar-roll','lumbar-yaw'};

jointi.pelvis.list  = 1; 
jointi.pelvis.rot   = 2; 
jointi.pelvis.tilt  = 3; 
jointi.pelvis.tx    = 4;
jointi.pelvis.ty    = 5;
jointi.pelvis.tz    = 6;
jointi.hip_flex.l   = 7;
jointi.hip_add.l    = 8;
jointi.hip_rot.l    = 9;
jointi.hip_flex.r   = 10;
jointi.hip_add.r    = 11;
jointi.hip_rot.r    = 12;
jointi.knee.l       = 13;
jointi.knee.r       = 14;
jointi.ankle.l      = 15;
jointi.ankle.r      = 16;
jointi.subt.l       = 17;
jointi.subt.r       = 18;
jointi.trunk.ext    = 19;
jointi.trunk.ben    = 20;
jointi.trunk.rot    = 21;

jointi.l = [jointi.hip_flex.l,jointi.hip_add.l,jointi.hip_rot.l,jointi.knee.l,jointi.ankle.l,jointi.subt.l];
jointi.r = [jointi.hip_flex.r,jointi.hip_add.r,jointi.hip_rot.r,jointi.knee.r,jointi.ankle.r,jointi.subt.r];

threshold = 20;
N = 100; % number of time points for interpolation
% Timing initial contacts & GRFs
for i = 1:length(subjects)
    count = 1;
    for j = 1:length(Trials.(subjects{i}))
        pathGRF = [pathOpenSimModel,'GRF/Gait/GRF_',Trials.(subjects{i}){j},'.mot'];
        GRF = getGRF(pathGRF);
        idx_IC.(subjects{i}){count}(1) = find(round(GRF.time,4)==round(timeICon.(subjects{i})(j),2));
        idx_IC.(subjects{i}){count}(2) = find(round(GRF.time,4)==round(timeICoff.(subjects{i})(j),2));
        time_IC.(subjects{i}){count}(1) = round(GRF.time(idx_IC.(subjects{i}){count}(1)),4);
        time_IC.(subjects{i}){count}(2) = round(GRF.time(idx_IC.(subjects{i}){count}(2)),4);         
        GRFsel.(subjects{i}){count} = GRF.val.all(idx_IC.(subjects{i}){count}(1):idx_IC.(subjects{i}){count}(2),:);                
        % interpolation over N points
        step = (GRFsel.(subjects{i}){count}(end,1)-GRFsel.(subjects{i}){count}(1,1))/(N-1);
        interval.(subjects{i}){count} = GRFsel.(subjects{i}){count}(1,1):step:GRFsel.(subjects{i}){count}(end,1);
        GRFinterp.(subjects{i})(:,:,count) = interp1(GRFsel.(subjects{i}){count}(:,1),GRFsel.(subjects{i}){count}(:,2:end),interval.(subjects{i}){count});
        GRFinterp.(subjects{i})(:,:,count) = GRFinterp.(subjects{i})(:,:,count)./(body_weight/100);
                
        % Kinematics: joint positions
        pathIK = [pathOpenSimModel,'IK/Gait/KS_',Trials.(subjects{i}){j},leg,'_MRI_extROM.mot'];
        Qs = getIK_MRI(pathIK,joints);            
        IKsel.(subjects{i}){count} = Qs.allfilt;
        IKinterp.(subjects{i})(:,:,count) = interp1(IKsel.(subjects{i}){count}(:,1),IKsel.(subjects{i}){count}(:,2:end),interval.(subjects{i}){count});                  
        % Adjust pelvis tx to start from 0       
        IKinterp.(subjects{i})(:,jointi.pelvis.tx,count) = IKinterp.(subjects{i})(:,jointi.pelvis.tx,count) - IKinterp.(subjects{i})(1,jointi.pelvis.tx,count);       
        % Convert into degrees 
        IKinterp.(subjects{i})(:,[1:3,7:end],count) = IKinterp.(subjects{i})(:,[1:3,7:end],count)*180/pi;
             
        % Kinematics: joint velocities and accelerations
        % Spline approximation
        for ii = 1:size(IKinterp.(subjects{i}),2)
            Qss.datafiltspline(ii) = ...
                spline(interval.(subjects{i}){count}',...
                IKinterp.(subjects{i})(:,ii,count));             
            [Qs_spline.(subjects{i})(:,ii,count),...
                Qdots_spline.(subjects{i})(:,ii,count),...
                Qdotdots_spline.(subjects{i})(:,ii,count)] = ...
                SplineEval_ppuval(Qss.datafiltspline(ii),...
                interval.(subjects{i}){count}',1);
        end    
        
        % Kinetics
        pathID = [pathOpenSimModel,'ID/Gait/ID_',Trials.(subjects{i}){j},leg,'_MRI_extROM.sto'];
        ID = getID_MRI(pathID,joints); 
        IDinterp.(subjects{i})(:,:,count) = interp1(ID.all(:,1),ID.all(:,2:end),interval.(subjects{i}){count}); 
        
        % EMG
        pathEMG = [pathOpenSimModel,'EMG/Gait/EMG_filt.mat'];
        load(pathEMG,'EMG_filt');
        EMG_s.all = EMG_filt.(trials.(subjects{i}){j}).data;
        EMGinterp.(subjects{i}).data(:,:,count) = interp1(EMG_s.all(:,1),EMG_s.all(:,2:end),interval.(subjects{i}){count});            
        EMGinterp.(subjects{i}).colheaders = EMG_filt.(trials.(subjects{i}){j}).colheaders(2:end);
        
        time_ICall.(subjects{i})(j,:) = time_IC.(subjects{i}){count};                      
        
        % Step width
        pathPK_r = [pathOpenSimModel,'PointKinematics/',Trials.(subjects{i}){j},leg,'/',Trials.(subjects{i}){j},leg,'_PointKinematics_origin_r_pos.sto'];
        pathPK_l = [pathOpenSimModel,'PointKinematics/',Trials.(subjects{i}){j},leg,'/',Trials.(subjects{i}){j},leg,'_PointKinematics_origin_l_pos.sto'];
        PK_r = importdata(pathPK_r);
        PK_l = importdata(pathPK_l);
        PK_2_rinterp.(subjects{i})(:,count) = interp1(PK_r.data(:,1),PK_r.data(:,strcmp(PK_r.colheaders,'state_2')),interval.(subjects{i}){count});    
        PK_2_linterp.(subjects{i})(:,count) = interp1(PK_l.data(:,1),PK_l.data(:,strcmp(PK_l.colheaders,'state_2')),interval.(subjects{i}){count});    
        PK_diff_2_interp.(subjects{i})(:,count) = mean(abs(PK_2_rinterp.(subjects{i})(:,count)-PK_2_linterp.(subjects{i})(:,count)));
        
        % Stride length    
        PK_0_rinterp.(subjects{i})(:,count) = interp1(PK_r.data(:,1),PK_r.data(:,strcmp(PK_r.colheaders,'state_0')),interval.(subjects{i}){count});    
        PK_1_rinterp.(subjects{i})(:,count) = interp1(PK_r.data(:,1),PK_r.data(:,strcmp(PK_r.colheaders,'state_1')),interval.(subjects{i}){count});  
        PK_2_rinterp.(subjects{i})(:,count) = interp1(PK_r.data(:,1),PK_r.data(:,strcmp(PK_r.colheaders,'state_2')),interval.(subjects{i}){count});    
        sLength_interp.(subjects{i})(count) = sqrt((PK_0_rinterp.(subjects{i})(end,count)-PK_0_rinterp.(subjects{i})(1,count))^2 + ...
            (PK_1_rinterp.(subjects{i})(end,count)-PK_1_rinterp.(subjects{i})(1,count))^2 + ...
            (PK_2_rinterp.(subjects{i})(end,count)-PK_2_rinterp.(subjects{i})(1,count))^2);   
            
        count =  count + 1;   
    end 
    
    GRF_ref.(subjects{i}).all = GRFinterp.(subjects{i});
    GRF_ref.(subjects{i}).mean = mean(GRFinterp.(subjects{i}),3);
    GRF_ref.(subjects{i}).std = std(GRFinterp.(subjects{i}),[],3);
    GRF_ref.(subjects{i}).colheaders = {'Fore-aft-r','Vertical-r','Lateral-r','Fore-aft-l','Vertical-l','Lateral-l'};    
    Q_ref.(subjects{i}).Qs.all = Qs_spline.(subjects{i});
    Q_ref.(subjects{i}).Qs.mean = mean(Qs_spline.(subjects{i}),3);
    Q_ref.(subjects{i}).Qs.std = std(Qs_spline.(subjects{i}),[],3);
    Q_ref.(subjects{i}).Qs.colheaders = joints;
    Q_ref.(subjects{i}).Qdots.all = Qdots_spline.(subjects{i});
    Q_ref.(subjects{i}).Qdots.mean = mean(Qdots_spline.(subjects{i}),3);
    Q_ref.(subjects{i}).Qdots.std = std(Qdots_spline.(subjects{i}),[],3);
    Q_ref.(subjects{i}).Qdots.colheaders = joints;
    Q_ref.(subjects{i}).Qdotdots.all = Qdotdots_spline.(subjects{i});
    Q_ref.(subjects{i}).Qdotdots.mean = mean(Qdotdots_spline.(subjects{i}),3);
    Q_ref.(subjects{i}).Qdotdots.std = std(Qdotdots_spline.(subjects{i}),[],3);   
    Q_ref.(subjects{i}).Qdotdots.colheaders = joints;
    Torques_ref.(subjects{i}).all = IDinterp.(subjects{i});
    Torques_ref.(subjects{i}).mean = mean(IDinterp.(subjects{i}),3);
    Torques_ref.(subjects{i}).std = std(IDinterp.(subjects{i}),[],3);
    Torques_ref.(subjects{i}).colheaders = joints;
    EMG_ref.(subjects{i}).all = EMGinterp.(subjects{i}).data;
    EMG_ref.(subjects{i}).mean = nanmean(EMGinterp.(subjects{i}).data,3);
    EMG_ref.(subjects{i}).std = nanstd(EMGinterp.(subjects{i}).data,[],3);
    EMG_ref.(subjects{i}).colheaders = EMGinterp.(subjects{i}).colheaders;
    StepWidth_ref.(subjects{i}).all = PK_diff_2_interp.(subjects{i});
    StepWidth_ref.(subjects{i}).mean = mean(PK_diff_2_interp.(subjects{i}),2);
    StepWidth_ref.(subjects{i}).std = std(PK_diff_2_interp.(subjects{i}),[],2);
    StrideLength_ref.(subjects{i}).all = sLength_interp.(subjects{i});
    StrideLength_ref.(subjects{i}).mean = mean(sLength_interp.(subjects{i}),2);
    StrideLength_ref.(subjects{i}).std = std(sLength_interp.(subjects{i}),[],2);    
    ExperimentalData.Q = Q_ref;
    ExperimentalData.GRFs = GRF_ref;
    ExperimentalData.Torques = Torques_ref;
    ExperimentalData.EMG = EMG_ref;
    ExperimentalData.StepWidth = StepWidth_ref;
    ExperimentalData.StrideLength = StrideLength_ref;   
    
    if saveExperimentalData
        pathExperimentalData = [pathPredictiveSimulations,'/ExperimentalData/'];
        save([pathExperimentalData,'ExperimentalData_',leg],'ExperimentalData');
    end              
    
    if plotGRF
        figure()
        if strcmp(leg,'r')        
        for c = 1:size(GRFinterp.(subjects{i}),3)
            plot(GRFinterp.(subjects{i})(:,1,c),'b','linewidth',0.5);
            hold on;
            plot(GRFinterp.(subjects{i})(:,2,c),'r','linewidth',0.5);
            plot(GRFinterp.(subjects{i})(:,3,c),'m','linewidth',0.5);        
        end
        plot(GRF_ref.(subjects{i}).mean(:,1),'b','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,1)+GRF_ref.(subjects{i}).std(:,1),'b--','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,1)-GRF_ref.(subjects{i}).std(:,1),'b--','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,2),'r','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,2)+GRF_ref.(subjects{i}).std(:,2),'r--','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,2)-GRF_ref.(subjects{i}).std(:,2),'r--','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,3),'m','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,3)+GRF_ref.(subjects{i}).std(:,3),'m--','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,3)-GRF_ref.(subjects{i}).std(:,3),'m--','linewidth',2)
        title('GRF-Right','Fontsize',20);  
        else
        for c = 1:size(GRFinterp.(subjects{i}),3)
            plot(GRFinterp.(subjects{i})(:,4,c),'b','linewidth',0.5);
            hold on;
            plot(GRFinterp.(subjects{i})(:,5,c),'r','linewidth',0.5);
            plot(GRFinterp.(subjects{i})(:,6,c),'m','linewidth',0.5);        
        end
        plot(GRF_ref.(subjects{i}).mean(:,4),'b','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,4)+GRF_ref.(subjects{i}).std(:,4),'b--','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,4)-GRF_ref.(subjects{i}).std(:,4),'b--','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,5),'r','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,5)+GRF_ref.(subjects{i}).std(:,5),'r--','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,5)-GRF_ref.(subjects{i}).std(:,5),'r--','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,6),'m','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,6)+GRF_ref.(subjects{i}).std(:,6),'m--','linewidth',2);
        plot(GRF_ref.(subjects{i}).mean(:,6)-GRF_ref.(subjects{i}).std(:,6),'m--','linewidth',2)
        title('GRF-Left','Fontsize',20);  
        end
    end

    if plotKinematics
        figure()
        for jj = 1:size(Q_ref.(subjects{i}).Qs.mean,2)
            subplot(4,6,jj)
            for c = 1:size(IKinterp.(subjects{i}),3)
                plot(Qs_spline.(subjects{i})(:,jj,c),'k','linewidth',0.5); hold on
            end
            plot(Q_ref.(subjects{i}).Qs.mean(:,jj),'k','linewidth',2); hold on
            plot(Q_ref.(subjects{i}).Qs.mean(:,jj)+Q_ref.(subjects{i}).Qs.std(:,jj),'k--','linewidth',2);
            plot(Q_ref.(subjects{i}).Qs.mean(:,jj)-Q_ref.(subjects{i}).Qs.std(:,jj),'k--','linewidth',2);
            title(joints_tit{jj},'Fontsize',16);
        end  
        s = suptitle('Joint positions');
        set(s,'Fontsize',20);
        
        figure()
        for jj = 1:size(Q_ref.(subjects{i}).Qdots.mean,2)
            subplot(4,6,jj)
            for c = 1:size(IKinterp.(subjects{i}),3)
                plot(Qdots_spline.(subjects{i})(:,jj,c),'k','linewidth',0.5); hold on
            end
            plot(Q_ref.(subjects{i}).Qdots.mean(:,jj),'k','linewidth',2); hold on
            plot(Q_ref.(subjects{i}).Qdots.mean(:,jj)+Q_ref.(subjects{i}).Qdots.std(:,jj),'k--','linewidth',2);
            plot(Q_ref.(subjects{i}).Qdots.mean(:,jj)-Q_ref.(subjects{i}).Qdots.std(:,jj),'k--','linewidth',2);
            title(joints_tit{jj},'Fontsize',16);
        end  
        s = suptitle('Joint velocities');
        set(s,'Fontsize',20);
        
        figure()
        for jj = 1:size(Q_ref.(subjects{i}).Qdotdots.mean,2)
            subplot(4,6,jj)
            for c = 1:size(IKinterp.(subjects{i}),3)
                plot(Qdotdots_spline.(subjects{i})(:,jj,c),'k','linewidth',0.5); hold on
            end
            plot(Q_ref.(subjects{i}).Qdotdots.mean(:,jj),'k','linewidth',2); hold on
            plot(Q_ref.(subjects{i}).Qdotdots.mean(:,jj)+Q_ref.(subjects{i}).Qdotdots.std(:,jj),'k--','linewidth',2);
            plot(Q_ref.(subjects{i}).Qdotdots.mean(:,jj)-Q_ref.(subjects{i}).Qdotdots.std(:,jj),'k--','linewidth',2);
            title(joints_tit{jj},'Fontsize',16);
        end  
        s = suptitle('Joint accelerations');
        set(s,'Fontsize',20);
    end
    
    if plotKinetics
        figure()
        for jj = 1:size(jointi.(leg),2)
            subplot(2,3,jj)
            for c = 1:size(IDinterp.(subjects{i}),3)
                plot(IDinterp.(subjects{i})(:,jointi.(leg)(jj),c),'k','linewidth',0.5); hold on
            end
            plot(Torques_ref.(subjects{i}).mean(:,jointi.(leg)(jj)),'k','linewidth',2); hold on
            plot(Torques_ref.(subjects{i}).mean(:,jointi.(leg)(jj))+Torques_ref.(subjects{i}).std(:,jointi.(leg)(jj)),'k--','linewidth',2);
            plot(Torques_ref.(subjects{i}).mean(:,jointi.(leg)(jj))-Torques_ref.(subjects{i}).std(:,jointi.(leg)(jj)),'k--','linewidth',2);
            title(joints_tit{jointi.(leg)(jj)},'Fontsize',16);
        end 
        s = suptitle('Joint kinetics');
        set(s,'Fontsize',20);
    end   
     
    if plotEMG
        figure()
        for jj = 1:size(EMG_ref.(subjects{i}).mean,2)
            subplot(4,4,jj)
            for c = 1:size(EMGinterp.(subjects{i}),3)
                plot(EMGinterp.(subjects{i}).data(:,jj,c),'k','linewidth',0.5); hold on
            end
            plot(EMG_ref.(subjects{i}).mean(:,jj),'r','linewidth',2); hold on
            plot(EMG_ref.(subjects{i}).mean(:,jj)+EMG_ref.(subjects{i}).std(:,jj),'b--','linewidth',2);
            plot(EMG_ref.(subjects{i}).mean(:,jj)-EMG_ref.(subjects{i}).std(:,jj),'b--','linewidth',2);
            title(EMG_ref.(subjects{i}).colheaders{jj},'Fontsize',14);
            ylim([0 inf])
        end 
        s = suptitle('EMG');
        set(s,'Fontsize',20);
    end   
    
    % Get average time
    time_ICalls.(subjects{i}).all(1,:) = (time_ICall.(subjects{i})(:,2) - time_ICall.(subjects{i})(:,1))';
    time_ICalls.(subjects{i}).mean = mean(time_ICalls.(subjects{i}).all);
    time_ICalls.(subjects{i}).std = std(time_ICalls.(subjects{i}).all);
    
    step = time_ICalls.(subjects{i}).mean/(N-1);
    time.(subjects{i}).mean = (0:step:time_ICalls.(subjects{i}).mean)';
    
    % Get average distance traveled (pelvis horizontal displacement)
    dist_travelled.(subjects{i}).all(1,:) = Qs_spline.(subjects{i})(end,4,:);
    dist_travelled.(subjects{i}).mean = mean(dist_travelled.(subjects{i}).all,2);
    dist_travelled.(subjects{i}).std = std(dist_travelled.(subjects{i}).all,[],2);
    
    % Get average speed
    speed.(subjects{i}).all(1,:) = dist_travelled.(subjects{i}).all./time_ICalls.(subjects{i}).all;
    speed.(subjects{i}).mean = mean(speed.(subjects{i}).all,2);
    speed.(subjects{i}).std = std(speed.(subjects{i}).all,[],2); 
    
    if createAverageMotion
        % Create fake motion from the average data
        JointAngle.labels = {'time'};
        JointAngle.labels = [JointAngle.labels,joints];
        JointAngle.data = [time.(subjects{i}).mean,Q_ref.(subjects{i}).Qs.mean];
        namescript = ['walking_average_',leg];
        filenameJointAngles = [pathOpenSimModel,'IK/Gait/KS_Gait_average_',leg,'.mot'];
        write_motionFile(JointAngle, filenameJointAngles)        
    end
end
