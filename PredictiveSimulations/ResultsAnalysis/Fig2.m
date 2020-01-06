% This script reproduces Fig 2
% Author: Antoine Falisse
% Date: 1/7/2019

clear all
close all
clc

%% User settings
% Selected trials
ww = [3,1];
subject = 'subject1';
ocp_path = 'PredSimOCP';

%% Paths
pathmain = pwd;
[pathPredictiveSimulations,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathPredictiveSimulations);
pathSettings = [pathPredictiveSimulations,'/Settings'];
addpath(genpath(pathSettings));
pathVariousFunctions = [pathRepo,'/VariousFunctions'];
addpath(genpath(pathVariousFunctions));

%% Fixed settings
switch subject
    case 'subject1'
        body_mass = 3.69727 + 4.64393 + 1.43323 + 0.01986 + 0.39058 + ...
            0.04303 + 4.64393 + 1.43323 + 0.01986 + 0.39058 + 0.04303 + ...
            16.25541;      
end
body_weight = body_mass*9.81;

%% Load results
% Pre-allocation structures
Qs_r            = struct('m',[]);
Acts_noSpas_r   = struct('m',[]);
Acts_r          = struct('m',[]);
Ts_r            = struct('m',[]);
GRFs_r          = struct('m',[]);
COT             = struct('m',[]);
Qs_toTrack      = struct('m',[]);
Ts_toTrack      = struct('m',[]);
GRFs_toTrack    = struct('m',[]);
StrideLength    = struct('m',[]);
corrCP          = struct('m',[]);
corrTD          = struct('m',[]);
corrCPAll       = struct('m',[]);
corrTDAll       = struct('m',[]);
NPhases         = struct('m',[]);
TO              = struct('r',[]);
HS              = struct('l',[]);
PredSimOCP_settings
pathresults = [pathPredictiveSimulations,'/Results/',ocp_path];
load([pathresults,'/ResultsPredSim.mat']);
threshold = 30;
for k = 1:length(ww) 
    N = settings(ww(k),3);   % number of mesh intervals
    NSyn = settings(ww(k),5);   % number of synergies per leg
    NPhases(ww(k)).m = 1;
    for mpi = 1:NPhases(ww(k)).m
        Qs_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Qs_r(mpi).mp;
        Ts_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Ts_r(mpi).mp;
        GRFs_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).GRFs_r(mpi).mp;
        Qs_toTrack(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Qs_toTrack(mpi).mp;
        Ts_toTrack(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Ts_toTrack(mpi).mp;
        GRFs_toTrack(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).GRFs_toTrack(mpi).mp;
        Acts_noSpas_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_noSpas_r(mpi).mp;
        Acts_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_r(mpi).mp;
        COT(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).COT(mpi).mp;         
        StrideLength(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).StrideLength(mpi).mp; 
        % Extract stance/swing
        TO(mpi).r(k) = find(GRFs_r(ww(k)).m(mpi).mp(:,2)<threshold,1,'first');
        TO(mpi).l(k) = find(GRFs_r(ww(k)).m(mpi).mp(:,5)<threshold,1,'first');
        tempDiff = false(N,1);    
        tempDiff(2:end) = diff(GRFs_r(ww(k)).m(mpi).mp(:,5)>threshold) & ...
            diff(GRFs_r(ww(k)).m(mpi).mp(:,5))>0;
        HS(mpi).l(k) = find(tempDiff);        
    end
end

%% Common settings for plots
label_fontsize  = 12;
line_linewidth  = 3;
% Colors
color_all(1,:) = [253,174,97]./255; % Orange
color_all(2,:) = [171,221,164]/255; % Green
color_all(3,:) = [72,133,237]/255;  % Blue

%% Experimental data TODO
pathExperimentalData = [pathPredictiveSimulations,'/ExperimentalData'];
load([pathExperimentalData,'/ExperimentalData_r.mat'],'ExperimentalData');

%% Plot joint angles
ref_Qs = ExperimentalData.Q;
pos_Qs = 1:5;
ylim_Qs = [-30,30;-30,30;-30,30;-2,2;0.5,1;-1,1;-10,50;-25,25;-25,25;-10,50;...
    -25,25;-25,25;0,90;0,90;-30,40;-30,40;-50,20;-50,20];
idx_Qs = [10,11,12,14,16]; 
refNames_Qs = {'lower-torso-RX','lower-torso-RY','lower-torso-RZ',...
    'lower-torso-TX','lower-torso-TY','lower-torso-TZ','Hip flexion',...
    'Hip adduction','Hip rotation','Hip flexion','Hip adduction',...
    'Hip rotation','Knee flexon','Knee flexion','Ankle flexion',...
    'Ankle flexion','subt-angle-l','subt-angle-r'};
names_Qs = {'lower_torso_RX','lower_torso_RY','lower_torso_RZ',...
    'lower_torso_TX','lower_torso_TY','lower_torso_TZ','hip_flex_l',...
    'hip_add_l','hip_rot_l','hip_flex_r','hip_add_r','hip_rot_r',...
    'knee_flex_l','knee_flex_r','ankle_flex_l','ankle_flex_r',...
    'subt_angle_l','subt_angle_r','lumbar_pitch','lumbar_roll','lumbar_yaw'};
NumTicks_Qs = 2;

figure()
for mpi = 1:NPhases(ww(1)).m
count = 1;
for i = 1:length(idx_Qs)
    p = gobjects(1,length(ww)+1);
    subplot(5,5,pos_Qs(i))      
    % Experimental data of the CP child.
    idx_jref = strcmp(ref_Qs.(subject).Qs.colheaders,names_Qs{idx_Qs(i)});
    meanPlusSTD = ref_Qs.(subject).Qs.mean(:,idx_jref) + 2*ref_Qs.(subject).Qs.std(:,idx_jref);
    meanMinusSTD = ref_Qs.(subject).Qs.mean(:,idx_jref) - 2*ref_Qs.(subject).Qs.std(:,idx_jref);          
    stepQ = (size(Qs_r(ww(1)).m(mpi).mp,1)-1)/(size(meanPlusSTD,1)-1);
    intervalQ = 1:stepQ:size(Qs_r(ww(1)).m(mpi).mp,1);
    sampleQ = 1:size(Qs_r(ww(1)).m(mpi).mp,1);
    meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
    meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);
    meanInterp = interp1(intervalQ,ref_Qs.(subject).Qs.mean(:,idx_jref),sampleQ);
    hold on
    x = Qs_toTrack(ww(1)).m(mpi).mp(:,1)';
    yy = fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'-','linewidth',2);
    yy.EdgeColor = [192,192,192]/255;
    yy.FaceColor = [192,192,192]/255;    
    p(length(ww)+2) = yy; 
    % Experimental data of the TD child.
    Qs_toTrack_deg = Qs_toTrack(ww(1)).m(mpi).mp;
    Qs_toTrack_deg(:,[2:4,8:end]) = Qs_toTrack_deg(:,[2:4,8:end])*180/pi;    
    p(length(ww)+1) = plot(Qs_toTrack_deg(:,1),Qs_toTrack_deg(:,idx_Qs(i)+1),'color','k','linewidth',line_linewidth);        
    for k = 1:length(ww)
        % Simulation results
        p(k) = plot(Qs_toTrack(ww(k)).m(mpi).mp(:,1),Qs_r(ww(k)).m(mpi).mp(:,idx_Qs(i)),'color',color_all(k,:),'linewidth',line_linewidth);    
        hold on;        
        [corrTD(ww(k)).m(mpi).mp.r2.angles.all(i),corrTD(ww(k)).m(mpi).mp.rmse.angles.all(i)] = ...
            rsquare(Qs_toTrack_deg(:,idx_Qs(i)+1),Qs_r(ww(k)).m(mpi).mp(:,idx_Qs(i))); 
        corrTDAll(mpi).mp.r2.angles.all(count,i) = corrTD(ww(k)).m(mpi).mp.r2.angles.all(i);
        corrTDAll(mpi).mp.rmse.angles.all(count,i) = corrTD(ww(k)).m(mpi).mp.rmse.angles.all(i);
        [corrCP(ww(k)).m(mpi).mp.r2.angles.all(i),corrCP(ww(k)).m(mpi).mp.rmse.angles.all(i)] = ...
            rsquare(meanInterp',Qs_r(ww(k)).m(mpi).mp(:,idx_Qs(i))); 
        corrCPAll(mpi).mp.r2.angles.all(count,i) = corrCP(ww(k)).m(mpi).mp.r2.angles.all(i);
        corrCPAll(mpi).mp.rmse.angles.all(count,i) = corrCP(ww(k)).m(mpi).mp.rmse.angles.all(i);
        count = count+1;
        % Stance/swing
        temp_name = names_Qs{idx_Qs(i)};
        if strcmp(temp_name(end-1:end),'_l')
            plot([Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1)],...
                [ylim_Qs(idx_Qs(i),1) ylim_Qs(idx_Qs(i),2)],...
                'color',color_all(k,:),'linewidth',1);
            plot([Qs_toTrack(ww(k)).m(mpi).mp(HS(mpi).l(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(HS(mpi).l(k),1)],...
                [ylim_Qs(idx_Qs(i),1) ylim_Qs(idx_Qs(i),2)],...
                'color',color_all(k,:),'linestyle','--','linewidth',1);            
        else
            plot([Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).r(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).r(k),1)],...
                [ylim_Qs(idx_Qs(i),1) ylim_Qs(idx_Qs(i),2)],...
                'color',color_all(k,:),'linewidth',1);
        end
    end 
    count = 1;
    % Plot settings 
    set(gca,'Fontsize',label_fontsize);    
    title(refNames_Qs{idx_Qs(i)},'Fontsize',label_fontsize);  
    % Y-axis
    ylim([ylim_Qs(idx_Qs(i),1) ylim_Qs(idx_Qs(i),2)]);
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks_Qs));       
    if i == 1
        ylabel('Angles (°)','Fontsize',label_fontsize);
    end      
    % X-axis
    xlim([Qs_toTrack_deg(1,1),Qs_toTrack_deg(end,1)])
    set(gca,'XTick',[]);
    box off;
end 
end
l = legend(p,{'Generic parameters','Personalized parameters',...
    'Reference TD child','Experimental CP child'});
set(l,'Fontsize',20)
% Analyse RMS
corrCPAll.mp.rmse.angles.diff = diff(corrCPAll.mp.rmse.angles.all,1,1);
corrCPAll.mp.rmse.angles.diff_per = corrCPAll.mp.rmse.angles.diff./...
    corrCPAll.mp.rmse.angles.all(2:end,:)*100;

%% Plot ground reaction forces
ref_GRF = ExperimentalData.GRFs;
pos_GRF = 21:23;
ylim_GRF = [-50,50;0 200;-25 25;-50,50;0 200;-25 25];
NumTicks_GRF = 2;
names_GRF = {'Fore-aft','Vertical','Lateral','Fore-aft L','Vertical L','Lateral L'};

for mpi = 1:NPhases(ww(1)).m
temp = GRFs_toTrack(ww(1)).m(mpi).mp(:,[3,6]);
temp(temp<0) = 0;
GRFs_toTrack(ww(1)).m(mpi).mp(:,[3,6])= temp;
count = 1;
for i = 1:length(names_GRF)/2
    subplot(5,5,pos_GRF(i))
    % Experimental data of the CP child.
    meanPlusSTD = ref_GRF.(subject).mean(:,i) + 2*ref_GRF.(subject).std(:,i);    
    meanMinusSTD = ref_GRF.(subject).mean(:,i) - 2*ref_GRF.(subject).std(:,i);   
    stepGRF = (size(GRFs_r(ww(1)).m(mpi).mp,1)-1)/(size(meanPlusSTD,1)-1);
    intervalGRF = 1:stepGRF:size(GRFs_r(ww(1)).m(mpi).mp,1);
    sampleGRF = 1:size(GRFs_r(ww(1)).m(mpi).mp,1);
    meanPlusSTD = interp1(intervalGRF,meanPlusSTD,sampleGRF);
    meanMinusSTD = interp1(intervalGRF,meanMinusSTD,sampleGRF);
    hold on
    x = GRFs_toTrack(ww(1)).m(mpi).mp(:,1)';
    yy = fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'-','linewidth',2);
    yy.EdgeColor = [192,192,192]/255;
    yy.FaceColor = [192,192,192]/255; 
    % Experimental data of the TD child.
    plot(GRFs_toTrack(ww(k)).m(mpi).mp(:,1),GRFs_toTrack(ww(1)).m(mpi).mp(:,i+1)./(body_weight/100),...
        'color','k','linewidth',line_linewidth);
    for k = 1:length(ww)
        % Simulation results
        plot(GRFs_toTrack(ww(k)).m(mpi).mp(:,1),GRFs_r(ww(k)).m(mpi).mp(:,i)./(body_weight/100),...
            'color',color_all(k,:),'linewidth',line_linewidth);
        hold on;   
        [corrTD(ww(k)).m(mpi).mp.r2.GRF.all(i),corrTD(ww(k)).m(mpi).mp.rmse.GRF.all(i)] = ...
            rsquare(GRFs_toTrack(ww(k)).m(mpi).mp(:,i+1)./(body_weight/100),GRFs_r(ww(k)).m(mpi).mp(:,i)./(body_weight/100));  
        corrTDAll(mpi).mp.r2.GRF.all(count,i) = corrTD(ww(k)).m(mpi).mp.r2.GRF.all(i);
        corrTDAll(mpi).mp.rmse.GRF.all(count,i) = corrTD(ww(k)).m(mpi).mp.rmse.GRF.all(i);
        count = count+1;
        % Stance/swing
        if i > 3
            plot([Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1)],...
                [ylim_GRF(i,1) ylim_GRF(i,2)],...
                'color',color_all(k,:),'linewidth',1);
            plot([Qs_toTrack(ww(k)).m(mpi).mp(HS(mpi).l(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(HS(mpi).l(k),1)],...
                [ylim_GRF(i,1) ylim_GRF(i,2)],...
                'color',color_all(k,:),'linestyle','--','linewidth',1);            
        else
            plot([Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).r(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).r(k),1)],...
                [ylim_GRF(i,1) ylim_GRF(i,2)],...
                'color',color_all(k,:),'linewidth',1);
        end
    end
    count = 1;    
    % Plot settings 
    set(gca,'Fontsize',label_fontsize);    
    title(names_GRF{i},'Fontsize',label_fontsize);  
    % Y-axis
    ylim([ylim_GRF(i,1) ylim_GRF(i,2)]);
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks_GRF));       
    if i == 1
        ylabel('GRFs (%BW)','Fontsize',label_fontsize);
    end      
    % X-axis
    xlim([Qs_toTrack_deg(1,1),Qs_toTrack_deg(end,1)])
    set(gca,'XTick',[]);
    box off;   
end
end

%% Plot joint torques
ref_Ts = ExperimentalData.Torques;
ylim_Ts = [0,0;0,0;0,0;0,0;0,0;0,0;-50,50;-30,30;-20,20;-50,50;-30,30;...
    -20,20;-60,20;-60,20;-50,50;-50,50;-20,20;-20,20;0,0;0,0];
pos_Ts = 6:10;
idx_Ts = [10:12,14,16,18];

for mpi = 1:NPhases(ww(1)).m
count=1;
for i = 1:length(idx_Qs)
    subplot(5,5,pos_Ts(i))
    if find(idx_Ts==idx_Qs(i))
    % Experimental data of the CP child.
    idx_jref = strcmp(ref_Ts.(subject).colheaders,names_Qs{idx_Qs(i)});
    meanPlusSTD = ref_Ts.(subject).mean(:,idx_jref) + 2*ref_Ts.(subject).std(:,idx_jref);
    meanMinusSTD = ref_Ts.(subject).mean(:,idx_jref) - 2*ref_Ts.(subject).std(:,idx_jref);  
    stepID = (size(Ts_r(ww(1)).m(mpi).mp,1)-1)/(size(meanPlusSTD,1)-1);
    intervalID = 1:stepID:size(Ts_r(ww(1)).m(mpi).mp,1);
    sampleID = 1:size(Ts_r(ww(1)).m(mpi).mp,1);
    meanPlusSTD = interp1(intervalID,meanPlusSTD,sampleID);
    meanMinusSTD = interp1(intervalID,meanMinusSTD,sampleID); 
    hold on
    yy = fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'-','linewidth',2);
    yy.EdgeColor = [192,192,192]/255;
    yy.FaceColor = [192,192,192]/255; 
    end
    % Experimental data of the TD child.
    plot(Ts_toTrack(ww(1)).m(mpi).mp(:,1),Ts_toTrack(ww(1)).m(mpi).mp(:,idx_Qs(i)+1),...
        'color','k','linewidth',line_linewidth);
    for k = 1:length(ww)        
        % Simulation results
        plot(Ts_toTrack(ww(k)).m(mpi).mp(:,1),Ts_r(ww(k)).m(mpi).mp(:,idx_Qs(i)),...
            'color',color_all(k,:),'linewidth',line_linewidth);
        hold on;            
        [corrTD(ww(k)).m(mpi).mp.r2.T.all(i),corrTD(ww(k)).m(mpi).mp.rmse.T.all(i)] = ...
            rsquare(Ts_toTrack(ww(k)).m(mpi).mp(:,idx_Qs(i)+1),Ts_r(ww(k)).m(mpi).mp(:,idx_Qs(i)));  
        corrTDAll(mpi).mp.r2.T.all(count,i) = corrTD(ww(k)).m(mpi).mp.r2.T.all(i);
        corrTDAll(mpi).mp.rmse.T.all(count,i) = corrTD(ww(k)).m(mpi).mp.rmse.T.all(i);
        count = count+1;
        % Stance/swing
        temp_name = names_Qs{idx_Qs(i)};
        if strcmp(temp_name(end-1:end),'_l')
            plot([Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1)],...
                [ylim_Ts(idx_Qs(i),1) ylim_Ts(idx_Qs(i),2)],...
                'color',color_all(k,:),'linewidth',1);
            plot([Qs_toTrack(ww(k)).m(mpi).mp(HS(mpi).l(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(HS(mpi).l(k),1)],...
                [ylim_Ts(idx_Qs(i),1) ylim_Ts(idx_Qs(i),2)],...
                'color',color_all(k,:),'linestyle','--','linewidth',1);
            
        else
            plot([Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).r(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).r(k),1)],...
                [ylim_Ts(idx_Qs(i),1) ylim_Ts(idx_Qs(i),2)],...
                'color',color_all(k,:),'linewidth',1);
        end
    end
    count=1;         
    % Plot settings 
    set(gca,'Fontsize',label_fontsize);    
    % Y-axis
    ylim([ylim_Ts(idx_Qs(i),1) ylim_Ts(idx_Qs(i),2)]);
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks_Qs));       
    if i == 1
        ylabel('Torques (Nm)','Fontsize',label_fontsize);
    end      
    % X-axis
    xlim([Qs_toTrack_deg(1,1),Qs_toTrack_deg(end,1)])
    set(gca,'XTick',[]);
    box off;
end    
end

%% Plot muscle activations
ref_EMG = ExperimentalData.EMG;
muscleNames = {'Glut max 1','Glut max 2','Glut max 3','Glut med 1',...
    'Glut med 2','Glut med 3','Glut min 1','Glut min 2',...
    'Glut min 3','Add long','Add brev','Add mag 1','Add mag 2',...
    'Add mag 3','Pectineus','Iliacus','Psoas','Quad fem',...
    'Gemellus','Piri','TFL','Gracilis','Semimem','Semiten',...
    'Bi fem lh','Bi fem sh','Sartorius','Rectus fem',...
    'Vas med','Vas int','Vas lat','Gas med',...
    'Gas lat','Soleus','Tib post','Tib ant','Ext dig',...
    'Ext hal','Flex dig','Flex hal','Per brev','Per long',...
    'Per tert'};
muscleNames_tit = {'Glut max 1','Glut max 2','Glut max 3','Glut med 1',...
    'Gluteus medius 2','Glut med 3','Glut min 1','Glut min 2',...
    'Glut min 3','Add long','Add brev','Add mag 1','Add mag 2',...
    'Add mag 3','Pectineus','Iliacus','Psoas','Quad fem',...
    'Gemellus','Piri','TFL','Gracilis','Semimembranosus','Semiten',...
    'Biceps femoris lh','Bi fem sh','Sartorius','Rectus femoris',...
    'Vas med','Vas int','Vastus lateralis','Gastrocnemius medialis',...
    'Gastrocnemius lateralis','Soleus','Tib post','Tibialis anterior','Ext dig',...
    'Ext hal','Flex dig','Flex hal','Per brev','Per long',...
    'Per tert'};
EMGchannel = {'no','no','no','no',...
    'GLU','no','no','no',...
    'no','no','no','no','no',...
    'no','no','no','no','no',...
    'no','no','no','no','MEH','MEH',...
    'BIF','no','no','no',...
    'no','no','no','GAS',...
    'GAS','SOL','no','TIA','no',...
    'no','no','no','no','no',...
    'no'};
EMGcol = ones(1,length(muscleNames));
for i = 1:length(muscleNames)
    if strcmp(EMGchannel{i},'no')
        EMGcol(i) = 99;
    end
end
pos_Acts = 11:20;
idx_Acts = [5,17,23,25,28,31,32,33,34,36]; 
NMuscles = size(Acts_noSpas_r(ww(k)).m(mpi).mp,2);
% Experimental data of the CP child.
load([pathRepo,'/OpenSimModel/',subject,'/EMG/Gait/EMG_filt.mat'],'EMG_filt');
for mpi = 1:NPhases(ww(1)).m
count = 1;
for i = 1:length(idx_Acts)
    subplot(5,5,pos_Acts(i))    
    if EMGcol(idx_Acts(i)) ~= 99
        x = 1:(100-1)/(size(Acts_noSpas_r(ww(1)).m(mpi).mp,1)-1):100;
        % Normalize peak EMG to peak muscle activation
        a_peak_vec = zeros(1,length(ww));
        for k = 1:length(ww)
        a_peak_vec(k) = max(Acts_r(ww(k)).m(mpi).mp(:,idx_Acts(i)+NMuscles/2));
        end
        a_peak = max(a_peak_vec);    
        emg_peak = zeros(1,size(ref_EMG.(subject).all,3));
        for j = 1:size(ref_EMG.(subject).all,3)
            emg_peak(j) = nanmax(ref_EMG.(subject).all(:,strcmp(ref_EMG.(subject).colheaders,['R',EMGchannel{idx_Acts(i)}]),j),[],1);
        end
        norm_f = a_peak./emg_peak;
        tempp(:,:) = ref_EMG.(subject).all(:,strcmp(ref_EMG.(subject).colheaders,['R',EMGchannel{idx_Acts(i)}]),:);
        emg1processednorm.all = tempp.*repmat(norm_f,size(tempp,1),1);   
        emg1processednorm.mean = mean(emg1processednorm.all,2);
        emg1processednorm.std = std(emg1processednorm.all,[],2);         
        meanPlusSTD = emg1processednorm.mean + 2*emg1processednorm.std;
        meanMinusSTD = emg1processednorm.mean - 2*emg1processednorm.std; 
        intervalInterp = 1:(size(tempp,1)-1)/(size(Acts_noSpas_r(ww(1)).m(mpi).mp,1)-1):size(tempp,1);
        meanPlusSTD = interp1(1:size(tempp,1),meanPlusSTD,intervalInterp);
        meanMinusSTD = interp1(1:size(tempp,1),meanMinusSTD,intervalInterp); 
        meanInterp = interp1(1:size(tempp,1),emg1processednorm.mean,intervalInterp); 
        yy = fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'-','linewidth',2);
        yy.EdgeColor = [192,192,192]/255;
        yy.FaceColor = [192,192,192]/255;
        hold on
    end
    p = gobjects(1,length(ww));   
    % Simulation results
    for k = 1:length(ww)
        x = 1:(100-1)/(size(Acts_noSpas_r(ww(k)).m(mpi).mp,1)-1):100;
        p(k) = plot(x,Acts_noSpas_r(ww(k)).m(mpi).mp(:,idx_Acts(i)+NMuscles/2),'color',...
            color_all(k,:),'linewidth',line_linewidth);
        hold on
        count = count+1;
        plot([x(TO(mpi).r(k)) x(TO(mpi).r(k))],[0 1],'color',color_all(k,:),'linewidth',1);
    end   
    count = 1;
    % Plot settings
    set(gca,'Fontsize',label_fontsize)
    title(muscleNames_tit{idx_Acts(i)},'Fontsize',label_fontsize);    
    % X-axis
    set(gca,'XTick',[]);
    % Y-axis
    ylim([0,1]);
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
    if i == 1 || i == 6
        ylabel('Activations (-)','Fontsize',label_fontsize);
    end  
    box off;
end
end

%% Plot metabolic cost of transport (COT)
COT_order = zeros(1,length(ww));
for k = 1:length(ww)
    COT_order(k) = COT(ww(k)).m(mpi).mp;
end
COT_order_mean = COT_order;
COT_order_mean = [COT_order_mean;zeros(1,length(COT_order_mean))];
COT_order_std = zeros(2,length(COT_order_mean));
% Plot
ylim_COT = [0,4];
xlim_COT = [0.4 1.6];
NumTicks_COT = 2;
subplot(5,5,24)
h_COT = barwitherr(COT_order_std,COT_order_mean);
% Plot settings 
set(gca,'Fontsize',label_fontsize);
title('COT','Fontsize',label_fontsize);
% X-axis
xlim([xlim_COT(1,1) xlim_COT(1,2)]);
set(gca,'XTick',[]);
% Y-axis
ylim([ylim_COT(1,1) ylim_COT(1,2)]);
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks_COT));
ylabel('(J kg-1 m-1)','Fontsize',label_fontsize);
for k = 1:length(ww)
    set(h_COT(k),'FaceColor',color_all(k,:));
end
box off;

%% Plot stride length
SLengthref.mean = 1.241;
SLengthref.std = 0; 
% Post-processing for barplot
SLength_order = zeros(1,length(ww));
for k = 1:length(ww)
    SLength_order(k) = StrideLength(ww(k)).m(mpi).mp;
end
SLength_order_mean = SLength_order;
SLength_order_mean = [SLength_order_mean,SLengthref.mean,ExperimentalData.StrideLength.(subject).mean];
SLength_order_mean = [SLength_order_mean;zeros(1,length(SLength_order_mean))];
SLength_order_std = zeros(2,length(SLength_order_mean));
SLength_order_std(1,end) = 2*ExperimentalData.StrideLength.(subject).std;
% Plot
ylim_SLength = [0,1.5];
xlim_SLength = [0.4 1.6];
NumTicks_SLength = 2;
subplot(5,5,25)
h_SLength = barwitherr(SLength_order_std,SLength_order_mean);
% Plot settings 
set(gca,'Fontsize',label_fontsize);
title('Stride length','Fontsize',label_fontsize);
% X-axis
xlim([xlim_SLength(1,1) xlim_SLength(1,2)]);
set(gca,'XTick',[]);
% Y-axis
ylim([ylim_SLength(1,1) ylim_SLength(1,2)]);
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks_SLength));
ylabel('(m)','Fontsize',label_fontsize);
% Colors
for k = 1:length(ww)
    set(h_SLength(k),'FaceColor',color_all(k,:));
end
set(h_SLength(k+1),'FaceColor',[0,0,0]);
set(h_SLength(k+2),'FaceColor',[192,192,192]/255);
box off;
