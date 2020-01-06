% This script reproduces Fig 6
% Author: Antoine Falisse
% Date: 1/6/2020

clear all
close all
clc

%% User settings
% Selected trials
ww = [1,5];
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

%% Load results
% Pre-allocation structures
Qs_r            = struct('m',[]);
Acts_noSpas_r   = struct('m',[]);
Acts_r          = struct('m',[]);
GRFs_r          = struct('m',[]);
lMtilde_r       = struct('m',[]);
Fpe_r           = struct('m',[]);
Qs_toTrack      = struct('m',[]);
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
NSyn = zeros(1,length(ww));
for k = 1:length(ww) 
    N = settings(ww(k),3);   % number of mesh intervals
    NSyn(k) = settings(ww(k),5);   % number of synergies per leg
    NPhases(ww(k)).m = 1;
    for mpi = 1:NPhases(ww(k)).m
        Qs_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Qs_r(mpi).mp;
        GRFs_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).GRFs_r(mpi).mp;
        Qs_toTrack(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Qs_toTrack(mpi).mp;
        Acts_noSpas_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_noSpas_r(mpi).mp;
        Acts_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_r(mpi).mp;
        Fpe_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Fpe_r(mpi).mp;
        lMtilde_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).lMtilde_r(mpi).mp;
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

%% Experimental data
pathExperimentalData = [pathPredictiveSimulations,'/ExperimentalData'];
load([pathExperimentalData,'/ExperimentalData_r.mat'],'ExperimentalData');

%% Plot joint kinematics
ref_Qs = ExperimentalData.Q;
pos_Qs = 1;
ylim_Qs = [-30,30;-30,30;-30,30;-2,2;0.5,1;-1,1;-10,50;-25,25;-25,25;-10,50;...
    -25,25;-25,25;0,90;0,90;-30,40;-30,40;-50,20;-50,20];
idx_Qs = 14; 
refNames_Qs = {'lower-torso-RX','lower-torso-RY','lower-torso-RZ',...
    'lower-torso-TX','lower-torso-TY','lower-torso-TZ','Hip flexion',...
    'Hip adduction','Hip rotation','Hip flexion','Hip adduction','Hip rotation',...
    'Knee flexion','Knee flexion','Ankle flexion','Ankle flexion',...
    'subt-angle-l','subt-angle-r'};
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
    p(length(ww)+1) = plot(Qs_toTrack_deg(:,1),Qs_toTrack_deg(:,idx_Qs(i)+1),...
        'color','k','linewidth',line_linewidth);        
    for k = 1:length(ww)
        % Simulation results
        p(k) = plot(Qs_toTrack(ww(k)).m(mpi).mp(:,1),Qs_r(ww(k)).m(mpi).mp(:,idx_Qs(i)),...
            'color',color_all(k,:),'linewidth',line_linewidth);    
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
% Analyse RMS
corrCPAll.mp.rmse.angles.diff = diff(corrCPAll.mp.rmse.angles.all,1,1);
corrCPAll.mp.rmse.angles.diff_per = ...
    corrCPAll.mp.rmse.angles.diff./corrCPAll.mp.rmse.angles.all(2:end,:)*100;

%% Plot muscle activations
EMGref = ExperimentalData.EMG;
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
muscleNames_tit = {'Glut max 1','Gluteus maximus','Glut max 3','Glut med 1',...
    'Gluteus medius','Glut med 3','Glut min 1','Glut min 2',...
    'Glut min 3','Add long','Add brev','Add mag 1','Add mag 2',...
    'Add mag 3','Pectineus','Iliacus','Psoas','Quad fem',...
    'Gemellus','Piri','TFL','Gracilis','Semimembranosus','Semiten',...
    'Biceps femoris lh','Biceps femoris sh','Sartorius','Rectus femoris',...
    'Vastus medialis','Vas int','Vastus lateralis','Gastrocnemius medialis',...
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
pos_As = [2,3];
idx_As = [2,5]; 
NMuscles = size(Acts_noSpas_r(ww(k)).m(mpi).mp,2);

% Experimental data of the CP child.
load([pathRepo,'/OpenSimModel/',subject,'/EMG/Gait/EMG_filt.mat'],'EMG_filt');  
for mpi = 1:NPhases(ww(1)).m
count = 1;
for i = 1:length(idx_As)
    subplot(5,5,pos_As(i)) 
    p = gobjects(1,length(ww));   
    % Simulation results
    for k = 1:length(ww)
        x = 1:(100-1)/(size(Acts_noSpas_r(ww(k)).m(mpi).mp,1)-1):100;
        p(k) = plot(x,Acts_noSpas_r(ww(k)).m(mpi).mp(:,idx_As(i)+NMuscles/2).^10,'color',...
            color_all(k,:),'linewidth',line_linewidth);
        hold on
        count = count+1;
        plot([x(TO(mpi).r(k)) x(TO(mpi).r(k))],[0 4e-4],'color',color_all(k,:),'linewidth',1);
    end   
    count = 1;
    % Plot settings
    set(gca,'Fontsize',label_fontsize)
    title(muscleNames_tit{idx_As(i)},'Fontsize',label_fontsize);    
    % X-axis
    L = get(gca,'XLim');
        set(gca,'XTick',[]);
    % Y-axis
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
    if i == 1
        ylabel('Activations10 (-)','Fontsize',label_fontsize);
    end  
    box off;
end
end

%% Plot passive forces
Fpe_opt_adj = struct('m',[]);
lMtilde_opt_adj = struct('m',[]);
NMuscle_noFLV_r = struct('m',[]);
musi_noFLV_r  = struct('m',[]);
musi_noFLV = struct('m',[]);
musi_FLV = struct('m',[]);
muscleNames_noFLV(ww(k)).m = {};
pathOpenSimModel = [pathRepo,'/OpenSimModel/',subject,'/'];
pathParameterEstimation = [pathOpenSimModel,'ParameterEstimation/'];   
load([pathParameterEstimation,'MTParameters_personalized.mat']);  
for k = 1:length(ww)
    NMuscle_noFLV_r(ww(k)).m  = length(muscleNames_noFLV(ww(k)).m);
    musi_noFLV_r(ww(k)).m = zeros(1,NMuscle_noFLV_r(ww(k)).m);
    for i = 1:NMuscle_noFLV_r(ww(k)).m
        musi_noFLV_r(ww(k)).m(i) = find(strcmp(muscleNames,muscleNames_noFLV(ww(k)).m{i}));
    end
    musi_noFLV(ww(k)).m = [musi_noFLV_r(ww(k)).m,musi_noFLV_r(ww(k)).m+NMuscles/2]; 
    musi = 1:NMuscles;
    musi_FLV(ww(k)).m = musi;
    musi_FLV(ww(k)).m(musi_noFLV(ww(k)).m) = [];    
    for mpi = 1:NPhases(ww(k)).m    
        Fpe_opt_adj(ww(k)).m(mpi).mp = zeros(size(Fpe_r(ww(k)).m(mpi).mp,1),NMuscles);
        Fpe_opt_adj(ww(k)).m(mpi).mp(:,musi_FLV(ww(k)).m) = Fpe_r(ww(k)).m(mpi).mp;        
        lMtilde_opt_adj(ww(k)).m(mpi).mp = zeros(size(lMtilde_r(ww(k)).m(mpi).mp,1),NMuscles);
        lMtilde_opt_adj(ww(k)).m(mpi).mp(:,musi_FLV(ww(k)).m) = lMtilde_r(ww(k)).m(mpi).mp;
    end    
end
for mpi = 1:NPhases(ww(k)).m
pos_pF = [4,5];
idx_pF = [17,26]; 
for i = 1:length(idx_pF)
    subplot(5,5,pos_pF(i))    
    for k = 1:length(ww)
        x = 1:(100-1)/(size(Acts_noSpas_r(ww(k)).m(mpi).mp,1)-1):100;            
        plot(x,Fpe_opt_adj(ww(k)).m(mpi).mp(:,idx_pF(i)+NMuscles/2)./MTParameters(1,idx_pF(i)+NMuscles/2),...
            'color',color_all(k,:),'linewidth',line_linewidth);
        hold on;
        plot([x(TO(mpi).r(k)) x(TO(mpi).r(k))],[0 1],'color',color_all(k,:),'linewidth',1);
    end
    set(gca,'Fontsize',label_fontsize)
    title(muscleNames_tit{idx_pF(i)},'Fontsize',label_fontsize);
    set(gca,'XTick',[]);
    ylim([0,0.5]);
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
    if i == 1
        ylabel('Passive forces (-)','Fontsize',label_fontsize);
    end  
    box off;
end
end
