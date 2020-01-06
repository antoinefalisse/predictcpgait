% This script reproduces Fig 4
% Author: Antoine Falisse
% Date: 1/6/2020

clear all
close all
clc

%% User settings
% Selected trials
ww = [1,7,22];
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
thesholdSynW = 0.5;

%% Load results
% Pre-allocation structures
Qs_r            = struct('m',[]);
Acts_noSpas_r   = struct('m',[]);
Acts_syn_r      = struct('m',[]);
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
Wsyn            = struct('m',[]);
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
        Ts_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Ts_r(mpi).mp;
        GRFs_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).GRFs_r(mpi).mp;
        Qs_toTrack(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Qs_toTrack(mpi).mp;
        Ts_toTrack(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Ts_toTrack(mpi).mp;
        GRFs_toTrack(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).GRFs_toTrack(mpi).mp;
        Acts_noSpas_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_noSpas_r(mpi).mp;
        Acts_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_r(mpi).mp;
        COT(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).COT(mpi).mp;         
        StrideLength(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).StrideLength(mpi).mp; 
        if NSyn(k) ~= 99
            Acts_syn_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_syn_r(mpi).mp;
            Wsyn(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).w_syn(mpi).mp;
        end
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

%% Plot joint angles
ref_Qs = ExperimentalData.Q;
pos_Qs = 1:3;
ylim_Qs = [-30,30;-30,30;-30,30;-2,2;0.5,1;-1,1;-10,50;-25,25;-25,25;-10,50;...
    -25,25;-25,25;0,90;0,90;-30,40;-30,40;-50,20;-50,20];
idx_Qs = [10,14,16]; 
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
l = legend(p,{'No synergies','4 synergies','3 synergies','Reference TD child','Experimental CP child'});
set(l,'Fontsize',20)

%% Plot joint torques
ref_Ts = ExperimentalData.Torques;
ylim_Ts = [0,0;0,0;0,0;0,0;0,0;0,0;-50,50;-30,30;-20,20;-50,50;-30,30;-20,20;...
    -60,20;-60,20;-50,50;-50,50;-20,20;-20,20;0,0;0,0];
pos_Ts = 6:8;
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
pos_As = [11,16];
idx_As = [32,33]; 
NMuscles = size(Acts_noSpas_r(ww(k)).m(mpi).mp,2);

% Experimental data of the CP child.
load([pathRepo,'/OpenSimModel/',subject,'/EMG/Gait/EMG_filt.mat'],'EMG_filt');
for mpi = 1:NPhases(ww(1)).m
count = 1;
for i = 1:length(idx_As)
    subplot(5,5,pos_As(i))    
    if EMGcol(idx_As(i)) ~= 99
        x = 1:(100-1)/(size(Acts_noSpas_r(ww(1)).m(mpi).mp,1)-1):100;
        % Normalize EMG
        % Normalize peak EMG to peak muscle activation
        a_peak_vec = zeros(1,length(ww));
        for k = 1:length(ww)
        a_peak_vec(k) = max(Acts_r(ww(k)).m(mpi).mp(:,idx_As(i)+NMuscles/2));
        end
        a_peak = max(a_peak_vec);                
        emg_peak = zeros(1,size(EMGref.(subject).all,3));
        for j = 1:size(EMGref.(subject).all,3)
            emg_peak(j) = nanmax(EMGref.(subject).all(:,strcmp(EMGref.(subject).colheaders,['R',EMGchannel{idx_As(i)}]),j),[],1);
        end
        norm_f = a_peak./emg_peak;
        tempp(:,:) = EMGref.(subject).all(:,strcmp(EMGref.(subject).colheaders,['R',EMGchannel{idx_As(i)}]),:);
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
        p(k) = plot(x,Acts_noSpas_r(ww(k)).m(mpi).mp(:,idx_As(i)+NMuscles/2),'color',...
            color_all(k,:),'linewidth',line_linewidth);
        hold on
        count = count+1;
        plot([x(TO(mpi).r(k)) x(TO(mpi).r(k))],[0 1],'color',color_all(k,:),'linewidth',1);
    end   
    count = 1;
    % Plot settings
    set(gca,'Fontsize',label_fontsize)
    title(muscleNames_tit{idx_As(i)},'Fontsize',label_fontsize);    
    % X-axis
    set(gca,'XTick',[]);
    % Y-axis
    ylim([0,1]);
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
    ylabel('Activations (-)','Fontsize',label_fontsize);
    box off;
end
end

%% Plot synergy activation profiles: 4 synergies
pos_Syn4.right = [10,25,20,15];
[~,pos_Syn4i.right] = sort(pos_Syn4.right);
pos_Syn4.left = [15,25,10,20];
[~,pos_Syn4i.left] = sort(pos_Syn4.left);
for mpi = 1:NPhases(ww(1)).m
for i = 1:4
    subplot(5,5,pos_Syn4.right(i))
    p = gobjects(1,length(ww));
    for k = 2
        if settings(ww(k),13) ~= 0
            if size(Acts_syn_r(ww(k)).m(mpi).mp,2)/2 >= i
            x = 1:(100-1)/(size(Acts_noSpas_r(ww(k)).m(mpi).mp,1)-1):100;
            p(k) = plot(x,Acts_syn_r(ww(k)).m(mpi).mp(:,i+NSyn(k)),'color',color_all(k,:),...
                'linewidth',line_linewidth);
            A_syn.four.right(:,i) = Acts_syn_r(ww(k)).m(mpi).mp(:,i+NSyn(k));
            A_syn.four.left(:,i) = Acts_syn_r(ww(k)).m(mpi).mp(:,i);
            hold on;
            plot([x(TO(mpi).r(k)) x(TO(mpi).r(k))],[0 1],'color',color_all(k,:),'linewidth',1);
            end
        end
    end
     % Plot settings
    set(gca,'Fontsize',label_fontsize)
    title(['4 synergies: ',num2str(pos_Syn4i.right(i))],'Fontsize',label_fontsize);    
    % X-axis
    set(gca,'XTick',[]);
    % Y-axis
    ylim([0,1]);
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
    ylabel('Activations (-)','Fontsize',label_fontsize);   
    box off;
end
end

%% Plot synergy weights: 4 synergies
pos_WSyn4.right = pos_Syn4.right-1;
[~,pos_WSyn4i.right] = sort(pos_WSyn4.right);
pos_WSyn4.left = pos_Syn4.left-1;
[~,pos_WSyn4i.left] = sort(pos_WSyn4.left);
for mpi = 1:NPhases(ww(1)).m
for i = 1:4
    subplot(5,5,pos_WSyn4.right(i))
    p = gobjects(1,length(ww));
    for k = 2
        if settings(ww(k),13) ~= 0
            if size(Acts_syn_r(ww(k)).m(mpi).mp,2)/2 >= i
                w_syn_resh.four.right = reshape(Wsyn(ww(k)).m(mpi).mp(:,2),NMuscles/2,NSyn(k)); 
                w_syn_resh.four.right = w_syn_resh.four.right./max(w_syn_resh.four.right);
                w_syn_resh.four.left = reshape(Wsyn(ww(k)).m(mpi).mp(:,1),NMuscles/2,NSyn(k)); 
                w_syn_resh.four.left = w_syn_resh.four.left./max(w_syn_resh.four.left);
                b = bar(w_syn_resh.four.right(:,i));
                b.FaceColor = color_all(k,:);
                b.EdgeColor = color_all(k,:);
                hold on;
                for j = 5:5:40
                    plot([j j],[0 1],'k:','linewidth',1);
                end            
            end
        end
    end
    % Plot settings
    set(gca,'Fontsize',label_fontsize)
    title(['4 synergies: ',num2str(pos_Syn4i.right(i))],'Fontsize',label_fontsize);   
    % X-axis
    set(gca,'XTick',[]);
    % Y-axis
    ylim([0,1]);
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
    ylabel('Weights (-)','Fontsize',label_fontsize);
    box off;
end
end
% Identify muscles in synegies
highWSynergies.four.four.names = muscleNames(w_syn_resh.four.right(:,pos_WSyn4i.right(1)) > thesholdSynW);
highWSynergies.four.four.values = w_syn_resh.four.right((w_syn_resh.four.right(:,pos_WSyn4i.right(1)) > thesholdSynW),pos_WSyn4i.right(1));
highWSynergies.four.one.names = muscleNames(w_syn_resh.four.right(:,pos_WSyn4i.right(2)) > thesholdSynW);
highWSynergies.four.one.values = w_syn_resh.four.right((w_syn_resh.four.right(:,pos_WSyn4i.right(2)) > thesholdSynW),pos_WSyn4i.right(2));
highWSynergies.four.two.names = muscleNames(w_syn_resh.four.right(:,pos_WSyn4i.right(3)) > thesholdSynW);
highWSynergies.four.two.values = w_syn_resh.four.right((w_syn_resh.four.right(:,pos_WSyn4i.right(3)) > thesholdSynW),pos_WSyn4i.right(3));
highWSynergies.four.three.names = muscleNames(w_syn_resh.four.right(:,pos_WSyn4i.right(4)) > thesholdSynW);
highWSynergies.four.three.values = w_syn_resh.four.right((w_syn_resh.four.right(:,pos_WSyn4i.right(4)) > thesholdSynW),pos_WSyn4i.right(4));

%% Plot synergy activation profiles: 3 synergies
% adjust order to match 4 synergies
pos_Syn3.right = [23,18,13];
[~,pos_Syn3i.right] = sort(pos_Syn3.right);
pos_Syn3.left = [23,18,13];
[~,pos_Syn3i.left] = sort(pos_Syn3.left);
for mpi = 1:NPhases(ww(1)).m
for i = 1:3
    subplot(5,5,pos_Syn3.right(i))
    p = gobjects(1,length(ww));
    for k = 3
        if settings(ww(k),13) ~= 0
            if size(Acts_syn_r(ww(k)).m(mpi).mp,2)/2 >= i
            x = 1:(100-1)/(size(Acts_noSpas_r(ww(k)).m(mpi).mp,1)-1):100;
            p(k) = plot(x,Acts_syn_r(ww(k)).m(mpi).mp(:,i+NSyn(k)),'color',color_all(k,:),...
                'linewidth',line_linewidth);
            A_syn.three.right(:,i) = Acts_syn_r(ww(k)).m(mpi).mp(:,i+NSyn(k));
            A_syn.three.left(:,i) = Acts_syn_r(ww(k)).m(mpi).mp(:,i);
            hold on;
            plot([x(TO(mpi).r(k)) x(TO(mpi).r(k))],[0 1],'color',color_all(k,:),'linewidth',1);
            end
        end
    end
     % Plot settings
    set(gca,'Fontsize',label_fontsize)
    title(['3 synergies: ',num2str(pos_Syn3i.right(i))],'Fontsize',label_fontsize);    
    % X-axis
    set(gca,'XTick',[]);
    % Y-axis
    ylim([0,1]);
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
    ylabel('Activations (-)','Fontsize',label_fontsize);
box off;
end
end

%% Plot synergy weights: 3 synergies
% adjust order to match 4 synergies
pos_WSyn3.right = pos_Syn3.right-1;
[~,pos_WSyn3i.right] = sort(pos_WSyn3.right);
pos_WSyn3.left = pos_Syn3.left-1;
[~,pos_WSyn3i.left] = sort(pos_WSyn3.left);
for mpi = 1:NPhases(ww(1)).m
for i = 1:3
    subplot(5,5,pos_WSyn3.right(i))
    p = gobjects(1,length(ww));
    for k = 3
        if settings(ww(k),13) ~= 0
            if size(Acts_syn_r(ww(k)).m(mpi).mp,2)/2 >= i
                w_syn_resh.three.right = reshape(Wsyn(ww(k)).m(mpi).mp(:,2),NMuscles/2,NSyn(k)); 
                w_syn_resh.three.right = w_syn_resh.three.right./max(w_syn_resh.three.right);
                w_syn_resh.three.left = reshape(Wsyn(ww(k)).m(mpi).mp(:,1),NMuscles/2,NSyn(k)); 
                w_syn_resh.three.left = w_syn_resh.three.left./max(w_syn_resh.three.left);
                b = bar(w_syn_resh.three.right(:,i));
                b.FaceColor = color_all(k,:);
                b.EdgeColor = color_all(k,:);
                hold on;
                for j = 5:5:40
                    plot([j j],[0 1],'k:','linewidth',1);
                end            
            end
        end
    end
     % Plot settings
    set(gca,'Fontsize',label_fontsize)
    title(['3 synergies: ',num2str(pos_WSyn3i.right(i))],'Fontsize',label_fontsize);    
    % X-axis
    set(gca,'XTick',[]);
    % Y-axis
    ylim([0,1]);
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
    ylabel('Weights (-)','Fontsize',label_fontsize);
    box off;
end
end
% Identify muscles in synegies
highWSynergies.three.one.names = muscleNames(w_syn_resh.three.right(:,pos_WSyn3i.right(1)) > thesholdSynW);
highWSynergies.three.one.values = w_syn_resh.three.right((w_syn_resh.three.right(:,pos_WSyn3i.right(1)) > thesholdSynW),pos_WSyn3i.right(1));
highWSynergies.three.two.names = muscleNames(w_syn_resh.three.right(:,pos_WSyn3i.right(2)) > thesholdSynW);
highWSynergies.three.two.values = w_syn_resh.three.right((w_syn_resh.three.right(:,pos_WSyn3i.right(2)) > thesholdSynW),pos_WSyn3i.right(2));
highWSynergies.three.three.names = muscleNames(w_syn_resh.three.right(:,pos_WSyn3i.right(3)) > thesholdSynW);
highWSynergies.three.three.values = w_syn_resh.three.right((w_syn_resh.three.right(:,pos_WSyn3i.right(3)) > thesholdSynW),pos_WSyn3i.right(3));

%% Plot metabolic cost of transport (COT)
COT_order = zeros(1,length(ww));
for k = 1:length(ww)
    COT_order(k) = COT(ww(k)).m(mpi).mp;
end
COT_order_mean = COT_order;
COT_order_mean = [COT_order_mean;zeros(1,length(COT_order_mean))];
COT_order_std = zeros(2,length(COT_order_mean));
% Plot
ylim_COT = [0,5];
xlim_COT = [0.4 1.6];
NumTicks_COT = 2;
subplot(5,5,4)
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

%% Stride lengths
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
subplot(5,5,5)
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

%% Expand synergy plot
figure()
subplot(5,5,[2:5])
k=2;
b = bar(w_syn_resh.four.right(:,2),'stacked');
hold on
for j = 5:5:40
    plot([j j],[0 1],'k:','linewidth',1);
end
b.FaceColor = color_all(k,:);
b.EdgeColor = color_all(k,:);
size_x = length(w_syn_resh.four.right(:,1))+1;
Xt = 0:size_x;
Xl = [0 size_x];
set(gca,'XTick',Xt,'XLim',Xl);
temp = cell(1,size_x+1);
temp(2:end-1) = muscleNames;
set(gca,'XTickLabel',temp);
xtickangle(90)
ylim([0,1]);
NumTicks = 2;
LY = get(gca,'YLim');
set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
ylabel('Weights (-)','Fontsize',label_fontsize);
box off;

%% Compare synergy composition
nameFields = fieldnames(highWSynergies.three);
for i=1:length(nameFields)
    count = 0;
    for j = 1:length(highWSynergies.three.(nameFields{i}).names)
        if find(strcmp(highWSynergies.three.(nameFields{i}).names{j},highWSynergies.four.(nameFields{i}).names))
            count = count + 1;
            highWSynergies.comp.(nameFields{i}).names{count} = highWSynergies.three.(nameFields{i}).names{j};
        end
    end
end

%% Compare synergy composition
% Pearson's coefficents
w_syn_resh.four.all = [w_syn_resh.four.left,w_syn_resh.four.right];
w_syn_resh.three.all = [w_syn_resh.three.left,w_syn_resh.three.right];
for i = 1:3
    [corrSyn.right.W.r(i),corrSyn.right.W.p(i)] = corr(w_syn_resh.four.right(:,pos_WSyn4i.right(i+1)),...
        w_syn_resh.three.right(:,pos_WSyn3i.right(i)));       
    [corrSyn.right.A.r2(i),corrSyn.right.A.rmse(i)] = ...
        rsquare(A_syn.four.right(:,pos_Syn4i.right(i+1)),A_syn.three.right(:,pos_Syn3i.right(i)));    
    
    [corrSyn.left.W.r(i),corrSyn.left.W.p(i)] = corr(w_syn_resh.four.left(:,pos_WSyn4i.left(i+1)),...
        w_syn_resh.three.left(:,pos_WSyn3i.left(i)));       
    [corrSyn.left.A.r2(i),corrSyn.left.A.rmse(i)] = ...
        rsquare(A_syn.four.left(:,pos_Syn4i.left(i+1)),A_syn.three.left(:,pos_Syn3i.left(i)));  
end

corrSyn.both.W.r.all = [corrSyn.left.W.r,corrSyn.right.W.r];
corrSyn.both.W.r.mean = round(mean(corrSyn.both.W.r.all),2);
corrSyn.both.W.r.std = round(std(corrSyn.both.W.r.all),2);
corrSyn.both.A.r2.all = [corrSyn.left.A.r2,corrSyn.right.A.r2];
corrSyn.both.A.r2.mean = round(mean(corrSyn.both.A.r2.all),2);
corrSyn.both.A.r2.std = round(std(corrSyn.both.A.r2.all),2);
sumResults.w_syn_resh = w_syn_resh;
sumResults.A_syn = A_syn;
sumResults.pos_WSyn4i = pos_WSyn4i;
sumResults.pos_WSyn3i = pos_WSyn3i;
% save('sumResults_personalized','sumResults')
