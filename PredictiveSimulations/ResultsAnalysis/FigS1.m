% This script reproduces Fig S1
% Author: Antoine Falisse
% Date: 1/6/2020

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
Qs_l            = struct('m',[]);
Acts_noSpas_l   = struct('m',[]);
Acts_l          = struct('m',[]);
Ts_l            = struct('m',[]);
GRFs_l          = struct('m',[]);
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
        Qs_l(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Qs_l(mpi).mp;
        Ts_l(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Ts_l(mpi).mp;
        GRFs_l(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).GRFs_l(mpi).mp;
        Qs_toTrack(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Qs_toTrack(mpi).mp;
        Ts_toTrack(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Ts_toTrack(mpi).mp;
        GRFs_toTrack(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).GRFs_toTrack(mpi).mp;
        Acts_noSpas_l(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_noSpas_l(mpi).mp;
        Acts_l(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_l(mpi).mp;
        COT(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).COT(mpi).mp;         
        StrideLength(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).StrideLength(mpi).mp; 
        % Extract stance/swing
        TO(mpi).r(k) = find(GRFs_l(ww(k)).m(mpi).mp(:,2)<threshold,1,'first');
        TO(mpi).l(k) = find(GRFs_l(ww(k)).m(mpi).mp(:,5)<threshold,1,'first');
        tempDiff = false(N,1);    
        tempDiff(2:end) = diff(GRFs_l(ww(k)).m(mpi).mp(:,2)>threshold) & ...
            diff(GRFs_l(ww(k)).m(mpi).mp(:,2))>0;
        HS(mpi).r(k) = find(tempDiff);       
    end
end
threshold = 30;
HS_TD = zeros(1,NPhases(ww(1)).m);
for mpi = 1:NPhases(ww(1)).m    
    HS_TD(mpi) = ...
        find(GRFs_toTrack(ww(1)).m(mpi).mp(:,6)<threshold,1,'last') + 1;     
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
load([pathExperimentalData,'/ExperimentalData_l.mat'],'ExperimentalData');

%% Plot joint kinematics
ref_Qs = ExperimentalData.Q;
pos_Qs = 1:5;
ylim_Qs = [-30,30;-30,30;-30,30;-2,2;0.5,1;-1,1;-10,50;-35,35;-25,25;-10,50;...
    -25,25;-25,25;0,90;0,90;-30,40;-30,40;-50,20;-50,20];
idx_Qs = [7,8,9,13,15]; 
refNames_Qs = {'lower-torso-RX','lower-torso-RY','lower-torso-RZ',...
    'lower-torso-TX','lower-torso-TY','lower-torso-TZ','Hip flexion',...
    'Hip adduction','Hip internal rotation','Hip flexion','Hip adduction',...
    'Hip internal rotation',...
    'Knee flexion','Knee flexion','Ankle dorsiflexion','Ankle dorsiflexion',...
    'subt-angle-l','subt-angle-r'};
names_Qs = {'lower_torso_RX','lower_torso_RY','lower_torso_RZ',...
    'lower_torso_TX','lower_torso_TY','lower_torso_TZ','hip_flex_l',...
    'hip_add_l','hip_rot_l','hip_flex_r','hip_add_r','hip_rot_r',...
    'knee_flex_l','knee_flex_r','ankle_flex_l','ankle_flex_r',...
    'subt_angle_l','subt_angle_r','lumbar_pitch','lumbar_roll','lumbar_yaw'};
NumTicks_Qs = 2;

figure()
for mpi = 1:NPhases(ww(1)).m
for i = 1:length(idx_Qs)
    p = gobjects(1,length(ww)+1);
    subplot(5,5,pos_Qs(i))      
    % Experimental data of the CP child.
    idx_jref = strcmp(ref_Qs.(subject).Qs.colheaders,names_Qs{idx_Qs(i)});
    meanPlusSTD = ref_Qs.(subject).Qs.mean(:,idx_jref) + 2*ref_Qs.(subject).Qs.std(:,idx_jref);
    meanMinusSTD = ref_Qs.(subject).Qs.mean(:,idx_jref) - 2*ref_Qs.(subject).Qs.std(:,idx_jref);          
    stepQ = (size(Qs_l(ww(1)).m(mpi).mp,1)-1)/(size(meanPlusSTD,1)-1);
    intervalQ = 1:stepQ:size(Qs_l(ww(1)).m(mpi).mp,1);
    sampleQ = 1:size(Qs_l(ww(1)).m(mpi).mp,1);
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
    % Part 1: from left heel strike to right heel strike
    p(length(ww)+1) = plot(Qs_toTrack_deg(1:N-HS_TD(mpi)+1,1),...
        Qs_toTrack_deg(HS_TD(mpi):end,idx_Qs(i)+1),...
        'color','k','linewidth',line_linewidth);
    % Part 2: from right heel strike to left heel strike
    plot(Qs_toTrack_deg(N-HS_TD(mpi)+2:N,1),...
        Qs_toTrack_deg(1:HS_TD(mpi)-1,idx_Qs(i)+1),...
        'color','k','linewidth',line_linewidth);
    % Indicate where it is discontinuous
    plot([Qs_toTrack_deg(N-HS_TD(mpi)+1,1),Qs_toTrack_deg(N-HS_TD(mpi)+1,1)],...
        [ylim_Qs(idx_Qs(i),1) ylim_Qs(idx_Qs(i),2)],'k--','linewidth',1);  
    for k = 1:length(ww)
        % Simulation results
        p(k) = plot(Qs_toTrack(ww(k)).m(mpi).mp(:,1),...
            Qs_l(ww(k)).m(mpi).mp(:,idx_Qs(i)),...
            'color',color_all(k,:),'linewidth',line_linewidth);    
        hold on;     
        % Stance/swing
        temp_name = names_Qs{idx_Qs(i)};
        if strcmp(temp_name(end-1:end),'_r')
            plot([Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).r(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).r(k),1)],...
                [ylim_Qs(idx_Qs(i),1) ylim_Qs(idx_Qs(i),2)],...
                'color',color_all(k,:),'linewidth',1);
            plot([Qs_toTrack(ww(k)).m(mpi).mp(HS(mpi).r(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(HS(mpi).r(k),1)],...
                [ylim_Qs(idx_Qs(i),1) ylim_Qs(idx_Qs(i),2)],...
                'color',color_all(k,:),'linestyle','--','linewidth',1);            
        else
            plot([Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1) ...
                Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1)],...
                [ylim_Qs(idx_Qs(i),1) ylim_Qs(idx_Qs(i),2)],...
                'color',color_all(k,:),'linewidth',1);
        end
    end 
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
l = legend(p,{'Generic parameters','Personalized parameters','Reference TD child','Experimental CP child'});
set(l,'Fontsize',20)

%% Plot ground reaction forces
ref_GRF = ExperimentalData.GRFs;
pos_GRF = 21:23;
ylim_GRF = [-50,50;0 200;-25 25;-50,50;0 200;-25 25];
NumTicks_GRF = 2;
names_GRF = {'Fore-aft','Vertical','Lateral','Fore-aft','Vertical','Lateral'};

for mpi = 1:NPhases(ww(1)).m
temp = GRFs_toTrack(ww(1)).m(mpi).mp(:,[3,6]);
temp(temp<0) = 0;
GRFs_toTrack(ww(1)).m(mpi).mp(:,[3,6])= temp;
% figure()
for i = 1:length(names_GRF)/2
    subplot(5,5,pos_GRF(i))
    % Experimental data of the CP child.
    meanPlusSTD = ref_GRF.(subject).mean(:,i+3) + 2*ref_GRF.(subject).std(:,i+3);    
    meanMinusSTD = ref_GRF.(subject).mean(:,i+3) - 2*ref_GRF.(subject).std(:,i+3);   
    stepGRF = (size(GRFs_l(ww(1)).m(mpi).mp,1)-1)/(size(meanPlusSTD,1)-1);
    intervalGRF = 1:stepGRF:size(GRFs_l(ww(1)).m(mpi).mp,1);
    sampleGRF = 1:size(GRFs_l(ww(1)).m(mpi).mp,1);
    meanPlusSTD = interp1(intervalGRF,meanPlusSTD,sampleGRF);
    meanMinusSTD = interp1(intervalGRF,meanMinusSTD,sampleGRF);
    hold on
    x = GRFs_toTrack(ww(1)).m(mpi).mp(:,1)';
    yy = fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'-','linewidth',2);
    yy.EdgeColor = [192,192,192]/255;
    yy.FaceColor = [192,192,192]/255; 
    % Experimental data of the CP child.
    % Part 1: from left heel strike to right heel strike
    plot(GRFs_toTrack(ww(k)).m(mpi).mp(1:N-HS_TD(mpi)+1,1),GRFs_toTrack(ww(1)).m(mpi).mp(HS_TD(mpi):end,i+1+3)./(body_weight/100),...
        'color','k','linewidth',line_linewidth);
    % Part 2: from right heel strike to left heel strike
    plot(GRFs_toTrack(ww(k)).m(mpi).mp(N-HS_TD(mpi)+2:N,1),GRFs_toTrack(ww(1)).m(mpi).mp(1:HS_TD(mpi)-1,i+1+3)./(body_weight/100),...
        'color','k','linewidth',line_linewidth);
    % Indicate where it is discontinuous
    plot([GRFs_toTrack(ww(k)).m(mpi).mp(N-HS_TD(mpi)+1,1),...
        GRFs_toTrack(ww(k)).m(mpi).mp(N-HS_TD(mpi)+1,1)],[ylim_GRF(i,1) ylim_GRF(i,2)],...
        'k--','linewidth',1);
    for k = 1:length(ww)
        % Simulation results
        plot(GRFs_toTrack(ww(k)).m(mpi).mp(:,1),GRFs_l(ww(k)).m(mpi).mp(:,i+3)./(body_weight/100),...
            'color',color_all(k,:),'linewidth',line_linewidth);
        hold on;   
        % Plot stance/swing
        plot([Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1) ...
            Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1)],...
            [ylim_GRF(i,1) ylim_GRF(i,2)],...
            'color',color_all(k,:),'linewidth',1);
    end
    % Plot settings 
    set(gca,'Fontsize',label_fontsize);    
    title(names_GRF{i},'Fontsize',label_fontsize);  
    % Y-axis
    ylim([ylim_GRF(i,1) ylim_GRF(i,2)]);
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks_GRF));       
    if i == 1
        ylabel('GRF (%BW)','Fontsize',label_fontsize);
    end      
    % X-axis
    xlim([Qs_toTrack_deg(1,1),Qs_toTrack_deg(end,1)])
    L = get(gca,'XLim');
    set(gca,'XTick',[]);
    box off;   
end
end

%% Plot joint torques
ref_Ts = ExperimentalData.Torques;
ylim_Ts = [0,0;0,0;0,0;0,0;0,0;0,0;-50,50;-35,35;
    -20,20;-50,50;-30,30;-20,20;-60,20;-60,20;-50,50;...
    -50,50;-20,20;-20,20;0,0;0,0];
pos_Ts = 6:10;
idx_Ts = [7:9,13,15,17];

for mpi = 1:NPhases(ww(1)).m
for i = 1:length(idx_Qs)
    subplot(5,5,pos_Ts(i))
    if find(idx_Ts==idx_Qs(i))
    % Experimental data of the CP child.
    idx_jref = strcmp(ref_Ts.(subject).colheaders,names_Qs{idx_Qs(i)});
    meanPlusSTD = ref_Ts.(subject).mean(:,idx_jref) + 2*ref_Ts.(subject).std(:,idx_jref);
    meanMinusSTD = ref_Ts.(subject).mean(:,idx_jref) - 2*ref_Ts.(subject).std(:,idx_jref);  
    stepID = (size(Ts_l(ww(1)).m(mpi).mp,1)-1)/(size(meanPlusSTD,1)-1);
    intervalID = 1:stepID:size(Ts_l(ww(1)).m(mpi).mp,1);
    sampleID = 1:size(Ts_l(ww(1)).m(mpi).mp,1);
    meanPlusSTD = interp1(intervalID,meanPlusSTD,sampleID);
    meanMinusSTD = interp1(intervalID,meanMinusSTD,sampleID); 
    hold on
    yy = fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'-','linewidth',2);
    yy.EdgeColor = [192,192,192]/255;
    yy.FaceColor = [192,192,192]/255; 
    end
    % Experimental data of the CP child.
    % Part 1: from left heel strike to right heel strike
    plot(Ts_toTrack(ww(1)).m(mpi).mp(1:N-HS_TD(mpi)+1,1),Ts_toTrack(ww(1)).m(mpi).mp(HS_TD(mpi):end,idx_Qs(i)+1),...
        'color','k','linewidth',line_linewidth);
    % Part 2: from right heel strike to left heel strike
    plot(Ts_toTrack(ww(1)).m(mpi).mp(N-HS_TD(mpi)+2:N,1),Ts_toTrack(ww(1)).m(mpi).mp(1:HS_TD(mpi)-1,idx_Qs(i)+1),...
        'color','k','linewidth',line_linewidth);
    % Indicate where it is discontinuous
    plot([Ts_toTrack(ww(1)).m(mpi).mp(N-HS_TD(mpi)+1,1),...
        Ts_toTrack(ww(1)).m(mpi).mp(N-HS_TD(mpi)+1,1)],[ylim_Ts(idx_Qs(i),1) ylim_Ts(idx_Qs(i),2)],...
        'k--','linewidth',1);  
    for k = 1:length(ww)        
        % Simulation results
        plot(Ts_toTrack(ww(k)).m(mpi).mp(:,1),Ts_l(ww(k)).m(mpi).mp(:,idx_Qs(i)),...
            'color',color_all(k,:),'linewidth',line_linewidth);
        hold on;       
        % Plot stance/swing
        temp_name = names_Qs{idx_Qs(i)};
        plot([Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1) ...
            Qs_toTrack(ww(k)).m(mpi).mp(TO(mpi).l(k),1)],...
            [ylim_Ts(idx_Qs(i),1) ylim_Ts(idx_Qs(i),2)],...
            'color',color_all(k,:),'linewidth',1);
    end   
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
    'Gluteus medius','Glut med 3','Glut min 1','Glut min 2',...
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
pos_As = 11:20;
idx_As = [5,17,23,25,28,31,32,33,34,36]; 
% Experimental data of the CP child.
load([pathRepo,'/OpenSimModel/',subject,'/EMG/Gait/EMG_filt.mat'],'EMG_filt');
for mpi = 1:NPhases(ww(1)).m 
count = 1;
for i = 1:length(idx_As)
    subplot(5,5,pos_As(i))    
    if EMGcol(idx_As(i)) ~= 99
        x = 1:(100-1)/(size(Acts_noSpas_l(ww(1)).m(mpi).mp,1)-1):100;
        % Normalize EMG
        % Normalize peak EMG to peak muscle activation
        a_peak_vec = zeros(1,length(ww));
        for k = 1:length(ww)
        a_peak_vec(k) = max(Acts_l(ww(k)).m(mpi).mp(:,idx_As(i)));
        end
        a_peak = max(a_peak_vec);    
        emg_peak = zeros(1,size(ref_EMG.(subject).all,3));
        for j = 1:size(ref_EMG.(subject).all,3)
            emg_peak(j) = nanmax(ref_EMG.(subject).all(:,strcmp(ref_EMG.(subject).colheaders,['L',EMGchannel{idx_As(i)}]),j),[],1);
        end
        norm_f = a_peak./emg_peak;
        tempp(:,:) = ref_EMG.(subject).all(:,strcmp(ref_EMG.(subject).colheaders,['L',EMGchannel{idx_As(i)}]),:);
        emg1processednorm.all = tempp.*repmat(norm_f,size(tempp,1),1);   
        emg1processednorm.mean = mean(emg1processednorm.all,2);
        emg1processednorm.std = std(emg1processednorm.all,[],2); 
        meanPlusSTD = emg1processednorm.mean + 2*emg1processednorm.std;
        meanMinusSTD = emg1processednorm.mean - 2*emg1processednorm.std; 
        intervalInterp = 1:(size(tempp,1)-1)/(size(Acts_noSpas_l(ww(1)).m(mpi).mp,1)-1):size(tempp,1);
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
        x = 1:(100-1)/(size(Acts_noSpas_l(ww(k)).m(mpi).mp,1)-1):100;
        p(k) = plot(x,Acts_noSpas_l(ww(k)).m(mpi).mp(:,idx_As(i)),'color',...
            color_all(k,:),'linewidth',line_linewidth);
        hold on
        count = count+1;
        plot([x(TO(mpi).l(k)) x(TO(mpi).l(k))],[0 1],'color',color_all(k,:),'linewidth',1);
    end   
    count = 1;
    % Plot settings
    set(gca,'Fontsize',label_fontsize)
    title(muscleNames_tit{idx_As(i)},'Fontsize',label_fontsize);    
    % X-axis
    L = get(gca,'XLim');
    set(gca,'XTick',[]);
    % Y-axis
    ylim([0,1]);
    NumTicks = 2;
    LY = get(gca,'YLim');    
    if i == 1 || i == 6
        ylabel('Activations (-)','Fontsize',label_fontsize);
        set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
    else
        set(gca,'YTick',[]);
    end  
    box off;
end
end
