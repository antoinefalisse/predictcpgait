% This script reproduces Fig 5
% Author: Antoine Falisse
% Date: 1/6/2020

clear all
close all
clc

%% User settings
% Selected trials
ww = [1,30,7,34];
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
GRFs_r          = struct('m',[]);
Acts_noSpas_r   = struct('m',[]);
Acts_Ff_r       = struct('m',[]);
Acts_dFf_r      = struct('m',[]);
Acts_r          = struct('m',[]);
NPhases         = struct('m',[]);
TO              = struct('r',[]);
HS              = struct('l',[]);
PredSimOCP_settings
pathresults = [pathPredictiveSimulations,'/Results/',ocp_path];
load([pathresults,'/ResultsPredSim.mat']);
threshold = 28;
NSyn = zeros(1,length(ww));
for k = 1:length(ww) 
    N = settings(ww(k),3);   % number of mesh intervals
    NPhases(ww(k)).m = 1;
    for mpi = 1:NPhases(ww(k)).m
        GRFs_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).GRFs_r(mpi).mp;
        Acts_noSpas_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_noSpas_r(mpi).mp;
        Acts_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_r(mpi).mp;
        Acts_Ff_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_Ff_r(mpi).mp;
        Acts_dFf_r(ww(k)).m(mpi).mp = ResultsPredSim.(['Case_',num2str(ww(k))]).Acts_dFf_r(mpi).mp;
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
    'Gemellus','Piri','TFL','Gracilis','Semimembranosus','Semitendinosus',...
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
% Load scaling factors  
side_lc = 'r';
muscleNames_Spas = {'Bi fem lh','Semimem','Semiten','Gas med','Gas lat'};
NMuscle_Spas_r  = length(muscleNames_Spas); % one side
musi_Spas_r = zeros(1,NMuscle_Spas_r);
for i = 1:NMuscle_Spas_r
    musi_Spas_r(i) = find(strcmp(muscleNames,muscleNames_Spas{i}));
end
NMuscles_spas = length(muscleNames_Spas);
pos_As = 1:5;
idx_As = [23,24,25,32,33]; 
NMuscles = size(Acts_noSpas_r(ww(k)).m(mpi).mp,2);

% Experimental data of the CP child.
load([pathRepo,'/OpenSimModel/',subject,'/EMG/Gait/EMG_filt.mat'],'EMG_filt');
for mpi = 1:NPhases(ww(1)).m
count = 1;
figure()
for i = 1:length(idx_As)
    subplot(5,5,pos_As(i))    
    if EMGcol(idx_As(i)) ~= 99
        x = 1:(100-1)/(size(Acts_noSpas_r(ww(1)).m(mpi).mp,1)-1):100;
        % Normalize EMG      
        % Normalize peak EMG to peak muscle activation
        a_peak_vec = zeros(1,2);
        for k = 1:2
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
    for k = 1:2
        x = 1:(100-1)/(size(Acts_noSpas_r(ww(k)).m(mpi).mp,1)-1):100;
        if settings(ww(k),11) ~= 1
            p(k) = plot(x,Acts_noSpas_r(ww(k)).m(mpi).mp(:,idx_As(i)+NMuscles/2),'color',...
                color_all(k,:),'linewidth',line_linewidth);
        end
        hold on
        if settings(ww(k),11) == 1
            p(k) = plot(x,Acts_noSpas_r(ww(k)).m(mpi).mp(:,idx_As(i)+NMuscles/2),'k','linestyle',':','linewidth',line_linewidth);          
            if find(musi_Spas_r==idx_As(i))
                idx_temp = find(musi_Spas_r==idx_As(i));
                plot(x,Acts_Ff_r(ww(k)).m(mpi).mp(:,idx_temp+NMuscles_spas) + ...
                    Acts_dFf_r(ww(k)).m(mpi).mp(:,idx_temp+NMuscles_spas),'k','linestyle','-','linewidth',2);
            end
            plot(x,Acts_r(ww(k)).m(mpi).mp(:,idx_As(i)+NMuscles/2),'color',...
            color_all(k,:),'linestyle','-','linewidth',2);
        end
        count = count+1;
        plot([x(TO(mpi).r(k)) x(TO(mpi).r(k))],[0 1],'color',color_all(k,:),'linewidth',1);
    end   
    count = 1;    
    % Plot settings
    set(gca,'Fontsize',label_fontsize)
    title(muscleNames_tit{idx_As(i)},'Fontsize',label_fontsize);    
    % X-axis
    L = get(gca,'XLim');
    set(gca,'XTick',[]);
    % Y-axis
    ylim([0,0.7]);
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
    if i == 1
        ylabel('Activations (-)','Fontsize',label_fontsize);
    end  
    box off;
end
end

for mpi = 1:NPhases(ww(1)).m  
pos_As = 6:10;
count = 1;
p = gobjects(1,4);
for i = 1:length(idx_As)
    subplot(5,5,pos_As(i))    
    if EMGcol(idx_As(i)) ~= 99
        x = 1:(100-1)/(size(Acts_noSpas_r(ww(1)).m(mpi).mp,1)-1):100;
        % Normalize EMG
        % Normalize peak EMG to peak muscle activation
        a_peak_vec = zeros(1,2);
        for k = 3:4
        a_peak_vec(k-2) = max(Acts_r(ww(k)).m(mpi).mp(:,idx_As(i)+NMuscles/2));
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
        p(5) = yy;
        hold on
    end       
    % Simulation results
    for k = 3:4
        x = 1:(100-1)/(size(Acts_noSpas_r(ww(k)).m(mpi).mp,1)-1):100;
        if settings(ww(k),11) ~= 1
            p(1) = plot(x,Acts_noSpas_r(ww(k)).m(mpi).mp(:,idx_As(i)+NMuscles/2),'color',...
                color_all(k-2,:),'linewidth',line_linewidth);
        end
        hold on
        if settings(ww(k),11) == 1
            p(3) = plot(x,Acts_noSpas_r(ww(k)).m(mpi).mp(:,idx_As(i)+NMuscles/2),'k','linestyle',':','linewidth',line_linewidth);          
            if find(musi_Spas_r==idx_As(i))
                idx_temp = find(musi_Spas_r==idx_As(i));
                p(4) = plot(x,Acts_Ff_r(ww(k)).m(mpi).mp(:,idx_temp+NMuscles_spas) + ...
                    Acts_dFf_r(ww(k)).m(mpi).mp(:,idx_temp+NMuscles_spas),'k','linestyle','-','linewidth',2);
            end
            p(2) = plot(x,Acts_r(ww(k)).m(mpi).mp(:,idx_As(i)+NMuscles/2),'color',...
            color_all(k-2,:),'linestyle','-','linewidth',2);
        end
        count = count+1;
        plot([x(TO(mpi).r(k)) x(TO(mpi).r(k))],[0 1],'color',color_all(k-2,:),'linewidth',1);
    end   
    count = 1;    
    % Plot settings
    set(gca,'Fontsize',label_fontsize)  
    % X-axis
    set(gca,'XTick',[]);
    % Y-axis
    ylim([0,0.7]);
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
    if i == 1
        ylabel('Activations (-)','Fontsize',label_fontsize);
    end  
    box off;
end
end
l = legend(p,{'No Spasticity','Spasticity: total activation','Spasticity: supra-spinal activation','Spasticity: spinal/spastic activation','Experimental CP child'});
set(l,'Fontsize',20)
