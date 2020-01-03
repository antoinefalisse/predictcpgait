% This function generates polynomials to approximate muscle-tendon lengths
% and moment arms. The code is from Wouter Aerts and was adapted for use 
% with CasADi.
%
% Author: Antoine Falisse
% Date: 03/04/2018
%
clear all
close all
clc

%% User inputs
runPolynomialfit = 1;
saveQdot = 0;
savePolynomials = 1;

%% Extract time and angles from dummy motion

pathmain = pwd;
[pathPredictiveSimulations,~,~] = fileparts(pathmain);
[pathRepo,~,~] = fileparts(pathPredictiveSimulations);
pathOpenSimModel = [pathRepo,'/OpenSimModel/'];

name_dummymotion = 'dummy_motion.mot';
subject = 'subject1';

pathPolynomials = [pathOpenSimModel,subject,'/Polynomials/'];
pathMuscleAnalysis = [pathPolynomials,'/MuscleAnalysis/'];

dummy_motion = importdata([pathPolynomials,name_dummymotion]);
% Order of dofs: 6*pelvis, hip flex r, hip add r, hip rot r, knee flex r, 
% ankle flex r, hip flex l, hip add l, hip rot l, knee flex l,
% ankle flex l, lumbar ext, lumbar bend, lumbar rot, subtalar r, subtalar l
% Although we should get different polynomials for both legs since the
% model is asymmetric (MRI-based model), we use the same Qs as inputs
q = dummy_motion.data(:,[8:12,21]).*(pi/180);

% Generate random numbers between -1000 and 1000 (°/s) 
if saveQdot
    a = -1000;
    b = 1000;
    r1 = (b-a).*rand(size(q,1),1) + a;
    r2 = (b-a).*rand(size(q,1),1) + a;
    r3 = (b-a).*rand(size(q,1),1) + a;
    r4 = (b-a).*rand(size(q,1),1) + a;
    r5 = (b-a).*rand(size(q,1),1) + a;
    r6 = (b-a).*rand(size(q,1),1) + a;
    r = [r1,r2,r3,r4,r5,r6];
    qdot = r.*(pi/180);
    dummy_qdot = qdot;
    save([pathPolynomials,'dummy_qdot.mat'],'dummy_qdot');
end
load([pathPolynomials,'dummy_qdot.mat']);
qdot = dummy_qdot(:,:);

% The nodel is asymmetric so we run twice the polynomial approximation,
% once for each leg 
sides = {'r','l'};
Qs_idx.r = [8:12,21];
Qs_idx.l = [13:17,22];
for s = 1:length(sides)
    side = sides{s};
    %% Import data
    % lMT
    lMT = importdata([pathMuscleAnalysis,'subject01_MuscleAnalysis_Length.sto']);
    % hip flexion
    MA.hip.flex = importdata([pathMuscleAnalysis,'subject01_MuscleAnalysis_MomentArm_hip_flex_',side,'.sto']);
    % hip adduction
    MA.hip.add = importdata([pathMuscleAnalysis,'subject01_MuscleAnalysis_MomentArm_hip_add_',side,'.sto']);
    % hip rotation
    MA.hip.rot = importdata([pathMuscleAnalysis,'subject01_MuscleAnalysis_MomentArm_hip_rot_',side,'.sto']);
    % knee flexion
    MA.knee.flex = importdata([pathMuscleAnalysis,'subject01_MuscleAnalysis_MomentArm_knee_flex_',side,'.sto']);
    % ankle flexion
    MA.ankle.flex = importdata([pathMuscleAnalysis,'subject01_MuscleAnalysis_MomentArm_ankle_flex_',side,'.sto']);
    % subtalar
    MA.sub = importdata([pathMuscleAnalysis,'subject01_MuscleAnalysis_MomentArm_subt_angle_',side,'.sto']);

    %% Organize MuscleData
    if runPolynomialfit
        MuscleData.dof_names = dummy_motion.colheaders(Qs_idx.(side)); 
        muscleNames = {'glut_max1','glut_max2','glut_max3','glut_med1',...
        'glut_med2','glut_med3','glut_min1','glut_min2',...
        'glut_min3','add_long','add_brev','add_mag1','add_mag2',...
        'add_mag3','pectineus','iliacus','psoas','quad_fem',...
        'gemellus','piri','TFL','gracilis','semimem','semiten',...
        'bi_fem_lh','bi_fem_sh','sartorius','rectus_fem',...
        'vas_med','vas_int','vas_lat','gas_med',...
        'gas_lat','soleus','tib_post','tib_ant','ext_dig',...
        'ext_hal','flex_dig','flex_hal','per_brev','per_long',...
        'per_tert'};        
        for m = 1:length(muscleNames)
            MuscleData.muscle_names{m} = [muscleNames{m},'_',side];
            MuscleData.lMT(:,m)     = lMT.data(:,strcmp(lMT.colheaders,[muscleNames{m},'_',side]));            % lMT    
            MuscleData.dM(:,m,1)    = MA.hip.flex.data(:,strcmp(lMT.colheaders,[muscleNames{m},'_',side]));    % hip_flex
            MuscleData.dM(:,m,2)    = MA.hip.add.data(:,strcmp(lMT.colheaders,[muscleNames{m},'_',side]));     % hip_add
            MuscleData.dM(:,m,3)    = MA.hip.rot.data(:,strcmp(lMT.colheaders,[muscleNames{m},'_',side]));     % hip_rot
            MuscleData.dM(:,m,4)    = MA.knee.flex.data(:,strcmp(lMT.colheaders,[muscleNames{m},'_',side]));   % knee
            MuscleData.dM(:,m,5)    = MA.ankle.flex.data(:,strcmp(lMT.colheaders,[muscleNames{m},'_',side]));  % ankle
            MuscleData.dM(:,m,6)    = MA.sub.data(:,strcmp(lMT.colheaders,[muscleNames{m},'_',side]));         % subt 
        end
        MuscleData.q = q;
        MuscleData.qdot = qdot;
    end

    %% Call PolynomialFit
    if runPolynomialfit
        if strcmp(side,'r')            
            [muscle_spanning_joint_INFO_r,MuscleInfo_r] = ...
                PolynomialFit(MuscleData);
            MuscleData_r = MuscleData;
        else
            [muscle_spanning_joint_INFO_l,MuscleInfo_l] = ...
                PolynomialFit(MuscleData);
            MuscleData_l = MuscleData;
        end
        if savePolynomials
            save([pathPolynomials,'MuscleData_',subject,'_',side],['MuscleData_',side]);
            save([pathPolynomials,'muscle_spanning_joint_INFO_',subject,'_',side],['muscle_spanning_joint_INFO_',side]);
            save([pathPolynomials,'MuscleInfo_',subject,'_',side],['MuscleInfo_',side]);
        end
    end

    %% Create CasADi functions
    import casadi.*
    % Order mobilities: hip_flex, hip_add, hip_rot, knee, ankle, subt
    load([pathPolynomials,'muscle_spanning_joint_INFO_',subject,'_',side,'.mat']);
    if strcmp(side,'r')  
        muscle_spanning_joint_INFO = muscle_spanning_joint_INFO_r;
        MuscleInfo = MuscleInfo_r;
    else
        muscle_spanning_joint_INFO = muscle_spanning_joint_INFO_l;
        MuscleInfo = MuscleInfo_l;
    end
    load([pathPolynomials,'MuscleInfo_',subject,'_',side,'.mat']);
    NMuscle = length(MuscleInfo.muscle);
    q_leg = 6;
    qin     = SX.sym('qin',1,q_leg);
    qdotin  = SX.sym('qdotin',1,q_leg);
    lMT     = SX(NMuscle,1);
    vMT     = SX(NMuscle,1);
    dM      = SX(NMuscle,q_leg);
    for i=1:NMuscle     
        index_dof_crossing  = find(muscle_spanning_joint_INFO(i,:)==1);
        order               = MuscleInfo.muscle(i).order;
        [mat,diff_mat_q]    = n_art_mat_3_cas_SX(qin(1,index_dof_crossing),order);
        lMT(i,1)            = mat*MuscleInfo.muscle(i).coeff;
        vMT(i,1)            = 0;
        dM(i,1:q_leg) = 0;
        nr_dof_crossing     = length(index_dof_crossing); 
        for dof_nr = 1:nr_dof_crossing
            dM(i,index_dof_crossing(dof_nr)) = (-(diff_mat_q(:,dof_nr)))'*MuscleInfo.muscle(i).coeff;
            vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing(dof_nr))*qdotin(1,index_dof_crossing(dof_nr)));
        end 
    end
    f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM});

    %% Check results
    load([pathPolynomials,'MuscleData_',subject,'_',side,'.mat']);
    if strcmp(side,'r')  
        MuscleData = MuscleData_r;
    else
        MuscleData = MuscleData_l;
    end
    lMT_out = zeros(size(q,1),NMuscle);
    vMT_out = zeros(size(q,1),NMuscle);
    dM_out = zeros(size(q,1),NMuscle,q_leg);
    for i = 1:size(q,1)
        [out1,out2,out3] = f_lMT_vMT_dM(MuscleData.q(i,:),MuscleData.qdot(i,:));
        lMT_out(i,:) = full(out1);
        vMT_out(i,:) = full(out2);
        dM_out(i,:,1) = full(out3(:,1));
        dM_out(i,:,2) = full(out3(:,2));
        dM_out(i,:,3) = full(out3(:,3));
        dM_out(i,:,4) = full(out3(:,4));
        dM_out(i,:,5) = full(out3(:,5));   
        dM_out(i,:,6) = full(out3(:,6));
    end

    %% lMT
    % right
    figure()
    subplot(4,4,1)
    scatter(MuscleData.q(:,4),lMT_out(:,26)); hold on;
    scatter(MuscleData.q(:,4),MuscleData.lMT(:,26));
    xlabel('q knee');
    title('BFSH');
    subplot(4,4,2)
    scatter(MuscleData.q(:,4),lMT_out(:,29)); hold on;
    scatter(MuscleData.q(:,4),MuscleData.lMT(:,29));
    xlabel('q knee');
    title('VM');
    subplot(4,4,3)
    scatter(MuscleData.q(:,4),lMT_out(:,30)); hold on;
    scatter(MuscleData.q(:,4),MuscleData.lMT(:,30));
    xlabel('q knee');
    title('VI');
    subplot(4,4,4)
    scatter(MuscleData.q(:,4),lMT_out(:,31)); hold on;
    scatter(MuscleData.q(:,4),MuscleData.lMT(:,31));
    xlabel('q knee');
    title('VL');
    subplot(4,4,5)
    scatter3(MuscleData.q(:,4),MuscleData.q(:,5),lMT_out(:,32)); hold on;
    scatter3(MuscleData.q(:,4),MuscleData.q(:,5),MuscleData.lMT(:,32));
    xlabel('q knee');
    ylabel('q ankle');
    title('GM');
    subplot(4,4,6)
    scatter3(MuscleData.q(:,4),MuscleData.q(:,5),lMT_out(:,33)); hold on;
    scatter3(MuscleData.q(:,4),MuscleData.q(:,5),MuscleData.lMT(:,33));
    xlabel('q knee');
    ylabel('q ankle');
    title('GL');
    subplot(4,4,7)
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,39)); hold on;
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,39));
    xlabel('q knee');
    ylabel('q ankle');
    title('Flex dig');
    subplot(4,4,8)
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,40)); hold on;
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,40));
    xlabel('q knee');
    ylabel('q ankle');
    title('Flex hal');
    subplot(4,4,9)
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,36)); hold on;
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,36));
    xlabel('q knee');
    ylabel('q ankle');
    title('TA');
    subplot(4,4,10)
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,41)); hold on;
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,41));
    xlabel('q knee');
    ylabel('q ankle');
    title('Per brev');
    subplot(4,4,11)
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,42)); hold on;
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,42));
    xlabel('q knee');
    ylabel('q ankle');
    title('Per long');
    subplot(4,4,12)
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,43)); hold on;
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,43));
    xlabel('q knee');
    ylabel('q ankle');
    title('Per ter');
    subplot(4,4,13)
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,37)); hold on;
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,37));
    xlabel('q knee');
    ylabel('q ankle');
    title('Ext dig');
    subplot(4,4,14)
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),lMT_out(:,38)); hold on;
    scatter3(MuscleData.q(:,5),MuscleData.q(:,6),MuscleData.lMT(:,38));
    xlabel('q knee');
    ylabel('q ankle');
    title('Ext hal');
    legend('Polynomial','Model');
    suptitle('lMT right');

    %% Assert results
    assertLMT = zeros(size(q,1),NMuscle);
    for i = 1:NMuscle  
        assertLMT(:,i) = abs(lMT_out(:,i) - MuscleData.lMT(:,i));
        assertdM.hip.flex(:,i) = abs(dM_out(:,i,1) - MuscleData.dM(:,i,1));
        assertdM.hip.add(:,i) = abs(dM_out(:,i,2) - MuscleData.dM(:,i,2));
        assertdM.hip.rot(:,i) = abs(dM_out(:,i,3) - MuscleData.dM(:,i,3));
        assertdM.knee(:,i) = abs(dM_out(:,i,4) - MuscleData.dM(:,i,4));
        assertdM.ankle(:,i) = abs(dM_out(:,i,5) - MuscleData.dM(:,i,5));
        assertdM.sub(:,i) = abs(dM_out(:,i,6) - MuscleData.dM(:,i,6));
    end

    assertLMTmax = max(max(assertLMT));
    assertdM.hip.flexmax = max(max(assertdM.hip.flex));
    assertdM.hip.addmax = max(max(assertdM.hip.add));
    assertdM.hip.rotmax = max(max(assertdM.hip.rot));
    assertdM.kneemax = max(max(assertdM.knee));
    assertdM.anklemax = max(max(assertdM.ankle));

end
 