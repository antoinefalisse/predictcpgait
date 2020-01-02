%%  Muscle-tendon parameter optimal estimation
% This script formulates and solves the optimal control problem underlying
% the parameter estimation
%
% Author: Antoine Falisse and Lorenzo Pitto
% Date: 1/2/2020
%
function [params,activations,residuals,forces,scale,lMtilde,time,...
    stats,extlMT]=...
    parameterEstimation_ocp(JointMom,MuscMomArm,Fmax,EMGopt,indEMG2misc,...
    ParametersInit,MuscTenLen,MuscTenVel,SelectedMusclesInfo,modelName,...
    muscleNames,sides,Fvparam,Fpparam,Faparam,lMtilde_ext,W,...
    clinicalReport,idx_spanning,NindEMG2misc,idx_EMGscale2misc,...
    SelectedMusclesInfoSelSpanning,SelectedMusclesInfoOptSpanning,...
    SelectedMusclesInfoSpanningSel,timeOpt,side_t,dev_p,EMGScale_min)

%% get muscle-tendon lengths at extremes of the range of motion
Ntr = length(JointMom);
for s = 1:length(sides)
    NMuscleAll = length(muscleNames.(sides{s}).s);
    LMTMAX.(sides{s}).s = zeros(1,NMuscleAll);
    LMTMIN.(sides{s}).s = 100*ones(1,NMuscleAll);
end
for t=1:Ntr
    % find the max/min MT-lengths to constrain the fiber lengths   
    LMTMAX.(side_t{t}).s(1,idx_spanning{t}==1) = ...
        max([LMTMAX.(side_t{t}).s(1,idx_spanning{t}==1);MuscTenLen{t}],[],1);
    LMTMIN.(side_t{t}).s(1,idx_spanning{t}==1) = ...
        min([LMTMIN.(side_t{t}).s(1,idx_spanning{t}==1);MuscTenLen{t}],[],1);
end
xlsF = clinicalReport.xlsF;
xlsD = clinicalReport.xlsD;  
for s = 1:length(sides)
    Side = 'L';
    if strcmp(sides{s},'r')
        Side = 'R';
    end
    Values = readClinicalExam([xlsD,xlsF],Side);
    [LMT,~,NmuscCE] = getROM(modelName,Values,[],[],0,Side);
    for q = 1:length(muscleNames.(sides{s}).s)
        Im(q) = find(strcmpi(muscleNames.(sides{s}).s{q},NmuscCE));
    end
    LMTce.(sides{s}).s = LMT(Im)';
end
for s = 1:length(sides)
    figure,hold all
    bar([LMTMAX.(sides{s}).s;LMTMIN.(sides{s}).s;LMTce.(sides{s}).s]')
    set(gca,'XTickLabels',strrep(muscleNames.(sides{s}).s,'_',' '))
    set(gca,'XTick',1:length(muscleNames.(sides{s}).s))
    set(gca,'XTickLabelRotation',60)
    title('Muscle-tendon lengths used for estmating the parameters')
    legend({'MAX from gait/ISA','MIN from gait/ISA','MAX from clinic exam'})
    % update the upper bounds with the values from the clinical exam
    LMTMAX.(sides{s}).s = max([LMTMAX.(sides{s}).s;LMTce.(sides{s}).s],[],1);
    idx_lMopt_max.(sides{s}).s = ...
        find((LMTce.(sides{s}).s>0) & SelectedMusclesInfo{1}.Bool);
end

%% Formulate problem    
% Values for bounds
ParMIN = 0.5;
ParMAX = 2;
ScaleMAX = 5;
lMtildeMax = lMtilde_ext.max;
lMtildeMin = lMtilde_ext.min;
% Set up CasADi/Opti
import casadi.*
opti = casadi.Opti(); % Create opti instance
% EMG scale factors
for s = 1:length(sides)
    scale.(sides{s}).s = opti.variable(NindEMG2misc);
    opti.subject_to(EMGScale_min < scale.(sides{s}).s < ScaleMAX); 
    opti.set_initial(scale.(sides{s}).s,0);
    % Tendon slack lengths
    lTs.(sides{s}).s = opti.variable(length(SelectedMusclesInfo{1}.Index));
    opti.subject_to(...
        ParametersInit.(sides{s}).TenSlackLen(SelectedMusclesInfo{1}.Index)' ...
        *ParMIN < lTs.(sides{s}).s < ...
        ParametersInit.(sides{s}).TenSlackLen(SelectedMusclesInfo{1}.Index)' ...
        *ParMAX);
    opti.set_initial(lTs.(sides{s}).s,...
        ParametersInit.(sides{s}).TenSlackLen(SelectedMusclesInfo{1}.Index)');
    % Optimal fiber lengths
    lMopt.(sides{s}).s = opti.variable(length(SelectedMusclesInfo{1}.Index));
    opti.subject_to(...
        ParametersInit.(sides{s}).OptFibLen(SelectedMusclesInfo{1}.Index)' ...
        *ParMIN < lMopt.(sides{s}).s < ...
        ParametersInit.(sides{s}).OptFibLen(SelectedMusclesInfo{1}.Index)' ...
        *ParMAX); 
    opti.set_initial(lMopt.(sides{s}).s,...
        ParametersInit.(sides{s}).OptFibLen(SelectedMusclesInfo{1}.Index)');
end
% The parameters from both legs cannot deviate by more than +/- dev_p/100%
if length(sides) > 1
    opti.subject_to(1 - dev_p/100 < lTs.(sides{1}).s./lTs.(sides{2}).s < ...
        1 + dev_p/100);
    opti.subject_to(1 - dev_p/100 < lMopt.(sides{1}).s./lMopt.(sides{2}).s < ...
        1 + dev_p/100);
end

% Collocation scheme
d = 3; % degree of interpolating polynomial
method = 'radau'; % other option is 'legendre' (collocation scheme)
[tau_root,C,D,~] = CollocationScheme(d,method);
% Problem bounds
a_min = 0; a_max = 1; % bounds on muscle activation
aT_min = -1; aT_max = 1; % bounds on reserve actuators
vA_min = -1/100; vA_max = 1/100; % bounds on derivative of muscle activation
F_min = 0; F_max = 5; % bounds on normalized tendon force
dF_min = -100; dF_max = 100; % bounds on derivative of normalized tendon force
% Initial components of cost
C1 = 0;
C2 = 0;
C3 = 0;
C4 = 0;
C5 = 0;
C6 = 0;
C7 = 0;
% Loop over trials
Ntr = length(JointMom);
for t=1:Ntr
    N = size(JointMom{t},1); % number of frames
    Nmusc = size(MuscMomArm{t},2); % number of spanning muscles
    % number of spanning muscles whose parameters are optimized
    NselMusc = length(SelectedMusclesInfoSpanningSel{t}.Index); 
    Nj = size(JointMom{t},2); % number of joints
    lMT = MuscTenLen{t};
    vMT = MuscTenVel{t};
    aTendon = 35*ones(length(SelectedMusclesInfoSpanningSel{t}.Index),1); 
    shift = zeros(length(SelectedMusclesInfoSpanningSel{t}.Index),1);    
    % States
    % activations
    a{t} = opti.variable(Nmusc,N+1);
    amesh{t} = opti.variable(Nmusc,d*N);   
    opti.subject_to(a_min <= a{t} <= a_max);
    opti.subject_to(a_min < amesh{t} < a_max);
    opti.set_initial(a{t},0.2);
    opti.set_initial(amesh{t},0.2);
    % Muscle-tendon forces
    FTtilde{t} = opti.variable(NselMusc,N+1);
    FTtildemesh{t} = opti.variable(NselMusc,d*N);
    opti.subject_to(F_min < FTtilde{t} < F_max);
    opti.subject_to(F_min < FTtildemesh{t} < F_max);
    opti.set_initial(FTtilde{t},0.2);
    opti.set_initial(FTtildemesh{t},0.2);    
    % Controls
    % time derivative of muscle activations (states)
    vA{t} = opti.variable(Nmusc,N);
    tact = 0.015; tdeact = 0.06;
    opti.subject_to(vA_min/tdeact < vA{t} < vA_max/tact);
    opti.set_initial(vA{t},0);
    scaling.vA = 100; % Scaling factor: derivative muscle activations
    % Reserve actuators
    aT{t} = opti.variable(Nj,N);
    opti.subject_to(aT_min < aT{t} < aT_max);
    opti.set_initial(aT{t},0);
    scaling.aT = 150;
    % time derivative of muscle-tendon forces (states)
    dFTtilde{t} = opti.variable(NselMusc,N);
    opti.subject_to(dF_min < dFTtilde{t} < dF_max);
    opti.set_initial(dFTtilde{t},0.01);
    scaling.dFTtilde = 10; % Scaling factor: derivative muscle-tendon forces
    % Loop over mesh points formulating NLP
    h = (timeOpt{t}(2)-timeOpt{t}(1))/N;
    for k=1:N
        % Variables within current mesh interval
        ak = a{t}(:,k); FTtildek = FTtilde{t}(:,k);
        ak_colloc = [ak amesh{t}(:,(k-1)*d+1:k*d)]; 
        FTtildek_colloc = [FTtildek FTtildemesh{t}(:,(k-1)*d+1:k*d)];
        dFTtildek = dFTtilde{t}(:,k); aTk = aT{t}(:,k); vAk = vA{t}(:,k);
        % Loop over collocation points
        for j=1:d
            % Expression of the state derivatives at the collocation points
            ap = ak_colloc*C(:,j+1);
            FTtildep = FTtildek_colloc*C(:,j+1);        
            % Append collocation equations
            % Activation dynamics (explicit formulation)
            opti.subject_to(h*vAk.*scaling.vA - ap == 0);
            % Contraction dynamics (implicit formulation)
            opti.subject_to(h*dFTtildek.*scaling.dFTtilde - FTtildep == 0)
            % Add contribution to the quadrature function
            C1 = C1 + sumsqr(ak_colloc(:,j+1))*h;
            C2 = C2 + sumsqr(aTk)*h;
            C3 = C3 + sumsqr(vAk)*h;
            C4 = C4 + sumsqr(dFTtildek)*h;
        end      
        % State continuity at mesh transition
        opti.subject_to(a{t}(:,k+1)== ak_colloc*D);
        opti.subject_to(FTtilde{t}(:,k+1) == FTtildek_colloc*D);
        % Constraints
        % Get muscle-tendon forces and derive Hill-equilibrium
        % For some muscles, muscle force is activation
        FT = MX(Nmusc,1);
        idx_noFLV = 1:Nmusc;
        for tt = 1:length(SelectedMusclesInfoSpanningSel{t}.Index)
           idx_noFLV(idx_noFLV==SelectedMusclesInfoSpanningSel{t}.Index(tt))=[];
        end         
        if ~isempty(idx_noFLV)            
            FT(idx_noFLV) = ak(idx_noFLV).*(Fmax{t}(idx_noFLV))';
        end 
        % For other muscles, muscle force depends on FLV curves
        for m = 1:length(SelectedMusclesInfoSpanningSel{t}.Index) 
            [Hilldiff, FT(SelectedMusclesInfoSpanningSel{t}.Index(m)),~] = ...
                parameterEstimation_ForceEquilibrium_FtildeState(...
                ak(SelectedMusclesInfoSpanningSel{t}.Index(m)),...                 
                FTtildek(m),...                    
                dFTtildek(m)*scaling.dFTtilde,...                 
                lMT(k,SelectedMusclesInfoSpanningSel{t}.Index(m)),...
                vMT(k,SelectedMusclesInfoSpanningSel{t}.Index(m)),...
                Fmax{t}(SelectedMusclesInfoSpanningSel{t}.Index(m)),...
                lMopt.(side_t{t}).s(...
                    SelectedMusclesInfoSelSpanning{t}.Index(m)),...
                lTs.(side_t{t}).s(...
                    SelectedMusclesInfoSelSpanning{t}.Index(m)),...
                ParametersInit.(side_t{t}).PennAng(...
                SelectedMusclesInfoOptSpanning{t}.Index(m)),...
                Fvparam,Fpparam,Faparam,aTendon(m),shift(m)); 
            % Hill-equilibrium constraint
            opti.subject_to(Hilldiff == 0);
        end 
        % Joint torque is sum of muscle torque and reserve actuator
        for j = 1:Nj
            IDdiff = JointMom{t}(k,j) - (sum(FT.*MuscMomArm{t}(k,:,j)') + ...
                aTk(j)*scaling.aT);
            opti.subject_to(IDdiff == 0);        
        end 
        % EMG-driven activations
        for m=1:size(indEMG2misc{t},1)
            % get activation from muscle for which we have EMGs
            Act = ak(indEMG2misc{t}(m,3));
            % get EMG        
            EMG = EMGopt{t}(k,indEMG2misc{t}(m,1));      
            % weighted squared difference EMG vs act
            actEMGdiff = Act-EMG.*scale.(side_t{t}).s(idx_EMGscale2misc{t}(m));
            C5 = C5 + sumsqr(actEMGdiff)*W.EMG.ind(idx_EMGscale2misc{t}(m));
            if indEMG2misc{t}(m,1) ~= 1 % not for REF
                opti.subject_to(-0.1 <= actEMGdiff <= 0.1);  
            end
        end  
        % Activation dynamics constraints
        act1 = vAk*scaling.vA + ak/tdeact;
        opti.subject_to(act1 >= 0);
        act2 = vAk*scaling.vA + ak/tact;
        opti.subject_to(act2 <= 1/tact);
    end  
end
% get muscle fiber lengths at the maximum muscle-tendon lengths assuming
% nul muscle-tendon velocities, nul muscle activations, and nul dfdt
NselMusc_all = length(SelectedMusclesInfo{1}.Index);
for s = 1:length(sides)    
    F_lMtilde_max{s} = opti.variable(NselMusc_all,1);
    opti.subject_to(F_min <= F_lMtilde_max{s} <= F_max);
    opti.set_initial(F_lMtilde_max{s},0.2);
    F_lMtilde_min{s} = opti.variable(NselMusc_all,1);
    opti.subject_to(F_min <= F_lMtilde_min{s} <= F_max);
    opti.set_initial(F_lMtilde_min{s},0.2);    
    idx_side = find(strcmp(side_t,sides{s}),1,'first');    
    LMTMAX_side_temp = LMTMAX.(sides{s}).s(:,SelectedMusclesInfo{1}.Index);
    LMTMIN_side_temp = LMTMIN.(sides{s}).s(:,SelectedMusclesInfo{1}.Index);
    lMopt_side_temp = lMopt.(sides{s}).s(...
        SelectedMusclesInfoSelSpanning{1}.Index);
    lTs_side_temp = lTs.(sides{s}).s(SelectedMusclesInfoSelSpanning{1}.Index);
    Fmax_side_temp = Fmax{idx_side}(SelectedMusclesInfoSpanningSel{1}.Index);
    PennAngle_side_temp = ParametersInit.(sides{s}).PennAng(...
        SelectedMusclesInfoOptSpanning{1}.Index);    
    Lfibmax.(sides{s}).s = MX(1,NMuscleAll);
    Lfibmin.(sides{s}).s = MX(1,NselMusc_all);
for m = 1:NselMusc_all
    % Upper extreme
    [Hilldiff_lMtilde_max, ~, ...
        Lfibmax.(sides{s}).s(1,SelectedMusclesInfo{1}.Index(m))] = ...
        ForceEquilibrium_FtildeState_ParamEst(0,F_lMtilde_max{s}(m),0,...
        LMTMAX_side_temp(m),0,Fmax_side_temp(m),lMopt_side_temp(m),...
        lTs_side_temp(m),PennAngle_side_temp(m),Fvparam,Fpparam,Faparam,35,0); 
    % Lower extreme
    [Hilldiff_lMtilde_min, ~, Lfibmin.(sides{s}).s(m)] = ...
        ForceEquilibrium_FtildeState_ParamEst(0,F_lMtilde_min{s}(m),0,...
        LMTMIN_side_temp(m),0,Fmax_side_temp(m),lMopt_side_temp(m),...
        lTs_side_temp(m),PennAngle_side_temp(m),Fvparam,Fpparam,Faparam,35,0);
    % Hill-equilibrium constraint
    opti.subject_to(Hilldiff_lMtilde_max == 0);
    opti.subject_to(Hilldiff_lMtilde_min == 0);
end 
% deviation max fiber length to lMtildeMax  
C6 = C6 + sumsqr(Lfibmax.(sides{s}).s(1,idx_lMopt_max.(sides{s}).s)-lMtildeMax); 
% Constraints
opti.subject_to(...
    Lfibmax.(sides{s}).s(1,SelectedMusclesInfo{1}.Index)-lMtildeMax < 0);
opti.subject_to(-Lfibmax.(sides{s}).s(1,SelectedMusclesInfo{1}.Index)+1< 0);
opti.subject_to(Lfibmin.(sides{s}).s-1< 0);
opti.subject_to(-Lfibmin.(sides{s}).s+lMtildeMin< 0);
% Penalize the deviation for the linearly-scaled parameters
C7 = C7 + sum(lMopt.(side_t{t}).s);
end
J = C1*W.a+C2*W.aT+C3*W.vA+C4*W.dF+C5*W.EMG.all+C6*W.lMopt_max+W.lMopt*C7;
opti.minimize(J);
% Settings
optionssol.ipopt.hessian_approximation = 'limited-memory';
optionssol.ipopt.mu_strategy = 'adaptive';
optionssol.ipopt.nlp_scaling_method = 'gradient-based';
optionssol.ipopt.linear_solver = 'mumps';
optionssol.ipopt.tol = 1e-6;
optionssol.ipopt.max_iter = 10000;
% Solver
opti.solver('ipopt',optionssol);
sol = opti.solve();
stats = opti.stats();
% Retrieve solutions
for s = 1:length(sides)
    lMopt_opt.(sides{s}) = sol.value(lMopt.(sides{s}).s);
    lTs_opt.(sides{s}) = sol.value(lTs.(sides{s}).s);
    params.(sides{s}) = ParametersInit.(sides{s});
    params.(sides{s}).OptFibLen(1,SelectedMusclesInfo{1}.Index) = ...
        lMopt_opt.(sides{s});
    params.(sides{s}).TenSlackLen(1,SelectedMusclesInfo{1}.Index) = ...
        lTs_opt.(sides{s});
    scale.(sides{s}) = sol.value(scale.(sides{s}).s);    
end
for t=1:Ntr
    activations{t} = sol.value(a{t})';
    residuals{t} = sol.value(aT{t})'*scaling.aT;
    tgrid = linspace(timeOpt{t}(1),timeOpt{t}(2),N+1);
    dtime = zeros(1,d+1);
    for i = 1:d+1
        dtime(i)=tau_root(i)*((timeOpt{t}(2)-timeOpt{t}(1))/N);
    end
    time{t} = tgrid;
    forces{t} = sol.value(FTtilde{t})'.*repmat(...
        Fmax{t}(SelectedMusclesInfoSpanningSel{t}.Index),...
        size(sol.value(FTtilde{t})',1),1);
    lMTinterp = MuscTenLen{t}(:,SelectedMusclesInfoSpanningSel{t}.Index);    
    params_temp = zeros(4,length(SelectedMusclesInfoSelSpanning{t}.Index));   
    params_temp(2,:) = ...
        lMopt_opt.(side_t{t})(SelectedMusclesInfoSelSpanning{t}.Index);
    params_temp(3,:) = ...
        lTs_opt.(side_t{t})(SelectedMusclesInfoSelSpanning{t}.Index);
    params_temp(4,:) = ParametersInit.(side_t{t}).PennAng(...
        SelectedMusclesInfoOptSpanning{t}.Index);
    aTendon = 35; % default
    shift = 0;       
    Ftemp = sol.value(FTtilde{t})';
    [~,lMtilde{t}] = FiberLength_TendonForce_tendon(Ftemp(1:N,:),...
        params_temp,lMTinterp,aTendon,shift);
end  
for s = 1:length(sides)
    extlMT.LMTMAX.(sides{s}) = LMTMAX.(sides{s}).s;
    extlMT.LMTMIN.(sides{s}) = LMTMIN.(sides{s}).s;
    extlMT.idx_lMopt_max.(sides{s}) = idx_lMopt_max.(sides{s}).s;
    extlMT.SelectedMusclesInfo.(sides{s}) = SelectedMusclesInfo{1}.Index;
    extlMT.F_lMtilde_max_opt.(sides{s}) = sol.value(F_lMtilde_max{s});
    extlMT.F_lMtilde_min_opt.(sides{s}) = sol.value(F_lMtilde_min{s});
    extlMT.LMTce.(sides{s}) = LMTce.(sides{s}).s;
end
end
