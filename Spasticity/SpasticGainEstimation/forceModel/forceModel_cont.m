function phaseout = forceModel_cont(input)

% Get phase-independent input data
Nsegment        = input.auxdata.Nsegment;
NMuscles_spas   = input.auxdata.NMuscles_spas;
Ff_td           = input.auxdata.Ff_td;
dFf_td          = input.auxdata.dFf_td;
threshold_Ff    = input.auxdata.threshold_Ff;
threshold_dFf   = input.auxdata.threshold_dFf;
bspas           = input.auxdata.bspas;
min_EMG_ordered_beg_mat = input.auxdata.min_EMG_ordered_beg_mat;

for ms = 1:Nsegment
    % Get phase-dependent input data 
    FTtildeSpline = input.auxdata.splinestruct(ms).MS.Ftilde;
    dFTtildeSpline = input.auxdata.splinestruct(ms).MS.dFtilde;
    EMGSpline = input.auxdata.splinestruct(ms).MS.EMG;  

    % Get controls
    FTtilde = input.phase(ms).control(:,1:NMuscles_spas);
    dFTtilde = input.phase(ms).control(:,NMuscles_spas+1:end);

    % Get states
    eFf = input.phase(ms).state(:,1:NMuscles_spas);
    edFf = input.phase(ms).state(:,NMuscles_spas+1:end);
    
    numColPoints = size(input.phase(ms).state,1);
    % Get parameters
    gFf = input.phase(ms).parameter(:,1:NMuscles_spas);
    gdFf = input.phase(ms).parameter(:,NMuscles_spas+1:2*NMuscles_spas);

    % PATH CONSTRAINTS
    FTtildediff = FTtilde-FTtildeSpline;
    dFTtildediff = dFTtilde-dFTtildeSpline;
    min_EMG_ordered_beg_mat_rep = ...
        repmat(min_EMG_ordered_beg_mat(ms,:),numColPoints,1);
    EMGdiff = min_EMG_ordered_beg_mat_rep + eFf + edFf - EMGSpline; 

    phaseout(ms).path = [FTtildediff,dFTtildediff,EMGdiff];

    % DYNAMIC CONSTRAINTS
    % Spindle dynamics: force component
    deFfdt = zeros(numColPoints,NMuscles_spas);
    threshold_Ff_N = repmat(threshold_Ff(ms,:),numColPoints,1);
    for m = 1:NMuscles_spas
        deFfdt(:,m) = forceFDynamics(eFf(:,m),FTtilde(:,m),Ff_td,...
            gFf(:,m),bspas,threshold_Ff_N(:,m));
    end
    
    % Spindle dynamics: derivative of force component
    dedFfdt = zeros(numColPoints,NMuscles_spas);
    threshold_dFf_N = repmat(threshold_dFf(ms,:),numColPoints,1);
    for m = 1:NMuscles_spas
        dedFfdt(:,m) = dFdtFDynamics(edFf(:,m),dFTtilde(:,m),dFf_td,...
            gdFf(:,m),bspas,threshold_dFf_N(:,m));
    end
    
    phaseout(ms).dynamics = [deFfdt dedFfdt];

    % OBJECTIVE FUNCTION
    EMGdiffweight = EMGdiff;    
    EMGdiffnorm = sum(EMGdiffweight.^2,2); 
    
    phaseout(ms).integrand = EMGdiffnorm; 

end
end







