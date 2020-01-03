function output = forceModel_end(input)

Nsegment = input.auxdata.Nsegment;
obj = 0;
for ms = 1:Nsegment
    q = input.phase(ms).integral;
    obj = obj + q;
end
output.objective = obj;

NMuscles_spas = input.auxdata.NMuscles_spas;
sum_eFf_edFf = zeros(size(input.phase(ms).finalstate(:,1),1),...
    Nsegment*NMuscles_spas);
sum_eFfi_edFfi = zeros(size(input.phase(ms).initialstate(:,1),1),...
    Nsegment*NMuscles_spas);

min_EMG_ordered_beg_mat = input.auxdata.min_EMG_ordered_beg_mat;

for ms = 1:Nsegment  
    
    bc = min_EMG_ordered_beg_mat(ms,:);
    
    eFf = input.phase(ms).finalstate(:,1:NMuscles_spas);
    edFf = input.phase(ms).finalstate(:,NMuscles_spas+1:2*NMuscles_spas);
    sum_eFf_edFf(:,(ms-1)*NMuscles_spas+1:ms*NMuscles_spas) = ...
        eFf + edFf + bc;
    
    eFfi = input.phase(ms).initialstate(:,1:NMuscles_spas);
    edFfi = input.phase(ms).initialstate(:,...
        NMuscles_spas+1:2*NMuscles_spas);
    sum_eFfi_edFfi(:,(ms-1)*NMuscles_spas+1:ms*NMuscles_spas) = ...
        eFfi + edFfi + bc;    
end
output.eventgroup.event = [sum_eFfi_edFfi,sum_eFf_edFf];
